package AaronDriver;

use strict;

use lib $ENV{'AARON'};

use Constants qw(:THEORY :PHYSICAL :OTHER_USEFUL :COMPARE);
use AaronInit qw(%arg_in %arg_parser %template_info $template_job $system grab_cata_coords);
use AaronOutput qw(print_message $parent close_logfile print_ee);
use AaronTools;

use File::Path qw(make_path);
use Cwd qw(getcwd cwd);
use Exporter qw(import);
use Math::Trig;
use File::Path qw(make_path);
use File::Copy qw(cp);
use List::MoreUtils qw(firstidx);
use Data::Dumper;

our @EXPORT = qw(get_sub make_directories get_status run_stepX print_status analyze_result count_time);


my $status;
my $sub_ref;
my $constraints;
my @cat_coords;
my $thermo = {};

my $launch_failed = 0;
my $RT = BOLTZMANN * $arg_in{temperature};

#################
#ref to some useful function
################
my $imitate;
$imitate = sub {
    my ($working_dir, $f_file, $f_dir) = @_;
    my $start_dir = cwd;
    chdir($working_dir) or die "Unable to enter dir $working_dir: $!\n";
    opendir(DIR, '.') or die "Unable to open dir $working_dir: $!\n";
    my @names = readdir(DIR) or die "Unable to read $working_dir: $!\n";
    closedir(DIR);
    my $find_folder;
    foreach my $name (@names) {
        next if ($name eq '.');
        next if ($name eq '..');

        if (-d $name) {
            &$imitate($name, $f_file, $f_dir);
            $find_folder = 1;
            next;
        }
        &$f_file($name);
    }
    $f_dir && &$f_dir($find_folder);
    chdir($start_dir);
};  #end of $imitate

#Check for .sub file, to see if the substitute is needed
#if using this method will return a hash ref
#$sub_ref->{3-CH3-2-OH} => {
#                              3 => 'CH3',
#                              2 => 'OH'}
sub get_sub {
    my $sub_file = $arg_in{jobname} . ".sub";
    my $msg;
    my $sub_got;

    sub examine_sub {
        my ($sub) = @_;
        my %sub = %$sub;
        my $msg;
        if ($arg_in{catalyst}) {
            my @wrong_sub = grep {$_ >= $template_info{first_cat_atom}} keys %sub;
            if (@wrong_sub) {
                delete @sub{@wrong_sub};
                $msg .= "Since you require new catalyst mapping, ";
                $msg .= "substituent on the original catalyst is not allowed. ";
                $msg .= "Substitution at @wrong_sub are removed.\n";
            }
        }
        return $msg;
    }

    if (-e $sub_file) {
        $sub_ref = {};
        print "Make substituents\n";
        open SUB, "<$sub_file" or die "Can't open substitue";

        my %cat_re;
        my %sub_re;
        while(<SUB>) {
            /^$/ && do { next; };

            /^&cat/i && do {
                my @array = split /\s+/;
                shift (@array);
                my %cat_temp = map { my ($atom, $sub) = split '=';
                                     ($atom + $template_info{first_cat_atom}, $sub) } @array;
                my $geometry = join('-', map { "$_-$cat_temp{$_}" } sort keys %cat_temp);
                $cat_re{$geometry} = { %cat_temp };
                next;
            };

            /^&sub/i && do {
                my @array = split /\s+/;
                shift (@array);
                my %sub_temp = map { my ($atom, $sub) = split '=';
                                     ($atom, $sub) } @array;

                $msg = &examine_sub(\%sub_temp, $msg);

                my $geometry = join('-', map { "$_-$sub_temp{$_}" } sort keys %sub_temp);
                $sub_re{$geometry} = { %sub_temp };
                next;
            };

            my @array = split /\s+/;
            my $cat_idx = firstidx {$_ =~ /cat/i} @array;
            my $sub_idx = firstidx {$_ =~ /sub/i} @array;

            my @cat;
            my @sub;
            if ($cat_idx < $sub_idx) {
                @cat = @array[$cat_idx+1..$sub_idx-1];
                @sub = @array[$sub_idx+1..$#array];
            }else {
                @sub = @array[$sub_idx+1..$cat_idx-1];
                @cat = @array[$cat_idx+1..$#array];
            }
            my %cat = map { my ($atom, $sub) = split '=';
                            ($atom + $template_info{first_cat_atom}, $sub) } @cat;
            my %sub = map { my ($atom, $sub) = split '=';
                            ($atom, $sub) } @sub;

            $msg = &examine_sub(\%sub, $msg);

            my $geometry = '';
            if (%cat) {
                $geometry .= 'Cat-' . join('-', map{"$_-$cat{$_}"} sort keys %cat);
            }
            if (%sub) {
                if (%cat) {
                    $geometry .= '-Sub-';
                }else{
                    $geometry .= 'Sub-';
                }
                $geometry .= join('-', map{"$_-$sub{$_}"} sort keys %sub);
            } 

            $sub_ref->{$geometry}->{sub} = { %cat, %sub };

            #map new catalyst we need also define rotated conformers
            my @rot_sub;
            foreach my $sub_atom (keys %{ $sub_ref->{$geometry}->{sub} }) {
                if (sub_geo($sub_ref->{$geometry}->{sub}->{$sub_atom}, 'conf')) {
                    push (@rot_sub, $sub_atom);
                }
            }
            $sub_ref->{$geometry}->{conf} = [ @rot_sub ];
        }


        foreach my $cat (sort keys %cat_re) {
            foreach my $sub (sort keys %sub_re) {
                my $geometry = 'Cat-' . $cat . '-Sub-' . $sub;
                $sub_ref->{$geometry}->{sub} = { %{ $cat_re{$cat} }, %{ $sub_re{$sub} } };
                my @rot_sub;
                foreach my $sub_atom (keys %{ $sub_ref->{$geometry}->{sub} }) {
                    if (sub_geo($sub_ref->{$geometry}->{sub}->{$sub_atom}, 'conf')) {
                        push (@rot_sub, $sub_atom);
                    }
                }
                $sub_ref->{$geometry}->{conf} = [ @rot_sub ];
            }
        }

        foreach my $geometry (keys %$sub_ref) {
        }

        close SUB;
        $sub_got = 1;
    }elsif($arg_in{catalyst}) {
        $msg .= "New catalyst is being mattached, no substitution are undergoing now...\n";
        @cat_coords = grab_cata_coords();
        $sub_got = 0;
    }else{
        $msg .= "Please use -map_cata in the command line to explore a new catalyst.\n";
        $msg .= "Or make a <jobname>.sub file containing substitues you want to make to the \n";
        $msg .= "template. Otherwise I don't know what you want to do.\n";
        $msg .= "Exit without launching AARON\n";
        print_message($msg);
        close_logfile();
        exit(1);
    }

    if ($msg) {
        print_message($msg);
    }
    return $sub_got;
}


#initiate a a directory tree and status hashes
sub make_directories {
    my ($dir_type) = @_;

    if ($dir_type eq 'sub') {
        for my $sub (keys %{ $sub_ref }) {
            if (!(-d $sub)) {
                mkdir "$sub";
                $thermo->{$sub} = {};
            }
            &dir_tree($arg_in{TS_path_sub}, $sub);
        }
    }
    if ($dir_type eq 'map') {
        if (!(-d 'mapping')) {
            mkdir 'mapping';
            $thermo->{mapping} = {};
        }
        &dir_tree($arg_in{TS_path_map}, 'mapping');
    }
}


#mapping a directory tree and do any operation related to the specific file
sub dir_tree {
    my @new_dir;
    #top directory to mimic; top directory to make;
    #constraints for each gemotry;
    #status for each geometry;
    my ($top_dir_tar, $top_dir_make) = @_;

    if (!(-d $top_dir_tar)) {
        return;
    }

    chdir($top_dir_tar) or die "Unable to enter dir $top_dir_tar: $!\n";

    my $f_file = sub {
        my ($name) = @_;
        if ($name =~ /(\S+).xyz/) {
            my $extend = $1;
            my $tempdir = cwd;
            if ($tempdir =~ /$top_dir_tar(\S+)/) {
                
                my $newdir = $1 . '/' . $extend;
                #make distance hashes for each geometry
                my $key = $newdir;
                my @model_TS_coords = grab_coords($name);

                foreach my $constraint (@{ $template_info{constraints} }) {
                    my $distance = distance(@{ $constraint }, \@model_TS_coords);

                    unless (defined $constraints->{$constraint}) {
                        $constraints->{$constraint} = {};
                    }

                    $constraints->{$constraint}->{$key} = $distance;
                }
                push (@new_dir, $newdir);
            }
        }
    };

    &$imitate($top_dir_tar, $f_file);
    chdir(LAUNCH_DIR);

    #make paths and initialize eevery geometry to be on step0 1st attempt
    foreach my $newdir (@new_dir) {
        my $geometry = $top_dir_make . $newdir;

        unless (-d $geometry) {
            make_path($geometry);
        }
        &init_status($geometry);
    }

}   #end of dir_tree


#get status for each directories, the hash ref status is initiated in the dirtree
sub init_status {
    my ($head) = @_;

    my $current_dir = cwd;

    my $f_file = sub { return; };

    my $f_dir = sub { 
        my ($dir_found) = @_;

        if (! $dir_found) {
            my $dir = cwd;
            if ($dir =~ /($head.*)/) {
                my $geometry = $1;
                if (! $status->{$geometry}) {
                    $status->{$geometry} = {};
                    $status->{$geometry}->{step} = 1;
                    #This count how many attempt has try for this step
                    $status->{$geometry}->{attempt} = 1;
                    #This count how many cycles has been done 
                    $status->{$geometry}->{cycle} = 1;
                    $status->{$geometry}->{status} = 'start';
                }
            }else{
                print "Cannot initiate status hashses for $dir\n";
            }
        }
    };

    &$imitate($head, $f_file, $f_dir);
    chdir($current_dir);
}


sub get_status {

    foreach my $geometry (sort keys %{ $status }) {
        my $head = &get_head($geometry);

        if ($status->{$geometry}->{status} eq 'killed' || $status->{$geometry}->{status} eq 'finished') {
            next;
        }
        &check_steps($geometry);
    }
}


sub run_stepX {
    foreach my $geometry (sort keys %{ $status }) {

        if ($status->{$geometry}->{status} eq 'killed' || 
            $status->{$geometry}->{status} eq 'finished') {
            next;
        }

        chdir($geometry) or die "Cannot change into $geometry\n";

        my $step;
        my $coords;
        my $subs_atoms = [];

        my $new_com;

        my ($quit, $com_file);

        my $error = 'none';
        if ($status->{$geometry}->{status} eq 'start') {
            $new_com = 1;
            my $head = &get_head($geometry);
            my $file = $geometry;
            $file =~ s/$head\///;

            if ($sub_ref->{$head}) {
                $coords = [ grab_coords("$arg_in{TS_path_sub}/$file.xyz") ];
                foreach my $atomn(sort keys %{ $sub_ref->{$head}->{sub} }) {
                    my $sub = $sub_ref->{$head}->{sub}->{$atomn};
                    ($coords, my $sub_atoms) = substitute($atomn, $sub, $coords, 0);
                    push (@$subs_atoms, @$sub_atoms);
                }
            }else {
                $coords = [ grab_coords("$arg_in{TS_path_map}/$file.xyz") ];
                $coords = [ map_catalyst($coords, \@cat_coords, $template_info{key_atoms},
                                       $arg_in{key_atoms_cata}, $arg_in{bonds_RMSD},
                                       $arg_in{bonds_LJ}, $template_info{first_cat_atom}) ];
            }
            $step = $status->{$geometry}->{step};

        }elsif ($status->{$geometry}->{status} eq 'done') {
            $new_com = 1;
            $step = $status->{$geometry}->{step};
            my $filename = &get_file_name($geometry);
            $coords = [ grab_coords("$filename.$step.log") ];
            $step ++;
            
            $status->{$geometry}->{step} = $step;
            $status->{$geometry}->{attempt} = 1;

        }elsif ($status->{$geometry}->{status} eq 'failed') {
            $new_com = 1;
            $step = $status->{$geometry}->{step};
            my $filename = &get_file_name($geometry);

            $coords = grab_coords("$filename.$step.log") ?
                        [ grab_coords("$filename.$step.log") ] :
                        [ grab_coords("$filename.$step.com") ];

            $error = $status->{$geometry}->{error};

            $status->{$geometry}->{attempt} ++;
            if ($status->{$geometry}->{attempt} > 5) {
                $status->{$geometry}->{status} = 'killed';
                $status->{$geometry}->{msg} = 'killed because of too many attempts. ';
                next;
            }

        }elsif ($status->{$geometry}->{status} eq  '2submit') {
            my $step = $status->{$geometry}->{step};
            $status->{$geometry}->{status} = 'pending';
            my $filename = &get_file_name($geometry);
            $com_file = "$filename.$step.com";
        }

        if ($new_com) {
            ($quit, $com_file) = &build_com($geometry, $coords, $step, $error, $subs_atoms);
        }

        if ($quit) {
            close_logfile();
            chdir(LAUNCH_DIR);
            close_logfile();
            exit 0;
        }

        if ($com_file) {
            if ($status->{$geometry}->{step} == 1) {
                &submit_job($com_file, 1);
            }else {
                &submit_job($com_file);
            }
        }
        chdir(LAUNCH_DIR);
    }
} #End run_stepX;


sub check_steps {
    my ($geometry) = @_;

    my $head = &get_head($geometry);

    my $filename = &get_file_name($geometry);
   
    my $step = MAXSTEP;

    my $not_found = 1;

    while($step >= 0 && $not_found) {
        my $current = cwd;
        if (-e "$geometry/$filename.$step.log") {
            if (finished("$geometry/$filename.$step.log")) {
                $status->{$geometry}->{step} = $step;
                $status->{$geometry}->{status} = 'done';
                $not_found = 0;
            }elsif (findJob($geometry)) {
                $status->{$geometry}->{step} = $step;
                $status->{$geometry}->{status} = 'running';
                $not_found = 0;
            }else {
                my $error = get_error("$geometry/$filename.$step.log");
                $status->{$geometry}->{step} = $step;
                $status->{$geometry}->{status} = 'failed';
                $status->{$geometry}->{error} = $error;
                $not_found = 0;
            }
        }elsif (-e "$geometry/$filename.$step.com") {
            $status->{$geometry}->{step} = $step;
            $not_found = 0;
            if (findJob($geometry)) {
                $status->{$geometry}->{status} = 'pending';
                $status->{$geometry}->{step} = $step;
                $status->{$geometry}->{msg} = "waiting in queue";
            }else {
                $status->{$geometry}->{status} = '2submit';
                $status->{$geometry}->{step} = $step;
            }
        }
        $step--;
    }

    if ($status->{$geometry}->{step} >= 2) {
        #check if the reaction goes the correct way

        my $step = $status->{$geometry}->{step};
        my $go_wrong = &check_active_center($geometry);

        if ($go_wrong) {
            return 0;
        }
        #check if the geometry is a unique one
        my $geo_to_kill = &check_unique($geometry);
        #kill the repeated geometry
        if ($geo_to_kill) {
            $status->{$geo_to_kill}->{status} = 'killed';
            &kill($geo_to_kill);
            system("rm -fr $geo_to_kill");
            $status->{$geo_to_kill}->{msg} = 'killed because of repeated conformers';
            sleep(2);
            #if the killed geometry is the current checking one, reutrn 0
            if ($geo_to_kill eq $geometry) {
                return 0;
            }
        }

        if ($status->{$geometry}->{status} eq 'done') {
            my $conf_made = &make_conformers($geometry);
            unless ($conf_made) {
                if ($status->{$geometry}->{step} == 3) {
                    #optimization is done, making xyz file
                    my @coords = grab_coords("$geometry/$filename.4.log");
                    my $comment = "optimized geometry of $geometry";
                    printXYZ(\@coords, $comment, "$geometry/$filename.xyz");
                    my $failed_to_build = &build_sub_lib($geometry);
                    if ($failed_to_build) {
                        my $time = &count_time(30);
                        my $msg = "AARON failed to build TS library for substitution ".
                                  "for new catalyst $geometry\n".
                                  "AARON think this is a fatal problem, please check that manually\n".
                                  "AARON will sleep for 30min and continue\n".
                                  "AARON will restart on $time\n".
                                  "Sleeping...\n";
                        print_message($msg);
                        sleep(1800);
                    }
                }elsif ($status->{$geometry}->{step} >= 4) {
                    $thermo->{$head}->{$geometry} = [ get_thermo("$geometry/$filename.4.log", 
                                                                    $arg_in{temperature}) ];
                    if ($status->{$geometry}->{step} == 5) {
                        my $energy_h = get_energy("$geometry/$filename.5.log");
                        my @thermo_h = map {$energy_h + $_ - 
                            $thermo->{$head}->{$geometry}->[0]} @{ $thermo->{$head}->{$geometry} };

                        $thermo->{$head}->{$geometry} = [ @{$thermo->{$head}->{$geometry}},
                                                          @thermo_h ];
                        $status->{$geometry}->{status} = 'finished';
                    }
                }
            }
        }
    }
}


sub check_active_center {

    my ($geometry) = @_;

    my ($energy, $coords) = &get_energy_coords($geometry);

    unless($energy && $coords) {
        return;
    }

    foreach my $constraint (@{ $template_info{constraints} }) {
        foreach my $geo_constraint (sort keys %{ $constraints->{$constraint} }) {

            if ($geometry =~ /$geo_constraint/) {
                my $dist_ref = $constraints->{$constraint}->{$geo_constraint};
                my $dist = distance(@{ $constraint }, $coords);
                if ($dist - $dist_ref > CUTOFF->{D_CUTOFF}) {
                    &revert_to_step2($geometry, $constraint, 'reduce');
                    return 1;
                }elsif ($dist_ref - $dist > CUTOFF->{D_CUTOFF}) {
                    &revert_to_step2($geometry, $constraint, 'increase');
                    return 1;
                }
                next;
            }
        }
    }
    return 0;    
} 


sub revert_to_step2 {
    my ($geometry, $constraint, $option) = @_;
    my $filename = &get_file_name($geometry);
    &kill($geometry);
    print_message("Reverting to step2...\n");

    my @coords = grab_coords("$geometry/$filename.2.log");
    
    unless(@coords) {
        my $message = "Cannot find coordinates from the step2 log file, check that\n";
        $status->{$geometry}->{status} = 'killed';
        $status->{$geometry}->{mesg} = $message;
        print_message($message);
        return 0
    }

    chdir("$parent/$geometry") or die "Can't change into directory $parent/$geometry\n";

    #remove evering thing more than step 2
    foreach my $later_step (2..MAXSTEP) {
        if (-e "$filename.$later_step.com") {
            print_message("Removing $filename.$later_step.com...\n");
            if ($arg_in{debug}) {
                system("mv $filename.$later_step.com $filename.$later_step.prec");
            }else {
                system("rm -fr $filename.$later_step.com");
            }
        }

        if (-e "$filename.$later_step.log") {
            print_message("Removing $filename.$later_step.log...\n");
            if ($arg_in{debug}) {
                system("mv $filename.$later_step.log $filename.$later_step.prel");
            }else {
                system("rm -fr $filename.$later_step.log");
            }
        }

        if (-e "$filename.$later_step.log") {
            print_message("Removing $filename.$later_step.log...\n");
            system("rm -fr $filename.$later_step.com");
        }
        if (-e "$filename.$later_step.job") {
            print_message("Removing .job files for $later_step step...\n");
            system("rm -fr $filename.$later_step.job*");
        }
    }

    my $over_cycle;

    $status->{$geometry}->{step} = 2;
    $status->{$geometry}->{attempt} = 1;
    $status->{$geometry}->{cycle} ++;

    if ($status->{$geometry}->{cycle} > 5) {
        $status->{$geometry}->{status} = 'killed';
        $status->{$geometry}->{msg} = "killed becuse of too many recycles. ";
        $over_cycle = 1;
    }

    unless ($over_cycle) {
        $status->{$geometry}->{status} = "pending";
        $status->{$geometry}->{msg} = "revert to step 2, now waiting in the queue";

        #fix geometry -- make small changes to distances of constrained bond lengths from step 2
        my $distance = distance(@{ $constraint }, \@coords);
        if ($option eq 'reduce') {
            print_message( "Shortening the distance between @{ $constraint } by 0.1A\n");
            $distance -= 0.1;
        }elsif ($option eq 'increase') {
            print_message("Lengthing the distance between @{ $constraint } by 0.1A\n");
            $distance += 0.1;
        }else {
            my $random = 0.1 * (2 * int(rand(2)) - 1);
            print_message("Randomly move the distance by $random A\n");
            $distance += $random;
        }
        change_distance(@{ $constraint }, $distance, \@coords);
        sleep(2);

        my $com_file = &build_com($geometry, \@coords, 2);

        &submit_job($com_file);
    }

    chdir("$parent") or die "Can't return to $parent!\n";
    return 0;
}


sub check_unique {
    my($geometry) = @_;

    my $filename = &get_file_name($geometry); 

    my ($energy, $coords) = &get_energy_coords($geometry);

    #if no energy or coords, give up checking unique
    unless($energy && $coords) {
        return;
    }

    my @geos_conf = &get_confs($geometry);

    my $geo_to_kill;
    foreach my $geo_ref (@geos_conf) {
        if ($geo_ref eq $geometry) {
            next;
        }
        my ($energy_ref, $coords_ref) = &get_energy_coords($geo_ref);

        unless($energy_ref && $coords_ref) {
            next;
        }
        if (abs($energy_ref - $energy)*HART_TO_KCAL < CUTOFF->{E_CUTOFF}) {

            if (RMSD_align($coords, $coords_ref, CUTOFF->{RMSD_CUTOFF})) {
                $geo_to_kill = &determine_priority($geometry, $geo_ref);
            }
        }
    }
    return $geo_to_kill;
}

sub determine_priority {
    my ($geo1, $geo2) = @_;

    my $slow_geo;
    if ($status->{$geo1}->{step} > $status->{$geo2}->{step}) {
        $slow_geo = $geo2;
    }elsif ($status->{$geo1}->{step} < $status->{$geo2}->{step}) {
        $slow_geo = $geo1;
    }elsif (findJob($geo2)) {
        $slow_geo = $geo1;
    }

    return $slow_geo;
}


sub get_confs {

    my ($geometry) = @_;
    my @geos_conf;

    foreach my $geo (sort keys %{ $status }) {
        if ($geo =~ /(\S+)\/Cf\d+/) {
            my $main = $1;
            if ($geometry =~ /^$main\/Cf\d+/) {
                push (@geos_conf, $geo);
             }
        }
    }

    return @geos_conf;
}


sub kill {
    
    my ($geometry) = @_;

    foreach my $job (findJob("$parent/$geometry")) {
        &kill_job($job, $system);
    }
}


sub kill_job {
    
    my ($job, $system) = @_;

    if ($queue_type eq 'LSF') {
        system("bkill $job");
    }elsif ($queue_type eq 'PBS') {
        system("qdel $job");
    }

} #end kill_job;

#get filename from $geometry
sub get_file_name {
    my ($filename) = @_;
    $filename =~ s/\//\./g;
    return $filename;
}

sub get_head {
    my ($geometry) = @_;
    my @geo_part = split('/', $geometry);
    my $head = shift(@geo_part);
    return $head;
}


#get energy and coords for the geometry from the latest .log file
sub get_energy_coords {
    my ($geometry) = @_;

    my $filename = &get_file_name($geometry);

    my $step;
    if ($status->{$geometry}->{status} eq 'pending' || 
        $status->{$geometry}->{status} eq '2submit') {
        if ($status->{$geometry}->{step} == 2) {
            return;
        }else {
            $step = $status->{$geometry}->{step} - 1;
        }
    }else {
        $step = $status->{$geometry}->{step};
    }
   
    my $energy;
    my $coords;

    if (-e "$geometry/$filename.$step.log") {
        $energy = get_energy("$geometry/$filename.$step.log");
        $coords = [ grab_coords("$geometry/$filename.$step.log") ];
    }

    return($energy, $coords);
}


sub build_com {
    my ($geometry, $coords, $step, $error, $subs_atoms) = @_;

    my $filename = get_file_name($geometry);

    my $route;
    my $footer;

    my $method;
    my $high_method;
    my $gen_basis;

    #subroutine reduce maxstep from maxstep 
    sub reduce_maxstep {
        my ($route_temp) = @_ ;
        my $step_change;
        if ($route_temp =~ /maxstep=(\d+)/) {
            my $new_max_step = $1 - 2;
            if ($new_max_step > 0) {
                $route_temp =~ s/maxstep=$1/max=$new_max_step/;
                $step_change = 1;
            }
        }else {
            $route_temp =~ s/opt=\(/opt=\(maxstep=5,/;
            $step_change = 1;
        }

        return ($route_temp, $step_change);
    }

    if (! $arg_in{gen}) {
        $method = "$arg_in{method}/$arg_in{basis}";
        $high_method = "$arg_in{high_method}/$arg_in{basis}";
    }else {
        $method = "$arg_in{method}/gen";
        $high_method = "$arg_in{high_method}/gen";
        $gen_basis = "\@$arg_in{gen}" . "$arg_in{basis}/N";
    }

    SWITCH: {
        if ($step == 1) { $route = "\%chk=$filename.chk\n";
                          $route .= "#$arg_in{lowmethod} opt=(maxcyc=5000) nosym";
                          #add constrats to substrate and old part of catalyst
                          $coords = &add_constraints($coords, $subs_atoms);
                          last SWITCH; }

        if ($step == 2) { $route = "\%chk=$filename.chk\n";
                          $route .= "#$method opt=(modredundant,maxcyc=1000)";

                          foreach my $constraint (@{ $template_info{constraints} }) {
                              my @con = map {$_ + 1} @$constraint;
                              $footer .= "B $con[0] $con[1] F\n"
                          }
                          $footer .= "\n";
                          last SWITCH; }

        if ($step == 3) { $route = "\%chk=$filename.chk\n";
                          $route .= "#$method opt=(readfc,ts,maxcyc=1000)";
                          last SWITCH; }

        if ($step == 4) { $route = "\%chk=$filename.chk\n";
                          $route .= "#$method freq=(hpmodes,noraman,temperature=$arg_in{temperature})";
                          last SWITCH; }

        if ($step == 5) { $route = "\%chk=$filename.chk\n";
                          $route .= "#$high_method";
                          last SWITCH; }
    }

    if ($arg_in{solvent} ne "gas" && $step > 1) {
        $route .= " scrf=($arg_in{pcm},solvent=$arg_in{solvent})";
    }
    if ($arg_in{gen} && $step != 1) {
        $footer .= $gen_basis;
    }

    if ($error) {
        ERROR: {    
        if ($error eq 'CONV') {my $scf_change = $route =~ /scf=xqc/ ?
                                                $route .= " scf=xqc" : 0;
            
                                my $message = " SCF convergence problems with $geometry ";
                                if ($scf_change) {
                                    $message .= "...scf=xqc is in use\n";
                                }

                                $status->{$geometry}->{msg} = $message;
                                $status->{$geometry}->{status} = 'restart';
                                last ERROR;
                              }

        if ($error eq 'EIGEN') { $route =~ s/opt=\(/opt=\(noeigen,/;
                                 my $message = "Worng number of negative eigenvalues for $geometry ";
                                 $message .= "...Adding noeigen\n";
                                 $status->{$geometry}->{msg} = $message;
                                 $status->{$geometry}->{status} = 'restart';
                                 last ERROR;
                               }

        if ($error eq 'QUOTA') { if ( ! $arg_parser{no_quota}) {
                                    my $msg = "\nAARON thinks you hit the disk quota. "  .
                                           "Make more space or have quota increased. Then add " .
                                           "\"-quota\" in the command line to restart\n";
                                    print_message($msg);
                                    return 1;
                                 }else{
                                     last ERROR;
                                 }
                               }

        if ($error eq "CHK") { if (-e "$filename.chk") {
                                   unlink "$filename.chk";
                                   my $msg = "Problem with check point file, using calcfc\n";
                                   $route =~ s/readfc/calcfc/;
                                   $status->{$geometry}->{msg} = $msg;
                               }else {
                                   my $msg = "Errors detected with the $geometry, " .
                                             "This is a hard case to address, try to modified the structure " .
                                             "around manually, AARON will retry this step\n";
                                   $status->{$geometry}->{msg} = $msg;
                                   $status->{$geometry}->{status} = 'restart';
                                   return 0;
                               }
                               last ERROR;
                             }
        if ($error eq "CLASH") { my $msg = "Atoms too crowded in $filename.$step.com\n" .
                                           "Aaron has skipped this. If you want to continue this job " .
                                           "Please fix the geometry and restart AARON\n";
                                 $status->{$geometry}->{msg} = $msg;
                                 $status->{$geometry}->{status} = 'killed';
                                 
                                 return 0;
                               }
        if ($error eq "CHARGEMULT") { my $msg = "The combination of multipicity is not correct " .
                                                "AARON believe this is a fatal error so AARON quit " .
                                                "at this point, after you fix the problem, restart AARON\n";
                                      print_message($msg);
                                      return 1;
                                    }

        if ($error eq "REDUND") { my $msg = "Bend failed for some angles for $geometry, " .
                                            "find more details in $geometry.$step.log file\n" .
                                            "AARON has skipped this one, " .
                                            "If you want to continue this geometry after you fix the problem, " .
                                            "please restart AARON\n";

                                  $status->{$geometry}->{msg} = $msg;
                                  $status->{$geometry}->{status} = 'killed';
                                  return 0;
                                }

        if ($error eq "UNKNOWN") { my $msg = "$geometry stopped by unknown reason, " .
                                             "AARON will retry the failed step by changing some setup. ";
                                   $msg .= "Please also check $filename.$step.log manually\n";

                                   $status->{$geometry}->{msg} = $msg;
                                   $status->{$geometry}->{status} = 'restart';
                                   my $step = $status->{$geometry}->{step};
                                   system("rm -fr $filename.$step.*");
                                 }
        }
    }

    if ($status->{$geometry}->{cycle} > 1) {
        my $fc_change = ($route =~ s/readfc/calcfc/);
        ($route, my $step_change) = &reduce_maxstep($route);
        if ($fc_change) {
            $status->{$geometry}->{msg} .= "calculate fc instead of read fc from .chk file.\n";
            system("rm -fr $geometry.chk");
        }
        if ($step_change) {
            $status->{$geometry}->{msg} .= "using smaller step.\n";
        }
    }

    if ($status->{$geometry}->{attempt} > 3) {
        ($route, my $step_change) = &reduce_maxstep($route);
        if ($step_change) {
            $status->{$geometry}->{msg} .= "using smaller step.\n";
        }
    }

    my $com_file = "$filename.$step.com";
    my $comment = "$geometry step $step (attempt $status->{$geometry}->{attempt})";
    write_com($route, $comment, $arg_in{charge}, $arg_in{mult}, $coords,
              $footer, 1, $com_file);

    return (0, $com_file);
} #End build new com file

sub make_conformers {
    my ($geometry) = @_;
    my $head = &get_head($geometry);
    my $filename = &get_file_name($geometry);
    my $conf_made;

    if (my $conf = $sub_ref->{$head}->{conf}) {
        my $count = () = $geometry =~ /Cf/g;
        if ($count < @{ $conf }) {

            $conf_made = 1;
            my $atomn = $conf->[$count];
            my $sub  = $sub_ref->{$head}->{sub}->{$atomn};

            chdir("$parent/$geometry");

            my $step = $status->{$geometry}->{step};
            my @coords = grab_coords("$filename.$step.log");

            foreach my $cf_n(1..sub_geo($sub, 'num')) {
                my $cf_dir = "Cf$cf_n";
                my $degree = ($cf_n - 1) * sub_geo($sub, 'deg');
                my @coords_temp = copy_coords(\@coords);
                @coords_temp = sub_rotate($atomn, $degree, \@coords_temp);

                my $geometry_temp = $geometry . "/Cf$cf_n";
                $status->{$geometry_temp} = {};
                $status->{$geometry_temp}->{attempt} = 1;
                $status->{$geometry_temp}->{cycle} = 1;

                unless ( -d $cf_dir ) {
                    mkdir ($cf_dir) or die "cannot make $geometry\n";

                    if ($cf_n == 1) {

                        opendir (DIR, "$parent/$geometry") or die "cannot open\n";
                        my @files = readdir(DIR);
                        my $pattern = $geometry;
                        $pattern =~ s/\//\./g;

                        for (@files) {
                            /^\.+$/ && do {next;};
                            -d && do {next;};
                            /^$pattern(\S+)/ && do {
                                my $new_file = "$cf_dir/$pattern.$cf_dir" . $1;
                                system ("mv $_ $new_file");
                            };
                        }
                        chdir ($cf_dir);
                        if ($status->{$geometry}->{step} < MAXSTEP) {
                            $status->{$geometry_temp}->{step} = $status->{$geometry}->{step} + 1;
                            my ($quit, $com_file) = &build_com($geometry_temp, \@coords_temp, 
                                                        $status->{$geometry_temp}->{step});
                            &submit_job($com_file);
                            $status->{$geometry_temp}->{msg} = "New conformer waiting in queue";
                        }else {
                            $status->{$geometry_temp}->{step} = $status->{$geometry}->{step};
                            $status->{$geometry_temp}->{msg} = "Job finished but pending for more conformers";
                        }
                    }else {
                        chdir ($cf_dir);
                        $status->{$geometry_temp}->{step} = 2;
                        my ($quit, $com_file) = &build_com($geometry_temp, \@coords_temp, 2);
                        &submit_job($com_file);
                        $status->{$geometry_temp}->{msg} = "New conformer waiting in queue";
                    }
                    $status->{$geometry_temp}->{status} = 'pending';


                    chdir("$parent/$geometry") or die "cannot change into directory $parent/$geometry\n";
                }
            }

            $status->{$geometry}->{status} = 'new_conformers';
            $status->{$geometry}->{msg} = "new comformers made from this old one\n";

            chdir("$parent");
        }
    }
    return $conf_made;
} #End of make conforemers



sub submit_job {
    my ($com_file, $short) = @_;

    my ($wall, $nprocs);
    if ($short) {
        $wall = $system->{SHORT_WALL};
        $nprocs = $system->{SHORT_PROCS};
    }else {
        $wall = $system->{WALL};
        $nprocs = $system->{N_PROCS};
    }

    my $queue_name = $system->{QUEUE_NAME};

    unless($arg_in{nosub}) {
        $launch_failed += submit($com_file, $wall, $nprocs, $template_job, $queue_name);
    }
    if ($launch_failed > MAX_LAUNCH_FAILED) {

        my $time = &count_time(60);

        my $msg = "AARON has failed to submit jobs to queue more than ";
        $msg .= MAX_LAUNCH_FAILED;
        $msg .= " times. AARON believe you may run out of hours or ";
        $msg .= "something wrong with the system. ";
        $msg .= "AARON will just sleep for one hour and continue.\n";
        $msg .= "AARON will restart at $time\n";
        $msg .= "sleeping...";

        print_message($msg);
        $launch_failed = 0;

        sleep(3600);
    }
}


sub add_constraints {
    my ($coords, $subs_atoms) = @_;

    my @constraint_atoms;

    if (@$subs_atoms) {
        my %in = map {$_ => 1} @$subs_atoms;
        @constraint_atoms = grep {not $in{$_}} (0..$#{ $coords });
    }else {
        @constraint_atoms = (0..$template_info{first_cat_atom}-1);
    }

    foreach my $atom (@constraint_atoms) {
        $coords->[$atom]->[1] = "-1";
    }

    return $coords;
}


sub build_sub_lib {
    my ($geometry) = @_;
    my $head = &get_head($geometry);
    my $file_name = &get_file_name($geometry);

    unless ($head =~ /^mapping/) {
        return 0;
    }

    my ($stereo, $ts);
    if ($geometry =~ /([RS])\/(ts\d+)/) {
        $stereo = $1;
        $ts = $2;
    }

    unless (-d "$arg_in{TS_path_sub}/$stereo") {
        print "arg_in{TS_path_sub}/$stereo";
        make_path("$arg_in{TS_path_sub}/$stereo") or return 1;
    }

    unless (-e "$arg_in{TS_path_sub}/$stereo/$ts.xyz") {
        cp("$geometry/$file_name.xyz", "$arg_in{TS_path_sub}/$stereo/$ts.xyz") or return 1;
        my $msg .= "TS library for substitution were built in $arg_in{TS_path_sub}\n";
        $status->{$geometry}->{msg} .= $msg;
    }

    return 0;
}

sub count_time {
    my ($sleep_time) = @_;

    my $time = localtime;
    
    my $sleep_hour = int($sleep_time/60);
    my $sleep_minute = $sleep_time%60;

    if ($time =~ /\s(\d+)\:(\d+)/) {
        my $hour = $1 + $sleep_hour;
        my $minute = $2 + $sleep_minute;
        
        if ($minute > 60) {
            $hour += 1;
            $minute += $minute%60;
        }

        $time =~ s/\s\d+\:\d+/ $hour:$minute/;
    }

    return $time;
}


sub print_status {

    print "\033[2J";                                              #clear the screen
    print "\033[0;0H";                                    #jump to 0,0
    my $date1=localtime;                          #current date
    
    my @start;
    my @running;
    my @done;
    my @finished;
    my @restart;
    my @pending;
    my @killed;

    foreach my $geometry (sort keys %{ $status }) {
        STATUS: {
            if ($status->{$geometry}->{status} eq 'start') { push(@start, $geometry);
                                                             last STATUS; }

            if ($status->{$geometry}->{status} eq 'restart') { push(@restart, $geometry);
                                                               last STATUS; }

            if ($status->{$geometry}->{status} eq 'pending') { push(@pending, $geometry);
                                                               last STATUS;}

            if ($status->{$geometry}->{status} eq 'done') { push(@done, $geometry);
                                                            last STATUS;}

            if ($status->{$geometry}->{status} eq 'finished') { push(@finished, $geometry);
                                                                last STATUS;}

            if ($status->{$geometry}->{status} eq 'running') {push(@running, $geometry);
                                                               last STATUS;}

            if ($status->{$geometry}->{status} eq 'killed' ||
                $status->{$geometry}->{status} eq 'new_conformers') {push(@killed, $geometry);
                                                                    last STATUS;}
        }
    }

    print "Status of all jobs...($date1)\n";
    my $msg;

    @start && do {$msg .= "The following jobs are going to start:\n";};
    for my $geometry(@start) {
        $msg .= "job $geometry is starting the AARON workflow using geometry from TS libraries\n";
    }

    @done && do {$msg .= "The following jobs are done:\n";};
    for my $geometry(@done) {
        my $step_done = $status->{$geometry}->{step} - 1;
        $msg .= "job $geometry step $step_done done\n";
    }

    @finished && do {$msg .= "The following AARON are finished: \n";};
    for my $geometry(@finished) {
        $msg .= "$geometry finished normally\n";
    }

    @running && do {$msg .= "The following jobs are running:\n";};
    for my $geometry(@running) {
        my $filename = &get_file_name($geometry);
        my $step = $status->{$geometry}->{step};
        $msg .= "job $geometry step $step attempt $status->{$geometry}->{attempt} " .
                "cycle $status->{$geometry}->{cycle} running:\n";
    }

    @pending && do {$msg .= "The following jobs are pending\n";};
    for my $geometry(@pending) {
        $msg .= "job $geometry step $status->{$geometry}->{step} attempt $status->{$geometry}->{attempt}: ";
        $msg .= $status->{$geometry}->{msg} ? "$status->{geometry}->{msg}\n" : "No msg recorded\n";
    }

    @restart && do {$msg .= "The following jobs are restarted by some reasons\n";};
    for my $geometry(@restart) {
        $msg .= "job $geometry step $status->{$geometry}->{step} " .
                "restarted by reason: \n";
        if ($status->{$geometry}->{msg}) {
            $msg .= $status->{$geometry}->{msg};
        }else {
            $msg .= "No restart message recorded\n";
        }
        $msg .= "Now at attempt $status->{$geometry}->{attempt}.\n";
    }

    @killed && do {$msg .= "The following jobs are stopped by reason:\n";};
    for my $geometry(@killed) {
        $msg .= "job $geometry\n step $status->{$geometry}->{step} attemp $status->{$geometry}->{attempt}: ";
        $msg .= $status->{$geometry}->{msg} ?
                "$status->{$geometry}->{msg}\n" : "No msg recorded";

        if ($status->{$geometry}->{status} eq 'new_conformers') {
            delete $status->{$geometry};
        }
    }

    print "$msg\n\n";

    unless (@done || @running || @pending || @restart) {
        return 0;
    }else {
        return 1;
    }
}


sub analyze_result {

    my $thermo_num = $arg_in{high_method} ? 8 : 4;
    
    foreach my $head (sort keys %{ $thermo }) {
        my $ee = [];
        my $rel_thermos = {};
        foreach my $i (0..$thermo_num-1) {
            my $R_sum;
            my $S_sum;
            my @geo = grep { $thermo->{$head}->{$_}->[$i] } keys %{ $thermo->{$head} };
            my @sorted_geo = sort {$thermo->{$head}->{$a}->[$i] 
                                   <=> $thermo->{$head}->{$b}->[$i]} @geo;

            my $lowest_geo = $sorted_geo[0];

            foreach my $geo (@sorted_geo) {
                my $rel_thermo = ($thermo->{$head}->{$geo}->[$i] 
                                  - $thermo->{$head}->{$lowest_geo}->[$i])
                                  * HART_TO_KCAL; 

                $R_sum += exp(-$rel_thermo/$RT) if ($geo =~ /\/R\//);
                $S_sum += exp(-$rel_thermo/$RT) if ($geo =~ /\/S\//);

                $geo =~ s/$head\///;
                unless ($rel_thermos->{$geo}) {
                    $rel_thermos->{$geo} = [];
                }
                $rel_thermos->{$geo}->[$i] = $rel_thermo;
            }

            if ($R_sum && $S_sum) {
                $ee->[$i] = 100 * ($R_sum - $S_sum) / ($R_sum + $S_sum);
            }
        }

        #print result to STANDOUT and .log file
        print_ee($head, $ee, $rel_thermos);
    }
}


1;
