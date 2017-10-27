package AaronInit;

use strict;
use lib $ENV{'AARON'};

use Getopt::Long;
use Constants qw(:INFORMATION :THEORY :PHYSICAL :SYSTEM :JOB_FILE);
use Exporter qw(import);
use AaronTools;

#default values for some argument

my @authors = @{ INFO->{AUTHORS} };
my $version = INFO->{VERSION};
my $lastupdate = INFO->{LASTUPDATE};

our @EXPORT = qw(%arg_in %arg_parser %template_info $template_job $system init_main grab_cata_coords);

my $helpMsg = "\nAARON: An Automated Reaction Optimizer for Non-metal-catalyzed reactions  - Computational toolkit to aid in optimizing transition states with the Gaussian09 quantum chemistry software.\n\nAuthors: @authors.\n\nLast Update: $lastupdate\n\nAARON automates the optimization of transition states for a wide variety of metal-free asymmetric reactions. Based on a library of TS structures previously computed using model catalysts, AARON replaces the model catalyst with a user-supplied catalyst and then performs a prescribed series of constrained and uncontrained optimizations to arrive at final predicted structures and energies for all transition states.  

Command Line Options
     \"-debug\"    - Runs all optimizations using the PM6 semi-empirical method with a wall time of 2 hours.
     \"-pcm\"      - Changes the solvent model used from the standard PCM to a user specified model. Prompts user for Gaussian
     \"-solvent\"  - select solvent use in calculation
     \"-nosub\"    - Builds .com files for the current step without submitting anything
     \"-method\"   - allows user to set method 
     \"-basis\"    - allows user to set basis
     \"-gen\"      - allows user to set basis from gen
     \"-noquota\"  - ignore \"Erroneous write\" errors and re-run jobs
     NOTE: some of these don't work yet in version $version\n
";

#arguments for AARON taking from command line
our %arg_parser;

#content of template job file
our $template_job = {};

#all arguments for AARON
our %arg_in = (
    jobname => '',
    debug => 0,
    nosub => 0,
    noquota => 0,
    nosubstitute => 0,
    TS_path => '',
    MaxRTS => 1,
    MaxSTS => 1,
    solvent => 'gas',
    temperature => ROOM_TEMPERATURE,
    reaction_type => '',
    method => METHOD,
    high_method => HIGH_METHOD,
    lowmethod => LOW_METHOD,
    basis => BASIS,
    high_basis => HIGH_BASIS,
    gen => '',
    pcm => PCM,
    catalyst => '',
    charge => '',
    mult => '',
    template => '',
    key_atoms_cata => [],
    bonds_RMSD => [],
    bonds_LJ => [],
);

#THe default system is Ada if not specified
#$system is a hash reference
our $system = ADA;

#hash stroing template information
our %template_info;

#read command line arguments
sub read_args{
    GetOptions(
        'debug' => \$arg_parser{debug},
        'wall=i' => \$arg_parser{wall},
        'shortwall=i' => \$arg_parser{shortwall},
        'nprocs=i' => \$arg_parser{nprocs},
        'shortprocs' => \$arg_parser{shortporcs},
        'nosub' => \$arg_parser{nosub},
        'noquota' => \$arg_parser{noquota},
        'help' => \$arg_parser{help},
        'jobname=s' => \$arg_parser{jobname},
        'map_cata' => \$arg_parser{map_cata},
        'restart' => \$arg_parser{restart},
    ) or die ("Error in command line arguments\n $helpMsg");

    if ($arg_parser{help}) {
        print $helpMsg;
        exit(1);
    }

    #validae arguments
    unless ($arg_parser{jobname}) {
        print "You must indicate a jobname to start Aaron\n";
        exit (1);
    }
}


#read arguments from input file
sub read_params {
    my $input_file = $arg_parser{jobname} . '.in';
    open my $in_h, "< $input_file" or die "Can't open job_info!";

    while(<$in_h>) {
        /TS_path=(\S+)/ && do {$arg_in{TS_path} = $1; next};
        /MaxRTS=(\S+)/ && do {$arg_in{MaxRTS} = $1; next};
        /MaxSTS=(\S+)/ && do {$arg_in{MaxSTS} = $1; next;};
        /[sS]olvent=(\S+)/ && do {$arg_in{solvent} = $1; next;};
        /[tT]emperature=(\S+)/ && do {$arg_in{temperature} = $1; next;};
        /[rR]eaction_type=(\S+)/ && do {$arg_in{reaction_type} = $1; next;};
        /^[mM]ethod=(\S+)/ && do {$arg_in{method} = $1; next;};
        /[hH]igh_method=(\S+)/ && do {$arg_in{high_method} = $1; next;};
        /^[bB]asis=(\S+)/ && do {$arg_in{basis} = $1; next;};
        /^[hH]igh_basis=(\S+)/ && do {$arg_in{high_basis} = $1; next;};
        /[gG]en=(\S+)/ && do {$arg_in{gen} = $1; next;};
        /[lL]owmethod=(\S+)/ && do {$arg_in{lowmethod} = $1; next;};
        /^[wW]all=(\S+)/ && do {$arg_in{wall} = $1; next;};
        /[Nn]procs=(\S+)/ && do {$arg_in{nprocs} = $1; next;};
        /[sS]hortwall=(\S+)/ && do {$arg_in{shortwall} = $1; next;};
        /[sS]hortprocs=(\S+)/ && do {$arg_in{shortprocs} = $1; next;};
        /[cC]atalyst=(\S+)/ && do {$arg_in{catalyst} = $1; next;};
        /[cC]harge=(\S+)/ && do {$arg_in{charge} = $1; next;};
        /[mM]ult=(\S+)/ && do {$arg_in{mult} = $1; next;}; 
        /[tT]emplate=(\S+)/ && do {$arg_in{template} = $1; next;};
        /[kK]ey_atoms_cata=(.*)/ && do {
                                        my @key_atoms_cata;
                                        my @atoms_temp = split(/;/, $1);
                                        while(@atoms_temp) {
                                            my @atom_temp = split(/,/, shift(@atoms_temp));
                                            push (@key_atoms_cata, \@atom_temp);
                                        }
                                        $arg_in{key_atoms_cata} = [ @key_atoms_cata ];
                                        next;
                                       };
        /[bB]onds_RMSD=(.*)/ && do {
                                     my @bonds_RMSD;
                                     my @bonds_temp = split(/;/, $1);
                                     while(@bonds_temp) {
                                         my @bond_temp = split(/,/, shift(@bonds_temp));
                                         push (@bonds_RMSD, \@bond_temp);
                                     }
                                     $arg_in{bonds_RMSD} = [ @bonds_RMSD ];
                                     next;
                                    };
        /[bB]onds_LJ=(.*)/ && do {
                                   my @bonds_LJ;
                                   my @bonds_temp = split(/;/, $1);
                                   while(@bonds_temp) {
                                       my @bond_temp = split(/,/, shift(@bonds_temp));
                                       push (@bonds_LJ, \@bond_temp);
                                   }
                                   $arg_in{bonds_LJ} = [ @bonds_LJ ];
                                   next;
                                 };
        /^\&[sS]ystem/ && do {
                              while(<$in_h>) {
                                  /^\&$/ && do {last;};
                                  /n_procs=(\d+)/ && do {$system->{N_PROCS}=$1; next;};
                                  /wall=(\d+)/ && do {$system->{WALL}=$1; next;};
                                  /short_procs=(\d+)/ && do {$system->{SHORT_PROCS}=$1; next;};
                                  /short_wall=(\d+)/ && do {$system->{SHORT_WALL}=$1; next;};
                                  /node_type=(\S+)/ && do {$system->{NODE_TYPE}=$1; next;};
                              }
                             };
    }
    close $in_h;
    
    #examine arguments;
    if  (! $arg_in{TS_path}) {
        #no TS_path defined in <filename>.in file, built_in TS_lib is in use
        $arg_in{TS_path} = $AARON . "/TS_geoms/";
    }
    $arg_in{TS_path_reac} = $arg_in{TS_path}.$arg_in{reaction_type}.'/';

    unless($arg_in{MaxRTS}) {
        print "Can't find MaxRTS in $arg_parser{jobname}.in to run Aaron\n";
        exit(1);
    }

    unless($arg_in{reaction_type}) {
        print "Can't find reaction_type in $arg_in{reaction_type}.in to run Aaron\n";
        exit(1);
    }

    unless($arg_in{temperature}) {
        print "Can't find temperature in $arg_in{temperature}.in to run Aaron\n";
        exit(1);
    }

    unless ($arg_in{template}) {
        print "A template must be figure out explicitly in the <jobname>.in file ";
        print "by template=catalyst. \n";
        print "If this catalyst contains different steps in a reaction, ";
        print "you should figure out the step too. e.g. template=catalyst/TS1. \n";
        print "Exit without calculation\n";
        exit(1);
    }

    
    #make TS_path_map and TS_path_sub
    if ($arg_parser{map_cata}) {
        unless ($arg_in{catalyst}) {
            print "You want to map a new catalyst, but you didn't figure out the name of catalyst.\n";
            print "If the transition state you want to explore is one of the steps, ";
            print "please write catalyst/TSx in the <jobname>.in. e.g. catalyst=catalyst/TS1.\n";
            print "Exit without calculation\n";
            exit(1);
        }

        $arg_in{TS_path_map} = $arg_in{TS_path_reac}.$arg_in{template};
        $arg_in{TS_path_sub} = $arg_in{TS_path_reac}.$arg_in{catalyst};
        &read_reaction_data($arg_in{TS_path_map});
    }else{
        $arg_in{TS_path_sub} = $arg_in{TS_path_reac}.$arg_in{template};
        &read_reaction_data($arg_in{TS_path_sub});
    }

    #assign each element in arg_parser to arg_in if exist
    foreach my $arg_parser_key (keys %arg_parser) {
        if (defined $arg_parser{$arg_parser_key} && defined $arg_in{$arg_parser_key}) {
            $arg_in{$arg_parser_key} = $arg_parser{$arg_parser_key};
        }
    }
}

#reads template information from reaction_data file and store in our %template_info.
#returns key_atom_description
#my $key_atom_description = read_reaction_data();
sub read_reaction_data {
    my ($TS_path_temp) = @_;
    my $key_atoms_desc;
    open INFILE, "<$TS_path_temp/reaction_data" or die "Can't open reaction_data";
    while (<INFILE>) {
        chomp;
        /bond_constraint=(.*)/           && do {
                                                    my @constraints;
                                                    my @constraints_temp = split(/;/, $1);
                                                    while(@constraints_temp) {
                                                        my @constraint_temp = split(/,/, shift(@constraints_temp));
                                                        push (@constraints, \@constraint_temp);
                                                    }
                                                    $template_info{constraints} = [ @constraints ];
                                                    next;
                                                  }; 
          
         /subdirs=(.*)/                 && do { $template_info{subdirs} = [ split(/\s+/, $1) ]; 
                                                    next;
                                                  };
         /first_cat_atom=(\d+)/             && do { $template_info{first_cat_atom} = $1; 
                                                    next;
                                                  };
         /key_atoms_description=(.+)/    && do { $template_info{key_atoms_desc} = $1; next;};

         /key_atoms=(.*)/           && do {
                                              my @key_atoms;
                                              my @atoms_temp = split(/;/, $1);
                                              while(@atoms_temp) {
                                                 my @atom_temp = split(/,/, shift(@atoms_temp));
                                                 push (@key_atoms, \@atom_temp);
                                              }
                                              $template_info{key_atoms} = [ @key_atoms ];
                                              next;
                                             };
    }
    close(INFILE);
}

sub write_params {
    print "Initializing job...\n";
    my $infile = $arg_in{jobname}.'.in';
    open my $in_h, ">$infile" or die "Can't open job_info!";
    print $in_h "Input parameters needed to start Aaron...See Aaron manul for details\n\n";

    if ($arg_in{TS_path} ne $AARON) {
        print $in_h "TS_path=$arg_in{TS_path}\n";
    }

    $arg_in{catalyst} && print $in_h "catalyst=$arg_in{catalyst}\n";

    $#{ $arg_in{key_atoms_cata} } >= 0 && do{
        print $in_h "key_atoms_cata=";
        my $temp_key_atoms;
        for(@{ $arg_in{key_atoms_cata} }) {
           my $pair = join(',', @{$_});
           $temp_key_atoms .= "$pair;";
         }
         print $in_h "$temp_key_atoms\n";
    };

    $#{ $arg_in{bonds_RMSD} } >= 0 && do{
        print $in_h "bonds_RMSD=";
        my $temp_bonds;
        for(@{ $arg_in{bonds_RMSD} }) {
          my $pair = join(',', @{$_});
          $temp_bonds .= "$pair;";
        }
        print $in_h  "$temp_bonds\n";
    };

    $#{ $arg_in{bonds_LJ} } >= 0 && do{
        print $in_h "bonds_LJ=";
        my $temp_bonds;
        for(@{ $arg_in{bonds_LJ} }) {
          my $pair = join(',', @{$_});
          $temp_bonds .= "$pair;";
        }
        print $in_h "$temp_bonds\n";
    };

    print $in_h "solvent=$arg_in{solvent}\n";
    print $in_h "template=$arg_in{template}\n";
    print $in_h "temperature=$arg_in{temperature}\n";
    print $in_h "reaction_type=$arg_in{reaction_type}\n";
    print $in_h "method=$arg_in{method}\n";
    $arg_in{high_method} && print $in_h "high_method=$arg_in{high_method}\n";
    print $in_h "pcm=$arg_in{pcm}\n";
    print $in_h "basis=$arg_in{basis}\n";
    $arg_in{gen} && print $in_h "gen=$arg_in{gen}\n";
    $arg_in{high_basis} && print $in_h "high_basis=$arg_in{high_basis}\n";
    $arg_in{lowmethod} && print $in_h "lowmethod=$arg_in{lowmethod}\n";
    print $in_h "charge=$arg_in{charge}\n";
    print $in_h "mult=$arg_in{mult}\n";
    
    if ($system) {
        print $in_h "&System\n";
        print $in_h "n_procs=$system->{N_PROCS}\n";
        print $in_h "wall=$system->{WALL}\n";
        print $in_h "short_procs=$system->{SHORT_PROCS}\n" if $system->{SHORT_PROCS};
        print $in_h "short_wall=$system->{SHORT_WALL}\n" if $system->{SHORT_WALL};
        print $in_h "node_type=$system->{NODE_TYPE}\n" if $system->{NODE_TYPE};
        print $in_h "&"
    }
    close $in_h;
}


sub get_job_template {
    if ( -e "$AARON/template.job") {
        my $job_invalid;
        my $template_pattern = TEMPLATE_JOB;
        $template_job->{job} = "$AARON/template.job";
        $template_job->{formula} = {};
        open JOB, "<$AARON/template.job";
        #get formulas
        while (<JOB>) {
            /&formula&/ && do { while (<JOB>) {
                                    /&formula&/ && last;
                                    /^(\S+)=(\S+)$/ && do {  my $formula = $2;
                                                             my @pattern = grep {$formula =~ 
                                                                    /\Q$_\E/} values %$template_pattern;

                                                             unless (@pattern) {
                                                                print "template.job in $AARON is invalid. " .
                                                                      "Formula expression is wrong. " .
                                                                      "Please see manual.\n";
                                                                $job_invalid = 1;
                                                                last;
                                                            }
                                                            $template_job->{formula}->{$1} = $2 };
                                }
                                last if $job_invalid;
                              }
        }

    }
}


sub grab_cata_coords {
    #change catalyst from catalyst/TS into catalyst_TS
    my @temp = split('/', $arg_in{catalyst});
    my $catalyst = shift(@temp);
    $catalyst =~ s/\//_/g;
    my $cata_xyz = $catalyst.".xyz";
    unless (-e $cata_xyz) {
        print "You want to map new catalyst, but you didn't provide geometry of the new catalyst.\n";
        print "Exist without calculation.\n";
        exit(1);
    }
    my @cat_coords = grab_coords($cata_xyz);

    return @cat_coords;
}


#main function to initiate Aaron job
sub init_main {
    print "Preparing to run transition state searches...\n";
    sleep(2);
    &read_args();
    &read_params();
    &write_params();
    &get_job_template();
}
1;
