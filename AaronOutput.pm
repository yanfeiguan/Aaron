package AaronOutput;

use strict;

use lib $ENV{'AARON'};

use Constants qw(:INFORMATION :THEORY :PHYSICAL);
use AaronInit qw(%arg_in %template_info $system);
use AaronTools;

use Cwd;
use Exporter qw(import);

our @EXPORT = qw(&init_log $parent header restart_header print_message print_params close_logfile clean_up print_ee);
our $parent = getcwd;

my $ol;
my $out_file;

sub init_log {
    $out_file = $parent . '/' . $arg_in{jobname} . "_aaron.log";
    open $ol, ">>$out_file" or die "Can't open $out_file\n";
}

#Prints Aaron header to $ol.  open STDOUT as $ol to write to screen
sub header {

    my $date = localtime;
    my $version = INFO->{VERSION};
    my @authors = @{ INFO->{AUTHORS} };
    my $year = INFO->{YEAR};

    &print_message("Aaron job started on $date\n\n");

    print $ol "                            Welcome to AARON!\n";
    print $ol "                                (v. $version)\n\n";
    print $ol "                               Written by\n";
    foreach my $author (@authors) {
      print $ol "                            $author\n";
    }
    print $ol "                          Texas A&M University\n";
    print $ol "                         September, 2013 - 2017\n\n";
    print $ol "                                          /\\          /\\ \n";
    print $ol "                                         ( \\\\        // )\n";
    print $ol "                                          \\ \\\\      // / \n";
    print $ol "                                           \\_\\\\||||//_/  \n";
    print $ol "                                            \\/ _  _ \\    \n";
    print $ol "                                           \\/|(O)(O)|    \n";
    print $ol "                                          \\/ |      |    \n";
    print $ol "                      ___________________\\/  \\      /    \n";
    print $ol "                     //  |          |  //     |____|     \n";
    print $ol "                    //   |A.A.R.O.N.| ||     /      \\    \n";
    print $ol "                   //|   |          | \\|     \\ 0  0 /    \n";
    print $ol "                  // \\   ||||||||||||  V    / \\____/     \n";
    print $ol "                 //   \\     /________(   _ /             \n";
    print $ol "                \"\"     \\   /    /    |  ||               \n";
    print $ol "                       /  /\\   /     |  ||               \n";
    print $ol "                      /  / /  /      \\  ||               \n";
    print $ol "                      | |  | |        | ||               \n";
    print $ol "                      |_|  |_|        |_||               \n";
    print $ol "                       \\_\\  \\_\\        \\_\\\\              \n\n\n";
    print $ol "                              Automated\n";
    print $ol "                              Alkylation\n";
    print $ol "                              Reaction\n";
    print $ol "                              Optimizer for\n";
    print $ol "                              New catalyst\n\n";
    print $ol "Citation:\n";
    print $ol "AARON, verson $version, S. E. Wheeler and B. J. Rooks, Texas A&M University, $year.\n\n";
    print $ol "B. J. Rooks, M. R. Haas, D. Sepulveda, T. Lu, and S. E. Wheeler, \"Prospects for the
Co  mputational Design of Bipyridine N,N\'-Dioxide Catalysts for Asymmetric Propargylations\" ACS Catalysis 
5,   272 (2015).\n\n";
    print $ol "The following should also be cited when AARON is used for bidentate Lewis-base catalyzed alkylation reactions:\n\n";
    print $ol "1. T. Lu, M. A. Porterfield, and S. E. Wheeler, \"Explaining the Disparate Stereoselectivities
of   N-Oxide Catalyzed Allylations and Propargylations of Aromatic Aldehydes\", Org. Lett. 14, 
53  10 (2012).\n\n";
    print $ol "2. T. Lu, R. Zhu, Y. An, and S. E. Wheeler, \"Origin of Enantioselectivity in the 
Pr  opargylation of Aromatic Aldehydes Catalyzed by Helical N-Oxides\", J. Am. Chem. Soc. 134, 
30  95 (2012).\n\n";
    print $ol "3. D. Sepulveda, T. Lu, and S. E. Wheeler, \"Performance of DFT Methods and Origin of
St  ereoselectivity in Bipyridine N,N\'-Dioxide Catalyzed Allylation and Propargylation Reactions\", 
12  , 8346 (2014).\n\n";
    print $ol "The development of AARON is sponsored in part by the National Science Foundation,\nGrant CHE-1266022.\n\n\n";
} #end sub header


sub restart_header {
    my $date=localtime;
    if (-e $out_file) {
        print $ol "\n---------------------------------------------------------\nAaron job restarted on $date\n\n";
    } else {
        print $ol "Aaron job restarted on $date\n\n";
        &header();
        &print_params();
    }
}


sub print_message {
    print "$_[0]";
    print $ol "$_[0]";
}


#print all job parameters to $ol
sub print_params {
    my $version = INFO->{VERSION};
    my $AARON_HOME = INFO->{AARON_HOME};

    print $ol "----------------------------------------------------------------------------------\n";
    print $ol "Parameters\n";
    print $ol "----------------------------------------------------------------------------------\n";
    print $ol " AARON_HOME          = $AARON_HOME\n";
    print $ol "  version            = $version\n";
    print $ol "\n Reaction parameters:\n";
    print $ol "  reaction_type      = $arg_in{reaction_type}\n";
    print $ol "  solvent            = $arg_in{solvent}\n";
    print $ol "  temperature        = $arg_in{temperature} K\n";
    print $ol "  MaxRTS             = $arg_in{MaxRTS}\n";
    print $ol "  MaxSTS             = $arg_in{MaxSTS}\n";
    print $ol "  TS_path            = $arg_in{TS_path}\n";
    $arg_in{TS_path_matp} && print $ol "  TS_path_map        = $arg_in{TS_path_map}\n";
    print $ol "  TS_path_sub        = $arg_in{TS_path_sub}\n";
    print $ol "\n Methods:\n";
    print $ol "  method = $arg_in{method}\n";
    print $ol "  high level method  = $arg_in{high_method}\n";
    if($arg_in{basis}) {
        print $ol "  basis set file     = $arg_in{basis}\n";
    }
    print $ol "  solvent model      = $arg_in{pcm}\n";
    print $ol "  low-level method   = $arg_in{lowmethod}\n";
    print $ol "\n Queue parameters:\n";
    print $ol "  wall               = $system->{WALL} hours\n";
    print $ol "  nprocs             = $system->{N_PROCS}\n";
    print $ol "  shortwall          = $system->{SHORT_WALL} hours\n" if $system->{SHORT_WALL};
    print $ol "  shortprocs         = $system->{SHORT_PROCS}\n" if $system->{SHORT_PROCS};
    print $ol "  queue_name         = $system->{QUEUE_NAME}\n" if $system->{QUEUE_NAME};

    if(@ARGV) {
        print $ol "\n command-line flags  = @ARGV\n";
    }
    print $ol "----------------------------------------------------------------------------------\n\n";
} #end sub print_params


sub print_ee {
    my ($head, $ee, $rel_thermos) = @_;

    my $data = "\n\n====================================================================================\n";                    
    $data = "$head:\n";

    my @geo_R = grep { $_ =~ /^R\// } sort keys %{ $rel_thermos };
    my @geo_S = grep { $_ =~ /^S\// } sort keys %{ $rel_thermos };

    my $print_thermo = sub {
        my ($geo, $data) = @_;
        my $first = 1;
        foreach my $thermo (@{ $rel_thermos->{$geo} }) {
            if (defined $thermo) {
                if ($first) {
                    $data .= sprintf "%10.1f", $thermo;
                    $first = 0;
                }else {
                    $data .= sprintf "%20.1f", $thermo;
                }
            } else {
                if ($first) {
                    $data .= sprintf "%10s", '';
                    $first = 0;
                }else {
                    $data .= sprintf "%20s", '';
                }
            }
        }
        $data .= "\n";
        return $data;
    };

    if (@geo_R && @geo_S) {
        $data .= sprintf "%22s%20s%20s%20s%20s%20s%20s%20s\n", 'E', 'H', 'G', 'G_Grimme',
                                                'E_High', 'H_High', 'G_High', 'G_Grimme_High';
        $data .= sprintf "%2s", 'ee';
        foreach my $e (@$ee) {
            if ($e) {
                $data .= sprintf "%19.1f%%", $e;
            }else {
                $data .= sprintf "%20s", '';
            }
        }
        $data .= "\n";

        $data .= "R relative energy\n";
        foreach my $geo (@geo_R) {
            $data .= sprintf "%-12s", $geo;

            $data = &$print_thermo($geo, $data);
        }

        $data .= "S relative energy\n";
        foreach my $geo (@geo_S) {
            $data .= sprintf "%-12s", $geo;

            $data = &$print_thermo($geo, $data);
        }

        &print_message($data);
    }
}

sub clean_up {
    my ($workingfile) = @_;
    my $startdir = cwd;
    chdir ($workingfile) or die "Unable to enter dir $workingfile:$!\n";
    opendir(DIR, ".") or die "Unable to open dir $workingfile:$!\n";
    my @names = readdir(DIR) or die "Unable to read $workingfile:$!\n";
    closedir(DIR);
    for (@names) {
        /^\.+$/ && do {next;};

        -d && do {
            &clean_up($_);
            next;
        };

        /\.xyz/ && do {
            for (@names) {
                /\.job/ && do { system("rm -f $_"); next;};
            }
            next;
        }
    }
    chdir ($startdir);
} #End clean_up


sub close_logfile {
    my $date = localtime;
    print $ol "AARON stopped $date\n";
    close ($ol);
}

1;
