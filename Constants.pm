package Constants;

use strict;
use Cwd qw(cwd);
use base 'Exporter';

our @EXPORT = ();
our @EXPORT_OK = ('INFO', 'EOS', 'ADA', 'METHOD', 'HIGH_METHOD', 'BASIS', 'HIGH_BASIS', 'LOW_METHOD', 
                  'PCM', 'BOLTZMANN', 'AMU_TO_KG', 'HART_TO_KCAL', 'ROOM_TEMPERATURE', 'LAUNCH_DIR',
                  'CUTOFF', 'MAXSTEP', 'MAX_LAUNCH_FAILED', 'TEMPLATE_JOB' );

our %EXPORT_TAGS = (
                     INFORMATION => [ 'INFO' ],
                     SYSTEM => [ 'EOS', 'ADA' ],
                     THEORY => [ 'METHOD', 'HIGH_METHOD', 'BASIS',
                                 'HIGH_BASIS', 'LOW_METHOD', 'PCM'],
                     PHYSICAL => [ 'BOLTZMANN', 'AMU_TO_KG', 'HART_TO_KCAL', 'ROOM_TEMPERATURE'],
                     OTHER_USEFUL => [ 'LAUNCH_DIR', 'MAXSTEP', 'MAX_LAUNCH_FAILED' ],
                     COMPARE => [ 'CUTOFF' ],
                     JOB_FILE => ['TEMPLATE_JOB'],
                   );

#Code Information
use constant INFO => {
    VERSION => 1.0,
    YEAR => 2017,
    LASTUPDATE => '9/07/17',
    AUTHORS => ["Yanfei Guan", "Benjamin J. Rooks", "Steven E. Wheeler"],
    AARON_HOME => "/home/einsteinguan/bin/Aaron"  		#absolute path to Aaron directory
};

#System depedent Variables
use constant EOS => {
    N_PROCS => 8,
    WALL => 48,
    SHORT_PROCS => 4,
    SHORT_WALL => 2,
    QUEUE => 'PBS'
};

use constant ADA => {
    N_PROCS => 20,
    WALL => 24,
    SHORT_PROCS => 20,
    SHORT_WALL => 12,
    NO_GET_QUOTA => 1,
    QUEUE => 'LSF'
};

#Theoretical methods
use constant {
    METHOD => 'B97D',
    HIGH_METHOD => 'wB97XD',
    BASIS => 'TZV2d2p',
    HIGH_BASIS => 'TZV2d2p',
    LOW_METHOD => 'PM6',
    PCM => 'PCM',
};

#RMSD constants
use constant CUTOFF => {
    E_CUTOFF => 0.2,
    RMSD_CUTOFF => 0.5,
    D_CUTOFF => 0.35,
};

#template.job constant
use constant TEMPLATE_JOB => {
    JOB_NAME => '$jobname',
    WALL_TIME => '$walltime',
    N_PROCS => '$numprocs',
    NODE_TYPE => '$nodetype',
};

#Physical constants
use constant BOLTZMANN => 0.001987204;
use constant AMU_TO_KG => 1.66053886E-27;
use constant HART_TO_KCAL => 627.5095;
use constant ROOM_TEMPERATURE => 298.15;

#Other useful constants
use constant LAUNCH_DIR => cwd;
use constant MAXSTEP => 5;
use constant MAX_LAUNCH_FAILED => 5;
