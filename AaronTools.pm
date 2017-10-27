#
#Aaron::AaronTools     Tools for working with molecular structures and G09 input and output files
#
#
# Copyright (c) 2015 Steven E. Wheeler. All rights reserved. This program
# is free software; you can redistribute it and/or modify it under the
# same terms as Perl itself.
#

package AaronTools;
#Perl module containing commonly used subroutines for building molecules, manipulating coordinates, parsing G09 log files, interacting with LSF and PBS queues, etc.

use strict;
use warnings;
use lib qw(/scratch/group/wheeler/perl/lib/perl5);

use Constants qw(:INFORMATION :SYSTEM :PHYSICAL :JOB_FILE);

use Cwd;
use Switch;
use Math::Trig;
use File::Basename;
use List::Util qw(min max);
use Math::Vector::Real;
use Math::MatrixReal;

use Exporter qw(import);
our @EXPORT = qw($AARON $queue_type $LSF %radii %vdw_radii %masses finished grab_coords mirror_coords same_structure copy_coords combine_coords printXYZ printPDB write_com num_imaginary grab_freqs print_freqs find_str findJob submit make_job_file change_distance distance angle dihedral change_dihedral set_dihedral get_quota coord_shift remove_atoms remove_fragment get_error get_connected get_connectivity check_connectivity rotate genrotate quat_rot center_genrotate center_ring substitute fused_ring shortest_path get_all_connected get_gradient get_energy get_homo_lumo get_thermo RMSD_align minimize_torsion LJ_energy calculate_ee match_molecules sub_geo map_catalyst fix_coords get_molecules sub_rotate RMSD random_string flatten read_IRC nt_builder);

#get system information from environement vironment
our $queue_type = $ENV{'QUEUE_TYPE'} ? 
                 $ENV{'QUEUE_TYPE'} : ADA->{QUEUE};						#flag for queue type
our $g09root = $ENV{'G09_ROOT'} ?
              $ENV{'G09_ROOT'} : "/software/lms/g09_D01";       #absolute path to root directory for Gaussian09
our $AARON = $ENV{'AARON'} ?
             $ENV{'AARON'} : "/home/einsteinguan/bin/AARON";

#Physical consants
my $boltzmann = 0.001987204;			#Boltzmann constant in kcal/mol
my $h = 6.62606957E-34; 			#J.s
my $kb = 1.380662E-23; 				#J/K
my $c = 29979245800; 				#cm/s
my $R = 1.987204; 				#cal/mol
my $kcal2hartree = 0.0015936;
my $amu2kg = 1.66053886E-27;
my $hart2kcal = 627.5095;
my $ang2bohr = 1.889725989;

#ATOM data
#covalent radii (from jmol source code)
our %radii = ( H=>0.32, He=>0.93, Li=>1.23, Be=>0.90, B=>0.82, C=>0.77, N=>0.75, O=>0.73, F=>0.72, Ne=>0.71, Na=>1.54, Mg=>1.36, Al=>1.18, Si=>1.11, P=>1.06, S=>1.02, Cl=>0.99, Ar=>0.98, K=>2.03, Ca=>1.74, Sc=>1.44, Ti=>1.32, V=>1.22, Cr=>1.18, Mn=>1.17, Fe=>1.17, Co=>1.16, Ni=>1.15, Cu=>1.17, Zn=>1.25, Ga=>1.26, Ge=>1.22, As=>1.20, Se=>1.16, Br=>1.14, Kr=>1.12, Rb=>2.16, Sr=>1.91, Y=>1.62, Zr=>1.45, Nb=>1.34, Mo=>1.30, Tc=>1.27, Ru=>1.25, Rh=>1.25, Pd=>1.28, Ag=>1.34, Cd=>1.48, In =>1.44, Sn=>1.41, Sb=>1.40, Te=>1.36, I=>1.33, Xe=>1.31, Cs=>2.35, Ba=>1.98, La=>1.69, Lu=>1.60, Hf=>1.44, Ta=>1.34, W=>1.30, Re=>1.28, Os=>1.26, Ir=>1.27, Pt=>1.30, Au=>1.34, Hg=>1.49, Tl=>1.48, Pb=>1.47, Bi=>1.46, X=>0);

our %vdw_radii = (H => 1.20, He => 1.40, Li => 1.82, Be => 1.3725, B => 0.795, C => 1.70, N => 1.55, O => 1.52, F => 1.47, Ne => 1.54, Na => 2.27, Mg => 1.73, Al => 1.7, Si => 2.10, P => 1.80, S => 1.80, Cl => 1.75, Ar => 1.88, K => 2.75, Ca => 2.45, Sc => 1.37, Ti => 1.37, V => 1.37, Cr => 1.37, Mn => 1.37, Fe => 1.456, Co => 0.88, Ni => 0.69, Cu => 0.72, Zn => 0.74, Ga => 1.37, Ge => 1.95, As => 1.85, Se => 1.90, Br => 1.85, Kr => 2.02, Rb => 1.58, Sr => 2.151, Y => 1.801, Zr => 1.602, Nb => 1.468, Mo => 1.526, Tc => 1.360, Ru => 1.339, Rh => 1.345, Pd => 1.376, Ag => 1.27, Cd => 1.424, In => 1.663, Sn => 2.10, Sb => 2.05, Te => 2.06, I => 1.98, Xe=>2.00, Cs=>1.84, Ba=>2.243, La=>1.877, Lu=>2.17, Hf=>1.580, Ta=>1.467, W=>1.534, Re=>1.375, Os=>1.353, Ir=>1.357, Pt=>1.75, Au=>1.66, Hg=>1.55, Tl=>1.96, Pb=>2.02, Bi=>2.15, X=>0);

our %masses = (X => '0.', H => '1.00782503207', He => '4.00260325415', Li => '7.016004548', Be => '9.012182201', B => '11.009305406', C => '12.0', N => '14.00307400478', O => '15.99491461956', F => '18.998403224', Ne => '19.99244017542', Na => '22.98976928087', Mg => '23.985041699', Al => '26.981538627', Si => '27.97692653246', P => '30.973761629', S => '31.972070999', Cl => '34.968852682', Ar => '39.96238312251', K => '38.963706679', Ca => '39.962590983', Sc => '44.955911909', Ti => '47.947946281', V => '50.943959507', Cr => '51.940507472', Mn => '54.938045141', Fe => '55.934937475', Co => '58.933195048', Ni => '57.935342907', Cu => '62.929597474', Zn => '63.929142222', Ga => '68.925573587', Ge => '73.921177767', As => '74.921596478', Se => '79.916521271', Br => '78.918337087', Kr => '85.910610729', Rb => '84.911789737', Sr => '87.905612124', Y => '88.905848295', Zr => '89.904704416', Nb => '92.906378058', Mo => '97.905408169', Tc => '98.906254747', Ru => '101.904349312', Rh => '102.905504292', Pd => '105.903485715', Ag => '106.90509682', Cd => '113.90335854', In => '114.903878484', Sn => '119.902194676', Sb => '120.903815686', Te => '129.906224399', I => '126.904472681', Xe => '131.904153457', Cs => '132.905451932', Ba => '137.905247237', La => '138.906353267', Lu => '174.940771819', Hf => '179.946549953', Ta => '180.947995763', W => '183.950931188', Re => '186.955753109', Os => '191.96148069', Ir => '192.96292643', Pt => '194.964791134', Au => '196.966568662', Hg => '201.970643011', Tl => '204.974427541', Pb => '207.976652071', Bi => '208.980398734');

#Lennard-Jones C12 and C6 parameters from Autodock (http://www.csb.yale.edu/userguides/datamanip/autodock/html/Using_AutoDock_305.a.html)
my %LJC12 = ( "CC" => '2516582.400', "CN" => '1198066.249', "CO" => '820711.722', "CS" => '2905899.052', "CH" => '29108.222', "NC" => '1198066.249', "NN" => '540675.281', "NO" => '357365.541', "NS" => '1383407.742', "NH" => '10581.989', "OC" => '820711.722', "ON" => '357365.541', "OO" => '230584.301', "OS" => '947676.268', "OH" => '6035.457', "SC" => '2905899.052', "SN" => '1383407.742', "SO" => '947676.268', "SS" => '3355443.200', "SH" => '33611.280', "HC" => '29108.222', "HN" => '10581.989', "HO" => '6035.457', "HS" => '33611.280', "HH" => '81.920', "X" => '0');

my %LJC6 = ("CC" => '1228.800000', "CN" => '861.634784', "CO" => '754.059521', "CS" => '1418.896022', "CH" => '79.857949', "NC" => '861.634784', "NN" => '588.245000', "NO" => '505.677729', "NS" => '994.930149', "NH" => '48.932922', "OC" => '754.059521', "ON" => '505.677729', "OO" => '429.496730', "OS" => '870.712934', "OH" => '39.075098', "SC" => '1418.896022', "SN" => '994.930149', "SO" => '870.712934', "SS" => '1638.400000', "SH" => '92.212017', "HC" => '79.857949', "HN" => '48.932922', "HO" => '39.075098', "HS" => '92.212017', "HH" => '2.560000', "X" => '0');

#Substituent geometry, if linear 0, nonlienar 1
my $sub_geo = { 
                Me => { conf => 0 },
                OMe => { conf => 1,
                         deg => pi,
                         num => 2,
                       },
                Ph => { conf => 1,
                        deg => pi/2,
                        num => 2,
                      },
                H => { conf => 0 },
              };

#checks if job has finished.  Returns 1 if finished, 0 else.  Accepts name of log file
sub finished {
  my $file = $_[0];
  my $norm="Normal termination";			#message for finished optimization
  open INFILE, "<$file" or die "Can't open $file!";
  while (<INFILE>) {
    if($_ =~ /$norm/) {
      close(INFILE);
      return 1;
    }
  }
  close(INFILE);
  return 0;
}  #end sub finished

# For .log files, extracts geometry of the last step in an optimization from a Gaussian LOG file.
# For .pdb, grabs coords, element name, residue name
# For .xyz, simply grabs coords
# In both cases, returns array of $coords: $coords[$atom][0] = atomic symbol, $coords[$atom][1] = 0, $coords[$atom][2-4] = x, y, z
# For PDB files, additional information in $coords: $coords[$atom][5] = ATOM|HETATM, $coords[$atom][6] = atom type, $coords[$atom][7] = residue name, $coods[$atom][8] = chain
# If file does not exist or can't find coords, returns 0
sub grab_coords {
  my $filename = $_[0];
  my @coords;
  my @bond_fix;
  if(-e $filename) {                                   #check to make sure file exists first
    open (INFILE, "<$filename") or die "Can't open $filename";
    my @elements=('Bq','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','X');

    if($filename =~ /\S+\.log/) {                       #G09 log file

    # Snatches last geometry from LOG file
      while (<INFILE>) {
        my $line=$_;
        if($line =~ / orientation:/) {
          @coords = ();
          $line = <INFILE>;
          $line = <INFILE>;
          $line = <INFILE>;
          $line = <INFILE>;
          $line = <INFILE>;
          do {
            if($line =~ /^\s+\d+\s+(\S+)\s+\S*\s+(\S+)\s+(\S+)\s+(\S+)/) {
              my @atom = ($elements[$1], 0, $2, $3, $4);
              push(@coords, \@atom);
              $line = <INFILE>;
            }
          } while(!($line =~ /--/));
        }
      }
    } elsif ($filename =~ /\S+\.xyz/) {         #XYZ file
      my $numatoms = <INFILE>;
      my $comment = <INFILE>;
      if($comment =~ /F:(\S+)/) {
        my @temp = split(/;/, $1);
        while (@temp) {
          my $bond = shift @temp;
          my @bond = split(/-/, $bond);
          push @bond_fix, \@bond;
        }
      }
      while(<INFILE>) {
        chomp;
        $_ =~ s/^\s+//;
        if($_ =~ /(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
          my @atom = ($1, 0, $2, $3, $4);
          push(@coords, \@atom);
        }
      }
    } elsif ($filename =~ /\S+\.pdb/) {		#PDB file
      while(<INFILE>) {
        #strip leading whitespace
        $_ =~ s/^\s+//;
	#Typical PDB file line: 
	#ATOM     13  CB  ASP A  23      -0.219   5.194 -16.219  1.00 30.89           C
        #****         **  ***            ******   *****  ******                       *
        #This is ugly, but allows us to skip over 'missing' entries in formatted PDB file!
	#              (ATOM)     13  (CB)  (ASP) A  (23)    (  -0.219)(   5.194)( -16.219)  1.00 30.89          ( C)
        if($_ =~ /(ATOM|HETATM).......(..)..(...)....(..)....(........)(........)(........)......................(..)/) {
          chomp;
          #strip leading whitespace off of everything
          my $type = $1;
          my $atom_type = $2;
	  my $res_name = $3;
	  my $chain = $4;
	  my $x = $5;
	  my $y = $6; 
          my $z = $7;
          my $element = $8;
          #strip all whitespace from element name, otherwise I can't use this as a hash key later!
          $element =~ s/\s+//;
          my @atom = ($element, 0, $x, $y, $z, $type, $atom_type, $res_name, $chain);
          push(@coords, \@atom);
        }
      }
    } elsif ($filename =~ /\S+\.com/) {	#G09 com file (Cartesian input!)
      while (<INFILE>) {
        chomp;
        $_ =~ s/^\s+//;

        if($_ =~ /^\s?([a-zA-Z]+)\s+(-?\d)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s*$/) {
          my @atom = ($1, $2, $3, $4, $5);
          push(@coords, \@atom);
        } elsif($_ =~ /^\s?([a-zA-Z]+)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s*$/) {
          my @atom = ($1, 0, $2, $3, $4);
          push(@coords, \@atom);
        }
      }
    }
    close(INFILE);
  }
  if (@bond_fix) {
    return (\@coords, \@bond_fix);
  } else {
    return @coords;
  }                          #This will be empty if either can't find coords or file does not exist
} #End sub grab_coords

#makes a copy of a @coords array, mirroring x-values.  Returns new @coords array
#mirror_coords($ref_to_coords);
sub mirror_coords {
  my ($coords_ref) = @_;

  my @new_coords = copy_coords($coords_ref);
  foreach my $atom (0..$#new_coords) {
    $new_coords[$atom][2] *= -1;
  }
   
  return @new_coords;
} #End sub mirror_coords

#compares RMSD of two structures (and their mirror images); returns 0 or 1 if structures are unique or different 
#runs RMSD_align several times, attempting to renumber atoms after each round!
sub same_structure {
  my ($coords1_ref, $coords2_ref, $cutoff) = @_;
  my @coords1 = @{$coords1_ref};
  my @coords2 = @{$coords2_ref};
  my $rmsd = 1e99;
  ($rmsd,@coords2) = RMSD_align(\@coords1, \@coords2);
  if($rmsd < $cutoff) {
    return 1;
  }
  foreach (0..2) {
    #attempt to renumber...
    my @temp_coords2;
    foreach my $atom (0..$#coords1) {	#for each atom in coords1, find closest atom of same element in coords2 and put into temp_compare_coords
      my $short_distance = 1e99;
      my $short_atom = 0;
      foreach my $atom2 (0..$#coords1) {
        if($coords2_ref->[$atom2][0] eq $coords1[$atom][0]) {
          my $distance = distance($atom2, $atom, $coords2_ref, \@coords1);
          if($distance < $short_distance) {
            $short_distance = $distance;
            $short_atom = $atom2;
          }
        }
      }
  #push closest atom onto @temp_compare_coords and remove from @compare_coords
      push(@temp_coords2, $coords2[$short_atom]);
    }
    @coords2 = @temp_coords2;
    ($rmsd,@coords2) = RMSD_align(\@coords1, \@coords2);
    if($rmsd < $cutoff) {
      return 1;
    }
  }
  return 0;
} #end same_structure


#makes a copy of a @coords array.  Returns new @coords array
#copy_coords($ref_to_coords);
sub copy_coords {
  my ($coords_ref) = @_;

  my @new_coords = map { [@$_] } @{$coords_ref};
   
  return @new_coords;
} #End sub copy_coords


#combines two sets of coords and returns new @coords array
#combine_coords($ref_to_coords1, $ref_to_coords2);
sub combine_coords {
  my ($coords1_ref, $coords2_ref) = @_;

  my @combined = (@{$coords1_ref}, @{$coords2_ref});
  my @new_coords = copy_coords(\@combined);

  return @new_coords;
} #End combine coords


#Prints coords in XYZ format from @coords
#printxyz(ref_to_coords, commentsub printXYZ {
sub printXYZ {
  my ($coords_ref, $comment, $filename) = @_;
  my @coords = @{$coords_ref};

  my $num_atoms = $#coords + 1;
  if($filename) {
    #print "Wrinting geometry to $filename\n";
    open OUTFILE, ">>$filename" or die "Can't open $filename";
    print OUTFILE "$num_atoms\n$comment\n";
    foreach my $atom (0..$#coords) {
      print OUTFILE "$coords[$atom][0] ";
      foreach my $j (2..4) {
        printf OUTFILE "%14.6f  ", $coords[$atom][$j];
      }
      print OUTFILE "\n";
    }
    print OUTFILE "\n\n";
    close(OUTFILE);
  } else {
    print "$num_atoms\n$comment\n";
    foreach my $atom (0..$#coords) {
      print "$coords[$atom][0] ";
      foreach my $j (2..4) {
        printf "%14.6f  ", $coords[$atom][$j];
      }
      print "\n";
    }
  }
} #End printXYZ


#Writes com file
#write_com(route, comment, charge, multiplicity, ref_to_coords, footer, flag)
#write_com(route, comment, charge, multiplicity, ref_to_coords, footer, flag, filename)
#where footer contains anything that goes after the coords (gen basis specification, mod redundant commands, etc)
#flag = 0 will print only elements and coords
#flag = 1 will print 0/-1 column as well
sub write_com {
  my ($route, $comment, $charge, $mult, $coords_ref, $footer, $flag, $filename) = @_;
  my @coords = @{$coords_ref};
  
  my $num_atoms = $#coords + 1;
  if($filename) {	#print to file
    open OUTFILE, ">$filename" or die "Can't open $filename";
    print OUTFILE "$route\n\n";
    print OUTFILE "$comment\n\n";
    print OUTFILE "$charge $mult\n";
    foreach my $atom (0..$#coords) {
      print OUTFILE "$coords[$atom][0] ";
      if($flag) {
        print OUTFILE " $coords[$atom][1] ";
      }
      foreach my $j (2..4) {
        printf OUTFILE "%14.6f  ", $coords[$atom][$j];
      }
      print OUTFILE "\n";
    }
    print OUTFILE "\n";
    if($footer) {
      print OUTFILE "$footer\n";
    }
    print OUTFILE "\n\n";
    close(OUTFILE);
  } else {		#print to STDOUT
    print "$route\n\n";
    print "$comment\n\n";
    print "$charge $mult\n";
    foreach my $atom (0..$#coords) {
      print "$coords[$atom][0] ";
      if($flag) {
        print " $coords[$atom][1] ";
      }
      foreach my $j (2..4) {
        printf "%14.6f  ", $coords[$atom][$j];
      }
      print "\n";
    }
    print "\n";
    if($footer) {
      print "$footer\n";
    }
    print "\n\n";
  }
} #End write_com

#Prints coords in PDB format from @coords.  If some information is missing (atom type, residue, etc), just prins placeholders
#printPDB(ref_to_coords) or
#printPDB(ref_to_coords, filename)
sub printPDB {
  my ($coords_ref, $filename) = @_;
  my @coords = @{$coords_ref};

  if($filename) {
    print "Writing PDB file to $filename\n";
    open OUTFILE, ">$filename" or die "Can't open $filename";
    foreach my $atom (0..$#coords) {
      printf OUTFILE "%-6s%5d  %-3s %3s  %4s    %8.3f%8.3f%8.3f                       %s\n", $coords[$atom][5], $atom + 1, $coords[$atom][6], $coords[$atom][7], $coords[$atom][8], $coords[$atom][2], $coords[$atom][3], $coords[$atom][4], $coords[$atom][0];
    }
    print OUTFILE "END\n";
  } else {
    foreach my $atom (0..$#coords) {
      printf "%-6s%5d  %-3s %3s          %8.3f%8.3f%8.3f                       %s\n", $coords[$atom][5], $atom + 1, $coords[$atom][6], $coords[$atom][7], $coords[$atom][2], $coords[$atom][3], $coords[$atom][4], $coords[$atom][0];
    }
    print "END\n";
  }
}


#Reads G09 log file and returns number of imaginary frequencies
sub num_imaginary {
  my ($logfile) = @_;

  open INFILE, "<$logfile" or die "Can't open $logfile";
  while (<INFILE>) {
    if($_ =~ /(\d+) imaginary frequencies \(negative Signs\)/) {
      close(INFILE);
      return $1;
    }
  }
  close(INFILE);
  return 0;
} #end num_imaginary


#Frequency data stored as an array of arrays of arrays:
#my @freqs (frequency values)
#my @vectors (array of arrays of arrays)
#Later, replace with frequency object!
#grab_freqs($filename);
#returns three references: frequencies, intensities, vectors (normal mode vectors)
sub grab_freqs {
  my ($file) = @_;

  open INFILE, "<$file" or die "Can't open $file";

  my @freqs;
  my @intensities;
  my @vectors = ();
  my $numatoms;

  my $freq_num = 0;
READFILE:
  while (<INFILE>) {
    my $line = $_;
    chomp $line;
    if($line =~ /NAtoms=\s+(\d+)/) {
      $numatoms = $1;
    }
    if($line =~ / Harmonic frequencies /) {
      foreach (0..4) {	#skip down to actual frequency data
        <INFILE>;
      }
      my $next_line = "";
      while($next_line !~ / Thermochemistry /) {
        $next_line = <INFILE>;
        chomp($next_line);
        $next_line =~ s/^\s+//;
        if($next_line =~ /Frequencies/) {
          my @array = split(/\s+/, $next_line);
          shift(@array); #discard "Frequencies"
          shift(@array); #discard "--"
          my $num_freq_this_row = $#array + 1;
          push(@freqs, @array);
          # Skip over Red. masses and Frc consts
          <INFILE>;
          <INFILE>;
          $next_line = <INFILE>;
          $next_line =~ s/^\s+//;
          @array = split(/\s+/, $next_line);
          shift(@array); #discared "IR";
          shift(@array); #discared "Inten";
          shift(@array); #discared "--";
          push(@intensities, @array);
          #Skip coordinate labels
          my $next_temp = <INFILE>;
          next READFILE unless $next_temp =~ /x\s+y\s+z/i;
          foreach my $atom (0..$numatoms) {
            $next_line = <INFILE>;
            if($next_line =~ /\s+(\d+)\s+\d+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
              #three frequencies 
              $vectors[$freq_num][$1-1][0] = $2;
              $vectors[$freq_num][$1-1][1] = $3;
              $vectors[$freq_num][$1-1][2] = $4;
              $vectors[$freq_num+1][$1-1][0] = $5;
              $vectors[$freq_num+1][$1-1][1] = $6;
              $vectors[$freq_num+1][$1-1][2] = $7;
              $vectors[$freq_num+2][$1-1][0] = $8;
              $vectors[$freq_num+2][$1-1][1] = $9;
              $vectors[$freq_num+2][$1-1][2] = $10;
            } elsif($next_line =~ /\s+(\d+)\s+\d+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
              #two frequencies
              $vectors[$freq_num][$1-1][0] = $2;
              $vectors[$freq_num][$1-1][1] = $3;
              $vectors[$freq_num][$1-1][2] = $4;
              $vectors[$freq_num+1][$1-1][0] = $5;
              $vectors[$freq_num+1][$1-1][1] = $6;
              $vectors[$freq_num+1][$1-1][2] = $7;
            } elsif($next_line =~ /\s+(\d+)\s+\d+\s+(\S+)\s+(\S+)\s+(\S+)/) {
              #one frequency
              $vectors[$freq_num][$1-1][0] = $2;
              $vectors[$freq_num][$1-1][1] = $3;
              $vectors[$freq_num][$1-1][2] = $4;
            }
          }
          $freq_num += $num_freq_this_row;
        }
      }
    }
  }
  return (\@freqs, \@intensities, \@vectors);
} #End sub grab_freqs

sub print_freqs {
  my ($freq_ref, $intens_ref, $vec_ref, $numatoms) = @_;
  my @frequencies = @{$freq_ref};
  my @intensities = @{$intens_ref};
  my @vectors2 = @{$vec_ref};

  foreach my $frequency (0..$#frequencies) {
    print "Freq = $frequencies[$frequency] ($intensities[$frequency])\n";
    print "X\tY\tZ\n";
    for my $atom (0..$numatoms) {
      printf "%3d %10.2f %10.2f %10.2f\n", $atom, $vectors2[$frequency][$atom][0],$vectors2[$frequency][$atom][1],$vectors2[$frequency][$atom][2];
    }
    print "\n";
  }
} #End sub grab_freq


#Finds any frequencies involving the concerted motion of a pair of atoms
#returns refences to arrays holding frequencies and intensities
#Checks to make sure normal mode vectors are in opposite directions and oriented along bond, and that change is at least X% of distance
#find_str($atom1, $atom2, \@freqs, \@intensities, \@vectors)
sub find_str {
  my ($atom1, $atom2, $freq_ref, $intens_ref, $vec_ref, $coords_ref, $debug) = @_;
  my @freqs = @{$freq_ref};
  my @intensities = @{$intens_ref};
  my @vectors = @{$vec_ref};
  my @coords = @{$coords_ref};

  my @str_freqs;
  my @str_intens;
  my $bond_vec = V($coords[$atom2][2] - $coords[$atom1][2],
                   $coords[$atom2][3] - $coords[$atom1][3],
                   $coords[$atom2][4] - $coords[$atom1][4]);
  $bond_vec /= $bond_vec->norm();

  foreach my $frequency (0..$#freqs) {
    my $vec1 = V(@{$vectors[$frequency][$atom1]});
    my $vec2 = V(@{$vectors[$frequency][$atom2]});

    my $change_vec = $vec1 - $vec2;
    my $percent_change = $change_vec->norm()/distance($atom1, $atom2, \@coords);
    if($vec1->norm() != 0 && $vec2->norm() != 0) {
      $vec1 /= $vec1->norm();
      $vec2 /= $vec2->norm();
  
      my $dot1 = $vec1 * $vec2;
      my $dot2 = $vec1 * $bond_vec;

      if($debug) {
        printf "%10.2f: %10.2f %10.2f %10.2f ", $freqs[$frequency],$dot1,$dot2,$percent_change;
      }
      if($dot1 < -0.7 && abs($dot2) > 0.7 && $percent_change > 0.2) {
        push(@str_freqs, $freqs[$frequency]);
        push(@str_intens, $intensities[$frequency]);
        if($debug) {
          print "YES";
        }
      }
      if($debug) {
        print "\n";
      }
    }
  }

  return (\@str_freqs, \@str_intens);
} #End find_str

#Returns jobIDs of all jobs (queued or running) in current directory, returns 0 if no jobs in queue match current directory
#This could be improved by searching for $Path more carefully!
#Works for PBS (default) or LSF
sub findJob {
  my $Path = $_[0];
  chomp($Path);

  #Strip leading directories off of $Path, to deal with different ways $HOME is treated
  $Path =~ s/^\S+$ENV{USER}//;
  
  my @jobIDs;

  if($queue_type eq 'LSF') {				#LSF queue
    my $bjobs=`bjobs -l 2> /dev/null`;
    #Combine into one line with no whitespace
    $bjobs =~ s/\s+//g;
    $bjobs =~ s/\r|\n//g;

    #First grab all jobs
    my @jobs = ($bjobs =~ m/(Job<\d+>.+?RUNLIMIT)/g);

    #parse each job looking for $Path
    foreach my $job (@jobs) {
      if ($job =~ /Job<(\d+)>\S+CWD<.+$Path>/) {
        push(@jobIDs,$1);
      }
    }
  }elsif ($queue_type eq 'PBS') {				#PBS (default)
    my $qstat = `qstat -fx`;
  
    #First grab all jobs
    my @jobs = ($qstat =~ m/<Job>(.+?)<\/Job>/g);
  
    #Grab jobIDs for all jobs matching $Path
    foreach my $job (@jobs) {
    	if ($job =~ m/<Job_Id>(\d+)\S+PBS_O_WORKDIR=\S+$Path</) {
    		push(@jobIDs, $1);
    	}
    }
  }elsif ($queue_type eq 'Slurm') {
    my @alljobs=`squeue -o %i_%Z`;
    foreach my $job (@alljobs) {
      if($job =~ /$Path/) {
        my @array = split(/_/, $job);
        push(@jobIDs,$array[0]);
      }
    }
  }
  
  if(@jobIDs) {
  	return @jobIDs;
  }
  return;
} #end sub findJob


#call as submit(com file, catalyst, walltime, nprocs, nocheck)
#Works for LSF or PBS (default)
sub submit {
  my ($com_file, $walltime, $numprocs, $template_job, $node) = @_;
  chomp(my $jobname=`basename $com_file .com`);
  chomp(my $Path = getcwd());

  #Don't submit if job is already running in current directory!
  my $job_running = 0;
  if(&findJob($Path)) {
    print "Job already running in directory of $jobname.com.  New job NOT submitted\n";
    $job_running = 1;
  }

  unless ($job_running) { 
    #if template_job were provide, use template_job to build job file
    if ($template_job->{job}) {

        my $job_file = $template_job->{job};

        my $job_found;

        my $template_pattern = TEMPLATE_JOB;

        open JOB_TEM, "<$job_file" or die "Cannot open $job_file:$!\n";
        my $job_content = do {local $/; <JOB_TEM>};
        close (JOB_TEM);

        $job_content =~ s/$template_pattern->{JOB_NAME}/$jobname/g && ($job_found = 1);
        $job_content =~ s/$template_pattern->{WALL_TIME}/$walltime/g;
        $job_content =~ s/$template_pattern->{N_PROCS}/$numprocs/g;
        $job_content =~ s/$template_pattern->{NODE_TYPE}/$node/g;
        #remove the formula part
        $job_content =~ s/&formula&\n(.*\n)*&formula&\n//g;

        for my $var (sort keys %{ $template_job->{formula} }) {
            my $var_value = eval($template_job->{formula}->{$var});
            $job_content =~ s/$var/$var_value/g;
        }

        if ($job_found) {
            open JOB, ">$jobname.job";
            print JOB $job_content;
            close (JOB);
        }
    }
    
    unless (-e "$jobname.job") {

      my $memory=$numprocs*120000000;
      my $mb;
      if($queue_type eq 'LSF') {
        $memory = 0.8*$numprocs*2*10**9/8;		#memory in words per job
        $mb = 2700;				#memory in MB per core + padding
      }

      open JOB, ">$jobname.job";

      if($queue_type eq 'LSF') {					#LSF queue
        print JOB "#BSUB -J $jobname\n";
        print JOB "#BSUB -o $jobname.job.%J\n";
        print JOB "#BSUB -L /bin/bash\n";
        print JOB "#BSUB -W $walltime:00\n";
        print JOB "#BSUB -M $mb\n";
        print JOB "#BSUB -R 'rusage[mem=$mb]'\n";
        print JOB "#BSUB -n $numprocs\n";
        print JOB "#BSUB -R 'span[ptile=$numprocs]'\n";
        print JOB "export g09root=$g09root\n";
        print JOB ". \$g09root/g09/bsd/g09.profile\n";
        print JOB "trap \"rm -r \$SCRATCH/\$LSB_JOBID\" 0 1 2 3 9 13 14 15\n";
        print JOB "mkdir \$SCRATCH/\$LSB_JOBID\n";
        print JOB "cd \$SCRATCH/\$LSB_JOBID\n";
        print JOB "echo -P- $numprocs > Default.Route\n";
        print JOB "echo -M- $memory >> Default.Route\n";
        print JOB "module purge\n";
        print JOB "env\n";
        print JOB "cp \$LS_SUBCWD/*.chk .\n";
        print JOB "g09  < \$LS_SUBCWD/$jobname.com  > \$LS_SUBCWD/$jobname.log\n";
        print JOB "cp *.chk \$LS_SUBCWD/\n";
        print JOB "exit\n";
      }elsif ($queue_type eq 'PBS') {					#PBS (default)
        print JOB "#!/bin/bash\n";
        print JOB "#PBS -l walltime=$walltime:00:00,mem=8gb,nodes=1:ppn=$numprocs\n\n\n";
        print JOB "export g09root=$g09root\n";
        print JOB ". \$g09root/g09/bsd/g09.profile\n\n";
        print JOB "trap \"\\rm -r \$TMPDIR/\$PBS_JOBID\" 0 1 2 3 9 13 14 15\n\n";
        print JOB "mkdir \$TMPDIR/\$PBS_JOBID\n";
        print JOB "cd \$TMPDIR/\$PBS_JOBID\n\n";
        print JOB "echo -P- $numprocs > Default.Route\n";
        print JOB "echo -M- $memory >> Default.Route\n\n";
        print JOB "module purge\n\n";
        print JOB "env\n\n";
        print JOB "cp \$PBS_O_WORKDIR/*.chk .\n";
        print JOB "g09  < \$PBS_O_WORKDIR/$jobname.com > \$PBS_O_WORKDIR/$jobname.log\n\n";
        print JOB "cp *.chk \$PBS_O_WORKDIR/\n\n";
        print JOB "exit";
      }elsif ($queue_type eq 'Slurm') {
        print JOB "#!/bin/bash\n";
        print JOB "#\n";
        print JOB "#SBATCH -J $jobname -e $jobname.job.e%j -o $jobname.job.o%j\n";
        if ($node) {
          print JOB "#SBATCH -p $node\n";
        }else {
          print JOB "#SBATCH -p medium\n";
        }
        print JOB "#SBATCH -t $walltime:00:00 -n $numprocs --mem=56G\n";
        print JOB "#SBATCH --ntasks-per-node=$numprocs\n";
        print JOB "\n";
        print JOB "\n";
        print JOB "cd \$TMPDIR\n";
        print JOB "\n";
        print JOB "df -h\n";
        print JOB "export g09root=/sw/group/lms/sw/g09_D01\n";
        print JOB ". $g09root/g09/bsd/g09.profile\n";
        print JOB "\n";
        print JOB "echo -P- 28 > Default.Route \n";
        print JOB "echo -M- 56GB >> Default.Route \n";
        print JOB "\n";
        print JOB "module purge \n";
        print JOB "\n";
        print JOB "env\n";
        print JOB "\n";
        print JOB "g09  <  \$SLURM_SUBMIT_DIR/$jobname.com  >  \$SLURM_SUBMIT_DIR/$jobname.log\n";
        print JOB "\n";
        print JOB "df -h\n";
        print JOB "\n";
        print JOB "ls -al\n";
        print JOB "\n";
        print JOB "\cp *.wf? *.47 \$LS_SUBCWD\n";
        print JOB "\n";
        print JOB "exit\n";
      }
      close(JOB);
    }

    #Alert user if qsub (or bsub) returns error
    #FIXME
    if($queue_type eq 'LSF') {
      if(system("bsub < $jobname.job >& /dev/null")) {
        print "Submission denied!\n";
        return 1;
      }
    } else {
      if(system("qsub $jobname.job -N $jobname >& /dev/null")) { 
        print "Submission denied!\n";
        return 1;
      }
    }
  }
  return 0;
} #end sub submit


#Makes job file for arbitrary number of com files (LSF)
sub make_job_file {
  my ($jobname, $walltime, $mb, $numprocs, $memory, @comfiles) = @_;

  open JOB, ">$jobname.job" or die "Can't open $jobname.job";
  print JOB "#BSUB -J $jobname\n";
  print JOB "#BSUB -o $jobname.job.%J\n";
  print JOB "#BSUB -L /bin/bash\n";
  print JOB "#BSUB -W $walltime:00\n";
  print JOB "#BSUB -M $mb\n";
  print JOB "#BSUB -R 'rusage[mem=$mb]'\n";
  print JOB "#BSUB -n $numprocs\n";
  print JOB "#BSUB -R 'span[ptile=$numprocs]'\n";
  print JOB "export g09root=/software/lms/g09_D01\n";
  print JOB ". \$g09root/g09/bsd/g09.profile\n";
  print JOB "trap \"rm -r \$SCRATCH/\$LSB_JOBID\" 0 1 2 3 9 13 14 15\n";
  print JOB "mkdir \$SCRATCH/\$LSB_JOBID\n";
  print JOB "cd \$SCRATCH/\$LSB_JOBID\n";
  print JOB "echo -P- $numprocs > Default.Route\n";
  print JOB "echo -M- $memory >> Default.Route\n";
  print JOB "module purge\n";
  print JOB "env\n";
  print JOB "cp \$LS_SUBCWD/*.chk .\n";
  foreach (@comfiles) {
    print JOB "g09  < \$LS_SUBCWD/$_.com  > \$LS_SUBCWD/$_.log\n";
  }
  print JOB "cp *.chk \$LS_SUBCWD/\n";
  print JOB "exit\n";

  close(JOB);
} #End make_job_file


#This just moves the two selected atoms!!!! Can be dangerous if large changes are made.
#TODO: write another version that moves two molecular fragments until two atoms are a selected distance apart
#if $flag is set, only moves $atom2 rather than moving atoms symmetrically!
sub change_distance {
  my ($atom1, $atom2, $distance, $coords_ref, $flag) = @_;
  my @coords = @{$coords_ref};

  my $current_distance = &distance($atom1, $atom2, \@coords);
  my $difference = $distance - $current_distance;
  my @displacement;
  foreach my $i (2..4) {
    if($flag) {
      $coords[$atom2][$i] -= $difference*($coords[$atom1][$i] - $coords[$atom2][$i])/($current_distance);
    } else {
      $coords[$atom1][$i] += $difference*($coords[$atom1][$i] - $coords[$atom2][$i])/(2*$current_distance);
      $coords[$atom2][$i] -= $difference*($coords[$atom1][$i] - $coords[$atom2][$i])/(2*$current_distance);
    }
  }
  return;
} #end sub change_distance

#Get distance between two atoms
#distance($atom1, $atom2, $coords_ref);
#distance($atom1, $atom2, $coords_ref, $coords_ref2);	#if two coords_refs given, pulls atom1 from coords_ref and atom2 from coords_ref2
sub distance {
  my ($atom1, $atom2, $coords_ref, $coords_ref2) = @_;
  my @coords1 = @{$coords_ref};
  my @coords2;
  if($coords_ref2) {
    @coords2 = @{$coords_ref2};
  } else {
    @coords2 = @{$coords_ref};
  }

  my $x1 = $coords1[$atom1][2];
  my $y1 = $coords1[$atom1][3];
  my $z1 = $coords1[$atom1][4];
  
  my $x2 = $coords2[$atom2][2];
  my $y2 = $coords2[$atom2][3];
  my $z2 = $coords2[$atom2][4];
  
  my $distance = sqrt(($x1-$x2)**2 + ($y1-$y2)**2 + ($z1-$z2)**2);
  return $distance;
}


#calculate angle angle for atoms 1-3
#angle(atom1, atom2, atom3, ref_to_coords)
sub angle {
  my ($atom1, $atom2, $atom3, $coords_ref) = @_;
  my @coords = @{$coords_ref};

  my $v21 = V($coords[$atom1][2] - $coords[$atom2][2], $coords[$atom1][3] - $coords[$atom2][3], $coords[$atom1][4] - $coords[$atom2][4]);
  my $v23 = V($coords[$atom3][2] - $coords[$atom2][2], $coords[$atom3][3] - $coords[$atom2][3], $coords[$atom3][4] - $coords[$atom2][4]);

  $v21 /= $v21->norm();
  $v23 /= $v23->norm();

  my $angle = acos($v21 * $v23);
  return (rad2deg($angle));
}

#calculate dihedral angle for atoms 1-4
#dihedral(atom1, atom2, atom3, atom4, ref_to_coords)
sub dihedral {
  my ($atom1, $atom2, $atom3, $atom4, $coords_ref) = @_;
  my @coords = @{$coords_ref};

  my $v12 = V($coords[$atom1][2] - $coords[$atom2][2], $coords[$atom1][3] - $coords[$atom2][3], $coords[$atom1][4] - $coords[$atom2][4]);
  my $v23 = V($coords[$atom2][2] - $coords[$atom3][2], $coords[$atom2][3] - $coords[$atom3][3], $coords[$atom2][4] - $coords[$atom3][4]);
  my $v34 = V($coords[$atom3][2] - $coords[$atom4][2], $coords[$atom3][3] - $coords[$atom4][3], $coords[$atom3][4] - $coords[$atom4][4]);

  my $dihedral = atan2((($v12 x $v23) x ($v23 x $v34)) * $v23/abs($v23), ($v12 x $v23) * ($v23 x $v34));

  return rad2deg($dihedral);
}

#changes dihedral about bond atom1-atom2 by an angle angle_change (in radians!)
#change_dihedral(atom1, atom2, angle_change, ref_to_coords
sub change_dihedral {
  my ($atom1, $atom2, $angle_change, $coords_ref) = @_;

  my @coords = @{$coords_ref};
  my @connected = get_connected(\@coords);
  my @connected_atoms = get_all_connected($atom1, $atom2, @connected);

  if($#connected_atoms < 0) {
    print "Cannot change this dihedral angle...";
    return 0;
  }

  my @axis = ($coords[$atom2][2] - $coords[$atom1][2],$coords[$atom2][3] - $coords[$atom1][3],$coords[$atom2][4] - $coords[$atom1][4]);
  center_genrotate($atom1, @axis, $angle_change, \@coords, @connected_atoms);

  return;
}  #End sub change_dihedral


#get quota and return summary message
sub get_quota {
  my $limit;
  my $used;
  my $quota;
  my $message;

  if($queue_type eq 'LSF') {
    $quota = `/software/tamusc/local/bin/showquota`;
  } else {
    $quota = `/usr/lpp/mmfs/bin/mmlsquota`;
  }
  if($queue_type eq 'LSF') {
    if($quota =~ /scratch\s+(\S+)\s+(\S+)/) {
      ($limit, $used) = ($2, $1);
    }
    $message = sprintf "scratch used = $used of $limit";
  } else {
    if($quota =~ /scratch\s+\S+\s+(\d+)\s+(\d+)/ ) {
      ($limit, $used) = ($2, $1);
      $used /= 1048576;
      $limit /= 1048576;
      my $percent = 100*$used/$limit;
      $message = sprintf "scratch used = %.0f of %.0f GB (%0.f%%)", $used, $limit, $percent;
      if($percent > 90) {
        $message .= " (Dangerously full!)";
      }
    }
  }

  return $message;
} #end sub get_quota


#coord_shift--shifts atoms first_atom through last_atom along X, Y, Z coords
sub coord_shift {
  my ($X, $Y, $Z, $coords_ref, @targets) = @_;
  my @coords = @{$coords_ref};
  
  #if no target atoms specified, rotate all atoms
  if($#targets == -1) {
    @targets = (0..$#coords);
  }

  foreach my $atom (@targets) {
    $coords[$atom][2] += $X;
    $coords[$atom][3] += $Y;
    $coords[$atom][4] += $Z;
  }
} #sub coord_shift

sub remove_atoms {                                #remove atoms from atom list and return the removed atoms
  my ($coords_ref, @targets) = @_;
  my @removed;
  #sort targets so that atoms can be removed from last to first
  @targets = sort {$b <=> $a} @targets;

  foreach my $target (@targets) {
   my @temp = splice(@{$coords_ref}, $target, 1);
   push @removed, @temp;
  }
  return @removed;
} #End remove_atoms

#Removes fragment by cutting $atom1-$atom2 bond and removing $atom2 and all connected atoms, replaces with hydrogen along original $atom1-$atom2 bond.
#remove_fragment(\@coords, $atom1, $atom2);
sub remove_fragment {
  my ($coords_ref, $atom1, $atom2) = @_;

  my @coords = @{$coords_ref};

  my @connected = get_connected($coords_ref);
  #get all atoms connected to $atom2, avoiding $atom1 (automatically checks if $atom2 and $atom1 are in same cyclic substructure)
  my @connected_atoms = get_all_connected($atom2, $atom1, @connected);

  #remove $atom2 from list of connected atoms (so it can be converted to a hydrogen
  my $index_of_atom2 = 0;
  foreach my $atom (0..$#connected_atoms) {
    if($connected_atoms[$atom] == $atom2) {
      $index_of_atom2 = $atom;
    }
  }
  #replace $atom1 with hydrogen and scale bond length to sum of covalent radii
  $coords[$atom2][0] = 'H';
  my $new_distance = $radii{$coords[$atom1][0]} + $radii{$coords[$atom2][0]};
  change_distance($atom1, $atom2, $new_distance, $coords_ref, 1);

  splice(@connected_atoms, $index_of_atom2, 1);
  remove_atoms($coords_ref, @connected_atoms);

  return;
} #end sub remove_fragment


#Comb through G09 log file and look for reason for failure
#returns name of error or "UNKOWN" if can't find anything obvious
sub get_error {
  my $filename = $_[0];

  my %errors = (
    "CHK"   => "NtrErr Called from FileIO",			#delete
    "EIGEN" => "Wrong number of Negative eigenvalues", 		#opt=noeigen
    "CONV"  => "Convergence failure -- run terminated.", 	#scf=xqc
    "QUOTA" => "Erroneous write", 				#check quota and alert user; REMOVE error from end of file!
    "CLASH" => "Atoms too close", 				#flag as CLASH
    "CHARGEMULT" => "The combination of multiplicity",		#die and alert user to check catalyst structure or fix reaction_data!
    "REDUND" => "Bend failed for angle"                         #Using opt=cartesian
  );

  if(-e $filename) {
    open INFILE, "$filename" or die "Can't open $filename!";
    while (<INFILE>) {
      foreach my $error (keys %errors) {
        if($_ =~ /$errors{$error}/) {
          close(INFILE);
          return $error;
        } 
      }
    }
    close(INFILE);
    return "UNKNOWN";
  } 
  return "NOFILE";
} #end sub get_error

#Compute connected atoms for each atom in @coords
#return 2-D array, where each row is a list of connected atoms [connected atom1, connected atom2, ...]
#NB: this is DIFFERENT from get_connectivity!!!
#This is used by shortest_path
sub get_connected {
  my ($coords_ref) = @_;
  my @coords = @{$coords_ref};
  my $tolerance = 0.2;

  my @connected;
  foreach my $atom1 (0..$#coords) {
    my @row;
    foreach my $atom2 (0..$#coords) {
      my $distance = &distance($atom1, $atom2, \@coords);
      my $cutoff = $radii{$coords[$atom1][0]} + $radii{$coords[$atom2][0]} + $tolerance;
      if($distance < $cutoff && $atom1 != $atom2) {
        push(@row, $atom2);
      }
    }
    push(@connected, \@row);
  }
  return @connected;
} #end of get_connected


#Compute connectivity for each atom (ie: number of atoms bound to that atom)!
#return 2-D array, where each row is [element, number of connected atoms]
sub get_connectivity {
  my ($coords_ref) = @_;
  my @coords = @{$coords_ref};
  my $tolerance = 0.2;

  my @connectivity; 
  foreach my $atom1 (0..$#coords) {
    my $connected = 0;
    foreach my $atom2 (0..$#coords) {
      my $distance = &distance($atom1, $atom2, \@coords);
      my $cutoff = $radii{$coords[$atom1][0]} + $radii{$coords[$atom2][0]} + $tolerance;
      if($distance < $cutoff && $atom1 != $atom2) {
        $connected++;
      } 
    }
    push(@connectivity, [$coords[$atom1][0], $connected]);
  }
  return @connectivity;
} #end of get_connectivity


#check connectivity array for hypervalent atoms
#needs @connectivity from get_connectivity
#return offending atom number if problems detected (if multiple problems, only returns first one found!
#return -1 if no problems
sub check_connectivity {
  my @connectivity = @_;
  #maximum connectivity
  my %max_connectivity = ( H => '1', B => '4', C => '4', N => '4', O => '2', F => '1', Si => '6', P => '4', S => '4', Cl => '1', I => '1', Br => '1', X => '1000');

  foreach my $atom1 (0..$#connectivity) {
    if($connectivity[$atom1][1] > $max_connectivity{$connectivity[$atom1][0]}) {
      return $atom1;
    }
  }
  return -1;
} #end of check_connectivity


#Simple rotations around Cartesian axes of all coords
#rotate($axis, $angle_in_radians, \@coords)
#Returns 0 if it fails
sub rotate {
  my ($axis, $angle, $coords_ref) = @_;
  my @coords = @{$coords_ref};
  switch ($axis) {
    case /[Xx]/ {		
      genrotate(1,0,0, $angle, \@coords);
      return 1;
    }
    case /[Yy]/ {	
      genrotate(0,1,0, $angle, \@coords);
      return 1;
    }
    case /[Zz]/ {	
      genrotate(0,0,1, $angle, \@coords);
      return 1;
    } else {
      die "Attempt at rotation other than around Cartesian axes!";
      return 0;
    }
  }
} #End sub rotate

#rotation about axis given by k_x, k_y, k_z (nb: this is now used by rotate)
#Uses quaternion rotation to perform rotation about arbitrary axis (vector fomulation of Euler-Rodrigues formula)
#genrotate(k_x, k_y, k_z, angle_in_radians, ref_to_coords, @list_of_atoms_to_rotate)
#If no atom list is given, applies rotation to all atoms
sub genrotate {
  my ($x, $y, $z, $angle, $coords_ref, @targets) = @_;
  my @coords = @{$coords_ref};
  my $numatoms = $#coords;

  #if no target atoms specified, rotate all atoms
  if($#targets == -1) {
    @targets = (0..$#coords);
  }

  my $a = cos($angle/2);
  my $w = V($x, $y, $z);

  #normalize rotation axis
  $w /= $w->norm();
  $w *= sin($angle/2);
  
  quat_rot($a, $w->[0], $w->[1], $w->[2], $coords_ref, @targets);
  return;
} #End sub genrotate


#Performs quaternion rotation; assumes quaternions already normalized!
#quat_rot(a, b, c, d, ref_to_coords, @atoms to be rotated)
sub quat_rot {
  my ($a, $b, $c, $d, $coords_ref, @targets) = @_;
  my @coords = @{$coords_ref};
  
  my $w = V($b, $c, $d);
  #rotate coords for each atom in target list
  foreach my $atom (@targets) {
    my $vec = V($coords[$atom][2], $coords[$atom][3], $coords[$atom][4]);
    my $wx = $w x $vec;
    my $new_vec = $vec + 2*$a*$wx + 2*($w x $wx);
    $coords[$atom][2] = $new_vec->[0];
    $coords[$atom][3] = $new_vec->[1];
    $coords[$atom][4] = $new_vec->[2];
  }
  
  return;
} #End quat_rot

#shifts atom to origin, performs genrotate, then shifts back to original origin
#center_genrotate(atom, k_x, k_y, k_z, angle_in_radians, ref_to_coords, @list_of_atoms_to_rotate)
sub center_genrotate {
  my ($atom, $x, $y, $z, $angle, $coords_ref, @targets) = @_;
  my @coords = @{$coords_ref};

  my @center = @{$coords[$atom]};
  coord_shift(-$center[2], -$center[3], -$center[4], $coords_ref);
  genrotate($x, $y, $z, $angle, \@coords, @targets);
  coord_shift($center[2], $center[3], $center[4], $coords_ref);

  return;
} #End sub center_genrotate




#centers molecule so that average of atom1, atom2, and atom3 is at origin, and orients molecule so that atom1 is on x-axis and atom2 and atom3 are as close as possible to XY-plane
#center_ring($atom1, $atom2, $atom3, $ref_to_coords);
sub center_ring {
  my ($atom1, $atom2, $atom3, $coords_ref) = @_;

  my @coords = @{$coords_ref};

  my @COM = (0, 0, 0);
  for my $i (0..2) {
    $COM[$i] += $coords[$atom1][$i+2];
    $COM[$i] += $coords[$atom2][$i+2];
    $COM[$i] += $coords[$atom3][$i+2];
  }
  for my $i (0..2) {
    $COM[$i] /= 3;
  }
  
  #shift geom to COM
  coord_shift(-$COM[0], -$COM[1], -$COM[2], \@coords);
  
  #Put atom1 along x-axis
  my $v1 = V($coords[$atom1][2], $coords[$atom1][3], $coords[$atom1][4]);
  my $vx = V(1,0,0);
  my $cross1 = $v1 x $vx;
  my $angle1 = atan2($v1, $vx);
  genrotate($cross1->[0], $cross1->[1], $cross1->[2], -$angle1, \@coords);
  printXYZ(\@coords);
  
  #Now put atom2 and atom3 in XY-plane (or as close as possible)
  my $v2 = V(0, $coords[$atom2][3], $coords[$atom2][4]);
  my $v3 = V(0, $coords[$atom3][3], $coords[$atom3][4]);
  my $vz = V(0,0,1);
  my $cross2 = $v2 x $vz;
  my $cross3 = $v3 x $vz;
  my $angle2;
  if($v2->norm() != 0) {
    $angle2 = asin($cross2->norm()/$v2->norm());
  } else {
    $angle2 = pi()/2;
  }
  my $angle3;
  if($v3->norm() != 0) {
    $angle3 = asin($cross3->norm()/$v3->norm());
  } else {
    $angle3 = pi()/2;
  }
  my $av_angle = pi()/2 - ($angle2*2)/2;
  my $test = 180*$av_angle/pi();

  genrotate(1, 0, 0, -$av_angle, \@coords);
  return;
}


#Replaces atom $target with substituent $sub (name of XYZ file)
#Returns coords in same orientation as given, with atom ordering preserved
#&substitute($target, $sub, \@coords, $no_min);
#TO DO: after adding substituent, rotation added atoms around new bond to maxmimize distance to all other atoms!
sub substitute {
  my ($target, $sub, $coords_ref, $degree, $no_min) = @_;
  my @coords = @{$coords_ref};
  my $nearest_neighbor;
  #get connected atoms and check that atom to be replaced is monovalent
  my @connected = &get_connected(\@coords);
  my @nearest = @{$connected[$target]};
  my @targets;
  if($#nearest > 0) {
    print "Target is a substituent group!\n";
    #get the nearest connected atoms in the main body
    my $min = 999;
    for (@nearest) {
      my @get_all = get_all_connected($target, $_, \@connected);
      if ($#get_all < $min) {
        $min = $#get_all;
        $nearest_neighbor = $_;
        @targets = @get_all;
      }
    }
  }else{ 
    $nearest_neighbor = $connected[$target][0];
    @targets = ($target);
  }
  #get nearest neighbor from @connected

  #grab substituent coords
  my @sub_coords = &grab_coords("$AARON/Subs/$sub.xyz");

  #Rotate to align along nearest_neighbor-target bond then shift sub_coords to nearest_neighbor position
  my @nearest_neighbor_position = ($coords[$nearest_neighbor][2], $coords[$nearest_neighbor][3], $coords[$nearest_neighbor][4]);
  my $bond_axis = V($coords[$target][2] - $coords[$nearest_neighbor][2], 
                    $coords[$target][3] - $coords[$nearest_neighbor][3], 
                    $coords[$target][4] - $coords[$nearest_neighbor][4]); 
  $bond_axis /= $bond_axis->norm();

  #sub_coords are aligned along x-axis, so find rotation axis that transforms x-axis to bond_axis
  my $v_x = V(1,0,0);
  my $cross = $v_x x $bond_axis;
  my $angle = atan2($bond_axis, $v_x);
  
  genrotate($cross->[0], $cross->[1], $cross->[2], $angle, \@sub_coords);
  coord_shift(@nearest_neighbor_position, \@sub_coords);
  
  #replace target with first atom of substituent
  splice(@coords, shift(@targets), 1, shift(@sub_coords));
  #if any element remaining in the @targets, remove them
  while (@targets) {
     splice(@coords,shift(@targets),1);
     foreach my $i (0..$#targets) {
        $targets[$i]--;
     }
  }
  #build list of substituent coords
  my @targets_tosub = ($target);
  foreach my $atom (0..$#sub_coords) {
    push(@targets_tosub, $#coords + $atom + 1);
  } 

  if($#sub_coords >= 1) {
     #Add remainder of substituent coords
     push(@coords, @sub_coords);
     if (not defined $no_min) {
        minimize_torsion($nearest_neighbor, $target, \@coords, @targets_tosub);
     }

     if ($degree) {
       my $axis = V($coords[$target][2] - $coords[$nearest_neighbor][2],
                    $coords[$target][3] - $coords[$nearest_neighbor][3],
                    $coords[$target][4] - $coords[$nearest_neighbor][4]);
       center_genrotate($nearest_neighbor, $axis->[0], $axis->[1], $axis->[2], $degree, \@coords, @targets_tosub); 
     }
  }
  return (\@coords, \@targets_tosub);
} #End sub substitute

#rotate substitute
#sub_rotate(target,angle,\@coords)
sub sub_rotate{
   my ($target, $angle, $coords) = @_;
   my @coords = @{$coords};
   my @connected = &get_connected(\@coords);
   my @nearst = @{$connected[$target]};
   #now we need to define which nearst atom is the target.
   my $min = 999;
   my $atom2;
   my @targets;
   for (@nearst) {
      my @get_all = get_all_connected($target, $_, \@connected);
      if (@get_all && $#get_all < $min) {
        $min = $#get_all;
        $atom2 = $_;
        @targets = @get_all;
      }
   }
   my $axis = V($coords[$target][2] - $coords[$atom2][2],
                $coords[$target][3] - $coords[$atom2][3],
                $coords[$target][4] - $coords[$atom2][4]);
   center_genrotate($target, $axis->[0], $axis->[1], $axis->[2], $angle, \@coords, @targets);
   return @coords;
    
} #End sub sub_rotate

#Replaces atoms $target1 and $target2 with six-membered fused ring
#Possible types: Ar (aromatic) LD_chair, LU_chair, D_boat, U_boat
#Returns coords in same orientation as given, with atom ordering preserved
#fused_ring($target1, $target2, $type, @coords);
sub fused_ring {
  my ($target1, $target2, $type, $coords_ref) = @_;
  my @coords = @{$coords_ref};

  #get connected atoms
  my @connected = &get_connected(\@coords);
  if($#{$connected[$target1]} > 0 || $#{$connected[$target2]} > 0) {
    print "Trying to substitute non-monovalent atom!\n";
    return 0;
  } 

  #Figure out what needs to be added to match fused ring
  my $path_length = &shortest_path($target1, $target2, @connected);
  if($path_length < 3 || $path_length > 5) {
    print "Can't figure out how to build fused ring connecting those two atoms...\n";
    return 0;
  }
  #check to make sure known type
  if($type !~ /a_pinene/ && $type !~ /LD_chair/ && $type !~ /LU_chair/ && $type !~ /D_boat/ && $type !~ /U_boat/ && $type !~ /Ar/) {
    print "Unknown ring type!\n";
    return 0;
  }

  my @ring_coords;

  switch($type) {
    case /Ar/ {
      @ring_coords = &grab_coords("$AARON/Ring_fragments/six.$path_length.xyz");
    }
    else {	#Chairs
      if($path_length != 3) {
        print "Can't figure out how to build fused ring connecting those two atoms...\n";
        return 0;
      } else {
        @ring_coords = &grab_coords("$AARON/Ring_fragments/Chairs/$type.xyz");
      }
    }
  }
  
  #get nearest neighbors from @connected
  my $nearest_neighbor1 = $connected[$target1][0];
  my $nearest_neighbor2 = $connected[$target2][0];

  #shift coords so that nearest_neighbor1 is at the origin
  my @origin = @{$coords[$nearest_neighbor1]};
  &coord_shift(-$origin[2], -$origin[3], -$origin[4], \@coords);

  #Orient so that nearest_neighbor2 is along x-axis
  my $v1 = V($coords[$nearest_neighbor2][2], $coords[$nearest_neighbor2][3], $coords[$nearest_neighbor2][4]);
  my $vx = V(1,0,0);
  my $cross1 = $v1 x $vx;
  my $angle1 = atan2($v1, $vx);
  genrotate($cross1->[0], $cross1->[1], $cross1->[2], $angle1, \@coords); #accounts for sign of dot product to get sign of angle right!

  #final rotation around x-axis to bring $target1 and $target2 into XY-plane
  my $chi1 = atan2($coords[$target1][4],$coords[$target1][3]);
  my $chi2 = atan2($coords[$target2][4],$coords[$target2][3]);
  my $chi = ($chi1 + $chi2)/2;
  &rotate('x',-$chi,\@coords);
  
  #replace target1 with 1st ring atom
  splice(@coords, $target1, 1, shift(@ring_coords));
  #replace target2 with 2nd ring atom
  splice(@coords, $target2, 1, shift(@ring_coords));
  #Add remainder of ring coords
  push(@coords, @ring_coords);
  
  #Return geometry to original orientation/position
  &rotate('x',$chi,\@coords);
  genrotate($cross1->[0], $cross1->[1], $cross1->[2], -$angle1, \@coords);
  &coord_shift($origin[2], $origin[3], $origin[4], \@coords);

  return @coords;
} #End sub fused_ring


#Find shortest path between atom1 and atom2
#Performs a breadth first search, returns length of shortest path
#returns -1 if no path found
sub shortest_path {
  my ($atom1, $atom2, @connected) = @_;
  my @positions = ($atom1);
  my %visited = ( $atom1 => '0') ;	#keys are numbers of the visited atoms, values are the corresponding depths

  #loop until @positions is empty
  while(@positions) {
    my $position = shift(@positions);
    if ($position eq $atom2) {
      return $visited{$position};
    }
    foreach my $atom (@{$connected[$position]}) { 	#if not found yet, grab all atoms connected to current atom and add to queue (unless already visited)
      if($atom =~ /\d/ && !exists $visited{$atom}) {	#skip over element, just add atom numbers
        push(@positions, $atom);
        $visited{$atom} = $visited{$position} + 1;
      }
    }
  }
  return -1;	#return -1 if no path found in BFS
} #end shortest_path


#Returns list of all atoms connected to a start_atom without passing through avoid_atom
#Returns empty list if start_atom and avoid_atom are part of a cyclic structure
#get_all_connected(start_atom, avoid_atom, @connected)
sub get_all_connected {
  my ($start_atom, $avoid_atom, $connected) = @_;
  my @connected = @{$connected};
  my @connected_temp = map { [@$_] } @connected;
  #remove avoid_atom from list of atoms connected to $start_atom
  my $avoid_atom_found = 0;
  foreach my $atom (0..$#{$connected_temp[$start_atom]}) {
    if($connected_temp[$start_atom][$atom] == $avoid_atom) {
      $connected_temp[$start_atom][$atom] = -1;
      $avoid_atom_found = 1;
    }
  }

  my @positions = ($start_atom);
  #I can probably use something simpler than a hash here (since I'm not really using the values)
  my %visited = ( $start_atom => '0') ;	#keys are numbers of the visited atoms, values are not used

  #loop until @positions is empty
  while(@positions) {
    my $position = shift(@positions);
    foreach my $atom (@{$connected_temp[$position]}) { 	#grab all atoms connected to current atom and add to queue (unless already visited)
      if($atom >= 0 && !exists $visited{$atom}) {
        push(@positions, $atom);
        $visited{$atom} = 0;
      }
    }
  }

  my @all_connected_atoms = keys %visited;
  #Check all_connected_atoms for avoid_atom
  foreach (@all_connected_atoms) {
    if($_ == $avoid_atom) {
      return ();
    }
  }
  #change the start_atom in the @all_connected_atoms to the first element
  if ($all_connected_atoms[1]) {
    @all_connected_atoms = grep {$_ != $start_atom} @all_connected_atoms;
    @all_connected_atoms = sort @all_connected_atoms;
    unshift @all_connected_atoms, $start_atom;
  }
  return @all_connected_atoms;
} #End sub get_all_connected


#Returns summary of the last optimization step of an optimization
#TODO: Modify to fill progress hash with numerical values between 0 and 100 based on progress towards four convergence criteria
sub get_gradient {
  my $file = $_[0];

  my $maxforce;
  my $rmsforce;
  my $maxdisp;
  my $rmsdisp;
  my @converged; 	#array to hold "YES/NO" for each of the four criteria

  if(-e $file) {
    open INFILE, "<$file";
    while (<INFILE>) {
      if($_ =~ /Threshold  Converged/) {
        my $line = <INFILE>;
        if($line =~ /Maximum Force\s+(\d\.\d+)\s+\d+\.\d+\s+(\S+)/) {
          $maxforce = $1;
          push(@converged, $2);
        }
        $line = <INFILE>;
        if($line =~ /RMS     Force\s+(\d\.\d+)\s+\d+\.\d+\s+(\S+)/) {
          $rmsforce = $1;
          push(@converged, $2);
        }
        $line = <INFILE>;
        if($line =~ /Maximum Displacement\s+(\d\.\d+)\s+\d+\.\d+\s+(\S+)/) {
          $maxdisp = $1;
          push(@converged, $2);
        }
        $line = <INFILE>;
        if($line =~ /RMS     Displacement\s+(\d\.\d+)\s+\d+\.\d+\s+(\S+)/) {
          $rmsdisp = $1;
          push(@converged, $2);
        }
      }
    }
    close INFILE;
  } else {
    return 0;		#file not found!
  }

  if(@converged) {
    return "Max Force: $maxforce ($converged[0]), RMS Force: $rmsforce ($converged[1]), Max Disp: $maxdisp ($converged[2]), RMS Disp: $rmsdisp ($converged[3])";
  } else {
    return "No steps yet...";
  }
} #end sub get_gradient

sub get_homo_lumo {
  my ($infile) = @_;
  open INFILE, $infile or die "Can't open $infile";

  my @occ;
  my @vir;

  while (<INFILE>) {
    chomp;
    if($_ =~ /The electronic state is/) {
      @occ = [];
      @vir = [];
      my $line = '';
      while ($line !~ /Condensed to atoms/) {
        $line = <INFILE>;
        chomp($line);
        if($line =~ /Alpha  occ\. eigenvalues --\s+(.+)$/) {
          push(@occ, split(/\s+/, $1));
        } elsif($line =~ /Alpha virt\. eigenvalues --\s+(.+)$/) {
          push(@vir, split(/\s+/, $1));
        }
      }
    }
  }
  return ($occ[-1], $vir[1]);
} #end get_homo_lumo

#Returns energy from G09 log file
sub get_energy {
  my ($infile) = @_;
  open INFILE, $infile or die "Can't open $infile";
  
  my $energy = 0;
  while (<INFILE>) {					#read in data from Gaussian log file
    chomp;
    if($_ =~ /SCF Done/o) {
    	my @array= split(/\s+/o,$_);
    	$energy = $array[5];
    }
  }
  close(INFILE);
  return $energy;
} #end get_energy

#returns energy, enthalpy, G, and grimme_G (grimme_G is calculated; the rest are read from teh log file)
#Based on Gaussian Thermochemistry white paper, calculate: E + ZPVE, H (0K), H(Temp), G(Temp) (http://www.gaussian.com/g_whitepap/thermo/thermo.pdf)
sub get_thermo {
  my ($infile,$T) = @_;

  my $v0 = 100; #cutoff for quasi-RRHO (in cm-1)
  my $P = 1*101317; #1 ATM
 
  my $energy; 
  my $mass; 						#molecular mass
  my $mult; 						#Spin multiplicity
  my @rottemps; 					#rotational temperatures (Kelvin)
  my @vibtemps; 					#vibrational temperatures (Kelvin)
  my @vibfreqs; 					#vibrational frequencies (cm-1)
  
  my $sigmar; 						#Rotational symmetry number
  my $ZPVE; 						#Zero point energy read from G09 output
  my $Elec_zpve; 					#Sum of electronic and zero-point Energies read from G09 output
  my $enthalpy;  					#Sum of electronic and thermal enthalpies read from G09 output
  my $G;  						#Sum of electronic and thermal Free Energies read from G09 output
  
  open INFILE, $infile or die "Can't open $infile";
  while (<INFILE>) {					#read in data from Gaussian log file
    my $line = $_;
    if($line =~ / Harmonic frequencies/) { 		#reset everything if new frequency calculation is found
      @vibtemps = ();
      @vibfreqs = ();
    }
    if($line =~ /^ Frequencies --/) {
            chomp;
      $line =~ s/^\s+//;
      my @array= split(/\s+/o,$line);
            foreach my $freq (2..$#array) {
        if($array[$freq] > 0 ) {
          push(@vibfreqs, $array[$freq]);
          push(@vibtemps, $array[$freq]*$c*$h/$kb);
        }
      }
    }
    if($_ =~ /SCF Done/o) {
    	my @array= split(/\s+/o,$_);
    	$energy = $array[5];
    }
    if($line =~ /^ Rotational constants \(GHZ\):\s+(\S+)\s+(\S+)\s+(\S+)/) {
      @rottemps = ($1, $2, $3);
      foreach my $rot (0..$#rottemps) {
        $rottemps[$rot] *= $h*(10**9)/($kb);
      }
    }
    if($line =~ /^ Molecular mass:\s+(\S+)/) {
      $mass = $1*$amu2kg;
    }
    if($line =~ /^ Sum of electronic and zero-point Energies=\s+(\S+)/) {
      $Elec_zpve = $1;
    } 
    if($line =~ /^ Sum of electronic and thermal Enthalpies=\s+(\S+)/) {
      $enthalpy = $1;
    }
    if($line =~ /^ Sum of electronic and thermal Free Energies=\s+(\S+)/) {
      $G = $1;
    }
    if($line =~ /^ Zero-point correction=\s+(\S+)/) {
      $ZPVE = $1;
    }
    if($line =~ / Multiplicity = (\d+)/) {
      $mult = $1;
    }
    if($line =~ / Rotational symmetry number\s+(\d+)/) {
      $sigmar = $1;
    }
  } #end read G09 log file
  close(INFILE);
  unless ($enthalpy) {
    return ($energy);
  }
  #Calculate average moment of inertia for Grimme's quasi-RRHO approach
  my $Bav = ($h**2/(24*pi**2*$kb))*(1/$rottemps[0] + 1/$rottemps[1] + 1/$rottemps[2]);
  
  if ($#rottemps!=2) {
    die "Problem reading Rotational constants";
  }
  
  #Translational component of Entropy
  my $qt = (2*pi*$mass*$kb*$T/($h*$h))**(3/2)*$kb*$T/$P;
  my $St = $R*(log($qt) + 5/2);
  
  #Translation component of Energy
  my $Et = 3*$R*$T/2;
  
  #Electronic component of Entropy
  my $Se = $R*(log($mult));
  
  #Rotational component of Entropy
  my $qr = (sqrt(pi)/$sigmar)*($T**(3/2)/sqrt($rottemps[0]*$rottemps[1]*$rottemps[2]));
  my $Sr = $R*(log($qr) + 3/2);
  
  #Rotational component of Energy
  my $Er = 3*$R*$T/2;
  
  #Vibrational component of Entropy and Energy
  my $Ev = 0;
  my $Sv = 0;
  my $Sv_quasiRRHO = 0;
  
  foreach my $i (0..$#vibtemps) {
    my $Sv_temp = $vibtemps[$i]/($T*(exp($vibtemps[$i]/$T)-1)) - log(1-exp(-$vibtemps[$i]/$T));
    
    $Sv += $Sv_temp;
      $Ev += $vibtemps[$i]*(1/2 + 1/(exp($vibtemps[$i]/$T) - 1));
    
    #calculate quasi-RRHO contribution to Sv
    my $mu = $h/(8*pi**2*$vibfreqs[$i]*$c);
    my $mu_prime = $mu*$Bav/($mu + $Bav);
    my $Sr = 1/2 + log(sqrt(8*pi**3*$mu_prime*$kb*$T/$h**2));
    
    my $weight = 1/(1+($v0/$vibfreqs[$i])**4);
    
    $Sv_quasiRRHO += $weight*$Sv_temp + (1-$weight)*$Sr;
  }
  
  $Sv *= $R;
  $Ev *= $R;
  $Sv_quasiRRHO *= $R;
  
  #Grab Electronic energy from $Elec_zpve and $ZPVE
  my $E_e = $Elec_zpve - $ZPVE;
  my $Etot = $Et + $Er + $Ev;
  my $Hcorr = $Etot + $R*$T;
  my $Stot = $St + $Sr + $Sv + $Se;
  my $Stot_quasiRRHO = $St + $Sr + $Sv_quasiRRHO + $Se;
  my $Gcorr = $Hcorr - $T*$Stot;
  my $Gcorr_quasiRRHO = $Hcorr - $T*$Stot_quasiRRHO;
  $Hcorr *= $kcal2hartree/1000;
  $Gcorr_quasiRRHO *= $kcal2hartree/1000;
  
  my $Grimme_G = $E_e + $Gcorr_quasiRRHO;
  return ($E_e, $enthalpy, $G, $Grimme_G);
} #end of sub get_thermo


#Aligns coords2 with coords1 to minimize RMSD for pairs of atoms from two lists
#my @new_coords = RMSD_align(\@coords, \@cat_coords, \@atoms1, \@atoms2);
#returns RMSD value and a new set of coords
sub RMSD_align {
  my ($coords1_ref, $coords2_ref, $atoms1_ref, $atoms2_ref) = @_;
  my @ref_coords = copy_coords($coords1_ref);
  my @compare_coords = copy_coords($coords2_ref);
  my @ref_atoms;
  my @compare_atoms;
  if($atoms1_ref && $atoms2_ref) {
    @ref_atoms = @{$atoms1_ref};
    @compare_atoms = @{$atoms2_ref};
  } else {
    #make lists of heavy atoms
    if($#ref_coords != $#compare_coords) {
      die "All-heavy-atom RMSD requires equal numbers of atoms!";
    } else {
      foreach my $atom (0..$#ref_coords) {
        if($ref_coords[$atom][0] ne "H") {
          push(@ref_atoms, $atom);
          push(@compare_atoms, $atom);
        }
      }
    }
  }

  if($#ref_atoms != $#compare_atoms) {
    print "Can't use RMSD align with different length atom lists!\n";
    print "ref_atoms = @ref_atoms\n";
    print "compare_atoms = @compare_atoms\n";
    return 0;
  }

  #shift coords to the COM of the specified three atoms for that structure
  my @centroid1 = (0,0,0,0,0);
  my @centroid2 = (0,0,0,0,0);
  foreach my $atom (0..$#ref_atoms) {
    foreach my $i (2..4) {
      $centroid1[$i] += $ref_coords[$ref_atoms[$atom]][$i];
      $centroid2[$i] += $compare_coords[$compare_atoms[$atom]][$i];
    }
  }
  foreach my $i (2..4) {
    $centroid1[$i] /= ($#ref_atoms+1);
    $centroid2[$i] /= ($#compare_atoms+1);
  }
  coord_shift(-$centroid1[2],-$centroid1[3],-$centroid1[4],\@ref_coords);
  coord_shift(-$centroid2[2],-$centroid2[3],-$centroid2[4],\@compare_coords);
  
  #Matrix to diagonalize
  my $matrix = new Math::MatrixReal(4,4);
  
  foreach my $atom (0..$#ref_atoms) {
    my $xm = $compare_coords[$compare_atoms[$atom]][2] - $ref_coords[$ref_atoms[$atom]][2];
    my $xp = $compare_coords[$compare_atoms[$atom]][2] + $ref_coords[$ref_atoms[$atom]][2];
    my $ym = $compare_coords[$compare_atoms[$atom]][3] - $ref_coords[$ref_atoms[$atom]][3];
    my $yp = $compare_coords[$compare_atoms[$atom]][3] + $ref_coords[$ref_atoms[$atom]][3];
    my $zm = $compare_coords[$compare_atoms[$atom]][4] - $ref_coords[$ref_atoms[$atom]][4];
    my $zp = $compare_coords[$compare_atoms[$atom]][4] + $ref_coords[$ref_atoms[$atom]][4];
  
    my $temp_matrix = Math::MatrixReal->new_from_rows( [[$xm*$xm + $ym*$ym + $zm*$zm, $yp*$zm - $ym*$zp,          $xm*$zp - $xp*$zm,           $xp*$ym - $xm*$yp],
                                                        [$yp*$zm - $ym*$zp,           $yp*$yp + $zp*$zp + $xm*$xm,$xm*$ym - $xp*$yp,           $xm*$zm - $xp*$zp],
                                                        [$xm*$zp - $xp*$zm,           $xm*$ym - $xp*$yp,          $xp*$xp + $zp*$zp + $ym*$ym, $ym*$zm - $yp*$zp], 
                                                        [$xp*$ym - $xm*$yp,           $xm*$zm - $xp*$zp,          $ym*$zm - $yp*$zp,           $xp*$xp + $yp*$yp + $zm*$zm]]);
  
    $matrix += $temp_matrix
  }
  
  my ($eigenvalues,$evectors) = $matrix->sym_diagonalize();
  
  #find smallest of four eigenvalues and save corresponding eigenvectors
  my $sd = 999;
  my $Q = new Math::MatrixReal(1,4);
  foreach my $i (1..4) {
    my $value = $eigenvalues->element($i,1);
    if($value < $sd) {
      $sd = $value;
      $Q = $evectors->column($i);
    }
  }
  
  my $rmsd = 0;
  if($sd > 0) { #to avoid very small negative numbers for sd (-1e-16, etc)
    $rmsd = sqrt($sd/($#ref_atoms+1));
  }

  quat_rot($Q->element(1,1), $Q->element(2,1), $Q->element(3,1), $Q->element(4,1), \@compare_coords, 0..$#compare_coords);
  #Now shift compare_coords back to original ref_coords COM
  coord_shift($centroid1[2],$centroid1[3],$centroid1[4],\@compare_coords);

  return ($rmsd, @compare_coords);
} #End RMSD_align


#finds minimum energy structure (based on LJ potential) by rotating list of target atoms around bond between $atom1 and $atom2
#TODO: if no list of @targets provided, rotate fragment that starts with atom1!
sub minimize_torsion {
  my ($atom1, $atom2, $coords_ref, @targets) = @_;
  my @coords = @{$coords_ref};
  my $increment = 5; #angle increment to search over

  my $E_min = 1E20;
  my $angle_min = 0;

  my @axis = ($coords[$atom2][2] - $coords[$atom1][2],$coords[$atom2][3] - $coords[$atom1][3],$coords[$atom2][4] - $coords[$atom1][4]);
  foreach my $count (0..360/$increment + 1) {
    my $angle = $count*$increment;
    unless (@targets) {
      my @connected = get_connected(\@coords);
      @targets = get_all_connected($atom1, $atom2, @connected);
    }
    center_genrotate($atom1, @axis, deg2rad($increment), \@coords, @targets);
    my $energy = LJ_energy(\@coords);
    if($energy < $E_min) {
      $angle_min = $angle;
      $E_min = $energy;
    }
  }
  
  center_genrotate($atom1, @axis, deg2rad($angle_min), \@coords, @targets);
}

#calculates  LJ-6-12 potential energy based on autodock Rij and Eij parameters
#simply ignores any atom pair involving elements for which parameters are missing (which shouldn't be anything!).
sub LJ_energy {
  my ($coords_ref) = @_;
  my @coords = @{$coords_ref};

  #Lennard-Jones Rij (sigma) and Eij (epsilon) parameters from Autodock4 (http://autodock.scripps.edu/local_files/scripts/lj4.py)
  my %Eij = ("CC" => 0.1500, "CN" => 0.1549, "NC" => 0.1549, "CO" => 0.1732, "OC" => 0.1732, "CP" => 0.1732, "PC" => 0.1732, "CS" => 0.1732,
    "SC" => 0.1732, "CH" => 0.0548, "HC" => 0.0548, "CFe" => 0.0387, "FeC" => 0.0387, "CF" => 0.1095, "FC" => 0.1095, "CCl" => 0.2035,
    "ClC" => 0.2035, "CBr" => 0.2416, "BrC" => 0.2416, "CI" => 0.2877, "IC" => 0.2877, "CMg" => 0.3623, "MgC" => 0.3623, "CZn" => 0.2872,
    "ZnC" => 0.2872, "CCa" => 0.2872, "CaC" => 0.2872, "NC" => 0.1549, "CN" => 0.1549, "NN" => 0.1600, "NO" => 0.1789, "ON" => 0.1789,
    "NP" => 0.1789, "PN" => 0.1789, "NS" => 0.1789, "SN" => 0.1789, "NH" => 0.0566, "HN" => 0.0566, "NFe" => 0.0400, "FeN" => 0.0400,
    "NF" => 0.1131, "FN" => 0.1131, "NCl" => 0.2101, "ClN" => 0.2101, "NBr" => 0.2495, "BrN" => 0.2495, "NI" => 0.2972, "IN" => 0.2972,
    "NMg" => 0.3742, "MgN" => 0.3742, "NZn" => 0.2966, "ZnN" => 0.2966, "NCa" => 0.2966, "CaN" => 0.2966, "OC" => 0.1732, "CO" => 0.1732,
    "ON" => 0.1789, "NO" => 0.1789, "OO" => 0.2000, "OP" => 0.2000, "PO" => 0.2000, "OS" => 0.2000, "SO" => 0.2000, "OH" => 0.0632,
    "HO" => 0.0632, "OFe" => 0.0447, "FeO" => 0.0447, "OF" => 0.1265, "FO" => 0.1265, "OCl" => 0.2349, "ClO" => 0.2349, "OBr" => 0.2789,
    "BrO" => 0.2789, "OI" => 0.3323, "IO" => 0.3323, "OMg" => 0.4183, "MgO" => 0.4183, "OZn" => 0.3317, "ZnO" => 0.3317, "OCa" => 0.3317,
    "CaO" => 0.3317, "PC" => 0.1732, "CP" => 0.1732, "PN" => 0.1789, "NP" => 0.1789, "PO" => 0.2000, "OP" => 0.2000, "PP" => 0.2000,
    "PS" => 0.2000, "SP" => 0.2000, "PH" => 0.0632, "HP" => 0.0632, "PFe" => 0.0447, "FeP" => 0.0447, "PF" => 0.1265, "FP" => 0.1265,
    "PCl" => 0.2349, "ClP" => 0.2349, "PBr" => 0.2789, "BrP" => 0.2789, "PI" => 0.3323, "IP" => 0.3323, "PMg" => 0.4183, "MgP" => 0.4183,
    "PZn" => 0.3317, "ZnP" => 0.3317, "PCa" => 0.3317, "CaP" => 0.3317, "SC" => 0.1732, "CS" => 0.1732, "SN" => 0.1789, "NS" => 0.1789,
    "SO" => 0.2000, "OS" => 0.2000, "SP" => 0.2000, "PS" => 0.2000, "SS" => 0.2000, "SH" => 0.0632, "HS" => 0.0632, "SFe" => 0.0447,
    "FeS" => 0.0447, "SF" => 0.1265, "FS" => 0.1265, "SCl" => 0.2349, "ClS" => 0.2349, "SBr" => 0.2789, "BrS" => 0.2789, "SI" => 0.3323,
    "IS" => 0.3323, "SMg" => 0.4183, "MgS" => 0.4183, "SZn" => 0.3317, "ZnS" => 0.3317, "SCa" => 0.3317, "CaS" => 0.3317, "HC" => 0.0548,
    "CH" => 0.0548, "HN" => 0.0566, "NH" => 0.0566, "HO" => 0.0632, "OH" => 0.0632, "HP" => 0.0632, "PH" => 0.0632, "HS" => 0.0632,
    "SH" => 0.0632, "HH" => 0.0200, "HFe" => 0.0141, "FeH" => 0.0141, "HF" => 0.0400, "FH" => 0.0400, "HCl" => 0.0743, "ClH" => 0.0743,
    "HBr" => 0.0882, "BrH" => 0.0882, "HI" => 0.1051, "IH" => 0.1051, "HMg" => 0.1323, "MgH" => 0.1323, "HZn" => 0.1049, "ZnH" => 0.1049,
    "HCa" => 0.1049, "CaH" => 0.1049, "FeC" => 0.0387, "CFe" => 0.0387, "FeN" => 0.0400, "NFe" => 0.0400, "FeO" => 0.0447, "OFe" => 0.0447,
    "FeP" => 0.0447, "PFe" => 0.0447, "FeS" => 0.0447, "SFe" => 0.0447, "FeH" => 0.0141, "HFe" => 0.0141, "FeFe" => 0.0100, "FeF" => 0.0283,
    "FFe" => 0.0283, "FeCl" => 0.0525, "ClFe" => 0.0525, "FeBr" => 0.0624, "BrFe" => 0.0624, "FeI" => 0.0743, "IFe" => 0.0743, "FeMg" => 0.0935,
    "MgFe" => 0.0935, "FeZn" => 0.0742, "ZnFe" => 0.0742, "FeCa" => 0.0742, "CaFe" => 0.0742, "FC" => 0.1095, "CF" => 0.1095, "FN" => 0.1131,
    "NF" => 0.1131, "FO" => 0.1265, "OF" => 0.1265, "FP" => 0.1265, "PF" => 0.1265, "FS" => 0.1265, "SF" => 0.1265, "FH" => 0.0400,
    "HF" => 0.0400, "FFe" => 0.0283, "FeF" => 0.0283, "FF" => 0.0800, "FCl" => 0.1486, "ClF" => 0.1486, "FBr" => 0.1764, "BrF" => 0.1764,
    "FI" => 0.2101, "IF" => 0.2101, "FMg" => 0.2646, "MgF" => 0.2646, "FZn" => 0.2098, "ZnF" => 0.2098, "FCa" => 0.2098, "CaF" => 0.2098,
    "ClC" => 0.2035, "CCl" => 0.2035, "ClN" => 0.2101, "NCl" => 0.2101, "ClO" => 0.2349, "OCl" => 0.2349, "ClP" => 0.2349, "PCl" => 0.2349,
    "ClS" => 0.2349, "SCl" => 0.2349, "ClH" => 0.0743, "HCl" => 0.0743, "ClFe" => 0.0525, "FeCl" => 0.0525, "ClF" => 0.1486, "FCl" => 0.1486,
    "ClCl" => 0.2760, "ClBr" => 0.3277, "BrCl" => 0.3277, "ClI" => 0.3903, "ICl" => 0.3903, "ClMg" => 0.4914, "MgCl" => 0.4914, "ClZn" => 0.3896,
    "ZnCl" => 0.3896, "ClCa" => 0.3896, "CaCl" => 0.3896, "BrC" => 0.2416, "CBr" => 0.2416, "BrN" => 0.2495, "NBr" => 0.2495, "BrO" => 0.2789,
    "OBr" => 0.2789, "BrP" => 0.2789, "PBr" => 0.2789, "BrS" => 0.2789, "SBr" => 0.2789, "BrH" => 0.0882, "HBr" => 0.0882, "BrFe" => 0.0624,
    "FeBr" => 0.0624, "BrF" => 0.1764, "FBr" => 0.1764, "BrCl" => 0.3277, "ClBr" => 0.3277, "BrBr" => 0.3890, "BrI" => 0.4634, "IBr" => 0.4634,
    "BrMg" => 0.5834, "MgBr" => 0.5834, "BrZn" => 0.4625, "ZnBr" => 0.4625, "BrCa" => 0.4625, "CaBr" => 0.4625, "IC" => 0.2877, "CI" => 0.2877,
    "IN" => 0.2972, "NI" => 0.2972, "IO" => 0.3323, "OI" => 0.3323, "IP" => 0.3323, "PI" => 0.3323, "IS" => 0.3323, "SI" => 0.3323,
    "IH" => 0.1051, "HI" => 0.1051, "IFe" => 0.0743, "FeI" => 0.0743, "IF" => 0.2101, "FI" => 0.2101, "ICl" => 0.3903, "ClI" => 0.3903,
    "IBr" => 0.4634, "BrI" => 0.4634, "II" => 0.5520, "IMg" => 0.6950, "MgI" => 0.6950, "IZn" => 0.5510, "ZnI" => 0.5510, "ICa" => 0.5510,
    "CaI" => 0.5510, "MgC" => 0.3623, "CMg" => 0.3623, "MgN" => 0.3742, "NMg" => 0.3742, "MgO" => 0.4183, "OMg" => 0.4183, "MgP" => 0.4183,
    "PMg" => 0.4183, "MgS" => 0.4183, "SMg" => 0.4183, "MgH" => 0.1323, "HMg" => 0.1323, "MgFe" => 0.0935, "FeMg" => 0.0935, "MgF" => 0.2646,
    "FMg" => 0.2646, "MgCl" => 0.4914, "ClMg" => 0.4914, "MgBr" => 0.5834, "BrMg" => 0.5834, "MgI" => 0.6950, "IMg" => 0.6950, "MgMg" => 0.8750,
    "MgZn" => 0.6937, "ZnMg" => 0.6937, "MgCa" => 0.6937, "CaMg" => 0.6937, "ZnC" => 0.2872, "CZn" => 0.2872, "ZnN" => 0.2966, "NZn" => 0.2966,
    "ZnO" => 0.3317, "OZn" => 0.3317, "ZnP" => 0.3317, "PZn" => 0.3317, "ZnS" => 0.3317, "SZn" => 0.3317, "ZnH" => 0.1049, "HZn" => 0.1049,
    "ZnFe" => 0.0742, "FeZn" => 0.0742, "ZnF" => 0.2098, "FZn" => 0.2098, "ZnCl" => 0.3896, "ClZn" => 0.3896, "ZnBr" => 0.4625, "BrZn" => 0.4625,
    "ZnI" => 0.5510, "IZn" => 0.5510, "ZnMg" => 0.6937, "MgZn" => 0.6937, "ZnZn" => 0.5500, "ZnCa" => 0.5500, "CaZn" => 0.5500, "CaC" => 0.2872,
    "CCa" => 0.2872, "CaN" => 0.2966, "NCa" => 0.2966, "CaO" => 0.3317, "OCa" => 0.3317, "CaP" => 0.3317, "PCa" => 0.3317, "CaS" => 0.3317,
    "SCa" => 0.3317, "CaH" => 0.1049, "HCa" => 0.1049, "CaFe" => 0.0742, "FeCa" => 0.0742, "CaF" => 0.2098, "FCa" => 0.2098, "CaCl" => 0.3896,
    "ClCa" => 0.3896, "CaBr" => 0.4625, "BrCa" => 0.4625, "CaI" => 0.5510, "ICa" => 0.5510, "CaMg" => 0.6937, "MgCa" => 0.6937, "CaZn" => 0.5500,
    "ZnCa" => 0.5500, "CaCa" => 0.5500, "X" => 0);
  
  my %Rij = ("CC" => 4.00, "CN" => 3.75, "NC" => 3.75, "CO" => 3.60, "OC" => 3.60, "CP" => 4.10, "PC" => 4.10, "CS" => 4.00,
    "SC" => 4.00, "CH" => 3.00, "HC" => 3.00, "CFe" => 2.65, "FeC" => 2.65, "CF" => 3.54, "FC" => 3.54, "CCl" => 4.04,
    "ClC" => 4.04, "CBr" => 4.17, "BrC" => 4.17, "CI" => 4.36, "IC" => 4.36, "CMg" => 2.65, "MgC" => 2.65, "CZn" => 2.74,
    "ZnC" => 2.74, "CCa" => 2.99, "CaC" => 2.99, "NC" => 3.75, "CN" => 3.75, "NN" => 3.50, "NO" => 3.35, "ON" => 3.35,
    "NP" => 3.85, "PN" => 3.85, "NS" => 3.75, "SN" => 3.75, "NH" => 2.75, "HN" => 2.75, "NFe" => 2.40, "FeN" => 2.40,
    "NF" => 3.29, "FN" => 3.29, "NCl" => 3.79, "ClN" => 3.79, "NBr" => 3.92, "BrN" => 3.92, "NI" => 4.11, "IN" => 4.11,
    "NMg" => 2.40, "MgN" => 2.40, "NZn" => 2.49, "ZnN" => 2.49, "NCa" => 2.74, "CaN" => 2.74, "OC" => 3.60, "CO" => 3.60,
    "ON" => 3.35, "NO" => 3.35, "OO" => 3.20, "OP" => 3.70, "PO" => 3.70, "OS" => 3.60, "SO" => 3.60, "OH" => 2.60,
    "HO" => 2.60, "OFe" => 2.25, "FeO" => 2.25, "OF" => 3.15, "FO" => 3.15, "OCl" => 3.65, "ClO" => 3.65, "OBr" => 3.77,    "BrO" => 3.77, "OI" => 3.96, "IO" => 3.96, "OMg" => 2.25, "MgO" => 2.25, "OZn" => 2.34, "ZnO" => 2.34, "OCa" => 2.59,
    "CaO" => 2.59, "PC" => 4.10, "CP" => 4.10, "PN" => 3.85, "NP" => 3.85, "PO" => 3.70, "OP" => 3.70, "PP" => 4.20,
    "PS" => 4.10, "SP" => 4.10, "PH" => 3.10, "HP" => 3.10, "PFe" => 2.75, "FeP" => 2.75, "PF" => 3.65, "FP" => 3.65,
    "PCl" => 4.14, "ClP" => 4.14, "PBr" => 4.27, "BrP" => 4.27, "PI" => 4.46, "IP" => 4.46, "PMg" => 2.75, "MgP" => 2.75,
    "PZn" => 2.84, "ZnP" => 2.84, "PCa" => 3.09, "CaP" => 3.09, "SC" => 4.00, "CS" => 4.00, "SN" => 3.75, "NS" => 3.75,
    "SO" => 3.60, "OS" => 3.60, "SP" => 4.10, "PS" => 4.10, "SS" => 4.00, "SH" => 3.00, "HS" => 3.00, "SFe" => 2.65,
    "FeS" => 2.65, "SF" => 3.54, "FS" => 3.54, "SCl" => 4.04, "ClS" => 4.04, "SBr" => 4.17, "BrS" => 4.17, "SI" => 4.36,
    "IS" => 4.36, "SMg" => 2.65, "MgS" => 2.65, "SZn" => 2.74, "ZnS" => 2.74, "SCa" => 2.99, "CaS" => 2.99, "HC" => 3.00,
    "CH" => 3.00, "HN" => 2.75, "NH" => 2.75, "HO" => 2.60, "OH" => 2.60, "HP" => 3.10, "PH" => 3.10, "HS" => 3.00,
    "SH" => 3.00, "HH" => 2.00, "HFe" => 1.65, "FeH" => 1.65, "HF" => 2.54, "FH" => 2.54, "HCl" => 3.04, "ClH" => 3.04,
    "HBr" => 3.17, "BrH" => 3.17, "HI" => 3.36, "IH" => 3.36, "HMg" => 1.65, "MgH" => 1.65, "HZn" => 1.74, "ZnH" => 1.74,
    "HCa" => 1.99, "CaH" => 1.99, "FeC" => 2.65, "CFe" => 2.65, "FeN" => 2.40, "NFe" => 2.40, "FeO" => 2.25, "OFe" => 2.25,
    "FeP" => 2.75, "PFe" => 2.75, "FeS" => 2.65, "SFe" => 2.65, "FeH" => 1.65, "HFe" => 1.65, "FeFe" => 1.30, "FeF" => 2.19,
    "FFe" => 2.19, "FeCl" => 2.69, "ClFe" => 2.69, "FeBr" => 2.81, "BrFe" => 2.81, "FeI" => 3.01, "IFe" => 3.01, "FeMg" => 1.30,
    "MgFe" => 1.30, "FeZn" => 1.39, "ZnFe" => 1.39, "FeCa" => 1.64, "CaFe" => 1.64, "FC" => 3.54, "CF" => 3.54, "FN" => 3.29,
    "NF" => 3.29, "FO" => 3.15, "OF" => 3.15, "FP" => 3.65, "PF" => 3.65, "FS" => 3.54, "SF" => 3.54, "FH" => 2.54,
    "HF" => 2.54, "FFe" => 2.19, "FeF" => 2.19, "FF" => 3.09, "FCl" => 3.59, "ClF" => 3.59, "FBr" => 3.71, "BrF" => 3.71,
    "FI" => 3.90, "IF" => 3.90, "FMg" => 2.19, "MgF" => 2.19, "FZn" => 2.29, "ZnF" => 2.29, "FCa" => 2.54, "CaF" => 2.54,
    "ClC" => 4.04, "CCl" => 4.04, "ClN" => 3.79, "NCl" => 3.79, "ClO" => 3.65, "OCl" => 3.65, "ClP" => 4.14, "PCl" => 4.14,
    "ClS" => 4.04, "SCl" => 4.04, "ClH" => 3.04, "HCl" => 3.04, "ClFe" => 2.69, "FeCl" => 2.69, "ClF" => 3.59, "FCl" => 3.59,
    "ClCl" => 4.09, "ClBr" => 4.21, "BrCl" => 4.21, "ClI" => 4.40, "ICl" => 4.40, "ClMg" => 2.69, "MgCl" => 2.69, "ClZn" => 2.79,
    "ZnCl" => 2.79, "ClCa" => 3.04, "CaCl" => 3.04, "BrC" => 4.17, "CBr" => 4.17, "BrN" => 3.92, "NBr" => 3.92, "BrO" => 3.77,
    "OBr" => 3.77, "BrP" => 4.27, "PBr" => 4.27, "BrS" => 4.17, "SBr" => 4.17, "BrH" => 3.17, "HBr" => 3.17, "BrFe" => 2.81,
    "FeBr" => 2.81, "BrF" => 3.71, "FBr" => 3.71, "BrCl" => 4.21, "ClBr" => 4.21, "BrBr" => 4.33, "BrI" => 4.53, "IBr" => 4.53,
    "BrMg" => 2.81, "MgBr" => 2.81, "BrZn" => 2.91, "ZnBr" => 2.91, "BrCa" => 3.16, "CaBr" => 3.16, "IC" => 4.36, "CI" => 4.36,
    "IN" => 4.11, "NI" => 4.11, "IO" => 3.96, "OI" => 3.96, "IP" => 4.46, "PI" => 4.46, "IS" => 4.36, "SI" => 4.36,
    "IH" => 3.36, "HI" => 3.36, "IFe" => 3.01, "FeI" => 3.01, "IF" => 3.90, "FI" => 3.90, "ICl" => 4.40, "ClI" => 4.40,
    "IBr" => 4.53, "BrI" => 4.53, "II" => 4.72, "IMg" => 3.01, "MgI" => 3.01, "IZn" => 3.10, "ZnI" => 3.10, "ICa" => 3.35,
    "CaI" => 3.35, "MgC" => 2.65, "CMg" => 2.65, "MgN" => 2.40, "NMg" => 2.40, "MgO" => 2.25, "OMg" => 2.25, "MgP" => 2.75,
    "PMg" => 2.75, "MgS" => 2.65, "SMg" => 2.65, "MgH" => 1.65, "HMg" => 1.65, "MgFe" => 1.30, "FeMg" => 1.30, "MgF" => 2.19,
    "FMg" => 2.19, "MgCl" => 2.69, "ClMg" => 2.69, "MgBr" => 2.81, "BrMg" => 2.81, "MgI" => 3.01, "IMg" => 3.01, "MgMg" => 1.30,
    "MgZn" => 1.39, "ZnMg" => 1.39, "MgCa" => 1.64, "CaMg" => 1.64, "ZnC" => 2.74, "CZn" => 2.74, "ZnN" => 2.49, "NZn" => 2.49,
    "ZnO" => 2.34, "OZn" => 2.34, "ZnP" => 2.84, "PZn" => 2.84, "ZnS" => 2.74, "SZn" => 2.74, "ZnH" => 1.74, "HZn" => 1.74,
    "ZnFe" => 1.39, "FeZn" => 1.39, "ZnF" => 2.29, "FZn" => 2.29, "ZnCl" => 2.79, "ClZn" => 2.79, "ZnBr" => 2.91, "BrZn" => 2.91,
    "ZnI" => 3.10, "IZn" => 3.10, "ZnMg" => 1.39, "MgZn" => 1.39, "ZnZn" => 1.48, "ZnCa" => 1.73, "CaZn" => 1.73, "CaC" => 2.99,
    "CCa" => 2.99, "CaN" => 2.74, "NCa" => 2.74, "CaO" => 2.59, "OCa" => 2.59, "CaP" => 3.09, "PCa" => 3.09, "CaS" => 2.99,
    "SCa" => 2.99, "CaH" => 1.99, "HCa" => 1.99, "CaFe" => 1.64, "FeCa" => 1.64, "CaF" => 2.54, "FCa" => 2.54, "CaCl" => 3.04,
    "ClCa" => 3.04, "CaBr" => 3.16, "BrCa" => 3.16, "CaI" => 3.35, "ICa" => 3.35, "CaMg" => 1.64, "MgCa" => 1.64, "CaZn" => 1.73,
    "ZnCa" => 1.73, "CaCa" => 1.98, "X" => 0);

  my $energy = 0;

  foreach my $atom1 (0..$#coords - 1) {
    foreach my $atom2 ($atom1+1..$#coords) {
      my $string = $coords[$atom1][0] . $coords[$atom2][0];
      if(my $sigma = $Rij{$string}) {
        my $epsilon = $Eij{$string};
        my $R = distance($atom1, $atom2, \@coords);
        $energy += $epsilon*(($sigma/$R)**12 - ($sigma/$R)**6);
      }
    }
  }
  return $energy;
}

#Given references to arrays holding (R) and (S) absolute energies, compute ee based on Boltzmann weighted sums
#my $ee = calculate_ee(\@R_energies, \@S_energies, $temp);
sub calculate_ee {
  my ($R_ref, $S_ref, $temp) = @_;

  my $RT = $boltzmann*$temp;
  my @R_vals = @{$R_ref};
  my @S_vals = @{$S_ref};
  my $R_sum = 0;
  my $S_sum = 0;
  foreach (@R_vals) {
    $R_sum += exp(-$hart2kcal*($_ - $R_vals[0])/$RT);
  }
  foreach (@S_vals) {
    $S_sum += exp(-$hart2kcal*($_ - $R_vals[0])/$RT);
  }
  my $ee = 0;
  if($R_sum + $S_sum != 0) {
    $ee = ($R_sum - $S_sum)/($R_sum + $S_sum);
  }
  return $ee;
} #end calculate_ee

#Prints random string
#random_string(16, 'a'..'z') will be a random string of 16 letters from a to z
sub random_string { 
  join'', @_[ map{ rand @_ } 1 .. shift ] 
}

#takes reference to a @coords array and converts it to a scalar with escaped line returns (for ALEXANDER mostly)
#my $geom = flatten(\@coords);
sub flatten {
  my ($coords_ref) = @_;
  my @coords = @{$coords_ref};

  my $num_atoms = $#coords + 1;
  my $geometry = "$num_atoms\\\\n\\\\n";
  foreach my $atom (0..$#coords) {
    $geometry = $geometry . "$coords[$atom][0]  $coords[$atom][2]    $coords[$atom][3]   $coords[$atom][4]\\\\n";
  }

  return $geometry;
} #end flatten


#
#@key_atoms = [1,2] #central atoms
#@other_atom = [4 1 3]
#@other_atom = [5 3 4]
#overlay to similar molecules together by rotating the lingking bond of key fragment of molecules. 
#my @coords = match_molecules(\@coords1, \@coords2, \@keyatoms1, \@keyatoms2, \@bonds1, \@bonds2)
#the @bonds1 and @bonds2 are array of arrays
#@bonds1 = (
#          [1,3],
#          [3,4]
#          [3,4]
#          )
#@keyatoms = (
#            [1,2,3],
#            [4],
#            [5],
#            [5],
#            )
# Since we need 3 keyatoms for the central part and one key atom for each fragment,the order of the peripherial fragment should match that of bonds
#note: that the atom number and the linking order of bonds must be exatly the same
#note: this is the basic version, it can only used for Yanfei's catalyst, or other similar system with three rigid fragment in the molecule.
sub match_molecules {
   my ($coords1, $coords2, $keyatoms1, $keyatoms2, $bonds2) = @_;
   my @coords1 = @{$coords1};
   my @coords2 = @{$coords2};
   my @keyatoms1 = @{$keyatoms1};
   my @keyatoms2 = @{$keyatoms2};
   my @RMSDatoms1 = @{$keyatoms1[0]};
   my @RMSDatoms2 = @{$keyatoms2[0]};
   my @bonds2 = @{$bonds2};
   my @connected1 = &get_connected($coords1);
   my @connected2 = &get_connected($coords2);
   (my $RMSD_temp1, @coords2) = &RMSD_align(\@coords1, \@coords2, \@RMSDatoms1, \@RMSDatoms2);
   my @coords2_final = @coords2;
   for my $i (0..$#bonds2) {
      my $r_atom1 = $bonds2[$i][0];
      my $r_atom2 = $bonds2[$i][1];
      my @fragment2 = &get_all_connected($r_atom1, $r_atom2, \@connected2);
      my @keyatom1 = @{$keyatoms1[$i+1]};
      my @keyatom2 = @{$keyatoms2[$i+1]};
      push(@RMSDatoms1, @keyatom1);
      push(@RMSDatoms2, @keyatom2);
      my $increment = 2;
      my $minRMSD = 9999;
      for my $n(0..360/$increment) { 
         my @axis = ($coords2[$r_atom1][2]-$coords2[$r_atom2][2], $coords2[$r_atom1][3]-$coords2[$r_atom2][3], $coords2[$r_atom1][4]-$coords2[$r_atom2][4]);
         &center_genrotate($r_atom1, @axis, deg2rad($increment), \@coords2, @fragment2);
         my ($RMSD_temp, @coords_temp) = &RMSD_align(\@coords1, \@coords2, \@RMSDatoms1, \@RMSDatoms2);
         if ($RMSD_temp < $minRMSD) {
           $minRMSD = $RMSD_temp;
           @coords2_final = @coords_temp;
         }
      }
   }
   return @coords2_final;
}  

sub sub_geo {
 my ($sub, $prop) = @_;
 return $sub_geo->{$sub}->{$prop};
}


sub map_catalyst {
    my ($coords, $cat_coords, $ts_key_atoms, $cata_key_atoms, $bonds_RMSD, $bonds_LJ, $first_cat_atom) = @_;
    print "New catalyst is being mapped to the template, this may take several minutes...\n";

    my @new_cata = &match_molecules($coords, $cat_coords, $ts_key_atoms, $cata_key_atoms, $bonds_RMSD);
    #remove old cat and then add new cat
    &remove_atoms($coords, $first_cat_atom..$#{ $coords });
    my @coords = &combine_coords($coords, \@new_cata);
    #if there are only two keyatoms, this match can not give the right position, so rotate the 
    #catalyst around the bonds between key atoms to minimize the LJ potential
    if ($#{ $cata_key_atoms->[0] } == 1) {
        my @bond = map {$_ + $first_cat_atom} @{ $cata_key_atoms->[0] };
        minimize_torsion(@bond, \@coords, $first_cat_atom..$#coords);
    }

    #now rotate the peripheral bonds to minimize the LJ potential of the system
    #First we need to remark the bonds_LJ since the catalyst has been added to the ts and the
    #number of atoms has changed
    #rotate the bonds in @bonds_LJnew to minimze the LJ potential
    for (@{ $bonds_LJ }) {
        #change label for each bond
        my @bond = map {$_ + $first_cat_atom} @{ $_ };
        
        my @connected = get_connected(\@coords);
        #filt out catalyst atoms
        my @target = grep { $_ >= $first_cat_atom} get_all_connected(@bond, \@connected);
        minimize_torsion(@bond, \@coords, @target);
    }
    return @coords;
} #End Map_catalyst

          
   
sub fix_coords {
  my ($coords, $subs, $specified) = @_;
  my @coords = @{$coords};
  my @subs = @{$subs};
  my @specified = @{$specified};
  my @connect = get_connected(\@coords);
  my @molecules = get_molecules(@connect);
  my @fix;
  for (@molecules) {
    my @molecule = @{$_};
    my $check;
    for my $i (@subs) { 
      if (grep { $i == $_ } @molecule) {
        $check = 1;
      }
    }
    if(!$check) {
      push(@fix, @molecule);
    } 
  }
  push(@fix, @specified);
  return @fix;
} 

                
sub get_molecules {
  my @connected = @_;
  my @molecules;
  my @molecules_count;
#  #remove avoid_atom from list of atoms connected to $start_atom
  my $start_atom = 0;
  for(;;) {
    my $i=0;
    my @positions = ($start_atom);
    #I can probably use something simpler than a hash here (since I'm not really using the values)
    my %visited = ( $start_atom => '0') ; #keys are numbers of the visited atoms, values are not used

    #loop until @positions is empty
    while(@positions) {
      my $position = shift(@positions);
      foreach my $atom (@{$connected[$position]}) {       #grab all atoms connected to current atom and add to queue (unless already visited)
        if($atom >= 0 && !exists $visited{$atom}) {
          push(@positions, $atom);
          $visited{$atom} = 0;
        }
      }
    }
    my @molecule = keys %visited;
    push (@molecules, \@molecule);
    push (@molecules_count, @molecule);
    while (grep { $i == $_ } @molecules_count) {
         $i++;
    }
    $start_atom= $i;
    $start_atom == $#connected+1 && last;
  }
  return (@molecules);
}


#return the RMSD of two coords
sub RMSD {
  my ($ref_coords_ref, $compare_coords_ref) = @_; 
  my @ref_coords = @$ref_coords_ref;
  my @compare_coords = @$compare_coords_ref;


  my $num_common_atoms;
  if($#ref_coords < $#compare_coords) {
  	$num_common_atoms = $#ref_coords;
  } else {
  	$num_common_atoms = $#compare_coords;
  }
  
  #shift coords to centroid
  my @centroid1 = (0,0,0,0,0);
  my @centroid2 = (0,0,0,0,0);
  foreach my $atom (0..$num_common_atoms) {
  	foreach my $i (2..4) {
  		$centroid1[$i] += $ref_coords[$atom][$i];
  		$centroid2[$i] += $compare_coords[$atom][$i];
  	}
  }
  foreach my $i (2..4) {
  	$centroid1[$i] /= ($num_common_atoms+1);
  	$centroid2[$i] /= ($num_common_atoms+1);
  }
  coord_shift(-$centroid1[2],-$centroid1[3],-$centroid1[4],\@ref_coords);
  coord_shift(-$centroid2[2],-$centroid2[3],-$centroid2[4],\@compare_coords);
  
  my @sd;
  
  #Run regular RMSD comparison as well as mirroring compare_coords across each plane
  #Report minimum value as the RMSD
  foreach my $reorder (0,1) {
    foreach my $mirrorx (0,1) {
      foreach my $mirrory (0,1) {
        foreach my $mirrorz (0,1) {
         if($reorder==0 && ($mirrory || $mirrorz)) {	#always re-order, except for non-mirrored and x-mirrored
          last;
         }
          if($mirrorx) {
            foreach my $atom (0..$num_common_atoms) {
              $compare_coords[$atom][2] = -$compare_coords[$atom][2];
            }
          }
          if($mirrory) {
            foreach my $atom (0..$num_common_atoms) {
              $compare_coords[$atom][3] = -$compare_coords[$atom][3];
            }
          }
          if($mirrorz) {
            foreach my $atom (0..$num_common_atoms) {
              $compare_coords[$atom][4] = -$compare_coords[$atom][4];
            }
          }
  
          if($reorder) {
            my @temp_compare_coords;
            foreach my $atom (0..$num_common_atoms) {	#for each atom in ref_coords, find closest atom of same element in compare_coords and put into temp_compare_coords
              my $short_distance = 999;
              my $short_atom = 0;
              foreach my $atom2 (0..$#compare_coords) {
                if($compare_coords[$atom2][0] eq $ref_coords[$atom][0]) {
                  my $distance;
                  foreach my $i (2..4) {
                    $distance += ($compare_coords[$atom2][$i] - $ref_coords[$atom][$i])**2;
                  }
                  $distance = sqrt($distance);
                  if($distance < $short_distance) {
                    $short_distance = $distance;
                    $short_atom = $atom2;
                  }
                }
              }
  	#push closest atom onto @temp_compare_coords and remove from @compare_coords
              push(@temp_compare_coords, splice(@compare_coords, $short_atom, 1));
            }
            @compare_coords = @temp_compare_coords;
          }
  
          #Matrix to diagonalize
          my $matrix = new Math::MatrixReal(4,4);
          
          foreach my $atom (0..$num_common_atoms) {
            if($compare_coords[$atom][0] ne $ref_coords[$atom][0]) {
              print "Something terribly wrong with atom ordering! (elements do not match)\n";
            }
        
            my $xm = $compare_coords[$atom][2] - $ref_coords[$atom][2];
            my $xp = $compare_coords[$atom][2] + $ref_coords[$atom][2];
            my $ym = $compare_coords[$atom][3] - $ref_coords[$atom][3];
            my $yp = $compare_coords[$atom][3] + $ref_coords[$atom][3];
            my $zm = $compare_coords[$atom][4] - $ref_coords[$atom][4];
            my $zp = $compare_coords[$atom][4] + $ref_coords[$atom][4];
          
            my $temp_matrix = Math::MatrixReal->new_from_rows( [[$xm*$xm + $ym*$ym + $zm*$zm, $yp*$zm - $ym*$zp,          $xm*$zp - $xp*$zm,           $xp*$ym - $xm*$yp],
                                                                [$yp*$zm - $ym*$zp,           $yp*$yp + $zp*$zp + $xm*$xm,$xm*$ym - $xp*$yp,           $xm*$zm - $xp*$zp],
                                                                [$xm*$zp - $xp*$zm,           $xm*$ym - $xp*$yp,          $xp*$xp + $zp*$zp + $ym*$ym, $ym*$zm - $yp*$zp], 
                                                                [$xp*$ym - $xm*$yp,           $xm*$zm - $xp*$zp,          $ym*$zm - $yp*$zp,           $xp*$xp + $yp*$yp + $zm*$zm]]);
          
            $matrix += $temp_matrix
          }
          
          my $eigenvalues = $matrix->sym_eigenvalues();
      	#Now unmirror!
          if($mirrorx) {
            foreach my $atom (0..$num_common_atoms) {
              $compare_coords[$atom][2] = -$compare_coords[$atom][2];
            }
          }
          if($mirrory) {
            foreach my $atom (0..$num_common_atoms) {
              $compare_coords[$atom][3] = -$compare_coords[$atom][3];
            }
          }
          if($mirrorz) {
            foreach my $atom (0..$num_common_atoms) {
              $compare_coords[$atom][4] = -$compare_coords[$atom][4];
            }
          }
          
          #find smallest of four eigenvalues
          my $sd = 999;
          foreach my $i (1..4) {
          	my $value = $eigenvalues->element($i,1);
          	if($value < $sd) {
          		$sd = $value;
          	}
          }
          push(@sd,$sd);	#compute RMSD
        }
      }
    }
  }
  
  @sd = sort {$a <=> $b} @sd;
  if($sd[0] < 0) {
    return sqrt(-$sd[0]/($num_common_atoms+1));
  } else {
    return sqrt($sd[0]/($num_common_atoms+1));
  }
} #End sub RMSD

#reads log file with IRC and and returns references to three arrays (@irc, @energies, @reaction_path)
#@irc is an array of @coords-like 2-D arrays
sub read_IRC {
  my $filename = $_[0];
  my @irc;
  my @energies;
  my @reaction_path;

  if(-e $filename) {                                   #check to make sure file exists first
    open (INFILE, "<$filename") or die "Can't open $filename";
    my @elements=('Bq','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe');

    if($filename !~ /\S+\.log/) { 
      print "Expecting .log file!\n";
      return 0;
    }
    #grab each geometry from log file
    while (<INFILE>) {
      my $line=$_;
      if($line =~ / orientation:/) {
        my @coords = ();
        $line = <INFILE>;
        $line = <INFILE>;
        $line = <INFILE>;
        $line = <INFILE>;
        $line = <INFILE>;
        do {
          if($line =~ /^\s+\d+\s+(\S+)\s+\S*\s+(\S+)\s+(\S+)\s+(\S+)/) {
            my @atom = ($elements[$1], 0, $2, $3, $4);
            push(@coords, \@atom);
            $line = <INFILE>;
          }
        } while(!($line =~ /--/));
        push(@irc, \@coords);
      } elsif ($line =~ /Summary of reaction path following/) {
        $line = <INFILE>;
        $line = <INFILE>;
        $line = <INFILE>;
        do {
          if($line =~ / \d+\s+(\S+)\s+(\S+)/) {
            push(@energies, $1*$hart2kcal);
            push(@reaction_path, $2/$ang2bohr);
            $line = <INFILE>;
          }
        } while(!($line =~ /--/));
      }
    }
  }
  return (\@irc, \@energies, \@reaction_path);
} #end read_irc


#sets a given dihedral angle
#set_dihedral(atom1, atom2, atom3, atom4, new_angle (in degrees!), \@coords)
sub set_dihedral {
  my ($atom1, $atom2, $atom3, $atom4, $new_tau, $coords_ref) = @_;
  my @coords = @{$coords_ref};

  my $tau = dihedral($atom1, $atom2, $atom3, $atom4, \@coords);
  change_dihedral($atom2, $atom3, deg2rad($new_tau - $tau), \@coords);
}


#Builds carbon nanotubes or nanotube fragment of a given size and radius of curvature
#For complete (N,N)-CNT: @coords = nt_builder($N, $length)
#For fragment: @coords = nt_builder($width, $length, $radius)
#or
#@coords = nt_builder($width, $length, $radius, $angular_offset) where angular_offset is number of rings rotated about z-axis
sub nt_builder {
  my ($width, $length, $radius, $angular_offset) = @_;
  #Standard C-C bond length and C-H bond length
  my $CC = 1.415;
  my $CH = 1.08;
  
  my $fragment = 0;
  my $angle_tally;
  my @coords;
  if(not defined($angular_offset)) {
    $angular_offset = 0;
  }
  
  if(defined($radius)) {
    $fragment = 1;
  } else {
    if($width < 2) {
      die("Can't build nanotube smaller than (2,2)!\n");
    } else {
      $radius = newton($CC, $width);  #Quick Newton-Raphson solver to get nanotube radius numerically
    }
  }

  
  my $CC_angle = 2*asin($CC/(2*$radius));
  my $CC_halfangle = 2*asin($CC/(4*$radius));
  my $CC_side = $CC*sqrt(3.0)/2.0;
  
  my @atom = ("C", 0, $radius, 0, 0);
  my @atom_coords;
  push(@atom_coords, \@atom);
  
  #Offset rotation to place center of fragment normal to x-axis
  rotate('z', -($width/2+$width-1)*$CC_angle - $angular_offset*2*($CC_angle+$CC_halfangle), \@atom_coords);
  #slide to center tube/fragment along z-axis
  #slide($CC_side*($length-1)/2);
  coord_shift(0, 0, $CC_side*($length-1)/2, \@atom_coords);
  $angle_tally = 0;
  
  for(my $row=0; $row<$length; $row++) {
    if($row%2==1) {
      if($row!=$length-1 || $fragment==0) {
        @coords = combine_coords(\@coords, \@atom_coords);
      }
      rotate('z', $CC_angle+2*$CC_halfangle, \@atom_coords);
      $angle_tally += $CC_angle+2*$CC_halfangle;
    } else {
      rotate('z', $CC_halfangle, \@atom_coords);
      $angle_tally += $CC_halfangle;
    }
    for (my $ring=0; $ring<$width; $ring++) {
      if($row!=$length-1 || $row%2!=1 || $ring != $width-1 || $fragment==0) {
        @coords = combine_coords(\@coords, \@atom_coords);
      }
      rotate('z', $CC_angle, \@atom_coords);
      $angle_tally += $CC_angle;
      if($row%2!=1 || $ring != $width-1) {
        @coords = combine_coords(\@coords, \@atom_coords);
      }
      rotate('z', $CC_angle+2*$CC_halfangle, \@atom_coords);
      $angle_tally += $CC_angle+2*$CC_halfangle;
    }
  
    #Reset and shift
  #  slide(-$CC_side);
    coord_shift(0, 0, -$CC_side, \@atom_coords); 
    rotate('z', -$angle_tally, \@atom_coords);
    $angle_tally = 0;
  }
  
  #Cap open valences
  my $numCs = $#coords;
  foreach my $atom1 (0..$numCs) {
    my @vector;
    my $neighbors = 0;
    foreach my $atom2 (0..$numCs) {
      if($atom1 != $atom2) {
        if(distance($atom1, $atom2, \@coords) < $CC+0.1) {
          $neighbors++;
          $vector[0] += $coords[$atom2][2]-$coords[$atom1][2];
          $vector[1] += $coords[$atom2][3]-$coords[$atom1][3];
          $vector[2] += $coords[$atom2][4]-$coords[$atom1][4];
        }
      }
    }
    if($neighbors < 3) {
      my $norm = sqrt($vector[0]**2+$vector[1]**2+$vector[2]**2);
      foreach (@vector) {
        $_ *= $CH/$norm
      }
    
      my @Hatom = ("H", 0, $coords[$atom1][2] - $vector[0], $coords[$atom1][3] - $vector[1], $coords[$atom1][4] - $vector[2]);
      push(@coords, [@Hatom]);
    }
    if($neighbors < 2) {
      die "Dangling carbon $atom1!";
    }
    if($neighbors > 4) {
      die "Too many neighbors for atom $atom1 (radius too small to accommodate $width rings)";
    }
  }

  return @coords;
}
  
  
#Dirty Newton solver to get radius of closed CNTs
sub newton {
  my ($CC, $width) = @_;
  #Threshold for Newton-Raphson solver to get radius
  my $THRESH = 1E-10;

  my $lastradius = 3*$CC*$width/(2*pi); #starting guess from conventional formula
  my $old_gap = get_CNT_gap($lastradius, $CC, $width);
  my $radius = $lastradius + 0.01; #arbitrary step to get process started
  my $gap = get_CNT_gap($radius, $CC, $width);

  #Simple Newton solver using very crude finite difference derivative
  while(abs($gap) > $THRESH) {
    my $newradius = $radius - $gap*($radius - $lastradius)/($gap - $old_gap);
    $old_gap = $gap;
    $gap = get_CNT_gap($newradius, $CC, $width);
    $lastradius = $radius;
    $radius = $newradius;
  }
  return $radius;
}

#Function to be minimized to get the radius of closed CNT
sub get_CNT_gap {
  my ($guess, $CC, $width) = @_;
  my $value = asin($CC/(2*$guess)) + asin($CC/(4*$guess)) - pi/(2*$width);
  return $value;
}
1;
