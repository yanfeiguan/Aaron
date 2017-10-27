#!/usr/bin/perl -w
#Orient Aarom TS geometry so that Si-bonded atoms (either O,O or N,O) are along axes
#First, drop atom 1 to xy plane
#Rotate 1 to x-axis by rotation about z-axis
#Rotate 2 to y-axis by rotation about x-axis
#ie: assumes that the oxygens are atoms 24 and 25!

use strict;
use Math::Trig;

my @coords;
my @geom;
my $numatoms;
my $comment;

#Read XYZ file
while (<>) {
        chomp;
        if(/[A-Z]? /) {
                push(@geom, $_);
        } elsif(/^\d/) {
                $numatoms = $_;
		$comment = (<>);
		chomp $comment;
        } 
}

#parse coordinates
foreach (@geom) {
	my @line = split;
	push (@coords, \@line);
}

#identify silicon
my $silicon;
foreach my $atom (0..$#coords) {
	if($coords[$atom][0] =~ /Si/) {
		$silicon = $atom;
	}
}
#shift Si to origin
my @origin = ($coords[$silicon][0], $coords[$silicon][1], $coords[$silicon][2], $coords[$silicon][3]);
foreach my $atom (0..$#coords) {
	foreach my $j (1..3) {
		$coords[$atom][$j] -= $origin[$j];
	}
}

my $phi = atan2($coords[23][3],$coords[23][1]);
rotate(1,-$phi);
$phi = atan2($coords[23][2],$coords[23][1]);
rotate(2,-$phi);
$phi = atan2($coords[24][3],$coords[24][2]);
rotate(0,-$phi);
printgeom(@coords);


sub rotate {
	my $axis = $_[0];
	my $angle = $_[1];
	
	if($axis==0) { #X-axis rotation)
		foreach my $atom (0..$#coords) {
			my $y = $coords[$atom][2];
			my $z = $coords[$atom][3];
			$coords[$atom][2] = $y*cos($angle)-$z*sin($angle);
			$coords[$atom][3] = $y*sin($angle)+$z*cos($angle);
		}
		return;
	} 
	if($axis==1) { #Y-axis rotation
		foreach my $atom (0..$#coords) {
			my $x = $coords[$atom][1];
			my $z = $coords[$atom][3];
			$coords[$atom][1] = $x*cos($angle)-$z*sin($angle);
			$coords[$atom][3] = $x*sin($angle)+$z*cos($angle);
		}
		return;
	} 
	if($axis==2) { #Z-axis rotation
		foreach my $atom (0..$#coords) {
			my $x = $coords[$atom][1];
			my $y = $coords[$atom][2];
			$coords[$atom][1] = $x*cos($angle)-$y*sin($angle);
			$coords[$atom][2] = $x*sin($angle)+$y*cos($angle);
		}
		return;
	} 
	die "Attempt at rotation beyond 3-D\n";
}

sub printgeom {
	my @coords = @_;
	print "$numatoms\n$comment\n";
	foreach my $atom (0..$#coords) {
		printf "%s ", $coords[$atom][0];
		foreach my $j (1..3) {
			printf "    %10.6f", $coords[$atom][$j];
		}
		print "\n";
	}
}


