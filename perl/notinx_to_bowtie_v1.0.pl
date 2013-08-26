#!/usr/bin/perl
use 5.8.8;
use strict;
use warnings;

use Benchmark;
my $t0 = Benchmark->new;

# Author: Matt LaFave
# Written: 4/4/12

# Converts the output of remove_nearby_v3.0.pl (which I refer to as a "not in X"
# file, where X is the number of bases around an integration in which no other
# integrations in the same group is allowed) to a format based on the output of
# Bowtie.

################################################################################

# The intent is that the bowtie file used will be _nofarpair.
unless ( @ARGV == 2 ) {
	die "Usage: notinx_to_bowtie_v1.0.pl not_in_file bowtie_file\n$!";
}

unless ( -e $ARGV[0] ) {
	die "$ARGV[0] doesn't exist: $!";
}

unless ( -e $ARGV[1] ) {
	die "$ARGV[1] doesn't exist: $!";
}

# Declare global variables
my %position_hash;


unless ( open NOTINX, "<", "$ARGV[0]" ) {
	die "Cannot open file $ARGV[0]: $!";
}

while (<NOTINX>) {
	
	chomp;
	my @line = split;
	
	# Put both the position and the chromosome in a hash.  Since this will be 
	# the key, and must therefore be unique, concatenate them with a hyphen.
	my $cat = "$line[2]-$line[1]";
	
	$position_hash{$cat} = 1;
	
}

close NOTINX;


unless ( open BOWTIEIN, "<", "$ARGV[1]" ) {
	die "Cannot open file $ARGV[1]: $!";
}

unless ( open BOWTIEOUT, ">", "$ARGV[1]_island" ) {
	die "Cannot open file $ARGV[1]: $!";
}

# Read through the bowtie file, presumably the output of 
# bowtie_examine_pairs_v1.4.pl of good reads without p-pairs represented.
while (<BOWTIEIN>) {
	
	chomp;
	my @bowline = split;
	
	# Make the same combination of position and chromosome, as above.
	my $bowcat = "$bowline[3]-$bowline[2]";
	
	if ( exists $position_hash{$bowcat} ) {
		
		# If the aligned read survived remove_nearby_v3.0.pl, print it back out
		# in bowtie-output-style format.
		print BOWTIEOUT "$_\n";
		
	} # if
	
} # while

close BOWTIEOUT;

close BOWTIEIN;

my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);
print "The code took:",timestr($td),"\n";

exit;
