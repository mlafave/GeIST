#!/usr/bin/perl
use 5.8.8;
use strict;
use warnings;

use Benchmark;
my $t0 = Benchmark->new;

# Author: Matt LaFave
# Written: 7/13/12

# Designed to append "barcode group" information to a bowtie-style file that
# already has barcode information in it (like an _island file).  Given the 
# barcode sequence, this script identifies which group the read belongs in.
# The version originally run in the main shell script broke reads into four
# groups, corresponding to the four original T75 flasks.

# This script assumes the bowtie file has 11 columns; group information will 
# become the 12th.  Barcode information is assumed to be in the 10th column.

################################################################################

unless ( @ARGV == 2 ) {
	die "Usage: bowtie_barcode-group_v1.0.pl barcode_groups bowtie_file\n$!";
}

unless ( -e $ARGV[0] ) {
	die "$ARGV[0] doesn't exist: $!";
}

unless ( -e $ARGV[1] ) {
	die "$ARGV[1] doesn't exist: $!";
}

my %group_hash;

unless ( open GROUP, "<", "$ARGV[0]" ) {
	die "Cannot open file $ARGV[0]: $!";
}

# Make a hash with barcodes as keys, and group numbers as values.

while (<GROUP>) {
	
	chomp;
	my @line = split;
	
	# The first column of the file is the barcode (key); the second is the group
	# number (value)
	$group_hash{$line[0]} = $line[1];
	
}

close GROUP;


unless ( open BOWTIE, "<", "$ARGV[1]" ) {
	die "Cannot open file $ARGV[1]: $!";
}

unless ( open OUTPUT, ">", "$ARGV[1]_grouped" ) {
	die "Cannot open file $ARGV[1]_grouped: $!";
}

while (<BOWTIE>) {
	
	chomp;
	my @line = split;
	
	# $line[9] should hold the barcode information.
	# At this point, every barcode in the bowtie file should be in the hash, but
	# it won't hurt to check.
	if (exists $group_hash{$line[15]}) {
		
		print OUTPUT "$_\t$group_hash{$line[15]}\n";
		
	} else {
		
		# If there's a barcode in the bowtie file that isn't in the hash, 
		# print NO_GROUP.
		print OUTPUT "$_\tNO_GROUP\n";
		
	}	
	
}



my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);
print "The code took:",timestr($td),"\n";

exit;
