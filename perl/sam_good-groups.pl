#!/usr/bin/perl
use 5.8.8;
#use strict;
#use warnings;

# Author: Matt LaFave
# Created: 12/5/12

# Goes through a bowtie-style output file that has barcode group information 
# within, and prints out those that are permissible as defined by a cutoff file.

# It's very similar to notinx_to_bowtie_v1.0.pl.

################################################################################

unless ( @ARGV == 2 ) {
	die "Usage: sam_good-groups.pl cutoff_file bowtie_with_groups\n$!";
}

unless ( -e $ARGV[0] ) {
	die "$ARGV[0] doesn't exist: $!";
}

unless ( -e $ARGV[1] ) {
	die "$ARGV[1] doesn't exist: $!";
}

# Declare global variables
my %position_hash;


unless ( open CUTOFF, "<", "$ARGV[0]" ) {
	die "Cannot open file $ARGV[0]: $!";
}

while (<CUTOFF>) {
	
	chomp;
	my @line = split;
	
	# Put the position, chromosome, and group in a hash.  Since this will be 
	# the key, and must therefore be unique, concatenate them with a hyphen.
	my $cat = "$line[2]-$line[1]-$line[3]";
	
	$position_hash{$cat} = 1;
	
}

close CUTOFF;

unless ( open BOWTIEIN, "<", "$ARGV[1]" ) {
	die "Cannot open file $ARGV[1]: $!";
}

# Read through the bowtie file.
while (<BOWTIEIN>) {
	
	chomp;
	my @bowline = split;
	
	# Make the same combination of position, chromosome, and group as above.
	my $bowcat = "$bowline[3]-$bowline[2]-$bowline[17]";
	
	if ( exists $position_hash{$bowcat} ) {
		
		# If the aligned read survived remove_nearby_v3.0.pl, print it back out
		# in bowtie-output-style format.
		print "$_\n";
		
	} # if
	
} # while

close BOWTIEIN;


exit;
