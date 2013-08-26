#!/usr/bin/perl
use 5.8.8;
use strict;
use warnings;

use Benchmark;
my $t0 = Benchmark->new;

# Author: Matt LaFave
# Written: 3/1/12

# This script looks through the second file and only keeps entries that have the
# same name as those in the first.  It's a way of only keeping the "overlap"
# of two fastq files.

################################################################################

my %names_hash;
my $name;
my $count;


unless (@ARGV == 2) {
	die "Usage: fastq_file_overlap_v1.0.pl fastq_to_compare fastq_to_reduce\n$!";
}

unless ( -e $ARGV[0] ) {
	die "$ARGV[0] doesn't exist: $!";
}

unless ( -e $ARGV[1] ) {
	die "$ARGV[1] doesn't exist: $!";
}

unless ( open FASTQ_COMP, '<', $ARGV[0]) {
	die "Can't open $ARGV[0]: $!";
}

# First, examine the "fastq_to_compare" file.  Put all the entry names into a
# hash.
while (<FASTQ_COMP>) {

	chomp;
	
	# If the line looks like a name, put it in the hash.  Otherwise, skip to the
	# next line.
	if (m%^(@(.*)/(1|2))$%) {
		
		$names_hash{$_} = 1;
		
		# Skip to what should be the next name line
		<FASTQ_COMP>;
		<FASTQ_COMP>;
		<FASTQ_COMP>;
		
	} else {
		
		next;
		
	}

}

close FASTQ_COMP;


unless ( open FASTQ_RED, '<', $ARGV[1]) {
	die "Can't open $ARGV[1]: $!";
}

unless ( open OUTPUT, '>', "$ARGV[1]_overlap.fastq") {
	die "Can't open $ARGV[1]_overlap.fastq: $!";
}

# Read through the other fastq file, and only keep entries that had an entry
# with the same name in the first fastq file.
while ( defined ($name=<FASTQ_RED>)) {
	
	chomp $name;

	if ( exists ($names_hash{$name})) {
	
		chomp (my $seq = <FASTQ_RED>);
		chomp (my $plus = <FASTQ_RED>);
		chomp (my $qual = <FASTQ_RED>);
		print OUTPUT "$name\n$seq\n$plus\n$qual\n";
		
		delete $names_hash{$name};
		$count = keys %names_hash;
		
		# There probably won't be an instance in which I exit early because of 
		# completely emptying the hash, but it won't hurt to put that 
		# functionality in.
		last if ($count == 0);
		
	}

}

close OUTPUT;

close FASTQ_RED;


my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);
print "The code took:",timestr($td),"\n";

exit;
