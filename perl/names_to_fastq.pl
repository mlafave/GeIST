#!/usr/bin/perl
use 5.8.8;
use strict;
use warnings;

use Benchmark;
my $t0 = Benchmark->new;

# Author: Matt LaFave
# Written: 8/17/12, as a modified version of fastq_file_overlap.pl.

my %names_hash;
my $name;
my $count;


unless (@ARGV == 2) {
	die "Usage: fastq_file_overlap.pl name_list fastq_to_reduce\n$!";
}

unless ( -e $ARGV[0] ) {
	die "$ARGV[0] doesn't exist: $!";
}

unless ( -e $ARGV[1] ) {
	die "$ARGV[1] doesn't exist: $!";
}

unless ( open NAMES, '<', $ARGV[0]) {
	die "Can't open $ARGV[0]: $!";
}


while (<NAMES>) {

	chomp;
	
	$names_hash{$_} = 1;

}

close NAMES;


unless ( open FASTQ, '<', $ARGV[1]) {
	die "Can't open $ARGV[1]: $!";
}

unless ( open OUTPUT, '>', "$ARGV[0].fastq") {
	die "Can't open $ARGV[0].fastq: $!";
}


while ( defined ($name=<FASTQ>)) {
	
	chomp $name;

	if ( exists ($names_hash{$name})) {
	
		chomp (my $seq = <FASTQ>);
		chomp (my $plus = <FASTQ>);
		chomp (my $qual = <FASTQ>);
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

close FASTQ;


my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);
print "The code took:",timestr($td),"\n";

exit;
