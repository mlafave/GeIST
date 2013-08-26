#!/usr/bin/perl
use 5.8.8;
use strict;
use warnings;

use Benchmark;
my $t0 = Benchmark->new;

# Author: Matt LaFave
# Written: 3/16/12

# This script is a slight modification of another script called 
# pair_in_other_file_v1.0.pl. It's designed to return the fastq entries corresponding
# to the pairs of entries of a given fastq file.  The two scripts are close, if
# not identical in function at this point (see next paragraph), but the 
# filehandles used in pair_in_other_file_v1.0.pl are a little confusing in 
# certain contexts, so I've kept this script around, as well.

# This version was originally designed to get the pairs of reads, and add "|p" 
# to the ends of those reads.  However, after discovering that Bowtie alignments
# can differ based on the name of the read, I removed the "|p"-adding 
# functionality and moved it downstream, to the "add_p_v1.0.pl script. 

################################################################################


my %pair_hash;
my $name;
my $count;


unless (@ARGV == 2) {
	die "Usage: fastq_get_pair_v1.0.pl trimmed_fastq original_fastq\n$!";
}

unless ( -e $ARGV[0] ) {
	die "$ARGV[0] doesn't exist: $!";
}

unless ( -e $ARGV[1] ) {
	die "$ARGV[1] doesn't exist: $!";
}

unless ( open INPUT, '<', $ARGV[0]) {
	die "Can't open $ARGV[0]: $!";
}

# First, put the names of the trimmed fastq entries in a hash.  It's not really 
# the true names, though;  it's the names of the pair.  So, if the input has /2,
# that means we want to fish out /1 from NIOBE.fastq.  So, we put /1 into the 
# hash now.

while (<INPUT>) {

	chomp;
	
	# If the line looks like a name, put its pair in the hash.  Otherwise, skip
	# to the next line.
	if (m%^(@(.*)/(1|2))$%) {
		
		my $pair_name = $_;
	
		if ( substr($pair_name, -1) == 1) {
			substr($pair_name, -1) = 2;
		} elsif (substr($pair_name, -1) == 2) {
			substr($pair_name, -1) = 1;
		} else {
			# If the name doesn't end with a 1 or a 2, it's weird, and shouldn't
			# have made it this far. Go to the next entry.
			next;
		}
		
		# Record the pair name in a hash
		$pair_hash{$pair_name} = 1;
		
		# Skip to what should be the next name line
		<INPUT>;
		<INPUT>;
		<INPUT>;
		
	} else {
		
		next;
		
	}

}

close INPUT;

unless ( open ORIGINAL, '<', $ARGV[1]) {
	die "Can't open $ARGV[1]: $!";
}

unless ( open OUTPUT, '>', "$ARGV[0]_pairs-with-p.fastq") {
	die "Can't open $ARGV[0]_pairs-with-p.fastq: $!";
}

# Next, look through the original fastq file, looking for names that are in the 
# hash.  
while ( defined ($name = <ORIGINAL>)) {
	
	chomp $name;
	
	# It will only match if the current line is a name, so there's no need to 
	# check that again.
	if ( exists ($pair_hash{$name})) {
		
		chomp (my $seq = <ORIGINAL>);
		chomp (my $plus = <ORIGINAL>);
		chomp (my $qual = <ORIGINAL>);
		print OUTPUT "$name\n$seq\n$plus\n$qual\n";
		
		# If the hash empties before you're done looking through the file, exit
		# the loop.
		delete $pair_hash{$name};
		$count = keys %pair_hash;
		last if ($count == 0);
		
		
	} else {
		
		next;
		
	}

}

close OUTPUT;

close ORIGINAL;


my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);
print "The code took:",timestr($td),"\n";

exit;
