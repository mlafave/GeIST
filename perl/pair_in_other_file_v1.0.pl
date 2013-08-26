#!/usr/bin/perl
use 5.8.8;
use strict;
use warnings;

use Benchmark;
my $t0 = Benchmark->new;

# Author: Matt LaFave
# Written: 2/29/12

# Used to reduce a fastq file to only the reads that have paired reads in 
# another fastq file.

# Seaches the linker_F file to make sure that the reads in LTR_F_long_nobadlink 
# have  forward-facing linker on their paired read.  If they don't, delete the 
# name from the LTR_F hash.

################################################################################


my %linker_F_pair_hash;
my $LTR_line;
my $count;


unless (@ARGV == 2) {
	die "Usage: pair_in_other_file_v1.0.pl linker_F_nodup LTR_F_long_nobadlink\n$!";
}

unless ( -e $ARGV[0] ) {
	die "$ARGV[0] doesn't exist: $!";
}

unless ( -e $ARGV[1] ) {
	die "$ARGV[1] doesn't exist: $!";
}

unless ( open LINKER_F, '<', $ARGV[0]) {
	die "Can't open $ARGV[0]: $!";
}

# First, put the names of the linker_F read in a hash.  It's not really the true
# names, though;  it's the names of the pair.  So, if linker_F has /2, that 
# means we want to fish out /1 from LTR_F (if LTR_F has it).  So, we put /1 into
# the hash now.

while (<LINKER_F>) {

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
		
		# Put the pair name in a hash as a key.
		$linker_F_pair_hash{$pair_name} = 1;
		
		# Skip to what should be the next name line
		<LINKER_F>;
		<LINKER_F>;
		<LINKER_F>;
		
	} else {
		
		next;
		
	}

}

close LINKER_F;

unless ( open LTR_F, '<', $ARGV[1]) {
	die "Can't open $ARGV[1]: $!";
}

unless ( open OUTPUT, '>', "$ARGV[1]_pairhaslink.fastq") {
	die "Can't open $ARGV[1]_pairhaslink.fastq: $!";
}

# Next, look through LTR_F, looking for names that are in the hash.  These are 
# the names of reads that have forward-facing linker on their paired read.

while ( defined ($LTR_line = <LTR_F>)) {
	
	chomp $LTR_line;
	
	# It will only match if the current line is a name, so there's no need to 
	# check that again.
	if ( exists ($linker_F_pair_hash{$LTR_line})) {
		
		chomp (my $seq = <LTR_F>);
		chomp (my $plus = <LTR_F>);
		chomp (my $qual = <LTR_F>);
		print OUTPUT "$LTR_line\n$seq\n$plus\n$qual\n";
		
		# It's very unlikely that I'll be able to exit early from completely 
		# emptying the hash, but just in case, I'll add in that functionality:
		delete $linker_F_pair_hash{$LTR_line};
		$count = keys %linker_F_pair_hash;
		last if ($count == 0);
		
		
	} else {
		
		next;
		
	}

}

close OUTPUT;

close LTR_F;


my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);
print "The code took:",timestr($td),"\n";

exit;
