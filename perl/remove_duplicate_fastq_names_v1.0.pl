#!/usr/bin/perl
use 5.8.8;
use strict;
use warnings;

use Benchmark;
my $t0 = Benchmark->new;

# Author: Matt LaFave
# Written: 2/19/12

# Removes entries with repeated names (2 or more times) from a fastq file.
# In the main shell script, it's used to remove reads that have inverted LTR or
# inverted linker (that is, two occurences of either sequence on a single read,
# in opposite orientations).

# General procedure: "while" through the fastq file, only looking at names.  Add
# each to a hash. If the name is already in the hash, add it to a "duplicate 
# names" hash.

################################################################################


my %initial_hash;
my %duplicate_hash;
my $name;
my $dup_count = 0;

unless ( @ARGV == 1 ) {
	die "Usage: remove_duplicate_fastq_names_v1.0.pl fastq_file\n$!";
}

unless ( -e $ARGV[0] ) {
	die "$ARGV[0] doesn't exist: $!";
}

unless ( open FASTQIN, "<", "$ARGV[0]" ) {
	die "Cannot open file $ARGV[0]: $!";
}

while ( defined ($name = <FASTQIN>) ) {
	
	chomp $name;
	
	# If it looks like a name, operate on it.  Otherwise, check the next line.
	if ( $name =~ m%^(@(.*)/(1|2))$%) {
		
		# If the name has been seen before, put it in the dulicate hash.  
		# Otherwise, put it in the initial hash.
		if ( exists $initial_hash{$name}) {
			
			$duplicate_hash{$name} = 1;
			
		} else {
			
			$initial_hash{$name} = 1;
			
		}
		
	} else {
	
		next;
	
	}
	
	# Skip the sequence line, plus line, and quality line.
	<FASTQIN>;
	<FASTQIN>;
	<FASTQIN>;
	
}

close FASTQIN;

# "while" through the fastq file again; if a name isn't in the duplicate names 
# hash, print it, and the rest of the entry, to output.  Otherwise, don't print 
# it, and skip to the next entry.

unless ( open FASTQIN, "<", "$ARGV[0]" ) {
	die "Cannot open file $ARGV[0]: $!";
}

unless ( open OUTPUT, ">", "$ARGV[0]_nodup.fastq" ) {
	die "Cannot open file $ARGV[0]_nodup.fastq: $!";
}

while (defined ($name = <FASTQIN>)) {
	
	chomp $name;
	
	# Make sure you're reading a name line; if not, go to the next line.
	if ( $name =~ m%^(@(.*)/(1|2))$%) {
	
		if ( exists $duplicate_hash{$name} ) {
			
			# Keep track of how many duplicate names are in the file.
			$dup_count++;
			<FASTQIN>;
			<FASTQIN>;
			<FASTQIN>;
			
		} else {
			
			# If it's not in the duplicate hash, print it to the output.
			chomp (my $seq = <FASTQIN>);
			chomp (my $plus = <FASTQIN>);
			chomp (my $qual = <FASTQIN>);
			print OUTPUT "$name\n$seq\n$plus\n$qual\n";
			
		}
		
	} else {
		
		next;
		
	}
	
}

close OUTPUT;

close FASTQIN;

print "$ARGV[0] had $dup_count duplicate names.\n";

my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);
print "The code took:",timestr($td),"\n";

exit;
