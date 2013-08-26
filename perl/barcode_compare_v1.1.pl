#!/usr/bin/perl
use 5.8.8;
use strict;
use warnings;

use Benchmark;
my $t0 = Benchmark->new;

# Author: Matt LaFave
# Written: 5/3/12
# Last revised: 6/11/12

# Designed to look through ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair,
# produced by bowtie_examine_p-pairs_v1.4.pl, and remove reads for which the
# barcodes are inconsistent between reads.

################################################################################

unless ( @ARGV == 3 ) {
	die "Usage: barcode_compare_v1.1.pl bowtie_file_nofarpair norm_barcode_cat p-pair_barcode_cat\n$!";
}

unless ( -e $ARGV[0] ) {
	die "$ARGV[0] doesn't exist: $!";
}

unless ( -e $ARGV[1] ) {
	die "$ARGV[1] doesn't exist: $!";
}

unless ( -e $ARGV[2] ) {
	die "$ARGV[2] doesn't exist: $!";
}

my %name_hash;
my %barcode_hash;

unless ( open BOWTIE, "<", "$ARGV[0]" ) {
	die "Cannot open file $ARGV[0]: $!";
}

# First, make a hash of names of entries for which barcodes should be obtained.
# Do this by reading through the "nofarpair" file, which indicates the 
# LTR-containing reads that are reasonably nearby their paired reads, and in
# correct orientation.  That file has a bowtie-output-style format, hence the
# filehandle.

while (<BOWTIE>) {
	
	chomp;
	my @line = split;
	
	# Detect the name, put it in the format of a fastq file
	my $name = '@'.$line[0];
	my $pair_name = $name;
	
	# Calculate the name of the read's "p-pair", the version of the paired read
	# that did not have low-quality bases removed during trimming.
	# Note that reads that already have "|p" should not be in the "nofarpair"
	# file, and thus don't need to be considered when finding the pair.
	if ( substr($pair_name, -1) eq '1' ) {
		substr($pair_name, -1) = '2|p';
	} elsif ( substr($pair_name, -1) eq '2' ) {
		substr($pair_name, -1) = '1|p';
	} else {
		next;
	}
	
	# Put the names of both the original read and the paired read in a hash.
	$name_hash{$name} = 1;
	$name_hash{$pair_name} = 1;
	
}

close BOWTIE;


# Next, go retrieve those barcode sequences.  Empty the name hash as you go.
unless ( open NORM_BARCODE, "<", "$ARGV[1]" ) {
	die "Cannot open file $ARGV[1]: $!";
}

# Read in the barcodes.  Note they're being read into a hash and stored as 
# values, with the names serving as keys.
while (<NORM_BARCODE>) {
	
	chomp;
	
	# Check to see if you're currently on a 'name' line.  If not, move on until
	# you reach a name line.
	if (m%^(@(.*)/(1|2))$%) {
		
		# Once you know you're on a name line, make sure it's a barcode that
		# needs to be recorded.  If not, skip to the next entry.
		if ( exists $name_hash{$_} ) {	
			
			my $name = $_;		
			chomp ( my $seq = <NORM_BARCODE>);
			<NORM_BARCODE>; # Skip the plus
			<NORM_BARCODE>; # Skip the quality line
			
			$barcode_hash{$name} = $seq;
			delete $name_hash{$_};
			
		} else {
			
			<NORM_BARCODE>; # Skip the sequence
			<NORM_BARCODE>; # Skip the plus
			<NORM_BARCODE>; # Skip the quality line
			
		}			
	
	} else {
	
		next;
		
	} # if

} # while

close NORM_BARCODE;


unless ( open P_BARCODE, "<", "$ARGV[2]" ) {
	die "Cannot open file $ARGV[2]: $!";
}

# Read in the p-pair barcodes.
while (<P_BARCODE>) {
	
	chomp;
	
	# Check to see if you're currently on a 'name' line.  If so, read in the
	# rest of the entry; if not, move on until you reach a name line.
	if (m%^(@(.*)/(1|2)(\|p)?)$%){
		
		if ( exists $name_hash{$_} ) { 
			
			my $name = $_;		
			chomp ( my $seq = <P_BARCODE>);
			<P_BARCODE>; # Skip the plus
			<P_BARCODE>; # Skip the quality line
			
			$barcode_hash{$name} = $seq;
			delete $name_hash{$_};
			
		} else {
			
			<P_BARCODE>; # Skip the sequence
			<P_BARCODE>; # Skip the plus
			<P_BARCODE>; # Skip the quality line
			
		}
	
	} else {
	
		next;
		
	} # if

} # while

close P_BARCODE;

# If there's anything left in the %name_hash, remove it.
%name_hash = ();


# Look through the bowtie file again.  For reads that have two barcodes, verify 
# that they're the same.

unless ( open BOWTIE, "<", "$ARGV[0]" ) {
	die "Cannot open file $ARGV[0]: $!";
}

unless ( open OUTPUT, ">", "$ARGV[0]_consistentbc" ) {
	die "Cannot open file $ARGV[0]_consistentbc: $!";
}

unless ( open BADBC, ">", "$ARGV[0]_inconsistentbc" ) {
	die "Cannot open file $ARGV[0]_inconsistentbc: $!";
}


while (<BOWTIE>) {
	
	chomp;
	my @line = split;
	
	my $name = '@'.$line[0];
	
	my $pair_name = $name;
	
	# The input doesn't have p-pairs in it, so I don't need to worry about
	# running across an entry like that.  If I did, it'd be skipped.
	if ( substr($pair_name, -1) eq '1' ) {
		substr($pair_name, -1) = '2|p';
	} elsif ( substr($pair_name, -1) eq '2' ) {
		substr($pair_name, -1) = '1|p';
	} else {
		next;
	}
	
	# Make sure a read has a barcode to compare to.
	if ( ( exists $barcode_hash{$name}) && ( exists $barcode_hash{$pair_name}) ) {
		
		if ( $barcode_hash{$name} eq $barcode_hash{$pair_name} ) {
			
			# If the barcodes match, good.  Print the current line, and may
			# as well delete both entries from the hash, to make it smaller.
			print OUTPUT "$_\n";
			delete $barcode_hash{$name};
			delete $barcode_hash{$pair_name};
			
		} else {
			
			# If the barcodes DON'T match, only print it to the 'bad' output.
			# Also print the offending barcodes.
			print BADBC "$_\tBarcode:$barcode_hash{$name}\tp-pair_barcode:$barcode_hash{$pair_name}\n";
			delete $barcode_hash{$name};
			delete $barcode_hash{$pair_name};

		}
		
	} else {
		
		# If a read has no paired barcode (as will be the case for long 
		# fragments), or doesn't have a barcode itself (the linker side of a 
		# long fragment), it can't have inconsistent barcodes.
		# Print it & remove it from the hash.
		print OUTPUT "$_\n";
		if ( exists $barcode_hash{$name} ) {
			delete $barcode_hash{$name};
		}
		if ( exists $barcode_hash{$pair_name} ) {
			delete $barcode_hash{$pair_name};
		} 
		
	}
	
}

close BADBC;

close OUTPUT;

close BOWTIE;



my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);
print "The code took:",timestr($td),"\n";

exit;
