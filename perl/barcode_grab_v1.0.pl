#!/usr/bin/perl
use 5.8.8;
use strict;
use warnings;

use Benchmark;
my $t0 = Benchmark->new;

# Author: Matt LaFave
# Written: 4/2/12

# This script records and removes barcodes from a fastq file.  HOWEVER, it
# expects the barcode to be the first seven bases from the RIGHT of the read. 
# That may mean you need to flip a file before providing it to the script, and
# flip it again when it's done - however, since the file would likely have been
# flipped for cutadapt anyway, you probably won't need to do it twice.

################################################################################

unless ( @ARGV == 2 ) {
	die "Usage: barcode_grab_v1.0.pl barcode_ref fastq\n$!";
}

unless ( -e $ARGV[0] ) {
	die "$ARGV[0] doesn't exist: $!";
}

unless ( -e $ARGV[1] ) {
	die "$ARGV[1] doesn't exist: $!";
}

my %barcode_hash;
my ($name, $seq, $plus, $qual);

# First, make a hash of acceptable barcodes.  Also uses a check such that,
# if you don't KNOW which barcodes are permissible, it just lets everything 
# through.

unless ( open BARCODEREF, "<", "$ARGV[0]" ) {
	die "Cannot open file $ARGV[0]: $!";
}

while (<BARCODEREF>) {
	
	chomp;
	# Use the value of the hash to keep track of how often that barcode was
	# detected, so initialize to 0.
	$barcode_hash{$_} = 0;
	
}

close BARCODEREF;



# Walk through the file, grabbing the name and the sequence line (and the 
# quality line, while you're at it - might be interesting to see, later).  
# Use substr to get the  7 bp and store it in a variable ($code), and use substr
# again to get the remaining bases & store that in another variable ($trimseq).
# Check the barcode against the list of acceptable barcodes, and if it passes,
# print the barcode to a barcode file, and the sequence to an output file.
# There's also a third file that shows how often each barcode was detected.

unless ( open FASTQ, "<", "$ARGV[1]" ) {
	die "Cannot open file $ARGV[1]: $!";
}

unless ( open BAROUT, ">", "$ARGV[1]_barcodes.fastq" ) {
	die "Cannot open file $ARGV[1]_barcodes.fastq: $!";
}

unless ( open TRIMMEDOUT, ">", "$ARGV[1]_bctrimmed.fastq" ) {
	die "Cannot open file $ARGV[1]_bctrimmed.fastq: $!";
}

unless ( open UNTRIMMED, ">", "$ARGV[1]_untrimmed.fastq" ) {
	die "Cannot open file $ARGV[1]_untrimmed.fastq: $!";
}

while (<FASTQ>) {
	
	chomp;
	
	# Check to see if you're currently on a 'name' line.  I included the
	# |p part because this needs to be flexible enough to be used in the
	# compare-p section of the main shell script.
	if (m%^(@(.*)/(1|2)(\|p)?)$%) {
		
		$name = $_;
		chomp ($seq = <FASTQ>);
		chomp ($plus = <FASTQ>);
		chomp ($qual = <FASTQ>);
	
	} else {
	
		next;
	
	} # if
	
	# Capture the barcode sequence
	my $code = substr($seq, -7);
	
	if ( exists $barcode_hash{$code} ) {
		
		# If the barcode is in the hash, good.  Print the barcode to one file,
		# the rest of the sequence to another, and increment the value in the
		# hash.
		my $code_qual = substr($qual, -7);
		
		print BAROUT "$name\n$code\n$plus\n$code_qual\n";
		
		my $trimseq = substr($seq, 0, (length($seq)-7));
		my $trimqual = substr($qual, 0, (length($qual)-7));
		
		print TRIMMEDOUT "$name\n$trimseq\n$plus\n$trimqual\n";
		
		$barcode_hash{$code}++;
		
	} else {
		
		# If the barcode is NOT in the hash, print those reads to a separate
		# file, unmodified. 
		
		print UNTRIMMED "$name\n$seq\n$plus\n$qual\n";
		
	} # if
	
	
} # while

close UNTRIMMED;

close TRIMMEDOUT;

close BAROUT;

close FASTQ;


unless ( open SUMMARY, ">", "$ARGV[1]_bcsummary" ) {
	die "Cannot open file $ARGV[1]_bcsummary: $!";
}

while ( my ($key, $value) = each %barcode_hash ) {
	
	# Print a summary of how many times each barcode was detected
	print SUMMARY "$key\t$value\n";

}

close SUMMARY;

################################################################################

my $t1 = Benchmark->new;
my $td = timediff($t1, $t0);
print "The code took:",timestr($td),"\n";

exit;
