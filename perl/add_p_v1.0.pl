#!/usr/bin/perl
use 5.8.8;

# Author: Matt LaFave
# Written: 5/22/12

# Designed to add the characters "|p" to the end of the name of each entry in
# a FASTQ file.  Assumes that the first line of the file will be the name of an
# entry.

# USAGE: cat ${fastq_file}_p-pair_barcodes.fastq | ./add_p_v1.0.pl > ${fastq_file}_p-pair_barcodes.fastq_withp

while (<>) {
	
	chomp ( my $name = $_);
	chomp ( my $seq = <>);
	chomp ( my $plus = <>);
	chomp ( my $qual = <>);
	print "$name|p\n$seq\n$plus\n$qual\n";
	
} 

exit;
