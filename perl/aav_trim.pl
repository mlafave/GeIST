#!/usr/bin/perl
use 5.8.8;
use strict;
use warnings;

# Author: Matt LaFave
# Written: 1/15/13
# Last modified: 3/1/13 (revision 27*)

# This script is used to trim the AAV sequence off the leftmost side of 
# sequencing reads.

################################################################################

# IF you decide to make this so it can deal with reads in which the AAV primer
# was not detected by cutadapt:
# Use a command-line flag to tell the script if the AAV primer had been directly
# adjacent to the sequence to trim, or not.
# Until then, assume that AAV had been detected an trimmed.

# Expect a FASTQ input.

unless ( @ARGV == 2 ) {
	die "Usage: aav_trim.pl fastq primer_detected?\n$!";
}

unless ( -e $ARGV[0] ) {
	die "$ARGV[0] doesn't exist: $!";
}

my $primer_flag;
if ( $ARGV[1] eq 'yes') {
	# If the primer WAS detected, set the flag.
	$primer_flag = 1;
} elsif ($ARGV[1] eq 'no') {
	# If the primer was not detected, the flag is 0.
	$primer_flag = 0;
} else {
	die "$ARGV[1] needs to be 'yes' or 'no': $!";
}

# Define global variables

my $itr = "TTGGCCACTCCCTCTCTGCGCGCTCGCTCGCTCACTGAGGCCGGGCGACCAAAGGTCGCCCGACGCCCGGGCTTTGCCCGGGCGGCCTCAGTGAG";
# The ITR is shown as the reverse complement, since incoming sequences will have
# had AAV on the right side.
# Beware of case-sensitivity when doing the match.

# Use a for loop to pull off bases one at a time, and see if they match the AAV
# sequence.

unless ( open FASTQ, "<", "$ARGV[0]" ) {
	die "Cannot open file $ARGV[0]: $!";
}

unless ( open TRIMOUT, ">", "$ARGV[0]_aavtrim.fastq" ) {
	die "Cannot open file $ARGV[0]_aavtrim.fastq: $!";
}

unless ( open AAVOUT, ">", "$ARGV[0]_itrbits.fastq" ) {
	die "Cannot open file $ARGV[0]_itrbits.fastq: $!";
}

if ($primer_flag){
	
	# Primer was detected; AAV and sequence read are both chewed from the right.

	while (<FASTQ>) {
		
		my ($name, $seq, $plus, $qual);
		
		chomp;
		
		# Make sure you're on a name line
		if (m%^(@(.*)/(1|2)(\|p)?)$%) {
			
			$name = $_;
			chomp ($seq = <FASTQ>);
			chomp ($plus = <FASTQ>);
			chomp ($qual = <FASTQ>);
		
		} else {
			
			next;
			
		}
		
		# Peel off bases from the rightmost side, one at a time
		for ( my $i=1; $i <= length($itr); $i++ ) {
			
			my $seq_tail = substr($seq, -${i});
			my $aav_tail = substr($itr, -${i});
			
			# If it matches, good, keep going
			if ( $seq_tail =~ m/^($aav_tail)$/i ) {
				next;
			} else {
				
				# If it doesn't match, you've hit the end.  Use i-1 for the cut.
				
				# First, make sure there's enough DNA left to map.  If there is less
				# than 11 bp, go to the next FASTQ entry.
				if ( (length($seq)-($i-1)) < 11 ) {
					last;								
				}
				
				# Remove i-1 bases from the right side of the sequence and quality
				# score
				my $trimmed_seq = substr($seq, 0, ( length($seq)-($i-1) ) );
				my $trimmed_qual = substr($qual, 0, ( length($qual)-($i-1) ) );
				
				print TRIMOUT "$name\n$trimmed_seq\n$plus\n$trimmed_qual\n";
				
				# Print the portion that was trimmed to a seperate file, IF anything
				# was trimmed.
				if ( ($i-1) > 0 ) {
					my $itr_seq = substr($seq, -($i-1));
					my $itr_qual = substr($qual, -($i-1));
					
					print AAVOUT "$name\n$itr_seq\n$plus\n$itr_qual\n";
					
				}
				
				# The current entry is done; move on to the next one.
				last;
							
			}
			
		}
		
		
	}

} else {
	
	# If the primer was not detected, the read gets chewed from the right, but
	# the AAV bit just needs to have contiguous bases.
	
	while (<FASTQ>) {
		
		my ($name, $seq, $plus, $qual);
		
		chomp;
		
		# Make sure you're on a name line
		if (m%^(@(.*)/(1|2)(\|p)?)$%) {
			
			$name = $_;
			chomp ($seq = <FASTQ>);
			chomp ($plus = <FASTQ>);
			chomp ($qual = <FASTQ>);
		
		} else {
			
			next;
			
		}
		
		# Peel off bases from the rightmost side of the READ only, one at a time
		for ( my $i=1; $i <= length($itr); $i++ ) {
			
			my $seq_tail = substr($seq, -${i});
						
			# If it matches, good, keep going
			if ( $itr =~ m/($seq_tail)/i ) {
				next;
			} else {
				
				# If it doesn't match, you've hit the end.  Use i-1 for the cut.
				
				# First, make sure there's enough DNA left to map.  If there is less
				# than 11 bp, go to the next FASTQ entry.
				if ( (length($seq)-($i-1)) < 11 ) {
					last;								
				}
				
				# Remove i-1 bases from the right side of the sequence and quality
				# score
				my $trimmed_seq = substr($seq, 0, ( length($seq)-($i-1) ) );
				my $trimmed_qual = substr($qual, 0, ( length($qual)-($i-1) ) );
				
				print TRIMOUT "$name\n$trimmed_seq\n$plus\n$trimmed_qual\n";
				
				# Print the portion that was trimmed to a seperate file, IF anything
				# was trimmed.
				if ( ($i-1) > 0 ) {
					my $itr_seq = substr($seq, -($i-1));
					my $itr_qual = substr($qual, -($i-1));
					
					print AAVOUT "$name\n$itr_seq\n$plus\n$itr_qual\n";
					
				}
				
				# The current entry is done; move on to the next one.
				last;
							
			}
			
		}
		
		
	}
	
}

close AAVOUT;

close TRIMOUT;

close FASTQ;

exit;
