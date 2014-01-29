#! /bin/bash

# Email address, and options, to whom reports should be sent.
# The -M is for the address, the -m is for when to send email.
# a=abort b=begin e=end s=suspend n=none
#$ -M matthew.lafave@nih.gov
#$ -m abes

# Redirect the STDOUT and STDERR files to the jobs directory
#$ -o $HOME/jobs/$JOB_ID_$JOB_NAME.stdout
#$ -e $HOME/jobs/$JOB_ID_$JOB_NAME.error

# Operate in the current directory
#$ -cwd

# Author: Matt LaFave
# Created: 8/9/13 (based on mlv_zf_int_driver.sh (as of 4/3/13, r31*), which was
# based on tol2_driver_4.0.sh, which was based on ds_integration_driver.sh; all 
# eventually owe ancestry to GeIST.sh)


# The script assumes it's operating on an SGE cluster, submitted via qsub, but
# it isn't required.

# This script is designed to identify MLV, HIV, or Foamy virus integrations, and
# report them as a BED file. A tarball containing files with more detailed 
# information on every read is also produced.  If access to the intermediate 
# files is required, entering "on", "yes", "y", or "true" on the command line as
# an optional sixth argument will turn on recovery of the informative 
# intermediate files.  Said files are produced as a separate tarball.

# Set up functions for file testing & error reporting.
function throw_error
{
  echo >&2 ERROR: $1
  exit 1
}

function test_file
{
 if 
   [ -f $1 ]
 then 
   echo "$1 produced."
 else  
   throw_error "$1 was not produced!"
 fi
}

# Verify that the necessary programs are installed (also hash commands for faster 
# lookup later)
hash cutadapt 2>/dev/null || throw_error "cutadapt not found"
hash bowtie 2>/dev/null || throw_error "bowtie not found"
hash samtools 2>/dev/null || throw_error "samtools not found"

# Verify that the programs to be called are of the correct version
ver=`cutadapt --version`
if 
   awk -v ver="$ver" 'BEGIN {exit !(ver >= 1.3)}' 
then 
   echo "cutadapt is version $ver"
else
   throw_error "cutadapt is not at least version 1.3."
fi
ver=0

ver=`bowtie --version | head -1 | awk '{print $3}'` 
echo "bowtie is version $ver ( >= v. 1.0.0 is preferred)"
ver=0

ver=`samtools 2>&1 | head -3 | awk '/^Version/ {print $2}' | sed 's/^\(.*\)-.*$/\1/'`
echo "samtools is version $ver ( >= v. 0.1.19 is preferred)"
ver=0


# Check for the correct number of arguments
if [ $# -lt 7 ]
then echo "usage: $0 fastq_file barcode_ref bowtie_index_path insert cutoff name group_number [intermediate_on?]\n"
  exit 1   # General error
fi

fastq_file=$1
barcode_ref=$2
bowtie_index=$3
insert=$4
cutoff=$5
name=$6
group_number=$7
inter=$8

# Make sure the fastq file exists & contains data:

if test ! -f "$fastq_file"
then echo "$fastq_file doesn't exist!"
  exit 1
elif test ! -s "$fastq_file"
then echo "$fastq_file is empty!"
  exit 1
else
  echo "$fastq_file detected."
fi

# Set parameters for the indicated insert

echo $insert | grep -i -q -w 'mlv' && insert=mlv
echo $insert | grep -i -q -w 'tol2' && insert=tol2
echo $insert | grep -i -q -w 'hiv' && insert=hiv
echo $insert | grep -i -q -w 'fv' && insert=fv
echo $insert | grep -i -q -w 'foamy' && insert=fv

if [ "$insert" = "mlv" ]
then 
   echo "MLV integrations will be mapped."
   LTR=TGAAAGACCCCCGCTGACGGGTAGTCAATCACTC
   footprint=4
elif [ "$insert" = "tol2" ]
then
   echo "Tol2 integrations will be mapped."
   LTR=CAGAGGTGTAAAAAGTACTCAAAAATTTTACTCAAGTGA
   footprint=8
elif [ "$insert" = "hiv" ]
then
   echo "HIV integrations will be mapped."
elif [ "$insert" = "fv" ]
then
   echo "Foamy Virus integrations will be mapped."
else
   throw_error "Insert $insert not recognized."
fi

# Check if intermediate files are to be deleted or not; default is to delete.

keep=off

echo $inter | grep -i -q -w 'on' && keep=on
echo $inter | grep -i -q -w 'yes' && keep=on
echo $inter | grep -i -q -w 'y' && keep=on
echo $inter | grep -i -q -w 'true' && keep=on

if [ "$keep" = "on" ]
then 
   echo "Retention of intermediate files has been turned on."
else
   echo "Intermediate files will be deleted."
fi

# Verify that the name does not have blanks

echo $name | grep -q [[:blank:]] && throw_error "'name' can't have blanks"

# Make a working directory for the output of this script.

workdir=$PWD/Workdir_${name}_${insert}_$JOB_ID

if [ -d $workdir ] ; then throw_error "$workdir already exists!"; fi

mkdir $workdir

# Need to move the file to $workdir because at least one perl script puts its 
# output file in the same directory as the input file.

echo ""
echo "Entering working directory ${workdir}..."
cd $workdir

# Create a file of just the barcode sequences
cut -f1 ../$barcode_ref > ${barcode_ref}_seq

# Adjust for relative paths
bowtie_index=`echo $bowtie_index | awk '{ if($1 ~ /^\//){print}else{print "../"$1} }'`

# Check to make sure certain assumptions are met: for example, that the bowtie 
# index you'll need actually exists, and is in the right place.

function verify_index
{
 if 
   test ! -e $1
 then 
   cd ..
   rmdir $workdir
   throw_error "Bowtie index $1 doesn't exist!"
 elif 
   test ! -s $1
 then 
   cd ..
   rmdir $workdir
   throw_error "Bowtie index $1 is empty!"
 else
   echo "Index $1 verified."
 fi
}

echo ""
verify_index ${bowtie_index}.1.ebwt
verify_index ${bowtie_index}.2.ebwt
verify_index ${bowtie_index}.rev.1.ebwt
verify_index ${bowtie_index}.rev.2.ebwt
# .3.ebwt and .4.ebwt aren't used for single-end alignment, so it doesn't matter
# if they're there.

################################################################################
# Step 1.
# Prior to running this script, use Bamtools (or similar) to convert a BAM file 
# into FASTQ format.  Read names must end in /1 or /2.

################################################################################
# Step 2.
# Make cutadapt files for LTR_F, LTR_R, linker_F, and LTRcat (which is a 
# concatenation of LTR F and R, then reduced with 
# remove_duplicate_fastq_names_v1.0.pl).  This step also makes linker_cat, which is 
# needed later to remove inverted linker.
# "F" refers to "forward"; "R" refers to "reverse".  Both are in reference to 
# the orientation of the sequence in question (proviral LTR or linker) on a 
# given read.

echo ""
echo "***Step 2: Make initial cutadapt files.***"
# First, we want to find "forward"-oriented LTR adapters that overlap the 5' end
# of the read, OR don't overlap either end; we want to keep the DNA downstream 
# of the LTR.

# A common theme throughout the script is that cutadapt can only be used to trim
# correctly if it is searching for the reverse complement of the target 
# (otherwise, it would trim the wrong side if the target did not overlap either 
# end).  As such, there are many points at which the fastq file needs to be
# reverse-complemented, trimmed, and then reverse-complemented back to normal.

# Make the reverse complement of the fastq file

echo ""
echo "Making reverse complement of input fastq file..."
cat ../$fastq_file | ../perl/rcFastq.pl > revcom_${fastq_file}

test_file revcom_${fastq_file}

# Use cutadapt to search for the reverse complement of the LTR sequence in the 
# reverse complement of the fastq file.  As explained above, the reason this is 
# necessary is because cutadapt doesn't have a dedicated function for finding 
# reads that occur where we want AND for getting the DNA from the side
# that we want (specifically, for KEEPING the DNA 3' of the LTR).

# First, remove reads with LTR+TTTG+GGGGCTC.  These reads are false positives;
# they are amplification events off the 5' LTR into the primer binding domain of
# the provirus, and do not provide accurate integration site information.  As
# such, they are removed.  Note that this is only the case if looking for MLV.

# '-e 0' indicates that no mismatches are allowed.

if [ "$insert" = "mlv" ]
then
    echo ""
    echo "Performing LTR_F+TTTG+GGGGCTC removal with cutadapt..."
    
    cutadapt -a GAGCCCCCAAATGAAAGACCCCCGCTGACGGGTAGTCAATCACTC -q 3 -m 11 -O 17 \
      --untrimmed-output=noLTR+GGGGCTCrevcom_${fastq_file} \
      -o /dev/null revcom_${fastq_file}
   
   echo ""
   test_file noLTR+GGGGCTCrevcom_${fastq_file}
   inputF=noLTR+GGGGCTCrevcom_${fastq_file}
else
   inputF=revcom_${fastq_file}
fi

# Set quality cutoff to 3 (to ignore reads if the LTR portion is of poor 
# quality; that is, if the FASTQ quality is #)
# Require that trimmed reads have at least 11 bp left
# Require an overlap of at least 20 bp between the read and the LTR (to be 
# conservative)
# We don't really care about the untrimmed output; this is just to separate it 
# from the trimmed output.  That's why it's sent to null.

echo ""
echo "Performing LTR_F detection with cutadapt..."
cutadapt -a ${LTR} -q 3 -m 11 -O 17 -e 0 \
 --untrimmed-output=/dev/null -o cut_LTRrevcom_${fastq_file} \
 ${inputF}

echo ""
test_file cut_LTRrevcom_${fastq_file}

# If the MLV intermediate was produced, remove it.
if [ -f noLTR+GGGGCTCrevcom_${fastq_file} ]
then
   rm noLTR+GGGGCTCrevcom_${fastq_file} 
fi

echo ""
echo "Flipping cutadapt output back to the original orientation..."
cat cut_LTRrevcom_${fastq_file} | ../perl/rcFastq.pl > cut_LTR_F_${fastq_file}

test_file cut_LTR_F_${fastq_file}

rm cut_LTRrevcom_${fastq_file}



# This part is next because it also makes use of revcom_${fastq_file}; when it's
# done, it can remove it.

# Make linker_F, a file in which forward-facing linker has been detected and 
# removed.

# Require an overlap of at least 7 bp between the read and the linker (this will
# end up being the same as requiring an overlap of 14 bp: 7 are the barcode, and
# 7 are the actual linker)
# Again, we don't really care about the untrimmed output; this is just to 
# separate it from the trimmed output.  That's why it's sent to null.

echo ""
echo "Performing linker_F detection with cutadapt..."
cutadapt \
 -a ATGCGCAGTCGACCACGC -q 3 -m 18 -O 10 -e 0 \
 --untrimmed-output=/dev/null -o cut_linkerrevcom_${fastq_file} revcom_${fastq_file}

test_file cut_linkerrevcom_${fastq_file}

rm revcom_${fastq_file}

echo ""
echo "Trimming linker_F barcodes.  Barcode grab.pl says:"
perl ../perl/barcode_grab_v1.0.pl ${barcode_ref}_seq cut_linkerrevcom_${fastq_file}

test_file cut_linkerrevcom_${fastq_file}_bctrimmed.fastq

# Note that cut_linkerrevcom_${fastq_file}_barcodes.fastq is also produced here,
# the list of barcodes from reads with linker on the "forward" end.
test_file cut_linkerrevcom_${fastq_file}_barcodes.fastq

# Remove superfluous files
rm cut_linkerrevcom_${fastq_file}
rm cut_linkerrevcom_${fastq_file}_untrimmed.fastq  # Lists untrimmed reads
rm cut_linkerrevcom_${fastq_file}_bcsummary  # Summary of how many times each
# barcode was detected

echo "Flipping cutadapt output back to the original orientation..."
cat cut_linkerrevcom_${fastq_file}_bctrimmed.fastq | \
 ../perl/rcFastq.pl > cut_linker_F_${fastq_file}

test_file cut_linker_F_${fastq_file}

rm cut_linkerrevcom_${fastq_file}_bctrimmed.fastq


# Make linker_R: this will be used to pare down duplicates from linker_F, later.
# The idea is that reads that contain both forward and reverse linker sequence
# have something wrong with them, and should be removed.

echo ""
echo "Performing linker_R detection with cutadapt..."
cutadapt \
 -a ATGCGCAGTCGACCACGC -q 3 -m 18 -O 10 -e 0 \
 --untrimmed-output=/dev/null -o cut_linker_R_withbc_${fastq_file} ../${fastq_file}

# "withbc" means it hasn't had the barcodes trimmed off yet.
test_file cut_linker_R_withbc_${fastq_file}

echo ""
echo "Trimming linker_R barcodes.  Barcode grab.pl says:"
perl ../perl/barcode_grab_v1.0.pl ${barcode_ref}_seq cut_linker_R_withbc_${fastq_file}

# Change file names
mv cut_linker_R_withbc_${fastq_file}_bctrimmed.fastq cut_linker_R_${fastq_file}
mv cut_linker_R_withbc_${fastq_file}_barcodes.fastq cut_linker_R_${fastq_file}_barcodes.fastq

echo ""
test_file cut_linker_R_${fastq_file}

test_file cut_linker_R_${fastq_file}_barcodes.fastq

# Remove superfluous files
rm cut_linker_R_withbc_${fastq_file}
rm cut_linker_R_withbc_${fastq_file}_untrimmed.fastq
rm cut_linker_R_withbc_${fastq_file}_bcsummary 


echo ""
echo "Combining linker_F and linker_R..."
cat cut_linker_F_${fastq_file} cut_linker_R_${fastq_file} > ${fastq_file}_linker_cat.fastq

test_file ${fastq_file}_linker_cat.fastq

rm cut_linker_R_${fastq_file} 

# IMPORTANT - there may be duplicate names in that final output file; that's why
# I run remove_duplicate_fastq_names_v1.0.pl, to remove reads with linker_F and _R.
# Output is ${fastq_file}_linker_cat.fastq_nodup.fastq

echo ""
echo "The perl script remove_duplicate_fastq_names_v1.0.pl says:"
perl ../perl/remove_duplicate_fastq_names_v1.0.pl ${fastq_file}_linker_cat.fastq

echo ""
test_file ${fastq_file}_linker_cat.fastq_nodup.fastq

rm ${fastq_file}_linker_cat.fastq

# Now, use the nodup.fastq file to pare down linker_F, such that it does not 
# contain any reads that have linker pointed both left and right in the same 
# read.

# The first file is the list of names that you want to keep; the second file is 
# the one you get the actual sequences from.

echo ""
echo "The perl script fastq_file_overlap_v1.0.pl is running to remove inverted linker" 
echo "from cut_linker_F_${fastq_file}."
echo "It says:"
perl ../perl/fastq_file_overlap_v1.0.pl ${fastq_file}_linker_cat.fastq_nodup.fastq \
 cut_linker_F_${fastq_file} 

echo ""
test_file cut_linker_F_${fastq_file}_overlap.fastq

rm cut_linker_F_${fastq_file}

# The output is cut_linker_F_${fastq_file}_overlap.fastq, but that's a 
# misleading name, so change it.

echo "Changing the name of cut_linker_F_${fastq_file}_overlap.fastq..."
mv cut_linker_F_${fastq_file}_overlap.fastq cut_linker_F_nodup_${fastq_file}

test_file cut_linker_F_nodup_${fastq_file}



# Next, look for the reverse complement of the LTR in the original "unflipped" 
# fastq file.  This makes the LTR_R file.

# As before, first remove reads with LTR_R+TTTG+GGGGCTC.

if [ "$insert" = "mlv" ]
then
   echo ""
   echo "Performing LTR_R+TTTG+GGGGCTC removal with cutadapt..."
   cutadapt -a GAGCCCCCAAATGAAAGACCCCCGCTGACGGGTAGTCAATCACTC -q 3 -m 11 -O 17 \
     --untrimmed-output=noLTR_R-GGGGCTC_${fastq_file} \
     -o /dev/null ../${fastq_file}
   
   test_file noLTR_R-GGGGCTC_${fastq_file}
   
   inputR=noLTR_R-GGGGCTC_${fastq_file}
else
   inputR=../${fastq_file}
fi

# Detect/trim actual LTR.

echo ""
echo "Performing LTR_R detection with cutadapt..."
cutadapt -a ${LTR} -q 3 -m 11 -O 17 -e 0 \
 --untrimmed-output=/dev/null -o cut_LTR_R_${fastq_file} \
 ${inputR}

test_file cut_LTR_R_${fastq_file}

if [ -f noLTR_R-GGGGCTC_${fastq_file} ]
then
   rm noLTR_R-GGGGCTC_${fastq_file}
fi

# Concatenate LTR F and R together

echo ""
echo "Combining LTR_F and LTR_R..."
cat cut_LTR_F_${fastq_file} cut_LTR_R_${fastq_file} > ${fastq_file}_LTR_cat.fastq

# IMPORTANT - there may be duplicate names in that final output file (reads that
# have both forward and reverse LTR); run remove_duplicate_fastq_names_v1.0.pl to 
# get rid of them.

echo ""
echo "Removing reads with both LTR_F and LTR_R."
echo "The perl script remove_duplicate_fastq_names_v1.0.pl says:"
perl ../perl/remove_duplicate_fastq_names_v1.0.pl ${fastq_file}_LTR_cat.fastq

echo ""
test_file ${fastq_file}_LTR_cat.fastq_nodup.fastq

rm ${fastq_file}_LTR_cat.fastq


# The only final output to all this should be cut_LTR_F_${fastq_file}, 
# cut_LTR_R_${fastq_file}, cut_linker_F_nodup_${fastq_file}, 
# ${fastq_file}_LTR_cat.fastq_nodup.fastq, 
# ${fastq_file}_linker_cat.fastq_nodup.fastq, 
# cut_linkerrevcom_${fastq_file}_barcodes.fastq, and 
# cut_linker_R_${fastq_file}_barcodes.fastq


################################################################################
# Step 3.
# Use split_LTR_to_long_and_short_v1.0.pl to look through LTR_F, and sort the reads 
# into long and short fragments (2 hashes). Short fragments are those in which 
# LTR could be detected by cutadapt in forward orientation in one read, and 
# reverse orientation on the paired read; long fragments are those in which only
# forward LTR could be detected.

echo ""
echo "***Step 3. Sort LTR_F into long and short fragments.***"

echo ""
echo "split_LTR_to_long_and_short_v1.0.pl says:"
perl ../perl/split_LTR_to_long_and_short_v1.0.pl cut_LTR_R_${fastq_file} \
 cut_LTR_F_${fastq_file} ${fastq_file}_LTR_cat.fastq_nodup.fastq

echo ""
test_file cut_LTR_F_${fastq_file}_short.fastq

test_file cut_LTR_F_${fastq_file}_long.fastq

rm cut_LTR_F_${fastq_file}

############################ Formatting LTR-F-long #############################
# Step 4.
# Use cutadapt to find forward linker in the LTR_F_long file.  If the read has 
# F linker in the same read that has F LTR, throw it out.

echo ""
echo "***Step 4. Begin formatting of LTR_F_long.***"

# This processes trims reads that have been designated as coming from long 
# fragments.

# First, remove any reads with forward-facing linker.
# In this case, I don't care what the barcode is. Further, it's MORE
# conservative to not set -e to 0, as you'll throw out more reads this way;
# the same goes for not setting -q.

echo ""
echo "Creating reverse complement of LTR_F_long..."
cat cut_LTR_F_${fastq_file}_long.fastq | ../perl/rcFastq.pl > \
 revcom_cut_LTR_F_${fastq_file}_long.fastq

test_file revcom_cut_LTR_F_${fastq_file}_long.fastq

echo ""
echo "Using cutadapt to remove reads containing linker_F..."
cutadapt -a ATGCGCAGTCGACCACGC -m 11 -O 17 \
 --untrimmed-output=noncut_revcom_cut_LTR_F_${fastq_file}_long.fastq \
 -o /dev/null revcom_cut_LTR_F_${fastq_file}_long.fastq

test_file noncut_revcom_cut_LTR_F_${fastq_file}_long.fastq

rm revcom_cut_LTR_F_${fastq_file}_long.fastq

echo "Flipping the cutadapt untrimmed output back to the original orientation..."
cat noncut_revcom_cut_LTR_F_${fastq_file}_long.fastq | ../perl/rcFastq.pl > \
 cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link

test_file cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link

rm noncut_revcom_cut_LTR_F_${fastq_file}_long.fastq



################################################################################
# Step 5.
# Retain only reads that have forward linker in the paired read.  This ensures
# that all long fragments will have LTR and linker that are oriented toward each
# other.

echo ""
echo "***Step 5. Only keep LTR_F_long reads that have linker_F in the paired read.***"

# Remove reads that do not have a paired read with forward-facing linker.  This 
# is the reason that the linker_F file needs to be included in this portion of 
# the script.

echo ""
echo "Running pair_in_other_file_v1.0.pl to print LTR_F_long reads that have linker_F"
echo "in their paired read.  It says:"
perl ../perl/pair_in_other_file_v1.0.pl \
 cut_linker_F_nodup_${fastq_file} cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link

# Output is cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink.fastq.

echo ""
test_file cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink.fastq

# Once I have _no-bad-F-link_pairhaslink.fastq, there's no need to keep 
# _no-bad-F-link around.
rm cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link

# if [ "$keep" = "off" ]; then rm cut_linker_F_nodup_${fastq_file}; fi

################################################################################
# Step 6.
# Use cutadapt to remove reads with inappropriate reverse linker sequences.  
# Note that reads that lack linker are NOT thrown out (--untrimmed is set, but 
# those results are later concatenated with the others).

echo ""
echo "***Step 6. Trim linker_R from LTR_F_long reads.***"

# Use cutadapt to remove reverse linker from the output of the perl script.  So,
# the final output should be the LTR-containing reads of long fragments, in 
# which any reverse linker has been trimmed, and there never was any forward 
# linker (on that read) in the first place.  The idea is that reverse linker IS 
# permissible on the read with forward LTR, but is not required.

# I use 'uo' to refer to untrimmed output of cutadapt.

echo ""
echo "Trimming linker_R using cutadapt..."
cutadapt \
 -a ATGCGCAGTCGACCACGC -q 3 -m 11 -O 1 -e 0 \
 --untrimmed-output=cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink_uo.fastq \
 -o cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink_rm-R-link_withbc.fastq \
 cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink.fastq

test_file cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink_rm-R-link_withbc.fastq
test_file cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink_uo.fastq

## Again, now that I have _no-bad-F-link_pairhaslink_rm-R-link.fastq, there's no 
## need to keep the previous file.
#rm cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink.fastq

# Make a list of the reads that look like they might have linker_R
echo "Extracting names that MIGHT have linker_R..."
awk 'NR%4==1' \
  cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink_rm-R-link_withbc.fastq \
  | sort \
  > cut_LTR_F_${fastq_file}_long.fastq_putative_linker_R

test_file cut_LTR_F_${fastq_file}_long.fastq_putative_linker_R


# Remove the barcodes; the relevant barcodes have already been recorded, so 
# these don't need to be kept.
echo ""
echo "Trimming linker_R barcodes.  Barcode grab.pl says:"
perl ../perl/barcode_grab_v1.0.pl ${barcode_ref}_seq \
 cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink_rm-R-link_withbc.fastq

test_file cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink_rm-R-link_withbc.fastq_bctrimmed.fastq

# Delete superfluous files
rm cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink_rm-R-link_withbc.fastq
rm cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink_rm-R-link_withbc.fastq_barcodes.fastq
rm cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink_rm-R-link_withbc.fastq_untrimmed.fastq
rm cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink_rm-R-link_withbc.fastq_bcsummary


# Extract the names of the putative linker_R reads that DID have barcode
echo "Extracting names of reads that DID have barcode..."
awk 'NR%4==1' \
  cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink_rm-R-link_withbc.fastq_bctrimmed.fastq \
  | sort \
  > cut_LTR_F_${fastq_file}_long.fastq_putative_linker_R_withbc

test_file cut_LTR_F_${fastq_file}_long.fastq_putative_linker_R_withbc


# Use join -v to identify the reads that did NOT actually have barcode; restore
# them to their pre-trimmed state.
echo ""
echo "Identifying putative linker_R reads that did NOT have barcode..."
join -v 1 \
  cut_LTR_F_${fastq_file}_long.fastq_putative_linker_R \
  cut_LTR_F_${fastq_file}_long.fastq_putative_linker_R_withbc \
  > cut_LTR_F_${fastq_file}_long.fastq_no_linker_R_barcode

test_file cut_LTR_F_${fastq_file}_long.fastq_no_linker_R_barcode

rm cut_LTR_F_${fastq_file}_long.fastq_putative_linker_R
rm cut_LTR_F_${fastq_file}_long.fastq_putative_linker_R_withbc


echo "Restoring non-linker_R reads to pre-trimmed state..."
perl ../perl/names_to_fastq.pl \
  cut_LTR_F_${fastq_file}_long.fastq_no_linker_R_barcode \
  cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink.fastq

test_file cut_LTR_F_${fastq_file}_long.fastq_no_linker_R_barcode.fastq

rm cut_LTR_F_${fastq_file}_long.fastq_no_linker_R_barcode
rm cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink.fastq

echo "Trimming the # bases from the pre-trimmed reads..."
cat cut_LTR_F_${fastq_file}_long.fastq_no_linker_R_barcode.fastq \
  | perl ../perl/fastq_lowqual_trim.pl \
  > cut_LTR_F_${fastq_file}_long.fastq_no_linker_R_barcode_qualtrimmed.fastq

test_file cut_LTR_F_${fastq_file}_long.fastq_no_linker_R_barcode_qualtrimmed.fastq

rm cut_LTR_F_${fastq_file}_long.fastq_no_linker_R_barcode.fastq


# Combine the reads without linker to those that had linker and barcode trimmed.
# Remember that "without linker" encompasses two files: those in which linker
# was never seen, and those where the linker was thought to be there, but 
# wasn't.
cat cut_LTR_F_${fastq_file}_long.fastq_no_linker_R_barcode_qualtrimmed.fastq \
 cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink_uo.fastq \
 cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink_rm-R-link_withbc.fastq_bctrimmed.fastq \
 > cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink_rm-R-link.fastq

test_file cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink_rm-R-link.fastq

rm cut_LTR_F_${fastq_file}_long.fastq_no_linker_R_barcode_qualtrimmed.fastq
rm cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink_uo.fastq
rm cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink_rm-R-link_withbc.fastq_bctrimmed.fastq

# Also, since I'm done processing LTR_F_long, I no longer need 
# cut_LTR_F_${fastq_file}_long.fastq.
rm cut_LTR_F_${fastq_file}_long.fastq

########################### Formatting LTR-F-short #############################
# Step 7.
# Remove reads with inverted linker from LTR-F-short.  Then, use cutadapt to 
# trim reverse linker from LTR_F reads from short fragments; only keep those 
# that had linker.

echo ""
echo "***Step 7. Trim linker_R from LTR_F_short.***"

# This step trims reverse linker from reads that were judged to be from short 
# fragments by split_LTR_to_long_and_short_v1.0.pl.

# First, remove any reads with inverted linker (forward linker AND reverse 
# linker).  Any truly good read from LTR_F_short will also be in 
# linker_cat (which already had inverted linkers removed), albeit trimmed 
# differently.  The perl script only looks at names, though, so the trimming 
# doesn't matter for this step.

echo ""
echo "Running fastq_file_overlap_v1.0.pl to only keep LTR_F_short reads that"
echo "do NOT have inverted linker.  It says:"
perl ../perl/fastq_file_overlap_v1.0.pl ${fastq_file}_linker_cat.fastq_nodup.fastq \
 cut_LTR_F_${fastq_file}_short.fastq

# Output is cut_LTR_F_${fastq_file}_short.fastq_overlap.fastq.

test_file cut_LTR_F_${fastq_file}_short.fastq_overlap.fastq

if [ "$keep" = "off" ]; then rm cut_LTR_F_${fastq_file}_short.fastq; fi

# Don't remove $linker_cat yet - I use it for step 8, too.
# Actually, no, it's better to use the linker_F file for that.

rm ${fastq_file}_linker_cat.fastq_nodup.fastq

# Trim reads of reverse linker.  As all short reads should have reverse linker, 
# --untrimmed is sent to null (meaning the only reads to be kept are those with
# reverse linker).  Low quality may obscure some reads that genuinely do have 
# linker, but for those, we just have to hope that the paired read is caught in 
# another step (for LTR_R reads).

echo ""
echo "Trimming linker_R off of LTR_F_short via cutadapt..."
cutadapt -a ATGCGCAGTCGACCACGC -q 3 -m 18 -O 1 -e 0 \
 --untrimmed-output=/dev/null \
 -o cut_LTR_F_${fastq_file}_short.fastq_trim-R-link_withbc.fastq \
 cut_LTR_F_${fastq_file}_short.fastq_overlap.fastq

echo ""
echo "Trimming linker_R barcodes.  Barcode grab.pl says:"
perl ../perl/barcode_grab_v1.0.pl ${barcode_ref}_seq cut_LTR_F_${fastq_file}_short.fastq_trim-R-link_withbc.fastq

test_file cut_LTR_F_${fastq_file}_short.fastq_trim-R-link_withbc.fastq_bctrimmed.fastq

rm cut_LTR_F_${fastq_file}_short.fastq_overlap.fastq

# Remove superfluous files
rm cut_LTR_F_${fastq_file}_short.fastq_trim-R-link_withbc.fastq_barcodes.fastq
rm cut_LTR_F_${fastq_file}_short.fastq_trim-R-link_withbc.fastq
rm cut_LTR_F_${fastq_file}_short.fastq_trim-R-link_withbc.fastq_untrimmed.fastq
rm cut_LTR_F_${fastq_file}_short.fastq_trim-R-link_withbc.fastq_bcsummary

# Change the name of the output file
mv cut_LTR_F_${fastq_file}_short.fastq_trim-R-link_withbc.fastq_bctrimmed.fastq \
 cut_LTR_F_${fastq_file}_short.fastq_trim-R-link.fastq
 
test_file cut_LTR_F_${fastq_file}_short.fastq_trim-R-link.fastq


############################### Formatting LTR-R ###############################
# Step 8.
# Use fastq_file_overlap twice - once to clean LTR-R of inverted LTR reads, and 
# once to clean it of inverted linker reads.  Then, use cutadapt on LTR-R, only 
# keeping reads that contain forward-facing linker.

echo ""
echo "***Step 8. Format LTR_R.***"

# Use fastq_file_overlap_v1.0.pl to eliminate entries in LTR_R that have inverted 
# LTR.  This IS a little different from split_LTR_to_long_and_short_v1.0.pl, in which
# the sequences were printed from LTR_cat; here, they're printed from LTR_R.  
# The result is the same, though.

echo ""
echo "Running fastq_file_overlap_v1.0.pl to remove entries from LTR_R that have inverted LTR."
echo "It says:"
perl ../perl/fastq_file_overlap_v1.0.pl ${fastq_file}_LTR_cat.fastq_nodup.fastq \
 cut_LTR_R_${fastq_file}

# Output is cut_LTR_R_${fastq_file}_overlap.fastq

echo ""
test_file cut_LTR_R_${fastq_file}_overlap.fastq

if [ "$keep" = "off" ]; then rm cut_LTR_R_${fastq_file}; fi
if [ "$keep" = "off" ]; then rm ${fastq_file}_LTR_cat.fastq_nodup.fastq; fi

# We also need to remove any reads that have BOTH forward and reverse linker.  
# Do so by comparing ${LTR_R}_overlap.fastq to the linker_cat file.  At this 
# point, we only want reads that have forward-facing linker anyway, so the good 
# reads won't be thrown out by doing this.

echo ""
echo "Running fastq_file_overlap_v1.0.pl to remove entries from LTR_R that have inverted linker."
echo "It says:"
# perl ../perl/fastq_file_overlap_v1.0.pl ${fastq_file}_linker_cat.fastq_nodup.fastq \
#  cut_LTR_R_${fastq_file}_overlap.fastq

perl ../perl/fastq_file_overlap_v1.0.pl cut_linker_F_nodup_${fastq_file} \
 cut_LTR_R_${fastq_file}_overlap.fastq

# Output is cut_LTR_R_${fastq_file}_overlap.fastq_overlap.fastq.

echo ""
test_file cut_LTR_R_${fastq_file}_overlap.fastq_overlap.fastq

rm cut_LTR_R_${fastq_file}_overlap.fastq

#rm ${fastq_file}_linker_cat.fastq_nodup.fastq

if [ "$keep" = "off" ]; then rm cut_linker_F_nodup_${fastq_file}; fi

# Then, use cutadapt to identify reads with linker_F on the same read.  Trim 
# these reads of linker, and keep only those reads.

# To do that, the fastq file must first be flipped.

echo ""
echo "Creating reverse complement of cut_LTR_R_${fastq_file}_overlap.fastq_overlap.fastq."
cat cut_LTR_R_${fastq_file}_overlap.fastq_overlap.fastq | ../perl/rcFastq.pl > \
 revcom_cut_LTR_R_${fastq_file}_doubleoverlap.fastq

test_file revcom_cut_LTR_R_${fastq_file}_doubleoverlap.fastq

rm cut_LTR_R_${fastq_file}_overlap.fastq_overlap.fastq

# To detect forward linker and trim it correctly, it's necessary to look for 
# reverse complement linker on reverse complement template.

echo ""
echo "Detecting linker_F on LTR_R reads (by finding linker_R on reverse"
echo "complement LTR_R reads)..."
cutadapt -a ATGCGCAGTCGACCACGC -q 3 -m 18 -O 10 -e 0 \
 --untrimmed-output=/dev/null \
 -o cut_revcom_cut_LTR_R_${fastq_file}_doubleoverlap_withbc.fastq \
 revcom_cut_LTR_R_${fastq_file}_doubleoverlap.fastq

test_file cut_revcom_cut_LTR_R_${fastq_file}_doubleoverlap_withbc.fastq

echo ""
echo "Trimming linker_F barcodes.  Barcode grab.pl says:"
perl ../perl/barcode_grab_v1.0.pl ${barcode_ref}_seq cut_revcom_cut_LTR_R_${fastq_file}_doubleoverlap_withbc.fastq

test_file cut_revcom_cut_LTR_R_${fastq_file}_doubleoverlap_withbc.fastq_bctrimmed.fastq

# Remove superfluous files
rm cut_revcom_cut_LTR_R_${fastq_file}_doubleoverlap_withbc.fastq_barcodes.fastq
rm cut_revcom_cut_LTR_R_${fastq_file}_doubleoverlap_withbc.fastq
rm cut_revcom_cut_LTR_R_${fastq_file}_doubleoverlap_withbc.fastq_untrimmed.fastq
rm cut_revcom_cut_LTR_R_${fastq_file}_doubleoverlap_withbc.fastq_bcsummary 

# Change the name of the trimmed file
mv cut_revcom_cut_LTR_R_${fastq_file}_doubleoverlap_withbc.fastq_bctrimmed.fastq \
 cut_revcom_cut_LTR_R_${fastq_file}_doubleoverlap.fastq


test_file cut_revcom_cut_LTR_R_${fastq_file}_doubleoverlap.fastq

rm revcom_cut_LTR_R_${fastq_file}_doubleoverlap.fastq

echo ""
echo "Flipping the cutadapt output back to the original orientation..."
cat cut_revcom_cut_LTR_R_${fastq_file}_doubleoverlap.fastq | ../perl/rcFastq.pl > \
 cut_LTR_R_${fastq_file}_overlap_trim-F-link.fastq

test_file cut_LTR_R_${fastq_file}_overlap_trim-F-link.fastq

rm cut_revcom_cut_LTR_R_${fastq_file}_doubleoverlap.fastq


########################### Begin compare-p ####################################
# Step 9.
# Make the |p pairs of the LTR_F reads.

# The idea of "|p pairs" is to take the untrimmed paired reads of all reads
# collected so far, and trim them WITHOUT automatically removing low-quality
# bases.  The reason this is useful is in cases in which the sequence of two
# paired reads are consistent with one another, but one of them has low 
# base-call quality scores that would otherwise lead to the fragment being 
# thrown out.  These "|p" reads (so named because of the |p I attach to the end
# of the read name in the fastq file) are ONLY kept long enough to be evaluated
# in the context of the normally-trimmed reads; |p reads are not reported in the
# end.

# This assumes that LTR_F_long, _short, and LTR_R were made with -q turned on.
# For these steps, -e should NOT be set to 0, since leaving it as 0.1 (the 
# default) is more conservative when removing reads that contain a given bad 
# sequence, and helps obtain a better trim when trimming reads that may or may 
# not have the sequence in question.

echo ""
echo "***Step 9. Make non-q pairs of LTR_F reads.***"

# "non-q" refers to not using the "-q 3" option in cutadapt, which trims 
# low-quality bases.

# Start with LTR_F long and short.
# Combine LTR_F_long and LTR_F_short into a single file. 

echo ""
echo "Concatenating LTR_F long and short..."
cat cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink_rm-R-link.fastq \
 cut_LTR_F_${fastq_file}_short.fastq_trim-R-link.fastq > ${fastq_file}_LTR_F_combo

# Grab all pairs, and but DO NOT YET add '|p' to the end of the read names.
# The reason for not changing the NAME of the read, even though it may be 
# trimmed differently, has to do with Bowtie.  In Bowtie version 0.12.7 
# (presumably, other versions, too), pseudo-random activity is sometimes called 
# upon in the alignment, and the name of the read is used to seed the 
# pseudo-randomness.  That means that reads with identical sequence and quality 
# scores, but different names, can have different alignment results.  To avoid 
# this, I keep the names constant until after the alignment.  The normal reads 
# can be distinguished from the |p reads by filename.

echo ""
echo "Grabbing pairs of LTR_F. fastq_get_pair_v1.0.pl says:"
perl ../perl/fastq_get_pair_v1.0.pl ${fastq_file}_LTR_F_combo \
 ../$fastq_file

test_file ${fastq_file}_LTR_F_combo_pairs-with-p.fastq

rm ${fastq_file}_LTR_F_combo



# Use cutadapt to remove any reads with linker_R (essentially to get rid of 
# reads with inverted linker; since the reads come from the original file, that
# would not have been filtered out).  Again, no -q.

echo ""
echo "Performing removal of reads that have linker_R with cutadapt..."
cutadapt -a ATGCGCAGTCGACCACGC -m 11 -O 10 \
 --untrimmed-output=${fastq_file}_LTR_F_combo_pairs-with-p_nolinkR.fastq \
  -o /dev/null ${fastq_file}_LTR_F_combo_pairs-with-p.fastq

# Note that, since this is only to throw out certain reads, there's no need
# to call barcode_grab_v1.0.pl.

test_file ${fastq_file}_LTR_F_combo_pairs-with-p_nolinkR.fastq

rm ${fastq_file}_LTR_F_combo_pairs-with-p.fastq

# Trim LTR_R from the reads.  If a read does not have LTR_R, allow it to pass
# through, anyway (i.e., don't set --untrimmed).  This is because some of these
# reads (all of which are pairs of LTR_F, remember) come from long fragments, 
# and would not be expected to have LTR_R.

# First, remove LTR_R+TTTG+GGGGCTC.

if [ "$insert" = "mlv" ]
then
   echo ""
   echo "Performing LTR_R+TTTG+GGGGCTC detection/removal with cutadapt..."
   cutadapt -a GAGCCCCCAAATGAAAGACCCCCGCTGACGGGTAGTCAATCACTC -m 11 -O 17 \
     --untrimmed-output=${fastq_file}_LTR_F_combo_pairs-with-p_nolinkR_noGGGGCTC.fastq \
     -o /dev/null ${fastq_file}_LTR_F_combo_pairs-with-p_nolinkR.fastq
   
   test_file ${fastq_file}_LTR_F_combo_pairs-with-p_nolinkR_noGGGGCTC.fastq
   
   pair_inputR=${fastq_file}_LTR_F_combo_pairs-with-p_nolinkR_noGGGGCTC.fastq
   
   #rm ${fastq_file}_LTR_F_combo_pairs-with-p_nolinkR.fastq
else
   pair_inputR=${fastq_file}_LTR_F_combo_pairs-with-p_nolinkR.fastq
fi

# I avoided setting -e to 0; everything makes it through anyway, so not setting
# it just means I'll get more accurate trims.

# The overlap values is unusually low here because it needs to accomodate paired
# reads that had as little as 8 bp of linker_R on the other read (7 bp barcode +
# 1 bp linker).

echo ""
echo "Performing LTR_R detection/trim with cutadapt..."
cutadapt -a ${LTR} -m 11 -O 5 \
 -o ${fastq_file}_LTR_F_combo_pairs-with-p_nolinkR_trimLTR-R.fastq \
 ${pair_inputR}

test_file ${fastq_file}_LTR_F_combo_pairs-with-p_nolinkR_trimLTR-R.fastq

rm ${fastq_file}_LTR_F_combo_pairs-with-p_nolinkR.fastq
# rm ${fastq_file}_LTR_F_combo_pairs-with-p_nolinkR_trimLTR-R+TTTG.fastq
if [ -f ${fastq_file}_LTR_F_combo_pairs-with-p_nolinkR_noGGGGCTC.fastq ]
then
   rm ${fastq_file}_LTR_F_combo_pairs-with-p_nolinkR_noGGGGCTC.fastq
fi

# Use cutadapt to remove any read with LTR_F;  do not use -q.
# I've lowered the -O (required overlap) from 20 to 14 to err on the side of 
# caution.

echo ""
echo "Making reverse complement of LTR_F |p fastq file..."
cat ${fastq_file}_LTR_F_combo_pairs-with-p_nolinkR_trimLTR-R.fastq | \
 ../perl/rcFastq.pl > \
 revcom_${fastq_file}_LTR_F_combo_pairs-with-p_nolinkR_trimLTR-R.fastq

test_file revcom_${fastq_file}_LTR_F_combo_pairs-with-p_nolinkR_trimLTR-R.fastq

rm ${fastq_file}_LTR_F_combo_pairs-with-p_nolinkR_trimLTR-R.fastq

# Remember, if sending the real output to --untrimmed, not setting -e to 0 is
# MORE conservative.

echo ""
echo "Performing removal of reads with LTR_F with cutadapt..."
cutadapt -a ${LTR} -m 11 -O 17 \
--untrimmed-output=revcom_${fastq_file}_LTR_F_combo_pairs-with-p_nolinkR_trimLTR-R_noLTRF.fastq \
 -o /dev/null \
 revcom_${fastq_file}_LTR_F_combo_pairs-with-p_nolinkR_trimLTR-R.fastq

test_file revcom_${fastq_file}_LTR_F_combo_pairs-with-p_nolinkR_trimLTR-R_noLTRF.fastq

rm revcom_${fastq_file}_LTR_F_combo_pairs-with-p_nolinkR_trimLTR-R.fastq


# Detect and trim linker_F from those reads WITHOUT -q.  Do use -m, though.
# For example, -m 11 ensures that any trimmed reads will have at least 11 bp 
# left.  Since linker trims are typically followed by a call to bowtie_grab.pl, 
# which removes another 7 bp without allowing mismatches, -m 18 is tantamount
# to -m 11 (but only if the linker is REQUIRED).

echo ""
echo "Performing linker_F detection/trim with cutadapt..."
cutadapt -a ATGCGCAGTCGACCACGC -m 11 -O 7 -e 0 \
 --untrimmed-output=revcom_${fastq_file}_LTR_F_combo_p_trim_uo.fastq \
 -o revcom_${fastq_file}_LTR_F_combo_p_trim_withbc.fastq \
 revcom_${fastq_file}_LTR_F_combo_pairs-with-p_nolinkR_trimLTR-R_noLTRF.fastq

test_file revcom_${fastq_file}_LTR_F_combo_p_trim_withbc.fastq
test_file revcom_${fastq_file}_LTR_F_combo_p_trim_uo.fastq

rm revcom_${fastq_file}_LTR_F_combo_pairs-with-p_nolinkR_trimLTR-R_noLTRF.fastq

echo ""
echo "Trimming linker_F barcodes.  Barcode grab.pl says:"
perl ../perl/barcode_grab_v1.0.pl ${barcode_ref}_seq \
 revcom_${fastq_file}_LTR_F_combo_p_trim_withbc.fastq

test_file revcom_${fastq_file}_LTR_F_combo_p_trim_withbc.fastq_bctrimmed.fastq

# Remove superfluous files
rm revcom_${fastq_file}_LTR_F_combo_p_trim_withbc.fastq
rm revcom_${fastq_file}_LTR_F_combo_p_trim_withbc.fastq_untrimmed.fastq
rm revcom_${fastq_file}_LTR_F_combo_p_trim_withbc.fastq_bcsummary

# Concatenate the trimmed and untrimmed files
cat revcom_${fastq_file}_LTR_F_combo_p_trim_uo.fastq \
 revcom_${fastq_file}_LTR_F_combo_p_trim_withbc.fastq_bctrimmed.fastq \
 > revcom_${fastq_file}_LTR_F_combo_p_trim.fastq

rm revcom_${fastq_file}_LTR_F_combo_p_trim_uo.fastq
rm revcom_${fastq_file}_LTR_F_combo_p_trim_withbc.fastq_bctrimmed.fastq


test_file revcom_${fastq_file}_LTR_F_combo_p_trim.fastq

echo ""
echo "Flipping cutadapt output back to the original orientation..."
cat revcom_${fastq_file}_LTR_F_combo_p_trim.fastq | ../perl/rcFastq.pl \
 > ${fastq_file}_LTR_F_combo_p_trim.fastq

test_file ${fastq_file}_LTR_F_combo_p_trim.fastq

rm revcom_${fastq_file}_LTR_F_combo_p_trim.fastq


################################################################################
# Step 10.
# Making p-pairs of LTR_R reads.

echo ""
echo "***Step 10. Make non-q pairs of LTR_R reads.***"

# Next, take the LTR_R input.
# Grab all the pairs, but DO NOT YET add '|p' to the end.

echo ""
echo "Grabbing pairs of LTR_R. fastq_get_pair_v1.0.pl says:"
perl ../perl/fastq_get_pair_v1.0.pl \
 cut_LTR_R_${fastq_file}_overlap_trim-F-link.fastq ../$fastq_file

test_file cut_LTR_R_${fastq_file}_overlap_trim-F-link.fastq_pairs-with-p.fastq

# Throw out any reads with LTR_R. Use cutadapt without -q.

echo ""
echo "Performing removal of p-reads with LTR_R with cutadapt..."
cutadapt -a ${LTR} -m 11 -O 17 \
 --untrimmed-output=${fastq_file}_LTR_R_pairs-with-p_noLTR-R.fastq \
 -o /dev/null \
  cut_LTR_R_${fastq_file}_overlap_trim-F-link.fastq_pairs-with-p.fastq

test_file ${fastq_file}_LTR_R_pairs-with-p_noLTR-R.fastq

rm cut_LTR_R_${fastq_file}_overlap_trim-F-link.fastq_pairs-with-p.fastq

# Throw out any reads with linker_F (again, removing inverted linker)

echo ""
echo "Making reverse complement of LTR_R |p fastq file..."
cat ${fastq_file}_LTR_R_pairs-with-p_noLTR-R.fastq | ../perl/rcFastq.pl > \
	revcom_${fastq_file}_LTR_R_pairs-with-p_noLTR-R.fastq

test_file revcom_${fastq_file}_LTR_R_pairs-with-p_noLTR-R.fastq

rm ${fastq_file}_LTR_R_pairs-with-p_noLTR-R.fastq

# Since this is to remove reads, there's no need for barcode_grab_v1.0.pl, or to set
# -e 0.

echo ""
echo "Performing removal of p-reads with linker_F with cutadapt..."
cutadapt \
 -a ATGCGCAGTCGACCACGC -m 11 -O 10 \
 --untrimmed-output=revcom_${fastq_file}_LTR_R_pairs-with-p_noLTR-R_nolinkF.fastq \
 -o /dev/null revcom_${fastq_file}_LTR_R_pairs-with-p_noLTR-R.fastq

test_file revcom_${fastq_file}_LTR_R_pairs-with-p_noLTR-R_nolinkF.fastq

rm revcom_${fastq_file}_LTR_R_pairs-with-p_noLTR-R.fastq

# Trim the remaining reads for LTR_F.

# First, remove reads with LTR_F+TTTG+GGGGCTC.

if [ "$insert" = "mlv" ]
then
   echo ""
   echo "Performing removal of LTR_F+TTTG+GGGGCTC with cutadapt..."
   cutadapt -a GAGCCCCCAAATGAAAGACCCCCGCTGACGGGTAGTCAATCACTC -m 11 -O 17 \
     --untrimmed-output=revcom_${fastq_file}_LTR_R_pairs-with-p_noLTR-R_nolinkF_noGGGGCTC.fastq \
     -o /dev/null revcom_${fastq_file}_LTR_R_pairs-with-p_noLTR-R_nolinkF.fastq
   
   test_file revcom_${fastq_file}_LTR_R_pairs-with-p_noLTR-R_nolinkF_noGGGGCTC.fastq
   
   pair_inputF=revcom_${fastq_file}_LTR_R_pairs-with-p_noLTR-R_nolinkF_noGGGGCTC.fastq
   
   # rm revcom_${fastq_file}_LTR_R_pairs-with-p_noLTR-R_nolinkF.fastq
else
   pair_inputF=revcom_${fastq_file}_LTR_R_pairs-with-p_noLTR-R_nolinkF.fastq
fi

echo ""
echo "Performing trim of LTR_F with cutadapt..."
cutadapt -a ${LTR} -m 11 -O 14 \
 -o revcom_${fastq_file}_LTR_R_pairs-with-p_noLTR-R_nolinkF_trimLTRF.fastq \
 ${pair_inputF}

test_file revcom_${fastq_file}_LTR_R_pairs-with-p_noLTR-R_nolinkF_trimLTRF.fastq

rm revcom_${fastq_file}_LTR_R_pairs-with-p_noLTR-R_nolinkF.fastq

if [ -f revcom_${fastq_file}_LTR_R_pairs-with-p_noLTR-R_nolinkF_noGGGGCTC.fastq ]
then
   rm revcom_${fastq_file}_LTR_R_pairs-with-p_noLTR-R_nolinkF_noGGGGCTC.fastq
fi

echo ""
echo "Flipping cutadapt output back to the original orientation..."
cat revcom_${fastq_file}_LTR_R_pairs-with-p_noLTR-R_nolinkF_trimLTRF.fastq \
 | ../perl/rcFastq.pl \
 > ${fastq_file}_LTR_R_pairs-with-p_noLTR-R_nolinkF_trimLTRF.fastq

rm revcom_${fastq_file}_LTR_R_pairs-with-p_noLTR-R_nolinkF_trimLTRF.fastq

# Detect/trim linker_R on these reads.  

echo ""
echo "Performing linker_R detection/trim with cutadapt..."
cutadapt \
 -a ATGCGCAGTCGACCACGC -m 11 -O 1 -e 0 \
 --untrimmed-output=${fastq_file}_LTR_R_p_trim_uo.fastq \
 -o ${fastq_file}_LTR_R_p_trim_withbc.fastq \
 ${fastq_file}_LTR_R_pairs-with-p_noLTR-R_nolinkF_trimLTRF.fastq

test_file ${fastq_file}_LTR_R_p_trim_withbc.fastq
test_file ${fastq_file}_LTR_R_p_trim_uo.fastq


echo "Extracting names that MIGHT have linker_R..."
awk 'NR%4==1' \
  ${fastq_file}_LTR_R_p_trim_withbc.fastq \
  | sort \
  > ${fastq_file}_LTR_R_p_trim_putative_linker_R
  
test_file ${fastq_file}_LTR_R_p_trim_putative_linker_R


echo ""
echo "Trimming linker_R barcodes.  Barcode grab.pl says:"
perl ../perl/barcode_grab_v1.0.pl ${barcode_ref}_seq \
 ${fastq_file}_LTR_R_p_trim_withbc.fastq

test_file ${fastq_file}_LTR_R_p_trim_withbc.fastq_bctrimmed.fastq

# Delete superfluous files
rm ${fastq_file}_LTR_R_p_trim_withbc.fastq
rm ${fastq_file}_LTR_R_p_trim_withbc.fastq_untrimmed.fastq
rm ${fastq_file}_LTR_R_p_trim_withbc.fastq_bcsummary

# Extract the names of the putative linker_R reads that DID have barcode
echo "Extracting names of reads that DID have barcode..."
awk 'NR%4==1' \
  ${fastq_file}_LTR_R_p_trim_withbc.fastq_bctrimmed.fastq \
  | sort \
  > ${fastq_file}_LTR_R_p_trim_putative_linker_R_withbc

test_file ${fastq_file}_LTR_R_p_trim_putative_linker_R_withbc


# Use join -v to identify the reads that did NOT actually have barcode; restore
# them to their pre-trimmed state.
echo ""
echo "Identifying putative linker_R reads that did NOT have barcode..."
join -v 1 \
  ${fastq_file}_LTR_R_p_trim_putative_linker_R \
  ${fastq_file}_LTR_R_p_trim_putative_linker_R_withbc \
  > ${fastq_file}_LTR_R_p_trim_no_linker_R_barcode

test_file ${fastq_file}_LTR_R_p_trim_no_linker_R_barcode

rm ${fastq_file}_LTR_R_p_trim_putative_linker_R
rm ${fastq_file}_LTR_R_p_trim_putative_linker_R_withbc


echo "Restoring non-linker_R reads to pre-trimmed state..."
perl ../perl/names_to_fastq.pl \
  ${fastq_file}_LTR_R_p_trim_no_linker_R_barcode \
  ${fastq_file}_LTR_R_pairs-with-p_noLTR-R_nolinkF_trimLTRF.fastq

test_file ${fastq_file}_LTR_R_p_trim_no_linker_R_barcode.fastq

rm ${fastq_file}_LTR_R_p_trim_no_linker_R_barcode
rm ${fastq_file}_LTR_R_pairs-with-p_noLTR-R_nolinkF_trimLTRF.fastq


# Concatenate the output files
cat ${fastq_file}_LTR_R_p_trim_no_linker_R_barcode.fastq \
 ${fastq_file}_LTR_R_p_trim_uo.fastq \
 ${fastq_file}_LTR_R_p_trim_withbc.fastq_bctrimmed.fastq \
 > ${fastq_file}_LTR_R_p_trim.fastq

test_file ${fastq_file}_LTR_R_p_trim.fastq

rm ${fastq_file}_LTR_R_p_trim_no_linker_R_barcode.fastq
rm ${fastq_file}_LTR_R_p_trim_uo.fastq
rm ${fastq_file}_LTR_R_p_trim_withbc.fastq_bctrimmed.fastq

################################################################################
# Step 11.
# Combine the two above output files with the three 'good' input files.

echo ""
echo "***Step 11. Combine -q and non-q pairs into their respective files.***"

echo ""
echo "Combining the good reads..."
cat cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink_rm-R-link.fastq \
 cut_LTR_F_${fastq_file}_short.fastq_trim-R-link.fastq \
 cut_LTR_R_${fastq_file}_overlap_trim-F-link.fastq > ${fastq_file}_non-p

test_file ${fastq_file}_non-p

echo ""
echo "Combining the p-pairs..."
cat ${fastq_file}_LTR_F_combo_p_trim.fastq \
 ${fastq_file}_LTR_R_p_trim.fastq > ${fastq_file}_p

test_file ${fastq_file}_p

rm ${fastq_file}_LTR_F_combo_p_trim.fastq
rm ${fastq_file}_LTR_R_p_trim.fastq

###############################################################################
# Step 12.
# Align remaining reads with bowtie, and sort.
# The sort is probably not strictly necessary, but it should be fast enough that
# it won't hurt.

echo ""
echo "***Step 12. Align -q and non-q reads with bowtie.***"

# A note on the options used:  -a indicates that all alignments will be 
# considered; --best means that the alignments will be ranked, giving preference
# to alignments with 1) fewer mismatches in the seed region, and 2) alignments 
# with a lower sum of quality scores in the mismatches they DO have; --strata 
# means that the alignments will be stratified by the number of reads they have
# in their seed region; -m 1 means that the only reads to be reported will be 
# those with only a single alignment in the best stratum.  Such alignments 
# aren't necessarily unique, but they ARE unambiguously the best alignment based
# on mismatches in the seed region.  --chunkmbs 256 requisitions additional
# memory for Bowtie, and -t reports timing statistics.

echo ""
echo "Aligning non-p reads to ${bowtie_index}..."
echo "Bowtie says:"
bowtie -t -a -m 1 --best --strata --chunkmbs 256 --sam ${bowtie_index} \
 ${fastq_file}_non-p ${fastq_file}_non-p_bowtie_strata.sam

echo ""
test_file ${fastq_file}_non-p_bowtie_strata.sam

if [ "$keep" = "off" ]; then rm ${fastq_file}_non-p; fi

# Remove the header; keep it for later.
cat ${fastq_file}_non-p_bowtie_strata.sam | awk '$1 ~ /^@/' > \
    ${fastq_file}_non-p_bowtie_strata.sam_header

# cut -f1-3 ${fastq_file}_non-p_bowtie_strata.sam_header > \
#     ${fastq_file}_non-p_bowtie_strata.sam_header_clean
# 
# test_file ${fastq_file}_non-p_bowtie_strata.sam_header_clean
# rm ${fastq_file}_non-p_bowtie_strata.sam_header

cat ${fastq_file}_non-p_bowtie_strata.sam | awk '$1 !~ /^@/' > \
    ${fastq_file}_non-p_bowtie_strata_aln.sam

test_file ${fastq_file}_non-p_bowtie_strata_aln.sam
if [ "$keep" = "off" ]; then rm ${fastq_file}_non-p_bowtie_strata.sam; fi

# Repeat for the p-pairs file.
echo ""
echo "Aligning p-reads to ${bowtie_index}..."
echo "Bowtie says:"
bowtie -t -a -m 1 --best --strata --chunkmbs 256 --sam ${bowtie_index} \
 ${fastq_file}_p ${fastq_file}_p_bowtie_strata.sam

echo ""
test_file ${fastq_file}_p_bowtie_strata.sam

if [ "$keep" = "off" ]; then rm ${fastq_file}_p; fi

cat ${fastq_file}_p_bowtie_strata.sam | awk '$1 !~ /^@/' > \
    ${fastq_file}_p_bowtie_strata_aln.sam

test_file ${fastq_file}_p_bowtie_strata_aln.sam
if [ "$keep" = "off" ]; then rm ${fastq_file}_p_bowtie_strata.sam; fi

# Add "|p" to the end of the names of all p-pair reads.
echo ""
echo "Adding |p to the end of p-pair alignment names..."
awk -v OFS="\t" '{ if(NF == 14){print $1"|p",$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}else{print $1"|p",$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12} }' \
    ${fastq_file}_p_bowtie_strata_aln.sam > ${fastq_file}_p_bowtie_strata_aln_rename.sam

test_file ${fastq_file}_p_bowtie_strata_aln_rename.sam
rm ${fastq_file}_p_bowtie_strata_aln.sam

# Concatenate the two output files together.  Sorting just makes it easier to 
# read.
echo ""
echo "Concatenating the non-p and p alignments together, "
echo " removing unmapped reads, then sorting by chromosome and alignment..."
cat ${fastq_file}_non-p_bowtie_strata_aln.sam \
    ${fastq_file}_p_bowtie_strata_aln_rename.sam | awk '$3 != "*"' | \
    sort -k3,3 -k4,4n > ${fastq_file}_p-pairs_bowtie_strata_sort.sam

test_file ${fastq_file}_p-pairs_bowtie_strata_sort.sam

rm ${fastq_file}_non-p_bowtie_strata_aln.sam
rm ${fastq_file}_p_bowtie_strata_aln_rename.sam

################################################################################
# Step 13.
# Use the 'examine' script to remove reads with suspicious pairs.
#
# This version has two 'good' output files.  The one named _nofarpair has only
# non-p reads in it; this one is used in future steps.  The other is
# _nofarpair_with-p, and is just there to provide some more information on 
# the main output file, if needed.

echo ""
echo "***Step 13. Remove reads with suspicious pairs, and only keep -q reads.***"


echo ""
echo "Removing reads that did not align near their |p pair, or aligned in the"
echo " wrong orientation.  bowtie_examine_p-pairs_v1.4.pl says:"
perl ../perl/bowtie_examine_p-pairs_sam.pl \
 ${fastq_file}_p-pairs_bowtie_strata_sort.sam

mv ${fastq_file}_p-pairs_bowtie_strata_sort.sam_nofarpair \
    ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair.sam

echo ""
test_file ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair.sam

if [ "$keep" = "off" ]; then rm ${fastq_file}_p-pairs_bowtie_strata_sort.sam; fi
if [ "$keep" = "off" ]; then rm ${fastq_file}_p-pairs_bowtie_strata_sort.sam_nofarpair_with-p; fi
if [ "$keep" = "off" ]; then rm ${fastq_file}_p-pairs_bowtie_strata_sort.sam_badorientation; fi
if [ "$keep" = "off" ]; then rm ${fastq_file}_p-pairs_bowtie_strata_sort.sam_toofar; fi
if [ "$keep" = "off" ]; then rm ${fastq_file}_p-pairs_bowtie_strata_sort.sam_difchrom; fi


cat ${fastq_file}_LTR_R_p_trim_withbc.fastq_barcodes.fastq \
 revcom_${fastq_file}_LTR_F_combo_p_trim_withbc.fastq_barcodes.fastq > \
 ${fastq_file}_p-pair_barcodes.fastq

echo ""
test_file ${fastq_file}_p-pair_barcodes.fastq

if [ "$keep" = "off" ]; then rm revcom_${fastq_file}_LTR_F_combo_p_trim_withbc.fastq_barcodes.fastq; fi
if [ "$keep" = "off" ]; then rm ${fastq_file}_LTR_R_p_trim_withbc.fastq_barcodes.fastq; fi

echo ""
echo "Adding |p to the p-pair barcode names..."
cat ${fastq_file}_p-pair_barcodes.fastq | ../perl/add_p_v1.0.pl > \
 ${fastq_file}_p-pair_barcodes.fastq_withp

test_file ${fastq_file}_p-pair_barcodes.fastq_withp

rm ${fastq_file}_p-pair_barcodes.fastq


echo ""
echo "Concatenating barcode reads.  The file is not to be used directly, as"
echo "duplicate reads have not been removed..."
cat cut_linkerrevcom_${fastq_file}_barcodes.fastq \
 cut_linker_R_${fastq_file}_barcodes.fastq > ${fastq_file}_barcodes_cat.fastq

echo ""
test_file ${fastq_file}_barcodes_cat.fastq

if [ "$keep" = "off" ]; then rm cut_linkerrevcom_${fastq_file}_barcodes.fastq; fi
if [ "$keep" = "off" ]; then rm cut_linker_R_${fastq_file}_barcodes.fastq; fi

# Remove reads that have a different barcode than their paired read
echo ""
echo "Removing reads with inconsistent barcodes..."
echo "barcode_compare.pl says:"
perl ../perl/barcode_compare_v1.1.pl ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair.sam \
 ${fastq_file}_barcodes_cat.fastq ${fastq_file}_p-pair_barcodes.fastq_withp 

mv ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair.sam_consistentbc \
    ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair_consistentbc.sam

echo ""
test_file ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair_consistentbc.sam

if [ "$keep" = "off" ]; then rm ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair.sam; fi
if [ "$keep" = "off" ]; then rm ${fastq_file}_p-pair_barcodes.fastq_withp; fi
if [ "$keep" = "off" ]; then rm ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair.sam_inconsistentbc; fi


################################################################################
# Step 14.
# Modify the bowtie output so the alignments reflect the bases directly adjacent
# to the LTR.

echo ""
echo "***Step 14. Identify bases adjacent to the LTR.***"

# Need to make an LTR_F file.  It's best to do this post-trim, so there won't be
# any overlap with LTR_R reads.

echo ""
echo "Concatenating LTR_F_long and _short..."
cat cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink_rm-R-link.fastq \
 cut_LTR_F_${fastq_file}_short.fastq_trim-R-link.fastq > LTR_F_temp

echo ""
test_file LTR_F_temp

if [ "$keep" = "off" ]; then rm cut_LTR_F_${fastq_file}_long.fastq_no-bad-F-link_pairhaslink_rm-R-link.fastq; fi
if [ "$keep" = "off" ]; then rm cut_LTR_F_${fastq_file}_short.fastq_trim-R-link.fastq; fi

echo ""
echo "Modifying the bowtie output to indicate the bases adjacent to the LTR."
echo "Note that the base reported will be the leftmost base of the N bases next"
echo "to the LTR (defined by the insert selected at runtime)."
echo "bowtie_LTR_adjacent.pl says:"
perl ../perl/bowtie_LTR_adjacent_sam.pl \
 LTR_F_temp cut_LTR_R_${fastq_file}_overlap_trim-F-link.fastq \
 ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair_consistentbc.sam ${footprint}

mv ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair_consistentbc.sam_LTRadjacent ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair_consistentbc_LTRadjacent

echo ""
test_file ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair_consistentbc_LTRadjacent
# It's not REALLY a SAM file right now, since the 15th column isn't formatted
# correctly - so I won't append ".sam" onto it.  Note it is STILL 1-based.


rm LTR_F_temp
if [ "$keep" = "off" ]; then rm cut_LTR_R_${fastq_file}_overlap_trim-F-link.fastq; fi
if [ "$keep" = "off" ]; then rm ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair_consistentbc.sam; fi

# Sort - this time, it IS necessary for the next step.

echo ""
echo "Sorting Bowtie alignment by name, then by chromosome, then by alignment..."
sort -k1,1 -k3,3 -k4,4n ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair_consistentbc_LTRadjacent \
 > ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair_consistentbc_LTRadjacent_sort

test_file ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair_consistentbc_LTRadjacent_sort

rm ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair_consistentbc_LTRadjacent

################################################################################
# Step 15.
# Pare down the bowtie reads so each fragment is represented by only one line.
# Note that the file needs to be sorted for this to work correctly.

echo ""
echo "***Step 15. Reduce bowtie output to one entry per fragment.***"

awk -v OFS="\t" '{base=substr($1, 1, length($1)-2); if(base != oldbase){print $1,$15,$3,$4}; oldbase=base}' \
   ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair_consistentbc_LTRadjacent_sort \
   > ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair_consistentbc_LTRadjacent_sort_singlefrag

echo ""
test_file ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair_consistentbc_LTRadjacent_sort_singlefrag


################################################################################
# Step 16.
# Split the fragments into groups, according their barcode. The file
# ${barcode_ref}_seq is used to accomplish this.

# The default files split the reads into four groups, corresponding to the four
# T75 flasks in which the cells were initially transfected with MLV.

echo ""
echo "***Step 16. Split fragments into barcode groups.***"


# use whittle_and_count on _singlefrag to add barcode sequence information

echo ""
echo "Adding barcode information."
echo "barcode_whittle_and_count_sam.pl says:"
perl ../perl/barcode_whittle_and_count_sam.pl \
 ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair_consistentbc_LTRadjacent_sort_singlefrag \
 ${fastq_file}_barcodes_cat.fastq ${barcode_ref}_seq

echo ""
test_file ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair_consistentbc_LTRadjacent_sort_singlefrag_barcodes

rm ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair_consistentbc_LTRadjacent_sort_singlefrag_barcodesum
rm ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair_consistentbc_LTRadjacent_sort_singlefrag

# Throw in dummy columns (just to make formatting easier down the road)

echo ""
echo "Reformatting the _singlefrag_barcodes file..."
cat ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair_consistentbc_LTRadjacent_sort_singlefrag_barcodes \
 | awk '{print $1"\t"$2"\t"$3"\t"$4"\t(-)\t(-)\t(-)\t(-)\t(-)\t(-)\t(-)\t(-)\t(-)\t(-)\t(-)\t"$5"\t"$6}' > \
 ${fastq_file}_combo_split

echo ""
test_file ${fastq_file}_combo_split

if [ "$keep" = "off" ]; then rm ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair_consistentbc_LTRadjacent_sort_singlefrag_barcodes; fi


# Run barcode-group

echo ""
echo "Adding barcode group information."
echo "bowtie_barcode-group_sam.pl says:"
perl ../perl/bowtie_barcode-group_sam.pl ../${barcode_ref} \
 ${fastq_file}_combo_split

echo ""
test_file ${fastq_file}_combo_split_grouped


rm ${fastq_file}_combo_split


# use awk to split into ${group_number} parts

echo ""
echo "Splitting into ${group_number} groups, and sorting..."
i=1
while
  [ $i -le ${group_number} ]
do
  cat ${fastq_file}_combo_split_grouped | \
  awk -v num="${i}" '$18 == num && $2 == "+" {print $3"\t"$4}' | \
  sort -k1,1 -k2,2n | uniq -c > ${fastq_file}_combo_split_grouped_plus_uniq_${i}
  
  cat ${fastq_file}_combo_split_grouped | \
  awk -v num="${i}" '$18 == num && $2 == "-" {print $3"\t"$4}' | \
  sort -k1,1 -k2,2n | uniq -c > ${fastq_file}_combo_split_grouped_minus_uniq_${i}
  
  test_file ${fastq_file}_combo_split_grouped_plus_uniq_${i}
  test_file ${fastq_file}_combo_split_grouped_minus_uniq_${i}
  
  i=$[ $i + 1 ]
done

if [ "$keep" = "off" ]; then rm ${fastq_file}_combo_split_grouped; fi


################################################################################
# Step 17.
# Pare down reads in which the leftmost base of the integration sites are 
# within 5 bp of each other.  Note this is done separately for each group; 
# if two fragments are from different flasks, they can't represent the same 
# event.

echo ""
echo "***Step 17. Pare down reads within 5 bp of each other.***"

echo ""
echo "Removing nearby reads for each group of fragments..."
i=1
while
  [ $i -le ${group_number} ]
do
  perl ../perl/remove_nearby_blacklist.pl \
    ${fastq_file}_combo_split_grouped_plus_uniq_${i} 5
  perl ../perl/remove_nearby_blacklist.pl \
    ${fastq_file}_combo_split_grouped_minus_uniq_${i} 5
  
  test_file ${fastq_file}_combo_split_grouped_plus_uniq_${i}.not_in_5
  test_file ${fastq_file}_combo_split_grouped_minus_uniq_${i}.not_in_5
  
  if [ "$keep" = "off" ]; then rm ${fastq_file}_combo_split_grouped_plus_uniq_${i}; fi
  if [ "$keep" = "off" ]; then rm ${fastq_file}_combo_split_grouped_minus_uniq_${i}; fi

  i=$[ $i + 1 ]
done



################################################################################
# Step 18.
# If a fragment cutoff was defined when the script was submitted to run, remove
# integrations with fewer fragments than the cutoff.

echo ""
echo "***Step 18. Remove integrations below fragment cutoff.***"

echo ""
echo "Removing integrations with fewer than ${cutoff} fragments..."
i=1
while
  [ $i -le ${group_number} ]
do
  awk -v cut="$cutoff" '$1 >= cut' ${fastq_file}_combo_split_grouped_plus_uniq_${i}.not_in_5 > \
    ${fastq_file}_combo_split_grouped_plus_uniq_${i}.not_in_5_cutoff${cutoff}
  awk -v cut="$cutoff" '$1 >= cut' ${fastq_file}_combo_split_grouped_minus_uniq_${i}.not_in_5 > \
    ${fastq_file}_combo_split_grouped_minus_uniq_${i}.not_in_5_cutoff${cutoff}
  
  test_file ${fastq_file}_combo_split_grouped_plus_uniq_${i}.not_in_5_cutoff${cutoff}
  test_file ${fastq_file}_combo_split_grouped_minus_uniq_${i}.not_in_5_cutoff${cutoff}

  if [ "$keep" = "off" ]; then rm ${fastq_file}_combo_split_grouped_plus_uniq_${i}.not_in_5; fi
  if [ "$keep" = "off" ]; then rm ${fastq_file}_combo_split_grouped_minus_uniq_${i}.not_in_5; fi
  
  i=$[ $i + 1 ]
done

################################################################################
# Step 19.
# Make a modified Bowtie-style output file for the integrations that remain. I
# refer to this as an "island" file, mostly because it's after integration sites
# have been made somewhat more discrete - nearby reads  with fever reads had 
# just been removed in the previous steps.

echo ""
echo "***Step 19. Make a SAM file of the remaining reads.***"


echo ""
echo "Combining the ${group_number} files into one, sorting..."
cat ${fastq_file}_combo_split_grouped_plus_uniq_*.not_in_5_cutoff${cutoff} \
   ${fastq_file}_combo_split_grouped_minus_uniq_*.not_in_5_cutoff${cutoff} | \
   sort -k2,2 -k3,3n > \
   ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}

test_file ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}


# Make island file.  Note that entries from both strands will be printed for
# a given site; this can be corrected a few lines later.

echo ""
echo "Making the 'island' file (modified Bowtie format)."
echo "notinx_to_bowtie_v1.0.pl says:"
perl ../perl/notinx_to_bowtie_v1.0.pl \
 ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff} \
 ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair_consistentbc_LTRadjacent_sort

echo ""
test_file ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair_consistentbc_LTRadjacent_sort_island

if [ "$keep" = "off" ]; then rm ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair_consistentbc_LTRadjacent_sort; fi


# Rename file

echo ""
echo "Renaming file..."
mv ${fastq_file}_p-pairs_bowtie_strata_sort_nofarpair_consistentbc_LTRadjacent_sort_island \
 ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island

test_file ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island


# Add barcodes

echo ""
echo "Adding barcode information."
echo "barcode_whittle_and_count_sam.pl says:"
perl ../perl/barcode_whittle_and_count_sam.pl \
 ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island \
 ${fastq_file}_barcodes_cat.fastq ${barcode_ref}_seq
 
echo ""
test_file ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes

if [ "$keep" = "off" ]; then rm ${fastq_file}_barcodes_cat.fastq; fi
rm ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island


# Sort, add barcode groups

echo ""
echo "Sorting, adding barcode group information..."
sort -k3,3 -k4,4n \
 ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes \
 > ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort

test_file ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort
rm ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes



echo "bowtie_barcode-group_sam.pl says:"
perl ../perl/bowtie_barcode-group_sam.pl \
 ../${barcode_ref} \
 ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort 

echo ""
test_file ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped

rm ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort



# Use of notinx_to_bowtie, above, guarantees that the sites selected pass the 
# cutoff in at least one barcode group.  Now, further reduce it to be a list 
# of ONLY the barcode groups (and strands) in which the site passes the cutoff.

# Split the island file by strand
echo ""
echo "Splitting the island file by strand..."
awk '$15 == "+"' \
   ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped \
   > ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_plus

awk '$15 == "-"' \
   ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped \
   > ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_minus

test_file ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_plus
test_file ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_minus


# Print the sites that passed the cutoff AND the group to which they belong.
# First, combine all cutoff files into two "stranded" files.
echo ""
echo "Combining cutoff files by appropriate strand..."
i=1
while 
  [ $i -le ${group_number} ]
do 
  awk -v i="$i" '{print $0"\t"i}' ${fastq_file}_combo_split_grouped_plus_uniq_${i}.not_in_5_cutoff${cutoff} \
   >> ${fastq_file}_combo_split_grouped_plus_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_group-info
  
  awk -v i="$i" '{print $0"\t"i}' ${fastq_file}_combo_split_grouped_minus_uniq_${i}.not_in_5_cutoff${cutoff} \
   >> ${fastq_file}_combo_split_grouped_minus_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_group-info
  
  i=$[ $i + 1 ] 
done

echo ""
test_file ${fastq_file}_combo_split_grouped_plus_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_group-info
test_file ${fastq_file}_combo_split_grouped_minus_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_group-info

if [ "$keep" = "off" ]; then rm ${fastq_file}_combo_split_grouped_plus_uniq_*.not_in_5_cutoff${cutoff}; fi
if [ "$keep" = "off" ]; then rm ${fastq_file}_combo_split_grouped_minus_uniq_*.not_in_5_cutoff${cutoff}; fi


# Only keep the sam file entries that are the correct position for that group
echo ""
echo "Removing entries that failed to pass the cutoff for the group/strand..."
perl ../perl/sam_good-groups.pl \
 ${fastq_file}_combo_split_grouped_plus_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_group-info \
 ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_plus \
 > ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_plus_good

perl ../perl/sam_good-groups.pl \
 ${fastq_file}_combo_split_grouped_minus_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_group-info \
 ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_minus \
 > ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_minus_good

echo ""
test_file ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_plus_good
test_file ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_minus_good


if [ "$keep" = "off" ]; then rm ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_group-info; fi
if [ "$keep" = "off" ]; then rm ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped; fi

echo ""
echo "Combining the two strand files together..."
cat ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_plus_good \
  ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_minus_good \
  | sort -k3,3 -k4,4n > ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_good

test_file ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_good

if [ "$keep" = "off" ]; then rm ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_plus_good; fi
if [ "$keep" = "off" ]; then rm ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_minus_good; fi


################################################################################
# Step 20.
# Convert to .BED format.

echo ""
echo "***Step 20. Convert to .BED format.***"

# Be careful, here - though I add a BED header appropriate for the UCSC browser,
# I plan on using the Ensembl index to do the actual mapping.

echo ""
echo "Converting output to .BED format..."

cat ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_good \
  | awk '{print $3"\t"$4"\t"$15"\t"$16 }' | sort | uniq | sort -k4,4 > \
  ${fastq_file}_1-${group_number}-grouped_not-in-5_cutoff${cutoff}_job${JOB_ID}_demibed

test_file ${fastq_file}_1-${group_number}-grouped_not-in-5_cutoff${cutoff}_job${JOB_ID}_demibed

sort -k1,1 ../${barcode_ref} > ${barcode_ref}_barsort

test_file ${barcode_ref}_barsort


join -1 4 -2 1 ${fastq_file}_1-${group_number}-grouped_not-in-5_cutoff${cutoff}_job${JOB_ID}_demibed \
  ${barcode_ref}_barsort \
  | awk -v OFS="\t" '{print $2,$3-1,$4,$6"_"$7"_"$8"_"$9}' | sort | uniq \
  | awk -v footprint="$footprint" '{print $1"\t"$2"\t"$2+footprint"\t"$4"\t0\t"$3"\t"$2"\t"$2+footprint}' \
  | sort -k1,1 -k2,2n \
  | awk -v description="${name}" -v insjob="${insert}${JOB_ID}" 'BEGIN{print "track name="insjob" type=bed description=\""description"\""}{print}' \
  > ${fastq_file}_1-${group_number}-grouped_not-in-5_cutoff${cutoff}_${insert}_job${JOB_ID}.bed


echo ""
test_file ${fastq_file}_1-${group_number}-grouped_not-in-5_cutoff${cutoff}_${insert}_job${JOB_ID}.bed

if [ "$keep" = "off" ]; then rm ${fastq_file}_1-${group_number}-grouped_not-in-5_cutoff${cutoff}_job${JOB_ID}_demibed; fi
if [ "$keep" = "off" ]; then rm ${barcode_ref}_barsort; fi


################################################################################

# Step 21.
# Make the output into a real SAM file, then to BAM.

echo ""
echo "***Step 21. Convert file to SAM/BAM format***"

echo ""
echo "Adding SAM formatting to extra columns, adding job number to name..."
cat ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_good \
    | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,"XO:A:"$15,"BC:Z:"$16,"XP:Z:"$17,"XG:i:"$18}' \
    > ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_job${JOB_ID}.sam

test_file ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_job${JOB_ID}.sam

echo ""
echo "Putting header on the SAM file..."
cat ${fastq_file}_non-p_bowtie_strata.sam_header \
    ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_job${JOB_ID}.sam \
    > ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_job${JOB_ID}_header.sam

test_file ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_job${JOB_ID}_header.sam
if [ "$keep" = "off" ]; then rm ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_job${JOB_ID}.sam; fi

echo ""
echo "Converting to BAM..."
samtools view -Sb \
    ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_job${JOB_ID}_header.sam \
    > ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_job${JOB_ID}_header.bam

test_file ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_job${JOB_ID}_header.bam
rm ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_job${JOB_ID}_header.sam


echo ""
echo "Adding insert and job number to the remaining 'island' file..."
mv ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodesum \
 ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodesum_${insert}_job${JOB_ID}


################################################################################

# Step 22.
# Clean up files.

# First, we move the input file and the final output file, so they won't be 
# included in the tarball.

echo ""
echo "***Step 22. Clean up files.***"

echo ""
echo "Making output directory..."
mkdir $PWD/Output_${name}_${insert}_$JOB_ID


function move_file
{
 if 
   [ -e Output_${name}_${insert}_$JOB_ID/$1 ]
 then 
   throw_error "Can't move ${1}; a file with that name already exists in parent directory!"
 else  
   mv $1 Output_${name}_${insert}_${JOB_ID}/ || throw_error "Didn't move ${1}!"
 fi
}

move_file ${fastq_file}_1-${group_number}-grouped_not-in-5_cutoff${cutoff}_${insert}_job${JOB_ID}.bed
move_file ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_job${JOB_ID}_header.bam
move_file ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodesum_${insert}_job${JOB_ID}

# rm ${fastq_file}_1-${group_number}-grouped_not-in-5_cutoff${cutoff}_${insert}_job${JOB_ID}.bed
# rm ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_job${JOB_ID}_header.bam
# rm ${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodesum_${insert}_job${JOB_ID}


echo ""
echo "Storing any remaining files in a tarball, moving it out of the working directory..."
if [ "$keep" = "on" ] 
then 
   tar czf ${fastq_file}_${name}_intermediatefiles_${insert}_$JOB_ID.tgz *${fastq_file}* 
   test_file ${fastq_file}_${name}_intermediatefiles_${insert}_$JOB_ID.tgz
   
   move_file ${fastq_file}_${name}_intermediatefiles_${insert}_$JOB_ID.tgz
   
fi

# Move the output directory out of the working directory
if 
   [ -e ../Output_${name}_${insert}_$JOB_ID ]
 then 
   throw_error "Can't move Output_${name}_${insert}_$JOB_ID; a file with that name already exists in parent directory!"
 else  
   mv Output_${name}_${insert}_${JOB_ID}/ .. || throw_error "Didn't move Output_${name}_${insert}_$JOB_ID!"
 fi 

# Remove the working directory, and the files remaining therein.

echo ""
echo "Removing working directory and files..."
cd ..
rm -r $workdir


echo ""
echo "$0 finished!"
exit
