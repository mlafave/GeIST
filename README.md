GeIST
=====

GeIST (Genomic Integration Site Tracker) provides the scripts used to identify
integration sites from the assay described in Varshney et al. 2013
(http://dx.doi.org/10.1101/gr.151464.112) and LaFave et al. 2014
(http://dx.doi.org/10.1093/nar/gkt1399). GeIST is currently capable of
identifying integrations from MLV, AAV, Tol2, or Ds-based vectors.

Note that the version of the assay used in LaFave et al. 2014 involved mixing
the barcodes in such a way as to be able to identify multiple integrations in
the same integration site.  The current version of GeIST does not automatically
assume that barcodes have been mixed in this way.  If you wish to apply this
script to an experiment in which the barcodes have been mixed, please refer to
the calculations used to print the BED file in the previous version of GeIST,
hosted at http://research.nhgri.nih.gov/software/GeIST/, and make the
appropriate modifications to your barcode file (see below).


Requirements
------------

* Perl (all scripts use version 5.8.8)
* cutadapt version 1.3
* bowtie version 1.0.0, along with the relevant index
* samtools version 0.1.19
* bamtools version 1.0.2, if starting with a BAM instead of a FASTQ file
* Sufficient RAM to load indexes/hold intermediate values in memory (many of
	the perl scripts store the contents of large files in hashes).  32 GB of
	RAM should be sufficient when analyzing a 50 GB FASTQ input file.
* Sufficient disk space to hold the output files.  The amount needed depends on
    the input; if the input is 50 GB, at least 150 GB free space will be needed.  
    If intermediate files are not kept (the default - see below), less space is 
    needed, but it will likely still be over 100 GB.  


Running the script
------------------

There is a single main shell script, `geist.sh`, that executes the workflow.  It
starts with either a FASTQ or BAM file of sequencing reads, and ends with a BAM
of aligned, trimmed reads and a BED file of integration locations.  

The script was developed to be run on an SGE cluster; as such, it can be 
submitted using 'qsub', but doing so it not required.  Output files will have 
the job number assigned by the queue included in their name; if qsub is not 
used, that part of the name will remain blank.

Run the program with the following command (should be run from within the 
GeIST directory):

qsub -l mem_free=32G geist.sh [ -v ] -r barcodes -b index -i insert [-c cutoff] [-n name] [ -k ] reads.[fastq | bam]


More information on the elements of the command:

* If 'qsub' is used, -l can be specified to request memory resources.  This
isn't required, but it may be helpful to request a large amount of memory when
working with a large FASTQ file.  For a 50 GB FASTQ file, 32 GB of memory is
suggested (if available).

* Any files specified in the command - such as the bowtie index, the barcode
file, or the input sequences - can be specified as either relative or absolute
paths.

* After specifying the name of the shell script, the order in which the
remaining files are entered is not important, with the exception that the input
FASTQ or BAM file must be listed last.

* If '-v' is used, GeIST will print the current version number to standard
output and exit.

* The file specified with the '-r' option is the file in which the barcodes are
listed.  

The file must be of the following format:

	GAAAAAA	1	plate1	A01	LINE1	A
	GAAAATT	1	plate1	A01	LINE1	A
	GAAAAGG	1	plate1	A02	LINE1	B
	GAAAACC	1	plate1	A02	LINE1	B
	GAAATAT	1	plate1	A03	LINE1	C
	GAAATTA	1	plate1	A03	LINE1	C
	GAAATGC	2	plate1	A04	LINE9	A
	GAAATCG	2	plate1	A04	LINE9	A

	* The first column contains the barcode sequence, starting with the 5' base
	as viewed on a read in which the linker is on the rightmost side of the
	read.
	* The second column is an integer that indicates the "group" to which the
	barcode belongs.  For example, the barcodes can be used to multiplex
	different samples, and each sample should have an integer associated with
	it so the program knows how to separate them back out.  If everything is
	from the same sample, just set this column to 1 on ever line.  The integers 
	in this column need to be sequential, and the first integer in the sequence 
	needs to be 1.
	* The third column is a string indicating which plate the sample came from.
	If they're all the same, just set them all to "plate1".
	* The fourth column is a string showing the well associated with the
	barcode.
	* The fifth column is a string that shows the name of the "group". This can
	be used to indicate flask number, organism/cell line, etc.  The values in
	this column correspond to those in column 2.
	* The fifth column is a string indicating the relevant subset of the group
	from the fourth column.  If using model organisms, this is a useful place
	to put a character indicating which individual this barcode was associated
	with.

* Use the -b option to enter the path to the Bowtie index to be used, ending
with the basename of the index.  The basename is the portion of the name common
to all six index files (the part that comes before the first period). 

* Use the -i option enter which type of integration you want to detect.  The
four options are mlv, aav, tol2, and ds.

* Use -c to enter a "cutoff" value.  This specifies the minimum number of
sequenced fragments that constitute an integration. The default is 0, which
allows all integration sites regardless of fragment count.  

* Use -n to enter a name without spaces.  This will be used in the name of the
working directory, the output directory, and the tarball of intermediate files.
It's a good way to tell different runs apart, especially if the job number won't
be printed (such as when the script is run without qsub).  The default is
"GeIST". 

* The -k entry is optional.  By default, the script will remove unnecessary
temporary files as it goes along.  If you would like to examine these temporary
files (for example, to see how the final output was derived), use -k.

* The last file should be either a BAM file or a FASTQ file of LM-PCR
paired-end sequencing data.  If a BAM file is submitted, it is automatically
converted to a FASTQ file by bamtools.  If starting from a FASTQ file, ensure
that the names of paired reads end in /1 and /2. Bamtools does such processing
automatically. 
GeIST accepts only a single file for sequence input, so if the user needs to
analyze multiple files (for example, if the paired-end sequence data was
returned in two separate files), those files should be concatenated prior to
running GeIST.


The script will create a working directory into which it will create the
intermediate files.  The fastq (or BAM) file and the barcode file are not
copied into the working directory. As such, they should remain in the same
directory for the duration of the run of the program.  When the script is
finished, relevant files are moved to an Output directory, and the working
directory is deleted.


Testing the script
------------------

	qsub -l mem_free=2.5G geist.sh -r BC1-9_barcodes+G_4_groups_example -b ~/indexes/hg19 -i mlv -n example -k script_test.fastq

`script_test.bam` can be used in place of `script_test.fastq`, and should produce
identical results.

The above assumes that you're using a human genome index with basename "hg19",
such as the one available from the bowtie website at
ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/hg19.ebwt.zip. The "indexes"
directory in this example is located in the `$HOME` directory, but that doesn't
have to be the case.

These settings will identify MLV integrations in the human genome, based on a
few reads in `script_test.fastq` or `script_test.bam`.  It uses a fragment cutoff of 0.  The expected
output can be found in the `Output_expected/` directory.


Understanding the output
------------------------

The main output file will have the form 
`${fastq_file}_1-${group_number}-grouped_not-in-5_cutoff${cutoff}_${insert}_job${JOB_ID}.bed`.  
The name indicates the file submitted, the number of barcode groups used, that
integrations within 5 bp of another integration with more reads were removed,
the cutoff used, the type of integration detected, and the job number (if qsub
was used to submit the script). The format is that of a typical BED file.  The
first line is a track header, useful if you want to examine the output on the
UCSC genome browser.  The columns are as follows:

	Column 1: Chromosome name
	Column 2: Start position (0-based)
	Column 3: End position (1-based)
	Column 4: The plate, well, line, and individual, according to barcode detected.  
	Column 5: "Score".  This field doesn't mean anything in this context, so I set everything to 0.
	Column 6: Orientation of the integration
	Column 7: thickStart.  Only used for visualization purposes.
	Column 8: thickEnd.  Only used for visualization purposes.


Information on the individual reads can be found in 
`${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodes_sort_grouped_job${JOB_ID}_header.bam`.  
Note that this file only contains the reads in which part of the integration
itself (the "LTR") was detected.  While not all of the elements detected by
GeIST have long terminal repeats, "LTR" is used throughout the code as
shorthand for the amplified portion of the end of the integrated element. 

This BAM file can be converted to SAM format by using `samtools view`. The
first fourteen columns are standard SAM entries. The final four columns contain
information specific to the assay:
	
	* XO: Indicates the orientation of the integration, relative to the end
	used for sequencing.  For a given fragment, both read /1 and read /2 will
	have the same entry in this field.  Note that this orientation is not
	necessarily the same as the orientation in column 2.  Column 2 is the
	orientation of the read, and column 15 is the orientation of the entire
	fragment.  Note that this is only informative if the primer used for insert
	amplification could only bind to one end of the insert.  For example, AAV
	has identical termini on both ends, in which case this column is not
	informative.
	* BC: The sequence of the barcode detected on this fragment.
	* XP: Indicates if the barcode was detected on the read itself
	(bc_from_read) or from the paired read (bc_from_pair).  Typically, short
	fragments will be bc_from_read and long fragments will be bc_from_pair, but
	sometimes quality scores make it such that this is not the case.
	* XG: The group number.


`${fastq_file}_combo_split_grouped_uniq_1-${group_number}.not_in_5_cutoff${cutoff}_island_barcodesum_${insert}_job${JOB_ID}`
simply indicates the number of fragments on which each barcode occurred.  This
file can be useful in determining which barcodes are over- or under-represented.


Adding additional integrated elements
-------------------------------------

If you would like to use GeIST to recover an element not among the four
currently supported, you can modify geist.sh to do so. In the section labeled
"Set parameters for the indicated insert", add the following line:

	echo $insert | grep -i -q -w '(some_name)' && insert=(some_name)

Add the following code between the `elif` statements near the top of the script
(around line 215):

	elif [ "$insert" = "(some_name)" ]
	then
	   echo "(some_optional_statement)"
	   LTR=(3'_end_of_element)
	   footprint=(integration_footprint_size)

The details for the placeholders are as follows:

* `(some_name)`: A string without whitespace, used to indicate your element.
This is the part that you would specify with the `-i` flag when running GeIST.
For example, if you were trying to recover HIV, this field could be `hiv`.
* `(some_optional_statement)`: A message printed as the script runs, reporting
which element was used.
* `(3'_end_of_element)`: The reverse complement of the sequence of the 3' end
of the integrated element. Typically, this consists of the sequence of the
primer used in the final round of LM-PCR, followed by any bases appearing
between the end of the primer and the end of the element. For example, the
primer used in our Ds experiment was `TTTACCGACCGTTACCGACCGTTTTCATC`, and the
sequence of the 3' end of the provirus was
`GAGGTATTTTACCGACCGTTACCGACCGTTTTCATCCCTA`. The combination of the primer and
the four additional bases between the end of the primer and the end of the
element is `TTTACCGACCGTTACCGACCGTTTTCATCCCTA`, the reverse complement of which
is `TAGGGATGAAAACGGTCGGTAACGGTCGGTAAA`, which is the sequenced used for
analyzing Ds integrations. Only use sequence that you know will be present in
every fragment. If using an element with a variable 3' end (like AAV), this
field should only contain the stable portion of the element.
* `(integration_footprint_size)`: The number of bases directly affected by the
integration. For example, MLV integrates into sets of four basepairs, creating
a 4 bp duplication on either side of the virus. See Fig. 1 of
http://dx.doi.org/10.1371/journal.ppat.0020060 for more examples. If your
element integrates between bases without creating a duplication, set this value
to `1`.


Citation
--------

If you use GeIST in your work, please use the following reference to cite it:

LaFave MC, Varshney GK, Gildea DE, Wolfsberg TG, Baxevanis AD, Burgess SM. MLV integration site selection is driven by strong enhancers and active promoters. 2014 Apr;42(7):4257-69.


Contact
-------

Matthew C. LaFave, Ph.D.,
Developmental Genomics Section, Translational and Functional Genomics Branch,
NHGRI, NIH

Email: matthew.lafave [at sign] nih.gov
