
README
======

For additional details, see:

   http://cappseq.stanford.edu/ides

LICENSE
-------
See LICENSE.txt included in the iDES package.

DESCRIPTION
-----------

iDES (Integrated Digital Error Suppression) is a method for the 
suppression of background artifacts in high-throughput sequencing data.
This package provides Perl implementations and for input file conversion,
background database construction, background polishing, and quality
control statistics.

REQUIREMENTS
------------
Unix operating system (Linux, Mac OS X, etc.).

Reference genome (FASTA format; e.g., hg19.fa or hg38.fa).

R (tested with version 3+), with fitdistrplus library.

Perl 5, with the following external dependencies:
	Statistics::Descriptive, Statistics::R, Proc::Fork.

SAMtools 0.1.20+
	copy/link/move to PATH (i.e., /usr/bin).

iDES STEPS
----------

    1. Format input files (required)
    2. Create background database (required)
    3. Perform background polishing
    4. Background error report (optional)

USAGE
=====

For options and output details, see http://cappseq.stanford.edu/ides

-----------------------------------------------------------
1. Convert input BAM files to frequency (FREQ) file format.
-----------------------------------------------------------

Command

    perl ides-bam2freq.pl [options] MySample.bam(s) genome.fa targets.bed

Input

    MySample.bam(s): Single position-sorted BAM file or directory of 
    position-sorted BAM files consisting of paired-end reads (with or
    without de-duplication). BAM index (BAI) file(s) should be present
    in the same directory.

    genome.fa is a reference genome in FASTA format (e.g., hg19.fa).

    targets.bed is a standard 3-column BED (chr start end) used to 
    restrict FREQ file(s) to genomic regions of interest.

Output

    MySample.[paired/allreads].Q(0-40+).freq.txt = Single FREQ file
    created for each input BAM file.

------------------------------------------------------
2. Create nucleotide substitution background database.
------------------------------------------------------

Command

    perl ides-makedb.pl [options] dir

Input

    dir: Directory of FREQ files for building background database. 

Output

    ides-bgdb.txt = background database.

--------------------------------
3. Perform background polishing.
--------------------------------

Command

    perl ides-polishbg.pl [options] MySample.freq(s) ides-bgdb.txt

Input

    MySample.freq(s): Single FREQ file or directory of FREQ files.

    ides-bgdb.txt: Background database.

Output

    MySample.[paired/allreads].Q(0-40+).freq.rmbg.txt = Background
    polished FREQ file(s).

------------------------------------
4. Generate background error report.
------------------------------------

Command

    perl ides-bgreport.pl [options] MySample.freq

Input

    MySample.freq: Single FREQ file (with or without background polishing).

Output

    STDOUT = background statistics.

