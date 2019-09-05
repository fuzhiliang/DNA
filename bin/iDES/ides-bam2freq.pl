#!/usr/bin/perl

#DESCRIPTION==================================================================
# Read in bam file(s) (single BAM or directory containing BAMs), genome (FASTA), and targeted regions (BED), and write FREQ file(s) containing position-specific nucleotide counts to disk

#LICENSE AND TERMS OF USE=====================================================
# Licensed under Docket S14-399. All rights are covered under the license.

# THE BOARD OF TRUSTEES OF THE LELAND STANFORD JUNIOR UNIVERSITY (“STANFORD”) provides iDES software (“Program”), including any accompanying information, materials, or manuals, free of charge for non-commercial use only and bound by the above license agreement. By accepting, receiving, or using the Program, you (“RECIPIENT”) agree to be bound by the terms of this agreement (“Agreement”). If you do not agree to the terms of this Agreement, then do not use the Program and promptly remove all copies of the program from your computer(s).

# RECIPIENT acknowledges that the Program is a research tool still in the development stage and that it is being supplied as is, without any accompanying services, support, or improvements from STANFORD. STANFORD makes no representations and extends no warranties of any kind, neither express nor implied.

# RECIPIENT shall not use the Program on behalf of any organization that is not a non-profit organization. RECIPIENT shall not use the Program for commercial advantage, or in the course of for-profit activities. This Agreement provides no license, express or implied, under any patent (i) to any organization that is not a non-profit research organization; (ii) to use the Program on behalf of any for-profit entity; or (iii) to use the Program for any commercial purpose. The use of the Program by, or on behalf of, any for-profit entity will terminate any and all licenses, express or implied, granted by STANFORD under this Agreement. Commercial entities wishing to use this Program should contact Stanford University’s Office of Technology Licensing and reference docket S14-399.

# RECIPIENT shall not distribute the Program or transfer it to any other person or organization without prior written permission from STANFORD or licensee as per the above docket.

# RECIPIENT shall not reverse engineer, reverse assemble, reverse compile, decompile, disassemble, or otherwise attempt to create the source code for the Program. RECIPIENT acknowledges that any programs created based on the Program will be considered a derivative of the Program and owned by STANFORD or licensee as per the above docket.

# RECIPIENT shall NOT make modifications to the Program. Title and copyright to the Program, any derivatives thereof, and any associated documentation will at all times remain with STANFORD or licensee as per the above docket, and RECIPIENT agrees to preserve same.

# RECIPIENT warrants that RECIPIENT will not remove or export any part of the Program from the United States except in full compliance with all United States and other applicable laws and regulations.

# RECIPIENT shall indemnify, hold harmless, and defend STANFORD against any claim of any kind arising out of, or related to, (i) the exercise of any rights granted under this Agreement, or (ii) the breach of this Agreement by RECIPIENT.

#DEPENDENCIES=================================================================
# SAMtools 1.2+ (samtools executable must be in PATH; e.g., /usr/bin)
# Perl and the following external library: Proc::Fork

use Proc::Fork;
use Getopt::Std;
use List::Util qw(max min);
use POSIX qw(:sys_wait_h);
use Cwd 'abs_path';
use File::Basename;
#==============================================================================

my $version = "1.0";

my $size = @ARGV;

my %opts = (q=>30, t=>1);
my %opts1;
getopts('o:q:t:ab', \%opts1);

die("
ides-bam2freq.pl version $version by Aaron M. Newman (amnewman\@stanford.edu)

Purpose:
Convert BAM file(s) to FREQ file(s).

Usage:
ides-bam2freq.pl [options] <bam(s)> <genome.fa> <targets.bed>

    <bam(s)>        Sorted input BAM file or directory of sorted BAM files
                    (*sorted.bam file extension required)
    <genome.fa>     Reference genome (FASTA format)
    <targets.bed>   Restrict output to targeted regions (BED format)

Options (defaults in parentheses):
    -o <dir>   output directory (same as input file directory)
    -q <int>   Phred quality score filter ($opts{q})
    -t <int>   number of CPUs to use for processing >1 input BAM ($opts{t})
    -a         output all reads (output properly paired reads only)
    -b         disable mpileup base alignment quality adjustment

For detailed help, see http://cappseq.stanford.edu/ides

") if($size == 0);

die("Less than 3 arguments detected. Abort!\n") if(@ARGV < 3);
my $INPUT = $ARGV[0]; #directory of FREQ files to build background database.
if(!(-e $INPUT)) {die("Input BAM file/directory not found. Abort!\n")};
my $GENOME = $ARGV[1]; #genome sequence in fasta format
my $TARGETS = $ARGV[2]; #bed coordinates to focus output

my $options='';
foreach my $opt(keys %opts1){
    $options.='|'.$opt.$opts1{$opt};
    $opts{$opt}=$opts1{$opt};
    if($opt eq "o") {
        die("\'$opts1{$opt}\' is not a directory. Abort!\n") if(!(-d $opts1{$opt}));
        print "Output directory: $opts1{$opt}\n";
    }
    if($opt eq "q") {
        die("Minimum Phred quality must be an integer. Abort!\n") if is_integer($opts1{$opt}) == 0;
        die("Minimum Phred quality must be >=0. Abort!\n") if $opts1{$opt}<0;
        print "Minimum Phred quality set to: $opts1{$opt}\n";
    }
    if($opt eq "t") {
        die("Number of CPUs must be an integer. Abort!\n") if is_integer($opts1{$opt}) == 0;
        die("Number of CPUs must be >0. Abort!\n") if $opts1{$opt}<1;
        print "Number of CPUs set to: $opts1{$opt}\n";
    }
    if($opt eq "a") {
        print "Disable requirement for properly paired reads (all reads allowed).\n";
    }
    if($opt eq "b") {
        print "Disable BAQ.\n";
    }
}

#PARAMETERS===================================================================

my $DOPAIRED = 1; #if 1, restrict depth/on-target statistics to properly paired reads
if(exists($opts{a})) {$DOPAIRED = 0;} #if 0 examine all reads for depth statistics; default is properly paired
my $BAQ = 1; #if 1, do BAQ
if(exists($opts{b})) {$BAQ = 0;} #if 0, disable BAQ
my $QUALITY = $opts{q}; # Phred quality cutoff (Q30 recommended)
my $FORKS = $opts{t};
my $OUTPUTDIR = dirname($INPUT); #directory to write output database
if(exists($opts{o})) {$OUTPUTDIR = $opts{o};}
#=============================================================================




my @files = (); #store BAM files

if(-d $INPUT){ #if input is directory, read BAM
    my $ABSP = abs_path($INPUT);
    opendir(DIR, $ABSP) or die $!;
    while (my $file = readdir(DIR)){
        chomp($file);
        next if !($file =~ m/sorted.bam$/); #check for sorted BAM file
        push(@files, "$ABSP/$file");
    }
}else{#else input is single file
    if($INPUT =~ m/sorted.bam$/){ #check for sorted BAM file
        push(@files, $INPUT);
    }
}

die("No sorted BAM files (*sorted.bam) detected. Abort!") if(length(@files) == 0);

forkpool(\@files, $FORKS);

sub run{

    my $file = $_[0];
    print ">Reading $file\n";
    my $insert = "-f 3";
    my $insert2 = "";
    if($DOPAIRED == 0) {$insert = ""; $insert2 = "-A";}
    my $insert3 = "";
    if($BAQ == 0) {$insert3 = "-B";}
    
    my $mpileup;

    open $mpileup, '-|', "samtools view -hb $insert $file | samtools mpileup -Q$QUALITY $insert2 $insert3 -d10000000 -f $GENOME -l $TARGETS -" or die "Failed to open pipe$!";

    $file = basename($file);

    if($DOPAIRED == 1) {
        my $replace = "freq.paired.Q$QUALITY";
        $file =~ s/bam/$replace/g;}
    else {
        my $replace = "freq.allreads.Q$QUALITY";
        $file =~ s/bam/$replace/g;
    }

    open(freq_output, ">$OUTPUTDIR/$file.txt");
    
    print freq_output "CHR\tPOS\tDEPTH\tREF\tR+\tR-\tA+\tA-\tC+\tC-\tT+\tT-\tG+\tG-\n";
    
    while(<$mpileup>) {
        my $line = $_;
        if($line =~ m/chr/){
            my @vars = split("\t",$line);
            my $chr = $vars[0];
            my $pos = $vars[1];
            my $ref = $vars[2];
            my $depth = $vars[3];
            my $base_string = $vars[4];
            my @freq = ();
            for (my $key = 0; $key < 10; $key++) {
                $freq[$key] = 0;
            }
	    my $indel = 0; #bit for indels
	    my $read_start = 0; #bit for start of read
            for (my $key = 0; $key < length($base_string); $key++) {
                my $sym = substr ($base_string,$key,1);
		if($indel == 1){
		    if($sym eq '.' || $sym eq ',') {$indel = 0;} 
		}
		next if $indel == 1;
		if($read_start == 1) {$read_start = 0; next;}
		if($sym eq '^') {$read_start = 1; next;}
		if($sym eq '+' || $sym eq '-') {$indel = 1;}
		else{
		    if($sym eq '.') {$freq[0]++;}
		    if($sym eq ',') {$freq[1]++;}
		    if($sym eq 'A') {$freq[2]++;}
		    if($sym eq 'a') {$freq[3]++;}
		    if($sym eq 'C') {$freq[4]++;}
		    if($sym eq 'c') {$freq[5]++;}
		    if($sym eq 'T') {$freq[6]++;}
		    if($sym eq 't') {$freq[7]++;}
		    if($sym eq 'G') {$freq[8]++;}
		    if($sym eq 'g') {$freq[9]++;}
		}
            }
	    my $depth_count = 0;
	    for(my $i = 0; $i < @freq; $i++){$depth_count += $freq[$i];}
            print freq_output "$chr\t$pos\t$depth_count\t$ref\t$freq[0]\t$freq[1]\t$freq[2]\t$freq[3]\t$freq[4]\t$freq[5]\t$freq[6]\t$freq[7]\t$freq[8]\t$freq[9]\n";
        }
    }
    
    close(freq_output);

}

#run each sample on a separate process using a fork pool-----------------------------------;                             

sub forkpool{
    print "Starting fork pool\n";
    my @cmds = @{$_[0]};
    my $nb_compute = @{$_[0]};
    #print "$nb_compute\n";                                                                                              
    my $max = min($_[1], $nb_compute);
    my %pids;
    my $itor = 0;
    print "$nb_compute\n";
    while ($itor < $nb_compute) {

    # max limit reached                                                                                                  
        while ($max == keys %pids) {
        # loop thru pid list until a child is released                                                                   
            for my $pid (keys %pids) {

                if (my $kid = waitpid($pid, WNOHANG)) {
                    delete $pids{ $kid };
                    print "finished child id: $pid\n";
                    last;
                }
            }
        }

        run_fork {
            parent {
                my $child = shift;
                print "start child id: $child\n";
                $pids{ $child } = 1;
            }
            child {
                print $cmds[$itor] . "\t$itor\n";
                run($cmds[$itor]);
		
                exit;
            }

        };
        $itor++;
    }

    while (keys %pids > 0) {
    # loop thru pid list until all children are finished                                                                 
        for my $pid (keys %pids) {

            if (my $kid = waitpid($pid, WNOHANG)) {
                delete $pids{ $kid };
		print "finished child id: $pid\n";
                last;
            }
        }
    }

    print "End of fork pool\n";
}


sub is_integer {
    local $_ = shift;
    # checks from perlfaq4
    return $] < 5.009002 unless defined;
    return 1 if (/^[+-]?\d+$/); # is a +/- integer
    0;
}

