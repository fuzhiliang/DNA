#!/usr/bin/perl

#DESCRIPTION==================================================================
# Calculate background database for polishing.

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
# Perl and the following external libraries:
# Statistics::Descriptive, Statistics::R
# R and the following external library: 'fitdistrplus'
use File::Spec;
use Getopt::Std;
use File::Basename;
use POSIX;
use Statistics::Descriptive;
use Statistics::R;
#=============================================================================

my $version = "1.1";

my $size = @ARGV;

my %opts = (m=>100, d=>20, n=>6);
my %opts1;
getopts('m:d:o:n:a:p', \%opts1);

die("
ides-makedb.pl version $version by Aaron M. Newman (amnewman\@stanford.edu)

Purpose:
Create a SNV background database from a set of sequencing libraries (FREQ file format).

Usage:
ides-makedb.pl [options] <DIR>

    <DIR>   Directory of FREQ files

Options (defaults in parentheses):
    -o <dir>   output directory (same as input FREQ file directory)
    -m <0-100> maximum allele frequency for training database ($opts{m})
    -d <int>   minimum total depth required for each genomic position ($opts{d})
    -n <int>   minimum number of input samples required ($opts{n})
    -a <str>   add custom name to background database
    -p         print Weibull parameter estimation errors

For detailed help, see http://cappseq.stanford.edu/ides

") if($size == 0);

my $dir = $ARGV[0]; #directory of FREQ files to build background database.

die("Not a directory. Abort!\n") if(!(-d $dir));

my $options='';
foreach my $opt(keys %opts1){
    $options.='|'.$opt.$opts1{$opt};
    $opts{$opt}=$opts1{$opt};
    if($opt eq "o") {
        die("\'$opts1{$opt}\' is not a directory. Abort!\n") if(!(-d $opts1{$opt}));
        print "Output directory: $opts1{$opt}\n";
    }
    if($opt eq "m") {
        die("Maximum allele frequency must be between 1 and 100. Abort!\n") if($opts1{$opt}<0 || $opts1{$opt}>100);
        print "Maximum allele frequency set to: $opts1{$opt}\n";
    }
    if($opt eq "d") {
        die("Minimum depth must be an integer. Abort!\n") if is_integer($opts1{$opt}) == 0;
        die("Minimum depth must be >0. Abort!\n") if $opts1{$opt}<1;
        print "Minimum depth for each genomic position set to: $opts1{$opt}\n";
    }
    if($opt eq "n") {
        die("Minimum number of samples must be an integer. Abort!\n") if is_integer($opts1{$opt}) == 0;
        die("Minimum number of samples must be >3 (>10 recommended). Abort!\n") if $opts1{$opt}<4;
        print "Minimum number of samples set to: $opts1{$opt}\n";
    }
    if($opt eq "a") {
        print "Set custom name to: $opts1{$opt}\n";
    }
    if($opt eq "p") {
        print "Print Weibull parameter estimation errors.\n";
    }
}

#PARAMETERS===================================================================
my $MAXAF = $opts{m}; #Note percentage, not fraction; only consider bg events with AF < $maxAF
my $MINDEPTH = $opts{d}; #Only consider genomic positions with depth >= $mindepth
my $MINSAMPLES = $opts{n}; #Change minimum number of input samples required
my $NAME = "";
if(exists($opts{a})) {$NAME = $opts{a} . "_";} #Add custom name to output database
my $OUTPUTDIR = dirname($dir); #directory to write output database
if(exists($opts{o})) {$OUTPUTDIR = $opts{o};}
my $PRINTERRORS = 0; #if 1, print Weibull parameter estimation errors
if(exists($opts{p})) {$PRINTERRORS = 1;}
my $R = Statistics::R->new();
$R->startR;
$R->run(q`library(fitdistrplus)`);
#=============================================================================

my $inputcount = 0;
my %store = (); #store FREQ files

opendir(DIR, $dir) or die$!;

#Read and store input FREQ files
while(my $file = readdir(DIR)){
    next if !($file =~ m/freq.(allreads|paired).Q.{1,2}.txt/); #check for freq files
    $inputcount++;
    print "Loading $dir/$file...\n";
    my $sample = $file;
    $sample =~ s/Sample_//g;
    $sample =~ s/.sorted.+//g;
    open(FILE, "$dir/$file") or die$!;

    while(<FILE>){
        my $line = $_;
        chomp($line);
        next if $line =~ m/DEPTH/;
        my @tokens = split("\t", $line);
        my $chr = $tokens[0];
        my $pos = $tokens[1];
        for(my $i = 6; $i < @tokens; $i+=2){
            my $varsupport = $tokens[$i] + $tokens[$i+1];
            my $var = getVar($i);
            my $key = "$chr:$pos:$var";
            $store{$key}->{$sample} = \@tokens;
        }
    }
}

die("$inputcount of $MINSAMPLES required FREQ files detected. Abort!") if($inputcount < $MINSAMPLES);

print "\nGenerating background database...\n";

open(output, ">$OUTPUTDIR/" . $NAME . "ides-bgdb.txt") or die $!;

print output "Chr\tPos\tRef\tVar\tNumPosSamples\tTotalSamples\tFracSamples\tFracBothStrands\tMeanReads\tMedianReads\tStdReads\tMeanAF\tMedianAF\tStdAF\tW_Shape\tW_Scale\tW_Corr\tW_Pval\n";

my @failed = (); #collect any position-bp pairs that fail Weibull parameter estimation

my $progress_itor = 0;
my $prev_prog = -1;
my $countasterisk = 0;
my $total_comp = keys %store; #count total events to be compared (for progressbar)
print "0%   20   40   60   80   100%\n";
print "|----|----|----|----|----|\n";


#Collect statistics
foreach my $key (keys %store){
    my $depth_stat = Statistics::Descriptive::Full->new();
    my $reads_stat = Statistics::Descriptive::Full->new();
    my $num_samples = 0;
    my $num_bothstrands = 0;
    my @tokens = split(":", $key);
    my $chr = $tokens[0];
    my $pos = $tokens[1];
    my $var = $tokens[2];
    my $ref = "";
    my $countfiles = 0;
    foreach my $sample (keys %{$store{$key}}){
        my @data = @{$store{$key}->{$sample}};
        $ref = uc($data[3]);
        my $total_depth = $data[2];
        my $index = getIndex($var);
        next if $index == -1;
        my $var_depth = $data[$index]+$data[$index+1];
        next if $total_depth < $MINDEPTH;
        my $AF = 100 * ($var_depth / $total_depth);
        next if $AF >= $MAXAF;
        $countfiles++;
        $depth_stat->add_data($AF);
        $reads_stat->add_data($var_depth);
        if($data[$index]>0 && $data[$index+1]>0) {$num_bothstrands++;}
        if($AF > 0) {$num_samples++;}
    }
    if($countfiles < 2) {$progress_itor++; next;};
    
    #remove maximum AF per position-bp pair (outlier removal)
    my $buffer = 0;
    if($countfiles > 2){
        $num_samples--;
        $countfiles--;
        $buffer = 1;
    }
    my $frac_samples = ($num_samples/$countfiles);
    my $frac_bothstrands = ($num_bothstrands/$countfiles);
    $depth_stat->sort_data();
    $reads_stat->sort_data();
    my $depth_stat_trimmax = Statistics::Descriptive::Full->new();
    my $reads_stat_trimmax = Statistics::Descriptive::Full->new();
    my @depthstats = $depth_stat->get_data();
    my @readstats = $reads_stat->get_data();
    for(my $j = 0; $j < @depthstats-$buffer; $j++){
        $depth_stat_trimmax->add_data($depthstats[$j]);
        $reads_stat_trimmax->add_data($readstats[$j]);
    }
    
    my $ave_AF = $depth_stat_trimmax->mean();
    if($ave_AF == 0) {$progress_itor++; next;} #skip current position-bp pair if no background
    
    #Weibull parameters
    my $shape = "NA";
    my $scale = "NA";
    my $cor = "NA";
    my $cpval = "NA";
    my @temp = $depth_stat_trimmax->get_data();
    $R->set('x', \@temp);
    $R->run(q`x <- x[which(x>0)]`);
    my $len = $R->get('length(unique(x))');
    if($len>2){
        $R->run(q`w <- NA`);
        $R->run(q`tryCatch({w<-suppressWarnings(summary(fitdist(x*100, "weibull"))$estimate)}, warning = function(war){ b<-1 }, error = function(err){ b<-1 }, finally = { print(w) })`);
        if($R->get('w[1]') ne "NA"){ #if parameter estimation successful
            $shape = $R->get('w[1]');
            $scale = $R->get('w[2]');
            $R->run(q`w <- qqplot(x*100, rweibull(100, w[1], w[2]))`);
            $R->run(q`co <- cor.test(w$x, w$y)`);
            $cor = $R->get('co$estimate');
            $cpval = $R->get('co$p.value');
        }else{#if parameter estimation unsuccessful
             push(@failed, "$chr:$pos:$ref>$var");
        }
    }
    my $median_AF =  $depth_stat_trimmax->median();
    my $std_AF = $depth_stat_trimmax->standard_deviation();
    my $ave_reads = $reads_stat_trimmax->mean();
    my $median_reads = $reads_stat_trimmax->median();
    my $std_reads = $reads_stat_trimmax->standard_deviation();
    print output "$chr\t$pos\t$ref\t$var\t$num_samples\t$countfiles\t$frac_samples\t$frac_bothstrands\t$ave_reads\t$median_reads\t$std_reads\t$ave_AF\t$median_AF\t$std_AF\t$shape\t$scale\t$cor\t$cpval\n";
    
    $progress_itor++;
    my $prog = int(100*$progress_itor/$total_comp);
    if($prog % 4 == 0 && $prog != $prev_prog){
        $| = 1;
        if($countasterisk < 26){print "*";}
        $countasterisk++;
        $prev_prog = $prog;
    }
}

while($countasterisk < 26){
    print "*";
    $countasterisk++;
}

close(output);

if(@failed > 0 && $PRINTERRORS == 1){
    print "\n" . @failed . " Weibull parameter estimation error(s):\n";
    foreach(@failed){
        my $line = $_;
        print "$line\n";
    }
}

print "\n\nBackground database complete! Written to: $OUTPUTDIR/" . $NAME . "ides-bgdb.txt\n\n";

if(-e "Rplots.pdf") {system("rm Rplots.pdf");}

sub getVar{
    my $i = $_[0];
    if($i == 6) {return 'A';}
    if($i == 8) {return 'C';}
    if($i == 10) {return 'T';}
    if($i == 12) {return 'G';}
    return 'N';
}

sub getIndex{
    my $var = $_[0];
    if($var eq "A") {return 6;}
    if($var eq "C") {return 8;}
    if($var eq "T") {return 10;}
    if($var eq "G") {return 12;}
    return -1;
}

sub is_integer {
    local $_ = shift;
    # checks from perlfaq4
    return $] < 5.009002 unless defined;
    return 1 if (/^[+-]?\d+$/); # is a +/- integer
    0;
}
