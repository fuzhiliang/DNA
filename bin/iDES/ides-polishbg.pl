#!/usr/bin/perl

#DESCRIPTION==================================================================
# Polish variants in FREQ file using ides-bgdb.txt
# Keep duplex supported variants if corresponding duplex FREQ files available

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
# Unix
# R
# Perl with the following external libraries: Statistics::R, Proc::Fork
use Statistics::R;
use Proc::Fork;
use List::Util qw(max min);
use POSIX qw(:sys_wait_h);
use POSIX;
use Cwd 'abs_path';
use File::Basename;
use Getopt::Std;
#==============================================================================

my $version = "1.0";

my $size = @ARGV;

my %opts = (f=>0.2, n=>2, p=>0.05, t=>1, m=>4, w=>5, a=>5, r=>10);
my %opts1;
getopts('o:f:n:p:t:w:m:a:r:d:b', \%opts1);

die("
ides-polishbg.pl version $version by Aaron M. Newman (amnewman\@stanford.edu)

Purpose:
Perform background polishing on FREQ file(s).

Usage:
ides-polishbg.pl [options] <freq(s)> <ides-bgdb.txt>

    <freq(s)>        Sorted input FREQ file or directory of FREQ files
    <ides-bgdb.txt>  Background database (output of ides-makedb.pl)

Options (defaults in parentheses):
    -o <dir>   output directory (same as input file directory)
    -f <0-1>   minimum fraction of non-zero background samples needed for polishing ($opts{f})
    -n <int>   minimum number of non-zero background samples needed for polishing ($opts{n})
    -m <int>   minimum number of total background samples ($opts{m})
    -w <int>   minimum number of non-zero samples needed for Weibull modeling ($opts{w})
    -a <0-100> maximum allele frequency cutoff for polishing ($opts{a})
    -r <int>   maximum number of supporting reads cutoff for polishing ($opts{r})
    -p <0-1>   nominal p-value threshold for background polishing ($opts{p})
    -t <int>   number of CPUs to use for processing >1 input FREQ file ($opts{t})
    -b         do background polishing on previously polished FREQ file(s).
    -d <dir>   directory of duplex-supported FREQ files
               (note: file names can contain \".duplex.\" tag, but must otherwise
                match input FREQ file names to ensure proper pairing)
    -d <file>  matching duplex-supported FREQ file

For detailed help, see http://cappseq.stanford.edu/ides

") if($size == 0);

die("Less than 2 arguments detected. Abort!\n") if(@ARGV < 2);
my $INPUT = $ARGV[0]; #single FREQ file or directory of FREQ files
if(!(-e $INPUT)) {die("Input FREQ file/directory not found. Abort!\n")};
my $DB = $ARGV[1]; #snvbg.txt (output of snvbg.pl)
if(!(-f $DB)) {die("Background database not found. Abort!\n")};

my $options='';
foreach my $opt(keys %opts1){
    $options.='|'.$opt.$opts1{$opt};
    $opts{$opt}=$opts1{$opt};
    if($opt eq "o") {
        die("\'$opts1{$opt}\' is not a directory. Abort!\n") if(!(-d $opts1{$opt}));
        print "Output directory: $opts1{$opt}\n";
    }
    if($opt eq "n") {
        die("Minimum number of non-zero background samples needed for polishing must be an integer. Abort!\n") if is_integer($opts1{$opt}) == 0;
        die("Minimum number of non-zero background samples needed for polishing must be >0. Abort!\n") if $opts1{$opt}<1;
        print "Minimum number of non-zero background samples needed for polishing set to: $opts1{$opt}\n";
    }
    if($opt eq "m") {
        die("Minimum number of total background samples must be an integer. Abort!\n") if is_integer($opts1{$opt}) == 0;
        die("Minimum number of total background samples must be >3. Abort!\n") if $opts1{$opt}<4;
        print "Minimum number of total background samples set to: $opts1{$opt}\n";
    }
    if($opt eq "w") {
        die("Minimum number of non-zero samples needed for Weibull modeling must be an integer. Abort!\n") if is_integer($opts1{$opt}) == 0;
        die("Minimum number of non-zero samples needed for Weibull modeling must be >3. Abort!\n") if $opts1{$opt}<4;
        print "Minimum number of non-zero samples needed for Weibull modeling set to: $opts1{$opt}\n";
    }
    if($opt eq "r") {
        die("Maximum number of supporting reads cutoff for polishing must be an integer. Abort!\n") if is_integer($opts1{$opt}) == 0;
        die("Maximum number of supporting reads cutoff for polishing must be >0. Abort!\n") if $opts1{$opt}<1;
        print "Maximum number of supporting reads cutoff for polishing set to: $opts1{$opt}\n";
    }
    if($opt eq "t") {
        die("Number of CPUs must be an integer. Abort!\n") if is_integer($opts1{$opt}) == 0;
        die("Number of CPUs must be >0. Abort!\n") if $opts1{$opt}<1;
        print "Number of CPUs set to: $opts1{$opt}\n";
    }
    if($opt eq "p") {
        die("Nominal p-value threshold for background polishing must be between 0 and 1. Abort!\n") if($opts1{$opt}<0 || $opts1{$opt}>1);
        print "Nominal p-value threshold for background polishing set to: $opts1{$opt}\n";
    }
    if($opt eq "f") {
        die("Minimum fraction of non-zero background samples needed for polishing must be between 0 and 1. Abort!\n") if($opts1{$opt}<0 || $opts1{$opt}>1);
        print "Minimum fraction of non-zero background samples needed for polishing set to: $opts1{$opt}\n";
    }
    if($opt eq "a") {
        die("Maximum allele frequency cutoff for polishing must be between 0 and 100. Abort!\n") if($opts1{$opt}<0 || $opts1{$opt}>100);
        print "Maximum allele frequency cutoff for polishing set to: $opts1{$opt}\n";
    }
    if($opt eq "b") {
        print "Run on background polished FREQ files (*rmbg.txt)\n";
    }
    if($opt eq "d") {
        die("\'$opts1{$opt}\' is not a directory or file. Abort!\n") if(!(-d $opts1{$opt}) && !(-f $opts1{$opt}));
        print "Duplex file/directory: $opts1{$opt}\n";
    }
    
}

#PARAMETERS===================================================================

my $MINFRAC = $opts{f}; #0.2; #remove variants present in at least $frac samples in the snvbg database
my $MINN = $opts{n}; #1; #minimum number of samples with a given variant
my $WPVAL = $opts{p}; #nominal pvalue for Weibull distribution and Z-score cutoff (will be adjusted with Bonferroni-correction)
my $FORKS = $opts{t};
my $OUTPUTDIR = dirname($INPUT); #directory to write output database
if(exists($opts{o})) {$OUTPUTDIR = $opts{o};}
my $MINSAMPLES = $opts{m}; #absolute minimum number of samples (with or without a variant)
my $MINWEIBULL = $opts{w}; #minimum number of samples with a given variant needed for Weibull distribution (otherwise Z test used)
my $MAXAFCUTOFF = $opts{a}; #maximum allele frequency cutoff for polishing
my $MAXREADSCUTOFF = $opts{r}; #maximum number of supporting reads cutoff for polishing
my $RUNBG = 0; #if 1, processess bg-polished freq files; if 0, original freq files
if(exists($opts{b})) {$RUNBG = 1;} #Run on background polished FREQ files
my $DUPLEX = dirname($INPUT) . "/duplex"; #directory to write output database
if(exists($opts{d})) {$DUPLEX = $opts{d};}
#=============================================================================


#LOAD BACKGROUND==============================================================
my %bg = ();
my $countbg = 0;
open(FILE, $DB) or die $!;
while(<FILE>){
    my $line = $_;
    chomp($line);
    next if $line =~ m/Pos/;
    my @tokens = split("\t", $line);
    my $chr = $tokens[0];
    my $pos = $tokens[1];
    my $var = $tokens[3];
    my $frac = $tokens[6];
    my $numsamples = $tokens[5];
    my $numerator = $tokens[4];
    next if $numsamples < $MINSAMPLES;
    my $meanreads = $tokens[8];
    my $stdreads = $tokens[10];
    my $meanAF = $tokens[11];
    my $stdAF = $tokens[13];
    my $shape = $tokens[14];
    my $scale = $tokens[15];
    my @data = ($meanreads, $stdreads, $meanAF, $stdAF, $frac, $shape, $scale, $numerator);
    if($frac >= $MINFRAC && $numerator >= $MINN) {$bg{"$chr:$pos:$var"} = \@data; $countbg++;}
}
close(FILE);
#---------------------------------------

#P-value cutoffs with Bonferroni correction
my $R = Statistics::R->new();
$R->start();
$R->set('z', $WPVAL / $countbg);
my $zscore = $R->get('qnorm(1 - (z))');
$WPVAL /= $countbg; #BC
$R->stop();

#read in FREQ file(s)

my @commands = ();

if(-d $INPUT){ #if input is directory, read FREQ files
    
    opendir (DIR, $INPUT) or die $!;
    while (my $file = readdir(DIR)){
        chomp($file);
        next if !($file =~ m/freq.(allreads|paired).Q.{1,2}/);
        if($RUNBG == 0){
            next if $file =~ m/rmbg.txt/;
        }else{
            next if !($file =~ m/rmbg.txt/);
        }
        push(@commands, "$INPUT/$file");
    }
    close(DIR);
}else{#else input is single file
    next if !($INPUT =~ m/freq.(allreads|paired).Q.{1,2}/);
    if($RUNBG == 0){
        next if $INPUT =~ m/rmbg.txt/;
    }else{
        next if !($INPUT =~ m/rmbg.txt/);
    }
    push(@commands, $INPUT);
}

forkpool(\@commands, $FORKS);


sub runfilter{

    my $R = Statistics::R->new();
    $R->start();
    
    my $freq = $_[0];
    
    print ">Processing $freq...\n";

    my %duplex = (); #store duplex data if they exist                                                                                 
    if(exists($opts{d})) {%duplex =  %{getduplexsupport($DUPLEX,$freq)};}

    open(FILE, $freq) or die $!;
    my $output = basename($freq);
    $output =~ s/txt/rmbg.txt/g;
    open(out, ">$OUTPUTDIR/$output") or die $!;

    while(<FILE>){
        my $line = $_;
        chomp($line);
        if($line =~ m/DEPTH/) {print out "$line\n"; next;}
        my @tokens = split("\t", $line);
        my $chr = $tokens[0];
        my $pos = $tokens[1];
        my $depth = $tokens[2];
        next if $depth == 0;
        my $ref = uc($tokens[3]);
        $totalbases += $depth;
        my $totalerrors = 0;
        my $outa = "";
        my $outb = "";
        my $minus = 0;
        for(my $i = 6; $i < @tokens-1; $i+=2){
            my $var = getVar($i);
                #if $ref eq $var;
            my $forw = $tokens[$i];
            my $rev = $tokens[$i+1];
            my $support = $forw + $rev;
            my $AF = (100*$support/$depth);
            if(defined($bg{"$chr:$pos:$var"})){
                my ($meanreads, $stdreads, $meanAF, $stdAF, $frac, $shape, $scale, $numsamples) = @{$bg{"$chr:$pos:$var"}};
                if($stdreads == 0){
                    $outa = $outa . "\t$tokens[$i]\t$tokens[$i+1]";
                    next;
                }
                my $readsZ = ($support - $meanreads) / $stdreads;
                my $AFZ = ($AF - $meanAF) / $stdAF;
                
                my $WP = 1;
                
                if($scale ne "NA" && $scale ne "" && $numsamples >= $MINWEIBULL){
                    $R->set('x', 100*$AF);
                    $R->set('shape', $shape);
                    $R->set('scale', $scale);
                    $WP = 1 - ((1-$frac) + ($frac * $R->get('pweibull(x,shape,scale,TRUE)')));
                }

                if((($AFZ < $zscore && ($numsamples < $MINWEIBULL ||$scale eq "NA" || $scale eq "")) ||
                ($scale ne "NA" && $scale ne "" && $numsamples >= $MINWEIBULL && $WP > $WPVAL)) && ($AF < $MAXAFCUTOFF || $support <= $MAXREADSCUTOFF)) {
                
                    #check for duplex support; if found, don't polish
                    if(defined $duplex{"$chr:$pos"}->{"$ref:$var"}){
                    }else{
                        $minus += $support;
                        $outa = $outa . "\t0\t0"; $outb = $outb .  "\t$tokens[$i]\t$tokens[$i+1]";
                        next;
                    }
                }
            }
            $outa = $outa . "\t$tokens[$i]\t$tokens[$i+1]";
            $totalerrors += $support;
        }
        print out "$chr\t$pos\t" . ($depth-$minus) . "\t$ref\t$tokens[4]\t$tokens[5]" . $outa . "\n";
    }
    close(FILE);
    close(out);
    $R->stop();
    
}

if(-e "Rplots.pdf") {system("rm Rplots.pdf");}

sub getVar{
    my $i = $_[0];
    if($i == 6) {return 'A';}
    if($i == 8) {return 'C';}
    if($i == 10) {return 'T';}
    if($i == 12) {return 'G';}
    return 'N';
}

#read duplex data for a given sample
sub getduplexsupport{
    my $duplexinput = $_[0];
    my $match = $_[1];
    $match =~ s/^\.\///g;
    my %duplex = ();
    my $sample = $match;

    $sample =~ s/^\.\///g;
    my $duplexfile = "";

    if(-d $duplexinput){ #if directory
        open my $dup, '-|', "ls $duplexinput" or die "Failed to open pipe$!";
        my $dupfile = "";
        while(<$dup>){
            my $line = $_;
            chomp($line);
            my $dfile = $line;
            $dfile =~ s/duplex.//g;
            if($dfile eq $match) {
                $duplexfile = "$duplexinput/$line";
                last;
            }
        }
        close($duplexinput);
    }else{
        $duplexfile = $duplexinput;
                print "$match\t$duplexfile\n";
    }
    
    
    open(FILE, $duplexfile) or die $!;
    while(<FILE>){
        my $line = $_;
        chomp($line);
	next if $line =~ m/DEPTH/;
        my @tokens = split("\t", "$line");
	my $chr = $tokens[0];
        my $pos = $tokens[1];
        my $depth = $tokens[2];
        next if $depth == 0;
        my $ref = uc($tokens[3]);
	for(my $i = 6; $i < @tokens-1; $i+=2){
            my $var = getVar($i);
	    my $forw = $tokens[$i];
            my $rev = $tokens[$i+1];
            my $support = $forw + $rev;
	    next if $support == 0 || (100*$support/$depth) > $MAXAFCUTOFF;
	    $duplex{"$chr:$pos"}->{"$ref:$var"} = $support;
	}
    }
    close(FILE);
    return \%duplex;
}

#run each sample on a separate process using a fork pool--------------------------------------------;
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
                runfilter($cmds[$itor]);
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
