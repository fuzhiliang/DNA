#!/usr/bin/perl

#DESCRIPTION==================================================================
# Collect error statistics from FREQ file

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
use List::Util qw(max min);
use Getopt::Std;
#=============================================================================


my $version = "1.0";

my $size = @ARGV;

my %opts = (f=>5, d=>200);
my %opts1;
getopts('d:f:', \%opts1);

die("
ides-bgreport.pl version $version by Aaron M. Newman (amnewman\@stanford.edu)

Purpose:
Collect background statistics from FREQ file.

Usage:
ides-polishbg.pl [options] <freq>

    <freq>   FREQ file

Options (defaults in parentheses):
    -f <0-100> Maximum AF% to consider ($opts{f}%)
    -d <int>   Minimum depth to consider ($opts{d})

For detailed help, see http://cappseq.stanford.edu/ides

") if($size == 0);

die("Input FREQ file not detected. Abort!\n") if(@ARGV < 1);

my $INPUT = $ARGV[0]; #.freq file input

my $options='';
foreach my $opt(keys %opts1){
    $options.='|'.$opt.$opts1{$opt};
    $opts{$opt}=$opts1{$opt};
    if($opt eq "f") {
        die("Maximum AF% to consider must be between 0 and 100. Abort!\n") if($opts1{$opt}<0 || $opts1{$opt}>100);
        print "Maximum AF% to consider set to: $opts1{$opt}\n";
    }
    if($opt eq "d") {
        die("Minimum depth to consider must be an integer. Abort!\n") if is_integer($opts1{$opt}) == 0;
        die("Minimum depth to consider must be >0. Abort!\n") if $opts1{$opt}<1;
        print "Minimum depth to consider set to: $opts1{$opt}\n";
    }
}

#PARAMETERS===================================================================
my $MAXAF = $opts{f}; # Maximum AF (%) to consider
my $MIND = $opts{d}; # Minimum depth to consider
#=============================================================================



open(FILE, $INPUT) or die $!;

my $zeros = 0;
my $bases = 0;
my $positions = 0;
my $all_errors = 0;
my %base_errors = ();

while(<FILE>){
    my $line = $_;
    chomp($line);
    my @tokens = split("\t", $line);
    my $depth = $tokens[2];
    my $chrpos = "$tokens[0]:$tokens[1]";
    next if $depth < $MIND;
    my $nuc = uc($tokens[3]);
    my $As = $tokens[6] + $tokens[7];
    my $Cs = $tokens[8] + $tokens[9];
    my $Ts = $tokens[10] + $tokens[11];
    my $Gs = $tokens[12] + $tokens[13];
    my $mx = max($As, $Cs, $Ts, $Gs);
    next if (100*$mx/$depth) > $MAXAF;
    $positions++;
    if($mx == 0) {$zeros++;}
    $bases += $depth;
    $all_errors += ($As + $Cs + $Ts + $Gs);
    my $perc_error = 100*($As + $Cs + $Ts + $Gs)/$depth;

    if($As > 0) {$base_errors{"$nuc>A"}->{$chrpos} += $As;}
    if($Cs > 0) {$base_errors{"$nuc>C"}->{$chrpos} += $Cs;}
    if($Ts > 0) {$base_errors{"$nuc>T"}->{$chrpos} += $Ts;}
    if($Gs > 0) {$base_errors{"$nuc>G"}->{$chrpos} += $Gs;}
}

die("There are no bases that satisfy current parameters. Abort!\n") if $positions == 0;

print "No. positions:\t$positions\n";
print "No. positions without errors:\t$zeros\n";
print "Percent positions without errors:\t" . (100*$zeros/$positions) . "\n";
print "No. bases sequenced:\t$bases\n";
print "No. errors:\t$all_errors\n";
print "Percent errors:\t" . (100*$all_errors/$bases)  . "\n";

print "Subst.\tPositions\tErrors\t%Errors\n";
foreach my $subs (sort keys %base_errors){
    my $count = keys %{$base_errors{$subs}};
    my $basecount = 0;
    foreach my $chrpos (keys %{$base_errors{$subs}}){
	$basecount += $base_errors{$subs}->{$chrpos};
    }
    print "$subs:\t$count\t$basecount\t" . (100*$basecount/$all_errors) . "\n";
}


sub is_integer {
    local $_ = shift;
    # checks from perlfaq4
    return $] < 5.009002 unless defined;
    return 1 if (/^[+-]?\d+$/); # is a +/- integer
    0;
}
