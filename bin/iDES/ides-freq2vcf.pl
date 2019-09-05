#!/usr/bin/perl -w
use strict;
use warnings;
use POSIX qw(strftime);
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin $Script);
use Cwd 'abs_path';
my $verbose ="v1.0";


my ($input,$output,$cutoff);
my %opts;
GetOptions(
    'verbose' => \$verbose,
	'input=s' => \$input,
	'output=s' => \$output, 
	'cutoff=f' => \$cutoff,
    "help|h" =>\&usage,
) or die $!;
unless(defined $input && defined $cutoff && defined $output ){&usage();exit 0;}
open I , $input or die $!;
open O ,">$output" or die $!;
my %hash;
my $head=<I>;
chomp $head;
print O "$head\tALT\tALT_freq\n";
while(<I>){
	chomp;
	my $col||=6;
	my @tmp=split "\t",$_;
	next if $tmp[2]==0;
	foreach my $base("A","C","T","G"){
		$hash{$base}=$tmp[$col]/$tmp[2];
		$col++;
		$hash{$base}+=$tmp[$col]/$tmp[2];
		$hash{$base}=0 if ($tmp[$col]*$tmp[$col-1]==0 && $hash{$base}<=3*$cutoff );## 正反向有为0的reads时，需要突变频率大于cutoff的3倍
		$col++;

	}
	foreach my $key(sort {$hash{$b}<=>$hash{$a}} keys %hash){
		next if  ($hash{$key}<=$cutoff);
		print O "$_\t$key\t$hash{$key}\n";
	}
} 
close I;
close O;


sub usage {
    die(
        qq!
Usage:根据iDES的结果，计算输出突变频率大于cutoff的所有位点。
	eg:
perl $0 -input /home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis_uid/output/v300013477_L03-1-24HD753_single.uid.sorted.freq.paired.Q30.txt -output /home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis_uid/output/v300013477_L03-1-24HD753_single.uid.sorted.freq.paired.Q30.txt_freq -cutoff 0.01
Command:	
	-input str   ides-bam2freq.pl or ides-polishbg.pl output file.
	-output str  output file.
	-cutoff	f 	print min cutoff.
Author:   Zhiliang Fu fuzl\@geneis.cn, QQ:594380908
Version:  v1.0
Update:   2019/9/1

\n!
    )
}




