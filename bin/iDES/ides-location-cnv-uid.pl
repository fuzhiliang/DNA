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


my ($input,$hotspot_vcf,$output,$cutoff);
my %opts;
GetOptions(
    'verbose' => \$verbose,
	'input=s' => \$input,  #uid bam 
	'hotspot_vcf=s' =>\$hotspot_vcf,
	'output=s' => \$output, 
#	'cutoff=f' => \$cutoff,
    "help|h" =>\&usage,
) or die $!;
unless(defined $input && defined $hotspot_vcf && defined $output ){&usage();exit 0;}
#bam 
#v300013477L4C001R0030646438     99      chr3    178935992       40      100M    =       178936049       157     TTACAGAGTAACAGACTAGCTAGAGACAATGAATTAAGGGAAAATGACAAAGAACAGCTCAAAGCAATTTCTACACGAGATCCTCTCTCTGAAATCACTG    GFGFFFEFFFFFGFFFGFFFFGGFFGFGGFFFFFGFFFF7@FFGGFG>GFFFBGFFFFGDGGFGFGGFFFEFEGBFFFGFGGFFFFGGFG>EFGFGGFGG    XA:Z:chr22,+17052912,100M,0;    MC:Z:100M       MD:Z:100        RG:Z:v300013477_L04-1-24HD753   NM:i:0  AS:i:100        XS:i:100        ID:Z:GTGTGGATG

#chr3    178936090       178936091       G_A_PIK3CA;PIK3CA:NM_006218:exon10:c.G1633A:p.E545K
#
my %mutation;
open V ,"$hotspot_vcf" or die $!;
while (<V>){
	chomp;
	my ($chr,$loc,$ref)=(split "\t",$_)[0,2,3];
	$mutation{$chr}{$loc}=$ref;

}

open I,"samtools view $input -L /home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/mutation/snv.contr.bed | " or die $!;
open O,">$output" or die $!;
while (<I>) {
	chomp;
	my @tmp=split "\t",$_;
	my $length=length($tmp[9]);
	foreach my $chr (keys %mutation){
		next unless (exists $mutation{$chr});
		foreach my $vcf_loc (keys %{$mutation{$chr}}){
			next if ($vcf_loc < $tmp[3]);
			next if ($vcf_loc-$tmp[3]>$length);
			my $reads_location=$vcf_loc-$tmp[3];
			my $base=substr($tmp[9],$reads_location,1);
			print O "$_\t$mutation{$chr}{$vcf_loc}\t$vcf_loc\t$base\n";
		}
	}
}

close V;
close I;
close O;
`awk \'\{print \$3\"\\t\"\$(NF-3)\"\\t\"\$(NF-2)\"\\t\"\$(NF-1)\"\\t\"\$NF\}\' $output |sort -k 2 >$output.uid`;
open U ,"$output.uid" or die $!;
open FO ,">$output.uid.make" or die $!;
my %UID;
my %UID_count;
while (<U>){
	chomp;
	my @tmp=split "\t",$_;
	my $u="$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]";
	$UID{$u}{$tmp[4]}++;
}
foreach my $key(sort {$a cmp $b}keys %UID){
	my @num=keys %{$UID{$key}};
	foreach my $base(keys %{$UID{$key}}){
		if ($#num>0){
			print FO "$key\t$base\t$UID{$key}{$base}\tfilter\n";
		}else{
			print FO "$key\t$base\t$UID{$key}{$base}\tpass\n";
			my ($chr,$ref,$loc)=(split "\t", $key)[0,2,3];
			my $l="$chr\t$ref\t$loc\t$base";
			$UID_count{$l}+=$UID{$key}{$base};
		}
	}
}

close U;
close FO;
open CO ,">$output.uid.make.count" or die $!;
foreach my $key(sort {$a cmp $b }keys %UID_count){
	print CO "$key\t$UID_count{$key}\n";
}

sub usage {
    die(
        qq!
Usage:根据iDES的结果，计算输出突变频率大于cutoff的所有位点。
	eg:
perl $0 -input /home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis_uid/output/v300013477_L03-1-24HD753_single.uid.sorted.freq.paired.Q30.txt -output /home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis_uid/output/v300013477_L03-1-24HD753_single.uid.sorted.freq.paired.Q30.txt_freq -cutoff 0.01
Command:	
	-input str   uid make format of bam.
	-output str  output file.
	-hotspot_vcf	f 	vcf of filter.
Author:   Zhiliang Fu fuzl\@geneis.cn, QQ:594380908
Version:  v1.0
Update:   2019/9/1

\n!
    )
}


