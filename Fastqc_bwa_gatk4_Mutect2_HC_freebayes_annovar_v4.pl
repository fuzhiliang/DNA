#!/usr/bin/perl -w
use strict;
use warnings;
use POSIX qw(strftime);
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin $Script);
use Cwd 'abs_path';
my $verbose ="v2.0";

###############################################################################
#Get the parameter and provide the usage.
###############################################################################
my ($fq1,$fq2,$outdir,$detail_cfg,$gatktype,$onlyprintcmd,$bed,$bdir,$step,$step_by_step,$fusion,$data_analyzer,$uid,$single,$cutoff);
my %opts;
GetOptions(
    'verbose' => \$verbose,
    'fq1=s' => \$fq1,
	'fq2=s' => \$fq2,
	'outdir=s' => \$outdir, 
	'detail_cfg=s' => \$detail_cfg,
	'bed=s' =>\$bed,
	'gatktype=i'=>\$gatktype,
	'bd:s'   =>\$bdir,
    'step:s' =>\$step,
    'sbs'  =>\$step_by_step,
	'onlyprintcmd' =>\$onlyprintcmd,
    'PL:s'   =>\$data_analyzer,
    'fusion' =>\$fusion,
    'uid'	=>\$uid,
    'single' =>\$single,
    'cutoff=s'=>\$cutoff,
    "help|h" =>\&usage,
) or die $!;
unless(defined $fq2 && defined $fq1 && defined $outdir ){&usage();exit 0;}

#######配置文件
$detail_cfg||="$Bin/config/ref_trans.detail.cfg";
$detail_cfg=abs_path($detail_cfg);
&log_current_time("$Script start...");
my (%detail_cfg);
&detail_cfg_read($detail_cfg,\%detail_cfg);

die "提供参数single时，必须有参数uid.\n" if (defined $single && !defined $uid);
########设置默认参数
$bed||="/home/fuzl/bed/Illumina_WES.bed";
$data_analyzer||='lip';
my $vf=$detail_cfg{gatk_vf};
my $threads=$detail_cfg{threads};
my $threads_bayes=$detail_cfg{threads_bayes};

########

########配置软件和数据库
my $SOFT="/home/fuzl/soft";
my $fastqc=$detail_cfg{fastqc};
my $qualimap=$detail_cfg{qualimap} ;
my $picard="java -jar $detail_cfg{picard} ";
my $samtools=$detail_cfg{samtools} ;
my $GATK=$detail_cfg{GATK_4};
my $bwa=$detail_cfg{bwa} ;
my $STAR_Fusion=$detail_cfg{STAR_Fusion}; 
my $Trimmomatic=$detail_cfg{Trimmomatic}; 
my $freebayesp=$detail_cfg{freebayesp}; 
my $freebayes=$detail_cfg{freebayes}; 
my $annovar=$detail_cfg{annovar}; 
my $split_multi=$detail_cfg{split_multi}; 

#################
my ($genome,$INDEX,$dbSNP,$phasel_1KG,$Mills_and_1KG,$Fusion_database);

my $hg19_root="$Bin/database_hg19";
$genome=$detail_cfg{genome};
$INDEX=$detail_cfg{INDEX};
$dbSNP=$detail_cfg{dbSNP};
$phasel_1KG=$detail_cfg{phasel_1KG};
$Mills_and_1KG=$detail_cfg{Mills_and_1KG};
$Fusion_database=$detail_cfg{Fusion_database};
my $humandb=$detail_cfg{humandb}; # 新版本注释数据库

########输出目录和备份目录
if (defined $bdir) {
   system "rm -r $bdir" if (-d $bdir);
   &MKDIR ($bdir);
   $bdir=abs_path($bdir);
}

$outdir=~'s/\/$//';
$outdir=abs_path($outdir);
&MKDIR ($outdir);

########
#my %step;
#$step ||= ($step_by_step) ? '1' : join ',',(1..7);
#&steps_process($step,$step_by_step,\%step);


my $start_time=time;
print strftime("Start time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "Input file is: \nfq1: $fq1\nfq2: $fq2\nOutput directory is: $outdir\nbed :$bed\n";
print "Database file is: $genome\n";

###############################################################################
my $cmd="";


#############################################
#QC
$cmd="";
$fq1=abs_path($fq1);
$fq2=abs_path($fq2);
die "Pleace check format of fastq filename, must end of '_1.fq' or '_1.fq.gz':\n$fq1\n" unless ($fq1=~/_1.f/) ;
my $fq_name=basename($fq1);
$fq_name=~s/_1.f(ast)?q(\.gz)?$//;

&MKDIR("$outdir/FASTQC");
$cmd .="$fastqc $fq1 -o $outdir/FASTQC/  && cd $outdir/FASTQC/ &&  unzip -o $outdir/FASTQC/${fq_name}_1_fastqc.zip \n";
$cmd .="$fastqc $fq2 -o $outdir/FASTQC/  && cd $outdir/FASTQC/ &&  unzip -o $outdir/FASTQC/${fq_name}_2_fastqc.zip \n ";

&MKDIR("$outdir/trim/");
$cmd .="java -jar $Trimmomatic PE -threads $threads -phred33 -trimlog $outdir/trim/logfile ";
$Trimmomatic=~m/(.*\/)/;
$cmd .="$fq1 $fq2 $outdir/trim/${fq_name}_1_paired.fq.gz  $outdir/trim/${fq_name}_1_unpaired.fq.gz  $outdir/trim/${fq_name}_2_paired.fq.gz $outdir/trim/${fq_name}_2_unpaired.fq.gz ";
$cmd .=" HEADCROP:4 " if ($single);
$cmd .="ILLUMINACLIP:$1/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \n";

$fq1="$outdir/trim/${fq_name}_1_paired.fq.gz";
$fq2="$outdir/trim/${fq_name}_2_paired.fq.gz";
$fq_name=basename($fq1);
$fq_name=~s/_1_paired.f(ast)?q(\.gz)?$//;
&MKDIR("$outdir/FASTQC_trim/");
$cmd .="$fastqc $fq1 -o $outdir/FASTQC_trim/  && cd $outdir/FASTQC_trim/ && unzip -o $outdir/FASTQC_trim/${fq_name}_1_paired_fastqc.zip & \n";
$cmd .="$fastqc $fq2 -o $outdir/FASTQC_trim/  && cd $outdir/FASTQC_trim/ && unzip -o $outdir/FASTQC_trim/${fq_name}_2_paired_fastqc.zip & \n ";

&runcmd("1 QC",$cmd);


#############################################################
#bwa
my $sample=$fq_name;
&MKDIR("$outdir/GATK/");
&MKDIR("$outdir/bwa");

$cmd="";
#$cmd .="mkdir $outdir/bwa \n" unless (-d "$outdir/bwa");
$cmd .="$bwa mem -M -t $threads -R '\@RG\\tID:$sample\\tSM:$sample\\tLB:$sample\\tPL:Illumina' $INDEX  $fq1 $fq2 >$outdir/bwa/$sample.sam  \n";
&runcmd("2 BWA",$cmd);

################################
$cmd ="";
$cmd .="$GATK --java-options '-Xmx${vf}G -Djava.io.tmpdir=./'  SortSam -SO coordinate  -I $outdir/bwa/$sample.sam  -O $outdir/GATK/$sample.bam  1>$outdir/GATK/log.sort 2>&1\n";
$cmd .="$samtools index $outdir/GATK/$sample.bam \n";
&runcmd("3 GATK_1 SortSam",$cmd);
my $bam;
if ($uid){
	#`awk 'NR==FNR{if (NR%4==1){a[$1]=\$NF}}NR>FNR{if (/^@/){print $0}else{if ("@"$1"/1" in a){print $0"\tID:Z:"a["@"$1"/1"]}}}' <(zcat $fq1 ) <(samtools view -h $outdir/GATK/$sample.bam) |awk  '{if(\$NF == CUST_ID ){$2+=1024;print $0 }else{if(\$NF != CUST_ID){CUST_ID=\$NF;print $0}}}' OFS="\t" |samtools view -Sb >$outdir/GATK/$sample.uid_make.bam && touch "3_GATK_2_UID_MarkDuplicates.sh.finish" ` unless (-f "3_GATK_2_UID_MarkDuplicates.sh.finish");
	$cmd="gunzip -c $fq1>$outdir/GATK/$sample.1.fq \n";
	$cmd .="samtools view -h $outdir/GATK/$sample.bam >$outdir/GATK/$sample.sam \n";
	$cmd .=" awk \'NR==FNR\{if \(NR%4==1\)\{a[\$1]=\$NF\}\}NR>FNR\{if \(\/\^\@\/\)\{print \$0\}else\{if \(\"\@\"\$1\"\/1\" in a\){print \$0\"\\tID:Z:\"a[\"\@\"\$1\"\/1\"]\}\}\}\' $outdir/GATK/$sample.1.fq $outdir/GATK/$sample.sam ";
	$cmd .="|awk  \'\{if\(\$NF == CUST_ID \)\{\$2+=1024;print \$0 \}else\{if\(\$NF != CUST_ID\)\{CUST_ID=\$NF;print \$0\}\}\}\' OFS=\"\\t\" |samtools view -Sb >$outdir/GATK/$sample.uid_make.sorted.bam \n ";
	$bam="$outdir/GATK/$sample.uid_make.sorted.bam ";
	$cmd .="samtools index $bam \n";
	$cmd .="rm $outdir/GATK/$sample.1.fq && rm $outdir/GATK/$sample.sam \n";
	#`$cmd && touch "3_GATK_2_UID_MarkDuplicates.sh.finish" ` unless (-f "3_GATK_2_UID_MarkDuplicates.sh.finish");
	&runcmd("3 GATK_2 UID_MarkDuplicates",$cmd);
	#$bam="$outdir/GATK/$sample.uid_make.sorted.bam ";


}else{
	#################################
	$cmd ="";
	$cmd .="$GATK  --java-options '-Xmx${vf}G -Djava.io.tmpdir=./'   MarkDuplicates  -I $outdir/GATK/$sample.bam -O $outdir/GATK/${sample}_marked.bam -M $outdir/GATK/${sample}.metrics 1>$outdir/GATK/log.mark  2>&1\n";
	$bam="$outdir/GATK/${sample}_marked.bam";
	$cmd .="$samtools index  $bam \n";
	&runcmd("3 GATK_2 MarkDuplicates",$cmd);

	#################################
	$cmd ="";
	$cmd .="$GATK  --java-options '-Xmx${vf}G -Djava.io.tmpdir=./'   FixMateInformation -I $bam -O $outdir/GATK/${sample}_marked_fixed.bam -SO coordinate   1>$outdir/GATK/log.fix 2>&1 \n";

	$bam="$outdir/GATK/${sample}_marked_fixed.bam";
	$cmd .="$samtools index  $bam \n";

	&runcmd("3 GATK_3 FixMateInformation",$cmd);
	##################################recal
	$cmd ="";
	$cmd .="$GATK  --java-options '-Xmx${vf}G -Djava.io.tmpdir=./'   BaseRecalibrator -R $genome  -I $bam -O $outdir/GATK/${sample}_recal_data.table --use-original-qualities ";
	$cmd .="--known-sites  $dbSNP " ;
	$cmd .="--known-sites  $phasel_1KG " ;
	$cmd .="--known-sites  $Mills_and_1KG " ;
	$cmd .=" 1>$outdir/GATK/log.recal1 2>&1\n";

	$cmd .="$GATK  --java-options '-Xmx${vf}G -Djava.io.tmpdir=./'  GatherBQSRReports -I  $outdir/GATK/${sample}_recal_data.table  -O  $outdir/GATK/${sample}_gether_recal_data.table 1>$outdir/GATK/log.recal2 2>&1\n"; #支持多个lane 一起矫正，-I (seq -I *table)

	$cmd .="$GATK  --java-options '-Xmx${vf}G -Djava.io.tmpdir=./' ApplyBQSR -R $genome  -I $bam -O $outdir/GATK/${sample}_recal.sorted.bam ";
	$cmd .="-bqsr  $outdir/GATK/${sample}_gether_recal_data.table  --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 --add-output-sam-program-record  --create-output-bam-md5  --use-original-qualities ";
	$cmd .=" 1>$outdir/GATK/log.recal3 2>&1\n";
	$bam=" $outdir/GATK/${sample}_recal.sorted.bam ";
	$cmd .="$samtools index  $bam \n";
	&runcmd("3 GATK_4 Recal",$cmd);
} 
############################
$cmd="";
$cmd .="$samtools flagstat $bam > $outdir/GATK/${sample}.alignment.flagstat & \n ";
$cmd .="$samtools stats -d $bam > $outdir/GATK/${sample}.alignment.stat  & \n";
$cmd .="awk \'\{print \$1\"\\t\"\$2\"\\t\"\$3\"\\ttmp\\t0\\t\+\"\}\'  $bed > $outdir/GATK/bed.m \n";
$bam=~s/\s+$//;
$cmd .="$qualimap  bamqc  -bam $bam -gff $outdir/GATK/bed.m   -outdir ${bam}_QC -outfile $sample -outformat PDF:HTML --java-mem-size=${vf}G -sd \n";

if ($uid) {
	&runcmd("4 Qualimap UID",$cmd) if ($bed);
}else{
	&runcmd("4 Qualimap",$cmd) if ($bed);
}
#################################
#Summary
open S ,">$outdir/Summary.txt" or die $! ;
#my $Summary="Total_Sequences_fq1\tTotal_Sequences_fq2\tTotal_Sequences_fq1_trim\tTotal_Sequences_fq2_trim\tReads_mapped\tOn_target\tcoverage50X\tMean_coverageData\tGC%\n";
my $Summary="RawData\tClearData\tReads_mapped\tBases_mapped\tOn_target_reads\tOn_target_bases\tCoverage50X\tMean_coverageData\tGC%\tDuplicate\tTarget_length\n";
my $Total_Sequences_fq1=`grep "Total Sequences"  $outdir/FASTQC/${sample}_1_fastqc/fastqc_data.txt |cut -f 2 `;chomp ($Total_Sequences_fq1);
my $Total_Sequences_fq2=`grep "Total Sequences"  $outdir/FASTQC/${sample}_2_fastqc/fastqc_data.txt |cut -f 2 `;chomp ($Total_Sequences_fq2);
my $Total_Sequences_fq1_trim=`grep "Total Sequences"  $outdir/FASTQC_trim/${sample}_1_paired_fastqc/fastqc_data.txt |cut -f 2 `;chomp ($Total_Sequences_fq1_trim);
my $Total_Sequences_fq2_trim=`grep "Total Sequences"  $outdir/FASTQC_trim/${sample}_2_paired_fastqc/fastqc_data.txt |cut -f 2 `;chomp ($Total_Sequences_fq2_trim);

my $Bases_mapped=`grep "bases mapped:" $outdir/GATK/${sample}.alignment.stat |cut -f 3`;chomp ($Bases_mapped);
my $Reads_mapped=`grep "number of mapped reads" ${bam}_QC/genome_results.txt |cut -d "=" -f 2 `;chomp ($Reads_mapped);
my $Bases_target_mapped=`grep "number of mapped bases" ${bam}_QC/genome_results.txt |cut -d "=" -f 2 `;chomp ($Bases_target_mapped);
$Bases_target_mapped=~s/\sbp//;$Bases_target_mapped=~s/,//g;
my ${coverage50X}=`grep "coverageData >= 50X" ${bam}_QC/genome_results.txt |awk '{print \$4}'`;chomp (${coverage50X});
my $mean_coverageData=`grep "mean coverageData" ${bam}_QC/genome_results.txt |awk '{print \$4}'|sed 's/X//'`;chomp ($mean_coverageData);
#my $On_target= `sed -n 173p  ${bam}_QC/qualimapReport.html|cut -d ">" -f 2|cut -d "<" -f 1|sed 's/\\s\\+\\/\\s\\+/ \\(/'|sed 's/\$/\\)/'`;chomp ($On_target);##qualimap 添加-sd参数，考虑dup
my $On_target= `sed -n 176p  ${bam}_QC/qualimapReport.html|cut -d ">" -f 2|cut -d "<" -f 1|sed 's/\\s\\+\\/\\s\\+/ \\(/'|sed 's/\$/\\)/'`;chomp ($On_target);
my $target=`grep -A 1  "Regions size/percentage of reference"  ${bam}_QC/qualimapReport.html|tail -n 1`; 
$target=~m/(\d*((,\d*)+)?)\s\/\s((\d+.)?\d+\%)/;
#print "$1\t$3";die;
my $target_length="$1($4)";
my $length=$1;$length=~s/,//g;
#my $target_bases=$length*$mean_coverageData;
my $target_bases_freq=$Bases_target_mapped/$Bases_mapped*100;$target_bases_freq=sprintf("%0.2f",$target_bases_freq);
#print "$length\t$mean_coverageData\t$Bases_mapped";
my $GC=`grep "GC percentage" ${bam}_QC/genome_results.txt |cut -d "=" -f 2 `;chomp ($GC);
my $Duplicate=`grep "number of duplicated reads" ${bam}_QC/genome_results.txt |cut -d "=" -f 2 `;chomp ($Duplicate);
my $RawData=$Total_Sequences_fq1+$Total_Sequences_fq2;
my $ClearData=$Total_Sequences_fq1_trim+$Total_Sequences_fq2_trim;
my $Reads_mapped_reads=(split " ", $Reads_mapped)[0];
$Duplicate=~s/,//g;
$Reads_mapped_reads=~s/,//g;
my $Duplicat_freq=100*$Duplicate/$Reads_mapped_reads;
$Duplicat_freq=sprintf("%4.2f" , $Duplicat_freq);
$Summary .="$RawData\t$ClearData\t$Reads_mapped\t$Bases_mapped\t$On_target\t$Bases_target_mapped($target_bases_freq%)\t$coverage50X\t$mean_coverageData\t$GC\t$Duplicate\($Duplicat_freq%\)\t$target_length\n";
$Summary=~s/\s\(/\(/g;
$Summary=~s/,//g;
print S "$Summary";
print "\nSummary done.\n";
######################################################

$cmd ="";
my $vcf="";
$outdir="$outdir/uid" if ($uid);
&MKDIR("$outdir/vcf");

$gatktype||=2;
if ($gatktype == "2"){
	$cmd .="$GATK  --java-options '-Xmx${vf}G -Djava.io.tmpdir=./'   Mutect2  -R $genome -I $bam ";
	#$cmd .=" --dbsnp $dbSNP "; 
	$cmd .=" -tumor $sample "; 
	$cmd .=" --intervals $bed ";
	#--heterozygosity
	$cmd .=" -O  $outdir/vcf/${sample}_mutect2.raw.vcf 1>$outdir/vcf/log.MT 2>&1\n";
	$vcf="$outdir/vcf/${sample}_mutect2.raw.vcf ";
}elsif($gatktype == "1"){
	$cmd .="$GATK  --java-options '-Xmx${vf}G -Djava.io.tmpdir=./'   HaplotypeCaller  -R $genome -I $bam ";
	$cmd .=" --dbsnp $dbSNP "; 
	$cmd .=" -O  $outdir/vcf/${sample}.HC.raw.vcf 1>$outdir/vcf/log.HC  2>&1 \n";
	$vcf="$outdir/vcf/${sample}.HC.raw.vcf ";
#	&filter($vcf,$outdir) ;
#	$vcf="$outdir/vcf/${sample}.filter.PASS.vcf ";
}elsif($gatktype == "4"){
	#$freebayes="$freebayesp  $threads_bayes	$freebayes ";
	system "rm $outdir/temp_bed/*bed" if (-d "$outdir/temp_bed");
	&MKDIR("$outdir/temp_bed");
	$threads_bayes=($threads_bayes>10?10:$threads_bayes); #默认最多拆10份
	my $bedtools=$detail_cfg{bedtools};
	#$cmd.="$bedtools  split -i $bed -n $threads_bayes -p $outdir/temp_bed/a ";
	`$bedtools  split -i $bed -n $threads_bayes -p $outdir/temp_bed/a `;
	$cmd="";
	my @bed=glob("$outdir/temp_bed/a*.bed");
	my @finish=glob("$outdir/temp_bed/*finish");
	print "finish number is: $#finish+1\nbed number is: $#bed+1\n";
	system "rm $outdir/shell/5.0_freebayes_split.sh.finish " if ($#finish!=$#bed);
	foreach my $bed(@bed){
		$cmd .="$bedtools sort -i $bed >$bed.m && mv $bed.m $bed \n";
		$cmd .="$freebayes -j -m 10 -q 30 -F 0.001 -C 1 -t $bed -i --no-indels -f $genome $bam > $bed.snp.freebayesp_q30.vcf && ";
		$cmd .="$freebayes -j -m 10 -q 20 -F 0.001 -C 1 -t $bed -I  -f $genome $bam > $bed.indel.freebayesp_q20.vcf && touch \"$bed.finish\" || touch \"$bed.error\" &\n";
	}
	&runsh("5.0 freebayes split",$cmd,"$outdir/temp_bed/");
	my @error;
	until($#finish==$#bed) {
		@finish=glob("$outdir/temp_bed/*finish");
		@error=glob("$outdir/temp_bed/*.error");
		exit(1) if ($#error>=0);
		sleep 2;
	}
   	system "touch $outdir/shell/5.0_freebayes_split.sh.finish ";

	print "split bed done.\n";
	$cmd ="cat $outdir/temp_bed/*.snp.freebayesp_q30.vcf|grep -v \"#\" |sort >$outdir/vcf/$sample.snp.freebayesp_q30.vcf \n";
	$cmd .="cat $outdir/temp_bed/*.indel.freebayesp_q20.vcf|grep -v \"#\" |sort >$outdir/vcf/$sample.indel.freebayesp_q20.vcf \n";
	$cmd .="cat $outdir/vcf/$sample.snp.freebayesp_q30.vcf $outdir/vcf/$sample.indel.freebayesp_q20.vcf  >$outdir/vcf/$sample.vcf \n";
	$cmd .= "perl $annovar/table_annovar.pl $outdir/vcf/$sample.vcf  $humandb -buildver hg19 "; 
	$cmd .= "-out  $outdir/vcf/${sample}  -remove -protocol refGene,cytoBand,dbnsfp30a,cosmic83,snp144,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,clinvar_20170905,esp6500siv2_all ";
	$cmd .= "-operation g,r,f,f,f,f,f,f,f,f,f,f,f --buildver hg19 --nastring . --vcfinput --otherinfo --dot2underline >$outdir/vcf/annovar.log \n";
	$cmd .= "perl $split_multi $outdir/vcf/${sample}.hg19_multianno.txt > $outdir/vcf/${sample}.hg19_multianno.txt.freq \n";
}else{
	die "Please select call mutation tools. \n"; 
}
if ($uid){
	&runcmd("5 Genotype Calling UID",$cmd);
}else{
	&runcmd("5 Genotype Calling",$cmd);
}

if($gatktype == "1"){
	&filter($vcf,$outdir);
	$vcf="$outdir/vcf/${sample}.filter.PASS.vcf ";
}
################################
#iDES 降噪
###############################
my $iDES=$detail_cfg{iDES};
&MKDIR("$outdir/iDES_analysis/input");
&MKDIR("$outdir/iDES_analysis/output");
my $iDES_output="$outdir/iDES_analysis";
$cmd="";
$cutoff||=0.001;
my $bamname=basename($bam);
$bamname=~s/.bam$//;
$cmd.="perl $iDES/ides-bam2freq.pl -o $iDES_output/input -t $threads $bam  $genome $bed \n";
$cmd.="perl $iDES/ides-freq2vcf.pl -input $iDES_output/input/$bamname.freq.paired.Q30.txt";
$cmd.=" -output $iDES_output/input/$bamname.freq.paired.Q30.cutoff$cutoff\_freq -cutoff $cutoff  \n";
#perl /home/fuzl/soft/iDES/ides-freq2vcf.pl -input $i -output $i.cutoff0.01_freq -cutoff 0.01 &
print "iDES output:\n$iDES_output/input/$bamname.freq.paired.Q30.txt\n";


$cmd.="perl $iDES/ides-polishbg.pl -o $iDES_output/output $iDES_output/input/$bamname.freq.paired.Q30.txt ";
$cmd.=" $iDES/519MGI_baseline/output/519MGI-db_ides-bgdb.txt \n";
print "iDES output rmbg:\n$iDES_output/output/$bamname.freq.paired.Q30.rmbg.txt\n";
$cmd.="perl $iDES/ides-freq2vcf.pl -input $iDES_output/output/$bamname.freq.paired.Q30.rmbg.txt";
$cmd.=" -output $iDES_output/output/$bamname.freq.paired.Q30.rmbg.cutoff$cutoff\_freq -cutoff $cutoff  \n";
#awk 'NR==FNR{a[$1"\t"$2"\t"$15]=$0}NR>FNR{if (FNR==1){print "Chr\ts\te\tref\talt\tgene\tp\t"a["CHR\tPOS\tALT"]};{if ($1"\t"$2"\t"$5 in a){print $0"\t"a[$1"\t"$2"\t"$5]}}}' \
#$cmd.="awk \'NR==FNR\{a[\$1\"\\t\"\$2\"\\t\"\$15]=\$0\}NR>FNR\{if (FNR==1)\{print \"Chr\\ts\\te\\tref\\talt\\tgene\\tp\\t\"a[\"CHR\\tPOS\\tALT\"]\};\{if (\$1\"\\t\"\$2\"\\t\"\$5 in a)\{print \$0\"\\t\"a[\$1\"\\t\"\$2\"\\t\"\$5]\}\}\}\' ";
$cmd.="awk \'NR==FNR\{a[\$1\"\\t\"\$2\"\\t\"\$15]=\$0\}NR>FNR\{if \(\$1\"\\t\"\$2\"\\t\"\$15 in a\)\{print a[\$1\"\\t\"\$2\"\\t\"\$15]\"\\t\"\$0\}\}\' ";
$cmd.=" $iDES_output/input/$bamname.freq.paired.Q30.cutoff$cutoff\_freq   $iDES_output/output/$bamname.freq.paired.Q30.rmbg.cutoff$cutoff\_freq " ;
$cmd.=" >$iDES_output/$bamname.Q30_vs_rmbg.txt \n";

&runcmd("6.0 iDES",$cmd) ;

####
#关注的区域，用uid 


#&MKDIR("$outdir/ANNOVAR");
$cmd ="";
$vcf=~s/\s+$//;
$cmd .= "perl  $annovar/convert2annovar.pl -format  vcf4  $vcf  >  $vcf.avinput \n ";
$cmd .= "perl $annovar/table_annovar.pl  $vcf.avinput  $humandb -buildver hg19 "; 
$cmd .= "-out  $vcf.annovar  -remove -protocol refGene,cytoBand,dbnsfp30a,cosmic83,snp144,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,clinvar_20170905,esp6500siv2_all ";
$cmd .= "-operation g,r,f,f,f,f,f,f,f,f,f,f,f --buildver hg19 --nastring . --vcfinput --otherinfo --dot2underline \n";
&runcmd("7 Annovar",$cmd) unless ($gatktype == "4");

$outdir=~s/\/uid// unless ($uid);
#my $size=`awk 'BEGIN{sum=0}{sum+=\$3-\$2+1}END{print sum}' $bed`;
#&MKDIR("$outdir/star_fusion");
&MKDIR("$outdir/fusion_factera");

my $SE=$detail_cfg{SE};
#my $fusion_bed ||="$SE/database/fusion_38cf.bed ";
my $fusion_bed=$detail_cfg{fusion_bed};
$cmd="perl /home/liq/hym_source/FACTERA-FUSION/factera.pl -o $outdir/fusion_factera -F $bam /home/liq/hym_source/FACTERA-FUSION/exons.bed  /home/liq/hym_source/FACTERA-FUSION/hg19.2bit  $bed & \n";
$cmd.="perl $SE/SEGF.pl -fq1 $fq1 -fq2 $fq2 -bed $fusion_bed -odir $outdir/fusion_SE -trim_len 10 -remain_len 35 -process $threads  1>/dev/null 2>&1 \n";
$cmd.="rm -r $outdir/fusion_SE/tmp $outdir/fusion_SE/Deal  \n";
#mkdir ${i%%GATK*}fusion_factera
#perl /home/liq/hym_source/FACTERA-FUSION/factera.pl -o ${i%%GATK*}fusion_factera -F $i /home/liq/hym_source/FACTERA-FUSION/exons.bed  /home/liq/hym_source/FACTERA-FUSION/hg19.2bit  /home/fuzl/project/boke519_test9sample/508p.sort.merge.bed &
&runcmd("8.0 Fusion_cmd",$cmd)if ($fusion) ;

#$cmd="";
#$cmd .="$STAR_Fusion --genome_lib_dir $Fusion_database --left_fq $fq1 --right_fq $fq2 --output_dir $outdir/star_fusion --CPU $threads threads\n";
#/home/fuzl/soft/Fusion/STAR-Fusion-v1.5.0/STAR-Fusion --genome_lib_dir /home/fuzl/soft/Fusion/GRCh37_v19_CTAT_lib_Feb092018/ctat_genome_lib_build_dir --left_fq /home/fuzl/project/demo/lib-HD-4_1.fq.gz --right_fq /home/fuzl/project/demo/lib-HD-4_2.fq.gz --output_dir /home/fuzl/project/demo/new_HD-4_star_fusion_outdir --CPU 36
#&runcmd("8.1 STAR Fusion",$cmd)if ($fusion);

`rm $outdir/bwa/$sample.sam ` if (-f "$outdir/bwa/$sample.sam");  
###############################################################################
my $duration_time=time-$start_time;
print strftime("End time is %Y-%m-%d %H:%M:%S\n", localtime(time));
print "This compute totally consumed $duration_time s\.\n";
#system "perl -e 'print \"Dear $data_analyzer:\\n\\t$project_id is done. The log is shown below:\\n\"; print \"\\nResult info:$flag \";' |mail -s \"$project_id $name Result\" $data_analyzer\@geneis.cn ";  

###############################################################################
sub usage {
    die(
        qq!
Usage:
	eg:
perl $0 -fq1 /home/fuzl/project/demo/demo_1.fq -fq2 /home/fuzl/project/demo/demo_2.fq -outdir /home/fuzl/project/demo/GATK4demo -gatktype 1 -vf 10
Function: Template for Perl FASTQC BWA GATK pipeline .
Command:	
	-fq1 str  		fq1 fastq format *_1.fq or *_1.fq.gz 
	-fq2 str  		fq2 
	-outdir	str 	outdir
	-detail_cfg str	detail cfg 	
	-gatktype INT 	gatk type, 1=HaplotypeCaller, 2=Mutect2 ,4=freebayes [2]
	-bed str 		bed file 
	-fusion     
	-uid	    对应MGI平台加UID分析，fq1中read id中必须把uid加到reads id后面，空格分隔
	-single		UID分析，用于切掉5’端4bp的uid
	-cutoff     iDES fliter min frequency          [0.001]
	-onlyprintcmd
	-help
Author:   Zhiliang Fu fuzl\@geneis.cn, QQ:594380908
Version:  v1.0
Update:   2018/8/8
Version:  v2.0
Update:   2019/2/27
Notes:   
添加 mutect2  
添加fusion 
添加检查是否执行sh，有*.sh.finish 时则不执行该shell
GATK 拆分多个shell
线程，内存，融合，参数外置
修改runcmd，添加报错中断机制
添加Summary
修改annovar注释版本
Version:  v3.0
Update:   2019/7/26
配置文件外置
添加freebayes
qualimap 加-sd ,考虑dup
Version:  v3.1
SE 的fusion_bed加都配置文件
SE 线程外置
Summary添加target bases
Version:  v4.0
添加UID分析
添加iDES分析
\n!
    )
}
sub detail_cfg_read {
    &log_current_time("detail config check:");
    my ($cfg_file, $detail_cfg) = @_;

    open (CFG,$cfg_file ) or die "$!: $cfg_file\n";
    while (<CFG>) {
        chomp;
        next if (/^\s+/ or /^#/);
        my ($key, $value) = (split /\s+/)[0,1];
         $detail_cfg->{$key} = $value;
        if ($key eq 'Project_name' or $key eq 'Principal' or $key eq 'Project_id' or $key eq 'Project_key') {
            $detail_cfg->{$key} = $value;
        }
        if ($key eq 'genome' or $key eq 'annovar') {
            die "$key: $value is not exist!\n" unless (-e $value);
            $detail_cfg->{$key} = $value;
        }
        print "$key: $value\n" if (exists $detail_cfg->{$key});
    }
    close CFG;
#    die "Must choose Queue_type1 and Queue_type2 !\n" unless (exists $detail_cfg->{Queue_type1} or exists $detail_cfg->{Queue_type2});
    &log_current_time("detail config check done.");
}

sub runcmd { # 
	my $name=shift @_;
	my $cmd=shift @_;
	&MKDIR("$outdir/shell/");
	my $start_time=time;
#	print strftime("Start $name analysis time is %Y-%m-%d %H:%M:%S\n", localtime(time));
	&log_current_time("Start $name analysis...\n");
#	print "Start $name analysis ... \n";
	my $n=$name;
	$n=~s/\s+/_/g ;
	my $log_file="$outdir/shell/$n.log";
	open CMD ,">$outdir/shell/$n.sh" or die $!;
	my $project_id=$detail_cfg{Project_id};
	print CMD "$cmd";
	close CMD;
	unless (defined $onlyprintcmd || -f "$outdir/shell/$n.sh.finish"){
		my $flag = system("sh $outdir/shell/$n.sh > $log_file") ;
	    if ($flag != 0 ){
	        &log_current_time("Error: command failed: $cmd");
	       	system "touch $outdir/shell/$n.sh.error";
#	       	system "perl -e 'print \"Dear $data_analyzer:\\n\\t$project_id of step $name is error. The log of error is shown below:\\n\"; print \"\\nError info:$flag \";' |mail -s \"$project_id $name Result\" $data_analyzer\@geneis.cn ";  
	        exit(1);
	    } else {
	        my $escaped_time = (time()-$start_time)."s";
	        &log_current_time("$name done, escaped time: $escaped_time.\n");
	       	system "rm $outdir/shell/$n.sh.error " if (-f "$outdir/shell/$n.sh.error");
	       	system "touch $outdir/shell/$n.sh.finish ";
	    }
	}else{
		&log_current_time("$name result has been finish, no duplicate run.");
	}
}
sub runsh { # 
	my $name=shift @_;
	my $cmd=shift @_;
	my $temp=shift @_;
	&MKDIR("$outdir/shell/");
	my $start_time=time;
	&log_current_time("Start $name analysis...\n");
	my $n=$name;
	$n=~s/\s+/_/g ;
	my $log_file="$outdir/shell/$n.log";
	open CMD ,">$outdir/shell/$n.sh" or die $!;
	my $project_id=$detail_cfg{Project_id};
	print CMD "$cmd";
	close CMD;
	unless (defined $onlyprintcmd || -f "$outdir/shell/$n.sh.finish"){
		my $flag = system("sh $outdir/shell/$n.sh > $log_file  & ") ;
	    if ($flag != 0 ){
	        &log_current_time("Error: command failed: $cmd");
	       	system "touch $outdir/shell/$n.sh.error";
	        exit(1);
	    } else {
	        my $escaped_time = (time()-$start_time)."s";
	        &log_current_time("$name done, escaped time: $escaped_time.\n");
	       	system "rm $outdir/shell/$n.sh.error " if (-f "$outdir/shell/$n.sh.error");
	       	system "touch $outdir/shell/$n.sh.finish ";
	    }
	}else{
		&log_current_time("$name result has been finish, no duplicate run.");
	}
}

sub MKDIR{
	my $dir=shift @_;
	system "mkdir  -p $dir " unless (-d "$dir");
}

sub log_current_time {
     # get parameter
     my ($info) = @_;

     # get current time with string
     my $curr_time = &date_time_format(localtime(time()));

     # print info with time
     print "[$curr_time] $info\n";
}

#############################################################################################################
sub date_time_format {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
sub filter{
	my $vcf=shift @_;
	my $outdit=shift @_;
	my $cmd="";
	$cmd .="$GATK --java-options '-Xmx20G -Djava.io.tmpdir=./' SelectVariants -R $genome -V $vcf -select-type SNP  -O $outdir/vcf/${sample}.raw.snp.vcf  \n";
	$cmd .="$GATK --java-options '-Xmx20G -Djava.io.tmpdir=./'  VariantFiltration -R $genome -V $outdir/vcf/${sample}.raw.snp.vcf   --filter-expression  \"QD < 2.0 || FS > 60.0 ||SOR > 3.0|| MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" --filter-name \"my_snp_filter\"  ";
	$cmd .=" -O $outdir/vcf/${sample}.filter.snp.vcf \n";
	#$cmd .="grep -w -v my_snp_filter $outdir/vcf/${sample}.filter.snp.vcf > $outdir/vcf/${sample}.filter.PASS.snp.vcf \n";
	$cmd .="$GATK --java-options '-Xmx20G -Djava.io.tmpdir=./' SelectVariants -R $genome -V $vcf -select-type INDEL  -O $outdir/vcf/${sample}.raw.InDel.vcf \n";
	$cmd .="$GATK --java-options '-Xmx20G -Djava.io.tmpdir=./'  VariantFiltration -R $genome -V $outdir/vcf/${sample}.raw.InDel.vcf   --filter-expression \"QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\"  --filter-name \"my_indel_filter\"  ";
	$cmd .="-O $outdir/vcf/${sample}.filter.InDel.vcf \n";
	$cmd .="$GATK MergeVcfs -I $outdir/vcf/${sample}.filter.InDel.vcf -I $outdir/vcf/${sample}.filter.snp.vcf -O $outdir/vcf/${sample}.filter.vcf \n";
	$cmd .="grep -w -v \"my_snp_filter\" $outdir/vcf/${sample}.filter.vcf|grep -w -v \"my_indel_filter\" > $outdir/vcf/${sample}.filter.PASS.vcf \n";
	#$cmd .="cat $outdir/vcf/${sample}.filter.PASS.snp.vcf $outdir/vcf/${sample}.filter.PASS.InDel.vcf >$outdir/vcf/${sample}.filter.PASS.vcf\n";
	&runcmd("5.1 filter haplotye vcf",$cmd);
	return ($cmd);
}

sub split_bed{
	my $bed=shift @_;
	my $number=shift @_;

}
