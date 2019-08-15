# DNA
DNA_mutation_analysis

Fastqc_trim_bwa_gatk4_mutect2_HC_freebayes_annovar_pipeline  
git init  
git add README.md  
git commit -m "DAN commit"  
git status  
git push -u origin master   
git branch my-branch  
git checkout my-granch   


#################   

perl Fastqc_bwa_gatk4_Mutect2_HC_freebayes_annovar_v3.2.pl  
  
Usage:  
        eg:  
perl Fastqc_bwa_gatk4_Mutect2_HC_freebayes_annovar_v3.2.pl -fq1 /home/fuzl/project/demo/demo_1.fq -fq2 /home/fuzl/project/demo/demo_2.fq -outdir /home/fuzl/project/demo/GATK4demo -gatktype 1 -vf 10  
Function: Template for Perl FASTQC BWA GATK pipeline .  
Command:  
        -fq1 str        fq1 fastq format *_1.fq or *_1.fq.gz  
        -fq2 str        fq2  
        -outdir str  outdir  
        -detail_cfg str detail cfg  [config/detail_cfg.txt]  
        -gatktype INT   gatk type, 1=HaplotypeCaller, 2=Mutect2 ,4=freebayes [4]  
        -bed str        bed file  
        -fusion      
        -onlyprintcmd  
        -help  
Author:   Zhiliang Fu fuzl@geneis.cn, QQ:594380908  
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
 
Version: 3.2  
Update:   2019/8/10  
freebayes bed 用bedtools 拆分，设置最大拆成10份  
计算mapping bases 用去重后的数据  

