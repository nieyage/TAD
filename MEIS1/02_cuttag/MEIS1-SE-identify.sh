
/public/home/nieyg/software/ROSE/annotation/hg19_refseq.ucsc

PATHTO=$HOME/software/ROSE-master
PYTHONPATH=$PATHTO/lib
export PYTHONPATH
export PATH=$PATH:$PATHTO/bin



/public/home/nieyg/software/ROSE/ROSE_main.py -g HG19 -i ${id}.narrowPeak.bed \
-r ${id}_H3K27AC.deduplicate.chr.bam \
-c ${id}_H3K27AC_INPUT.deduplicate.chr.bam \
-o ./${id}/ \
-s 12500 -t 2500 2>${id}.log




[optional: -s STITCHING_DISTANCE -t TSS_EXCLUSION_ZONE_SIZE -c CONTROL_BAM]
参数解释
-g refseq参考基因组
-i 输入gff文件
-r 排序后的bam文件，同时需为bam添加index
-o 输出文件目录
可选参数
-s STITCHING_DISTANCE，合并两个region的最大距离，默认值为12.5kb
-t TSS_EXCLUSION_ZONE_SIZE，排除TSS区域大小，排除与TSS前后某距离内的区域，以排除启动子偏差（默认值：0;推荐值：2500）。如果设置该值为0，将不会查找基因。
-c CONTROL_BAM，control样本的bam文件

-i后面的gff文件可以直接用MACS2的narrowPeak结果，但是文件名称要以bed结尾。

conda activate python27
cd /public/home/nieyg/software/ROSE
python ROSE_main.py -g HG19 -i ~/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/6_SE/MKac-1_peaks.narrowPeak.bed \
-r ~/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/2_bam/MKac-1_sorted_rmDup_mapped_rmbl.bam \
-c ~/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/2_bam/D5-IG_sorted_rmDup_mapped_rmbl.bam \
-o ~/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/6_SE/MKac-1/ \
-s 12500 -t 2500

python ROSE_main.py -g HG19 -i ~/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/6_SE/MKac-2_peaks.narrowPeak.bed \
-r ~/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/2_bam/MKac-2_sorted_rmDup_mapped_rmbl.bam \
-c ~/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/2_bam/D5-IG_sorted_rmDup_mapped_rmbl.bam \
-o ~/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/6_SE/MKac-2/ \
-s 12500 -t 2500

python ROSE_main.py -g HG19 -i ~/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/6_SE/WTac-1_peaks.narrowPeak.bed \
-r ~/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/2_bam/WTac-1_sorted_rmDup_mapped_rmbl.bam \
-c ~/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/2_bam/D5-IG_sorted_rmDup_mapped_rmbl.bam \
-o ~/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/6_SE/WTac-1/ \
-s 12500 -t 2500

python ROSE_main.py -g HG19 -i ~/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/6_SE/WTac-2_peaks.narrowPeak.bed \
-r ~/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/2_bam/WTac-2_sorted_rmDup_mapped_rmbl.bam \
-c ~/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/2_bam/D5-IG_sorted_rmDup_mapped_rmbl.bam \
-o ~/project/TAD/MEIS1/cuttag/20221015-MEIS1-H3K27ac-cuttag/6_SE/WTac-2/ \
-s 12500 -t 2500



python ROSE_geneMapper.py [options] -g [GENOME] -i [INPUT_ENHANCER_FILE]

-i INPUT, --i=INPUT   Enter a ROSE ranked enhancer or super-enhancer file
  -g GENOME, --genome=GENOME
                        Enter the genome build (MM9,MM8,HG18,HG19)
  -l GENELIST, --list=GENELIST
                        Enter a gene list to filter through
  -o OUT, --out=OUT     Enter an output folder. Default will be same folder as
                        input file
  -w WINDOW, --window=WINDOW
                        Enter a search distance for genes. Default is 50,000bp
  -f, --format          If flagged, maintains original formatting of input
                        table(chipseq) 

ROSE_bamToGFF.py: calculates density of .bam sequencing reads in .gff regions
ROSE_callSuper.R: ranks regions by their densities, creates a cutoff to separate super-enhancers from typical enhancers




