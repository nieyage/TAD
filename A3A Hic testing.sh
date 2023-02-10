# A3A Hic testing 

# random select some reads for testing 
seqtk sample -s100 /md01/nieyg/project/A3A/01.RawData/A76/A76_clean/A76/A76clean_R1.fastq.gz 10000 > /md01/nieyg/project/A3A/01.RawData/A76/test_clean/test/test_R1.fq
seqtk sample -s100 /md01/nieyg/project/A3A/01.RawData/A76/A76_clean/A76/A76clean_R2.fastq.gz 10000 > /md01/nieyg/project/A3A/01.RawData/A76/test_clean/test/test_R2.fq

# Hic Pro testing 
conda activate /data/R02/nieyg/software/hic-pro
time /data/R02/nieyg/software/HiC-Pro_3.1.0/bin/HiC-Pro -c config_test_latest.txt -i /md01/nieyg/project/A3A/01.RawData/A76/test_clean/ -o hicpro_latest_test
