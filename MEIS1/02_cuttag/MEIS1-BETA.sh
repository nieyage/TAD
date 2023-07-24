conda create -y -n beta_chip python=2.7.15
conda install -y -c hcc beta 
conda install -y libiconv

conda activate beta_chip
BETA plus \
    -p MEIS1-peak.bed \
    -e MS_H1_Day5.txt  \
    -k O \
    --info 1,3,7\
    -g hg19 \
    --gs /public/home/nieyg/reference/genome/hg19/hg19.fa \
    --pn 77870 \
    --gname2 \
    -n MEIS1-Day5 \
    --df 0.05 \
    -o BETA_MEIS1-Day_plus
