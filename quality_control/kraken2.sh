# /storage/zhangyanxiaoLab/suzhuojie/software/kraken2/kraken2-build --standard --threads 24 --db /storage/zhangyanxiaoLab/suzhuojie/software/kraken2/database/
nohup kraken2 --quick --paired --db /storage/zhangyanxiaoLab/suzhuojie/software/kraken2/database/ \
    --classified-out cseqs#.fq \
    /storage/zhangyanxiaoLab/suzhuojie/projects/axon_degradation/data/raw_data/20250512_wanrui_RNA/fastq/Soma_1_S1_L004_R1_001.fastq.gz \
    /storage/zhangyanxiaoLab/suzhuojie/projects/axon_degradation/data/raw_data/20250512_wanrui_RNA/fastq/Soma_1_S1_L004_R2_001.fastq.gz \
    --output /storage/zhangyanxiaoLab/suzhuojie/projects/axon_degradation/data/raw_data/20250512_wanrui_RNA/Soma_1_report.txt \
    --report /storage/zhangyanxiaoLab/suzhuojie/projects/axon_degradation/data/raw_data/20250512_wanrui_RNA/Soma_1_species_report.txt