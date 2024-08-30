# /storage/zhangyanxiaoLab/suzhuojie/software/kraken2/kraken2-build --standard --threads 24 --db /storage/zhangyanxiaoLab/suzhuojie/software/kraken2/database/
nohup kraken2 --quick --paired --db /storage/zhangyanxiaoLab/suzhuojie/software/kraken2/database/ \
    --classified-out cseqs#.fq \
    /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/CKJ006_deep_RNA/raw_data/CKJ006_R1.fq.gz \
    /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/CKJ006_deep_RNA/raw_data/CKJ006_R2.fq.gz  \
    --output /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/CKJ006_deep_RNA/kraken2/CKJ006_report.txt \
    --report /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/CKJ006_deep_RNA/kraken2/CKJ006_species_report.txt 2>&1>/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/logs/CKJ006_deep_kraken.log &