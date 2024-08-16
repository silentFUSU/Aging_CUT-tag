/storage/zhangyanxiaoLab/suzhuojie/software/kraken2/kraken2-build --standard --threads 24 --db /storage/zhangyanxiaoLab/suzhuojie/software/kraken2/database/
nohup kraken2 --quick --paired --db /storage/zhangyanxiaoLab/suzhuojie/software/kraken2/database/ \
    --classified-out cseqs#.fq \
    /storage/zhangyanxiaoLab/fastq/2024/2024-06-28-Lianchuan-CKJ/CKJ010/CKJ010_R1.fq.gz \
    /storage/zhangyanxiaoLab/fastq/2024/2024-06-28-Lianchuan-CKJ/CKJ010/CKJ010_R2.fq.gz \
    --output /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/test/kraken2/CKJ010_report.txt \
    --report /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/test/kraken2/CKJ010_species_report.txt 2>&1>/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/logs/CKJ010_kraken.log &