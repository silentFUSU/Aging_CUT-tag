data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/NIH_Roadmap_BJ/
# samples=(SRR13274560 SRR13274561)
samples=(ENCFF745TNC ENCFF265JIG)
mkdir -p ${data_path}peaks/
mkdir -p ${data_path}peaks/macs_narrowpeak
for sample in ${samples[@]}
do
    macs2 callpeak -B --SPMR --broad --nomodel --nolambda -t ${data_path}bam/${sample}.bam  -f BAM -n ${sample} --outdir ${data_path}peaks/macs_narrowpeak -g hs -q 0.0001  --keep-dup all &
    # macs2 callpeak -B --SPMR --nomodel -t ${data_path}bam/${sample}.bam  -f BAMPE -n ${sample} --outdir ${data_path}peaks/macs_narrowpeak -g hs -q 0.0001  --keep-dup all &
done

wait
for sample in ${samples[@]}
do
    python /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/public_WANG_HiC/bedtobins.py -i ${data_path}peaks/macs_narrowpeak/${sample}_treat_pileup.bdg -o ${data_path}peaks/macs_narrowpeak/${sample}_10kb.txt &
done