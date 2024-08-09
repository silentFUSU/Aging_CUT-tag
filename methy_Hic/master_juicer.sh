data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/raw_data/20240804_WGBS/juicer_meth_1M_trim/
juicer_path=/storage/zhangyanxiaoLab/suzhuojie/software/juicer/
samples=(XX315)

ref=mm10_lambda
site=DpnII
for sample in ${samples[@]}
do
    nohup ${juicer_path}scripts/juicer.sh -m  -g ${ref} -d ${data_path}${sample} -s ${site} \
        -a ${sample} -p ${juicer_path}references/${ref}/${ref}.chrom.sizes \
        -y ${juicer_path}restriction_sites/${ref}_${site}.txt -z ${juicer_path}references/${ref}/${ref}.fa \
        -D ${juicer_path} -b GATCGATC -t 10 2>&1>${data_path}${sample}/${sample}_juicer_meth.log &
done