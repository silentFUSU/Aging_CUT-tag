data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/public_data/m6a_monkey_CRA005942/CRA005942/
samples=($(basename -a ${data_path}CRR*/))
input_samples=()  
IP_samples=()  
  
for ((i=0; i<${#samples[@]}; i++)); do  
    if ((i % 2 == 0)); then  
        input_samples+=(${samples[i]})  
    else  
        IP_samples+=(${samples[i]})  
    fi  
done  

echo ${input_samples[@]}
echo ${IP_samples[@]}
if [ ! -d ${data_path}peaks ]; then  
    mkdir ${data_path}peaks
fi

if [ ! -d ${data_path}bed ]; then  
    mkdir ${data_path}bed
fi
for ((i=0; i<${#IP_samples[@]}; i++)); do  
    echo ${IP_samples[i]} ${input_samples[i]}
    macs2 callpeak -t ${data_path}${IP_samples[i]}/${IP_samples[i]}.bam -c ${data_path}${input_samples[i]}/${input_samples[i]}.bam -n ${IP_samples[i]} --outdir ${data_path}peaks  --nomodel -q 0.00001  --keep-dup all 
    awk '{OFS="\t"} {print $1, $2, $3, "peaks"NR}' ${data_path}peaks/${IP_samples[i]}_peaks.narrowPeak |Rscript /storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/code/keep_regular_chroms.r > ${data_path}bed/${IP_samples[i]}_peaks.bed
done  

if [ ! -d ${data_path}bw ]; then  
    mkdir ${data_path}bw
fi

for ((i=0; i<${#IP_samples[@]}; i++)); do  
    echo ${IP_samples[i]} ${input_samples[i]}
    samtools index ${data_path}${IP_samples[i]}/${IP_samples[i]}.bam 
    samtools index ${data_path}${input_samples[i]}/${input_samples[i]}.bam
    bamCompare -b1 ${data_path}${IP_samples[i]}/${IP_samples[i]}.bam -b2 ${data_path}${input_samples[i]}/${input_samples[i]}.bam \
        -bs 50 -p 16 -o ${data_path}bw/${IP_samples[i]}.bw 
done  