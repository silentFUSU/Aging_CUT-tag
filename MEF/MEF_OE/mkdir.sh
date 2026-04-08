data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/MEF_OE/
antibodys=(H3K27me3 H2AK119ub1)
conditions=(MEF_Vector MEF_Bmi1 MEF_Cbx2 MEF_Cbx7)
dirs=(bam bw bed)
for antibody in ${antibodys[@]}
do
    mkdir -p ${data_path}${antibody}
    for condition in ${conditions[@]}
    do
        mkdir -p ${data_path}${antibody}/${condition}
        for dir in ${dirs[@]}
        do
            mkdir -p ${data_path}${antibody}/${condition}/${dir}
        done
    done
done