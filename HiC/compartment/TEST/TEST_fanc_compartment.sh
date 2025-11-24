data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/HiC/
tissue=lung
resolution=1mb
young=${data_path}${tissue}/juicer/young_combined.allValidPairs.hic
old=${data_path}${tissue}/juicer/old_combined.allValidPairs.hic
ref=/storage/zhangyanxiaoLab/suzhuojie/ref_data/for_normal_mapping/mm10/mm10.fa
mkdir ${data_path}${tissue}/compartment

fanc compartments -g ${ref} \
                  -d ${data_path}${tissue}/compartment/fanc_young_${resolution}.compartment.bed \
                  -v ${data_path}${tissue}/compartment/fanc_young_${resolution}.ev.txt \
                  ${young}@${resolution}  ${data_path}${tissue}/compartment/fanc_young_${resolution}.ab &

fanc compartments -g ${ref} \
                  -d ${data_path}${tissue}/compartment/fanc_old_${resolution}.compartment.bed \
                  -v ${data_path}${tissue}/compartment/fanc_old_${resolution}.ev.txt \
                  ${old}@${resolution}  ${data_path}${tissue}/compartment/fanc_old_${resolution}.ab &

wait
echo all done