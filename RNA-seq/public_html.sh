data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/RNA/
tissues=(jejunum ileum colon cecum)
mkdir ~/public_html/Aging_CUT_Tag/RNA/
for tissue in ${tissues[@]}
do
    mkdir ~/public_html/Aging_CUT_Tag/RNA/${tissue}
    mkdir ~/public_html/Aging_CUT_Tag/RNA/${tissue}/bw
    ln -s ${data_path}${tissue}/bw/* ~/public_html/Aging_CUT_Tag/RNA/${tissue}/bw

done