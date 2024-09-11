data_path=/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/samples/
tissue=CB
antibody=H3K4me3

bedtools intersect -a ${data_path}${tissue}/${antibody}/bed/${antibody}_1kb_bins_diff_after_remove_batch_effect_down.bed\
            -b ~/projects/Aging_CUT_Tag/result/all/ChromHMM/until_ovary/15_until_ovary/state_transfer/E5_to_E9/bed/${tissue}.bed  > ~/projects/Aging_CUT_Tag/result/all/ChromHMM/until_ovary/15_until_ovary/state_transfer/E5_to_E9/bed/${tissue}_overlap_${antibody}_down.bed 