import sys
import glob 
genome = 'mm10'
baseurl = 'http://10.28.0.127/share/suzhuojie/Aging_CUT_Tag/HiC/'
tissues = ["brain","CB", "kidney", "liver", "lung", "bonemarrow", "colon", "heart", "Hip", "mammarygland", "stomach", "thymus"]

outfile = '/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/IGV_session/summarizing_IGV_session/HiC_comparment.xml'
with open(outfile,'w') as f:
    f.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?> <Session genome="'+genome+'" hasGeneTrack="true" hasSequenceTrack="true" version="8">\n')
    f.write('  <Resources>\n')
    for tissue in tissues:
        for sfile in sorted(glob.glob(tissue+"/compartment/homer_compartment/PC1/*_50000.PC1.bedGraph")):
            f.write('    <Resource path="'+baseurl+"/" + sfile + '"/>\n')
    f.write('  </Resources>\n')
    f.write('</Session>\n')