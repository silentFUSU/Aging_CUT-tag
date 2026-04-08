import sys
import glob 

genome = 'mm10'
baseurl = 'http://10.28.0.127/share/suzhuojie/Aging_CUT_Tag/'
tissues = ['aorta','BAT','bladder','bonemarrow','brain','CB','cecum','colon','heart','Hip','ileum','iWAT','jejunum','kidney','liver','lung','mammarygland','muscle','ovary','pancreas','skin','spleen','stomach','testis','thymus','tongue','uterus','iWAT']
antibodys = ['H3K27me3',"H3K9me3","H3K36me3"]
for antibody in antibodys:
    outfile = '/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/IGV_session/summarizing_IGV_session/'+antibody+'_log2FC.xml'
    with open(outfile,'w') as f:
        f.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?> <Session genome="'+genome+'" hasGeneTrack="true" hasSequenceTrack="true" version="8">\n')
        f.write('  <Resources>\n')
        for tissue in tissues:
            for sfile in sorted(glob.glob(tissue+"/"+antibody+"/bw/old_young_log2ratio.bw")):
                f.write('    <Resource path="'+baseurl+"/" + sfile + '"/>\n')
        f.write('  </Resources>\n')
        f.write('</Session>\n')


antibodys = ["ATAC"]
baseurl = 'http://10.28.0.127/share/suzhuojie/Aging_CUT_Tag/ATAC/'
for antibody in antibodys:
    outfile = '/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/IGV_session/summarizing_IGV_session/'+antibody+'.xml'
    with open(outfile,'w') as f:
        f.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?> <Session genome="'+genome+'" hasGeneTrack="true" hasSequenceTrack="true" version="8">\n')
        f.write('  <Resources>\n')
        for tissue in tissues:
            for sfile in sorted(glob.glob(tissue+"/"+antibody+"/bw/*nodup.bw")):
                f.write('    <Resource path="'+baseurl+"/" + sfile + '"/>\n')
        f.write('  </Resources>\n')
        f.write('</Session>\n')

outfile = '/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/IGV_session/summarizing_IGV_session/H3K9me3_H3K4me3.xml'
with open(outfile,'w') as f:
    f.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?> <Session genome="'+genome+'" hasGeneTrack="true" hasSequenceTrack="true" version="8">\n')
    f.write('  <Resources>\n')
    for tissue in tissues:
        for sfile in sorted(glob.glob(tissue+"/H3K9me3/bw/*bs1000.bw")):
            f.write('    <Resource path="'+baseurl+"/" + sfile + '"/>\n')
        for sfile in sorted(glob.glob(tissue+"/H3K4me3/bw/*.nodup.bw")):
            f.write('    <Resource path="'+baseurl+"/" + sfile + '"/>\n')
    f.write('  </Resources>\n')
    f.write('</Session>\n')

antibodys = ["RNA"]
baseurl = 'http://10.28.0.127/share/suzhuojie/Aging_CUT_Tag/RNA/'
for antibody in antibodys:
    outfile = '/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/IGV_session/summarizing_IGV_session/'+antibody+'.xml'
    with open(outfile,'w') as f:
        f.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?> <Session genome="'+genome+'" hasGeneTrack="true" hasSequenceTrack="true" version="8">\n')
        f.write('  <Resources>\n')
        for tissue in tissues:
            for sfile in sorted(glob.glob(tissue+"/bw/*sorted.bw")):
                f.write('    <Resource path="'+baseurl+"/" + sfile + '"/>\n')
        f.write('  </Resources>\n')
        f.write('</Session>\n')