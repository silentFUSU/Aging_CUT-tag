import sys
import glob 

genome = 'mm10'
baseurl = 'http://10.28.0.127/share/suzhuojie/Aging_CUT_Tag/'
tissues = ['aorta','BAT','bladder','bonemarrow','brain','CB','cecum','colon','heart','Hip','ileum','iWAT','jejunum','kidney','liver','lung','mammarygland','muscle','ovary','pancreas','skin','spleen','stomach','testis','thymus','tongue','uterus','iWAT']
for antibody in antibodys:
    outfile = '/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/IGV_session/summarizing_IGV_session/WGBS.xml'
    with open(outfile,'w') as f:
        for tissue in tissues:
            f.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?> <Session genome="'+genome+'" hasGeneTrack="true" hasSequenceTrack="true" version="8">\n')
            f.write('  <Resources>\n')
            for sfile in sorted(glob.glob(tissue+"/bw/DYQ*CpG.bw")):
                f.write('    <Resource path="'+baseurl+"/" + sfile + '"/>\n')
            f.write('  </Resources>\n')
            f.write('</Session>\n')