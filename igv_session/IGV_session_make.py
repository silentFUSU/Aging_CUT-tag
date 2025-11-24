import sys
import glob 

genome = 'mm10'
baseurl = 'http://10.28.0.127/share/suzhuojie/Aging_CUT_Tag/'
# tissues=['aorta','BAT','bladder','bonemarrow','brain','CB','cecum','colon','heart','Hip','ileum','iWAT','jejunum','kidney','liver','lung','mammarygland','muscle','ovary','pancreas','skin','spleen','stomach','testis','thymus','tongue','uterus']
tissues=['lung']
for tissue in tissues:
  outfile = '/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/IGV_session/summarizing_IGV_session/'+tissue+'_HiC_WGBS.xml'
  with open(outfile,'w') as f:
    f.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?> <Session genome="'+genome+'" hasGeneTrack="true" hasSequenceTrack="true" version="8">\n')
    f.write('  <Resources>\n')
    ## files 
    # for sfile in sorted(glob.glob(tissue+"/bw/*CpG.bw")):
    #   f.write('    <Resource path="'+baseurl+"/" + sfile + '"/>\n')
    for sfile in sorted(glob.glob(tissue+"/H3K9me3/bw/*bs1000.bw")):
      f.write('    <Resource path="'+baseurl+"/" + sfile + '"/>\n')
    for sfile in sorted(glob.glob(tissue+"/H3K36me3/bw/*bs1000.bw")):
      f.write('    <Resource path="'+baseurl+"/" + sfile + '"/>\n')
    for sfile in sorted(glob.glob(tissue+"/H3K27me3/bw/*bs1000.bw")):
      f.write('    <Resource path="'+baseurl+"/" + sfile + '"/>\n')
    for sfile in sorted(glob.glob(tissue+"/H3K27ac/bw/*.nodup.bw")):
      f.write('    <Resource path="'+baseurl+"/" + sfile + '"/>\n')
    for sfile in sorted(glob.glob(tissue+"/H3K4me3/bw/*.nodup.bw")):
      f.write('    <Resource path="'+baseurl+"/" + sfile + '"/>\n')
    for sfile in sorted(glob.glob(tissue+"/H3K4me1/bw/*.nodup.bw")):
      f.write('    <Resource path="'+baseurl+"/" + sfile + '"/>\n')
    for sfile in sorted(glob.glob("ATAC/"+tissue+"/ATAC/bw/*.nodup.bw")):
      f.write('    <Resource path="'+baseurl+"/" + sfile + '"/>\n')
    for sfile in sorted(glob.glob("WGBS/"+tissue+"/bw/*_CpG.bw")):
      f.write('    <Resource path="'+baseurl+"/" + sfile + '"/>\n')
    for sfile in sorted(glob.glob("HiC/"+tissue+"/compartment/homer_compartment/PC1/*_50000.PC1.bedGraph")):
      f.write('    <Resource path="'+baseurl+"/" + sfile + '"/>\n')
    # for sfile in sorted(glob.glob("*_enhancer.bed")):
    #   f.write('    <Resource path="'+baseurl+"/" + sfile + '"/>\n')
    f.write('  </Resources>\n')
    f.write('</Session>\n')