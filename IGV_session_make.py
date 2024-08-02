import sys
import glob 

argv = sys.argv

if len(argv)< 3: 
  print("Usage: python script.py genome baseurl outfile")
  exit(1)
genome = argv[1]
baseurl = argv[2]
outfile = argv[3]
genome = 'mm10'
# genome = 'GCF_000364345.1'
# genome = 'hg38'
# baseurl = 'http://10.28.0.127/share/suzhuojie/Aging_CUT_Tag/'
baseurl = 'http://10.28.0.127/share/suzhuojie/Aging_CUT_Tag/muscle/'
outfile = '/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/IGV_session/muscle_all_markers.xml'
# baseurl = 'http://10.28.0.127/share/suzhuojie/cellular_aging_GSE133292'
# outfile = '/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/IGV_session/cellular_aging_GSE133292.xml'
with open(outfile,'w') as f:
  f.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?> <Session genome="'+genome+'" hasGeneTrack="true" hasSequenceTrack="true" version="8">\n')
  f.write('  <Resources>\n')
  ## files 
  for sfile in sorted(glob.glob("H3K27me3/bw/*.nodup.bw")):
    f.write('    <Resource path="'+baseurl+"/" + sfile + '"/>\n')
  for sfile in sorted(glob.glob("H3K9me3/bw/*.nodup.bw")):
    f.write('    <Resource path="'+baseurl+"/" + sfile + '"/>\n')
  for sfile in sorted(glob.glob("H3K27ac/bw/*.nodup.bw")):
    f.write('    <Resource path="'+baseurl+"/" + sfile + '"/>\n')
  for sfile in sorted(glob.glob("H3K36me3/bw/*.nodup.bw")):
    f.write('    <Resource path="'+baseurl+"/" + sfile + '"/>\n')
  for sfile in sorted(glob.glob("H3K4me3/bw/*.nodup.bw")):
    f.write('    <Resource path="'+baseurl+"/" + sfile + '"/>\n')
  for sfile in sorted(glob.glob("H3K4me1/bw/*.nodup.bw")):
    f.write('    <Resource path="'+baseurl+"/" + sfile + '"/>\n')
  # for sfile in sorted(glob.glob("*_enhancer.bed")):
  #   f.write('    <Resource path="'+baseurl+"/" + sfile + '"/>\n')
  f.write('  </Resources>\n')
  f.write('</Session>\n')