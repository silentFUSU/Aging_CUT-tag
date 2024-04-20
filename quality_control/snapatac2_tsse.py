import snapatac2 as snap
import glob 
import os  
import anndata as ad
import gzip
data_path = '/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/data/'
result_path = '/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/'
antibodys =['ATAC','H3K27ac','H3K4me3']
# tissues = ["brain","liver", "testis", "colon", "kidney", "lung", "spleen", "muscle", "pancreas"]
# tissues = ["brain"]
# tissues = ["testis", "colon", "kidney", "lung", "spleen", "muscle", "pancreas"]
tissues = ["skin"]


for antibody in antibodys:
    if not os.path.exists("{}all/tsse/{}".format(result_path,antibody)):
        os.makedirs("{}all/tsse/{}".format(result_path,antibody))
    outfile = "{}all/tsse/{}/tsse.txt".format(result_path,antibody)
    for tissue in tissues:
        if antibody == 'ATAC':
            glob_path = "{}samples/ATAC/{}/{}/bam/*.bam".format(data_path,tissue,antibody)
        else:
            glob_path = "{}samples/{}/{}/bam/*.bam".format(data_path,tissue,antibody)
        for f in glob.glob(glob_path):
            f_name = f.rsplit('/',1)[-1].split('.')[0]
            input_fragment_file="{}all/tsse/{}/{}_{}.tsv".format(result_path,antibody,tissue,f_name)
            snap.pp.make_fragment_file(f,input_fragment_file,barcode_regex = r"^(.*?):\d+.*$")
            input_chrom_size={"chr1": 195471971, "chr2": 182113224, "chr3": 160039680, "chr4": 156508116, "chr5": 151834684, "chr6": 149736546, "chr7": 145441459, "chr8": 129401213, "chr9": 124595110, "chr10": 130694993, "chr11": 122082543, "chr12": 120129022, "chr13": 120421639, "chr14": 124902244, "chr15": 104043685, "chr16": 98207768, "chr17": 94987271, "chr18": 90702639, "chr19": 61431566, "chrX": 171031299, "chrY": 91744698} 
            gzip_fragment_file="{}all/tsse/{}/{}_{}.tsv.gz".format(result_path,antibody,tissue,f_name)
            with open(input_fragment_file, 'rb') as f_in, gzip.open(gzip_fragment_file, 'wb') as f_out:  
                f_out.writelines(f_in) 
            data = snap.pp.import_data(
                fragment_file=gzip_fragment_file,
                chrom_sizes=snap.genome.mm10,
                sorted_by_barcode=False
            )
            snap.metrics.tsse(data, snap.genome.mm10)
            with open(outfile,'a') as out_f:
                out_f.write(tissue+' '+f_name+' ')
            tsse_column = data.obs['tsse']
            tsse_column.to_csv(outfile, index=False, mode='a', header=False)
