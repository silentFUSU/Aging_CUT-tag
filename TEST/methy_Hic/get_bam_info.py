import pysam  

def extract_bam_info(input_bam, output_txt):  
    bamfile = pysam.AlignmentFile(input_bam, "rb")      
    with open(output_txt, 'w') as outfile:  
        outfile.write("ReadName\tReadLength\tChromosome\tPosition\tCIGAR\tStrand\n")  
        for read in bamfile:  
            if not read.is_unmapped:  # 过滤未对齐的reads  
                read_name = read.query_name  
                chromosome = bamfile.get_reference_name(read.reference_id)  
                read_length = read.query_length  
                position = read.reference_start + 1  # 转化为1-based坐标  
                cigar = read.cigarstring  
                strand = '-' if read.is_reverse else '+'  
                outfile.write(f"{read_name}\t{read_length}\t{chromosome}\t{position}\t{cigar}\t{strand}\n")
    bamfile.close()  

bam_file_path = "DYQ035_mapped_unique.bam"  
output_txt = "DYQ035_mapped_unique.txt"  
extract_bam_info(bam_file_path,output_txt)
