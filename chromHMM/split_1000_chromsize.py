import sys  
# input_file = sys.argv[1]  
# output_file = input_file.replace('.bed', '_1k.bed') 
import os  
import shutil  
  
# 定义输入文件夹和输出文件夹的路径  
input_folder = "/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/all/ChromHMM/15_until_skin/"  
output_folder = "/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/all/ChromHMM/15_until_skin/split_1k/"  
if not os.path.exists(output_folder):  
    os.makedirs(output_folder)  

for root, dirs, files in os.walk(input_folder):  
    for file in files:  
        if file.endswith("segments.bed"):  
            input_file_path = os.path.join(root, file)  
            output_file_path = os.path.join(output_folder, file.replace(".bed", "_1k.bed"))  
  
            with open(input_file_path, 'r') as file, open(output_file_path, 'w') as output:  
                for line in file:  
                    parts = line.split('\t')  
                    start = int(parts[1])  
                    end = int(parts[2])  
  
                    if end - start != 1000:  
                        num_lines = (end - start) // 1000  
                        for i in range(num_lines):  
                            new_start = start + i * 1000  
                            new_end = start + (i + 1) * 1000
                            output.write(f"{parts[0]}\t{new_start}\t{new_end}\t{parts[3]}")  
                    else:  
                        output.write(line) 