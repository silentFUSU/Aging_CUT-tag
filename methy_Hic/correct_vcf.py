import sys 
# 读取VCF文件  
input_vcf_file = sys.argv[1]  

with open(input_vcf_file, 'r') as file:  
    lines = file.readlines()  

# 处理行拆分  
corrected_lines = []  
  
for line in lines:  
    if not line.strip(): 
        continue 
    if line.startswith('#'):  
        corrected_lines.append(line)  
        continue  
    data = line.split('\t')
    if len(data) > 10:  # 如果第十个元素存在  
        elements = data[9].split('chr', 1)  # 以"chr"为分隔符拆分第十个元素          
        if len(elements) == 2:  # 如果成功拆分成两部分  
            # new_line = '\t'.join(data[:9] + [f"chr{elements[1]}"] + elements[0].split('chr', 1)[1:])  # 组合新的一行  
            new_line = '\t'.join(data[:9]  + elements[0].split('chr', 1)[:1])
            new_line += '\n'
            corrected_lines.append(new_line)  
            new_line= '\t'.join([f"chr{elements[1]}"] + data[10:])
            corrected_lines.append(new_line)  
    else:
        corrected_lines.append(line)  

last_part = input_vcf_file.rsplit('/', 1)[-1]  
  
corrected_last_part = "corrected_" + last_part  
  
first_part = input_vcf_file.rsplit('/', 1)[0]  
  
output_file = first_part + '/' + corrected_last_part 

with open(output_file, 'w') as file:  
    file.writelines(corrected_lines)  