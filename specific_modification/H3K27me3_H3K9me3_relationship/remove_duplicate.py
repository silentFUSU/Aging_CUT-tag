import sys
input_file = sys.argv[1]

with open(input_file, 'r') as file:  
    lines = file.readlines()  
    unique_lines = set(lines)

output_file = input_file.split('.')[0] + "_unique.bed"  
with open(output_file, 'w') as file:  
    file.writelines(unique_lines) 
