import os

def count_lines_in_bed_files(directory):
    # 初始化结果列表
    results = []

    # 遍历目录下所有文件
    for filename in os.listdir(directory):
        # 检查文件扩展名是否为.bed
        if filename.endswith(".bed"):
            # 构建完整的文件路径
            filepath = os.path.join(directory, filename)
            # 计算文件行数
            with open(filepath, 'r') as file:
                line_count = sum(1 for line in file)
            # 将文件名和行数以元组形式添加到列表中
            results.append((filename, line_count))

    # 将结果保存到文本文件
    with open(directory+"bed_files_line_counts.txt", "w") as output_file:
        for filename, line_count in results:
            output_file.write(f"{filename}\t{line_count}\n")

    print("Line counts have been written to bed_files_count.txt")

# 使用函数
directory = '/storage/zhangyanxiaoLab/suzhuojie/projects/Aging_CUT_Tag/result/all/ChromHMM/15_until_skin/state_transfer/E13_to_E6/bed/'  # 替换为你的目录路径
count_lines_in_bed_files(directory)