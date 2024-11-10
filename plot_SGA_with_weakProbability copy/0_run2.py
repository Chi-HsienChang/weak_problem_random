import sys
import os
import matplotlib.pyplot as plt
import math
import datetime

if __name__ == '__main__':

    # ell = int(sys.argv[1])
    # total_time = int(sys.argv[2])
    # subtasks = int(sys.argv[3])

    # total_time = int(sys.argv[1])
    # subtasks = int(sys.argv[2])
    # ell = int(sys.argv[3])
    # nInitial = int(sys.argv[4])       
    # selectionPressure = int(sys.argv[5])
    # pc = float(sys.argv[6])
    # pm = float(sys.argv[7])
    # maxGen = int(sys.argv[8])
    # maxFe = int(sys.argv[9])


    ell = int(sys.argv[1])
    total_time = int(sys.argv[2])
    subtasks = int(sys.argv[3])
    nInitial = int(sys.argv[4])       
    selectionPressure = int(sys.argv[5])
    pc = float(sys.argv[6])
    pm = float(sys.argv[7])
    maxGen = int(sys.argv[8])
    maxFe = int(sys.argv[9])

    total_sums = []  # 用于存储每行加总后的结果

    for s in range(0, subtasks):
        dir_name = f'./SGA_result/syn_{s}/'  # 构建目录名
        file_path = os.path.join(dir_name, 'a_total_weak.txt')  # 构建文件的完整路径
        
        try:
            with open(file_path, 'r') as file:  # 尝试打开文件进行读取
                for index, line in enumerate(file):  # 逐行读取
                    values = line.split()  # 将每行分割成数值列表
                    if index >= len(total_sums):  # 如果total_sums列表中的行数不够，添加一个新行
                        total_sums.append([0] * len(values))  # 使用0初始化新行
                    
                    for i, value in enumerate(values):  # 遍历数值列表
                        comb_ = math.comb(ell-1, i+1)
                        if comb_ != 0:
                            total_sums[index][i] += float(value)/(comb_*total_time) # 将每个数值加到对应的总和中
        except FileNotFoundError:
            print(f"File not found: {file_path}")  # 如果文件不存在，打印错误信息
        except ValueError:
            print("Error converting value to integer.")  # 如果无法将值转换为整数，打印错误信息

    col_count = len(total_sums[0]) if total_sums else 0  # 获取列数

    # 創建一個新的圖表
    plt.figure(figsize=(10, 6))

    # 找到所有数据中的最小值和最大值
    all_values = [item for sublist in total_sums for item in sublist]
    min_value = min(all_values)
    max_value = max(all_values)

    # 為每一列（除了第一列）繪製一條數據線
    for col in range(1, col_count):  # 從第二列開始迭代，跳過第一列
        col_values = [row[col] for row in total_sums]  # 获取每一行的当前列的值
        plt.plot(col_values, label=f'Weak {col+1}', marker='.')  # 为每一条线添加标签，跳过第一周

    # plt.title('Weak Epistasis over Generations', fontsize=20)
    plt.xlabel('Generation', fontsize=20)
    plt.ylabel('Probability', fontsize=20)
    # plt.ylim([min_value - 1.2, 12])
    plt.grid(True, linestyle='--', linewidth=0.5, color='gray')
    # plt.grid(True, linestyle='--')
    plt.legend(fontsize=15)  # 顯示圖例
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)

    plt.tight_layout()  # 调整布局
    plt.savefig(f'./png_result/SGA_ell_{ell}_run_{total_time}_date_{datetime.date.today()}.png')  # 保存圖表到文件
