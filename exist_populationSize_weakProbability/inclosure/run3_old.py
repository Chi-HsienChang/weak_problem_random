import matplotlib.pyplot as plt
import math
import sys
import datetime

ell = int(sys.argv[1])
total_time = int(sys.argv[2])
subtasks = int(sys.argv[3])
debug = int(sys.argv[4])       
GHC = int(sys.argv[5])
sample_size = int(sys.argv[6])

data_points = [[] for _ in range(ell-1)]  # [weak1, weak2, ..., weak(ell-1)]

# 循环读取每个文件
for i in range(sample_size):  # 从1到12
    filename = f"../txt_result/total_GHC_{GHC}_sample_{i+1}_ell_{ell}.txt"
    with open(filename, 'r') as file:
        line = file.readline()
        numbers = [int(num) for num in line.split()]
        # 将每个数据点的值添加到对应的列表中
        for j in range(len(numbers)):
            data_points[i].append(float(numbers[j]))


# 设置画布大小
plt.figure(figsize=(10, 8))  # 修改尺寸以适应单个图表

# 在同一张图上绘制所有数据点的折线图
for i in range(len(data_points)-1, 0): # i = 1 對應到 weak 2 所以是 C 的 ell-1 取 i+1
    if data_points[i]:  # 确保数据点列表不为空
        comb_ = math.comb(ell-1, i+1) # C的幾取幾、
        # print(data_points[i])
        for j in range(len(data_points[i])):
            data_points[i][j] = float(data_points[i][j])/(comb_*total_time)
        # print(data_points[i])
        plt.plot(range(ell-1, 0, -1), data_points[i], marker='o', linestyle='-', label=f'Weak {i+1}')

# plt.title('Weak Number by Epistasis Size and Sample Size', fontsize = 20)
plt.xlabel('Population Size', fontsize = 20)
plt.ylabel('Probability', fontsize = 20)
plt.xticks(range(1, ell+1), labels=[str(2**x) for x in range(ell, 0, -1)], fontsize=20)
plt.yticks(fontsize = 20)
# plt.ylim(0, 0.12)  # 设置y轴范围为0到max_value+1
plt.legend(fontsize = 18)  # 显示图例
plt.grid(True, linestyle='--', linewidth=0.5, color='gray')

plt.tight_layout()  # 调整布局
plt.savefig(f'../png_result/GHC_{GHC}_ell_{ell}_run_{total_time}_date_{datetime.date.today()}.png')  # 保存图像




