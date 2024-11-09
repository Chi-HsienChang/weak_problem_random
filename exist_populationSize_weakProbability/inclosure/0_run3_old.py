from math import comb
import os
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

data = {}
for i in range(sample_size):
    filename = f"../txt_result/total_GHC_{GHC}_sample_{i+1}_ell_{ell}.txt"
    if os.path.exists(filename):
        with open(filename, 'r') as file:
            line = file.read().strip()
            values = list(map(int, line.split()[:]))  # 轉換每個值為整數，忽略第一個數據（例如24）
            data[i+1] = values
    else:
        print(f"File not found: {filename}")


# 打印結果以確認
for sample, values in data.items():
    print(f"Sample {sample}: {values}")

populations = []
# proportions = {2: [], 3: [], 4: [], 5: [], 6: []}
proportions = {i: [] for i in range(2, ell)}

# 計算比例
for sample, counts in data.items():
    pop_size = 2 ** sample
    populations.append(pop_size)
    print(f"Population Size 2^{sample}: {pop_size}, Counts: {counts}")

    for i in range(1, ell-1):  # 處理 weak2 至 weak6, i 從 1 到 5
        # print("Hi")
        if i + 1 <= len(counts) and i + 1 <= ell - 1:
            comb_value = comb(ell - 1, i + 1)
            # print("comb_value: ", comb_value)
            if comb_value > 0 and counts[i] >= 0:
                # print(" (comb_value * total_time): ",  (comb_value * total_time))
                proportion = counts[i] / (comb_value * total_time)
                proportions[i + 1].append(proportion)
                print(f"Weak{i+1} - Counts: {counts[i]}, Comb({ell-1},{i+1}) = {comb_value}, Proportion: {proportion}")
            else:
                proportions[i + 1].append(0)  # 處理數據為零的情況
                print(f"Weak{i+1} - Counts: {counts[i]}, Comb({ell-1},{i+1}) = {comb_value}, Proportion: 0")
        else:
            proportions[i + 1].append(0)  # 若組合數無法計算，設為 0
            print(f"Weak{i+1} - Counts: {counts[i]}, Comb({ell-1},{i+1}) = N/A, Proportion: 0")

# 繪圖
fig, ax = plt.subplots()
for weak in range(2, ell):
    if proportions[weak]:  # 確保比例列表不是空的
        ax.plot(populations, proportions[weak], marker='o', label=f'Weak{weak}')
    else:
        print(f"No valid data for Weak{weak}, skipping.")

ax.set_xscale('log', base=2)
# plt.xticks(range(1, ell+1), labels=[str(2**x) for x in range(ell, 0, -1)], fontsize=20)
# ax.set_xticks([no for no in range(len(populations))])
ax.set_xlabel('Log Population Size')
ax.set_ylabel('Proportion')
ax.set_title('Proportions of Weak Categories vs. Population Size')
ax.legend(loc='upper right')
plt.gca().invert_xaxis()  # 顛倒 X 軸以匹配題目要求 (2^5 到 2^1)
# plt.legend(fontsize = 18)  # 显示图例
plt.grid(True, linestyle='--', linewidth=0.5, color='gray')

plt.tight_layout()  # 调整布局
plt.savefig(f'../png_result/GHC_{GHC}_ell_{ell}_run_{total_time}_date_{datetime.date.today()}.png')  # 保存图像
print(f"Saved plot to ../png_result/GHC_{GHC}_ell_{ell}_run_{total_time}_date_{datetime.date.today()}.png")
