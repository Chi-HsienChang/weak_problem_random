from math import comb
import os
import pandas as pd
import sys
import matplotlib.pyplot as plt
import datetime

# 參數設定
ell = int(sys.argv[1])
total_time = int(sys.argv[2])
subtasks = int(sys.argv[3])
debug = int(sys.argv[4])       
GHC = int(sys.argv[5])
sample_size = int(sys.argv[6])

# 建立空的 DataFrame
data = pd.DataFrame()

# 讀取每個樣本檔案並加入 DataFrame
for i in range(1, sample_size + 1):
    filename = f"../txt_result/total_GHC_{GHC}_sample_{i}_ell_{ell}.txt"
    if not os.path.isfile(filename):
        print(f"File {filename} not found, skipping...")
        continue

    # 讀取檔案內容
    with open(filename, 'r') as file:
        counts = list(map(int, file.read().split()))

    # data[f'popSize_2^{i}'] = counts

    # 計算比例並加入 DataFrame
    proportions = [
        # float(count) / (comb(ell - 1, k + 1) * total_time) if comb(ell - 1, k + 1) != 0 else 0
        float(count) / 1 if 1 != 0 else 0
        for k, count in enumerate(counts)
    ]
    data[f'{2**i}'] = proportions

# 檢查結果 DataFrame
print(data)



import matplotlib.pyplot as plt

# 排除索引為 0 的行數據
data_filtered = data.iloc[1:]  # 從索引 1 開始的所有行

# 繪製折線圖
plt.figure(figsize=(10, 6))
for idx, row_data in data_filtered.iterrows():
    plt.plot(row_data.index, row_data.values, marker='o', label=f'weak {idx+1}')

# 設定圖表屬性
plt.xlabel('Population Size', fontsize = 20)
plt.ylabel('Probability', fontsize = 20)
# plt.title('Data Across Different Population Sizes (Excluding Index 0)')
plt.legend()
plt.grid(True, linestyle='--', linewidth=0.5, color='gray')

# 顯示圖表
plt.gca().invert_xaxis()

plt.legend(fontsize=15)  # 顯示圖例
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.tight_layout()


# # 顯示圖表
# plt.tight_layout()
plt.savefig(f'../png_result/GHC_{GHC}_ell_{ell}_run_{total_time}_date_{datetime.date.today()}.png') 
