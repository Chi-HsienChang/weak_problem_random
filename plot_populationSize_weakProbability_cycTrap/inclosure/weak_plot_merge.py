import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# 定義原始數據
data1 = [1401389, 801530, 27609, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
data2 = [1398594, 808113, 29372, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
data3 = [1375879, 803391, 30634, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# 設定 X 軸的刻度範圍
x_labels_adjusted = list(range(2, 16))

# 創建圖表並設定大小
plt.figure(figsize=(14, 8))

# 繪製折線圖，使用藍色系顏色
plt.plot(x_labels_adjusted, data1[0:], label='Without GHC', marker='o', markersize=8, linewidth=2, color='lightblue')
plt.plot(x_labels_adjusted, data2[0:], label='One-bit GHC', marker='s', markersize=8, linewidth=2, color='royalblue')
plt.plot(x_labels_adjusted, data3[0:], label='Two-bit GHC', marker='^', markersize=8, linewidth=2, color='navy')

# 添加圖例
plt.legend(fontsize=20)

# 設定圖表標題和軸標籤
plt.title('Two-bit GHC', fontsize=20)  # 增大標題字體大小
plt.xlabel('Bits', fontsize=20)  # 增大X軸標籤字體大小
plt.ylabel('Weak', fontsize=20)  # 增大Y軸標籤字體大小
plt.xticks(fontsize=18)  # 增大X軸刻度數字大小
plt.yticks(fontsize=18)  # 增大Y軸刻度數字大小
plt.grid(True, linestyle='--')  # 加入網格
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) 

# plt.ylim(0, 1.6e6)  # 設定Y軸上限
plt.gca().yaxis.get_offset_text().set_fontsize(18)

# 顯示圖表
plt.savefig('merge2.png')
