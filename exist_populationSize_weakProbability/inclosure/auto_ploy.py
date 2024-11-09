import matplotlib.pyplot as plt

# 檔案路徑
file_paths = ['../MSO_0_GHC/result.txt', '../MSO_1_GHC/result.txt', '../MSO_2_GHC/result.txt']

# 初始化數據列表
data_list = []

# 從檔案中讀取數據
for path in file_paths:
    with open(path, 'r') as file:
        data = file.read().split()  # 將每行數據拆分成列表
        data_list.append([int(num) for num in data])  # 轉換成整數並添加到列表中

# 設定圖表的整體大小和標題
plt.figure(figsize=(20, 18))
plt.suptitle('Weak Epistasis', fontsize=20)

# 為每個數據集繪製子圖
for i, data in enumerate(data_list, start=1):
    ax = plt.subplot(3, 1, i)  # 創建子圖
    new_labels = [str(j+2) for j in range(len(data))]
    bars = ax.bar(new_labels, data, color='skyblue')
    ax.set_title(f'GHC {i-1}', fontsize=20)
    ax.set_xlabel('Bits', fontsize=20)
    ax.set_ylabel('Weak', fontsize=20)
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax.grid(True, linestyle='--')
    ax.set_ylim(0, 1.6e6)
    ax.yaxis.get_offset_text().set_fontsize(18)
    
    # 在柱狀上顯示數據
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, yval, f'{yval:.0f}', ha='center', va='bottom', fontsize=16, fontweight='bold')

# 調整子圖間的間距
plt.tight_layout(rect=[0, 0.03, 1, 0.95])

# 保存圖片
plt.savefig('auto_plot.png')
plt.show()
