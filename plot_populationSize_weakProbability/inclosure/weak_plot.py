import matplotlib.pyplot as plt

# 提供的數據
data = [1375879, 803391, 30634, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# 更新數據點的標籤為2到15
new_labels = [str(i+1) for i in range(len(data))]

# 重新繪製圖表，將網格改為虛線並增大X軸、Y軸的數字大小以及標題的大小
plt.figure(figsize=(20, 6))
bars = plt.bar(new_labels, data, color='skyblue')
plt.title('Two-bit GHC', fontsize=20)  # 增大標題字體大小
plt.xlabel('Bits', fontsize=20)  # 增大X軸標籤字體大小
plt.ylabel('Weak', fontsize=20)  # 增大Y軸標籤字體大小
plt.xticks(fontsize=18)  # 增大X軸刻度數字大小
plt.yticks(fontsize=18)  # 增大Y軸刻度數字大小
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))  # 使用科學記號
plt.grid(True, linestyle='--')  # 加入網格
plt.ylim(0, 1.6e6)  # 設定Y軸上限


plt.gca().yaxis.get_offset_text().set_fontsize(18)

# 在柱狀上顯示數據
for bar in bars:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, yval, f'{yval:.0f}', ha='center', va='bottom', fontsize=16, fontweight='bold')


# 再次保存圖片
plt.savefig('Twobit_GHC_16bits_64000times.png')