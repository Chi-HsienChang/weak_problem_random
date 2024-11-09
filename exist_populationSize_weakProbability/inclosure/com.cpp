#include <iostream>
#include <vector>

std::vector<std::vector<int>> generateBinarySequences(int n) {
    int totalSequences = 1 << n;  // 2^n
    std::vector<std::vector<int>> allSequences;

    for (int i = 0; i < totalSequences; ++i) {
        std::vector<int> sequence(n);

        for (int j = 0; j < n; ++j) {
            // 檢查第j位是否為1
            sequence[n - 1 - j] = (i >> j) & 1;
        }

        // 將序列添加到所有序列的vector中
        allSequences.push_back(sequence);
    }

    return allSequences;
}

int main() {
    int n = 3;  // 序列的長度
    auto sequences = generateBinarySequences(n);

    // 打印所有序列
    for (const auto& sequence : sequences) {
        for (int bit : sequence) {
            std::cout << bit;
        }
        std::cout << std::endl;
    }

    return 0;
}



// #include <iostream>
// #include <vector>
// #include <algorithm>

// // 修改函數以返回組合的vector
// std::vector<std::vector<int>> generateCombinations(int n, int k) {
//     std::vector<int> bitmask(k, 1);  // 創建k個1
//     bitmask.resize(n, 0);  // 後面填充n-k個0

//     std::vector<std::vector<int>> combinations;

//     do {
//         std::vector<int> currentCombination;
//         for (int i = 0; i < n; ++i) {
//             if (bitmask[i]) {
//                 currentCombination.push_back(i + 1);  // 輸出選中的元素
//             }
//         }
//         combinations.push_back(currentCombination);
//     } while (std::prev_permutation(bitmask.begin(), bitmask.end()));  // 生成下一個排列

//     return combinations;
// }

// int main() {
//     int n = 6;  // 從1到n
//     int k = 3;  // 選擇k個元素
//     auto combinations = generateCombinations(n, k);

//     // 印出所有組合
//     for (const auto& combination : combinations) {
//         for (int num : combination) {
//             std::cout << num << " ";
//         }
//         std::cout << std::endl;
//     }

//     return 0;
// }
