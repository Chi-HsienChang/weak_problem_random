/***************************************************************************
 *   Copyright (C) 2004 by Tian-Li Yu                                      *
 *   tianliyu@cc.ee.ntu.edu.tw                                             *
 *                                                                         *
 *   You can redistribute it and/or modify it as you like                  *
 ***************************************************************************/

#include <cmath>
#include <cstdio>
#include <iostream>
#include <cstdlib>

#include "statistics.h"
#include "ga.h"
#include "chromosome.h"
#include "global.h"

using namespace std;

int main (int argc, char *argv[])
{

    if (argc != 10) {
        printf ("GA ell nInitial selectionPressure pc pm maxGen maxFe repeat subtask\n");
        return -1;
    }

    int ell = atoi (argv[1]);    // problem size
                                 // initial population size
    int nInitial = atoi (argv[2]);
                                 // selection pressure
    int selectionPressure = atoi (argv[3]);
    double pc = atof (argv[4]);  // pc
    double pm = atof (argv[5]);  // pm
    int maxGen = atoi (argv[6]); // max generation
    int maxFe = atoi (argv[7]);  // max fe
    int repeat = atoi (argv[8]); // how many time to repeat
    int subtask = atoi (argv[9]);

    std::cout<<"----------subtask: "<<subtask<<std::endl;

    int i;

    Statistics stGenS, stGenF;
    int usedGen;

    int failNum = 0;

    for (i = 0; i < repeat; i++) {

        GA ga (ell, nInitial, selectionPressure, pc, pm, maxGen, maxFe, subtask);

        if (repeat == 1)
            usedGen = ga.doIt (true, nInitial, i, subtask);
        else
            usedGen = ga.doIt (false, nInitial, i, subtask);

        Chromosome ch(ell);
        // std::cout<<"ga: "<<ga.stFitness.getMax()<<std::endl;
        // std::cout<<"ch: "<<ch.getMaxFitness()<<std::endl;
 
        if (ga.stFitness.getMax() < ch.getMaxFitness()) {
            printf ("-");
            failNum++;
            stGenF.record (usedGen);
        }
        else {
            // std::cout<<"ga.stFitness.getMax(): "<<ga.stFitness.getMax()<<std::endl;
            // std::cout<<"ch.getMaxFitness(): "<<ch.getMaxFitness()<<std::endl;
            printf ("+");
            stGenS.record (usedGen);
        }

        fflush (NULL);

    }

    int fileCount = repeat; // 假设有3个文件
    std::vector<std::vector<int>> sum; // 用于存储加和结果的二维向量

    for (int i = 0; i < fileCount; ++i) {
        std::ostringstream fileName;
        fileName << "./SGA_result/syn_" << subtask << "/a_" << i << "_weak.txt"; // 构造文件名

        std::ifstream inFile(fileName.str());
        if (!inFile.is_open()) {
            std::cerr << "Unable to open file: " << fileName.str() << std::endl;
            continue;
        }

        std::string line;
        int lineCount = 0;
        while (getline(inFile, line)) {
            std::istringstream lineStream(line);
            int value;
            int columnCount = 0;

            // 扩展sum向量以容纳更多行
            if (lineCount >= sum.size()) {
                sum.push_back(std::vector<int>());
            }

            while (lineStream >> value) {
                // 扩展sum向量的内部向量以容纳更多列
                if (columnCount >= sum[lineCount].size()) {
                    sum[lineCount].push_back(0);
                }

                sum[lineCount][columnCount] += value; // 将值加到sum向量对应位置
                ++columnCount;
            }
            ++lineCount;
        }

        inFile.close();
    }

    // 将加和结果写入新文件
    std::ostringstream filePath; // 创建一个ostringstream对象
    filePath << "./SGA_result/syn_" << subtask << "/a_total_weak.txt"; // 使用<<操作符构建文件路径

std::ofstream outFile(filePath.str()); // 将构造好的字符串转换为std::string并传递给ofstream
    if (!outFile.is_open()) {
        std::cerr << "Unable to open file for writing." << std::endl;
        return -1;
    }

    for (const auto &row : sum) {
        for (const auto &elem : row) {
            outFile << elem << " ";
        }
        outFile << std::endl; // 每行结束后换行
    }

    outFile.close();    

    printf ("\nAverage Gen of Success: %f\n", stGenS.getMean());
    printf ("Average Gen of Failure: %f\n", stGenF.getMean());
    printf ("Total number of Failure: %d\n", failNum);

    return EXIT_SUCCESS;
}
