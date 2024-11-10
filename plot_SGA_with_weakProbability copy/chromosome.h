/***************************************************************************
 *   Copyright (C) 2004 by Tian-Li Yu                                      *
 *   tianliyu@cc.ee.ntu.edu.tw                                             *
 ***************************************************************************/

#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "global.h"
#include "statistics.h"
#include "myrand.h"
#include "ga.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <random>
#include <bitset>
#include <iomanip>
#include <utility>
#include <memory>
#include <cmath>
#include <cstdlib>

#ifndef _CHROMOSOME_H
#define _CHROMOSOME_H


class Chromosome
{
    public:
        Chromosome ();
        Chromosome (int n_ell);

        ~Chromosome ();

        Chromosome& operator= (const Chromosome & c);
        Chromosome(const Chromosome& other);

        void init (int n_ell);

        int getVal (int index) const;
        void setVal (int index, int val);

        double getFitness ();

        /** real evaluator */
        double evaluate ();

        double oneMax () const;

        bool isEvaluated () const;

        void printf () const;

        int getLength () const;

        double getMaxFitness () const;

        void set_fitness(double f);

        static std::vector<double> fitnessValues;

        static double maxFitnessValue; // 用于存储最大适应度值

        static void initializeFitnessValues(int ell) {
            std::default_random_engine generator;
            std::uniform_real_distribution<double> distribution(0.0, 1.0); // 假设适应度值在0.0到1.0之间

            fitnessValues.resize(std::pow(2, ell));
            maxFitnessValue = 0.0; // 初始化最大适应度值
            int theMaxIndex = -1;

            for(unsigned long i = 0; i < std::pow(2, ell); i++) {
                fitnessValues[i] = distribution(generator); // 为每个基因型随机分配适应度值

               
        //         printBinary(i, ell);
        //         std::cout<<" ";
           
        //         std::cout<<": " << std::setw(5) << std::setprecision(5) 
        // << std::fixed <<fitnessValues[i]<<std::endl;
                if (fitnessValues[i] > maxFitnessValue) {
                    maxFitnessValue = fitnessValues[i]; // 更新最大适应度值
                    theMaxIndex = i;
                }
            }

        //     std::cout<<"The global optimal: ";
        //     printBinary(theMaxIndex, ell);
        //     std::cout<<" ";
        //         std::cout<<": " << std::setw(5) << std::setprecision(5) 
        // << std::fixed <<fitnessValues[theMaxIndex]<<std::endl;

        // std::cout<<"-----------------------------\n";

        }

        static void printBinary(unsigned long num, int ell) {
            const int bits = ell; // 计算变量的位数，每个字节8位
            for (int i = bits - 1; i >= 0; i--) {
                unsigned long bit = (num >> i) & 1UL; // 将num右移i位，然后与1进行按位与操作，以获取第i位的值
                std::cout << bit;
            }
            // std::cout << std::endl; // 换行
        }


        bool *gene;
        unsigned long getGenotypeIndex();
        void printGene();



    protected:
        
        int length;
        double fitness;
        bool evaluated;
        
};
#endif
