/***************************************************************************
 *   Copyright (C) 2004 by Tian-Li Yu                                      *
 *   tianliyu@cc.ee.ntu.edu.tw                                             *
 ***************************************************************************/

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
#include <fstream> // 引入文件I/O库
#include <sstream> // 引入用于字符串流的库

GA::GA ()
{
    ell = 0;
    nInitial = 0;
    nCurrent = 0;
    fe = 0;

    nNextGeneration = 0;
    maxGen = -1;
    maxFe = -1;

    subtask = -1;
    
    
    detectPopulation = NULL;
    offspring = NULL;
    population = NULL;
    offspring = NULL;
    selectionIndex = NULL;

}


GA::GA (int n_ell, int n_nInitial, int n_selectionPressure, double n_pc, double n_pm, int n_maxGen, int n_maxFe, int subtask)
{
    init (n_ell, n_nInitial, n_selectionPressure, n_pc, n_pm, n_maxGen, n_maxFe, subtask);
}


GA::~GA ()
{
    delete[]detectPopulation;
    delete[]population;
    delete[]offspring;
    delete[]selectionIndex;
}


void
GA::init (int n_ell, int n_nInitial, int n_selectionPressure, double n_pc,
double n_pm, int n_maxGen, int n_maxFe, int subtask)
{
    int i;

    ell = n_ell;
    nInitial = n_nInitial;
    nCurrent = nInitial;
    selectionPressure = n_selectionPressure;
    pc = n_pc;
    pm = n_pm;
    maxGen = n_maxGen;
    maxFe = n_maxFe;
    fe = 0;
    subtask = subtask;

    detectPopulation  = new Chromosome[nInitial];
    population = new Chromosome[nInitial];
    offspring = new Chromosome[nInitial];
    selectionIndex = new int[nInitial];

    for (i = 0; i < nInitial; i++) {
        population[i].init (ell);
        offspring[i].init (ell);
    }

    initializePopulation ();
}


void GA::initializePopulation ()
{
    int i, j;
    double p = 0.5;

    for (i = 0; i < nInitial; i++)
        for (j = 0; j < ell; j++)
            if (myRand.uniform () > p)
                population[i].setVal (j, 1);
            else
                population[i].setVal (j, 0);

}

// For now, assuming fixed population size
int GA::getNextPopulation ()
{
    return nCurrent;
}

void GA::selection ()
{
    //rwSelection ();
    tournamentSelection ();
}

// Roulette wheel selection
// This is a O(n^2) implementation
// You can achieve O(nlogn) by using binary search
void GA::rwSelection ()
{
    int i, j;

    // Adjusting population size 
    nNextGeneration = getNextPopulation ();

    if (nNextGeneration != nCurrent) {
        delete[]selectionIndex;
        delete[]offspring;
        selectionIndex = new int[nNextGeneration];
        offspring = new Chromosome[nNextGeneration];

        for (i = 0; i < nNextGeneration; i++)
            offspring[i].init (ell);
    }

    double totalFitness = 0.0;
    for (i = 0; i < nCurrent; i++) 
	totalFitness += population[i].getFitness();

    for (i = 0; i < nNextGeneration; i++) {
	double pointer = totalFitness * myRand.uniform();
	int index = -1;
	double partialSum = 0.0;
	for (j = 0; j < nCurrent; j++) {
	    partialSum += population[j].getFitness();
            if (partialSum >= pointer) {
                index = j;
                break;
            }
	}
	if (index == -1) index = nCurrent - 1;

	selectionIndex[i] = index;
    }

}

// tournamentSelection without replacement
void GA::tournamentSelection ()
{
    int i, j;

    // Adjusting population size 
    nNextGeneration = getNextPopulation ();

    if (nNextGeneration != nCurrent) {
        delete[]selectionIndex;
        delete[]offspring;
        selectionIndex = new int[nNextGeneration];
        offspring = new Chromosome[nNextGeneration];

        for (i = 0; i < nNextGeneration; i++)
            offspring[i].init (ell);
    }

    int randArray[selectionPressure * nNextGeneration];

    int q = (selectionPressure * nNextGeneration) / nCurrent;
    int r = (selectionPressure * nNextGeneration) % nCurrent;

    for (i = 0; i < q; i++)
        myRand.uniformArray (randArray + (i * nCurrent), nCurrent, 0, nCurrent - 1);

    myRand.uniformArray (randArray + (q * nCurrent), r, 0, nCurrent - 1);

    for (i = 0; i < nNextGeneration; i++) {

        int winner = 0;
        double winnerFitness = -DBL_MAX;

        for (j = 0; j < selectionPressure; j++) {
            int challenger = randArray[selectionPressure * i + j];
            double challengerFitness = population[challenger].getFitness ();

            if (challengerFitness > winnerFitness) {
                winner = challenger;
                winnerFitness = challengerFitness;
            }

        }
        selectionIndex[i] = winner;
    }
}


void GA::crossover ()
{
    int i;

    if ((nNextGeneration & 0x1) == 0) { 
    	// nNextGeneration is even
    	
        for (i = 0; i < nNextGeneration; i += 2)
            pairwiseXO (population[selectionIndex[i]], population[selectionIndex[i + 1]],
                offspring[i], offspring[i + 1]);

    }
    else {
        for (i = 0; i < nNextGeneration - 1; i += 2) {
            pairwiseXO (population[selectionIndex[i]], population[selectionIndex[i + 1]],
                offspring[i], offspring[i + 1]);
        }
        offspring[nNextGeneration - 1] =
            population[selectionIndex[nNextGeneration - 1]];
    }

}


void GA::pairwiseXO (const Chromosome & p1, const Chromosome & p2, Chromosome & c1, Chromosome & c2)
{
    if (myRand.uniform () < pc) {
	onePointXO (p1, p2, c1, c2);
//      uniformXO (p1, p2, c1, c2, 0.5);
    }
    else {
        c1 = p1;
        c2 = p2;
    }
}

void GA::onePointXO (const Chromosome & p1, const Chromosome & p2, Chromosome & c1, Chromosome & c2)
{
    int i;
    int crossSite = myRand.uniformInt(1, ell-1);

    for (i = 0; i < crossSite; i++) {
            c1.setVal (i, p1.getVal(i));
            c2.setVal (i, p2.getVal(i));
    }

    for (i = crossSite; i < ell; i++) {
            c1.setVal (i, p2.getVal(i));
            c2.setVal (i, p1.getVal(i));
    }
}

void GA::uniformXO (const Chromosome & p1, const Chromosome & p2, Chromosome & c1, Chromosome & c2, double prob)
{
    int i;

    for (i = 0; i < ell; i++) {
        if (myRand.flip (prob)) {
            c1.setVal (i, p1.getVal(i));
            c2.setVal (i, p2.getVal(i));
        }
        else {
            c1.setVal (i, p2.getVal(i));
            c2.setVal (i, p1.getVal(i));
        }
    }
}

void GA::mutation ()
{
    //simpleMutation ();
    mutationClock ();
}


void GA::simpleMutation ()
{
    int i, j;

    for (i = 0; i < nNextGeneration; i++)
        for (j = 0; j< ell; j++)
            if (myRand.flip(pm)) {
                int val = offspring[i].getVal(j);
                offspring[i].setVal(j, 1-val);
            }
}

void GA::mutationClock ()
{
    if (pm <= 1e-6) return; // can't deal with too small pm

    int pointer = (int) (log(1-myRand.uniform()) / log(1-pm) + 1);

    while (pointer < nNextGeneration * ell) {

	int q = pointer / ell;
	int r = pointer % ell;

        int val = offspring[q].getVal(r);
        offspring[q].setVal(r, 1-val);

	// Compute next mutation clock
	pointer += (int) (log(1-myRand.uniform()) / log(1-pm) + 1);
    };
}


void GA::showStatistics ()
{

    printf ("Gen:%d  Fitness:(Max/Mean/Min):%f/%f/%f Chromsome Length:%d\n",
        generation, stFitness.getMax (), stFitness.getMean (),
        stFitness.getMin (), population[0].getLength ());
    printf ("best chromosome:");
    population[bestIndex].printf ();
    printf ("\n");
}


void GA::replacePopulation ()
{
    int i;

    if (nNextGeneration != nCurrent) {
        delete[]population;
        population = new Chromosome[nNextGeneration];
    }

    for (i = 0; i < nNextGeneration; i++)
        population[i] = offspring[i];

    nCurrent = nNextGeneration;
}


void GA::oneRun (bool output, int n_nInitial, int repeatNumber, int subtask)
{
    if (first)
    {
        // std::cout<<"------- Enumeration -------"<<std::endl;
        Chromosome::initializeFitnessValues(ell);   

        int new_size = buildDetectPopulation();
        std::vector<int> one_weak = detectWeak(new_size);
        
        first = false;

        // 打开一个文件用于写入，如果文件不存在则创建它
        // std::ofstream outFile("a_weak.txt");
        std::ostringstream fileName; // 创建一个字符串流
        fileName << "./SGA_result/syn_" << subtask << "/a_" << repeatNumber << "_weak.txt"; // 将文本和变量i拼接成文件名
        std::ofstream outFile(fileName.str()); // 使用字符串流的str()方法获取构造的字符串，并打开文件

        // 检查文件是否成功打开
        if (outFile.is_open()) {
            // 遍历one_weak向量并将每个值写入文件

            for(int i=1; i<ell; i++){
                outFile << one_weak[i] << " "; // 将值写入文件，每个值后换行
            }
            outFile <<std::endl;
            outFile.close(); // 完成写入后关闭文件
        } else {
            std::cerr << "Unable to open file for writing." << std::endl; // 如果文件无法打开，输出错误信息
        }

        

    }


    

    // std::cout<<"------- Population (Before evolution) -------"<<std::endl;
    // for (int i = 0; i < nCurrent; ++i) {
    //     population[i].printGene();
    //     std::cout<<": "<<population[i].getFitness()<<std::endl;
    // }
    

    

    int i;

    selection ();
    crossover ();
    mutation ();
    replacePopulation ();

    // std::cout<<"------- After evolution -------"<<std::endl;

    double max = -DBL_MAX;
    stFitness.reset ();
    for (i = 0; i < n_nInitial; i++) {
        double fitness = population[i].getFitness ();
        // population[i].printGene();
        // std::cout<<" : "<<fitness<<" "<<std::endl;
        if (fitness > max) {
            max = fitness;
            bestIndex = i;
        }
        stFitness.record (fitness);
    }



      

    int new_size = buildDetectPopulation();
    std::vector<int> one_weak = detectWeak(new_size);
    // dumpbuffer();

    // 打开一个文件用于写入，如果文件不存在则创建它
    // std::ofstream outFile("a_weak.txt", std::ios::app);

    std::ostringstream fileName; // 创建一个字符串流
    // std::cout<<"subtask: "<<subtask<<std::endl;
    fileName << "./SGA_result/syn_" << subtask << "/a_" << repeatNumber << "_weak.txt"; // 将文本和变量i拼接成文件名
    std::ofstream outFile(fileName.str(), std::ios::app); // 使用字符串流的str()方法获取构造的字符串，并打开文件

    // 检查文件是否成功打开
    if (outFile.is_open()) {
        // 遍历one_weak向量并将每个值写入文件

        for(int i=1; i<ell; i++){
            outFile << one_weak[i] << " "; // 将值写入文件，每个值后换行
        }
        outFile <<std::endl;
        outFile.close(); // 完成写入后关闭文件
    } else {
        std::cerr << "Unable to open file for writing." << std::endl; // 如果文件无法打开，输出错误信息
    }

    // first += 1;

    if (output)
        showStatistics ();

    generation++;
}


int GA::doIt (bool output, int n_nInitial, int i, int subtask)
{
    generation = 0;

    while (!shouldTerminate ()) {
        oneRun (output, n_nInitial, i, subtask);
    }
    return generation;
}


bool GA::shouldTerminate ()
{

    // std::cout<<"shouldTerminate\n";
    bool termination = false;

    //  std::cout<<"maxFe: "<<maxFe<<std::endl;
    //  std::cout<<"fe: "<<fe<<std::endl;

    // Reach maximal # of function evaluations
    if (maxFe != -1) {
        if (fe > maxFe){
            // std::cout<<"termination: "<<termination<<std::endl;
            termination = true;
        }
    }

    // Reach maximal # of generations
    if (maxGen != -1) {
        if (generation >= maxGen)
            termination = true;
    }

    // Found a satisfactory solution
    //if (stFitness.getMax() >= population[0].getMaxFitness())
    //    termination = true;

    // The population loses diversity
    // if (stFitness.getMax()-1e-6 < stFitness.getMean())
	// termination = true;

    return termination;

}


int GA::buildDetectPopulation ()
{

    // for (int i = 0; i < nCurrent; ++i) {
    //     population[i].printGene();
    //     std::cout<<": "<<population[i].getFitness()<<std::endl;
    // }

    // std::cout<<"########## buildDetectPopulation ##########"<<std::endl;
    int detectSize = 0; // 初始化detectPopulation的大小
    
    // std::cout<<"nCurrent: "<<nCurrent<<std::endl;
    
    for (int i = 0; i < nCurrent; ++i) {


        // double fitness = population[i].getFitness ();

        // population[i].printGene();
        // std::cout<<": "<<population[i].getFitness()<<std::endl;
        // std::cout<<"population["<<i<<"]: "<<fitness<<std::endl;

        bool isUnique = true;
        for (int j = 0; j < detectSize; ++j) {
            // if (i == 0){
            //     break;
            // }
            if (areChromosomesEqual(population[i], detectPopulation[j])) {
                isUnique = false;
                break;
            }
        }

        if (isUnique) {
            // std::cout<<"isUnique success: ";

            Chromosome c(ell);
            // c.init (ell);
            for (int j = 0; j < ell; ++j) {
                c.setVal(j,  population[i].getVal(j));
            }

            // std::cout<<"population["<<i<<"]: ";
            
            // population[i].printGene();
            // std::cout<<": "<<population[i].getFitness()<<std::endl;

            // std::cout<<"Chromosome c(ell): ";
            // c.printGene();
            // std::cout<<std::endl;
            // std::cout<<"getGenotypeIndex: "<<c.getGenotypeIndex()<<std::endl;

            detectPopulation[detectSize] = c; // 存儲相異的染色體


            // double detectFitness = detectPopulation[detectSize].getFitness ();
            // std::cout<<"detectPopulation fitness: "<<detectFitness<<std::endl;

            // for (int i = 0; i < population[i].getLength(); ++i) {
            //     detectPopulation[detectSize].setVal(i,  population[i].getVal(i));
            // }

            detectSize++;
        }
    }  


    // std::cout<<"########## buildDetectPopulation END ##########"<<std::endl;

    return detectSize;  
}

bool GA::areChromosomesEqual(Chromosome &a, Chromosome &b) {

    for (int i = 0; i < a.getLength(); ++i) {
        if (a.getVal(i) != b.getVal(i)) return false; // 逐位比較基因序列
    }

    return true; // 所有基因位都相同，則兩個染色體相等
}


std::vector<int> GA::detectWeak(int new_size){
    int debug = 0;
    if(debug) std::cout<<"-------------------- detectWeak start --------------------"<<std::endl;
    

    std::vector<std::vector<std::vector<int>>> weak_epi(ell);
    std::vector<int> one_weak(ell, 0); 
    std::vector<double> fitness(new_size, 0); 

    buffer.clear();

    if(debug) std::cout<<"new_size: "<<new_size<<std::endl;



    for (int i = 0; i < new_size; ++i) {
        Chromosome c(detectPopulation[i]);
        buffer.push_back(c); 
    }

    for(int i=0; i < new_size; i++){
         double fit =  buffer[i].getFitness();
         buffer[i].set_fitness(fit);
         fitness[i] = fit;
    }   

    
    std::sort(buffer.begin(), buffer.end(), GA::customSortfitness);

    if(debug){
        for (int i = 0; i < new_size; ++i) {
            buffer[i].printGene();
            std::cout<<": "<<buffer[i].getFitness()<<std::endl;
        }
    }


    // dumpbuffer();
    

    bool **graph = new bool*[ell];
    for (int i=0; i<ell; i++) graph[i]=new bool[ell];

    for (int i=0; i<ell; i++){
        for (int j=0; j<ell; j++){
            graph[i][j] = false;
        }
    }

    
    for (int u=1; u<ell; u++) 
    {
            
        int val = 1 - buffer[0].getVal(u); // flip bit

        int w = 1;

        // if (debug) std::cout<<"w index: "<<w<<std::endl;

        if(w > new_size-1) break;

        while (buffer[w].getVal(u) != val){
            // if (debug) std::cout<<"while: "<<w<<std::endl;
            w++;
            if(w > new_size-1) break;
        }

        // if (debug) std::cout<<"w index: "<<w<<std::endl;
        // if (debug) std::cout<<"======="<<std::endl;
        if(w > new_size-1) continue;

        
        // if (debug) std::cout<<"buffer[w].getVal(0): "<<detectPopulation[w].getVal(0)<<std::endl;
        // if (debug) std::cout<<"buffer[b].getVal(0): "<<detectPopulation[0].getVal(0)<<std::endl;
        if (buffer[w].getVal(0) != buffer[0].getVal(0))
        { // if change bit in position u, change bit in position v too
            if(debug) printf("%u -> 0\n", u);
            std::vector<int> u_vector(1, u);
            weak_epi[1].push_back(u_vector);
            one_weak[1]++;
            graph[u][0] = true;
        }
        
    }
    

    // 遍歷所有組合可能
    for (int epi_size = 2; epi_size < ell; epi_size++)
    {  
        auto combinations = generateCombinations(ell-1, epi_size); // combinations = { [1, 2], [1, 3], [2, 3] }
        
        for (auto& combination : combinations) // combination = [1, 2]
        { 

            bool not_find_smaller_epi = true;

            int smaller_epi_size = epi_size;
            bool is_subset;
            while(not_find_smaller_epi && smaller_epi_size>=2)
            {
                for (auto& previous : weak_epi[smaller_epi_size-1]) 
                {
                    is_subset = isSubset(previous, combination);
                    not_find_smaller_epi = !is_subset;
                    if(is_subset) break;
                }
                smaller_epi_size--;
            }

            
            // std::cout<<"epi_size "<<epi_size<<std::endl;
            // std::cout<<"not_find_smaller_epi "<<not_find_smaller_epi<<std::endl;
            // 如果沒有找到更小的epi，則進入if >>> 這裡使用 weak 機制排除
            if(not_find_smaller_epi)
            {
                // std::cout<<"enter if "<<std::endl;
                auto enumerations = generateBinarySequences(epi_size); // enumerations = { [0, 0], [0, 1], [1, 0], [1, 1] }    


                // 遍歷所有的 pattern 不用 remove best_sequence 的 pattern
                for (auto& enumeration : enumerations) // combination = [1, 2] and enumeration_without_best = { [0, 0], [0, 1], [1, 0]}
                { 
                    auto enumeration_original = enumeration;

                    // b is 主角pattern 出現 的 index
                    int b_original = 0;
                    int b = b_original;

                    int hit_original = 0; // 這個是用來計算有幾個相同的值
                    while(1)
                    {  
                        if(b > (new_size-1)) break;

                        hit_original = 0;
                        // combination 的角色是 index 的位置
                        for(int i = 0; i < enumeration.size(); i++){
                            if(buffer[b].getVal(combination[i]) == enumeration[i]) hit_original++; 
                        }
                        
                        if(hit_original == epi_size) break;
                        if(hit_original < epi_size) b++;    // b 不動 chrom_index chrom_index 往後移動
                    }
                    
                    if(b > (new_size-1)) break;

                    // 以當下 enumeration 為主角，去翻每一個 bit >>> 主角的 omega[0] != 去翻每一個omega[0] 這樣就有 epi >>> 就可以 break
                    //// int chrom_index = b+1; // 從fitness第二高的chromosome開始找


                    // chrom_index is 主角pattern flip 一個bit後 omega 的 index
                    int chrom_index = b_original; // 從fitness第二高的chromosome開始找

                    if(chrom_index > (new_size-1)) break;     

                    // while(buffer[chrom_index].getVal(0) <= 0){
                    //     chrom_index++;
                    //     if(chrom_index > (new_size-1)) break;
                    // }

                    // if(debug) printf("############## DEBUG ##############\n");
                    
                    int condition_holds = 0;
                    for (int condition_index = 0; condition_index < enumeration.size(); condition_index++)
                    {
                        enumeration = enumeration_original;
                        enumeration[condition_index] = 1 - enumeration[condition_index]; // flip bit 
                        chrom_index = b_original; // 從頭開始找


                        int hit = 0; // 這個是用來計算有幾個相同的值
                        while(1)
                        {  
                            if(chrom_index > (new_size-1)) break;

                            hit = 0;
                            // combination 的角色是 index 的位置
                            for(int i = 0; i < enumeration.size(); i++){
                                if(buffer[chrom_index].getVal(combination[i]) == enumeration[i]) hit++; 
                            }
                            
                            if(hit == epi_size) break;
                            if(hit < epi_size) chrom_index++;    // b 不動 chrom_index chrom_index 往後移動
                        }
                        
                        if(chrom_index > (new_size-1)) break;

                        if((buffer[chrom_index].getVal(0) != buffer[b].getVal(0)) && buffer[chrom_index].getVal(0) > 0)
                        { 
                            condition_holds += 1;
                        }else{

                            break; // 先不處理fitness相同的情況
                            // Continue searching and put chromosomes with the same fitness as buffer[chrom_index].getVal(0) into a vector
                            std::vector<int> same_fitness_indices;
                            int next_index = chrom_index + 1;
                            
                            // same_fitness_indices.push_back(chrom_index);
                            while (next_index < std::pow(2, ell) && buffer[next_index].getFitness() == buffer[chrom_index].getFitness()) {
                                if (buffer[next_index].getVal(0) > 0) {
                                    same_fitness_indices.push_back(next_index);
                                }
                                next_index++;
                            }

                            int go_to_next_enumeration = 1;
                            for (int idx : same_fitness_indices) {
                                if (buffer[idx].getVal(0) != buffer[b].getVal(0)) {
                                    condition_holds += 1;
                                    go_to_next_enumeration = 0;
                                    break;
                                }
                            }

                            if (go_to_next_enumeration) {
                                break;
                            };
                            
                        }
                    
                    }

                    if(condition_holds == epi_size)
                    {
                        one_weak[epi_size]++;
                        weak_epi[epi_size].push_back(combination);


                        
                        if(debug){

                            std::cout << "condition_holds is: "<< condition_holds << std::endl;
                            


                            std::cout << "!!!the pattern is: ";
                            for(int p = enumeration_original.size()-1; p >= 0; p--){
                                std::cout << enumeration_original[p] << " ";
                            }
                            std::cout << std::endl;
                            // for(auto val : enumeration_original){
                            //     std::cout << val << " ";
                            // }
                            // std::cout << std::endl;


                            std::cout << "!!!the fitenss is: ";
                            std::cout << buffer[b].getFitness() <<std::endl;




                            printf("{ ");
                            for (const auto& vertex : combination)
                            {
                                printf("%u ", vertex);
                            }
                            printf("} -> 0\n");
                        }

                        break;
                    }

                }
            }

        }
    }



    if(debug) std::cout<<"-------------------- detectWeak end haha--------------------"<<std::endl;

    // std::cout<<"END"<<std::endl;
    delete[] graph;

    return one_weak;
}


// bool GA::compareByFitness(Chromosome & a, Chromosome & b) {
//     return a.getFitness() > b.getFitness();
// }


bool GA::isSubset(const std::vector<int>& subset, const std::vector<int>& set) {
    for (const auto& elem : subset) {
        if (std::find(set.begin(), set.end(), elem) == set.end()) {
            // 如果找不到元素，則返回 false
            return false;
        }
    }
    // 所有元素都被找到，返回 true
    return true;
}

std::vector<std::vector<int>> GA::generateBinarySequences(int n) {
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

std::vector<std::vector<int>> GA::generateCombinations(int n, int k) {
    std::vector<int> bitmask(k, 1);  // 創建k個1
    bitmask.resize(n, 0);  // 後面填充n-k個0

    std::vector<std::vector<int>> combinations;

    do {
        std::vector<int> currentCombination;
        for (int i = 0; i < n; ++i) {
            if (bitmask[i]) {
                currentCombination.push_back(i + 1);  // 輸出選中的元素
            }
        }
        combinations.push_back(currentCombination);
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));  // 生成下一個排列

    return combinations;
}


std::vector<std::vector<int>> GA::remove_the_same_with_best(const std::vector<std::vector<int>>& enumerations, const std::vector<int>& best_sequence) {
    std::vector<std::vector<int>> result = enumerations;  

    auto newEnd = std::remove_if(result.begin(), result.end(),
                                 [&best_sequence](const std::vector<int>& sequence) {
                                     return sequence == best_sequence;
                                 });
    result.erase(newEnd, result.end());

    return result;
}

bool GA::customSortfitness (Chromosome &a, Chromosome &b) {
    return a.getFitness() > b.getFitness();
}






    
