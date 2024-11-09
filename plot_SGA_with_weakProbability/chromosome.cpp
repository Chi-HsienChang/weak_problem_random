/***************************************************************************
 *   Copyright (C) 2004 by Tian-Li Yu                                      *
 *   tianliyu@cc.ee.ntu.edu.tw                                             *
 ***************************************************************************/

#include <stdio.h>
#include "global.h"
#include "chromosome.h"
#include <vector>
#include <cmath>
#include <random>
#include <iostream>


Chromosome::Chromosome ()
{
    length = 0;
    gene = NULL;
    evaluated = false;
}


Chromosome::Chromosome (int n_length)
{
    gene = NULL;
    init (n_length);
}

Chromosome::Chromosome(const Chromosome& other) 
{
    length = other.length;
    evaluated = other.evaluated;
    fitness = other.fitness;
    if (other.gene) {
        gene = new bool[length];
        for (int i = 0; i < length; i++) {
            gene[i] = other.gene[i];
        }
    } else {
        gene = nullptr;
    }
}


Chromosome::~Chromosome ()
{
    delete[]gene;
}


void Chromosome::init (int n_length)
{
    length = n_length;

    if (gene != NULL)
        delete[]gene;

    gene = new bool[length];
    evaluated = false;
}

int Chromosome::getVal (int index) const
{
    if (index < 0 || index > length)
        outputErrMsg ("Index overrange in Chromosome::operator[]");

    return (gene[index])? 1:0;
}


void Chromosome::setVal (int index, int val)
{
    if (index < 0 || index > length)
        outputErrMsg ("Index overrange in Chromosome::operator[]");

    gene[index] = (val==1)? true:false;
    evaluated = false;
}


double Chromosome::getFitness ()
{
    if (evaluated)
        return fitness;
    else
        return (fitness = evaluate ());
}


bool Chromosome::isEvaluated () const
{
    return evaluated;
}

// TODO: random assign fitness
double Chromosome::maxFitnessValue = 0.0;

unsigned long Chromosome::getGenotypeIndex(){

    unsigned long genotypeIndex = 0;
    for (int i = 0; i < length; i++) {
        if (gene[i]) { // 如果gene[i]為true，則將對應的位設置為1
            genotypeIndex |= 1UL << (length - 1 - i);
        }
    }

    return genotypeIndex;
}

void Chromosome::printGene(){

    // std::cout<<"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT: \n";
    for (int i = 0; i < length; i++) {
        std::cout<<gene[i];
    }
    // std::cout<<"\n";
    // std::cout<<"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT: \n";
    // for (int i = 0; i < length; i++) {
    //     std::cout<<getVal(i);
    // }

    // std::cout<<std::endl;
}

std::vector<double> Chromosome::fitnessValues;

double Chromosome::evaluate() {
    if (!evaluated) {
        // 假设你有某种方式将基因数组转换为对应的基因型索引

        unsigned long genotypeIndex = 0;
        for (int i = 0; i < length; i++) {
            if (gene[i]) { // 如果gene[i]為true，則將對應的位設置為1
                genotypeIndex |= 1UL << (length - 1 - i);
            }
        }
        // 从预先生成的适应度值中检索适应度
        fitness = fitnessValues[genotypeIndex];
        evaluated = true;
    }
    return fitness;
}

// double Chromosome::evaluate ()
// {
//     evaluated = true;
//     return oneMax ();
// }


// OneMax
double Chromosome::oneMax () const
{
    int i;
    double result = 0;

    for (i = 0; i < length; i++)
        result += gene[i];

    return result;
}


Chromosome & Chromosome::operator= (const Chromosome & c)
{
    int i;

    if (length != c.length) {
        length = c.length;
        delete[]gene;
        init (length);
    }

    evaluated = c.evaluated;
    fitness = c.fitness;

    for (i = 0; i < length; i++)
        gene[i] = c.gene[i];

    return *this;
}


void Chromosome::printf () const
{
    int i;
    for (i = 0; i < length; i++)
        ::printf ("%d", gene[i]);
}


int Chromosome::getLength () const
{
    return length;
}


double Chromosome::getMaxFitness () const
{
    // For OneMax
    return ((double)maxFitnessValue-1e-6);
}


void Chromosome::set_fitness(double f){
    evaluated = true;
    fitness = f;
}




