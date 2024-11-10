/***************************************************************************
 *   Copyright (C) 2004 by Tian-Li Yu                                      *
 *   tianliyu@cc.ee.ntu.edu.tw                                             *
 ***************************************************************************/

#ifndef GA_H
#define GA_H

#include "chromosome.h"
#include "statistics.h"
#include "myrand.h"
#include <iostream>
#include <vector>

class GA
{

    public:

        int ell;                 // chromosome length
        int nInitial;            // initial population size
        int nCurrent;            // current population size
        int nNextGeneration;     // population size for the next generation
        int selectionPressure;

        double pc;               // prob of XO
        double pm;               // prob of Mutation

        std::vector<Chromosome> buffer;
        Chromosome *population;
        Chromosome *detectPopulation;
        Chromosome *offspring;
        int *selectionIndex;
        int maxGen;
        int maxFe;
        int repeat;
        int fe;
        int generation;
        int bestIndex;
        bool first = true;

        int subtask;

        GA ();
        GA (int n_ell, int n_nInitial, int n_selectionPressure, double n_pc,
            double n_pm, int n_maxGen, int n_maxFe, int subtask);

        ~GA ();

        void init (int n_ell, int n_nInitial, int n_selectionPressure, double n_pc,
            double n_pm, int n_maxGen, int n_maxFe, int subtask);

        void initializePopulation ();
        void evaluate ();

        void selection ();

        /** tournament selection without replacement */
        void tournamentSelection ();

	/** Roulette wheel selection */
	void rwSelection ();

        void crossover ();
        void pairwiseXO (const Chromosome &, const Chromosome &, Chromosome &, Chromosome &);
	void onePointXO (const Chromosome &, const Chromosome &, Chromosome &, Chromosome &);
        void uniformXO (const Chromosome &, const Chromosome &, Chromosome &, Chromosome &, double);

        void mutation ();
        void simpleMutation ();
	void mutationClock ();

        void replacePopulation ();

        void showStatistics ();
        void oneRun (bool output, int n_nInitial, int repeatNumber, int subtask);
        int doIt (bool output, int n_nInitial, int i, int subtask);

        bool shouldTerminate ();
        int getNextPopulation ();

        Statistics stFitness;

        int buildDetectPopulation ();
        bool areChromosomesEqual(Chromosome &a, Chromosome &b) ;
        std::vector<int> detectWeak(int new_size);
        bool isSubset(const std::vector<int>& subset, const std::vector<int>& set);
        std::vector<std::vector<int>> generateBinarySequences(int n);
        std::vector<std::vector<int>> generateCombinations(int n, int k);
        std::vector<std::vector<int>> remove_the_same_with_best(const std::vector<std::vector<int>>& enumerations, const std::vector<int>& best_sequence);
        static bool customSortfitness (Chromosome &a, Chromosome &b);
        void dumpbuffer();
        // std::vector<std::vector<int>> generateCombinations(int n, int k);
        // static bool compareByFitness(Chromosome& a, Chromosome& b);
        

    protected:

        // int ell;                 // chromosome length
        // int nInitial;            // initial population size
        // int nCurrent;            // current population size
        // int nNextGeneration;     // population size for the next generation
        // int selectionPressure;

        // double pc;               // prob of XO
        // double pm;               // prob of Mutation
        // Chromosome *population;
        // Chromosome *detectPopulation;
        // Chromosome *offspring;
        // int *selectionIndex;
        // int maxGen;
        // int maxFe;
        // int repeat;
        // int fe;
        // int generation;
        // int bestIndex;

};
#endif
