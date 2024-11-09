# weak_problem_random

## No weak? 後來發現有!
Consider the cycTrap, BB size = 4, overlapping is 1.

1111 --> 6

Others remains unchanged.

0000 -> 3

1000 -> 2

1100 -> 1

1110 -> 0


labadmin@TeilCluster:~/salima/weak_problem_random/plot_populationSize_weakProbability_cycTrap/inclosure$ ./inclosure 9 1 0 1 -1 9 > debug.txt


在chromosome.cpp中

// cycTrap
double Chromosome::oneMax () const {
    double fHigh = 6;
    double fLow = 3;
    int i, j;
    int u;
    int TRAP_M = length / (TRAP_K-1);
    if (length % (TRAP_K-1) != 0)
        outputErrMsg ("TRAP_k doesn't divide length for Cyclic Setting");
    double result = 0;
   
    for (i = 0; i < TRAP_M; i++) {
        u = 0;
        int idx = i * TRAP_K - i;
        for (j = 0; j < TRAP_K; j++) {
            int pos = idx + j;
            if (pos == length)
                pos = 0;
            else if (pos > length)
                outputErrMsg ("CYCLIC BUG");
            //
            u += getVal(pos);
        }
        result += trap (u, fHigh, fLow, TRAP_K);
    }

    

    return result;
}


double Chromosome::trap (int unitary, double fHigh, double fLow, int trapK) const {
    // if (unitary > trapK)
    //     return 0;

    // cout << "resultOOOOOO: " << endl;

    if (unitary == 4)
        return 6;
    if (unitary == 0)
        return 3;
    if (unitary == 1)
        return 2;
    if (unitary == 2)
        return 1;
    if (unitary == 3)
        return 0;
    else
        return 0;


    // else
    //     return fLow - unitary * fLow / (trapK-1);
}

以及在 inclosure.cpp 中:


    int test_problem_debug = 1;

    if (test_problem_debug){


        //trap debug
        // ############################
        unsigned long enum_size = std::pow(2, ell);

        for(unsigned long i=0; i< std::pow(2, ell); i++) { 
            // double fit = dis(gen);
            // double fit = buffer[i].evaluate_inclosure2(2);
            double fit = buffer[i].evaluate_inclosure2(0);
            // cout << "fit: " << fit << endl;
            // cout << "buffer: " << buffer[i].getGene(0) << endl;

            buffer[i].set_fitness(fit);
            fitness[i] = fit;
            count_chrom[buffer[i].getGene(0)]++;
            seed++;
        }  

        std::sort(buffer.begin(), buffer.end(), Inclosure::customSortfitness); 
        dumpbuffer();
        // ############################
    }else{

        std::random_shuffle(buffer.begin(), buffer.end());


        // sample t chromosome
        int i = 0;
        for (auto& c: buffer)
        {
            if(i == t) break;
            count_chrom[c.getGene(0)]++; // count_chrom 紀錄被sample到的chromosome
            i++;
        }

    }
