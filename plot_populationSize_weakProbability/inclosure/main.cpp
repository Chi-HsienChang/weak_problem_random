/***************************************************************************
 *   Copyright (C) 2015 by Tian-Li Yu                                      *
 *   tianliyu@ntu.edu.tw                                                   *
 ***************************************************************************/

#include <math.h>
#include <iostream>
#include <cstdlib>
#include <fstream>

#include "statistics.h"
#include "dsmga2.h"
#include "global.h"
#include "chromosome.h"

using namespace std;

int main (int argc, char *argv[]) {

    bool wantLog = false;
    ofstream outputFile("log.txt");
    streambuf* coutBuffer;
    if (outputFile.is_open()) {
        if(wantLog){
            cout.rdbuf(); 
            cout.rdbuf(outputFile.rdbuf()); 
        }

    if (argc != 9) {
        printf ("MSO ell nInitial function maxGen maxFe repeat display rand_seed\n");
        printf ("function: \n");
        printf ("     ONEMAX:  0\n");
        printf ("     MK    :  1\n");
        printf ("     FTRAP :  2\n");
        printf ("     CYC   :  3\n");
        printf ("     NK    :  4\n");
        printf ("     SPIN  :  5\n");
        printf ("     SAT   :  6\n");

        return -1;
    }

    int ell = atoi (argv[1]); // problem size
    int nInitial = atoi (argv[2]); // initial population size
    int fffff = atoi (argv[3]); // function
    int maxGen = atoi (argv[4]); // max generation
    int maxFe = atoi (argv[5]); // max fe
    int repeat = atoi (argv[6]); // how many time to repeat
    int display = atoi (argv[7]); // display each generation or not
    int rand_seed = atoi (argv[8]);  // rand seed

    if (fffff == 4) {

        char filename[200];
        sprintf(filename, "./NK_Instance/pnk%d_%d_%d_%d", ell, 4, 4, 0);

        printf("Loading: %s\n", filename);
        FILE *fp = fopen(filename, "r");
        loadNKWAProblem(fp, &nkwa);
        fclose(fp);
    }

    if (rand_seed != -1)  // time
        myRand.seed((unsigned long)rand_seed);

    int i;

    Statistics stGen, stFE, stLSFE;
    int usedGen;

    int failNum = 0;

    Chromosome::function = (Chromosome::Function)fffff;
    Chromosome::nfe = 0;
    Chromosome::lsnfe = 0;
    Chromosome::hitnfe = 0;
    Chromosome::hit = false;

    mso.init(ell);

	for (int i=0; i<nInitial; ++i) {
		Chromosome test;
		test.initR(ell);
		test.GHC();    
        if(i%100==0) printf("still in GHC\n");    
	}
    int count = 0;
    while (!mso.isReady()) {
        Chromosome test;
        test.initR(ell);
        test.getFitness();
        if(count++%100==0) printf("not ready\n");    
    }

    mso.printInfo();

    for (int i=0; i<100; ++i) {
        //if(i%10==0) 
            printf("Now in loop %d\n", i*10);    

		if (Chromosome::hit) {
            mso.printInfo();
			break;
		}

		while (mso.isDirty()) {
            cout << "dirty\n";
            mso.update();

            mso.flipMSO();
            mso.printInfo();

			if (Chromosome::hit) 
				break;

        };

		
        cout << "Dump all good\n"; // show the not best ones
		for (int i=0; i<ell; ++i)  
			cout << mso.good[i][1-mso.best->getVal(i)] << endl;

		if (Chromosome::hit) {
            mso.printInfo();
			break;
		}
    }

        cout.rdbuf(coutBuffer);
        outputFile.close();
    } 
    return EXIT_SUCCESS;
}
