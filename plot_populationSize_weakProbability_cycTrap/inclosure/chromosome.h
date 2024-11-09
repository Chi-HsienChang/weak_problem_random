/***************************************************************************
 *   Copyright (C) 2011 by TEIL                                        *
 *                                                                         *
 ***************************************************************************/

#pragma once

#include <cassert>
#include <iostream>
#include <vector>
#include <unordered_map>
#include "nk-wa.h"
#include "global.h"

using namespace std;


class Chromosome {

public:

    static enum Function {
        ONEMAX=0,
        MKTRAP=1,
        FTRAP=2,
        CYCTRAP=3,
        NK=4,
        SPINGLASS=5,
        SAT=6,
        CNIAH=7
    } function;


    Chromosome ();
	Chromosome (const Chromosome&);
    Chromosome (int n_ell);

    ~Chromosome ();

	// not const because we may update fitness
    bool isinP();

    bool hasSeen() const;

    bool GHC();
    void steepestDescent();

    void init (int _ell);
    void init0 (int _ell);
    void initR (int _ell);

    bool tryFlipping (int index);


    int getVal (int index) const;

    void setVal (int index, int val);

    unsigned long getKey () const;

    void setValF (int index, int val);

    void flip (int index);
    void flip_gene0 (int index);

    /** real evaluator */
    void evaluate ();

    bool isEvaluated () const;

    bool operator== (const Chromosome & c) const;
    Chromosome & operator= (const Chromosome & c);

    double getFitness ();
    double getFitness () const;
    double trap (int u, double high, double low, int trapK) const;
    double niah (int u, double high, double low, int trapK) const;
    double oneMax () const;
    double mkTrap (double high, double low) const;
    double cniah (double high, double low) const;
    double cycTrap(double fHigh, double fLow) const;
    double fTrap () const;
    double spinGlass () const;
    double nkFitness() const;
    double satFitness() const;


	friend ostream& operator<< (ostream& o, const Chromosome& c);
    void printOut () const;
    void shortPrintOut () const;

    int getLength () const;
    void setLength ();

    double getMaxFitness () const;
    void setGene(int index, unsigned long val);
    unsigned long getGene(int index) const;
    bool operator<(const Chromosome& other) const;

    // this part is used only in inclosure.cpp
    void evaluate_inclosure(int function);
    double evaluate_inclosure2(int function);
    void set_fitness(double f);

public:
    static int nfe;
    static int lsnfe;
    static int hitnfe;
    static bool hit;
    static unordered_map<unsigned long, double> cache;
    static vector<Chromosome> population;
    int count;

protected: 

    unsigned long *gene;
    int length;
    int lengthLong;
    double fitness;
    bool evaluated;
    unsigned long key;
    

};


