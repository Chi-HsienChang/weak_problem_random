
#pragma once

#include <list>
#include "chromosome.h"

using namespace std;

class Chromosome;

class MSO {
public:
    MSO();

    MSO(int _ell);

    void init(int _ell);

    void updateBest(const Chromosome& c);

    void update();

	void clearGraph();

    void printMSO() const {
        for (int i=0; i<ell; ++i) {
            cout << i << ": [ ";
            for (int j=0; j<ell; ++j)
                if (mso[i][j])
                    cout << j << " ";
            cout << "]" << endl;
        }
    }
    void printGraph() const {
        for (int i=0; i<ell; ++i) {
            for (int j=0; j<ell; ++j)
                cout << graph[i][j] << " ";
            cout << endl;
        }
    }

	void findClique(int, list<int>&);

    void pasteMSO();
    void flipMSO();
    void flipMSO(const Chromosome&);

    bool isReady() const;

    bool isDirty() const {
        return dirty;
    }

    void printInfo() const;

    ~MSO();

    Chromosome *best;

    bool **mso;
    Chromosome **good;

protected:
    void updateGraph();
	void updateGraph(const Chromosome&);
    void updateMSO();
    void updateMSO(int);
    void updateMSO(int, int, int);

    void genOrderELL();

    int ell;
    int uncertain;
    bool dirty;
    int **graph;
    int *orderELL;

};
