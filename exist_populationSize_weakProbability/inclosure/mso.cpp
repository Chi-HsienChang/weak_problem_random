
#include "chromosome.h"
#include "doublelinkedlistarray.h"
#include "mso.h"

MSO::MSO(): graph(nullptr), best(nullptr), good(nullptr), dirty(true), mso(nullptr), orderELL(nullptr) {}

MSO::MSO(int _ell): ell(_ell) {
    init(ell);
}

void MSO::init(int _ell) {

    ell = _ell;
    uncertain = 2*ell;
    dirty = true;

    best = new Chromosome(ell);

    good = new Chromosome*[ell];
    for (int i=0; i<ell; ++i)
        good[i] = new Chromosome[2];

    for (int i=0; i<ell; ++i)
        for (int j=0; j<2; ++j)
            good[i][j].init(ell);

    graph = new int*[ell];
    for (int i=0; i<ell; ++i)
        graph[i] = new int[ell];

    mso = new bool*[ell];
    for (int i=0; i<ell; ++i)
        mso[i] = new bool[ell];

    orderELL = new int[ell];
}

bool MSO::isReady() const {
    return (uncertain == 0);
}

void MSO::update() {
	if (!dirty) return;
    printf("updating");
	clearGraph();
    updateGraph();
    dirty = false;
}

void MSO::updateBest(const Chromosome& c) {

    for (int i=0; i<ell; ++i) {
        int val = c.getVal(i);

        if (!good[i][val].isEvaluated()) {
            good[i][val] = c;
            --uncertain;
        }

        else {
            if (good[i][val].getFitness() < c.getFitness()) {
                good[i][val] = c;
                dirty = true;
            }
        }
        
    }

    if (!best->isEvaluated()) 
        *best = c;
    else {
        if (best->getFitness() < c.getFitness()) {
            *best = c;
            dirty = true;
        }
    }
}

void MSO::updateGraph() {

	// should I do all?
	
	for (auto c:Chromosome::population)
		updateGraph(c);
	
		
	//updateGraph(*best);
}

void MSO::clearGraph() {
	for (int u=0; u<ell; ++u)
		for (int v=0; v<ell; ++v)
			graph[u][v] = 0;
}

void MSO::updateGraph(const Chromosome& c) {

    // u -> v ?
    for (int u=0; u<ell; ++u) {
        for (int v=0; v<ell; ++v) {

            int val = 1 - c.getVal(u); // flip bit
            if (good[u][val].getVal(v) != c.getVal(v)) // if change bit in position u, change bit in position v too
                ++graph[u][v];
        }

    }

}

// recursive
void MSO::updateMSO(int start, int now, int d) {

    /*
    if (start == 0) {
    	for (int i=0; i<4*d; ++i)
    		cout << " ";

    	cout << "-------(" << start << "," << now << ")" << endl;
    }
    */

    for (int i=0; i<ell; ++i) {
        if (graph[i][now] == 1) {
            if (!mso[start][i]) {
                mso[start][i] = true;
                updateMSO(start, i, d+1);
            }

        }
    }

}

void MSO::updateMSO(int start) {

    for (int i=0; i<ell; ++i) {
        if (graph[i][start] == 1) {
            mso[start][i] = true;
            updateMSO(start, i, 0);
        }
    }

}

void MSO::updateMSO() {

    for (int i=0; i<ell; ++i)
        for (int j=0; j<ell; ++j)
            mso[i][j] = false;

    for (int i=0; i<ell; ++i)
        updateMSO(i);

}

void MSO::flipMSO(const Chromosome& c) {

    for (int i=0; i<ell; ++i) {
		list<int> mask;
		findClique(i, mask);

		for (int size = 0; size < ell/2; ++size) {
	    
			Chromosome trial(*best);
			int x = 0;
			for (auto it = mask.begin(); x < size; ++x, ++it) 
				trial.flip(*it);
			trial.getFitness();

		}

	}
}

void MSO::flipMSO() {

	flipMSO(*best);

	for (int i=0; i<ell; ++i)
		flipMSO(good[i][1-best->getVal(i)]);
	/*
	size_t size = Chromosome::population.size();

	for (size_t i=0; i<size; ++i) 
		flipMSO(Chromosome::population[i]);
		*/
}

// paste MSO from best to good
void MSO::pasteMSO() {

    // from 0 to ell-1
    // should be from short to long mso
    for (int i=0; i<ell; ++i) {

        for (int j=0; j<ell; ++j) {
            if (!mso[i][j]) {
                Chromosome trial(ell);
                trial = *best;
                for (int k=0; k<ell; ++k) {
                    if (!mso[i][k])
                        trial.setVal(k, good[j][1-best->getVal(j)].getVal(k));
                }
                trial.getFitness();
            }

        }
    }

}

void MSO::findClique(int startNode, list<int>& result) {
	
    result.clear();

    DLLA rest(ell);
    
    genOrderELL();

    for (int i=0; i<ell; ++i) {
        if (orderELL[i]==startNode)
            result.push_back(orderELL[i]);
        else
            rest.insert(orderELL[i]); 
    }   

    double *connection = new double[ell];

    for (DLLA::iterator iter = rest.begin(); iter != rest.end(); ++iter)
        connection[*iter] = graph[*iter][startNode];

    while (!rest.isEmpty()) {

        double max = -INF;
        int index = -1; 
        for (DLLA::iterator iter = rest.begin(); iter != rest.end(); ++iter) {
            if (max < connection[*iter]) {
                max = connection[*iter];
                index = *iter;
            }
        }

        rest.erase(index);
        result.push_back(index);

        for (DLLA::iterator iter = rest.begin(); iter != rest.end(); ++iter)
            connection[*iter] += graph[*iter][index];
    }


    delete []connection;

}

inline void MSO::genOrderELL() {
    myRand.uniformArray(orderELL, ell, 0, ell-1);
}

void MSO::printInfo() const{
	cout << "Best Fitness: " << best->getFitness() << endl;
	cout << "Best Chromosome: " << endl;
	cout << *(best) << endl;
	cout << "out of " << Chromosome::population.size() << " in total." << endl;
}

MSO::~MSO() {
    if (best != NULL)
        delete best;

    if (good != NULL) {
        for (int i=0; i<ell; ++i)
            delete []good[i];
        delete []good;
    }

    if (graph != NULL) {
        for (int i=0; i<ell; ++i)
            delete []graph[i];
        delete []graph;
    }
    if (mso != NULL) {
        for (int i=0; i<ell; ++i)
            delete []mso[i];
        delete []mso;
    }
    if (orderELL != NULL){
        delete []orderELL;
    }
}
