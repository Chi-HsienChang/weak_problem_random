#include "inclosure.h"

int ell;
int mode; // mode 0: for enumeration, mode 1 for random assignment provided given n, mode 2 for random assign fitness, mode 3 for function
unsigned long n;
int function;
int randomseed =2024;
std::vector<std::vector<unsigned long>> epistasis;

std::ofstream outputFile("inclosure.txt");
std::default_random_engine rng(randomseed);

void dumpbuffer(std::vector<Chromosome> buffer);
void findInclosure(std::vector<int> &inclosure, bool *graph[], bool visited[], int index);
void printGraph(bool *graph[]);
bool isSO(std::vector<int> &inclosure, unsigned long rank[], std::vector<Chromosome> &buffer, bool isinP[]);
bool isSO(std::vector<int> &inclosure, double fitness[], std::vector<Chromosome> &buffer, bool isinP[]);
bool isMSO(std::vector<int> &inclosure, unsigned long rank[], std::vector<Chromosome> &buffer, bool isinP[]);
bool isMSO(std::vector<int> &inclosure, double fitness[], std::vector<Chromosome> &buffer, bool isinP[]);
bool customSortfitness (const Chromosome &a, const Chromosome &b);
void GHC_inclosure(std::vector<Chromosome> &buffer, double fitness[]);
void findEpistasis(std::vector<Chromosome> &buffer);
bool isWeak(std::vector<Chromosome> &buffer);
void printEpistasis(int v, int k);

unsigned long long factorial(int n);

int main(int argc, char *argv[]){

    std::vector<Chromosome> buffer;

    if(argc<3) {
        printf("Please check the input format\n:  mode 0: for enumeration, mode 1:random rank given n, mode 2: random assign fitness, mode 3: problem given\n");
        printf ("     ONEMAX:  0\n");
        printf ("     MK    :  1\n");
        printf ("     FTRAP :  2\n");
        printf ("     CYC   :  3\n");
        printf ("     NK    :  4\n");
        printf ("     SPIN  :  5\n");
        printf ("     SAT   :  6\n");
        return -1;
    }

    ell = atoi (argv[1]);
    mode = atoi (argv[2]);
    GHC = atoi(argv[3]);
    if(mode == 1 || mode == 2) n = atoi (argv[4]);
    else if(mode == 3) function = atoi(argv[4]);
    for(int i=0; i<ell; i++) epistasis.push_back(std::vector<unsigned long>());
    if (outputFile.is_open()) outputFile.flush();
    else printf("some thing wrong");

    // build EG for each rank
    bool *graph[ell];
    for (int i=0; i<ell; i++){
        graph[i] = new bool[ell];
    }
    unsigned long maxInDegree = 0;

    int counter = 0;

    for(unsigned long i=0; i < std::pow(2, ell); i++){
        Chromosome c;
        c.init(ell);
        c.setGene(0, i);
        buffer.push_back(c);
    }


    unsigned long poschromosome = std::pow(2, ell);
    double fitness[poschromosome];
    double a = 0;
    double b = 1;
    int notSO=0;
    int MSOnum=0;
    bool isSOnonWeak = false;
    bool notMSOnonWeak = false;

    for(unsigned long i=0; i<n; i++){
        buffer.clear();
        for(unsigned long i=0; i < std::pow(2, ell); i++){
            Chromosome c;
            c.init(ell);
            c.setGene(0, i);
            buffer.push_back(c);
        }
        for(unsigned long i=0; i<poschromosome; i++) buffer[i].set_fitness(myRand.uniform(a, b));
        for(unsigned long i=0; i<poschromosome; i++) fitness[i] = buffer[i].getFitness();
        if(GHC==1) GHC_inclosure(buffer, fitness);
        std::sort(buffer.begin(), buffer.end(), customSortfitness);
        bool isinP[poschromosome];
        for(auto &i: isinP) i=false;
        for(auto &i: buffer) isinP[i.getGene(0)] = true;

        if(counter == 100) {
            counter = 0;
            printf("Keep working\n");
        }

        for (int i=0; i<ell; i++){
            for (int j=0; j<ell; j++){
                graph[i][j] = false;
            }
        }

        // u -> v ?
        for (int u=0; u<ell; ++u) {
            for (int v=0; v<ell; ++v) {

                int val = 1 - buffer[0].getVal(u); // flip bit
                
                int w = 1;
                while (buffer[w].getVal(u) != val){
                    w++;
                }

                if (buffer[w].getVal(v) != buffer[0].getVal(v)) // if change bit in position u, change bit in position v too
                    graph[u][v] = true;
            }

        }

        std::vector<int> inclosure;
        bool visited[ell];
        maxInDegree = 0;
        bool isweak = false;
        //if(isweak) outputFile << "Weak exists in random assignment: " << i << "\n";

        for (int u=0; u<ell; ++u) {
            // find inclosure for each vertex
            inclosure.clear();
            for(int i=0; i<ell; i++) visited[i] = false;

            findInclosure(inclosure, graph, visited, u);
            if(inclosure.size() > maxInDegree) {
                // printf("\nIN(%u): ", u);
                // printf("\nNew max closure: ");
                // for(auto i: inclosure) std::cout << i;
                maxInDegree = inclosure.size();  
            }
            // for(auto i: inclosure) std::cout << i;
            // std::cout << "\n";
            if(!isSO(inclosure, fitness, buffer, isinP)) {
                notSO++;
                outputFile << "\nThis is not SO\n";
                outputFile << "n=" << i <<", IN(" << u << ")" <<"\n";
                for(auto i: inclosure) outputFile << i << ", ";
                outputFile << "indegree=" << inclosure.size() << "\n";
                //dumpbuffer(buffer);
                //printGraph(graph);
                // isweak = isWeak(buffer);
                // if(!isweak) {
                //     outputFile << "Inclosure doesn't equal to SO and non weak\n";
                //     isSOnonWeak = true;
                // }
            }
            else if(isMSO(inclosure, fitness, buffer, isinP)) {
                // printf("it's MSO\n");
                // printf("n= %lu, IN(%u): ",i, u);
                // for(auto i: inclosure) std::cout << i << ", ";
                // printf("indegree=%lu\n\n", inclosure.size());

                MSOnum++;
            }
            else{
                isweak = isWeak(buffer);
                if(!isweak) {
                    outputFile << "Inclosure doesn't equal to SO and non weak\n";
                    notMSOnonWeak = true;
                }
            }
            // printf("IN(%u): ", u);
            // for(auto i: inclosure) std::cout << i;
            // printf(", indegree=%lu\n ", inclosure.size());
        }
        
        // if(isweak) dumpbuffer(buffer);
        // printf("\nmax inclosure size: %lu\n", maxInDegree);
        counter++;
    }
    if(isSOnonWeak) printf("Inclosure doesn't equal to SO case exists if under non weak condition\n");
    else printf("Inclosure doesn't equal to SO case doesn't exist if under non weak condition\n");

    if(notMSOnonWeak) printf("Inclosure doesn't equal to MSO case exists if under non weak condition\n");
    else printf("Inclosure doesn't equal to MSO case doesn't exist if under non weak condition\n");
    printf("not SO rate: %f%%\n", 100.0*notSO/(n*ell));
    printf("Average SO but not MSO rate: %f%%\n", 100.0*(n*ell-notSO-MSOnum)/(n*ell));
    printf("Average MSO rate: %f%%\n", 100.0*MSOnum/(n*ell));
    
    
    printf("\nmax inclosure size: %ld\n", maxInDegree);
    outputFile.close();
    return 0;
}

void dumpbuffer(std::vector<Chromosome> buffer){
    outputFile << "\nbuffer content:\n";
    std::bitset<sizeof(unsigned long) * 8> bits_best(buffer[0].getGene(0));
    for(const Chromosome &c: buffer) {
        std::bitset<sizeof(unsigned long) * 8> bits(c.getGene(0));
        bits = ~(bits^bits_best);
        std::string lastBits = bits.to_string().substr(bits.size() - ell);
        outputFile << std::setw(5) << c.getGene(0) << ": " << std::setw(5) << std::setprecision(5) 
        << std::fixed << c.getFitness() << ", " << lastBits <<"\n";
    }
    outputFile << "=========================\n\n";

}

void findInclosure(std::vector<int> &inclosure, bool *graph[], bool visited[], int index){
    for(int i=0; i<ell; i++){
        if(!visited[i] && graph[i][index]){
            visited[i] = true;
            inclosure.push_back(i);
            findInclosure(inclosure, graph, visited, i);
        } 
    }
}

void printGraph(bool *graph[]){
    printf("Graph:\n");
    for (int u=0; u<ell; ++u) {
        for (int v=0; v<ell; ++v) {
            if(graph[u][v]) printf("O");
            else printf("X");
        }
        printf("\n");
    }
}

bool isSO(std::vector<int> &inclosure, unsigned long rank[], std::vector<Chromosome> &buffer, bool isinP[]){    
    sort(inclosure.begin(), inclosure.end());
    std::vector<int> notinclosure; 

    // seperate not inclosure and inclosure
    unsigned long index = 0;
    for(int i=0; i<ell; i++){
        if(i==inclosure[index]) index++;
        else notinclosure.push_back(i);
    }
    // printf("\n\nindegree=%lu\n", inclosure.size());
    // for(auto i: inclosure) std::cout << i;
    // printf("\n");
    // for(auto i: notinclosure) std::cout << i;

    unsigned long hash;
    for(int i=0; i<pow(2, notinclosure.size()); i++){
        hash = 0;
        unsigned long index = i;
        for(unsigned long  j=0; j<notinclosure.size(); j++){
            if(index%2==1) hash += pow(2, notinclosure[j]);
            index >>= 1;
        }// get the hash of no inclosure part

        unsigned long best=-1;
        unsigned long besthash=hash;
        for(unsigned long  j=0; j<inclosure.size(); j++){
            if(buffer[0].getVal(inclosure[j])==1) besthash += pow(2, inclosure[j]);
        }
        //printf("best hash: %u, best hash rank: %u\n", besthash, best);

        best = rank[besthash];

        for(unsigned long  j=0; j<pow(2, inclosure.size()); j++){
            int hash_ = hash;
            index = j;
            for(unsigned long  k=0; k<inclosure.size(); k++){
                if(index%2==1) hash_ += pow(2, inclosure[k]);
                index >>= 1;
            }

            //if(isinP[hash_]){
                if( best > rank[hash_]) {
                    //printf("This is not SO\n");
                    return false;
                }
            //}
        }
    }
    return true;
}

bool isSO(std::vector<int> &inclosure, double fitness[], std::vector<Chromosome> &buffer, bool isinP[]){    
    sort(inclosure.begin(), inclosure.end());
    std::vector<int> notinclosure; 

    // seperate not inclosure and inclosure
    unsigned long index = 0;
    for(int i=0; i<ell; i++){
        if(i==inclosure[index]) index++;
        else notinclosure.push_back(i);
    }
    // printf("\n\nindegree=%lu\n", inclosure.size());
    // for(auto i: inclosure) std::cout << i;
    // printf("\n");
    // for(auto i: notinclosure) std::cout << i;

    unsigned long hash;
    for(int i=0; i<pow(2, notinclosure.size()); i++){
        hash = 0;
        unsigned long index = i;
        for(unsigned long  j=0; j<notinclosure.size(); j++){
            if(index%2==1) hash += pow(2, notinclosure[j]);
            index >>= 1;
        }// get the hash of no inclosure part

        double best=-1.0;
        unsigned long besthash=hash;
        for(unsigned long  j=0; j<inclosure.size(); j++){
            if(buffer[0].getVal(inclosure[j])==1) besthash += pow(2, inclosure[j]);
        }

        best = fitness[besthash];

        for(unsigned long  j=0; j<pow(2, inclosure.size()); j++){
            int hash_ = hash;
            index = j;
            for(unsigned long  k=0; k<inclosure.size(); k++){
                if(index%2==1) hash_ += pow(2, inclosure[k]);
                index >>= 1;
            }

            if(best < fitness[hash_]) {
                //printf("This is not SO\n");
                return false;
            }
        }
    }
    return true;
}

bool customSortfitness (const Chromosome &a, const Chromosome &b) {
    return a.getFitness() > b.getFitness();
}

unsigned long long factorial(int n){
    if(n<0) return -1;
    else{
        if(n==1) return 1;
        else return n*factorial(n-1);
    }
}

bool isMSO(std::vector<int> &inclosure, unsigned long rank[], std::vector<Chromosome> &buffer, bool isinP[]){
    if(inclosure.empty()) return true;
    
    for(int k=1; k<std::pow(2, inclosure.size())-1; k++){
        int index = k;
        std::vector<int> inclosure_;
        for(size_t i=0; i<inclosure.size(); i++){
            if(index%2==1) inclosure_.push_back(inclosure[i]);
            index>>=1; 
        }
        if(isSO(inclosure_, rank, buffer, isinP)) {
            printf("Not MSO, SO of length=%ld exists\n\n", inclosure_.size());
            return false;
        }
    }

    return true;

}

bool isMSO(std::vector<int> &inclosure, double fitness[], std::vector<Chromosome> &buffer, bool isinP[]){
    if(inclosure.empty()) return true;
    
    for(int k=1; k<std::pow(2, inclosure.size())-1; k++){
        int index = k;
        std::vector<int> inclosure_;
        for(size_t i=0; i<inclosure.size(); i++){
            if(index%2==1) inclosure_.push_back(inclosure[i]);
            index>>=1; 
        }
        if(isSO(inclosure_, fitness, buffer, isinP)) {
            printf("Not MSO, SO of length=%ld exists\n\n", inclosure_.size());
            return false;
        }
    }

    return true;

}

void GHC_inclosure(std::vector<Chromosome> &buffer, double fitness[]){
    std::shuffle(buffer.begin(), buffer.end(), rng);

    for(Chromosome &c:buffer){
        int* order = new int [ell];
        myRand.uniformArray(order, ell, 0, ell-1);

        for (int i=0; i<ell; ++i) {
            printf("%lu", c.getGene(0));
            double oldF = fitness[(int)c.getGene(0)];
            //unsigned long oldGene = c.getGene(0);
            c.flip(order[i]);
            if (fitness[c.getGene(0)] <= oldF) {
                c.flip(order[i]);
            } else {
                c.set_fitness(fitness[c.getGene(0)]);
                //fitness[oldGene] = fitness[c.getGene(0)];
            }
        }

        delete []order;
    }
}

void findEpitasis(std::vector<Chromosome> &buffer){
    for(unsigned long i=1; i<std::pow(2, ell)-1; i++){ // i is a mask
        bool found = false;
        for(unsigned long k=1; k<buffer.size(); k++){
            if(found) break;
            unsigned long mask = i;                       
            if(((buffer[0].getGene(0)^buffer[k].getGene(0)) & mask)==mask){
                mask = ~mask;
                for(int j=0; j<ell; j++){
                    if(mask%2==1 && ((buffer[0].getVal(j)^buffer[k].getVal(j))==1)){   
                        epistasis[j].push_back(i);
                    }
                    mask >>= 1;
                }
                found = true;
            }
        }
    }
}

bool isWeak(std::vector<Chromosome> &buffer){
    findEpitasis(buffer);
    // for(int i=0; i<ell; i++){
    //     for(auto &j: epistasis[i]){
    //         printEpistasis(j, i);
    //     }
    // }

    bool isweak = false;
    for(int k=0; k<ell; k++){        
        for(unsigned long v=0; v<epistasis[k].size(); v++){
            bool flag = false; // flag represent if a subset points to the same point
            if ((epistasis[k][v] & (epistasis[k][v]-1))!=0){
                for(unsigned long r=0; r<epistasis[k].size(); r++){
                    //unsigned long r = epistasis[k][r_];
                    if(flag) break; 

                    if(epistasis[k][r]!=epistasis[k][v]){
                        //exist a subset points to the same point
                        if((epistasis[k][r]&~epistasis[k][v])==0UL) {
                            flag = true;
                            isweak = true;
                        }
                    }
                }
                if(!flag) printEpistasis(epistasis[k][v], k);
            }
        }
    }
    return isweak;
}

void printEpistasis(int v, int k){
    int index = v;
    //outputFile << "This is weak epistatsis:"; 
    outputFile << "{";
    for(int i=0; i<ell; i++){
        if(index%2 ==1) outputFile << i << " ";
        index >>= 1;
    }
    outputFile << "} -> {" << k << "}\n";
}
