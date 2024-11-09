#include "inclosure.h"
/**
 * @param ell_ The problem size
 * @param n_ The number of examinations
 * @param GHC_ How many bits GHC(greedy hill climbing) is needed.(if GHC is 0 implies no ghc in the process)
 * @param seed_ random seed of the program
 * In the Inclosure class, the examination of SO and MSO in the population is implemented.
*/

Inclosure::Inclosure(int ell_, int n_, int GHC_, int debug_, int seed_, int sample_size_):
ell(ell_), n(n_), GHC(GHC_), debug(debug_), seed(seed_), sample_size(sample_size_)
{
    unsigned long poschromosome = std::pow(2, ell);
    graph = new bool*[ell];
    for (int i=0; i<ell; i++) graph[i]=new bool[ell];
    fitness = new double[poschromosome];
    count_chrom = new int[poschromosome];
    visited = new bool[ell];

}

Inclosure::~Inclosure(){
    for (int i=0; i<ell; i++) delete[] graph[i];
    delete[] graph;
    delete[] fitness;
    delete[] visited;
}

void Inclosure::execute(){

    for (unsigned long i=0; i<n; i++){
        run();
    }

}

std::vector<int> Inclosure::run(){
    //set up buffer
    buffer.clear();
    for(unsigned long i=0; i < std::pow(2, ell); i++){
        Chromosome c;
        c.init(ell);
        c.setGene(0, i);
        buffer.push_back(c);
    }
    std::mt19937 gen(seed); 
    std::uniform_real_distribution<double> dis(fitness_lower, fitness_upper);
    
    for(unsigned long i=0; i< std::pow(2, ell); i++) { 
        // double fit = dis(gen);
        double fit = buffer[i].evaluate_inclosure2(7);
        buffer[i].set_fitness(fit);
        fitness[i] = fit;
        seed++;
    }

    // for(unsigned long i=0; i< std::pow(2, ell); i++) { 
    //     double fit = dis(gen);
    //     buffer[i].set_fitness(fit);
    //     fitness[i] = fit;
    //     seed++;
    // }

    // ghc();

    std::sort(buffer.begin(), buffer.end(), Inclosure::customSortfitness);

    // set up the graph
    for (int i=0; i<ell; i++){
        for (int j=0; j<ell; j++){
            graph[i][j] = false;
        }
    }

    for (int u=0; u<ell; ++u) 
    {
        int first = 1;
        double my_fitness;
        int val = 1 - buffer[0].getVal(u); // flip bit
        int w = 1;
        // std::cout<<"-----------------"<<std::endl;
        do 
        {
            while (buffer[w].getVal(u) != val){
                if (w == std::pow(2, ell)-1) break;
                w++;
            }

            if(first == 1){
               my_fitness = fitness[buffer[w].getGene(0)];
            }

            if(fitness[buffer[w].getGene(0)] == my_fitness){
                for (int v=0; v<ell; ++v) {
                    if (buffer[w].getVal(v) != buffer[0].getVal(v)){ // if change bit in position u, change bit in position v too
                        printf("%u -> %u\n", u, v);
                        graph[u][v] = true;
                    }
                }
            }

            first++;

            if(w < std::pow(2, ell)-2) w++;
            // std::cout<<fitness[buffer[w-1].getGene(0)]<<std::endl;
            // std::cout<<fitness[buffer[w].getGene(0)]<<std::endl;
        } while (fitness[buffer[w-1].getGene(0)] == fitness[buffer[w].getGene(0)] && w < std::pow(2, ell)-2);
       
        // while (buffer[w].getVal(u) != val) w++;
        // for (int v=0; v<ell; ++v) {
        //     if (buffer[w].getVal(v) != buffer[0].getVal(v)){ // if change bit in position u, change bit in position v too
        //         printf("%u -> %u\n", u, v);
        //         graph[u][v] = true;
        //     }
        // }

    }

    dumpbuffer();

    printf("Graph:\n");
    for (int u=0; u<ell; ++u) {
        for (int v=0; v<ell; ++v) {
            if(graph[u][v]) printf("O");
            else printf("X");
        }
        printf("\n");
    }

    std::vector<int> one_weak(ell, 0); 

    return one_weak;

}


bool Inclosure::customSortfitness (const Chromosome &a, const Chromosome &b) {
    return a.getFitness() > b.getFitness();
}


void Inclosure::dumpbuffer(){
    std::cout << "\nbuffer content:\n";
    unsigned long bits_best = buffer[0].getGene(0);
    for(const Chromosome &c: buffer) {
        unsigned long bits = c.getGene(0);
        bits = ~(bits^bits_best);
        std::cout << std::setw(5) << c.getGene(0) << ": " << std::setw(5) << std::setprecision(5) 
        << std::fixed << c.getFitness() << ", ";
        std::vector<char> b;
        for(int i=0; i<ell; i++){
            if(bits & 1UL == 1) b.push_back('1');
            else b.push_back('0');
            bits >>= 1;
        }
        std::reverse(b.begin(), b.end());
        for(auto &c: b) std::cout << c; 
        std::cout <<", ";
        b.clear();
        unsigned long gene = c.getGene(0);
        for(int i=0; i<ell; i++){
            if(gene & 1UL == 1) b.push_back('1');
            else b.push_back('0');
            gene >>= 1;
        }
        std::reverse(b.begin(), b.end());
        for(auto &c: b) std::cout << c; 
        std::cout <<"\n";
    }
    // std::cout << "=========================\n\n";

    // printGraph();
}



// void Inclosure::printGraph(){
//     printf("Graph:\n");
//     for (int u=0; u<ell; ++u) {
//         for (int v=0; v<ell; ++v) {
//             if(graph[u][v]) printf("O");
//             else printf("X");
//         }
//         printf("\n");
//     }
// }



