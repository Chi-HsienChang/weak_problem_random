#include "chromosome.h"
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>

int ell;
int mode; // mode 0: stop when all graph is supplied, mode 1: statistic of different type of graph
int posgraph;
int existgraph;

std::ofstream outputFile("ontoCheck.txt");

void printGraphState(bool* graphstate[]);
void printGraphStateCheck(int graphstatecheck[]);
void outGraphStateCheck(int graphstatecheck[]);
void printGraph(bool *graph[]);

int main(int argc, char *argv[]){
    
    if(!outputFile.is_open()) printf("File open error");

    std::vector<Chromosome> buffer;

    if(argc!=3) std::cout << "please check the input format";
    ell = atoi (argv[1]);
    mode = atoi (argv[2]);
    for(unsigned long i=0; i < std::pow(2, ell); i++){
        Chromosome c;
        c.init(ell);
        c.setGene(0, i);
        buffer.push_back(c);
    }

    bool *graph[ell];
    for (int i=0; i<ell; i++){
        graph[i] = new bool[ell];
    }

    bool *graphstate[4]; // 4 states for u x v , u -> v, u <- v, u <-> v 
    for (int i=0; i<4; i++){
        graphstate[i] = new bool[ell*(ell-1)/2];
        for (int j=0; j<ell*(ell-1)/2; j++) graphstate[i][j] =false;
    }
    posgraph = std::pow(2, ell*(ell-1));
    int graphstatecheck[posgraph];
    for (int i=0; i<posgraph; i++){
        graphstatecheck[i] = 0;
    }

    if(mode==1){
        int counter = 0;
        do{
            if(counter == 1000000) {
                counter = 0;
                printf("Keep working\n");
            }

            for (int i=0; i<ell; i++){
                for (int j=0; j<ell; j++){
                    graph[i][j] = false;
                }
            }
            for (int i=0; i<4; i++){
                for (int j=0; j<ell*(ell-1)/2; j++){
                    graphstate[i][j] = false;
                }
            }

            for (int u=0; u<ell; u++) {
                for (int v=0; v<ell; v++) {

                    int val = 1 - buffer[0].getVal(u); // flip bit
                    
                    int w = 1;
                    while (buffer[w].getVal(u) != val){
                        w++;
                    }

                    if (buffer[w].getVal(v) != buffer[0].getVal(v)) {// if change bit in position u, change bit in position v too
                        graph[u][v] = true;
                        
                    }
                }
            }

            for (int u=0; u<ell; u++) {
                for (int v=u+1; v<ell; v++) {
                    int index = u*(2*ell-u-1)/2+v-u-1;
                    //printf("index = %u\n", index);

                    if(graph[u][v]){
                        if(graph[v][u]) graphstate[3][index] = true;
                        else graphstate[1][index] = true;
                    }
                    else{
                        if(graph[v][u]) graphstate[2][index] = true;
                        else graphstate[0][index] = true;
                    }
                }
            }
            //printGraph(graph);
            //printGraphState(graphstate);
            unsigned long hash = 0;
            unsigned long shift = 1;
            for(int j=0; j<ell*(ell-1)/2; j++){
                for (int i=0; i<4; i++){
                    if(graphstate[i][j]){
                        hash += i * shift;
                    }
                }
                shift *= 4;
            }
            graphstatecheck[hash] ++;
            counter++;
        } while (std::next_permutation(buffer.begin(), buffer.end()));
    }        
    outGraphStateCheck(graphstatecheck);
    outputFile.close();

    return 0;
}

void printGraphState(bool* graphstate[]){
    for(int i=0; i<4; i++){
        for(int j=0; j<ell*(ell-1)/2; j++){
            if(graphstate[i][j]) printf("O");
            else printf("X");
        }
        printf("\n");
    }
    printf("==========================\n");
}

void printGraphStateCheck(int graphstatecheck[]){
    printf("\nGraphStateCheck: \n");
    for(int i=0; i<posgraph; i++) printf("%u", graphstatecheck[i]);
    printf("\n");
}

void outGraphStateCheck(int graphstatecheck[]){
    existgraph = 0;
    outputFile << "\nGraphStateCheck: \n";
    for(int i=0; i<posgraph; i++) {
        outputFile << graphstatecheck[i];
        if(graphstatecheck[i]!=0) existgraph++;
    }
    outputFile << "\nTotal graph number: " << posgraph <<"\n";
    outputFile << "Possible graph number: " << existgraph <<"\n";

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