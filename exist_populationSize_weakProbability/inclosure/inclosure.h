#pragma once
#include "chromosome.h"
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <random>
#include <bitset>
#include <iomanip>
#include <utility>
#include <memory>

class Inclosure{
public:
    /**
     * @brief Constructor of Inclosure class
     * @param ell_ The problem size
     * @param n_ The number of examinations
     * @param GHC_ How many bits GHC(greedy hill climbing) is needed.(if GHC is 0 implies no ghc in the process)
     * @param seed_ random seed of the program
     * Initialize constant and vector size.
     */
    Inclosure(int ell, int n, int GHC, int debug, int seed, int sample_size);

    /**
     * @brief Destructor of Inclosure class
     * 
     * Clean up the memory used in program.
     */
    ~Inclosure();

    /**
     * @brief Execute one random fitness assignment examination
     * 
     * 1. random fitness assignment to all chromosome in buffer
     * 2. pass the buffer through GHC 
     * 3. Find the epitasis graph of the buffer
     * 4. Examine if buffer weak
     * 5. If non-weak, examine if evrey inclosure for each loci is SO/
     * 6. If inclosure is SO, examine if inclousure MSO. 
     */
    std::vector<pair<double, double>> run();

    /**
     * @brief Execute run() for n times, and calculate the total number of weak, SO, MSO, etc..
     * @ref void Inclosure::run()
    */
    void execute();

private:
    bool **graph; /** < @private The 2-D array of size ell*ell stores the epistais graph for the buffer*/
    double *fitness; /**The array of size 2^ell stores the fitness of any genes*/
    int *count_chrom;
    bool *visited;/**A array of size ell used in findInclosure() to check if a loci is visited*/
    bool isweak; 
    unsigned long poschromosome; /**Constant 2^ell*/
    int seed; /**Random seed*/
    int debug;
    int b;
    int sample_size;

    std::vector<Chromosome> buffer; /**The buffer stores the chromosomes */

    const int ell;
    const unsigned long n;
    const int GHC;
    
    const double fitness_upper = 1.0;
    const double fitness_lower = 0.0;

    /**
     * @brief Provide a custom compare function for std::sort to sort out the buffer
     */
    static bool customSortfitness (const Chromosome &a, const Chromosome &b);

    /**
     * @brief Find the inclosure IN*(v) for given index by recursively called itself.
     */
    unsigned long findInclosure(int index);

    /**
     * @brief Examine if a given inclosure is a SO (stationary optimum)
     * 1. Find every assignment which is not specified by the inclosure.
     * 2. For each assigment, enumereate all pattern in the inclosure.
     * 3. If the patterns of inclosure are consistent under all assignment, the inclosure is a SO.
     *    Otherwise, if the patterns are inconsistent, then it's not a SO. 
    */
    bool isSO(unsigned long inclosure_);

    /**
     * @brief Examine if a given inclosure is a MSO (minimal stationary optimum)
     * @param inclosure_ The target inclosure to be examined
     * @param u The loci of which target inclosure points to 
     * For all subsets of the inclosure, if a true subset includes u is an SO, then the inclosure is not a MSO.
     */
    bool isMSO(unsigned long inclosure_, int u);

    /**
     * @brief Find the epitasis relation in the buffer
     * 1. Given a loci i, which is pointed by a number of sets.
     * 2. For specific i, enumrate combination of all other locus(mask in the program).
     * 3. Examine if the set covered by maks points to i.
     */
    void findEpitasis();

    /**
     * @brief Examine if the buffer weak
     * 1. Given the epistasis graph of buffer, examine if every epistasis realation weak.
     * 2. For a specific epistasis relation, examine if the sets that point to the same loci do not belong to this epistatis.
     * 3. If it's the case for any epistasis relation, then the buffer is weak.
     * 4. If it's not the case for all epistasis relations, then the buffer is not weak.
     */
    bool isWeak();

    /**
     * @brief Print out the content of buffer.
    */
    void dumpbuffer();

    /**
     * @brief Print out the epistatis relations of given index.
    */
    void printEpistasis(int v);

    /**
     * @brief Print out the epistatis graph.
    */
    void printGraph();
    bool isSubset(const std::vector<int>& subset, const std::vector<int>& set);
    std::vector<std::vector<int>> remove_the_same_with_best(const std::vector<std::vector<int>>& enumerations, const std::vector<int>& best_sequence);
    std::vector<std::vector<int>> generateCombinations(int n, int k);
    std::vector<std::vector<int>> generateBinarySequences(int n);
    // void printGraph_index_zero();

    /**
     * @brief For a given incloure, print out inclosure in the format readble.
    */
    void printInclosure(int inclosure_);

    /**
     * @brief Execute GHC in different bits. defalut = 0
    */
    void ghc();
};




