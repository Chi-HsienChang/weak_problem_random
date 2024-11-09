#include "inclosure.h"



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
    delete[] count_chrom;
    delete[] visited;
}

void Inclosure::execute(){

    std::vector<int> total_weak(ell, 0);
    std::vector<int> onetime_weak(ell, 0);

    for (int i=0; i<n; i++){
        if(debug) printf("----------------New run----------------\n");
        onetime_weak = run();
        for(int j=1; j<ell; j++){
            if(debug) printf("Onetime epi_sizw %u # = %u\n", j, onetime_weak[j]);
            total_weak[j] += onetime_weak[j];
        }
    }

    if(debug) printf("----------------Total runs----------------\n");
    if(debug) {
        for(int epi_size=1; epi_size<ell; epi_size++)
        {
            printf("Total epi_sizw %u # = %u\n", epi_size, total_weak[epi_size]);
        }
    }else{
        for(int epi_size=1; epi_size<ell; epi_size++)
        {
            printf("%u ", total_weak[epi_size]);
        }
    }

    // std::cout<<"EEEEEENNNNNNDDDDDD"<<std::endl;
}


std::vector<int> Inclosure::run(){

    int t = std::pow(2, sample_size);
    if (debug) std::cout<<"sample "<<t<<" chromosome"<<std::endl;

    if(debug){
        switch(GHC){
            case 0:
                std::cout << "No GHC\n";
                break;
            case 1:
                std::cout << "One-bit GHC\n";
                break;
            case 2:
                std::cout << "Two-bit GHC\n";
                break;
            default:
                break;
        }
    }



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


    // set up the fitness by random
    for(unsigned long i=0; i< std::pow(2, ell); i++) { 
        double fit = dis(gen);
        buffer[i].set_fitness(fit);
        fitness[i] = fit;
        count_chrom[i] = 0;
        seed++;
    }

    int test_problem_debug = 0;

    if (test_problem_debug){


        //trap debug
        // one max with one epi debug
        unsigned long enum_size = std::pow(2, ell);

        for(unsigned long i=0; i< std::pow(2, ell); i++) { 
            // double fit = dis(gen);
            double fit = buffer[i].evaluate_inclosure2(1);
            buffer[i].set_fitness(fit);
            fitness[i] = fit;
            count_chrom[buffer[i].getGene(0)]++;
            seed++;
        }  

        std::sort(buffer.begin(), buffer.end(), Inclosure::customSortfitness); 
        dumpbuffer();




        // one max with one epi debug
        // unsigned long enum_size = std::pow(2, ell);

        // for(unsigned long i=0; i< std::pow(2, ell); i++) { 
        //     // double fit = dis(gen);
        //     double fit = buffer[i].evaluate_inclosure2(0);
        //     buffer[i].set_fitness(fit);
        //     fitness[i] = fit;
        //     count_chrom[buffer[i].getGene(0)]++;
        //     seed++;
        // }  

        // std::sort(buffer.begin(), buffer.end(), Inclosure::customSortfitness); 

        // buffer[enum_size-1].set_fitness(1);
        // std::cout<<" buffer[enum_size-1].set_fitness(1): "<< fitness[enum_size-1] <<std::endl;
        // fitness[enum_size-1] = 1;
        // std::cout<<" buffer[enum_size-1].set_fitness(1): "<< fitness[enum_size-1] <<std::endl;
        // buffer[enum_size-2].set_fitness(0);
        // std::cout<<" buffer[enum_size-1].set_fitness(0): "<< fitness[enum_size-2] <<std::endl;
        // fitness[enum_size-2] = 0;
        // std::cout<<" buffer[enum_size-1].set_fitness(0): "<< fitness[enum_size-2] <<std::endl;
        // dumpbuffer();

    }else{

        std::random_shuffle(buffer.begin(), buffer.end());

        // sample t chromosome
        int i = 0;
        for (auto& c: buffer)
        {
            if(i == t) break;
            count_chrom[c.getGene(0)]++;
            i++;
        }

    }

    // count_chrom[buffer[0].getGene(0)]++;
    // count_chrom[buffer[7].getGene(0)]++;

    // std::sort(buffer.begin(), buffer.end(), Inclosure::customSortfitness);

    // dumpbuffer();

    // dumpbuffer();

    if(debug) std::sort(buffer.begin(), buffer.end(), Inclosure::customSortfitness);

    

    if(debug) dumpbuffer();

    ghc();

    std::sort(buffer.begin(), buffer.end(), Inclosure::customSortfitness);

    if(GHC && debug) dumpbuffer();


    


    // set up the graph
    for (int i=0; i<ell; i++){
        for (int j=0; j<ell; j++){
            graph[i][j] = false;
        }
    }

    std::vector<std::vector<std::vector<int>>> weak_epi(ell);
    std::vector<int> one_weak(ell, 0); 

    for (b=0; b< std::pow(2, ell);b++) {
        if(count_chrom[buffer[b].getGene(0)]) break;
    }

    
    if (debug) std::cout<<"best index: "<<b<<std::endl;
    
    if (debug){
        for (int u=ell-1; u>=0; u--) {
            
            std::cout<<buffer[b].getVal(u)<<" "; // flip bit
        }
        std::cout<<std::endl;
    }
     

    for (int u=1; u<ell; u++) {
        int val = 1 - buffer[b].getVal(u); // flip bit

        if(t==1){
            // std::cout<<"t=1"<<std::endl;
            break;
        }    
        int w = b+1;

        if (debug) std::cout<<"w index: "<<w<<std::endl;

        if(w > std::pow(2, ell)-1) break;

        while (buffer[w].getVal(u) != val || count_chrom[buffer[w].getGene(0)] == 0){
            if (debug) std::cout<<"while: "<<w<<std::endl;
            w++;
            if(w > std::pow(2, ell)-1) break;
        }

        if (debug) std::cout<<"w index: "<<w<<std::endl;
        if (debug) std::cout<<"======="<<std::endl;
        if(w > std::pow(2, ell)-1) continue;

        
        if (debug) std::cout<<"buffer[w].getVal(0): "<<buffer[w].getVal(0)<<std::endl;
        if (debug) std::cout<<"buffer[b].getVal(0): "<<buffer[b].getVal(0)<<std::endl;
        if (debug) std::cout<<"count_chrom[buffer[w].getGene(0)] "<<count_chrom[buffer[w].getGene(0)]<<std::endl;
        if (buffer[w].getVal(0) != buffer[b].getVal(0) && count_chrom[buffer[w].getGene(0)] > 0)
        { // if change bit in position u, change bit in position v too
            if(debug) printf("%u -> 0\n", u);
            std::vector<int> u_vector(1, u);
            weak_epi[1].push_back(u_vector);
            one_weak[1]++;
            

            graph[u][0] = true;
        }
        
    }


    // if(debug){
    //     printf("Graph:\n");
    //     for (int u=0; u<ell; u++) {
    //         for (int v=0; v<ell; v++) {
    //             if(graph[u][v]){
    //                 printf("O");
    //             }
    //             else{
    //                 printf("X");
    //             }
    //         }
    //         printf("\n");
    
    //     }       
    // }

    // for (auto& node : weak_epi[1]){
    //     std::cout<<"node: "<<node<<" ";
    // }
    // std::cout<<std::endl;


    
    
    for (int epi_size = 2; epi_size < ell; epi_size++)
    {  
        auto combinations = generateCombinations(ell-1, epi_size); // combinations = { [1, 2], [1, 3], [2, 3] }
        
        for (auto& combination : combinations) // combination = [1, 2]
        { 
            // std::cout<<"combination: "<<std::endl;
            // for (auto& combination : combinations){
            //      for (auto& index : combination){
            //         std::cout<<index<<" ";
            //      }
            //     std::cout<<std::endl;
            //     std::cout<<"----------------"<<std::endl;
            // }
            //  std::cout<<"+++++++++++++++++++"<<std::endl;



            bool not_find_smaller_epi = true;

            int smaller_epi_size = epi_size;
            bool is_subset;
            while(not_find_smaller_epi && smaller_epi_size>=2)
            {
                for (auto& previous : weak_epi[smaller_epi_size-1]) 
                {
                    is_subset = isSubset(previous, combination);
                    not_find_smaller_epi = !is_subset;
                    if(is_subset) break;
                }
                smaller_epi_size--;
            }

            

            
            // std::cout<<"epi_size "<<epi_size<<std::endl;
            // std::cout<<"not_find_smaller_epi "<<not_find_smaller_epi<<std::endl;
            if(not_find_smaller_epi)
            {
                // std::cout<<"enter if "<<std::endl;
                auto enumerations = generateBinarySequences(epi_size); // enumerations = { [0, 0], [0, 1], [1, 0], [1, 1] }
                std::vector<int> best_sequence(epi_size);
                
                int j = 0;
                for (auto num : combination)
                { // combination = [1, 2]
                    // std::cout << "b "<< b<<std:: endl;
                    best_sequence[j] = buffer[b].getVal(num);  // best_sequence = [g[1], g[2] ]
                    j++;
                }
                
                auto enumeration_without_best = remove_the_same_with_best(enumerations, best_sequence); // enumeration_without_best = { [0, 0], [0, 1], [1, 0]}

                // std::cout <<"-------------"<< std:: endl;

                

                for (auto& enumeration : enumeration_without_best) // combination = [1, 2] and enumeration_without_best = { [0, 0], [0, 1], [1, 0]}
                { 

                    int chrom_index = b+1;
                    if(chrom_index > (std::pow(2, ell)-1)) break;     

                    while(count_chrom[buffer[chrom_index].getGene(0)] <= 0){
                        chrom_index++;
                        if(chrom_index > (std::pow(2, ell)-1)) break;
                    }

                    // if(debug) printf("############## DEBUG ##############\n");

                    int hit = 0;
                    while(1)
                    {  
                        if(chrom_index > (std::pow(2, ell)-1)) break;

                        hit = 0;
                        for(int i = 0; i < enumeration.size(); i++){
                            if(buffer[chrom_index].getVal(combination[i]) == enumeration[i]) hit++;
                        }
                        
                        if(hit == epi_size) break;
                        if(hit < epi_size) chrom_index++;    
                    }
                    
                    if(chrom_index > (std::pow(2, ell)-1)) continue;

                    if((buffer[chrom_index].getVal(0) != buffer[b].getVal(0)) && count_chrom[buffer[chrom_index].getGene(0)] > 0)
                    { 
                        // std::cout << "index: ";
                        // for(int i = 0; i < enumeration.size(); i++){
                        // std::cout << combination[i]<<" ";
                        // }

                        // std::cout << "allel: ";
                        // for(int i = 0; i < enumeration.size(); i++){
                        // std::cout << enumeration[i]<< " ";
                        // }
                        // std::cout << std:: endl;

                        chrom_index += 1;
                        if(chrom_index > (std::pow(2, ell)-1)) continue;     

                        while(count_chrom[buffer[chrom_index].getGene(0)] <= 0){
                            chrom_index++;
                            if(chrom_index > (std::pow(2, ell)-1)) break;
                        }                  

                        int hit = 0;
                        while(1)
                        {  
                            if(chrom_index > (std::pow(2, ell)-1)) break;

                            hit = 0;
                            for(int i = 0; i < enumeration.size(); i++){
                                if(buffer[chrom_index].getVal(combination[i]) == enumeration[i]) hit++;
                            }
                            
                            if(hit == epi_size) break;
                            if(hit < epi_size) chrom_index++;    
                        }

                        if(chrom_index > (std::pow(2, ell)-1)) continue;   


                        if((buffer[chrom_index].getVal(0) == buffer[b].getVal(0)) && count_chrom[buffer[chrom_index].getGene(0)] > 0)
                        {
                            one_weak[epi_size]++;
                            weak_epi[epi_size].push_back(combination);

                            if(debug){
                                printf("{ ");
                                for (const auto& vertex : combination)
                                {
                                    printf("%u ", vertex);
                                }
                                printf("} -> 0\n");
                            }
                            break;
                        }
                    };

                }
            }

        }
    }

    // std::cout<<"END"<<std::endl;

    return one_weak;
}


bool Inclosure::isSubset(const std::vector<int>& subset, const std::vector<int>& set) {
    for (const auto& elem : subset) {
        if (std::find(set.begin(), set.end(), elem) == set.end()) {
            // 如果找不到元素，則返回 false
            return false;
        }
    }
    // 所有元素都被找到，返回 true
    return true;
}

std::vector<std::vector<int>> Inclosure::remove_the_same_with_best(const std::vector<std::vector<int>>& enumerations, const std::vector<int>& best_sequence) {
    std::vector<std::vector<int>> result = enumerations;  

    auto newEnd = std::remove_if(result.begin(), result.end(),
                                 [&best_sequence](const std::vector<int>& sequence) {
                                     return sequence == best_sequence;
                                 });
    result.erase(newEnd, result.end());

    return result;
}

std::vector<std::vector<int>> Inclosure::generateBinarySequences(int n) {
    int totalSequences = 1 << n;  // 2^n
    std::vector<std::vector<int>> allSequences;

    for (int i = 0; i < totalSequences; ++i) {
        std::vector<int> sequence(n);

        for (int j = 0; j < n; ++j) {
            // 檢查第j位是否為1
            sequence[n - 1 - j] = (i >> j) & 1;
        }

        // 將序列添加到所有序列的vector中
        allSequences.push_back(sequence);
    }

    return allSequences;
}

std::vector<std::vector<int>> Inclosure::generateCombinations(int n, int k) {
    std::vector<int> bitmask(k, 1);  // 創建k個1
    bitmask.resize(n, 0);  // 後面填充n-k個0

    std::vector<std::vector<int>> combinations;

    do {
        std::vector<int> currentCombination;
        for (int i = 0; i < n; ++i) {
            if (bitmask[i]) {
                currentCombination.push_back(i + 1);  // 輸出選中的元素
            }
        }
        combinations.push_back(currentCombination);
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));  // 生成下一個排列

    return combinations;
}

bool Inclosure::customSortfitness (const Chromosome &a, const Chromosome &b) {
    return a.getFitness() > b.getFitness();
}

void Inclosure::dumpbuffer(){
    std::cout << "====================================";
    std::cout << "\nbuffer content:\n";
    unsigned long bits_best = buffer[0].getGene(0);
    for(const Chromosome &c: buffer) {
        unsigned long bits = c.getGene(0);
        bits = ~(bits^bits_best);
        // std::cout << c.getGene(0) << ": " << std::setw(5) << std::setprecision(5) 
        // << std::fixed << c.getFitness() << ", ";

        std::cout << std::setw(5) << std::setprecision(5) 
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
        std::cout << " count: "<< count_chrom[c.getGene(0)];
        std::cout <<"\n";
    }
    std::cout << "====================================\n";
    
    // printGraph();
}


void Inclosure::ghc()
{
    switch (Inclosure::GHC)
    {
        case 0:
            // printf("no GHC\n");
            break;
        case 1:
            if(debug) printf("After one-bit GHC\n");
            std::random_shuffle(buffer.begin(), buffer.end());
            for (auto c: buffer)
            {
                if(count_chrom[c.getGene(0)] > 0)
                {
                    int* order = new int [ell];
                    myRand.uniformArray(order, ell, 0, ell-1);

                    if(debug){
                        std::cout<<"GHC order: ";
                        for (int i=0; i<ell; ++i){
                            std::cout<<order[i]<<" ";
                        }
                        std::cout<<std::endl;

                        std::cout<<"The chromosome:";
                        for (int i=ell-1; i>=0; i--){
                            std::cout<<c.getVal(i)<<" ";
                        }
                        std::cout<<std::endl;
                    }

                    for (int i=0; i<ell; ++i) 
                    {    
                        if(count_chrom[c.getGene(0)] > 0)
                        {
                            double oldF = fitness[c.getGene(0)];
                            c.flip_gene0(order[i]);   
                            if(debug) {
                                std::cout<<"Flip Chromosome:";
                                for (int i=ell-1; i>=0; i--){
                                    std::cout<<c.getVal(i)<<" ";
                                }
                                std::cout<<std::endl;
                            }
                            if (fitness[c.getGene(0)] > oldF) {
                                count_chrom[c.getGene(0)]++;
                                c.flip_gene0(order[i]);   
                                count_chrom[c.getGene(0)]--;
                                c.flip_gene0(order[i]);  
                            }else{
                                c.flip_gene0(order[i]);
                            }
                        }
                }

                delete[] order;

                }

            }
            break;
        case 2:
            if(debug) printf("After two-bit GHC\n");
            std::random_shuffle(buffer.begin(), buffer.end());
            for (auto c: buffer)
            {
                if(count_chrom[c.getGene(0)] > 0)
                {   if(debug) std::cout<<"------New c-------"<<std::endl;
                    auto combinations = generateCombinations(ell, 2);

                    for (auto& combination : combinations) 
                    {
                        for (auto& index : combination) 
                        {
                            index--;
                            if(debug) std::cout<<index<<" ";
                        }
                        if(debug) std::cout<<std::endl;
                    }

                    std::random_shuffle(combinations.begin(), combinations.end());

                    // std::cout<<"------After shuffle-------"<<std::endl;
                    // for (auto& combination : combinations) 
                    // {
                    //     for (auto& index : combination) 
                    //     {
                    //         std::cout<<index<<" ";
                    //     }
                    //     std::cout<<std::endl;
                    // }

                    double f1 = 0;
                    double f2 = 0;;
                    double f3 = 0;;                    
                    double f4 = 0;;


                    if(debug){
                        std::cout<<"GHC order: \n";

                        for (auto& combination : combinations) 
                        {
                            for (auto& index : combination) 
                            {
                                std::cout<<index<<" ";
                            }
                            std::cout<<std::endl;
                        }

                        std::cout<<"The chromosome:";
                        for (int i=ell-1; i>=0; i--){
                            std::cout<<c.getVal(i)<<" ";
                        }
                        std::cout<<std::endl;
                    }                   

                    for (auto& combination : combinations) // combination = [1, 2]
                    {   
                        if(count_chrom[c.getGene(0)] > 0)
                        {
                            f1 = fitness[c.getGene(0)];
                            // std::cout<<"f1 = "<< f1 <<std::endl;

                            c.flip_gene0(combination[0]);
                            f2 = fitness[c.getGene(0)];
                            // std::cout<<"f2 = "<< f2 <<std::endl;
                            c.flip_gene0(combination[0]);

                            c.flip_gene0(combination[1]);
                            f3 = fitness[c.getGene(0)];
                            // std::cout<<"f3 = "<< f3 <<std::endl;
                            c.flip_gene0(combination[1]);

                            c.flip_gene0(combination[0]);
                            c.flip_gene0(combination[1]);
                            f4 = fitness[c.getGene(0)];
                            // std::cout<<"f4 = "<< f4 <<std::endl;
                            c.flip_gene0(combination[0]);
                            c.flip_gene0(combination[1]);

                            // 初始化最大值為f1
                            double max = f1;
                            int max_index = 1;

                            // 比較f2，並更新最大值
                            if (f2 > max) 
                            {
                                max = f2;
                                max_index = 2;
                            }

                            // 比較f3，並更新最大值
                            if (f3 > max) 
                            {
                                max = f3;
                                max_index = 3;
                            }

                            // 比較f4，並更新最大值
                            if (f4 > max) 
                            {
                                max = f4;
                                max_index = 4;
                            }

                            // std::cout<<"max = "<< max <<" is f"<<max_index<<std::endl;


                            if(max_index == 1){
                                if(debug){
                                    std::cout<<"Flip Chromosome:";
                                    for (int i=ell-1; i>=0; i--){
                                        std::cout<<c.getVal(i)<<" ";
                                    }   
                                    std::cout<<std::endl;  
                                }                         
                            }else if (max_index == 2)
                            {
                                count_chrom[c.getGene(0)]--;
                                c.flip_gene0(combination[0]); 
                                count_chrom[c.getGene(0)]++;
                                if(debug){
                                    std::cout<<"Flip Chromosome:";
                                    for (int i=ell-1; i>=0; i--){
                                        std::cout<<c.getVal(i)<<" ";
                                    }     
                                    std::cout<<std::endl;
                                }
                                // c.flip_gene0(combination[0]); 
                            }else if (max_index == 3)
                            {
                                count_chrom[c.getGene(0)]--;
                                c.flip_gene0(combination[1]); 
                                count_chrom[c.getGene(0)]++;
                                if(debug){
                                    std::cout<<"Flip Chromosome:";
                                    for (int i=ell-1; i>=0; i--){
                                        std::cout<<c.getVal(i)<<" ";
                                    }   
                                    std::cout<<std::endl;
                                 }
                                // c.flip_gene0(combination[1]); 
                            }else
                            {
                                count_chrom[c.getGene(0)]--;
                                c.flip_gene0(combination[0]); 
                                c.flip_gene0(combination[1]); 
                                count_chrom[c.getGene(0)]++;
                                if(debug){
                                    std::cout<<"Flip Chromosome:";
                                    for (int i=ell-1; i>=0; i--){
                                        std::cout<<c.getVal(i)<<" ";
                                    }   
                                    std::cout<<std::endl;
                                }
                                // c.flip_gene0(combination[0]); 
                                // c.flip_gene0(combination[1]); 
                            }
                        }
                    }
                }

            }
            break;
        default:
            break;
    }
}
