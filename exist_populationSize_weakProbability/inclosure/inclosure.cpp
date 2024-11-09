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

    int t = std::pow(2, sample_size); // t是sample出來的chromosome數量
    if (debug) std::cout<<"sample "<<t<<" chromosome!!"<<std::endl;

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




        // one max with one epi debug
        // ############################

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

    // count_chrom[buffer[0].getGene(0)]++;
    // count_chrom[buffer[7].getGene(0)]++;

    // std::sort(buffer.begin(), buffer.end(), Inclosure::customSortfitness);

    // dumpbuffer();

    // dumpbuffer();

    if(debug) std::sort(buffer.begin(), buffer.end(), Inclosure::customSortfitness);

    

    if(debug) dumpbuffer();

    ghc();

    std::sort(buffer.begin(), buffer.end(), Inclosure::customSortfitness); // 已照fitness排序   

    if(GHC && debug) dumpbuffer();


    // 初始化EG 都捨定成false
    // set up the graph 
    for (int i=0; i<ell; i++){
        for (int j=0; j<ell; j++){
            graph[i][j] = false;
        }
    }

    // 創立 weak_epi 的 vector >>> 裡面有 ell 個 vector >>> vector 裡面有 int 來統計 weak epi 的數量
    std::vector<std::vector<std::vector<int>>> weak_epi(ell);
    std::vector<int> one_weak(ell, 0); 

    
    // 遍歷特定基因的組合，直到找到一個存在的基因為止
    // buffer[b].getGene(0) 是一個整數 例如 13
    for (b=0; b< std::pow(2, ell);b++) {
        if(count_chrom[buffer[b].getGene(0)]) break; // count_chrom 紀錄被sample到的chromosome
    }


    int b_original = b; // b_original 是 被 sample 到 fitness 最高的那一條
    
    if (debug) std::cout<<"best index: "<<b<<std::endl;
    
    if (debug){
        for (int u=ell-1; u>=0; u--) {
            
            std::cout<<buffer[b].getVal(u)<<" "; // flip bit
        }
        std::cout<<std::endl;
    }
     

    // EG 矩陣建立
    for (int u=1; u<ell; u++) { // 因為這裡是討論 指向 index 0 的情況的 epi，所以從 1 開始 
        int val = 1 - buffer[b].getVal(u); // flip bit

        if(t==1){  // t是sample出來的chromosome數量
            // std::cout<<"t=1"<<std::endl;
            break;
        }    
        int w = b+1;
        // if (debug) std::cout<<"w index: "<<w<<std::endl;

        if(w > std::pow(2, ell)-1) break;

        while (buffer[w].getVal(u) != val || count_chrom[buffer[w].getGene(0)] == 0){
            // if (debug) std::cout<<"while: "<<w<<std::endl;
            w++;
            if(w > std::pow(2, ell)-1) break;
        }

        // if (debug) std::cout<<"w index: "<<w<<std::endl;
        // if (debug) std::cout<<"======="<<std::endl;
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




    
    // 遍歷所有組合可能
    for (int epi_size = 2; epi_size < ell; epi_size++)
    {  
        auto combinations = generateCombinations(ell-1, epi_size); // combinations = { [1, 2], [1, 3], [2, 3] }
        
        for (auto& combination : combinations) // combination = [1, 2]
        { 

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
            // 如果沒有找到更小的epi，則進入if >>> 這裡使用 weak 機制排除
            if(not_find_smaller_epi)
            {
                // std::cout<<"enter if "<<std::endl;
                auto enumerations = generateBinarySequences(epi_size); // enumerations = { [0, 0], [0, 1], [1, 0], [1, 1] }    

                // 遍歷所有的 pattern 不用 remove best_sequence 的 pattern
                for (auto& enumeration : enumerations) // combination = [1, 2] and enumeration_without_best = { [0, 0], [0, 1], [1, 0]}
                { 
                    auto enumeration_original = enumeration;

                    // b is 主角pattern 出現 的 index
                    b = b_original;

                    int hit_original = 0; // 這個是用來計算有幾個相同的值
                    while(1)
                    {  
                        if(b > (std::pow(2, ell)-1)) break;

                        hit_original = 0;
                        // combination 的角色是 index 的位置
                        for(int i = 0; i < enumeration.size(); i++){
                            if(buffer[b].getVal(combination[i]) == enumeration[i]) hit_original++; 
                        }
                        
                        if(hit_original == epi_size) break;
                        if(hit_original < epi_size) b++;    // b 不動 chrom_index chrom_index 往後移動
                    }
                    
                    if(b > (std::pow(2, ell)-1)) break;

                    if(count_chrom[buffer[b].getGene(0)] <= 0){
                        // b is 主角pattern 出現 的 index
                        // b = b_original;
                        b = b + 1;
                        int hit_original = 0; // 這個是用來計算有幾個相同的值
                        while(1)
                        {  
                            if(b > (std::pow(2, ell)-1)) break;

                            hit_original = 0;
                            // combination 的角色是 index 的位置
                            for(int i = 0; i < enumeration.size(); i++){
                                if(buffer[b].getVal(combination[i]) == enumeration[i]) hit_original++; 
                            }
                            
                            if(hit_original == epi_size) break;
                            if(hit_original < epi_size) b++;    // b 不動 chrom_index chrom_index 往後移動
                        }
                        
                        if(b > (std::pow(2, ell)-1)) break;
                        if(count_chrom[buffer[b].getGene(0)] <= 0) break;
                    }


                    // 以當下 enumeration 為主角，去翻每一個 bit >>> 主角的 omega[0] != 去翻每一個omega[0] 這樣就有 epi >>> 就可以 break
                    //// int chrom_index = b+1; // 從fitness第二高的chromosome開始找


                    // chrom_index is 主角pattern flip 一個bit後 omega 的 index
                    int chrom_index = b_original; // 從被sample到的最高fitness的chromosome開始找

                    // 1102 應該不用
                    // if(chrom_index > (std::pow(2, ell)-1)) break;     
                    // while(count_chrom[buffer[chrom_index].getGene(0)] <= 0){
                    //     chrom_index++;
                    //     if(chrom_index > (std::pow(2, ell)-1)) break;
                    // }


                    // if(debug) printf("############## DEBUG ##############\n");
                    
                    int condition_holds = 0;
                    for (int condition_index = 0; condition_index < enumeration.size(); condition_index++)
                    {
                        enumeration = enumeration_original;
                        enumeration[condition_index] = 1 - enumeration[condition_index]; // flip bit 
                        chrom_index = b_original; // 從頭開始找


                        int hit = 0; // 這個是用來計算有幾個相同的值
                        while(1)
                        {  
                            if(chrom_index > (std::pow(2, ell)-1)) break;

                            hit = 0;
                            // combination 的角色是 index 的位置
                            for(int i = 0; i < enumeration.size(); i++){
                                if(buffer[chrom_index].getVal(combination[i]) == enumeration[i]) hit++; 
                            }
                            
                            if(hit == epi_size) break;
                            if(hit < epi_size) chrom_index++;    // b 不動 chrom_index chrom_index 往後移動
                        }
                        
                        if(chrom_index > (std::pow(2, ell)-1)) break;
                        if(count_chrom[buffer[chrom_index].getGene(0)] <= 0) {
                            // b is 主角pattern 出現 的 index
                            // b = b_original;
                            chrom_index = chrom_index + 1;
                            int hit_original = 0; // 這個是用來計算有幾個相同的值
                            while(1)
                            {  
                                if(chrom_index > (std::pow(2, ell)-1)) break;

                                hit_original = 0;
                                // combination 的角色是 index 的位置
                                for(int i = 0; i < enumeration.size(); i++){
                                    if(buffer[chrom_index].getVal(combination[i]) == enumeration[i]) hit_original++; 
                                }
                                
                                if(hit_original == epi_size) break;
                                if(hit_original < epi_size) chrom_index++;    
                            }
                            
                            if(chrom_index > (std::pow(2, ell)-1)) break;
                            if(count_chrom[buffer[chrom_index].getGene(0)] <= 0) break;                            
                        }


                        if((buffer[chrom_index].getVal(0) != buffer[b].getVal(0)) && count_chrom[buffer[chrom_index].getGene(0)] > 0)
                        { 
                            condition_holds += 1;
                        }else{

                            break; // 先不處理fitness相同的情況
                            // Continue searching and put chromosomes with the same fitness as buffer[chrom_index].getVal(0) into a vector
                            std::vector<int> same_fitness_indices;
                            int next_index = chrom_index + 1;
                            
                            // same_fitness_indices.push_back(chrom_index);
                            while (next_index < std::pow(2, ell) && buffer[next_index].getFitness() == buffer[chrom_index].getFitness()) {
                                if (count_chrom[buffer[next_index].getGene(0)] > 0) {
                                    same_fitness_indices.push_back(next_index);
                                }
                                next_index++;
                            }

                            int go_to_next_enumeration = 1;
                            for (int idx : same_fitness_indices) {
                                if (buffer[idx].getVal(0) != buffer[b].getVal(0)) {
                                    condition_holds += 1;
                                    go_to_next_enumeration = 0;
                                    break;
                                }
                            }

                            if (go_to_next_enumeration) {
                                break;
                            };

                            // std::cout << "debug pattern is: ";
                            // for(int p = enumeration.size()-1; p >= 0; p--){
                            //     std::cout << enumeration[p] << " ";
                            // }
                            // std::cout << std::endl;
                            
                            // std::cout << buffer[chrom_index].getFitness() <<std::endl;
                            
                        }
                    
                    }

                    if(condition_holds == epi_size)
                    {
                        one_weak[epi_size]++;
                        weak_epi[epi_size].push_back(combination);


                        
                        if(debug){

                            std::cout << "condition_holds is: "<< condition_holds << std::endl;
                            


                            std::cout << "!!!the pattern is: ";
                            for(int p = enumeration_original.size()-1; p >= 0; p--){
                                std::cout << enumeration_original[p] << " ";
                            }
                            std::cout << std::endl;
                            // for(auto val : enumeration_original){
                            //     std::cout << val << " ";
                            // }
                            // std::cout << std::endl;


                            std::cout << "!!!the fitenss is: ";
                            std::cout << buffer[b].getFitness() <<std::endl;




                            printf("{ ");
                            for (const auto& vertex : combination)
                            {
                                printf("%u ", vertex);
                            }
                            printf("} -> 0\n");
                        }

                        break;
                    }

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
    for(const Chromosome &c: buffer) 
    {
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
    std::cout << "\nbuffer content with count > 0:\n";
    unsigned long bits_best2 = buffer[0].getGene(0);
    for(const Chromosome &c: buffer) {
        if (count_chrom[c.getGene(0)] > 0)
        {
            unsigned long bits = c.getGene(0);
            bits = ~(bits^bits_best2);
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

    }
    
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
