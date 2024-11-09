import sys
import os
import matplotlib.pyplot as plt

if __name__ == '__main__':

    ell = int(sys.argv[1])
    total_time = int(sys.argv[2])
    subtasks = int(sys.argv[3])
    debug = int(sys.argv[4])       
    GHC = int(sys.argv[5])
    sample_size = int(sys.argv[6])


    
    flag = False
    # print(total_weak)

    if GHC == 0:
        ################## GHC 0 ####################

        for sample_size_i in range(1, sample_size+1):
            total_weak = [0]*(ell-1)
            flag = False
            complete_time_GHC_0 = 0
            dir_name = f'../weak_stat/MSO_0_{sample_size_i}/'
            for i in range(subtasks):
                    file_path = f'{dir_name}log_{i}.txt'
                    with open(file_path, 'r') as file:
                        for line in file:
                            count = [int(item) for item in line.split(" ") if item.strip().isdigit()]
                            if (line.split(" ")[-1] == "Complete"):
                                complete_time_GHC_0 += 1 
                            else:
                                print(f"Error processing line: {line}")
                                raise RuntimeError("Error processing line: {line}")      
                            total_weak = [a + b for a, b in zip(total_weak, count)]

            print("complete_time_GHC_0", complete_time_GHC_0)
            dir_name = '../txt_result/' 
            outfile_path = dir_name + f'total_GHC_{GHC}_sample_{sample_size_i}_ell_{ell}.txt'
            with open(outfile_path, 'w') as outfile:
                for i in range(0, ell-1):
                    outfile.write(f'{total_weak[i]} ')

            print("total_time: ", total_time)
   




