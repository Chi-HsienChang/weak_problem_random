import sys
import os
import time
from ipdb import set_trace
import time

"""
In this python file, the program of running inclosure on cluster is implemented.
:param ell the problem size
:param total_time the total time of random assignment
:param substasks the susbtasks would like to divided into (recommand = 64)
:param GHC how many bits GHC is used (GHC=0 if no GHC)
"""
if __name__ == '__main__':

    ell = int(sys.argv[1])
    total_time = int(sys.argv[2])
    subtasks = int(sys.argv[3])
    nInitial = int(sys.argv[4])       
    selectionPressure = int(sys.argv[5])
    pc = float(sys.argv[6])
    pm = float(sys.argv[7])
    maxGen = int(sys.argv[8])
    maxFe = int(sys.argv[9])

    if(total_time%subtasks != 0):
        print("Please check total time is divisable by subtasks.")
          
    one_time_repeat = int(total_time/subtasks)
    os.system(f'rm ./output_log_SGA_result/*')
    os.system(f'rm ./SGA_result/* -r')

    for s in range(0, subtasks):
        os.system(f'rm ./SGA_result/syn_{s}/*')
        dir_name = f'./SGA_result/syn_{s}/'
        os.system("mkdir -p "+dir_name)

        # cmd = f'nohup ./SGA {ell} 200 2 0.8 0.01 100 1000000 {one_time}'
        # cmd = f'nohup ./SGA {ell} 200 2 0.8 0.01 100 1000000 {one_time} {s}>> ' + dir_name + f'log_{s}.txt 2>&1 &'
        
        dir_name_log = f'./output_log_SGA_result/'
        os.system("mkdir -p "+dir_name_log)
        cmd = f'nohup ./SGA {ell} {nInitial} {selectionPressure} {pc} {pm} {maxGen} {maxFe} {one_time_repeat} {s} >> ' + dir_name_log + f'log_subtask_{s}.txt 2>&1 &'

        time.sleep(1)
        print("total_time: ", total_time)
        print(cmd)
        os.system(cmd)



        