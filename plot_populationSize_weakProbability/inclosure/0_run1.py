import sys
import os
import time
from ipdb import set_trace
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
    debug = int(sys.argv[4])       
    GHC = int(sys.argv[5])
    sample_size = int(sys.argv[6])

    if(total_time%subtasks != 0):
        print("Please check total time is divisable by subtasks.")
        
    os.system(f'rm -r ../weak_stat/*')

    for s in range(1, sample_size+1):
        one_time = int(total_time/subtasks)
        if GHC == 0:
            dir_name = f'../weak_stat/MSO_{GHC}_{s}/'
            # 檢查資料夾是否存在
            if not os.path.exists(dir_name):
                # 如果資料夾不存在，則建立資料夾
                os.makedirs(dir_name)
            
            for i in range(subtasks):
                # cmd = "cd /home/salima/MSO_GHC/inclosure"
                # os.system(cmd)
                cmd = f'nohup ./inclosure {ell} {one_time} {GHC} {debug} {time.time()+i} {s} >> ' + dir_name + f'log_{i}.txt 2>&1 &'
                print("total_time: ", total_time)
                time.sleep(1)
                print(cmd)
                os.system(cmd)


           