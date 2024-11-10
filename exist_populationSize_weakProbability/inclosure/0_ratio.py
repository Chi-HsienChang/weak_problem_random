import sys
import os
import time
import subprocess
from ipdb import set_trace
from tqdm import tqdm
"""
In this python file, the program of running inclosure on cluster is implemented.
:param ell the problem size
:param total_time the total time of random assignment
:param substasks the susbtasks would like to divided into (recommand = 64)
:param GHC how many bits GHC is used (GHC=0 if no GHC)
"""
if __name__ == '__main__':

    ell = int(sys.argv[1])
    runs = int(sys.argv[2])
    GHC = int(sys.argv[3])
    debug = int(sys.argv[4]) 
    seed = int(sys.argv[5]) 
    sample_size = int(sys.argv[6])

    print(f"ell: {ell}, runs: {runs}, GHC: {GHC}, debug: {debug}, seed: {seed}, sample_size: {sample_size}")

    os.system(f'rm -r ../weak_stat/*')
    
    # for s in range(1, sample_size+1):
    for s in tqdm(range(1, sample_size+1)):
        
        if GHC == 0:
            dir_name = f'../weak_stat/MSO_{GHC}_{s}/'
            # 檢查資料夾是否存在
            if not os.path.exists(dir_name):
                # 如果資料夾不存在，則建立資料夾
                os.makedirs(dir_name)
            

            command = f'./inclosure {ell} {runs} {GHC} {debug} -1 {s}'

            # 執行指令並捕獲輸出
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()

            # 將標準輸出轉為字串並存入變數
            result = stdout.decode()

            # 檢查是否有錯誤
            if stderr:
                print("錯誤訊息:", stderr.decode())

            # 輸出結果
            print(f"2^{s}: {result}")   
            time.sleep(1)
     


           