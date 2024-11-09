#/bin/bash
rand_seed=2024
# usage: ./MSO problem_size initial_population benchmark max_generation max_Function_evaluation repeat display $rand_seed
# usage: ./sweep problem_size numConvergence benchmark

make 
#ELL=("10" "20" "30")
ELL=("100")
truncate -s 0 fe.txt
file_path="./log/"

for ell in "${ELL[@]}"; do
    num=$(( $ell / 5 ))
    size=$(( $num + 500 ))
    ./MSO $ell $size 1 100 1000 5 1 $rand_seed
    echo "Ell: $ell done"
    echo -n $ell" " >> fe.txt
    tail -n 1 log.txt | sed -n 's/[^0-9]//gp' >> fe.txt

    file_path_="${file_path}${ell}_log.txt"
    if [ ! -e "$file_path_" ]; then
        cp log.txt "$file_path_"
    fi

done
gprof ./MSO gmon.out | gprof2dot -w | dot -Tpng -o profile.png

python plot.py
#./sweep 30 100 1