#/bin/bash
ell=30
numConvergence=3
fffff=1

./sweep $ell $numConvergence $fffff | tee "./sweep_log/${ell}_${numConvergence}_${fffff}.txt"