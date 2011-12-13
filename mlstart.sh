addpath(pwd)
cd core/
addpath(pwd)
cd ../exp
addpath(pwd)
cd ../solve
addpath(pwd)
cd ../misc
addpath(pwd)
cd ../tests

maxNumCompThreads(2)

test_steps3

cd ../
save mlsave.mat
