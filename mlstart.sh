addpath(pwd)
cd core/
addpath(pwd)
cd ../exp
addpath(pwd)
cd ../solve
addpath(pwd)
cd ../tests

maxNumCompThreads(2)

test_steps

cd ../
save mlsave.mat
