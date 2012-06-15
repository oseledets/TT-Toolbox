#!/bin/bash

if [ "$#" -lt "3" ]; then
	echo "Usage: genmexmakefile.sh source_file matlab_command compiler_command additional_compiler_flags"
	exit;
fi

fullname="$1.c";
if [ ! -f "$fullname" ]; then
	fullname="$1.cpp";
fi
if [ ! -f "$fullname" ]; then
	fullname="$1.cxx";
fi
if [ ! -f "$fullname" ]; then
	fullname="$1.cc";
fi
if [ ! -f "$fullname" ]; then
	fullname="$1.C";
fi
if [ ! -f "$fullname" ]; then
	fullname="$1.CC";
fi
if [ ! -f "$fullname" ]; then
	fullname="$1.CPP";
fi
if [ ! -f "$fullname" ]; then
	fullname="$1.F";
fi
if [ ! -f "$fullname" ]; then
	fullname="$1.f90";
fi
if [ ! -f "$fullname" ]; then
	fullname="$1.F90";
fi
if [ ! -f "$fullname" ]; then
	fullname="$1.f";
fi

echo "fprintf('');" > mexrun.log
echo "mex -n -largeArrayDims -lmwlapack -lmwblas $fullname" >> mexrun.log

$2 -nodesktop -nosplash -nodisplay < mexrun.log -logfile mexcommands.log

grep "$1" mexcommands.log > mexrun.log
sed '/^$/d' mexrun.log > mexcommands.log

curcc=`awk '{print $2}' mexcommands.log | tail -n 1`
garbage=`awk '{print $1}' mexcommands.log | tail -n 1`

sed "s/$garbage/\t/g" mexcommands.log > mexrun.log
sed "s/$curcc/$3 $4/g" mexrun.log > mexcommands.log

sed "s/-ansi//g" mexcommands.log > mexrun.log
sed "s/-DNDEBUG//g" mexrun.log > mexcommands.log
sed "s/-O/-O3/g" mexcommands.log > mexrun.log

cat mexrun.log

echo "all:" > makefile
cat mexrun.log >> makefile

rm -f mexrun.log
rm -f mexcommands.log
