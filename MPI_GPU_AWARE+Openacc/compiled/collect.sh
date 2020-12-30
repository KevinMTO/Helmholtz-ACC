#!/bin/bash


for file in ./*.txt
do
	name=${file##*/}
	base=${name%.txt}
	grep MFlops $file >> "${base}_Mflops.txt"
	grep Elapsed $file >> "${base}_Elapsed.txt"

done
