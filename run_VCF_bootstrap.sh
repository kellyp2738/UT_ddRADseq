#!/bin/bash

# positional arguments used to define extent of resampling
# $1 = minimum number of ticks to sample
# $2 = maximum number of ticks to sample
# $3 = number of repetitions of each sample size

# positional arguments to the bash script passed to python:
# $4 = -i; input ID file to permute
# $5 = -v; input VCF file to sample
# $6 = -o; output data file
# $7 = -m; max_missing for VCFtools

# check that output file does not already exist
# if it does exist, should it be overwritten?

if [ -f $6 ];
then
	#echo 	
	read -p ["Overwrite existing output file $6 ?"] -r
	echo    # move to a new line
	if [[ $REPLY =~ ^[Yy]$ ]];
	then
		echo "Removing $6."    		
		rm $6
	else
		echo "Please chooose a different output file name."
		exit 0
	fi
fi

echo "Sampling $1 to $2 individuals for $3 iterations and writing output to $6."
for j in `seq $1 $2`; # for $n between $1 and $2
do
	echo "Sampling $j individuals."
	for k in `seq 1 $3`; # sample $n $3 times
	do
		#echo $j
		python sample_IDs.py -n $j -i $4 -v $5 -o $6 -m $7
	done
done


