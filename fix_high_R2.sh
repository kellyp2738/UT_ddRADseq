#!/bin/bash

highR2=/scratch/02540/kellypie/Rsq/high*

for file in $highR2
do 
	#echo $file
	out1=$file.fixed_temp
	out2=$file.fixed_temp2
	out3=$file.fixed
	#echo $out1, $out2

	sed -n '1~2!p' $file > $out1
	paste -s -d' \n' $out1 > $out2
	awk '{print $2,$3,$4,$5,$7}' $out2 > $out3
	rm $out1
	rm $out2
done
