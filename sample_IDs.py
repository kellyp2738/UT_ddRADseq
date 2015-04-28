#!/bin/bash/python

import numpy.random as rand
import numpy as np
from subprocess import call, Popen, PIPE
import subprocess


if __name__ == '__main__':
    
    import argparse
    
    parser = argparse.ArgumentParser(description = 'Generates "keep_files" for boostrapping tick SNP data.')
    parser.add_argument('-n', '--subsample_size', help = 'Number of subsamples of tick data set')
 
    opts = parser.parse_args()

id=open('/home/antolinlab/Desktop/IDs.txt', 'r')

n=int(opts.subsample_size) # total number of ticks

id_names=id.readlines() #this creates a vector that has sample names and newline characters
id_permute=rand.permutation(id_names) # permute the array randomly to generate new combos of ticks to sample
id_sample=id_permute[0:n] # get the first n elements of the array

f=open('/home/antolinlab/Desktop/tmp.txt', 'w') # open a temp file
for i in range(0, n):
	#print id_sample[i]
	f.writelines(id_sample[i]) # add the sample name to the file. newlines are already embedded in the name.
f.close()

call(["vcftools", "--vcf ~/Desktop/D_variabilis_Pseudoref/MasterPseudoRefVCF_Copy/pseudoref_mapped_genotypes.vcf", "--keep ~/Desktop/tmp.txt", "--max-missing 0.75", "--out ~/Desktop/vcf_tmp"], shell=True, )

## this works!!!
my_env=os.environ.copy()
p = Popen(["vcftools", "--vcf", "/home/antolinlab/Desktop/D_variabilis_Pseudoref/MasterPseudoRefVCF_Copy/pseudoref_mapped_genotypes.vcf", "--keep", "/home/antolinlab/Desktop/tmp.txt", "--max-missing", "0.75", "--out", "/home/antolinlab/Desktop/vcf_tmp"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=my_env)
stdout, stderr = p.communicate()
