#!/usr/bin/python

import numpy.random as rand
import numpy as np
from subprocess import call, Popen, PIPE
import subprocess
import os, os.path
import sys


if __name__ == '__main__':
    
	import argparse

	parser = argparse.ArgumentParser(description = 'Generates "keep_files" for boostrapping tick SNP data.')
	parser.add_argument('-n', '--subsample_size', help = 'Number of subsamples of tick data set')
	parser.add_argument('-i', '--input_ID_file', help = 'Path to input ID file.')
	parser.add_argument('-v', '--input_VCF_file', help = 'Path to input VCF file.')
	parser.add_argument('-o', '--output', help = 'Path to output SNP count file.')
	parser.add_argument('-m', '--max_missing', help = 'Max-missing parameter for VCFtools.')

	opts = parser.parse_args()



n=int(opts.subsample_size) # total number of ticks
infile_ID = opts.input_ID_file
infile_VCF = opts.input_VCF_file
outfile = opts.output
max_missing = opts.max_missing

# open the file of tick IDs that will be sampled
id_file=open(infile_ID, 'r')

# delete the pre-existing output file?
#if os.path.isfile(outfile):
#	q = 'Overwrite existing data file?'
#	prompt = '[Y/n]'
#	valid = {"yes":True, "y":True, "Y":True, "Yes":True}
#	sys.stdout.write(q + prompt)
#	choice = raw_input().lower()
#	if choice in valid:
#		os.remove(outfile)
#	else:
#		print 'Please choose a different file name.'
#		quit()

# open a file to which sample data can be appended
final_data=open(outfile, 'a')

# permute the tick IDs for resampling/refiltering the VCF file
id_names=id_file.readlines() #this creates a vector that has sample names and newline characters
id_permute=rand.permutation(id_names) # permute the array randomly to generate new combos of ticks to sample
id_sample=id_permute[0:n] # get the first n elements of the array

# create a temporary file for holding the tick sample IDs
f=open('/home/antolinlab/Desktop/tmp.txt', 'w') # open a temp file
for i in range(0, n):
	#print id_sample[i]
	f.writelines(id_sample[i]) # add the sample name to the file. newlines are already embedded in the name.
f.close()

# call VCFtools with tmp.txt used as the "--keep" file
my_env=os.environ.copy()
# the arguments need to be fully separated for this to work (["--vcf", "/file/path/"], not ["--vcf /file/path/"])
p = Popen(["vcftools", "--vcf", infile_VCF, "--keep", "/home/antolinlab/Desktop/tmp.txt", "--min-meanDP", "20", "--minGQ", "25", "--maf", "0.05", "--max-missing", max_missing, "--recode", "--out", '/home/antolinlab/Desktop/vcf_tmp'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=my_env)

# extract the number of SNPs retained
stdout, stderr = p.communicate()
split_stderr = str.splitlines(stderr)
snp_line = split_stderr[15]
snp_count = [int(s) for s in snp_line.split() if s.isdigit()][0] # split snp_line, search for integers, and save the first one found

# convert the VCF file to a PLINK file
subprocess.call("vcftools --vcf /home/antolinlab/Desktop/vcf_tmp.recode.vcf --plink --out /home/antolinlab/Desktop/vcf_tmp_plink", shell=True)

# calculate LD on the plink file
# subprocess.call("plink --file vcf_tmp_plink --matrix", shell=True) # option A: pairwise matrix
#subprocess.call("plink --file /home/antolinlab/Desktop/vcf_tmp_plink --ld --inter-chr --allow-no-sex --noweb", shell=True) # option B: long-format list
subprocess.call("plink --file /home/antolinlab/Desktop/vcf_tmp_plink --r2 --inter-chr --noweb", shell=True) # option B: long-format list

# PLINK files have extra spaces to make them human readable... remove those spaces for R import
subprocess.call("sed -E 's|  *| |g' plink.ld > plink_edit-2.ld", shell=True)
r = Popen(["Rscript", "plink_LD_sig.r", "plink_edit-2.ld"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
outs, errs = r.communicate()
#split_stderr_r = str.splitlines(stderr_r)
print outs, errs
snp_pairs = split_stderr_r[1]
sig_pairs = split_stderr_r[2]

# save the sample size and snp count
write_line = [str(n), ",", str(snp_count), ",", str(snp_pairs), ",", str(sig_pairs), "\n"]
final_data.writelines(write_line)
final_data.close()



