#!python -OO

import numpy.random as rand
import numpy as np
from subprocess import call, Popen, PIPE
import subprocess
import os, os.path
import sys
import io
from scipy.stats import chisqprob
from operator import itemgetter
import csv

# line for profiling system usage
# python -m cProfile sample_IDs.py -min 10 -max 10 -r 1 -i /home/antolinlab/Desktop/IDs.txt -v /home/antolinlab/Desktop/D_variabilis_Pseudoref/MasterPseudoRefVCF_Copy/pseudoref_mapped_genotypes.vcf -os /home/antolinlab/Desktop/snp_bootstrap_test2.txt -or /home/antolinlab/Desktop/UT_ddRADseq/Rsq/high_R2_fraction.csv -m 0.75

if __name__ == '__main__':
    
    import argparse

    parser = argparse.ArgumentParser(description = 'Generates "keep_files" for boostrapping tick SNP data.')
    parser.add_argument('-min', '--min_sample', help = 'Minimum number of ticks to sample.')
    parser.add_argument('-max', '--max_sample', help = 'Maximum number of ticks to sample.')
    parser.add_argument('-r', '--repetitions', help = 'Number of repetitions for each sample size.')    
    parser.add_argument('-i', '--input_ID_file', help = 'Path to input ID file.')
    parser.add_argument('-v', '--input_VCF_file', help = 'Path to input VCF file.')
    parser.add_argument('-os', '--output_snp', help = 'Path to output SNP count file.')
    parser.add_argument('-or', '--output_R2', help = 'Path to output fraction high R2 file.')
    parser.add_argument('-m', '--max_missing', help = 'Max-missing parameter for VCFtools.')

    opts = parser.parse_args()



min_ticks=int(opts.min_sample)
max_ticks=int(opts.max_sample)
r=int(opts.repetitions)
infile_ID = opts.input_ID_file
infile_VCF = opts.input_VCF_file
outfile1 = opts.output_snp
outfile2 = opts.output_R2
max_missing = opts.max_missing



# delete the pre-existing output file?
if (os.path.isfile(outfile1) or os.path.isfile(outfile2)):
	print outfile1
	print outfile2
	q = 'Overwrite existing data files named above?'
	prompt = '[Y/n]'
	valid = {"yes":True, "y":True, "Y":True, "Yes":True}
	sys.stdout.write(q + prompt)
	choice = raw_input().lower()
	if choice in valid:
		if os.path.isfile(outfile1):
			os.remove(outfile1)
		if os.path.isfile(outfile2):
			os.remove(outfile2)
	else:
		print 'Please choose a different file name.'
		quit()
	
for n in range(min_ticks, (max_ticks+1)):
    for j in range(1, (r+1)):
	print 'repetition number', j        
	# open a file to which sample data can be appended
        final_data=open(outfile1, 'a')

	# open the file of tick IDs that will be sampled
	id_file=open(infile_ID, 'r')

        # permute the tick IDs for resampling/refiltering the VCF file
        id_names=id_file.readlines() #this creates a vector that has sample names and newline characters
        id_permute=rand.permutation(id_names) # permute the array randomly to generate new combos of ticks to sample
        id_sample=id_permute[0:n] # get the first n elements of the array
	print id_sample

        # create a temporary file for holding the tick sample IDs
        f=open('/home/antolinlab/Desktop/tmp.txt', 'w') # open a temp file
        for i in range(0, n):
	    print 'i', i
            print id_sample[i]
            f.writelines(id_sample[i]) # add the sample name to the file. newlines are already embedded in the name.
        f.close()

        # call VCFtools with tmp.txt used as the "--keep" file
        my_env=os.environ.copy()
        # the arguments need to be fully separated for this to work (["--vcf", "/file/path/"], not ["--vcf /file/path/"])
        p = Popen(["vcftools", "--vcf", infile_VCF, "--keep", "/home/antolinlab/Desktop/tmp.txt", "--min-meanDP", "20", "--minGQ", "25", "--maf", "0.05", "--max-missing", max_missing, "--recode", "--out", '/home/antolinlab/Desktop/vcf_tmp'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=my_env)

        # extract the number of SNPs retained
        stdout, stderr = p.communicate()
        split_stderr = str.splitlines(stderr)
        snp_line = split_stderr[17]
        snp_count = [int(s) for s in snp_line.split() if s.isdigit()][0] # split snp_line, search for integers, and save the first one found

	'''
        vcf_call = "vcftools --vcf /home/antolinlab/Desktop/vcf_chrom_rename_final.vcf --plink --out /home/antolinlab/Desktop/vcf_tmp_plink"
        processList = []
        for i in xrange(0,10):
            proc = Popen(args(i))
            processList.append(proc)
        for proc in processList:
            stdout, stderr = proc.communicate()
            with file as open(file,'w'):
            	file.write(stdout + os.endl)
	'''

        # Plink expects human data and doesn't like chromosome numbers > 22. We don't know the chromosome structure, so replace all the chromosome names with "1"
        # first get rid of the text in the fragment ID
        subprocess.call("sed 's|_pseudoreference_pe_concatenated_without_rev_complement||g' /home/antolinlab/Desktop/vcf_tmp.recode.vcf > /home/antolinlab/Desktop/vcf_chrom_rename.vcf", shell=True)
        # second, find the fragment number and replace it with a 1
        subprocess.call("sed -r 's/^[0-9]+/1/' /home/antolinlab/Desktop/vcf_chrom_rename.vcf > /home/antolinlab/Desktop/vcf_chrom_rename_2.vcf", shell=True)
        # now the SNP IDs are not unique... fix it with more sed and some R
        subprocess.call("sed 's|#CHROM|CHROM|' /home/antolinlab/Desktop/vcf_chrom_rename_2.vcf > /home/antolinlab/Desktop/vcf_chrom_rename_3.vcf", shell=True)
        subprocess.call("Rscript fix_vcf_pos.r /home/antolinlab/Desktop/vcf_chrom_rename_3.vcf /home/antolinlab/Desktop/vcf_chrom_rename_final.vcf", shell = True) # this will replace the "POS" column in the VCF file with consecutive numbers

        # make a "snps" file with 500 unique SNPs for LD stats

        # 1. get SNP IDs from vcf file
        # convert the VCF file to a PLINK file
        subprocess.call("vcftools --vcf /home/antolinlab/Desktop/vcf_chrom_rename_final.vcf --plink --out /home/antolinlab/Desktop/vcf_tmp_plink", shell=True)

        # calculate LD on the plink file
        subprocess.call("plink --file /home/antolinlab/Desktop/vcf_tmp_plink --r2 --matrix --noweb", shell=True) # option A: pairwise matrix
        #subprocess.call("plink --file /home/antolinlab/Desktop/vcf_tmp_plink --r2 --inter-chr --allow-no-sex --ld-window-r2 0 --noweb", shell=True) # option B: long-format list

	'''
	# Old scripts for resampling LD files and inputting into q-value calculation
	# the q-value calculation was never fully implemented in R...

        # for option B (long form) only: PLINK long-form files have extra spaces to make them human readable... remove those spaces for R import
        #subprocess.call("sed -E 's|  *| |g' plink.ld > plink_edit-2.ld", shell=True)

        # for option A (matrix), the plink file is often too big to be handled by R.
        # get the upper triangle from the matrix, convert to a vector, and take random subsamples of the data
        # convert to long form without loading the whole matrix into memory... the resulting long form should fit in memory

	# some explanation using a toy symmetric matrix
	# data...	index w/respect to full upper triangle... 	actual row indexes... 	row indexes when lower triangle and diagonal are removed
	#[[1, x, x],	[[-, 0, 1],					[[0, 1, 2],		[0,1]
	# [x, 1, x],	 [-, -, 2],					 [0, 1, 2],		[0]
	# [x, x, 1]]	 [-, -, -]]					 [0, 1, 2]]
	
	# the code chunk below relates the index with respect to the full upper triangle 
	# (corresponding to data we want) with the indexes assigned each row as it is read into memory
        
	save_name = '/home/antolinlab/Desktop/' + str(n) + "_" + 'retained_plink_R2.csv'   
	# if we want to overwrite the output file each time because we don't need it...
	os.remove(save_name)	
	print 'Resampling LD stats'
	for ldResample in range(10):
		idxSave = np.random.choice(range(snp_count), size=30, replace=False)
		i = 1 # start line counter
		s = 0 # start SNP counter
		with open('plink.ld') as f:		
			for line in f: # read file line by line
				rowData = (line.split())[i:snp_count] # get the upper triangle
				trueIdx = range(s, s+len(rowData)-1) # get the index w/respect to complete upper triangle (not just the current row loaded)					keepTrueIdx = list(set(trueIdx) & set(idxSave)) # find the intersection between what we have loaded and what we'd like to keep
				keepTrueIdx = list(set(trueIdx) & set(idxSave)) # find the intersection between what we have loaded and what we'd like to keep
				i += 1 # increase the line counter
				s = s + len(rowData) # increase the SNP counter
				if len(keepTrueIdx) > 0:
					for ti in keepTrueIdx:				
						keepRowIdx = trueIdx.index(ti)	# the index of trueIdx will correspond to the index of rowData. look up the index of the keepTrueIndexes (not very pythonic, but so goes it)			
						saveR2 = rowData[keepRowIdx]
						write_str = str(ldResample) + ',' + str(saveR2)					
						with open(save_name, 'a') as c:
							ldWriter = csv.writer(c, delimiter=',')
							ldWriter.writerow([str(ldResample)] + [str(saveR2)])
	'''	

	r = Popen(["Rscript", "save_R2_hist.r", str(n)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
	outs, errs = r.communicate()
	out_split = str.splitlines(outs) #out contains the print output from r (fraction 'significant', which in this context is simply greater than some R^2 threshold
	outss = str.split(out_split[0])
	outss[0] = n # replace the R print output line number, which we don't need anyway, with the sample size	
	with open(outfile2, 'a') as c: # open in append mode
		r2Writer = csv.writer(c, delimiter=",")
		r2Writer.writerow(outss)	

# save the sample size and snp count
write_line = [str(n), ",", str(snp_count), "\n"]
final_data.writelines(write_line)
final_data.close()



