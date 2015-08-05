#!python -OO

import numpy.random as rand
import numpy as np
from operator import itemgetter
import csv

save_name = '/home/antolinlab/Desktop/' + str(n) + "_" + 'retained_plink_R2.csv'   
print 'Resampling LD stats'
for ldResample in range(10):
	snp_count = 300
	n = 10
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

