#!/usr/bin/env python

import Bioinfo
import gzip
import matplotlib as mp
import matplotlib.pyplot as plt

file1 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
file2 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
file3 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
file4 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"
'''
fl1 = gzip.open(file1, "rt")
'''
fl2 = gzip.open(file2, "rt")
fl3 = gzip.open(file3, "rt")
'''
fl4 = gzip.open(file4, "rt")
'''
#found using /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | sed -n "2~4p" | wc -l
record_count = 363246735

#read1_score_array = np.full((101, record_count),0, dtype=int)
#read2_score_array = np.full((101, record_count),0, dtype=int)
#index1_score_array = np.full((8, record_count),0, dtype=int)
#index2_score_array = np.full((8, record_count),0, dtype=int)
index1_score_sums = [0 for f in range(8)]
index2_score_sums = [0 for f in range(8)]
'''
read1_score_sums = [0 for f in range(101)]
read2_score_sums = [0 for f in range(101)]
'''
#just need to have the sums of the quality scores, not all the scores themselves
#this block continuously sums of quality scores for each base position until the end of the file
while True:
    '''
    #run through lines in files until quality score is reached (every 4th line)
    read1inp = fl1.readline()
    read1inp = fl1.readline()
    read1inp = fl1.readline()
    #this one will be the quality score
    read1inp = fl1.readline().strip()
    if read1inp == "":
        break
    
    
    read2inp = fl4.readline()
    read2inp = fl4.readline()
    read2inp = fl4.readline()
    #this one will be the quality score
    read2inp = fl4.readline().strip()
    '''
    index1inp = fl2.readline()
    index1inp = fl2.readline()
    index1inp = fl2.readline()
    #this one will be the quality score
    index1inp = fl2.readline().strip()
    
    index2inp = fl3.readline()
    index2inp = fl3.readline()
    index2inp = fl3.readline()
    #this one will be the quality score
    index2inp = fl3.readline().strip()
    if index1inp == "":
        break
    '''
    for h in range(101):
        #read1_score_sums[h] += Bioinfo.convert_phred(read1inp[h])
        read2_score_sums[h] += Bioinfo.convert_phred(read2inp[h])
    
    '''
    
    for i in range(8):
        index1_score_sums[i] += Bioinfo.convert_phred(index1inp[i])
        index2_score_sums[i] += Bioinfo.convert_phred(index2inp[i])  

    '''
'''
#Closing the opened files since we are done with them now

#fl1.close()
fl2.close()
fl3.close()
#fl4.close()
'''
#mean the sums
read1_mean = []
for number in read1_score_sums:
    read1_mean.append(number/record_count)


read2_mean = []
for number in read2_score_sums:
    read2_mean.append(number/record_count)
'''
index1_mean = []
for number in index1_score_sums:
    index1_mean.append(number/record_count)
index2_mean = []
for number in index2_score_sums:
    index2_mean.append(number/record_count)

#create lists of read and index base numbers
read_base_num = [num for num in range(101)]
index_base_num = [num for num in range(8)]

#create a histogram plot for read 1, read 2, index 1, and index 2 base position quality means
'''
plt.rcParams['figure.figsize'] = [16, 6]
plt.bar(x=read_base_num,height=read1_mean, width=0.5)
plt.title("Mean Quality Score of Base Pairs At Each Position of Read 1")
plt.xlabel("# Base Pair")
plt.ylabel("Mean Quality Score")
plt.ylim(0,45)
plt.savefig("/projects/bgmp/wrosales/bioinfo/Bi622/Demultiplex/Assignment-the-first/Read1_hist.png")
plt.clf()

plt.rcParams['figure.figsize'] = [16, 6]
plt.bar(x=read_base_num,height=read2_mean, width=0.5)
plt.title("Mean Quality Score of Base Pairs At Each Position of Read 2")
plt.xlabel("# Base Pair")
plt.ylabel("Mean Quality Score")
plt.ylim(0,45)
plt.savefig("/projects/bgmp/wrosales/bioinfo/Bi622/Demultiplex/Assignment-the-first/Read2_hist.png")
plt.clf()
'''
plt.rcParams['figure.figsize'] = [16, 6]
plt.bar(x=index_base_num,height=index1_mean, width=0.5)
plt.title("Mean Quality Score of Base Pairs At Each Position of Index 1")
plt.xlabel("# Base Pair")
plt.ylabel("Mean Quality Score")
plt.ylim(0,45)
plt.savefig("/projects/bgmp/wrosales/bioinfo/Bi622/Demultiplex/Assignment-the-first/Index1_hist.png")
plt.clf()

plt.rcParams['figure.figsize'] = [16, 6]
plt.bar(x=index_base_num,height=index2_mean, width=0.5)
plt.title("Mean Quality Score of Base Pairs At Each Position of Index 2")
plt.xlabel("# Base Pair")
plt.ylabel("Mean Quality Score")
plt.ylim(0,45)
plt.savefig("/projects/bgmp/wrosales/bioinfo/Bi622/Demultiplex/Assignment-the-first/Index2_hist.png")
