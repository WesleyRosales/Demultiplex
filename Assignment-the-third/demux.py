#!/usr/bin/env python

import argparse
import gzip
import Bioinfo

'''
This script demultiplexed multiplexed sequencing reads. It requires read 1 zipped file, read 2 zipped file, index 1 zipped file, index 2 zipped file, barcode file, and output directory.

Poor quality index threshold is optional and is 20 by default because a default of 20 corresponds to a 99% base call certainty. 
The threshold is in reference to each base call in the indices and if ANY base call is below threshold, the whole read pair gets binned in the poor quality bin.

The output is 2 files (one for read 1 and another for read 2) for each unique dual-matched barcode pair, 2 files for index swapped barcode pairs, and 2 files for poor quality index.

Finally, counts of index swapped, counts of poor quality read-pairs, count of total read-pairs, and pecentage of each dual-matched barcode
'''
#Incorporate argparse parameters to make script useable for general demultiplexing
#7 parameters: read 1 file, read 2 file, index 1 file, index 2 file, barcode file, output directory, minimum quality score (optional)
def get_args():
    my_parser = argparse.ArgumentParser(description='Demultiplex reads')

    my_parser.add_argument('-r1',
                        '--read1',
                        action='store',
                        help='enter read 1 file',
                        required=True)
    
    my_parser.add_argument('-i1',
                        '--index1',
                        action='store',
                        help='enter index 1 file',
                        required=True)

    my_parser.add_argument('-i2',
                        '--index2',
                        action='store',
                        help='enter index 2 file',
                        required=True)

    my_parser.add_argument('-r2',
                        '--read2',
                        action='store',
                        help='enter read 2 file',
                        required=True)

    my_parser.add_argument('-bc',
                        '--barcode',
                        action='store',
                        help='enter barcode file',
                        required=True)

    my_parser.add_argument('-o',
                        '--output',
                        action='store',
                        help='designate output directory',
                        required=True)

    my_parser.add_argument('-q',
                        '--minqual',
                        action='store',
                        type=int,
                        help='poor quality index base call threshold, any reads with an index that has a base call below this number will be pooled in the poor_quality output files',
                        default=20,
                        required=False)
    return my_parser.parse_args()

args = get_args()

#1. open barcodes.txt file and save barcodes in a dictionary with barcode sequence as key and barcode name as value
barcode_dict = {}
with open(args.barcode, "r") as bc:
    header = bc.readline()
    while True:
        line = bc.readline().split()
        if line == []:
            break
        barcode_dict[line[4]] = line[3]

#2. Open each file as read
read1_file = gzip.open(args.read1, "rt")
read2_file = gzip.open(args.read2, "rt")
index1_file = gzip.open(args.index1, "rt")
index2_file = gzip.open(args.index2, "rt")

#3. open all necessary files to write (24 barcode files, 2 index swapped files, 2 poor quality files) as append to continuously add to them
#open index swapped and poor quality files manually, then open barcode files via a for loop
indswap_read1 = open(args.output+"/index_swap_R1.fastq", "w")
indswap_read2 = open(args.output+"/index_swap_R2.fastq", "w")
poorqual_read1 = open(args.output+"/poor_quality_R1.fastq", "w")
poorqual_read2 = open(args.output+"/poor_quality_R2.fastq", "w")
missed_read1 = open(args.output+"/missed_reads_R1.fastq", "w")
missed_read2 = open(args.output+"/missed_reads_R2.fastq", "w")

#create dictionary to hold write file variables as keys and the open function as the values
variable_names = {}
for v in barcode_dict.values():
    variable_names[v + "_read1"] = open(args.output+"/"+v+"_read1.fastq", "w")
    variable_names[v + "_read2"] = open(args.output+"/"+v+"_read2.fastq", "w")

#in preparation for next step, make a reverse complement function
def rev_comp(DNA: str) -> str:
    '''Takes any string of DNA (including "N") and outputs the reverse complement'''
    rev_comp = ""
    #make rev_comp dict with keys being DNA and values being reverse complement
    rev_comp_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    lenDNA = len(DNA)
    for i in reversed(range(lenDNA)):
        #work through DNA string in reverse and append that to 
        rev_comp += rev_comp_dict[DNA[i]]
    return rev_comp

#initialize counts to keep track of read-pair metrics
total_count = 0
barcode_count = {bc: 0  for bc in barcode_dict.values()}
poorqual_count = 0
indswap_count = 0
missed_reads = 0
#4. IN A WHILE LOOP:
while True:        
    #for each record, save each line as a variables.
    read1_head = read1_file.readline().strip()
    read1_seq = read1_file.readline().strip()
    read1_plus = read1_file.readline().strip()
    read1_qual = read1_file.readline().strip()

    read2_head = read2_file.readline().strip()
    read2_seq = read2_file.readline().strip()
    read2_plus = read2_file.readline().strip()
    read2_qual = read2_file.readline().strip()

    index1_head = index1_file.readline().strip()
    index1_seq = index1_file.readline().strip()
    index1_plus = index1_file.readline().strip()
    index1_qual = index1_file.readline().strip()

    index2_head = index2_file.readline().strip()
    index2_seq = index2_file.readline().strip()
    index2_plus = index2_file.readline().strip()
    index2_qual = index2_file.readline().strip()
    #reverse complement index 2 sequences
    index2_seq_rc = rev_comp(index2_seq)

    #check for "" in any variables to indicate end of file
    if read1_head == "":
        break

    #Create lists with boolean for qualisty scores in index 1 and 2 being less than quality threshold
    #If any score is < entered threshold (default 20), record goes into poor quality bin
    #This is a conservative filter because if any index is ambiguous, we cannot confidently tell if it is the index it is called as
    index1_scores = [Bioinfo.convert_phred(x) < args.minqual for x in index1_qual]
    index2_scores = [Bioinfo.convert_phred(x) < args.minqual for x in index2_qual]

    #if index 1 or index 2 have qscore < 20 OR they have "N" in sequence OR they are not in the list of barcodes:
        #make a fastq record for read 1 and read 2 with header having "index1-rev comp index2" appended to it
        #add that record to appropriate poor quality file to be written
    if "N" in index1_seq or "N" in index2_seq or any(index1_scores) or any(index2_scores) or index1_seq not in barcode_dict.keys() or index2_seq_rc not in barcode_dict.keys():
        poorqual_count += 1
        read1_new_head = read1_head + " " + index1_seq + "-" + index2_seq_rc
        read2_new_head = read2_head + " " + index1_seq + "-" + index2_seq_rc

        poorqual_read1.writelines(read1_new_head + "\n")
        poorqual_read1.writelines(read1_seq + "\n")
        poorqual_read1.writelines(read1_plus + "\n")
        poorqual_read1.writelines(read1_qual + "\n")

        poorqual_read2.writelines(read2_new_head + "\n")
        poorqual_read2.writelines(read2_seq + "\n")
        poorqual_read2.writelines(read2_plus + "\n")
        poorqual_read2.writelines(read2_qual + "\n")

    #elif index 1 != rev comp index 2 AND index 1 and rev comp index 2 are in barcodes dictionary:
        #make a fastq record for read 1 and read 2 with header having "index1-rev comp index2" appended to it
        #add that record to appropriate swapped file to be written
    elif index1_seq != index2_seq_rc and index1_seq in barcode_dict.keys() and index2_seq_rc in barcode_dict.keys():
        indswap_count += 1
        read1_new_head = read1_head + " " + index1_seq + "-" + index2_seq_rc
        read2_new_head = read2_head + " " + index1_seq + "-" + index2_seq_rc

        indswap_read1.writelines(read1_new_head + "\n")
        indswap_read1.writelines(read1_seq + "\n")
        indswap_read1.writelines(read1_plus + "\n")
        indswap_read1.writelines(read1_qual + "\n")

        indswap_read2.writelines(read2_new_head + "\n")
        indswap_read2.writelines(read2_seq + "\n")
        indswap_read2.writelines(read2_plus + "\n")
        indswap_read2.writelines(read2_qual + "\n")

    #if index 1 == rev comp index 2 AND index 1 is barcodes dictionary:
        #make a fastq record for read 1 and read 2 with header having "index1-rev comp index2" appended to it
        #add that record to file to be written in proper barcode file
    elif index1_seq == index2_seq_rc and index1_seq in barcode_dict.keys():
        barcode_count[barcode_dict[index1_seq]] += 1
        read1_new_head = read1_head + " " + index1_seq + "-" + index2_seq_rc
        read2_new_head = read2_head + " " + index1_seq + "-" + index2_seq_rc
        bcname1 = barcode_dict[index1_seq] + "_read1"
        bcname2 = barcode_dict[index1_seq] + "_read2"

        variable_names[bcname1].writelines(read1_new_head + "\n")
        variable_names[bcname1].writelines(read1_seq + "\n")
        variable_names[bcname1].writelines(read1_plus + "\n")
        variable_names[bcname1].writelines(read1_qual + "\n")

        variable_names[bcname2].writelines(read2_new_head + "\n")
        variable_names[bcname2].writelines(read2_seq + "\n")
        variable_names[bcname2].writelines(read2_plus + "\n")
        variable_names[bcname2].writelines(read2_qual + "\n")

    #This is a catch-all for any reads that I missed. If everything works properly, it should be empty.
    else:
        missed_reads += 1
        read1_new_head = read1_head + " " + index1_seq + "-" + index2_seq_rc
        read2_new_head = read2_head + " " + index1_seq + "-" + index2_seq_rc

        missed_read1.writelines(read1_new_head + "\n")
        missed_read1.writelines(read1_seq + "\n")
        missed_read1.writelines(read1_plus + "\n")
        missed_read1.writelines(read1_qual + "\n")

        missed_read2.writelines(read2_new_head + "\n")
        missed_read2.writelines(read2_seq + "\n")
        missed_read2.writelines(read2_plus + "\n")
        missed_read2.writelines(read2_qual + "\n")        



#close all opened read and write files
read1_file.close()
read2_file.close()
index1_file.close()
index2_file.close()

indswap_read1.close()
indswap_read2.close()
poorqual_read1.close()
poorqual_read2.close()
missed_read1.close()
missed_read2.close()

for v in variable_names.keys():
    variable_names[v].close()

#output demux metrics into txt file
with open(args.output+"/demux_stats.txt", "w") as dout:
    total_read_counts = sum(barcode_count.values())
    total_read_counts += indswap_count
    total_read_counts += poorqual_count
    dout.writelines("Total read-pairs: " + str(total_read_counts) + "\n")
    for k, v in barcode_count.items():
        dout.writelines("Percentage "+ str(k) + " dual-matched barcodes: " + str(v*100/total_read_counts) + "% of total read-pairs (raw count: " + str(v) + "\n")
    dout.writelines("Index swapped read-pairs: " + str(indswap_count) + "\n")
    dout.writelines("Poor quality index read-pairs: " + str(poorqual_count) + "\n")
    dout.writelines("Missed read-pairs: " + str(missed_reads))

