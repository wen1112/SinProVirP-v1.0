import csv
import sys
import numpy as np
#import re
from Bio import SeqIO

bam_stat = sys.argv[1]
protein_VC = sys.argv[2]
marker_ratio = sys.argv[3]
out = sys.argv[4]

bamstatfile = open(bam_stat, "rt")
proteinVCfile = open(protein_VC, "rt")
outfile = open(out, "w")
marker_ratio = float(marker_ratio)


##step1 protein --> VC ID
proteinid_VC = {}
vc_candimarkernum = {}
bamstatfile_head = bamstatfile.readline()

for line in proteinVCfile:
    line = line.rstrip('\n')
    array = line.split("\t")
    VC = array[5]
    proteinid = array[0]
    proteinid_VC[proteinid] = VC
    vc_candimarkernum[VC] = array[17]



##step2 VC --> markers stats

vc_markerlist = {}
unmapped_read_num = 0
mapped_reads_num = 0
vc_marker_readsnum = {}
for line in bamstatfile:
    line = line.rstrip('\n')
    array = line.split("\t")
    proteinid = array[0]
    proteinlength = int(array[1])
    protein_mapped_read_num = int(array[2])
    if array[0] == "*":
        unmapped_read_num = array[3]
    else:
        relative_reads_num = protein_mapped_read_num/proteinlength
        mapped_reads_num += protein_mapped_read_num
        if proteinid in proteinid_VC:
            vcid = proteinid_VC[proteinid]
            vc_markerlist.setdefault(vcid,{})[proteinid] = relative_reads_num
            vc_marker_readsnum.setdefault(vcid,{})[proteinid] = protein_mapped_read_num
        else:
            print(f'''Error:{proteinid} not exists!''')

##
total_reads_num = int(unmapped_read_num) + int(mapped_reads_num)
if total_reads_num != 0:
    mapping_ratio = round(mapped_reads_num/total_reads_num,4)
    print(f'''Reads mapping ratio: {mapping_ratio}''')
else:
    print("Error: Total reads num is 0!")


sum_vc_mean_abun = 0
sum_vc_reads_num = 0

for key,value in vc_markerlist.items():
    vcid = key
    ## VC marker stat
    vc_markerprotein_num = int(len(value.keys()))
   
    ## VC abudancne stat 
    markerprotein_ids = value.keys()
    markerprotein_relab = list(value.values())
    markerprotein_relab = sorted(markerprotein_relab,key=float)
    markerprotein_relab_valued = sum(relab>0 for relab in markerprotein_relab)
    marker_covered_ratio = float(markerprotein_relab_valued/len(markerprotein_relab))
    if marker_covered_ratio >= marker_ratio:

    ##VC reads counts
        markerprotein_reads = list(vc_marker_readsnum[vcid].values())
        markerprotein_reads = sorted(markerprotein_reads,key=int)

        low_index = int(vc_markerprotein_num * 0.1)
        up_index = int(vc_markerprotein_num * 0.9)+1
        if up_index>low_index:
            subset_abun = markerprotein_relab[low_index:up_index]  
            mean_abundance = np.mean(subset_abun)
            subset_reads = markerprotein_reads[low_index:up_index]
            sum_reads_num = np.sum(subset_reads)
            sum_vc_mean_abun  += mean_abundance
            sum_vc_reads_num  += sum_reads_num
        else:
            continue
    


outfile.write("VCID\tVC_markerprotein_num\tVC_markerprotein_valued\tSelected_marker_num\tSelected_marker_withabundance\tMean_abundance_measure\tRelative_abudannce\tTotal_reads_num\tRelative_reads_num\n")
for key,value in vc_markerlist.items():
    vcid = key
    ## VC marker stat
    vc_markerprotein_num = int(len(value.keys()))

    ## VC abudancne stat
    markerprotein_ids = value.keys()
    markerprotein_relab = list(value.values())
    markerprotein_relab = sorted(markerprotein_relab,key=float)
    markerprotein_relab_valued = sum(relab>0 for relab in markerprotein_relab)
    marker_covered_ratio = float(markerprotein_relab_valued/len(markerprotein_relab))
    if marker_covered_ratio >= marker_ratio:
#    if markerprotein_relab_valued > 0:
    ##VC reads counts
        markerprotein_reads = list(vc_marker_readsnum[vcid].values())
        markerprotein_reads = sorted(markerprotein_reads,key=int)
#    print(f'''{vcid}\t{markerprotein_reads}''')
        low_index = int(vc_markerprotein_num * 0.1)
        up_index = int(vc_markerprotein_num * 0.9)+1
        if up_index>low_index:
            subset_abun = markerprotein_relab[low_index:up_index]
            subset_abun_valued = sum(abun>0 for abun in subset_abun)
#            subset_abun_new = [new for new in subset_abun if new > 0]
            mean_abundance = np.mean(subset_abun) #get mean abudannce in middle 80% markers 
            subset_reads = markerprotein_reads[low_index:up_index]
            sum_reads_num = np.sum(subset_reads)
            if mean_abundance > 0:
            #    print(f'''{vcid}\t{markerprotein_reads}''')
                real_marker_num = int(up_index-low_index)
                relative_abudannce = mean_abundance/sum_vc_mean_abun
                relative_reads_num = sum_reads_num/sum_vc_reads_num
                if relative_abudannce > 0:
                    outfile.write(f'''{vcid}\t{vc_markerprotein_num}\t{markerprotein_relab_valued}\t{real_marker_num}\t{subset_abun_valued}\t{mean_abundance}\t{relative_abudannce}\t{sum_reads_num}\t{relative_reads_num}\n''')
                else:
                    continue
            else:
                continue
        else:
            continue
    else:
        continue
