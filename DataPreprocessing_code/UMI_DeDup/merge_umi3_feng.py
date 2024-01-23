#!/opt/conda/envs/condaFeng3.7/bin/python3 ##running on feng laptop P52s
#
### #!/home/ffeng/miniconda3/envs/FengPython3.7/bin/python3  #running on theano
#
#   #!/opt/conda/bin/python3.8
#!/usr/bin/python3
###
## #!/usr/bin/python3
#!/opt/conda/bin/python3.8   #please check and modify this one to use the python in the running machine. this is for P52s debian feng
## #!/usr/bin/python3   #this is for runninng on thoen
####   

import sys
import csv
import os
import operator
import subprocess
import time
import operator
from datetime import datetime
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from Bio import SeqIO
from Bio import Phylo
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from io import StringIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator 
from Bio.SeqIO.FastaIO import SimpleFastaParser
from itertools import combinations
import matplotlib  #these two lines (including one line below) are necessary for ssh which doesn't have X11 server
matplotlib.use('Agg')  #so who do plotting please remember to save them.

import matplotlib.pyplot as plt
from multiprocessing import Pool
import multiprocessing
import numpy as np
import argparse
import threading

#=================================
#for doing pythonnet, .NET calling
import clr
clr.AddReference('System.Collections')
from System.Collections.Generic import Dictionary
from System import String
from System.Collections.Generic import List
from System import Array
from System import Int32
from System import Byte

#assembly_path=os.path.join(os.path.dirname(os.path.abspath(__file__)),'packages/CloanalystDLL')
assembly_path=os.path.join(os.path.dirname(__file__),'packages/CloanalystDLL')
#assembly_path="./packages/CloanalystDLL"
#print(assembly_path)
sys.path.insert(0,assembly_path)
clr.AddReference("CloanalystDLL")


from CloanalystDLL import Interface
from CloanalystDLL import FileManager
from System.IO import FileStream
from System.IO import FileMode
#sys.path.append("/home/feng/Projects/BULQI/CloanalystCore/CloanalystCore/packages/Troschuetz.Random.4.1.3/lib/net46/")
clr.AddReference("Troschuetz.Random")
#sys.path.append("/home/feng/Projects/BULQI/CloanalystCore/CloanalystCore/packages/Thrower.4.0.6/lib/net46/")
clr.AddReference("PommaLabs.Thrower")

#for the modules coming together with the deployment
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'packages'))
#sys.path.insert(0, "./packages")
#sys.path.insert(0, "../packages")
from vote_consensus import VoteConsensus
import matplotlib.pyplot as plt
from simple_graph import Graph
import pFunction
import umis.network as network

import ctypes as ct
import AlignUtility as au

import subprocess
from project import Project

#gettign project range.
project = Project()
sample_list = project.sample_list
print("Sample list:", sample_list)

PDEBUG=project.DEBUG
run_Cloanalyst=project.Cloanalyst   #<---Note: this is calling .NET to do cloanlyst on UMI groups, it will slow down the run significantly.
#becareful to call it. By default this is not run, but you can modify the .config file to revise it.
#also this one is mutually exclusive with the sliding-window check on the umi.

#this code tries to separate reads by umi - and then 
#assemble each umi into one or more "compatible" consensus

#steps;
#1. run through trimming primers---
#2. do consensus by voting...
#3.....

###############################
# setup command line options.
###############################
parser = argparse.ArgumentParser(description='process fastq file and assemble reads by umi ')
#parser.add_argument('-c', '--consensus', action='store_false',  help='use consensus, instead of best read as output, default false')
#not that if not trimming primers, we are still triimming low quals, so the file generating sequence is still the the same.
args = parser.parse_args()

remove_singleton = project.Remove_Singleton
max_base_diff = 5      #for the first 100 bases, that is 95% similarity
#max_base_diff_reverse=
save_assembly = project.save_assembly
#not used anymore
right_align_trim = 4    #was used to discard 3'end of alignments, but since we are trimming primers off. so this is not needed for now.
                                #we now set it to be 4, since the read 2 has 4 random nt. updated by Feng 2020 May.
                                #we don't get rid of isotypes primers now. This right_align_trim will be used by process_file function to get rid of the 4 random nt.
max_N = 20              #maximum number of Ns

#if using dedup
dedup_option = project.Dedup_Option        #0 - don't do deup; 1- do dedup before assembly. 2- dedup afeter assembly.
dedup_prefix_len = 15

#level 0 at project level, 1 within project
parallel_level = project.Parallel
print("parallel_level:", parallel_level)
print("run_Cloanalyst:", run_Cloanalyst)

#species
species =project.species

#paired ?? and read2
Paired =project.isPaired
Joined=project.Joined

if (parallel_level==1) :
        print("  it is 1");
        
if (not parallel_level==1) and run_Cloanalyst:
        #this is not available for now.
        print("Currently running in parallel and calling to do Cloanalyst are not possible!!")
        print("\t please modify the \".config\" to do only one option.")
        sys.exit(-1)

#return index (0-base) range, note the range is [x0, x1)
def get_aligned_range(seq):
        start1 = 0
        for i in range(len(seq)):
                if seq[i] != '-':
                        start1 = i
                        break
        stop1 = len(seq)
        for i in reversed(range(len(seq))):
                if seq[i] != '-':
                        stop1 = i + 1
                        break
        return (start1, stop1)

##@profile
def count_base_differ(seq1, seq2):
        if len(seq1) != len(seq2):
                raise ValueError("sequences are of different length")
                sys.exit(1)
        #only check aligned part
        (start1, stop1) = get_aligned_range(seq1)
        (start2, stop2) = get_aligned_range(seq2)
        start = max(start1, start2)
        stop = min(stop1, stop2)

        ndiff = 0
        for i in range(start, stop):
                if seq1[i] != seq2[i]:
                        ndiff += 1
                        #if (ndiff > 5):return 100
        return ndiff

        #ndiff2 = sum(1 for a, b in zip(seq1[start:stop], seq2[start:stop]) if a != b)
        #return ndiff2

def score_alignment(seq1, seq2):
        if len(seq1) != len(seq2):
                raise ValueError("sequences are of different length")
        #only check aligned part
        (start1, stop1) = get_aligned_range(seq1)
        (start2, stop2) = get_aligned_range(seq2)
        start = max(start1, start2)
        stop = min(stop1, stop2) #range is NOT right inclusive
        ndiff = 0
        for i in range(start, stop):
                if seq1[i] != seq2[i]:
                        ndiff += 1
        return (stop - start) - ndiff * 4

#the partition_alignments was too slow
#this try to speed it up.
#read_dict contains key - readname, value- tuple of seq. qual...?
#we will convert the seq 
## some strange "low qual at right end???"
##@profile
def partition_alignments(read_dict, right_trim = right_align_trim):
        keylist = list(read_dict.keys())
        if (len(keylist) > 1000):
                print("making seq dictionary...", file=sys.stderr)
        #create a numpy array to speed up array comparison.
        #since we are using index, we should not need dictionary.
        #array should be fine.
        seq_array = []          #holds readname -> numpy char array
        aln_range = []          #holds aligned range of the the read
        id_read_dict = {}       #same as read_dict, but with int id

        for idx, item in enumerate(keylist):
                seq1 = read_dict[item][0]
                seq_array.append(np.array(list(seq1)))
                aln_range.append(get_aligned_range(seq1))
                id_read_dict[idx] = read_dict[item]

        #computing differences...
        if (len(keylist) > 1000):
                print("done.")
                print("finding read differences...", file=sys.stderr)
        pair_score = [] #{}
        aln_len = len(read_dict[keylist[0]][0])

        show_progress = len(keylist) >= 1000
        #if show_progress:
        #       start_time = time.time() 
        show_unit = len(keylist)//100
        keysize = len(keylist)
        for idx in range(0, len(keylist)-1):
                seqidx = seq_array[idx]
                (start1, stop1) = aln_range[idx]
                if show_progress:
                        if idx%show_unit == 0:
                                print(".", end='', file=sys.stderr)
                                #print(pid, ",", end='', file=sys.stderr, flush=True)
                for idy in range(idx+1, keysize):
                        (start2, stop2) = aln_range[idy]
                        start = max(start1, start2)
                        stop = min(stop1, stop2) - right_trim
                        #score = np.sum(seq_array[idx][start:stop] == seq_array[idy][start:stop])
                        score = np.sum(seqidx[start:stop] == seq_array[idy][start:stop])
                        ndiff = stop - start - score
                        if ndiff > max_base_diff:
                                continue
                        score -= ndiff * 2
                        pair_score.append(((idx, idy), score))

                        # seq1 = read_dict[keylist[idx]][0]
                        # seq2 = read_dict[keylist[idy]][0]
                        # ndiff2 = count_base_differ(seq1, seq2)
                        # if (ndiff != ndiff2):
                        #         print("diff differs", ndiff, ndiff2)


                
        #build a network and find number of "connected graph"
        #given the pair_score
        graph = Graph(len(keylist))
        for pair in pair_score:
                (idx, idy) = pair[0]
                graph.addEdge(idx, idy)
        cc = graph.connectedComponents()
                
        if show_progress:
                print("done", file=sys.stderr)
                print("sorting scores:", file=sys.stderr, flush=True)
        sorted_pairs = sorted(pair_score, key=operator.itemgetter(1), reverse= True)
        if show_progress:
                print("done", file=sys.stderr)
                print("finding groups...", file=sys.stderr, flush=True)

        if show_progress:
                #end_time = time.time()
                #print("timecost=", end_time - start_time, "seconds")
                print("pair score finding done.")

        groups_con = {}         #consensus dictionary
        item_group = {}         #read -> consensus key        
        group_id = 1    #new group id.
        for pair in sorted_pairs:
                (idx, idy) = pair[0]
                name1 = keylist[idx]
                name2 = keylist[idy]
                score = pair[1]
                grp1 = item_group.get(idx, -1)
                grp2 = item_group.get(idy, -1)
                if (grp1 != -1 and grp2 != -1):
                        if grp1 == grp2: continue
                        cons1 = groups_con[grp1]
                        cons2 = groups_con[grp2]
                        #check the difference between the two consensus
                        (start1, stop1) = cons1.aln_range
                        (start2, stop2) = cons2.aln_range
                        start = max(start1, start2)
                        stop = min(stop1, stop2) - right_trim
                        score = np.sum(cons1.consarray[start:stop] == cons2.consarray[start:stop])
                        ndiff = stop - start - score
                        if ndiff > max_base_diff:
                                continue
                        for item2 in cons2.read_list:
                                item_group[item2] = grp1
                        cons1.merge(cons2)
                        del groups_con[grp2]
                        #print("group", grp2, "merged to", grp1)
                elif grp1 != -1:
                        cons1 = groups_con[grp1]
                        (start1, stop1) = cons1.aln_range
                        (start2, stop2) = aln_range[idy]
                        start = max(start1, start2)
                        stop = min(stop1, stop2) - right_trim
                        score = np.sum(cons1.consarray[start:stop] == seq_array[idy][start:stop])
                        ndiff = stop - start - score
                        if ndiff > max_base_diff:
                                continue
                        cons1.add_read(idy)
                        item_group[idy] = grp1                        
                elif grp2 != -1:
                        cons2 = groups_con[grp2]
                        (start1, stop1) = aln_range[idx]
                        (start2, stop2) = cons2.aln_range
                        start = max(start1, start2)
                        stop = min(stop1, stop2) - right_trim
                        score = np.sum(seq_array[idx][start:stop] == cons2.consarray[start:stop])
                        ndiff = stop - start - score
                        if ndiff > max_base_diff:
                                continue
                        cons2.add_read(idx)
                        item_group[idx] = grp2                        
                else:
                        item_group[idx] = group_id 
                        item_group[idy] = group_id
                        groups_con[group_id] = VoteConsensus(id_read_dict, [idx, idy])
                        group_id += 1
        if show_progress:
                print("done.", file=sys.stderr, flush=True)
        id_groups = []
        for item in groups_con:
                cons1 = groups_con[item]
                idlist = [keylist[i] for i in cons1.read_list]
                id_groups.append(idlist)
        if not remove_singleton:
                for idx, item in enumerate(keylist):
                        if item_group.get(idx, -1) == -1:
                                id_groups.append([keylist[idx]])
        if len(id_groups) > 3:
                check_point_1 = 1

        if len(id_groups) != len(cc):
                check_point = 2

        return id_groups


def GetAlignMap(alignseq):
        index = 0
        align_map = [0] * len(alignseq)
        for i, c in enumerate(alignseq):
                if c == '-':
                        align_map[i] = -1
                else:
                        align_map[i] = index
                        index += 1

        for i, val in enumerate(align_map):
                if val >= 0:break
                
        return align_map


#make pmf consenus for readlist, 
#foreach id, the read_dict contains a tuple of (paddedSeq, unpadded qual)
#it seems that cons_len has not been used. 
def make_pmf_consensus(read_dict, readlist, cons_len, umi_code = ""):
        try:
            aut = au.AlignUtil()
            #aut.PrintAlignment()
            if (umi_code != ""):
                    umi_code = umi_code + "_"
            for id in readlist:
                    (padseq, qstr)  = read_dict[id]
                    seq = padseq.replace("-", "")
                    if len(seq) != len(qstr):
                            raise Exception("Read sequence and qual has to be same length")
                    qlist = [int(q)-33 for q in qstr.encode()]
                    qual = bytes(qlist)        
                    aut.AddQtSequence(umi_code + id, seq, qual, GetAlignMap(padseq))
                    #aut.PrintAlignment()
                    
            #aut.PrintAlignment()
            cons = aut.GetConsensus()
            if PDEBUG:
                print("--inside make_pmf_consensu.........ready to call printing cons", flush=True)
                aut.PrintQtSequence(cons)
                print("\tcheckpoint", flush=True)        
            #sys.exit()
            #clean up memory , is this right?
            #aut.DeleteQtAlignment()
            return cons
        finally:
            aut.DeleteQtAlignment()
            del aut
            
## find connected graph- nodes are connect if they overlap length >= hash_length
## since the reads have umi in the front and cut off, so the start should almost always 
## Feng, this function to further parse/partition the umi groups. since identical umi_code alone doesn't 
## guarantee the identical original RNA molecules. we check the first 100 nts (allowing 5 mutatoins)
## and decide the groupings. (first 100nt is not value now, 9/4/2020 Feng. we use the whole sequence )
## the input 
#                       read_dict : the umi -group sequences and they have been multi-aligned!!!
###                         input has a format of {readID : (align seq , quality score)}
#                       umi_code: the umi code for this current group 
#                       hash_length =200 , now this is not used (as of 9/4/2020, since we use the whole sequence to be hash)
#                                                                           the hash is the first step to group "identical" sequences within each UMI group.
#                       read_dict _r1: the umi-group sequences for read 1 in this case. they have been aligned. The rationale is
#                                   if there are single read and joined paired read, then this is None. If we do unjoined paired reads,
#                                    read_dict above will be read2, since we take care read 2 mainly for grouping or ungrouping. this 
#                                       read_dict_r1 is for read1 and need to do sliding window for checking the overall errors/mismatches.
#                                    in any case, we do read 1 through this input!!!
## output: return a list of read ID list [[readID1, readID3], [readiD2, readID4], [readID5].....]
#
#---update, now make the check sequence from the back(R2, rev'comp'ed), of course we remove the last 4 nt
#     as well as the primer , for mouse it is 25nt long
#       but constant region is variable based on the isotype, so we check the last 40nt for <=1 mismatch.
#       update000: now we don't use hash_length, but the whole length (aligned length)
def partition_alignments_3(read_dict, umi_code, hash_length= 200, read_dict_r1=None):
        len_rdNt=0 #for the right random nucleotide. here we set this to be zero. because this has been done by right_align_trim and not read in now. 
        len_mouse_prime=25  # this is not used now. we need to conisder the primer and the whole constant region.
        #hash_length=hash_length<100 ? 100:hash_length
        if PDEBUG:
                print("+++++++Inside partition alignments_3()+++++++ with UmI:", umi_code)
        keylist = list(read_dict.keys())
        if (len(keylist) > 1000):
                print("making seq dictionary...", file=sys.stderr)
        
        #check the last 100 bases. NOTE: aln_len is the sequence length, not the alingment length 
        #WE MIGHT NEED TO modify this !!!
        aln_len = len(read_dict[keylist[0]][0])-len_rdNt  #all sequence in the group should have the identical aln_length. they are multiple-aligned
        hash_length=aln_len
        #check_len = min(aln_len, hash_length)
        check_len=aln_len
        if(aln_len<100 ):
            print("aln_len too short, please check!! (return whole list without parsing )")
            #system.exit(-1)
            print("ID:", keylist)
            return [keylist]
        alpha=0.05
        #max_base_diff_2=0.05*check_len;#define the max possible errors
        #print("aln_len :", aln_len, "check_len:", check_len, "(aln_len-check_len):", (aln_len-check_len))
        #print("aligned seq:", read_dict[keylist[0]][0])
        #print("seqhash:" ,read_dict[keylist[0]][0][(aln_len-check_len):aln_len])
        
        #create a numpy array to speed up array comparison.
        #since we are using index, we should not need dictionary.
        #array should be fine.
        seq_array = []          #holds aligned sequence -> numpy char array
        id_read_dict = {}       #same as read_dict, but with int id, read_dict has readID (string) as the key
        seq_dict = defaultdict(list)    #<seqhash, list of ids> sequence read_indexes. and then only work on seq_dict
        for idx, item in enumerate(keylist):  #idx is the index id of the item, and item is the readID for each read
                #if PDEBUG:
                #        print("inside partition:", flush=True)
                #        print("idx:", idx, flush=True)
                #        print("item:", item, flush=True)
                
                seq1 = read_dict[item][0]  #sequence, aligned sequence  
                #seq1 is the aligned part of the sequence length, need to remove "-" 
                #seqhash =  seq1[0: check_len]    #using the first 100(?) (including
                seqhash =seq1; #seq1[(aln_len-check_len):aln_len];  #now checking form the 3' end, either joined reads or R2 (revcomp'ed)
                if Paired and not Joined:
                        try:
                                seqhash = read_dict_r1[item][0]+seqhash
                        except KeyError:
                                print("Error: KeyError discovered\n")
                                print("\t key length:", len(read_dict.keys()), " and ", len(read_dict_r1.keys()))
                                print("\t keys for R1:", list(read_dict_r1.keys()))
                                print("\t keys for R2:", list(read_dict.keys()))
                                sys.exit(-2);
                        except:
                                print("Error: unexpect error:", sys.exc_info()[0])
                                sys.exit(-1)
		#seqhash = read_dict_r1[item][0]+seqhash
                seq_dict[seqhash].append(idx)          
                seq_array.append(seqhash)                       #np.array(list(seqhash)))               
                id_read_dict[idx] = read_dict[item]  # read ids
        #done for the first step, put them in order based on the first (actually the last ) 100 nts
        if len(seq_dict)<=1:   #<-modified 9/3/2020, to only allow one single hash array to return early, 
        #if mismatch_count == 0:
                #then every reads are in one group
                #don't build the graph
                if PDEBUG:
                        print("len of seq_dict ==0 and every reads are in one group, returning to the caller")
                return [keylist]
        #computing differences...
        if (len(keylist) > 1000):
                print("done.")                

        #we now work on seq_dict and only convert things in seq_dict to numpy array. 
        
        #note, hashlist is the list of unique entries of hash sequnece
        #so even though there are 1000 sequences
        #the hashlist might be only 1
        #that means, every sequence has the first 100 bases.
        #that is why we see a lot of sequence
        #but we don't show_progress, because it is kind of a small number mostly.
        hashlist = list(seq_dict.keys())
        if PDEBUG:
                print("after hasing in partition, the hashlist (",len(hashlist), ":", hashlist, flush=True)
                for idk, hseq in enumerate(hashlist):
                        print("***hseq:",idk, "==", hseq)
                        for idz in seq_dict [hseq]:
                                print("\t --keylist: ", keylist[idz]);
        hash_array = []         #numpy char array
        for i in range(0, len(hashlist)):#seqhash in hashlist:  #now sequence hash is the whole sequence in numpy array!!!
                hash_array.append(np.array(list(hashlist[i].upper())))
        if PDEBUG:
                print("after hasing in partition, the hash_array:", hash_array, flush=True);
        show_progress = len(hashlist) >= 500
        if show_progress:
                print("finding read differences...", file=sys.stderr)
                #if show_progress:
                #       start_time = time.time() 
                show_unit = len(hashlist)//100
        mismatch_count = 0
        #start_time = time.time()
        dot_count = 0
        #max_base_diff_2=len(hashlist[0].replace("-","") )*0.05;
        if PDEBUG:
                print("alpha:", alpha, " and ", len(hashlist[0].replace("-","") ), flush=True)
        pair_score = [] #{}   in this array, we hold the pairs which are very similar (ndiff< max_base_diff), excluding the ones that are very different 
        #now do pairwise difference
        if PDEBUG:
                print("start  doing the pairwise comparison between hash_sequences:")
        #pair_score are using the index of hashSeq in the hashList as identity
        
        for idx in range(0, len(hashlist)-1):
                seqidx = hash_array[idx]
                if show_progress:
                        if idx%show_unit == 0:
                                print(".", end='', file=sys.stderr, flush=True)
                                dot_count += 1
                for idy in range(idx+1, len(hashlist)):
                        seqidy=hash_array[idy]
                        ns=check_sequence_similarity(hashlist[idx], hashlist[idy], seqidx, seqidy, alpha)
                        if(ns<0):
                                mismatch_count += 1
                                continue
                        
                        pair_score.append(((idx, idy), ns))
        if PDEBUG:
                print("pair_score:",pair_score, "\n",flush=True); 
                print("mismatch :" , mismatch_count, "\n")
        #end_time = time.time()
        #print("--- %s seconds ---" %(end_time - start_time))
        if show_progress:
                print("done", file=sys.stderr, flush=True)
        
        #   #<-modified 9/3/2020, to only allow one single hash array to return early, 
        #if mismatch_count == 0:
        #        #then every reads are in one group
        #        #don't build the graph
        #        if PDEBUG==False:
        #                print("mismatch_count==0 and every reads are in one group, returning to the caller")
        #        return [keylist]
        
        if PDEBUG:
                print("mismatch_count>0 and building the graph for possible merging")
        #build a network and find number of "connected graph"
        #given the pair_score, 
        #so the graph and cc (connected components) coming from the pair_score.
        #therefore they both use the index of hashseq in the hashlist as identity
        graph = Graph(len(hashlist))
        for pair in pair_score:
                (idx, idy) = pair[0]
                graph.addEdge(idx, idy)
        cc = graph.connectedComponents()
        #note: each connectedComponent could contain many entry they are similar by meeting a cutoff
        #print("connected graph:", cc)
        if PDEBUG:
                print("show connectiong cc of hash/umis:",cc, flush=True);
        #<----need to reenable this???
        #if len(cc) <= 1:
                #print("all connected, returning to caller");
        #        return [keylist]
        
        if PDEBUG:
                print("not all connected, prepare for building concensus and try one last time to compress/merge the concensus")
        #the idea is that we build consensus, so the conesuses could be connect or merged or compressed into bigger groups.
        #construct groups:
        id_groups = []  
        #this following for loop combine nodes and turn index groups into read id groups 
        
        #notation explaination:
        #group and cc: are using the index of hashseq in the hashlist as identity. indicating the hashSeq and hash group
        #hashlist : list of hashSeq 
        # seq_dict : a dictionary keyed by hashSeq and valued by list of index, which is index the read_id in the keylist
        #   keylist is the list of read_ids 
        #print("Before calling combineConnectedNodeGroup_________________")
        for group in cc:  #for each hash group in the connected graph 
                #print(group)                
                #---->
                #idlist = []
                #for index in group: #for  each index in the connected group 
                #        key = hashlist[index]   # this is the hash
                #        read_ids = seq_dict[key]   #this is the read information {idx , readIDs} that for each hash 
                #        items = [keylist[i] for i in read_ids] #
                 #       idlist.extend(items)
                ##if len(idlist) > 1:
                #id_groups.append(idlist)
                #<-------
                id_groups=combineConnectedNodeGroup(group, hashlist, seq_dict, keylist, graph, id_groups)
                #print(">>>>>>>>>idlist:",idlist)
                #id_groups.append(idlist)
                #print(">>>>>>>>>id_groups:", id_groups)
        #sys.exit(-1)
        if PDEBUG:
                print("id_groups by connecting graph components:",id_groups, flush=True)
        #testing compress(merging) groups...
        #what is the compressing??? and why we need to do compressing??
        if len(id_groups) > 1:  #will alsways have more than one group, since we check a couple lines ahead
                #double check if objects got deallocated.
                if PDEBUG:
                        print("there are more than one id_groups, so need to do concensus and compressing")
                
                #first check the constant region, using 35nt long right end sequence for testing first 
                compressor = au.AlignUtil()
                if PDEBUG:
                        print("checking groups before doing compressing  Constant REGION:", flush=True)
                min_len=10000;
                #read_dict_short={}
                #print(id_groups)
                for group in id_groups:
                        if PDEBUG:
                                print("group in comprssing block:", group, flush=True)
                        
                        #first make a shortened sequence
                        #read_dict_short=create_short_read_dict(read_dict, group, aln_len-35, aln_len) #again 35 is for isotype constant region len IgM, IgG is a little longer, but we just use it now.
                        cons = make_pmf_consensus(read_dict, group, check_len, umi_code); #check_len is not used, just leave it as it is.
                        if min_len>compressor.GetQtSeqLength(cons):
                                min_len=compressor.GetQtSeqLength(cons)
                        #compressor.PrintQtSequence(cons)
                        compressor.AddQtSequenceObject(cons)
                        if PDEBUG:
                                #print(cons, flush=True)
                                compressor.PrintAlignment();
                if PDEBUG:
                        print("##########min_len:", min_len, "; align length", aln_len)        
                        print("()()() partion alignment 3 after adding QtSequenceObject consensus:", flush=True)
                        compressor.PrintAlignment();
                compressor.AddCompressCriterion(0, aln_len, aln_len*0.95*1.30 ); #first criterion 0.95 match
                compressor.AddCompressCriterion(aln_len-25-20, aln_len-25, 19*1.3); #first criterion 0.95 match
                compressor.AddCompressCriterion(aln_len-90, aln_len-40, 50*0.95*1.30 ); #first criterion 0.95 match
                
                merge_groups = compressor.Compress(1.0) #we don't 1.0 criterion anyway, since the criterion is fed in above functions
                                        #0.90 is a relaxed version of 95/100, 1.38~=-LogPrior (log(1/4))
                                        #note about this setting. when doing compress(), my understanding is that we compress/merge the consensus.
                                        #so the rationale is that the concensus might be similar, so relax the 0.95 and 1.38 is the perfect match with very good quality
                #merge_groups=[]
                compressor.DeleteQtAlignment();
                
                if len(merge_groups) <=0:
                        if PDEBUG:
                                print("compressing...........:", merge_groups, flush=True)
                                print("\nNumber of compression:", len(id_groups), "=>", len(merge_groups), " ", umi_code)#, file=sys.stderr)
                        #---------------here, disabled on 8/10/2020, might need to enable later.
                        #print("\nNumber of compression:", len(id_groups), "=>", len(merge_groups), " ", umi_code, file=sys.stderr)
                        return id_groups
                #print("\nNumber of compression:", len(id_groups), "=>", len(merge_groups), " ", umi_code, file=sys.stderr)
                #print("\n merged groups:", merge_groups)
                #next, we need to do it again, now with "middle part of the sequence" or more precisely the CDR3 region 
                #compressor = au.AlignUtil()
                #if PDEBUG==False:
                #       print("checking CDR3 length groups before doing compressing :", flush=True)
                #min_len=10000;
                #read_dict_short={}
                #for group in id_groups:
                        #if PDEBUG:
                        #        print(group, flush=True)
                        
                        #first make a shortened sequence
                #        read_dict_short=create_short_read_dict(read_dict, group, aln_len-100, aln_len-50) #again 50 is for isotype constant region len IgM, IgG is a little longer, but we just use it now.
                #        cons = make_pmf_consensus(read_dict_short, group, check_len, umi_code)
                #        if min_len>compressor.GetQtSeqLength(cons):
                #                min_len=compressor.GetQtSeqLength(cons)
                #        if PDEBUG:
                #                print(cons, flush=True)
                #        compressor.AddQtSequenceObject(cons)
                
                
                #if PDEBUG==False:
                #        print("##########min_len:", min_len)        
                #        print("()()() partion alignment 3 after adding QtSequenceObject consensus:", flush=True)
                #        compressor.PrintAlignment();
                #merge_groups = compressor.Compress(min_len*0.95*1.30) #0.90 is a relaxed version of 95/100, 1.38~=-LogPrior (log(1/4))
                                        #note about this setting. when doing compress(), my understanding is that we compress/merge the consensus.
                                        #so the rationale is that the concensus might be similar, so relax the 0.95 and 1.38 is the perfect match with very good quality
                #after merging by Compress(), we got index of (groups) in the id_group s. (each group in id_groups could have many seqs 
                #print("\nNumber of compression:", len(id_groups), "=>", len(merge_groups), " ", umi_code, file=sys.stderr)
                #print("\n merged groups:", merge_groups)
                #now in here we turn group index (output from compress()) into sequences ids in order to send back to the caller 
                if len(merge_groups) > 0:
                        #check_point = 2
                        for items in merge_groups:
                                new_group = []
                                for index in items:
                                        new_group.extend(id_groups[index])
                                        id_groups[index].clear()
                                id_groups.append(new_group)
                        id_groups = [x for x in id_groups if x != []]
        return id_groups

#in this function, we take the connected node group and try to 
#combine them or break the connections based on the directional 
#network algorithm. Tom Smith UMItools genome research
#we  sort the connected nodes based on the size of the nodes (each node here  is one hash group)
# and check each nodes by turn 
# when we check each nodes we follow its edges and find its connected largest node
# and then check to see whether we combine them (small into big, when size of small <=1/2 of big) or break them 
# if we combine them, we update the position of them in the to_do_list
# else we break them,  write the small to the output list
# for either case, we will break the edges come out from this node
# and then loop through again to next node
#  input:   nodes: connected node group, contain indexes referring to the index of hashSeq in the hashlist, 
#                   hash_list: the list contains the hashseq, so to look up the hash sequence for each index in nodes
#                   hash_dict: the hash seq dictionary keyed by hash sequence (elements of hash_list )  and valued 
#                                   by the list of index to the keylist (which is a list of read ids)(so that we can get the size of each nodes)
#                   ReadID_list : list of read ids, so that can be indexed by hash_dict value list .(contain ids, so to look up the read ids to add to output id_list)
#                                                                               this is called keylist in above
#                   graph : the original connected group, in order to extract edges and vertexes, etc.
#  output: list of list of read ids 
def combineConnectedNodeGroup(nodes,   hash_list, hash_dict, ReadID_list, graph, idlist):
        templist=idlist
        if len(nodes)==0:
            return templist
        if len(nodes)==1:
            #print("***single node group --nodes[0]:",nodes[0])
            #print("hash_list[nodes[0]]:", hash_list[nodes[0]])
            #print("hash_dict[hash_list[nodes[0]]]:", hash_dict[hash_list[nodes[0]]])
            templist.append([ReadID_list[x] for x in hash_dict[hash_list[nodes[0]]]])
            return templist 
        #get the size of each node and be ready to do the sorting next
        size_nodes=[len(hash_dict[hash_list[i]]) for i in nodes]
        #print("size_nodes:", size_nodes)
        #print("nodes:", nodes)
        #sort the nodes based on the size of node (hash group)
        node_size_pairs=tuple(zip(size_nodes, nodes))
        sorted_size_nodes, sorted_nodes=zip(*[(y, x) for y,x in sorted(node_size_pairs)])
        sorted_size_nodes=list(sorted_size_nodes)
        sorted_nodes=list(sorted_nodes)
        #print("sorted_size_nodes:", sorted_size_nodes)
        #print("sorted_nodes:", sorted_nodes)
        #loop through from the smallest to the largest
        for i in range (len(sorted_nodes)):
                #print("\tDoing nodes i:", i, "-----------")
                #follow the edges out of this nodes find the largest node
                n,s, ind =graph.findConnectedNode_Largest(sorted_nodes[i], sorted_nodes, sorted_size_nodes)
                #print("\tn:", n , "; s:", s, "; ind:", ind )
                if n==-1 and s==-1: #no connected edges, because breaking of edges in the previous step
                        templist.append([ReadID_list[y] for y in hash_dict[hash_list[sorted_nodes[i]]]])
                        #print("\ttemplist:",templist)
                        continue
                        #sys.exit(-1)
                #find the largest nearest neighbor (direct link neighbor)
                #check to merge or break. to merge we need the a sizeable group with large # of members to
                # be only half size of the bigger one, or the group is to small (<3 or how small is too small??)
                if sorted_size_nodes[i]<=0.5*s or sorted_size_nodes[i]<=3 :
                        #do combining, update has_dict, sorted_size_nodes, and sorted_nodes. we need to update them in order
                        hash_dict[hash_list[n]].extend(hash_dict[hash_list[sorted_nodes[i]]])
                        sorted_size_nodes[ind]=sorted_size_nodes[ind]+sorted_size_nodes[i]
                        #now update the order of the sorted nodes and sorted nodes_size
                        graph.updateSortedConnectedNodeGroup(ind, sorted_nodes, sorted_size_nodes)
                        #print("\tafter updaing: sorted size:",sorted_size_nodes)
                        #print("\tafter updateing: sorted node:",sorted_nodes)
                else:
                        #in this case, we write it to the output
                        templist.append([ReadID_list[y] for y in hash_dict[hash_list[sorted_nodes[i]]]])
                        #print("after breaking templist:", templist)
                #remove edges out of this nodes
                graph.removeEdgesOfVertex(sorted_nodes[i])
                #done for loop 
        #print("done for loop and ready to return")        
        return templist        

#check the extra conditions for score a pair of sequence in order to do the connect group
#in partition alignment3 
#for now we only do for mouse. We do 3 things, 
#1)overall sequence similarity, total length, <0.05 errors
#2)constrant region. this is complicated. if we have KGTGTG, then only allow zero error. others we allow 1 error  
#3) CDR3 region, 0.05 region, can not all crowed in here.
#input the aligned seq1 and aligned seq2
#           numpy array seq1 and numpy array seq2.
#           assume all the input sequences are full length (aligned), no hashing
#           alpha : the rate for error , 0.05. previously it was max_base_diff, set up by the caller, usually 0.05*lenOfSeqnece(without "-")
#Output, ndiff, -1 means the sequences did not pass the filters and are not similar
#           otherwise it is positive value of ndiff (max n of difference of the whole sequence )
def check_sequence_similarity(seq1, seq2, np_seq1, np_seq2, alpha):
        #print("comparing ", idx, " and ", idy , "\n")
        #ndiff_max = np.sum(np_seq1 != np_seq2)
        #print("\t ndiff:", ndiff, ";")
        #we switch to do sliding window
        windowSize=100
        slideStep=50
        alpha=0.05
        sstart=0;
        send=len(seq2)-25
        ndiff_max, passed =pFunction.compareSeqSliding(np_seq1, np_seq2, sstart, send, windowSize, slideStep, alpha)
        #print("after sliding ndiff:", ndiff_max, " ; and passed :", passed)
        if not passed:  #ndiff_max > max_base_diff:
                #mismatch_count += 1
                return -1;
         #=============================No good section (new strategy below)===========================================       
        #check the read2 end for isotypes, for the last 6~7nt (IgM Ch1 by primer1, 30+/-, IgG ch1 primer 45+/-)
        # so we check for the last 6~7nt, and only allow no more than 1 mismatch
        # KGTGTG allow zero error and others allow for no more 1 mismatch
        #we also need to peel off the gaps at the end , find the first non-'-' char position 
        #=========================================================================
        #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        #Note: updated 8/27/2020, we change the strategy now!!
        #   The concern is now that the primer regoin is not good for demuxing. So now we change to distinguish ig subtypes.
        #   Rule of thumb: ignore the 25 very last nts and then check the next 20nts for at most two more mismatches.
        #           two more things to note: peel off the common gaps as above and this can not distinguish IgG2a from IgG2c.
        #   This so far only works for mouse. Need to revise for human cases 
        endIndex1=len(seq1.rstrip("-"))
        endIndex2=len(seq2.rstrip("-"))
        endIndex1=endIndex2 if endIndex1 < endIndex2 else endIndex1
        #ndiff=np.sum (np_seq1[(endIndex1-6):endIndex1]!=np_seq2[(endIndex1-6):endIndex1])
        ndiff=np.sum (np_seq1[(endIndex1-25-20):(endIndex1-25)]!=np_seq2[(endIndex1-25-20):(endIndex1-25)])
        #print("\t last 40 sequence -- ndiff:", ndiff, ";")
        if ndiff > 1:    #now very strict rule, hope the merge in below will merge the similar ones
                #no more than 1 mismatch<-----previously
                #now we need 2 mismatches for distinction
                return -1;
        #now check for a special case to distiguish between "GGGTGTG" vs "GTGTGTG", this was not used since 8/27/2020
        #if(species=="Mouse"):
        #    #print("Checking species ************************")
        #    #print("\tnp array:",seqidx[(hash_length-7):])
        #    #print("\tnp array:",hash_array[idy][(hash_length-7):])
        #    #still not perfect, did not take care of gap in the last six nts. 
        #    if(np.sum(np_seq1[(endIndex1-6):endIndex1]!=np.array(list("GGTGTG")))==0 or np.sum(np_seq1[(endIndex1-6):endIndex1]!=np.array(list("TGTGTG")))==0):
        #            if(ndiff>0):#skip don't go forward, since they are not supposed to be in the same sequence.
        #                    return -1
        #    if(np.sum(np_seq2[(endIndex1-6):endIndex1]!=np.array(list("TGTGTG")))==0 or np.sum(np_seq2[(endIndex1-6):endIndex1]!=np.array(list("GGTGTG")))==0):
        #            #skip don't go forward, since they are not supposed to be in the same sequence.
        #            if(ndiff>0):
        #                    return -1
        #check the shorten seq, for condensed errors, assuming we are reaching the ***CDR3*** part.
        #might be neccessary. will check later
        ndiff = np.sum(np_seq1[(endIndex1-80):(endIndex1-30)]!=np_seq2[(endIndex1-80):(endIndex1-30)])
        #ndiff = np.sum(np_seq1[(endIndex1-75):(endIndex1-25)]!=np_seq2[(endIndex1-75):(endIndex1-25)])
        #ndiff = np.sum(seqidx[0:int(hash_length*0.5)]!=hash_array[idy][0:int(hash_length*0.5)])
        #print("\t ndiff:", ndiff, ";\n")
        if ndiff > 2.5:#max_base_diff_2/2: 50*0.5
                #mismatch_count+=1
                return -1
        #if we are here, we will return a good score 
        return ndiff_max
        
'''
this is the function to rewrite the read_dict data struction and make it ready 
to put send to make_pmf_consensus. The important point is that we need to 
only take the certain region of specific part of the original sequence for doing the concensus
and then doing the merging. the hardest part is that we need make sure the quality part and 
sequence part are of the same length.
read_dict is the dictionary [key=readid, value=(seq, qal)]
ids is the [readids] for being processed.
start is the starting index for seq and qal
end is the ending index (exclusive). 
Note: the range is [start, end)
return new dictionary with short seq and short quality and only holds the ids records.
'''
def create_short_read_dict(read_dict, ids, start, end):
        short_read_dict={}
        for i in ids:
                (seq, qal)=read_dict[i]
                sseq=seq[start:end]
                ln=len(sseq.replace("-",""))
                #figure out the starting point
                qstart=start
                if not start==0:
                        qstart=len(seq[0:(start-1)].replace("-",""))
                sqal=qal[qstart:(qstart+ln)]
                short_read_dict[i]=(sseq,sqal)
        return short_read_dict

## traverse the phylo tree and build consensus from tips toward root.
## the tree object contains gloal informaiton about the tee, such as whether its roted
## it has one root clade, and under thant, it's nested lists of clades all the way down
## to tips (ref: https://biopython-cn.readthedocs.io/zh_CN/latest/en/chr13.html)
def find_consensus(root):
        read_dict = root.read_dict
        #print(len(read_dict))
        if len(root.clades) == 1:
                root.cons = VoteConsensus(read_dict, [root.name])
                return
        for clade in root.clades:
                clade.read_dict = read_dict
                find_consensus(clade)
        #collect allreadys(?) to build consenseuns?


#find the consensus -
#for a,c,g,t,n,-
baseIdx = {}
baseIdx['-'] = 0
baseIdx['a'] = 1
baseIdx['c'] = 2
baseIdx['g'] = 3
baseIdx['t'] = 4
baseIdx['n'] = 5
idxBase = ['-', 'a', 'c', 'g', 't', 'n']

#given read_dict and readlist
#make a consensus from it.
#this makes using majority vote, without consider quality.
def make_consensus(read_dict, readlist):
        aligned_len = len(next(iter(read_dict.values()))[0])
        freqlist = np.zeros([6, aligned_len])
        #keep track of "confirmed" range.
        min1, min2 = float('inf'), float('inf')
        max1, max2 = float('-inf'), float('-inf')
        
        for item in readlist:
                seq = read_dict[item][0]
                #print("read_len=", len(seq))
                x0, x1 = get_aligned_range(seq)

                #find second min, max
                if x0 <= min1: min1, min2 = x0, min1
                elif x0 < min2: min2 = x0
                if x1 >= max1: max1, max2 = x1, max1
                elif x1 > max2: max2 = x1

                for j in range(x0, x1):
                        i = baseIdx[seq[j]]
                        freqlist[i, j] += 1
        #we want double read range as consenuss...to avoid one longer junk read
        #alwasy get selected.....when range is doubled, then it is kind of "confirmed" because
        #we only allow 5 differences.
        maxcols = np.argmax(freqlist, axis = 0)

        #exlcude single read ends for consensus
        if min2 > 0:
                maxcols[0:min2] = 0
        if max2 < aligned_len - 1:
                maxcols[max2:] = 0
        cons = "".join([idxBase[x] for x in maxcols])

        return cons



#seq may have gaps, create a numpy arary
#for the qual- gapqual is the average of the bounding bases
def make_aligned_qual(seq, qstr):
        if len(seq) == len(qstr):
                qual = np.array([(ord(c) - 33) for c in qstr])
                return qual
        #insert gaps...
        qual = np.zeros(len(seq))
        q_idx = 0       #next base qual
        qrange = range(1, len(qstr))
        for i, c in enumerate(seq):
                if c != '-':
                        qual[i] = ord(qstr[q_idx]) - 33
                        q_idx += 1
                elif q_idx in qrange:
                        qual[i] = (ord(qstr[q_idx-1]) + ord(qstr[q_idx]) - 66)/2
        return qual

#this compute consensus such that quality is also considered.
#i.e. vote with some quality.
#if confirm_trim = true. then trim consensus with coverage of 1
def make_consensus_2(read_dict, readlist, confirm_trim = False):
        aligned_len = len(next(iter(read_dict.values()))[0])
        freqlist = np.zeros([6, aligned_len])
        #keep track of "confirmed" range.
        min1, min2 = float('inf'), float('inf')
        max1, max2 = float('-inf'), float('-inf')
        
        for item in readlist:
                seq = read_dict[item][0]
                qual = make_aligned_qual(seq, read_dict[item][1])

                x0, x1 = get_aligned_range(seq)

                #find second min, max
                if confirm_trim:
                        if x0 <= min1: min1, min2 = x0, min1
                        elif x0 < min2: min2 = x0
                        if x1 >= max1: max1, max2 = x1, max1
                        elif x1 > max2: max2 = x1

                for j in range(x0, x1):
                        i = baseIdx[seq[j]]
                        q = qual[j]
                        freqlist[i, j] += q
        #we want double read range as consenuss...to avoid one longer junk read
        #alwasy get selected.....when range is doubled, then it is kind of "confirmed" because
        #we only allow 5 differences.
        maxcols = np.argmax(freqlist, axis = 0)

        #exlcude single read ends for consensus
        if confirm_trim:
                if min2 > 0:
                        maxcols[0:min2] = 0
                if max2 < aligned_len - 1:
                        maxcols[max2:] = 0
        cons = "".join([idxBase[x] for x in maxcols])

        return cons


#for one umi_code group, we want to further perform alignment of the actual sequences to see 
#whether they are similar enough to coming from the same ancestor. 
#
#given list of seq record (who have the "identical umi" - perform alignment
#make and return conensus
#here record is a list of tuples, of the format(id, seq, qual), all of them string.
#input: (umi_code, seqlist), umi_code is the single umi_code for this group .
#           seqlist, the list of all sequences belong to the same umi, include (read id, seq, qual) 
#output:  [umi_code, ret_list, asm_str, umi_read_dict, ret_list_R2, asm_str_2]
#                  umi_code, umi_code for this group, this is identical to the input umi for this code of sequences. Nothing changed. will be 
#                           ignored by the caller.
#                   ret_list, this is the sequnce to be kept and will be returned to the called as a fasta format. For singleton, we will have it as it is, for multiple 
#                           sequences, we will have concensus. these will be written as output
#                   asm_str, to write out the assembly file. 
#                   umi_read_dict, the dictionary with new umi_code (appended with number indicateing different sub-group within the umi code group 
#                              { umi_code#: [seq, quality]}. if no sub-group formed, we will not append # and just use the input one. Need to keep read1 and 2 both.
#                               this contains R2 information seq and qual if available
#                   ret_list_R2 and asem_str_R2 are newly added output for R2 information 
#--updated by Feng May 31 2020, add to output umi freq dist, for the purpose to show stats in the outer caller. 
#@profile

def find_candidate_sequences(umi_seqlist):
        
        #print("I'm process", os.getpid())        
        (umi_code, seqlist) = umi_seqlist
        #if len(seqlist)==3:
        #print("****UMI_code:", umi_code, "; len:", len(seqlist), "\n")
                #print("seq list:", seqlist,"\n");
        
        #if only one entry, then do nothing rand return...actually we don't pass in the singletons
        if len(seqlist) == 1:
                entry = "".join([">", seqlist[0][0], "\n", seqlist[0][1], "\n"])
                if Paired and not Joined:
                        entry_R2=SeqRecord(Seq(seqlist[0][3]), id=seqlist[0][0])
                        entry_R2=entry_R2.reverse_complement(id=entry.id)
                        entry_R2="".join([">", entry.id, "\n", str(entry_R2.seq), "\n"])
                return (umi_code, [entry], None, None, [entry_R2], None)  # <------- we probably will need to do this a little carefully to make it 100% right.
                                        #the thing is that we might not need it anyway. since if everything right, we shold not pass singleton to here.

        read_dict = {}   #this is the aligned pair, id : [seq, quality], for this whole umi group (umi = umi_code) 
        read_dict_R2={}
        #aligned_len = -1

        #calling mafft
        #create "fasta" text input and make them ready for calling mafft first
        #first let's do for Read 1 or Joined Reads 
        out = StringIO()
        #out_R2=StringIO()
        for rd in seqlist:
                print(">" + rd[0], file=out)
                print(rd[1], file=out)
                #if Paired and not Joined: #case for R1 and R2 separed 
                #        print(rd[3], file=out)
                #else: #joined or single read 
                        
        fasta_bytes = out.getvalue().encode('utf-8')
        
        #to save the pre_aligned sequences to file.
        #with open("pre-align.fasta", "w") as paf:
        #        paf.write(out.getvalue())

        if len(seqlist) > 500:
                print("\ncalling mafft(", len(seqlist), ")...", sep='', flush=True)
        mafft_cmd = ['mafft', '--quiet']

        # do we need parallel here???
        if len(seqlist) > 1000:
                 mafft_cmd.extend(['--thread', '8'])
        mafft_cmd.append('-')             
        #p = subprocess.Popen(['mafft', '-'], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        p = subprocess.Popen(mafft_cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE)

        (stdout_str, stderr_str) = p.communicate(input=fasta_bytes) #communicate return tuple (stdout, stderr)
        if p.returncode != 0:
                print("mafft alignemnt failed!")
                sys.exit(1)
        #the alignment output
        align_out = stdout_str.decode()
        if len(seqlist) > 500:
                print("returned from mafft", flush=True)
        #if saving alingment
        save_align = False
        if save_align:
               with open("tmp_align.fasta", "w") as ff:
                        ff.write(align_out)

        #read the mafft output into the run 
        # the output of mafft is fasta alignment
        #then read into a read_dict, which has a format of {readID : (align seq , quality score)}
        #note: it is weird to have quality score, since the quality score is not "align". will have trouble to use it (need to be parsed/aligned)
        #            note to note: this not aligned quality score has been taken care of in the later code, by removing "-" gaps.
        fasta_io = StringIO(align_out)
        idx = 0
        #prepare the aligned sequences and quality, and make them ready for doing partition_alignments
        for id, seq in SimpleFastaParser(fasta_io):
                si = seqlist[idx][0]  #seqlist is a list, so we can check it according to it position/index
                if (si != id):
                        print("id error, align output order different?")
                        sys.exit(1)
                read_dict[id] = (seq, seqlist[idx][2])    
                #if Paired and not Joined:  #case where R1 and R2 separated 
                #        read_dict[id] = (seq, seqlist[idx][4])
                #else :
                #        read_dict[id] = (seq, seqlist[idx][2])
                idx += 1
        if idx==0 or len(read_dict[id])==0:
                print("Error: did not read any entry from the alignment (Read 1). something wrong with alignment. Quit!!!", flush=True) 
                sys.exit(-1)

        #if len(seqlist)==3:
        #print("Now after alignment, we have :", len(read_dict))
        if Paired and not Joined: # we have separated R2 need to do alignment
                out = StringIO()
                #out_R2=StringIO()
                for rd in seqlist:
                        print(">" + rd[0], file=out)
                        print(rd[3], file=out)
                        #if Paired and not Joined: #case for R1 and R2 separed 
                        #        print(rd[3], file=out)
                        #else: #joined or single read 
                                
                fasta_bytes = out.getvalue().encode('utf-8')
                if len(seqlist) > 500:
                        print("\ncalling mafft(", len(seqlist), ")...", sep='', flush=True)
                mafft_cmd = ['mafft', '--quiet']

                # do we need parallel here???
                if len(seqlist) > 1000:
                         mafft_cmd.extend(['--thread', '8'])
                mafft_cmd.append('-')             
                #p = subprocess.Popen(['mafft', '-'], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
                p = subprocess.Popen(mafft_cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE)

                (stdout_str, stderr_str) = p.communicate(input=fasta_bytes) #communicate return tuple (stdout, stderr)
                if p.returncode != 0:
                        print("mafft alignemnt failed!")
                        sys.exit(1)
                #the alignment output
                align_out = stdout_str.decode()
                if len(seqlist) > 500:
                        print("returned from mafft", flush=True)
                #if saving alingment
                save_align = False
                if save_align:
                       with open("tmp_align_R2.fasta", "w") as ff:
                                ff.write(align_out)
                fasta_io = StringIO(align_out)
                idx = 0
                #prepare the aligned sequences and quality, and make them ready for doing partition_alignments
                for id, seq in SimpleFastaParser(fasta_io):
                        si = seqlist[idx][0]  #seqlist is a list, so we can check it according to it position/index
                        if (si != id):
                                print("id error, align output order different?")
                                sys.exit(1)
                        read_dict_R2[id] = (seq, seqlist[idx][4])#R2 quality    
                        #if Paired and not Joined:  #case where R1 and R2 separated 
                        #        read_dict[id] = (seq, seqlist[idx][4])
                        #else :
                        #        read_dict[id] = (seq, seqlist[idx][2])
                        idx += 1
                if idx==0 or len(read_dict_R2)==0:
                        print("Error: did not read any entries from the alignment (Read 2). something wrong and quit!!!", flush=True)
                        sys.exit(1)

        #--------------done with alignment.
        #partition alignment into groups
        #aligned_groups = partition_alignments(read_dict)
        if PDEBUG:
                print ("after alignment and ready to do partion and the input aligned sequences for this umi is:")
                #print(read_dict, flush=True)
                print("**********:",umi_code, flush=True)
                print(read_dict_R2, flush=True)
        #note, we don't change the partition_alignment_3 for R1 and R2 separated case !! sine they don't have to know 
        #the difference. the take in the aligned sequence/quality, and partition alignments3 and return the groups (list of id groups (list)) 
        aligned_groups =[]
        if Paired and not Joined:   #case for R2
                #print("passing read_dict_R2 to do partition alignment...........")
                aligned_groups=partition_alignments_3(read_dict_R2, umi_code, 100, read_dict)
        else:
                aligned_groups=partition_alignments_3(read_dict, umi_code)
        #if len(seqlist)==3:
        #print("aligned_group:", aligned_groups, "\n")
        if PDEBUG:
                print("output of the partition_alignments_3:",aligned_groups, flush=True)
        #if aligned_groups is None: return       
        #print("-------in caller: aligned groups:", aligned_groups)
        
        #----------------------------------------------
        #now all the job finished, now working on output 
        #now make a consensus for each group, select the one sequence that 
        #has mimimum number of differences with the consensu as the one entry to keep for each group.
        ret_list = []   #returning tuple of (consensus, new_umi_code, count), the latter 2 is used for "post" dedup
        ret_list_R2=[]
        if save_assembly:
                asem_out = StringIO()
                asem_out_R2=StringIO()
        #print("num groups = ", len(aligned_groups))
        #debug
        tmp_item = next(iter(read_dict.values()))
        debug_len = len(tmp_item[0])
        umi_read_dict=defaultdict(list)
        umi_run_count=0;
        umi_run_code=umi_code
        for group in aligned_groups:
                if len(group) == 1:    #singleton 
                        umi_run_count=umi_run_count +1;
                        if umi_run_count >1:
                            umi_run_code ="".join([umi_code, "_", str(umi_run_count)])
                        tseq = read_dict[group[0]][0]
                        if '-' in tseq:  #get rid of '-'
                                tseq = "".join(c for c in tseq if c != '-')
                        if Paired and not Joined: #R2 separated
                                tseq_R2=read_dict_R2[group[0]][0]
                                if '-' in tseq_R2:  #get rid of '-'
                                        tseq_R2 = "".join(c for c in tseq_R2 if c != '-')
                                        
                                umi_read_dict[umi_run_code ].append((group[0], tseq, read_dict[group[0]][1], tseq_R2, read_dict_R2[group[0]][1]))
                        else:
                                umi_read_dict[umi_run_code ].append((group[0],tseq, read_dict[group[0]][1]))
                        #umi_read_dict[umi_run_code]
                        if save_assembly:
                                #print("\t save a singleton\n")
                                print(">", group[0], "\n", read_dict[group[0]][0], sep="", file=asem_out)
                                print(">", file = asem_out)
                                if Paired and not Joined:
                                        print(">", group[0], "\n", read_dict_R2[group[0]][0], sep="", file=asem_out_R2)
                                        print(">", file = asem_out_R2)
                                
                        if remove_singleton:
                                continue
                        #we are here, mean we need to save singleton (not remove_singleton)
                        entry = "".join([">", group[0], "\n", tseq, "\n"])
                        ret_list.append(entry)
                        if Paired and not Joined:
                                entry=SeqRecord(Seq(tseq_R2), id=group[0])
                                entry=entry.reverse_complement(id=entry.id)
                                entry = "".join([">", entry.id, "\n", str(entry.seq), "\n"])
                                ret_list_R2.append(entry)  
                        continue
                #if we are here, mean we have more items in the group, need to do cloanalyst for partition
                # and also importantly if inside the cloanalystPartition function to check for whether we 
                # want to cloanalyst, if not, we simply return back to the caller the groups
                clone_groups=[]
                if Paired and not Joined:
                        clone_groups=cloanalystPartition(read_dict_R2, group)
                else:
                        clone_groups=cloanalystPartition(read_dict, group) 
                #print("# of groups after cloanalyst: ", len(clone_groups))
                #if len(seqlist)==3:
                #       print("clone_analyst parsing length:", clone_groups, "\n")
                for  clone_group_items in clone_groups:
                        #save it first, no matter whether it is singleton or not
                        if save_assembly:
                                for item in clone_group_items:
                                        #print("\tsave ",item, "\n")
                                        print(">", item, "\n", read_dict[item][0], sep="", file=asem_out)
                                        if Paired and not Joined:
                                                print(">", item, "\n", read_dict_R2[item][0], sep="", file=asem_out_R2)

                        #for cloanalyst left with singleton, we need to take care of this 
                        if len(clone_group_items) == 1:
                                umi_run_count=umi_run_count +1;
                                if umi_run_count >1:
                                        umi_run_code ="".join([umi_code, "_", str(umi_run_count)])
                                #save the singleton to the outpt dictionary for further stats, note we have write the singlton to the keep_list below
                                tseq = read_dict[clone_group_items[0]][0]
                                if '-' in tseq:
                                        tseq = "".join(c for c in tseq if c != '-') 
                                if Paired and not Joined:
                                        tseq_R2=read_dict_R2[clone_group_items[0]][0]
                                        if '-' in tseq_R2:
                                                tseq_R2 = "".join(c for c in tseq_R2 if c != '-') 
                                        umi_read_dict[umi_run_code ].append((clone_group_items[0],tseq, read_dict[clone_group_items[0]][1], tseq_R2, read_dict_R2[clone_group_items[0]][1]))
                                else:
                                        umi_read_dict[umi_run_code ].append((clone_group_items[0],tseq, read_dict[clone_group_items[0]][1]))
                                if save_assembly:
                                        print(">", file=asem_out) #separator between different conseus.
                                        if Paired and not Joined :
                                                print(">", file=asem_out_R2) #separator between different conseus.
        
                                if remove_singleton:
                                        continue
                                entry = "".join([">", clone_group_items[0], "\n", tseq, "\n"])
                                ret_list.append(entry)
                                if Paired and not Joined:
                                        entry=SeqRecord(Seq(tseq_R2), id=clone_group_items[0])
                                        entry=entry.reverse_complement(id=entry.id)
                                        entry = "".join([">", entry.id, "\n", str(entry.seq), "\n"])
                                        ret_list_R2.append(entry)
                                #continue, don't have to make consensus
                                
                                continue
                                
                        #in this group, we have more than 1 items, meaning no singleton
                        #now write the umi group (with each sequences/quality) to the umi_read_dict for output
                        umi_run_count = umi_run_count +1;
                        if umi_run_count >1:
                                umi_run_code ="".join([umi_code, "_", str(umi_run_count)])
                        for sub_group in clone_group_items:
                                tseq = read_dict[sub_group][0]
                                if '-' in tseq:
                                        tseq = "".join(c for c in tseq if c != '-') 
                                                                #umi_count_array[umi_run_code ]=len(clone_group_items)
                                if Paired and not Joined:
                                        tseq_R2 = read_dict_R2[sub_group][0]
                                        if '-' in tseq_R2:
                                                tseq_R2 = "".join(c for c in tseq_R2 if c != '-') 
                                        umi_read_dict[umi_run_code ].append((sub_group,tseq, read_dict[sub_group][1], tseq_R2, read_dict_R2[sub_group][1]))
                                else:
                                        umi_read_dict[umi_run_code ].append((sub_group,tseq, read_dict[sub_group][1]))

                        cons = make_consensus_2(read_dict, clone_group_items, False) 
                        cons_R2=None
                        if Paired and not Joined:
                                cons_R2=make_consensus_2(read_dict_R2, clone_group_items, False)
                        if len(cons) != debug_len:
                                print("fatal error, cons worng length!")
                                print("debug_len = ", debug_len)
                                print("conlen = ", len(cons))

                                print(">cons")
                                print(cons)
                                print(">entry")
                                print(seqlist[0][1])
                                sys.exit(1)
                        if save_assembly:
                                print(">consensus\n", cons, "\n", sep="", file=asem_out)
                                if Paired and not Joined:
                                        print(">consensus\n", cons_R2, "\n", sep="", file=asem_out_R2)
                                

                        if "-" in cons:
                                cons = ''.join(c for c in cons if c != '-')
                        if Paired and not Joined:
                                if "-" in cons_R2:
                                        cons_R2 = ''.join(c for c in cons_R2 if c != '-')
                                
                        #when working with consensus, we are taking any name sa the name
                        entry = "".join([">", clone_group_items[0], "(", str(len(clone_group_items)), ")" "\n", cons, "\n"])
                        ret_list.append(entry)
                        if Paired and not Joined:
                                #reverse comlement 
                                entry=SeqRecord(Seq(cons_R2), id="".join([ clone_group_items[0], "(", str(len(clone_group_items)), ")"]))
                                entry=entry.reverse_complement(id=entry.id)
                                ret_list_R2.append("".join([">", entry.id, "\n", str(entry.seq),"\n"]))
                        if save_assembly:
                                print(">", file=asem_out) #separator between different conseus.
                                if Paired and not Joined:
                                        print(">", file=asem_out_R2);
        #each + is one umi done.
        #print("+", end="", file=sys.stderr, flush=True) 
        #if PDEBUG:
        #       print("+", end="", flush=True) 
        asem_str = None
        asem_str_R2=None
        if save_assembly:
                asem_str = asem_out.getvalue()
                if Paired and not Joined:
                        asem_str_R2=asem_out_R2.getvalue()
        return (umi_code, ret_list, asem_str, umi_read_dict, ret_list_R2, asem_str_R2)

                         
def cloanalystPartition(rdict, g):
        #will not do cloanlyst , it is slow to do it.
        #print("g is:",g);
        if not run_Cloanalyst:
            return [g]
            
        #in here we call cloanalyst .dll to do cloanlyst
        #first, we need to reformat the reads
        rd_net=Dictionary[String,Array[String]]()
        for item in g:
                #get rid of '-' due to alignment
                tseq = rdict[item][0]
                if '-' in tseq:
                        tseq = "".join(c for c in tseq if c != '-')
                q=rdict[item][1]
                rd_net[item]=Array[String]([tseq,q])
                
        #now ready to call it
        id_groups=Interface.RunCloanalyst(rd_net)
        return id_groups
        
        
def find_umi(seq):
        if "CGATAGGGG" not in seq: return None
        umi_index = seq.index("CGATAGGGG")
        if umi_index == -1 or umi_index < 16:
                return None        
        umi_code = seq[umi_index-16:umi_index]
        return umi_code      


#reformat the dictionary, so that the umi is taking the first 15 bases of the read!
def reformat_read_dict(read_dict, seqlen = 15):
        newdict = defaultdict(list)
        newdict_R2=defaultdict(list)
        for umi in read_dict:
                readlist = read_dict[umi]
                for item in readlist:
                        prefix = item[1][:seqlen]
                        newkey = umi + prefix.lower()
                        newdict[newkey].append(item)
                        #if Paired and not Joined:
                        #        read_dict_R2[umi][
                        #for item in readlist:
                        #        newdict
                                
        return newdict

#the diciton contains {umi_key: list of tuple(read_name, seq, qual)}
#for dedup before assembly, return new read_dict
def umi_tools_dedup1(read_dict_src, seqext = 15):
        read_dict = reformat_read_dict(read_dict_src, seqext)
        bundles = defaultdict(dict)
        for umi in read_dict:
                #print("umi:", umi)
                umi_bytes_key = bytes(umi, 'utf-8') #umi_tool expect bytes.
                bundles[umi_bytes_key]["count"] = len(read_dict[umi])
                bundles[umi_bytes_key]["read"] = read_dict[umi][0]                
        threshold = 1
        if seqext > 0:
                threshold += 1
        
        #call dedup
        #processor = network.ReadDeduplicator()
        processor =network.UMIClusterer(cluster_method="directional") #directional, adjacent , cluster, percentile???need to confirm 
        #umi_clusters = processor.make_clusters(bundle=bundles, threshold= threshold)
        #reads,umis, umi_counts = processor(bundle=bundles, threshold= threshold)
        umis = bundles.keys()
        counts = {umi: bundles[umi]["count"] for umi in umis}

        umi_clusters = processor( counts, threshold)
        outdict = defaultdict(list)
        for cluster in umi_clusters:
                cluster_key = None
                master_umi = None
                for umi_bytes in cluster:
                        umi_key = umi_bytes.decode('utf-8')
                        front_key = umi_key[:16]
                        if cluster_key == None:
                                cluster_key = umi_key
                                master_umi = front_key
                                outdict[cluster_key].extend(read_dict[umi_key])
                        else:
                                #outdict[cluster_key].extend(read_dict[umi_key])
                                if Paired and not Joined:
                                        for (name, seq, qual, seq_r2, qual_r2) in read_dict[umi_key]:
                                                #mark read indicating this is from a separte umi
                                                if front_key != master_umi:
                                                        name += "_" + front_key
                                                outdict[cluster_key].append((name, seq, qual, seq_r2, qual_r2))

                                else:
                                        for (name, seq, qual) in read_dict[umi_key]:
                                                #mark read indicating this is from a separte umi
                                                if front_key != master_umi:
                                                        name += "_" + front_key
                                                
                                                outdict[cluster_key].append((name, seq, qual))


        #debug find whihc group key1 is in
        # target_group = -1
        # for i, cluster in enumerate(umi_clusters):
        #         for m in cluster:
        #                 if m == key1:
        #                         target_group = i
        #                         print("target group = ", i)
        #                         break
                

        # for item in umi_clusters[target_group]:
        #         key = item.decode()
        #         cnt = len(read_dict[key])
        #         print(key, cnt)

        show_info = True
        #dedup_log=open("dedup.log", "w")
        if show_info:
                for cluster in umi_clusters:
                        for umi_bytes in cluster:
                                umi_key = umi_bytes.decode('utf-8')
                                #find one read.
                                sample_read = read_dict[umi_key][0] #any read                
                                dedup_log.write(("umi="+umi_key+"\t{0}\t"+  sample_read[0]+ "\t"+ sample_read[1][:30]+"\n").format(
                                            len(read_dict[umi_key])))
                        dedup_log.write("-----------------------\n")
                #the returned Umis is supposedly the ones kept.
                dedup_log.write("read_dict_src umis: {0}\n".format( len(read_dict_src)))
                dedup_log.write("total number of clusters:{0}\n".format( len(umi_clusters)))

                #check splits? i.e. orinally one group and now becomes 2 or more gorup, i.e. in different cluster?
                splits = defaultdict(list)
                i = 0
                for cluster in umi_clusters:
                        for umi_bytes in cluster:
                                umi_key = umi_bytes.decode('utf-8')
                                umi5half = umi_key[:16]
                                splits[umi5half].append( (i, umi_key))
                        i += 1
                dedup_log.write("splits:\n")
                for item in splits:
                        index_set = set()
                        if len(splits[item]) == 1: 
                                continue
                        for index, key in splits[item]:
                                index_set.add(index)
                        if len(index_set) < 2:
                                continue
                        itemlist = splits[item]
                        newlist = itemlist.copy()
                        index, head = newlist.pop(0)
                        dedup_log.write(""+head+"{0}\t{1}\n".format( len(read_dict[head]), index))
                        for index, u in newlist:
                                cnt = len(read_dict[u])
                                du = ""
                                for i, c in enumerate(u):
                                        if c == head[i]:du += "."
                                        else: du += c
                                dedup_log.write("{0} \t{1}\t{2}\n".format(du, cnt, index))
                        dedup_log.write("\n")
        #dedup_log.close()                
        #sys.exit(1)
        return outdict        


#for dedup after assembly
#input consists of list of (umi_code, count, item)
#return list of items to keep
def umi_tools_dedup2(dedup_items):
        item_to_keep = []
        item_dict = {}
        bundles = defaultdict(dict)
        for umi, count, item in dedup_items:
                if umi in item_dict:
                        #not allowing duplicte umi to go through dedup
                        item_to_keep.append(item)
                        continue
                umi_bytes_key = bytes(umi, 'utf-8') #umi_tool expect bytes.
                bundles[umi_bytes_key]["count"] = count
                item_dict[umi] = item
                bundles[umi_bytes_key]["read"] = item

        threshold = 1
        if dedup_prefix_len > 0:
                threshold += 1
        
        #call dedup
        processor = network.ReadDeduplicator()
        umi_clusters = processor.make_clusters(bundle=bundles, threshold= threshold)
        print("output cluster count=", len(umi_clusters))

        #exclude singlets that "merging" into a master group.
        for cluster in umi_clusters:
                maxcount = 0
                for umi_bytes in cluster:
                       ct = bundles[umi_bytes]["count"]
                       if ct > maxcount:
                               maxcount = ct
                mincount = 2 if maxcount  > 1 else 1
                for umi_bytes in cluster:
                        ct = bundles[umi_bytes]["count"]
                        if ct >= mincount:
                                item_to_keep.append(bundles[umi_bytes]["read"])
        print("total items kept:", len(item_to_keep))
        return item_to_keep

'''
This is the function used to go through the input array, read_dict to get the umi_code freq array, meaning that
for each umi, how many sequence are associated with it.
Input : read_dict , dictionary, key=UMI, and value = (readID, seq, qual)
output: a umi_freq_dict or  umi_count_dict, 
'''
def generate_umi_counts_dictionary(read_dict):
    umi_freq_dict={}
    for umi in read_dict :
            umi_freq_dict[umi]=len(read_dict[umi])
    return umi_freq_dict
    
'''
    going through umi_freq_dict (which holds for each umi, how many sequences are associated with it)
    to generate umi_freq distribution, like how many umis are having one seq (singletons), how many have 2 seqs, etc.
    Input, umi_freq_dict,  umi: number of sequences 
    Output, dictionary, holding information about frequence of each umi length group,  (numberof umi sequence, frequency)
'''
def generate_umi_dist (umi_freq_dict):
        umi_dist={}
        for umi in umi_freq_dict:
            c = umi_freq_dict[umi]
            if c in umi_dist:
                        umi_dist[c] += 1
            else:
                        umi_dist[c] = 1
        return umi_dist
        
#inthe format of "read_name(5)"
def get_readcount_from_id(idstr):
        if idstr.find("(") == -1 or idstr.find(")") == -1:
                return 1
        ctstr = idstr[idstr.find("(")+1:idstr.find(")")]
        try:
                return int(ctstr)
        except ValueError:
                return 1
#this is to process each fasta file (most likely for each sample)
#scan file, for every entry, find umi tag,
#and the sequence, if sequence identical - keep one. 
#or drop the shorting ones?
#---update :  add Read2 files for no joining paired reads for processing.
# input : fastq file name for R1 and R2.
# output: processing log message.
#--- assume read1 and read2 are in same order.
def process_file(fastq_file, fastq_file_R2=None):
        read_dict = defaultdict(list)
        #read_dict_R2=  defaultdict(list)
        
        read_count = 0
        #item_to_keep = []
        print("parsing  Read 1:", fastq_file, flush=True)
        if Paired and not Joined:
                print("parsing Read2 :", fastq_file_R2, flush=True)
        processing_log=StringIO()
        
        #start reading the R2 files into a pool , so to process and draw it as we do for read1 
        #the reason we do this is that read2 doesn't have umi. we need to "wait" for read1 to get umi. so to speak "synchronize"
        temp_records_R2=defaultdict(list)
        #with open(fastq_file_R2, "r") as handle:
        if Paired and not Joined :
                for record in SeqIO.parse(fastq_file_R2,'fastq'):
                        read_count+=1
                        rc=record.reverse_complement(id=record.id)
                        qual=bytearray(list(map(operator.add, rc.letter_annotations['phred_quality'], [33]*len(rc.letter_annotations['phred_quality']))))
                        qual=qual.decode('ASCII')
                        seq_len=len(rc.seq)-right_align_trim
                        temp_records_R2[record.id.split(" ", 1)[0]]= (str(rc.seq)[:seq_len], qual[:seq_len])
        #for read 1 to read and also add information of R2 
        read_count=0
        with open(fastq_file, "r") as handle:
                for title, seq, qual in FastqGeneralIterator(handle):
                        read_count += 1
                        #print("read #:", read_count)
                        #print(title)
                        #print(seq)
                        #print(qual)
                        
                        #does this contains molecular bar code?
                        if "CGATAGGGG" not in seq: continue
                        umi_index = seq.index("CGATAGGGG")
                        
                        if umi_index == -1 or umi_index < 16:continue
                        umi_code = seq[umi_index-16:umi_index]
                       
                        #for debug
                        # if umi_code != "TCACCTAAAGGTTAGC":
                        #         continue
                        #print("umi_code = ", umi_code)
                        # if 'debug_item' in locals() and debug_item in title:
                        #         print("umi_code = ", umi_code)
                        #         sys.exit(1)

                        #count number of 'N's, if more than 1, discard.
                        nn = seq.count('N', umi_index+9)
                        if nn > max_N:
                                print("nN=", nn)
                                continue
                        #print("nN:", nn)
                        seq_len=len(seq)
                        if not (Paired and not Joined):
                                seq_len=len(seq)-right_align_trim
                        
                        myrecord=(title.split(" ", 1)[0], seq[umi_index+9:seq_len], qual[umi_index+9:seq_len])
                        if Paired and not Joined:
                                if title.split(" ",1)[0] not in temp_records_R2.keys():
                                        print("Error, the read ID is not in the read2 file")
                                        exit(-1)
                                myrecord=(title.split(" ", 1)[0], seq[umi_index+9:seq_len], qual[umi_index+9:seq_len],temp_records_R2[title.split(" ",1)[0]][0], temp_records_R2[title.split(" ",1)[0]][1]) 
                        
                        read_dict[umi_code].append(myrecord)
                        
                        #if Paired and not Joined:
                        #       read_dict_R2[umi_code].append((title.split(" ",1)[0],temp_records_R2[title.split(" ",1)[0]]))
                        #print(read_dict)
        print("length(read_dict):",len(read_dict))               
        #print("length(read_dict_r2):",len(read_dict_R2))               
        print("\tdone with reading files and initial umi processing!!!!")
        
        #print item to keep to file and extract to another file.
        #print(fastq_file, "number of read:", read_count, " read kept:", len(item_to_keep), " duplicate count: ", dup_count)
        print("\t\t",fastq_file, "number of read:", read_count, " number of (umi) bins kept:", len(read_dict))
        
        if len(read_dict) == 0:
                print("Warning: No reads kept for " + fastq_file)
                return 0
        
        #find umi entry with maximum read count, informal only...
        print("Start counting the reads in each umi bin(", fastq_file, ").......", flush=True)
        umi_count = {}
        i = 1
        for umi in read_dict:
                #print (">item" + str(i), "count=", len(read_dict[umi]))
                i += 1
                #print(umi)
                umi_count[umi] = len(read_dict[umi])
        if len(umi_count) == 0:
                print("Warning: No umi_count ?? kept for " + fastq_file)
                return 0
        maximum = max(umi_count, key=umi_count.get)
        print("\tmaximum umi count(", fastq_file, "):", maximum, umi_count[maximum], flush =True)
        
        #check single read UMIs which is one edit distance away from a multi-reads umi.
        #checking very similar umis
        check_info = False
        if check_info:
                hb_UMIs = []    #high aboundance umis.
                for item in umi_count:
                        if umi_count[item] > 100:
                                hb_UMIs.append(item)
                #print("high aboundance UMIs:", hb_UMIs)

                umi_err_list = []      #item with 1 difference from a high aboundance one.
                for umi in umi_count:
                        for  hb_umi in hb_UMIs:
                                diff = sum(1 for a, b in zip(umi, hb_umi) if a != b)
                                if diff == 1:
                                        umi_err_list.append(umi)
                                        break
                
                print("Total number of UMIs", len(umi_count))
                print("Single read UMIs similar to a high aboundance ones:", len(umi_err_list))
                print("if remove those single read UMI, we will have left:", len(umi_count) - len(umi_err_list))
                #print(umi_err_list)
                
        #stats - one histogram we can do is 
        # count how many UMIs has one read, how many has 2 reads.
        # i.e. make a dictiony of (int, int) -> (read_count, num_of_umis)
        print("------------------------\n")
        print("UMI count stats -key umi_count, value- the umi have this many reads(", fastq_file, ").\n")
        umi_dist = {}           #dictionary of <int, int>, say we have an entry as (3, 10)                                 #means there are 10 different UMIs, each has 3 reads.
        for umi in umi_count:
                c = umi_count[umi]
                if c in umi_dist:
                        umi_dist[c] += 1
                else:
                        umi_dist[c] = 1
        print("\t(", fastq_file, ")",umi_dist, flush=True)
        #print("
        
        show_histogram = False
        if (show_histogram):
                xax = [ str(i) for i in list(umi_dist.keys())]
                plt.bar(xax, umi_dist.values(), color='g')
                plt.show()
        
        processing_log.write("{0}\tbefore UMI dedup\t{1}\n".format(fastq_file,umi_dist))
        #exit(-1)
        #if perform dedup before assembly
        
        if dedup_option == 1:
                print("Start doing UMI dedup(", fastq_file, ").........", flush=True)
                read_dict = umi_tools_dedup1(read_dict, 20) #in here, we add extra 20 nt to the umi to take care of the <------
                read_freq_dict_dedup=generate_umi_counts_dictionary(read_dict)
                read_count_dist_dedup=generate_umi_dist(read_freq_dict_dedup)
                print("\tDedup Summary(", fastq_file, "):", read_count_dist_dedup)
                print("\tdone with dedup(", fastq_file, ")!!", flush=True)
        #if PDEBUG:
        #        print("dedup option:", dedup_option)
        #        print("remove_singleton:", remove_singleton)
        #the output
        item_to_keep = [] 
        item_to_keep_R2=[]
        singleton_to_keep=[]
        singleton_to_keep_R2=[]
        #if "post" dedup needed.
        dedup_items = []  ###dedup _option ==2 not implemented !!! be careful
        
        #print to log the umi distribution
        umi_count={}    #<-not necessary, but make it clear
        umi_dist={}   #<-not neccessary, but make it clear 
        umi_count =generate_umi_counts_dictionary(read_dict)
        umi_dist=generate_umi_dist(umi_count)
        processing_log.write("{0}\tafter UMI dedup\t{1}\n".format(fastq_file,umi_dist))
        #process singleton."1", 
        
        #pick out singltons so that we can process umis of non-singlton the next step 
        #print("before picking out singlton, len of read_dict:", len(read_dict))
        for umi in read_dict:
                entry = read_dict[umi]
                if len(entry) != 1: continue  #for non singleton entries, we do nothing in here. we only 
                                                                                                #go ahead to take care of singletns
                
                #print("entry:", entry)
                #print("entry[0]:", entry[0])
                #for singletons, based one remove_singleton to do different things
                item = "".join([">", entry[0][0], "\n", entry[0][1], "\n"])  #seem to rewrite it as an a fasta item with name and seq
                if Paired and not Joined: #take care of read 2
                        entry=SeqRecord(Seq(entry[0][3]), id=entry[0][0])
                        entry=entry.reverse_complement(id=entry.id)
                        item_R2="".join([">", entry.id, "\n", str(entry.seq), "\n"])
                        
                if not remove_singleton: #pick the singleton and write directly to the "output" array in fasta format 
                        item_to_keep.append(item)
                        if Paired and not Joined:
                                item_to_keep_R2.append(item_R2)
                        if dedup_option == 2: #this is the option do the cloanalyst first and then dedup. so we want to put singleton in for dedup, but no cloanalyst
                                ext_umi_code = umi + entry[0][1][0:dedup_prefix_len].lower()
                                dedup_items.append((ext_umi_code, 1, item)) #umi, count, cons.
                else :  #for removing singltone, we and them to a different buffer. 
                        singleton_to_keep.append(item)
                        if Paired and not Joined:
                                singleton_to_keep_R2.append(item_R2)
                    
        #remove singleton which already handled in the previous step                       
        read_dict = {k : v for k, v in read_dict.items() if len(v) > 1 }
        #print("after picking out singlton, len of read_dict:", len(read_dict))
        #print("**************items to keep  length", len(item_to_keep), "\n")
        #print("**************items to keep  length", len(singleton_to_keep), "\n", flush=True)
        if save_assembly:  #get directory ready for new writing.
                tname = os.path.basename(fastq_file)
                pname = os.path.splitext(tname)[0]
                if not os.path.exists(pname):
                        os.mkdir(pname)
                else:   #if directory exists, then remove the fasta files inside.
                        filelist = [f for f in os.listdir(pname) if f.endswith(".fasta")]
                        for f in filelist:
                                os.remove(os.path.join(pname, f))        

        print("Start doing processing umi groups(", fastq_file, ").........")
        print("\t total number of umi groups  (non-singltons)(", fastq_file, "):", len(read_dict), flush=True)
        dedup_src = []
        parallelize = False #parallel_level == 1
        xx = 0
        umi_count={}
        umi_read_dict_all=defaultdict(list) # for doing the overall umi distribution and singleton to write 
        len_singleton=len(item_to_keep) 
        if remove_singleton:
            len_singleton=len(singleton_to_keep)
        #print("umi_read_dict:", umi_read_dict_all)
        
        if not parallelize:
                #now choose sequences from each umi...
                for umi in read_dict:
                        if PDEBUG:
                                print("#########finding the candidate sequences for umi:", umi , " for ", len(read_dict[umi]), " sequences",flush=True) 
                        #if umi != "GTAATTACCAGATTCG":
                        #        continue
                        seqlist = read_dict[umi]
                        #if len(seqlist)< 2000 or len(seqlist) > 5000:continue
                        if (len(seqlist) > 1000):
                                print("finding candiate seq. for ", len(seqlist), " sequences", file=sys.stderr)
                        #6/1/2020, now we add more output for R2 as output. Note: umi_read_dict itself contains more about R2, 
                        #we don't need to have separated output for read 2.
                        (_, keep_list, asm_txt, umi_read_dict, keep_list_R2, asm_txt_R2) = find_candidate_sequences((umi, seqlist))
                        umi_read_dict_all.update(umi_read_dict)
                        #print("**umi_read_dict:",umi_read_dict_all, flush=True)
                        if PDEBUG:
                                print("after find_candidate_sequences, show keep_list:", keep_list)
                        if keep_list is None: continue
                        
                        if dedup_option == 2:
                                for item in keep_list:
                                        #each item is a faster entry...
                                        fasta_io = StringIO(item)
                                        for id, seq in SimpleFastaParser(fasta_io):
                                                ext_umi_code = umi + seq[:dedup_prefix_len].lower()
                                                cnt = get_readcount_from_id(id)
                                                dedup_items.append((ext_umi_code, cnt, item))                        
                        if save_assembly and asm_txt is not None:  #this means asm_txt_R2 is not None.  asm_txt and asm_txt_R2 have same elements.
                                asem_path = os.path.join(pname, umi+".fasta")
                                with open(asem_path, "w") as asmf:
                                        asmf.write(asm_txt)
                                if Paired and not Joined and asm_txt_R2 is not None:  #<-----------none???
                                        asem_path = os.path.join(pname, umi+"_R2.fasta")
                                        with open(asem_path, "w") as asmf:
                                                asmf.write(asm_txt_R2)
                        item_to_keep.extend(keep_list)
                        if Paired and not Joined:
                                item_to_keep_R2.extend(keep_list_R2)
                        xx += 1
                        if xx %2500==0: print("Progress (", fastq_file, "): ", xx ,"/", len(read_dict),"....", flush=True)
                        if PDEBUG:
                                print("*****Done for this round !!!!!!!!!!!", flush=True)
                
        else: #do in parallel
                pool = Pool()
                result_list = pool.map(find_candidate_sequences, read_dict.items())
                for umi_code, keep_list, asm_txt , umi_read_dict, keep_list_R2, asm_txt_R2 in result_list:
                        item_to_keep.extend(keep_list)
                        if Paired and not Joined:
                                item_to_keep_R2.extend(keep_list_R2)
                        umi_read_dict_all.update(umi_read_dict)
                        if dedup_option == 2:
                                for item in keep_list:
                                        #each item is a fasta entry...
                                        fasta_io = StringIO(item)
                                        for id, seq in SimpleFastaParser(fasta_io):
                                                ext_umi_code = umi_code + seq[:dedup_prefix_len].lower()
                                                cnt = get_readcount_from_id(id)
                                                dedup_items.append((ext_umi_code, cnt, item))  
                        if save_assembly:
                                asem_path = os.path.join(pname, umi_code+".fasta")
                                with open(asem_path, "w") as asmf:
                                        asmf.write(asm_txt)
                                if Paired and not Joined:
                                        asem_path = os.path.join(pname, umi+"_R2.fasta")
                                        with open(asem_path, "w") as asmf:
                                                asmf.write(asm_txt_R2)
                pool.close()
                pool.join()

        
        if dedup_option == 2:
                item_to_keep = umi_tools_dedup2(dedup_items)
        umi_count=generate_umi_counts_dictionary(umi_read_dict_all)
        umi_dist=generate_umi_dist(umi_count)
        #print("-----------------checking the distribution ----------------")
        #print(umi_read_dict_all)
        if 1 in umi_dist:
            umi_dist[1]=umi_dist[1]+len_singleton
        else:
            umi_dist[1]=len_singleton
            
        print("\numi_dist(", fastq_file, "):", umi_dist)
        processing_log.write("{0}\tafter partition\t{1}\n".format(fastq_file,umi_dist))
        
        #<--------------check for remove_singleton and then write the singltone to the file.
        #
        for   umi_item in umi_read_dict_all:
                if len(umi_read_dict_all[umi_item])==1 : #singleton again
                        if remove_singleton:
                                entry=umi_read_dict_all[umi_item]
                                item = "".join([">", entry[0][0], "\n", entry[0][1], "\n"])  #seem to make it a fasta item with name and seq
                                singleton_to_keep.append(item)
                                if Paired and not Joined : #R2 case
                                        entry=SeqRecord(Seq(entry[0][3], generic_dna),id=entry[0][0])
                                        entry=entry.reverse_complement(id=entry.id)
                                
                                        item ="".join([">", entry.id, "\n", str(entry.seq), "\n"])
                                        singleton_to_keep_R2.append(item)
                        #in this case, we don't do anything for not remove_singleton, since it has been done in find_candidate_sequences() function
                
        print("\n", flush=True)
        print("**processing summary:")
        print("****file name:\"",fastq_file)
        if Paired and not Joined:
                 print("\t:\"", fastq_file_R2,"\"")
        print ("****total number of read processed:", read_count)
        print("****number of bins kept (without singleton):", len(read_dict))
        tempStr="****# of reads kept [2]"
        if not remove_singleton :
                tempStr= tempStr +" (including singletons):"
        else :
                tempStr=tempStr+" (without singletons):"
        print(tempStr, len(item_to_keep))
        print("******Note: the number of reads kept [2] might not be the same as the total number of bin kept with singleton[1],")
        print("******\tsince the number of bin is identical to the number of unique barcodes and each barcode might be split")
        print("******\tinto different groups of read. When we do umi_tool dedup, things are more complicated.");
        #exit(-1)

        outfile = ""
        #updated by feng, always save at the current folder.
        fasta_file_name = fastq_file.replace("fastq", "fasta")
        fasta_file_name = os.path.basename(fasta_file_name)
        #if fastq_file.startswith("../"): 
        #        fasta_file_name = fastq_file.replace("fastq", "fasta")[1:]
        outfile = fasta_file_name
        with open(fasta_file_name, "w") as faf:
                faf.writelines(item_to_keep)
        if Paired and not Joined: #R2
                fasta_file_name = fastq_file_R2.replace("fastq", "fasta")
                fasta_file_name = os.path.basename(fasta_file_name)
                #if fastq_file.startswith("../"): 
                #        fasta_file_name = fastq_file.replace("fastq", "fasta")[1:]
                outfile = fasta_file_name
                with open(fasta_file_name, "w") as faf:
                        faf.writelines(item_to_keep_R2)

        if remove_singleton:
                singleton_file = fastq_file.replace("fastq", "singleton.fasta")
                singleton_file= os.path.basename(singleton_file)
                with open(singleton_file, "w") as faf:
                        faf.writelines(singleton_to_keep)
                if Paired and not Joined:
                        singleton_file = fastq_file_R2.replace("fastq", "singleton.fasta")
                        singleton_file= os.path.basename(singleton_file)
                        with open(singleton_file, "w") as faf:
                                faf.writelines(singleton_to_keep_R2)
        print("Done!!! ", flush=True)                    
        return  processing_log.getvalue() #not running cloanalyst for testing.
        cmd2 = "run_cloanalyst_mouse_umi.sh " + outfile #fastq or fasta file depending on if using consensus
        return_code = subprocess.check_call(cmd2, shell=True)
        if return_code != 0:
                print("calling cloanalyst failed!")


###############################################################                
###this is the main function, normal running of the program
dedup_log=open("dedup"+(datetime.now().strftime("%Y%m%d%H%M%S"))+".log","w")
processing_log_file=open("processing"+(datetime.now().strftime("%Y%m%d%H%M%S"))+".log","w")
processing_log_file.write("sample\tdataSource\tUMI_dist\n")
if parallel_level == 1:
        (files, files_R2)=project.get_input_files(sample_list)
        #print("file list: ",files)
        #sys.exit(0)
        
        for i in range(0, len(files)): #range(13, 25): #79
                #if i < 13:continue
                #sample_name = "Sample" + str(i).rjust(2, '0')
                #infile = "./" + sample_name + "_trimmed.fastq"
                infile=files[i]
                
                infile_R2=None
                if Paired and not Joined :
                        infile_R2=files_R2[i]
                        #print("Read2 file :", infile_R2)
                        if not os.path.exists(infile_R2):
                                print("file ", infile_R2, "missing")
                                continue
                #print(infile);
                if not os.path.exists(infile):
                        print("file ", infile, "missing")
                        continue
                #print("checking file", infile)
                #sys.exit(-1)
                dedup_log.write("writing log for {0}......\n".format(infile))
                #exit(-1);
                
                start_time = time.time()
                pl=process_file(infile, infile_R2)
                end_time = time.time()
                #index=index+1
                processing_log_file.write(pl)
                print("done processing", infile)
                print("--- took %s seconds ---" %(end_time - start_time))
                dedup_log.write("\nDONE!\n")
                dedup_log.write("-------------------------------------------------------------------------\n");
else:   #do this only if runnning all projects
        file_list = []
        file_list_R2=[]
        
        (files, file_R2)=project.get_input_files(sample_list)
        #print(files, flush = True )
        #sys.exit(0)
        if len(files)==0 :
                print("Error: can not find the specified files, please double check!!!")
                sys.exit(-1)
        if Paired and not Joined :
                if len(file_R2)==0 :
                        print("Error: can not find the specified Read 2 files, please double check!!!")
                        sys.exit(-1)
                if len(files)!=len(file_R2):
                        print("ERROR: Read 1 and Read 2 files are paired. please double check ")
                        sys.exit(-1)    

                        
        for i in range(0,len(files)): #range(13, 25):ls run
                #if i < 4 or i >= 13 : continue
                #this is the joined/mapped/trimmed file.
                #sample_name = "Sample" + str(i) #str(i).rjust(2, '0')
                #sample_name = project.get_input_files(i)
                #infile = "./" + sample_name + ".fastq"
                infile =files[i]
                if not os.path.exists(infile):
                        print("file ", infile, "missing")
                        continue
                
                infile_R2=None
                if Paired and not Joined :
                        infile_R2=file_R2[i]
                        if not os.path.exists(infile_R2):
                                print("file ", infile, "missing")
                                continue
                        file_list_R2.append(infile_R2)
                file_list.append(infile)
                
                #now build a input, input is a list or list of tuple of two 
                inputs=file_list
                if Paired and not Joined:
                        inputs=list(map(lambda x, y:(x,y), file_list, file_list_R2)) 
                
        #print(file_list);
        #print(file_list_R2);
        print(inputs)
        
        pool = Pool()
        if Paired and not Joined:
                pl_list=pool.starmap(process_file, inputs)       
        else:
                pl_list=pool.map(process_file, inputs)
                
        for pl in pl_list:
                processing_log_file.write(pl)
dedup_log.close()
processing_log_file.close()
