#!/usr/bin/env python

##this script is to generate the summary of the te in the genome
##te_family_count te_super_family_count proportion size fusion or not

import re
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt

#####################################################
##initiate a function to get the summary of the table
##te family count/size/proportion
##te family fusion/single/Unclear count/size/proportion
##generate the score distribution table/figure
##      fam1        fam2      fam3
##0-10  number/pro
##20-30
##30-40




def store_all_te (input_te_final_file):
    all_te_dic = {}
    with open (input_te_final_file,'r') as ipt_all_te:
        for eachline in ipt_all_te:
            all_te_dic[eachline] = 1
    return (all_te_dic)

##these following two are used for the te_summary function
def store_type (fam_nm,size,type_dic):
    if fam_nm in type_dic.keys():
        type_dic[fam_nm]['count'] += 1
        type_dic[fam_nm]['size'] += size
    else:
        type_dic[fam_nm] = {'count': 1, 'size': size}

def generate_final_line (store_te_dic,total_size):
    summary_te_dic = {}
    for eachfam in store_te_dic:
        te_pro = store_te_dic[eachfam]['size']/int(total_size)
        final_line = eachfam + '\t' + str(store_te_dic[eachfam]['count']) + '\t' + \
                     str(store_te_dic[eachfam]['size']) + '\t' + str(te_pro)
        summary_te_dic[final_line] = 1
    return (summary_te_dic)

def te_summary (all_te_dic,input_genome_fas):

    ##initiate a dic to store count and size
    store_all_te_dic ={}
    store_fusion_te_dic = {}
    store_single_te_dic = {}
    store_unclear_te_dic = {}

    for eachline in all_te_dic:
        eachline = eachline.strip('\n')
        col = eachline.strip().split()

        ##get the family name and size information
        mt = re.match('.+\/\/(.+)_TE.+',col[3])
        fam_nm = mt.group(1)
        size = int(col[2]) - int(col[1])
        type = col[8]

        store_type(fam_nm,size,store_all_te_dic)

        if 'Unclear' in type:
            store_type(fam_nm,size,store_unclear_te_dic)

        if 'Fusion' in type:
            store_type(fam_nm,size,store_fusion_te_dic)

        if type == 'Single':
            store_type(fam_nm,size,store_single_te_dic)


    ##calcualte the bp in the input_genome_fas
    total_size = 0
    for seq_record in SeqIO.parse(input_genome_fas,'fasta'):
        chr_size = len(str(seq_record.seq))
        total_size = total_size + chr_size

    ##calculate the proporiton of the genome
    all_summary_te_dic = generate_final_line(store_all_te_dic, total_size)
    ##fusion
    fusion_summary_te_dic = generate_final_line(store_fusion_te_dic, total_size)
    ##single
    single_summary_te_dic = generate_final_line(store_single_te_dic, total_size)
    ##unclear
    unclear_summary_te_dic = generate_final_line(store_unclear_te_dic, total_size)

    return (all_summary_te_dic,fusion_summary_te_dic,single_summary_te_dic,unclear_summary_te_dic)

##draw score distribution, generate histgrame figure
def draw_score_distribution (input_te_final_file,output_dir):
    df = pd.read_table(input_te_final_file, header=None, delimiter=r"\s+")
    df.hist(column=14,bins=10,color='#86bf91',grid=False,figsize=(10,5), layout=(1,1),sharex=True,zorder=0, rwidth=0.7)
    plt.xlabel("Score (0-1)")
    plt.ylabel("TE number\n")
    plt.title('Score Distribution')
    plt.savefig(output_dir + '/opt_score_distribution.png')






































