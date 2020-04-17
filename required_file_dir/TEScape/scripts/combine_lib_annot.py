#!/usr/bin/env python

##This script will annotate no mite TE in the TE libraries
##This step will use the output dic from the previous step combine_lib

import subprocess
from Bio import SeqIO
from Bio.SeqUtils import six_frame_translations
import re
import sys
import os
import pandas as pd


##define a function to translate each sequence
##intial a translation function
##only translate the all the six reading frames for the no mite DNA TEs

##define a reverse complement functions
def getRC(te_seq_rc):
    out = te_seq_rc.reverse_complement()
    return out

##define a function to translate
def translate (te_seq):
    ##initial a dictionary to store the protein sequence for the six reading frame
    pro_seq_dic = {}
    ##Runs a six-frame translation for each gene in the string ts2seq
    f6 = six_frame_translations(te_seq)
    ##Splits each gene in the above mentioned translation with a new line
    f6tmp = str(f6).split('\n')
    num_m = re.match('(.+)5$', str(len(f6tmp)))
    num_cnt = num_m.group(1)
    ####################################################################
    ##the following is the first three reading frame of this te sequence
    ##get the location information of sequences
    l1 = 5
    l2 = 6
    l3 = 7
    ##Initializes a hand full of lists we will later fill with the genomic data to analyze
    list1 = []
    list2 = []
    list3 = []
    ##Sets range that the genomic data we are looking at can be based on the 6-frame translation
    for i in range(0, (int(num_cnt))):
        l1 = 5 + i * 10
        list1.append(l1)
        l2 = 6 + i * 10
        list2.append(l2)
        l3 = 7 + i * 10
        list3.append(l3)
    ##extract the sequence information from the location informationa
    ##the following is the first three list
    allist_f = [list1, list2, list3]
    ##create a string to be the index for the key of the dictionary
    fr_count = 0
    for eachlist in allist_f:
        fr_count = fr_count + 1
        ##extrac the sequence information
        ##For each gene, adds the word 'list' for each additional gene after the start + end codon
        framenm = 'f' + str(fr_count)
        ##initial a list to store all the protein sequences
        f6aa = []
        # Setting the frame protein for each translation
        for i in eachlist:
            pseq = f6tmp[i]
            f6aa.append(pseq)
        ##combination of the sequence
        comptlist = ' '.join(f6aa)
        ##split the '' between each aa and then replace it with ,
        comptstr = comptlist.replace(' ', '')
        ##store sequence information on the dictionary
        pro_seq_dic[framenm] = comptstr

    ###############################################################
    ##the following is last three reading frame of this te sequence

    ##Inputs the data to get the reversed compliment
    seq_reverse = getRC(te_seq)
    f6 = six_frame_translations(seq_reverse)
    f6tmp = str(f6).split('\n')
    num_m = re.match('(.+)5$', str(len(f6tmp)))
    num_cnt = num_m.group(1)
    l4 = 5
    l5 = 6
    l6 = 7
    list4 = []
    list5 = []
    list6 = []
    for i in range(0, (int(num_cnt))):
        l4 = 5 + i * 10
        list4.append(l4)
        l5 = 6 + i * 10
        list5.append(l5)
        l6 = 7 + i * 10
        list6.append(l6)
    allist_f = [list4, list5, list6]
    re_count = 3
    for eachlist in allist_f:
        re_count = re_count + 1
        framenm = 'f' + str(re_count)
        f6aa = []
        for i in eachlist:
            pseq = f6tmp[i]
            f6aa.append(pseq)
        comptlist = ' '.join(f6aa)
        comptstr = comptlist.replace(' ', '')
        pro_seq_dic[framenm] = comptstr

    return (pro_seq_dic)


##initial a function to identify the domain for the six reading frame
##define a function to identify te domain and generate the domain table
##import the pro_seq_dic generated from the translate function
def iden_domain_six_rdframe (pro_seq_dic,seqid,hmmscan_exe,working_dir):

    ##initial a dictionary to store the results
    te_domain_six_rdframe_dic = {}

    ##initial a dictionary to store the name of six reading frame
    six_frame_nm_dic = {}
    for eachframe in pro_seq_dic:
        te_frame_nm = seqid + '_' + eachframe
        te_opt_nm = eachframe + '.seq'
        te_pro_seq = open(working_dir + '/' + te_opt_nm, 'w')
        ##this will generate each reading frame file
        te_pro_seq.write('>' + te_frame_nm + '\n' + str(pro_seq_dic[eachframe]))
        te_pro_seq.close()  ##it is important to close the te_pro_seq file
        six_frame_nm_dic[eachframe] = 1

    ##use the hmm to detect domain information for each frame
    for eachnm in six_frame_nm_dic:
        tefm_seq = eachnm + '.seq'
        opt_tefm = eachnm + '.out'
        cmd = hmmscan_exe + ' ' + working_dir + '/minifam ' + working_dir + '/' + tefm_seq + ' > ' + working_dir + '/' + opt_tefm
        subprocess.call(cmd, shell=True)

        ##open the output file from hmmscan
        domain_te_file = open(working_dir + '/' + opt_tefm, 'r')
        for line in domain_te_file:
            line.strip()
            if re.match('.+\d+\s+PF\d+_.+', line):
                mt = re.match('(.+PF\d+_[A-Z]+).+', line)
                new_line = mt.group(1)
                te_line = seqid + '_' + eachnm + '\t' + new_line
                te_domain_six_rdframe_dic[te_line] = 1

    return (te_domain_six_rdframe_dic)


##define a function to store the te name and its domain number information in the dictionary
##the following funcitons is to filter the best reading frame
def cal_dm_num (te_nm,domain_dtf): ##note: the domain_dtf is generated by the next function.
    ##creat a dictionary to store number information for each domain
    domain_num_dic = {}
    ##collect the target te
    ##filter the different domain name for the dna_tr
    te_dt = domain_dtf.loc[domain_dtf[domain_dtf.columns[0]] == te_nm]
    te_dna_tr_dt = te_dt[te_dt[te_dt.columns[9]].str.contains('TR', flags=re.I, na=False)]
    domain_num_dic[te_nm] = {'TR': len(te_dna_tr_dt)}

    return domain_num_dic

##define a function to calculate the number of te_type
##this function will be used in the compre_rdframe
def cal_num (te_type):
    if len(str(te_type)) != 0:
        num = te_type
    else:
        num = 0
    return (num)

##intial a function to compare each reading frame and choose the best one
##the second output of the compare_rdframe function is to store the te domain number information
##import the pro_seq_dic from the translate function
def compare_rdframe (te_domain_six_rdframe_dic,working_dir):
    ##the import dic is the te_domain_six_rdframe_dic from the previous iden_domain_six_rdframe

    #########################################################################
    ##Step 1: write out because the next step will use the output of this dic
    with open(working_dir + '/opt_temp_six_rdframe_domain_hmm.tb', 'w+') as opt_te_ltr:
        for eachline in te_domain_six_rdframe_dic:
            opt_te_ltr.write(eachline + '\n')
    ##read the temp_domain_hmm table
    domain_df = pd.read_table(working_dir + '/opt_temp_six_rdframe_domain_hmm.tb', header=None, delimiter=r"\s+")
    te_nm = domain_df[domain_df.columns[0]].unique() ##contain information for six reading frame name
    ##the df dataframe only store te domain information for 0 or 1 that are used for filter the longest reading frame
    df = pd.DataFrame(index=te_nm, columns=['TR_DNA'])
    ##transfer the domain_df to the domain_dtf
    domain_dtf = pd.DataFrame(domain_df)

    ##get the reduced name
    reduced_nm = ''
    ################################################
    ##Step 2: create the table of domain information
    for eachte in te_nm:
        ##call the domain dictionary from function cal_dm_num
        te_domain_num_dic = cal_dm_num(eachte,domain_dtf)
        ##save the number of domain to the table

        ##for the nomit####################
        if te_domain_num_dic[eachte]['TR'] != 0:
            df.loc[eachte,'TR_DNA'] = 1
        else:
            df.loc[eachte, 'TR_DNA'] = 0

        ##get the reduced_nm
        mt = re.match('(.+)_f\d', eachte)
        reduced_nm = mt.group(1)


    ###############################################################################
    ##step 3: analyze each te with all reading frames and get the filtered te table
    ##create an empty dataframe to store the filtered information
    #df_ft = pd.DataFrame(index=set(te_nm_list), columns=['RT_LTR', 'IN_LTR', 'GAG_LTR', 'RH_LTR'])
    ##because df only contains six dataframe so do not need to create a new dt to store only dt for one dataframe
    max_te_nm = ''

    dna_tr = ''

    df['sum'] = df[['TR_DNA']].sum(axis=1)
    ##only select the best te
    te_dna_tr_dt = df[df['TR_DNA'] != 0]
    if (len(te_dna_tr_dt) != 0):
        max_te_nm = te_dna_tr_dt['sum'].idxmax()
        dna_tr = df.loc[max_te_nm, 'TR_DNA']
    else:
        max_te_nm = df['sum'].idxmax()
        dna_tr = df.loc[max_te_nm, 'TR_DNA']

    #print ('the dna_tr is ' + str(dna_tr))
    #########################################################
    ##step 4: store the number information for each max_te_nm
    ##get the reduced nm
    mt = re.match('(.+)_f\d',max_te_nm)
    domain_reduced_nm = mt.group(1)
    domain_num_dic = {} ##this dic store the number for the max_te_nm
    domain_num_dic[domain_reduced_nm] = {'TR': cal_num(dna_tr)}

    return(domain_num_dic)










