#!/usr/bin/env python

##this script is to identify if the RT is complete or not, and what the type of this LTR: RLG, RLC, RLX
##predict the time of LTRs
##generate a complete table to summary all the information

import sys
import re
from Bio import SeqIO
import subprocess
import math

####################
##set the input file
####################

######################################################################################
##initial a function to generate the complete table without divergence time prediction
def generate_compt_tb (input_ltr_pat_file,input_ltr_loc_file):

    ##return the dic storing all the information
    compt_tb_dic = {}

    ##store information of input_ltr_pat infor
    ##classify the ltr types
    ##this dic contains pat and type two reference key
    ltr_pat_dic = {}
    with open (input_ltr_pat_file,'r') as ipt_pat:
        for eachline in ipt_pat:
            col = eachline.strip().split()
            ##check eachline and filter out eachline containing one item
            tenm = col[0]

            ##make sure col[1] can work
            if len(col) > 1:

                pat = col[1]

                ##classify the te type
                rt_mt = re.compile('RT')
                rt_mt_list = rt_mt.findall(pat)

                ##only contain one RT
                ##this situation will divide several types:
                ##1: RT/RH - INT
                ##2:
                ##3:
                if len(rt_mt_list) == 1:

                    in_mt = re.compile('IN')
                    in_mt_list = in_mt.findall(pat)

                    ##if the IN number is not equal to the 0
                    if (in_mt_list) != 0:

                        rh_mt = re.compile('RH')
                        rh_mt_list = rh_mt.findall(pat)
                        if len(rh_mt_list) != 0:
                            ##if match the RLG_C: RTRH+IN+
                            if re.match('^RT[RH]+[IN]*IN$',pat):
                                ltr_pat_dic[tenm] = {'pattern':pat,'type':'RLG_C'}

                            ##if match the RLG_C: RH+RTIN+
                            if re.match('^[RH]+RT[IN]*IN$',pat):
                                ltr_pat_dic[tenm] = {'pattern': pat, 'type': 'RLG_C'}

                            ##if match the RLC_C: IN+RTRH+
                            if re.match('^[IN]+RT[RH]*RH$',pat):
                                ltr_pat_dic[tenm] = {'pattern': pat, 'type': 'RLC_C'}

                            ##if match the RLC_C: IN+RH+RT
                            if re.match('^[IN]+[RH]+RT$',pat):
                                ltr_pat_dic[tenm] = {'pattern': pat, 'type': 'RLC_C'}

                        if len(rh_mt_list) == 0:
                            ##if match the RLG_N: RTIN+
                            if re.match('^RT[IN]*IN$',pat):
                                ltr_pat_dic[tenm] = {'pattern': pat, 'type': 'RLG_N'}

                            ##if match the RLC_N:
                            if re.match('^[IN]+RT$',pat):
                                ltr_pat_dic[tenm] = {'pattern': pat, 'type': 'RLC_N'}

                    ##if the IN number is equal to the 0, it will be confused the type of LTR and marked as N_C
                    if (in_mt_list) == 0:
                        ltr_pat_dic[tenm] = {'pattern': pat, 'type': 'NC'}

                ##if there is not RT, it will be regarded as the N_C
                if len(rt_mt_list) == 0:
                    ltr_pat_dic[tenm] = {'pattern': pat, 'type': 'NC'}

                ##if there are more than two RTs, it will be regarded as the 'FU' (fused)
                if len(rt_mt_list) > 1:
                    ltr_pat_dic[tenm] = {'pattern': pat, 'type': 'FU'}

            if len(col) == 1:
                ltr_pat_dic[tenm] = {'pattern': 'NA', 'type': 'RLX'}


    ##open the input_ltr_loc_file to generate the compplete table
    with open (input_ltr_loc_file,'r') as ipt_loc:
        for eachline in ipt_loc:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            ltrnm = col[1]
            if ltrnm in ltr_pat_dic.keys():
                final_line = eachline + '\t' + ltr_pat_dic[ltrnm]['pattern'] + '\t' + ltr_pat_dic[ltrnm]['type']
                compt_tb_dic[final_line] = 1
            else:
                final_line = eachline + '\t' + 'NA' + '\t' + 'RLX'
                compt_tb_dic[final_line] = 1


    return(compt_tb_dic)


#####################################################################
##initial a function to extract target sequence and generate the file
def extr_ltr_seq_dect_time (compete_tb_dic,input_genome_fas_file,bedtools_exe,muscle_exe,ltr_dir):
    ##import the compete_tb_dic to get the ltr location information
    ##use the location information to extract the sequence from ltr_fas_file
    ##create the bed file for each ltrname
    ##initial a dic to store the compete_tb_dic
    compete_time_tb_dic = {}
    count = 0
    for eachltr in compete_tb_dic:
        count = count + 1
        print ('the analyzed te number is ' + str(count))
        col = eachltr.strip().split()
        ##the bedline contains two line the left and right ltr
        bedline = col[0] + '\t' + col[5] + '\t' + col[6] + '\t' + col[1] + '_1' + '\t' + '1' + '\t' + col[4] + \
                  '\n' + col[0] + '\t' + col[7] + '\t' + col[8] + '\t' + col[1] + '_2' + '\t' + '1' + '\t' + col[4]

        ##write out the bedline to a file
        with open (ltr_dir + '/opt_temp.bed','w') as opt_bed:
            opt_bed.write(bedline)

        ##use the bedtools to extract the sequence
        cmd = bedtools_exe + ' getfasta -fi ' + input_genome_fas_file + ' -bed ' + ltr_dir + '/opt_temp.bed' + \
              ' -s -name -fo ' + ltr_dir + '/opt_temp_te_ltr.fasta'
        subprocess.call(cmd, shell=True)

        ##use the muscle to align the sequence
        cmd = muscle_exe + ' -in ' + ltr_dir + '/opt_temp_te_ltr.fasta' + ' -out ' + ltr_dir + '/opt_temp_align.fasta'
        subprocess.call(cmd, shell=True)

        ##predict the time
        ##store the id and seq information in a dic
        id_seq_dic = {}
        length = 0
        for seq_record in SeqIO.parse(ltr_dir + '/opt_temp_align.fasta', 'fasta'):
            id_seq_dic[seq_record.id] = seq_record.seq
            length = len(seq_record.seq)

        tag_all = 0
        tag_div = 0

        for i in range(length):
            #print('the i is ' + str(i))
            letter1 = id_seq_dic[list(id_seq_dic.keys())[0]][i]
            if letter1 == '-' or letter1 == 'N':
                continue
            letter2 = id_seq_dic[list(id_seq_dic.keys())[1]][i]
            if letter2 == '-' or letter2 == 'N':
                continue
            tag_all += 1
            if letter1 != letter2:
                tag_div += 1

        print('the tag_all is ' + str(tag_all))
        print('the ta_div is ' + str(tag_div))

        ##calculate the time for each ltr name
        div_per = int(tag_div)/int(tag_all)
        print('the div_per is ' + str(div_per))
        if div_per < 0.75:
            #continue
        #else:
            time1 = -0.75*math.log(1 - (4 * div_per / 3))
            time2 = time1 * 100 / 2.6

            target_line = eachltr + '\t' + str(time2)
            print('the target_line is ' + target_line)
            compete_time_tb_dic[target_line] = 1
            #compete_time_tb_dic[eachltr] = str(compete_tb_dic[eachltr]) + '\t' + str(time2)

    return (compete_time_tb_dic)


















