#!/usr/bin/env python3.7

##update: 5.14 ##add the store dir for the output

import sys
import re
import glob
import os
import subprocess

input_accu_dir = sys.argv[1]
type_list = sys.argv[2]

input_store_dir = sys.argv[3] ##add the store dir for the output

##define a function to store all the te information
def store_all_accu_infor (input_accu_dir):
    ##initiate a dic to store all the information
    ##the key is each sample species, the value is the dic that contains type(key) and number (value)
    all_infor_dic = {}

    accu_file_dir_list = glob.glob(input_accu_dir + '/*')

    for eachfl_dir in accu_file_dir_list:

        type_vl_dic = {}

        ##get the name of fl_dir
        mt = re.match('.+/(\w+)$', eachfl_dir)
        sp_nm = mt.group(1)

        ##glob each file in the eachfl dir
        accu_fl_list = glob.glob(eachfl_dir + '/*')

        for eachfl in accu_fl_list:

            if 'opt_evaluate_tb.txt' in eachfl or 'opt_tmpt_evaluate_tb.txt' in eachfl:

                ##store information for this file
                with open (eachfl,'r') as ipt:

                    for eachline in ipt:
                        eachline = eachline.strip('\n')
                        col = eachline.strip().split('\t')
                        type = col[0]
                        vl = col[1]
                        type_vl_dic[type] = vl

                        ##store them in the all dic
                        all_infor_dic[sp_nm] = type_vl_dic

    return (all_infor_dic)


##define a function to write out results
def write_tb (all_infor_dic,type_list):

    ##initiate a dic to store the final_line
    final_line_dic = {}


    with open (type_list,'r') as ipt:
        for eachtp in ipt:

            eachtp =eachtp.strip('\n')

            vl_string = ''

            for eachnm in all_infor_dic:

                if eachtp not in all_infor_dic[eachnm]:
                    t_value = 'na'
                else:
                    t_value = all_infor_dic[eachnm][eachtp]

                vl_string =  vl_string + '\t' + t_value

            final_line = eachtp + vl_string

            final_line_dic[final_line] = 1

    return (final_line_dic)


##write out results
all_infor_dic = store_all_accu_infor (input_accu_dir)
final_line_dic = write_tb (all_infor_dic,type_list)

name_string = ''
for eachnm in all_infor_dic:
    name_string = name_string + '\t' + eachnm

with open (input_store_dir + '/opt_final_combine_all_genus.txt','w+') as opt:
    opt.write('Type' + name_string + '\n')
    for eachline in final_line_dic:
        opt.write(eachline + '\n')




