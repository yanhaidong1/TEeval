#!/usr/bin/env python3.7


##update 5.13 add the specific evaluation script on the landscape evaluation

import sys
import re
import glob
import os
import subprocess

##input files
input_te_scape_opt = sys.argv[1]
input_raw_data = sys.argv[2]
input_accu_dir = sys.argv[3] ##this input file means the analyzed opt will be stored in this accuracy dir

input_specific_script_path = sys.argv[4]


def batch_evaluation (input_te_scape_opt,input_raw_data,input_accu_dir,input_specific_script_path):

    te_loc_list = glob.glob(input_te_scape_opt + '/*')

    for eachte_dir in te_loc_list:

        ##this target opt fl will be used as ipt for accuracy detection
        target_opt_fl = ''

        file_list = glob.glob(eachte_dir + '/*')

        for eachfile in file_list:

            if 'output_dir' in eachfile:

                target_file_list = glob.glob(eachfile + '/*')

                for eachfile_2 in target_file_list:

                    if 'opt_final_tb.txt' in eachfile_2:

                        target_opt_fl = eachfile_2
        ##get the name of species
        mt = re.match('.+/(\w+)$',eachte_dir)
        sp_nm = mt.group(1)


        ##make empty dir for the accuracy
        sp_acu_dir = input_accu_dir + '/' + sp_nm
        if not os.path.exists(sp_acu_dir):  ##if dir not exit, it will create the dir
            os.makedirs(sp_acu_dir)

        ##conduct the analysis
        cmd = 'python3.7 ' + input_specific_script_path + ' ' + input_raw_data + '/opt_insert_te_location.txt ' + target_opt_fl + ' ' + sp_acu_dir + \
            ' ' + input_raw_data + '/opt_final.fasta'
        print (cmd)
        subprocess.call(cmd, shell=True)



batch_evaluation (input_te_scape_opt,input_raw_data,input_accu_dir,input_specific_script_path)

