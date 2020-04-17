#!/usr/bin/env python3.7

import sys
import re
import glob
import os
import subprocess



##input file
input_raw_landscape_dir = sys.argv[1]
input_opt_landscape_dir = sys.argv[2]
input_ipt_reduced_lib_dir = sys.argv[3]
input_genome_fas = sys.argv[4]

def batch_landscape (input_raw_landscape_dir,input_opt_landscape_dir,input_ipt_reduced_lib_dir):

    ##conduct analysis for each lib
    lib_fl_list = glob.glob(input_ipt_reduced_lib_dir + '/*')

    ##extract the name of the lib file
    for eachfl in lib_fl_list:
        mt = re.match('.+/(.+)\.lib',eachfl)
        lib_nm = mt.group(1)

        ##make the dir in the landscape dir
        lib_dir = input_opt_landscape_dir + '/' + lib_nm
        if not os.path.exists(lib_dir):
            os.makedirs(lib_dir)

        ##make the working_dir and output_dir in the lib_dir
        working_dir = lib_dir + '/working_dir'
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)

        output_dir = lib_dir + '/output_dir'
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        ##run the landscape
        cmd = 'python3.7 ' + input_raw_landscape_dir + '/TEScape/TEScape_up_nocdhit.py -d ' + working_dir + ' -o ' + output_dir \
              + ' -fas ' + input_genome_fas + ' -i ' + eachfl + ' -s ' + input_raw_landscape_dir + '/TEScape/supfile_dir/' \
              + ' --cdhit ~/software/cdhit-master/psi-cd-hit/psi-cd-hit.pl' \
              + ' --repeatmasker /data/software/RepeatMasker/RepeatMasker' \
              + ' --bedtools /usr/bin/bedtools' \
              + ' --hmmscan /usr/bin/hmmscan'

        print(cmd)
        subprocess.call(cmd,shell=True)


batch_landscape (input_raw_landscape_dir,input_opt_landscape_dir,input_ipt_reduced_lib_dir)