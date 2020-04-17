#!/usr/bin/env python3.7

import sys
import re
import glob
import os
import subprocess



##input file
input_raw_landscape_dir = sys.argv[1]
input_opt_landscape_dir = sys.argv[2]
input_genome_fas = sys.argv[3]
ltr_software_type = sys.argv[4]
file_name = sys.argv[5]

def batch_landscape (file_name,input_raw_landscape_dir,input_opt_landscape_dir, input_genome_fas,ltr_software_type):

    ##each_rep_fl contain each replicate information

    ##make the dir in the landscape dir
    lib_dir = input_opt_landscape_dir + '/' + file_name
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
    if ltr_software_type == 'double':

        cmd = 'python3.7 ' + input_raw_landscape_dir + '/TEScape/TEScape_LTR.py -d ' + working_dir + ' -o ' + output_dir \
              + ' -fas ' + input_genome_fas + ' -s ' + input_raw_landscape_dir + '/TEScape/supfile_dir/' \
              + ' --cdhit ~/software/cdhit-master/psi-cd-hit/psi-cd-hit.pl' \
              + ' --repeatmasker /data/software/RepeatMasker/RepeatMasker' \
              + ' --bedtools /usr/bin/bedtools' \
              + ' --hmmscan /usr/bin/hmmscan' \
              + ' --muscle /usr/bin/muscle' \
              + ' --ltrfdr /data/haidongy_projects/software/LTR_Finder/source/ltr_finder' \
              + ' --gt /data/haidongy_projects/software/bin/gt/bin/gt' \
              + ' -para ' + input_raw_landscape_dir + '/TEScape/test_input_dir/ltr_parameter.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)

    if ltr_software_type == 'ltr_finder':

        cmd ='python3.7 ' + input_raw_landscape_dir + '/TEScape/TEScape_LTR.py -d ' + working_dir + ' -o ' + output_dir \
              + ' -fas ' + input_genome_fas + ' -s ' + input_raw_landscape_dir + '/TEScape/supfile_dir/' \
              + ' --cdhit ~/software/cdhit-master/psi-cd-hit/psi-cd-hit.pl' \
              + ' --repeatmasker /data/software/RepeatMasker/RepeatMasker' \
              + ' --bedtools /usr/bin/bedtools' \
              + ' --hmmscan /usr/bin/hmmscan' \
              + ' --muscle /usr/bin/muscle' \
              + ' --ltrfdropt ' + input_opt_landscape_dir + '/ltr_finder_harvest/output_dir/ltr_finder.opt' \
              + ' -para ' + input_raw_landscape_dir + '/TEScape/test_input_dir/ltr_parameter.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)


    if ltr_software_type == 'ltr_harvest':

        cmd = 'python3.7 ' + input_raw_landscape_dir + '/TEScape/TEScape_LTR.py -d ' + working_dir + ' -o ' + output_dir \
              + ' -fas ' + input_genome_fas + ' -s ' + input_raw_landscape_dir + '/TEScape/supfile_dir/' \
              + ' --cdhit ~/software/cdhit-master/psi-cd-hit/psi-cd-hit.pl' \
              + ' --repeatmasker /data/software/RepeatMasker/RepeatMasker' \
              + ' --bedtools /usr/bin/bedtools' \
              + ' --hmmscan /usr/bin/hmmscan' \
              + ' --muscle /usr/bin/muscle' \
              + ' --ltrhvtopt ' + input_opt_landscape_dir + '/ltr_finder_harvest/output_dir/ltr_harvest.opt ' \
              + ' --ltrltrditopt ' + input_opt_landscape_dir + '/ltr_finder_harvest/output_dir/ltr_digest.opt ' \
              + ' -para ' + input_raw_landscape_dir + '/TEScape/test_input_dir/ltr_parameter.txt'

        print(cmd)
        subprocess.call(cmd,shell=True)



batch_landscape (file_name,input_raw_landscape_dir,input_opt_landscape_dir,input_genome_fas,ltr_software_type)