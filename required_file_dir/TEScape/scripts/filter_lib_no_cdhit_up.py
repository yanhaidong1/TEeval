#!/usr/bin/env python

##update: 1.25 with out cd-hit
##Because all the papers do not use cd-hit to filter the libraries


from Bio import SeqIO
import re
import os
import subprocess
import glob


#######################################################################
##Step 1: divide into different superfamily which will be ran by cd-hit
def check_double_lib_source_cd_hit (file,cdhit_exe):
    #initiate a dic store the number of seq
    sr_dic = {}
    for seq_record in SeqIO.parse(file,'fasta'):
        sr_dic[seq_record.id] = 1
    #calculate the key number in the dictionary
    sr_key_nm = len(sr_dic.keys())
    if sr_key_nm >= 2:
        ##create a file dictionary
        dir_nm = file + '_cdhit' + '/'
        os.makedirs(dir_nm)
        ##run cd-hit for this file
        ##create a file name as a optname
        mt = re.match('(.+)\/(.+)\.lib',file)
        new_fnm = mt.group(2)
        iptdir = mt.group(1)
        cmd = cdhit_exe + " -i " + file + " -o "+ dir_nm + "/" + new_fnm + \
              ".lib -c 0.8 -G 1 -g 1 -prog blastn -circle 0 -exec local"
        subprocess.call(cmd, shell=True)

        new_fnm_lib = dir_nm + new_fnm + '.lib'
        print ('the new_fnm_lib is ' + new_fnm_lib)
        file_list = glob.glob(dir_nm + "/*")
        print ('the file_list is')
        print (file_list)
        print ('the dir_nm is ' + dir_nm)
        if new_fnm_lib in file_list:
            ##copy the output to the dictionary to store all libs
            cmd = "cp " + dir_nm + "/" + new_fnm + ".lib " + iptdir + "/all_filtered_file_dir/"
            subprocess.call(cmd, shell=True)
            print(cmd)
        else: ##if there is no output from the cd-hit, use the lib information outside the dir
            cmd = "cp " + file + " " + iptdir + "/all_filtered_file_dir/"
            subprocess.call(cmd, shell=True)
            print(cmd)

    if sr_key_nm == 1:
        mt = re.match('(.+)\/(.+)\.lib', file)
        iptdir = mt.group(1)
        cmd = "cp " + file + " " + iptdir + "/all_filtered_file_dir/"
        subprocess.call(cmd, shell=True)


def run_cd_hit_sep (ipt_dir,lib_construction_dir):
    ##create a dictionary to store the filtered file after cd-hit
    dir_nm_allst = ipt_dir + '/all_filtered_file_dir/'
    if not os.path.exists(dir_nm_allst):  ##if dir not exit, it will create the dir
        os.makedirs(dir_nm_allst)
    file_list = glob.glob(ipt_dir + '/*.lib')

    #######################
    ##attention: comment these commands
    for eachfl in file_list:

        #check_double_lib_source_cd_hit(eachfl,cdhit_exe)
        cmd = 'cp ' + eachfl + " " + ipt_dir + "/all_filtered_file_dir/"
        subprocess.call(cmd, shell=True)

    ##combination of all the libraries in the dir storing the libs
    lib_list = ''
    lib_file_list = glob.glob(ipt_dir + '/all_filtered_file_dir/*.lib')
    for eachlib in lib_file_list:
        lib_list = eachlib + ' ' + lib_list
    print(str(lib_list))
    ##run the cat to get the final libaray
    cmd = "cat " + lib_list + " > " + lib_construction_dir + "/final.lib"
    print ('the cmd is ' + cmd)
    subprocess.call(cmd, shell=True)