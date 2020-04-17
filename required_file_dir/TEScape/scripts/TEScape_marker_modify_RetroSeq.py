#!/usr/bin/env python3.7

##this script is to modify the locationlist in the reference of the micclintock.
##delete the file without bed or without fasta
##import modules
import re
import glob
import sys
import subprocess

#input_location_dir = sys.argv[1]

#input_location_list = sys.argv[2]


def modify_location (input_location_dir,input_location_lis):

    ##initiate a dic to store all the te
    te_store_dic = {}

    ##initiate a dic to store the te with more thant 1 in the te_store_dic
    te_final_store_dic = {}

    ##initiate a dic to store the final line of location list
    final_line_loc_dic = {}

    all_file_list = glob.glob(input_location_dir + '/*')

    for eachfile in all_file_list:

        ##extract the name of TEs
        mt = re.match("(.+)/(.+)\..+",eachfile)
        te_nm = mt.group(2)
        path = mt.group(1)

        if te_nm in te_store_dic:
            te_store_dic[te_nm] += 1

            ##if the te_nm has been stored in the te_store_dic
            ##cp these two files to the output_location_dir
            #cmd = 'cp ' + path + '/' + te_nm + '.* ' + output_location_dir
            #subprocess.call(cmd,shell=True)

            te_final_store_dic[te_nm] = 1

        else:
            te_store_dic[te_nm] = 1

    ##remove the te without bed files in the location dir
    with open (input_location_lis,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            ##extract the name of TEs
            mt = re.match(".+/(.+)\..+", col[1])
            te_fl_nm = mt.group(1)

            ##remove the wrong te
            if te_fl_nm in te_final_store_dic.keys():

                final_line_loc_dic[eachline] = 1

    return (final_line_loc_dic)



#final_line_loc_dic = modify_location (input_location_dir,input_location_list)

##write out results
#with open ('opt_new_locationlist','w+') as opt:
#    for eachline in final_line_loc_dic:
#        opt.write(eachline + '\n')












