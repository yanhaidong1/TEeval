#!/usr/bin/env python

##update1.10 set the same TE connected each other that will be considered as fusion

##this script has been updated, and is aming to modify original detected results
##this updating will use a new output file, which has been added the results from ltrfinder and ltrharvest from the 05_2_modify_gff_rpmsk
##this updating has also changed the comparison by using the int function

##this script is to modify TEs (connect close te)
##this script is to create the match rate table
##this script is to extract the complete TEs sequences ready to the hmm

##opt files
##opt1: store all the information in the table including connection and match rate
##opt2: store all the information as bed files
##opt3: store all the domain te as bed files

##import functions
import sys
import re
import pandas as pd


#########################################
##step 1: store all information in a list
#########################################

###########################################
##initiate a function to order the te table
def order_location (te_gff_file):

    ltr_fl = pd.read_table(te_gff_file, header=None,low_memory=False)  ##delete delimiter=r"\s+"
    ##sort chr and start region
    ltr_fl.columns = ['Chr','Begin','End','Direct','Name','Lib_bg','Lib_end','Lib_left']
    ltr_sort_fl = ltr_fl.sort_values(by=['Chr','Begin'])
    return(ltr_sort_fl)

#####################################################
##initiate a function to add code for each te in the gff file
def code_te (te_sort_tb_file):
    final_dic = {}
    with open(te_sort_tb_file,'r') as ipt_sort_tb:
        id = 0
        for eachline in ipt_sort_tb:
            id += 1
            eachline = eachline.strip('\n')
            final_line = eachline + '\t' + str(id)
            final_dic[final_line] = 1
    return(final_dic)

#################################################
##initiate a function to store the te information
def store_te_infor (te_sort_id_tb_file):

    ##initial a list contain all the TE situations
    dic_list = []
    ##initial a dictionary to store the target lines
    dic_te = {}

    ##get the last id
    last_id = ''
    with open(te_sort_id_tb_file, 'r') as ipt_rmk_out:
        last_line = ipt_rmk_out.readlines()[-1]
        last_col = last_line.strip().split()
        last_id = last_col[8]

    ##store the te infor
    with open(te_sort_id_tb_file, 'r') as ipt_rmk_out:

        for line in ipt_rmk_out:
            col = line.strip().split()

            chr = col[0]
            te_nm = col[4]
            id = col[-1]
            bg = col[1]
            ed = col[2]
            dir = col[3]
            lib_bg = col[5]
            lib_ed = col[6]
            lib_left = col[7]

            if id == str(1):  ##if the id is 1, it will directly store in the dic_te
                dic_te[id] = {'chr': chr, 'begin': bg, 'end': ed, 'direct': dir, 'name': te_nm,
                              'lib_bg': lib_bg, 'lib_ed': lib_ed, 'lib_left': lib_left}

            else: ##if the id is over 1
                ##if the chr is the same as the previous one
                if dic_te[str(int(id)-1)]['chr'] == chr:

                    if dic_te[str(int(id)-1)]['name'] == te_nm: ##if the name is same
                        dic_te[id] = {'chr':chr,'begin':bg,'end':ed,'direct':dir,'name':te_nm,
                                      'lib_bg':lib_bg,'lib_ed':lib_ed,'lib_left':lib_left}

                    else:
                        if int(bg) <= int(dic_te[str(int(id)-1)]['end']):  ##if the begin of the current is smaller than the previous end, this means there is fusion for these two tes
                            dic_te[id] = {'chr':chr,'begin':bg,'end':ed,'direct':dir,'name':te_nm,
                                          'lib_bg':lib_bg,'lib_ed':lib_ed,'lib_left':lib_left}

                        else:
                            dic_list.append(dic_te)
                            dic_te = {}
                            dic_te[id] = {'chr':chr,'begin':bg,'end':ed,'direct':dir,'name':te_nm,'lib_bg':lib_bg,
                                          'lib_ed':lib_ed,'lib_left':lib_left}

                            if id == last_id:  ##should be out of the previous loop
                                dic_list.append(dic_te)

                else: ##if the chr is not the same as the previous one
                    ##store the dic_te which has been stored the te informations
                    dic_list.append(dic_te)
                    dic_te = {}
                    dic_te[id] = {'chr': chr, 'begin': bg, 'end': ed, 'direct': dir, 'name': te_nm, 'lib_bg': lib_bg,
                                  'lib_ed': lib_ed, 'lib_left': lib_left}

                    if id == last_id:  ##should be out of the previous loop
                        dic_list.append(dic_te)

    #print(dic_list)
    return (dic_list)


######################################################################
##step 2: define a function to store the information for the small dic
######################################################################
##define a function to store the connection information
def transfer (id_order,te_comp_dic,te_dic,type):
    te_comp_dic[id_order] = {'chr': te_dic[id_order]['chr'],'begin': te_dic[id_order]['begin'],
                             'end': te_dic[id_order]['end'],'direct': te_dic[id_order]['direct'],
                             'name': te_dic[id_order]['name'],'lib_bg': te_dic[id_order]['lib_bg'],
                             'lib_ed': te_dic[id_order]['lib_ed'],'lib_left': te_dic[id_order]['lib_left'],
                             'type': type}

    if id_order != list(te_dic)[0]:

        if type == 'Unclear':
            te_comp_dic[list(te_comp_dic)[-2]]['type'] = te_comp_dic[list(te_comp_dic)[-2]]['type'] + '/Unclear'

        if type == 'Fusion':
            te_comp_dic[list(te_comp_dic)[-2]]['type'] = te_comp_dic[list(te_comp_dic)[-2]]['type'] + '/Fusion'

##define a function to analyze each dic in the dic_list and store the information to the te_line_dic
def iden_te_type (te_dic,te_line_dic):


    ##the connect type rule is as follows:
    ##example: Single/Fusion: te has no relation to the previous one but has fusion relation with the next one

    ##if the key of the te_dic is equal to 1
    if len(te_dic.keys()) == 1:
        for key in te_dic:
            ##store all the information
            tar_te_line = te_dic[key]['chr'] + '\t' + te_dic[key]['begin'] + '\t' + te_dic[key]['end'] + '\t' + \
                          te_dic[key]['name'] + '\t' + te_dic[key]['direct'] + '\t' + te_dic[key]['lib_bg'] + '\t' + \
                          te_dic[key]['lib_ed'] + '\t' + te_dic[key]['lib_left'] + '\t' + 'Single'
            te_line_dic[tar_te_line] = 1

    ##if the key number of the te_dic is not equal to 1
    else:
        ##create a another te_dic to store the comparision results from this te_dic
        te_comp_dic = {}
        for eachid in te_dic:

            ##note: the location of te has been ordered
            ##store the first te to the comp dic
            ##any dic should allow the first te id as the single, because it has no relation to the previous one

            ##if the id is equal to the first id, store this id as single
            if eachid == list(te_dic)[0]:
                transfer(eachid, te_comp_dic, te_dic, 'Single')

            ######################
            ##Analyze the rest TEs
            ######################
            else:
                ##1: if previous one is larger than the current one (cover)
                if int(te_dic[eachid]['end']) <= int(te_comp_dic[list(te_comp_dic)[-1]]['end']):
                    ##reserve the previous one, no matter the name and the direct
                    continue
                ##2: if previous one is not larger than the current one (not cover)
                else:
                    ##update 1.10
                    ##if they have the same name
                    if te_dic[eachid]['name'] == te_comp_dic[list(te_comp_dic)[-1]]['name']:
                        ##if they have same direction
                        if te_dic[eachid]['direct'] == te_comp_dic[list(te_comp_dic)[-1]]['direct']:

                            ##if the direction is the +
                            if te_dic[eachid]['direct'] == '+':

                                ##if the tes are partially covered
                                if int(te_dic[eachid]['begin']) <= int(te_comp_dic[list(te_comp_dic)[-1]]['end']):

                                    ##if the lib location is covered
                                    if int(te_dic[eachid]['lib_ed']) >= int(te_comp_dic[list(te_comp_dic)[-1]]['lib_bg']) and \
                                            int(te_dic[eachid]['lib_bg']) < int(te_comp_dic[list(te_comp_dic)[-1]]['lib_ed']):
                                        transfer(eachid, te_comp_dic, te_dic, 'Fusion')
                                    ##if the lib location is not covered
                                    else:
                                        # print('the third one is the this situation')
                                        transfer(eachid, te_comp_dic, te_dic, 'Unclear')

                                ##if the tes are not covered
                                else:
                                    ##if the lib location is covered
                                    if int(te_dic[eachid]['lib_ed']) >= int(te_comp_dic[list(te_comp_dic)[-1]]['lib_bg']) and \
                                            int(te_dic[eachid]['lib_bg']) < int(te_comp_dic[list(te_comp_dic)[-1]]['lib_ed']):
                                        transfer(eachid, te_comp_dic, te_dic, 'Single')
                                    ##if the lib location is not covered
                                    else:
                                        # print('the third one is the this situation')
                                        transfer(eachid, te_comp_dic, te_dic, 'Unclear')

                            ##if the direction is the C
                            else:
                                ##if the tes are partially covered
                                if int(te_dic[eachid]['begin']) <= int(te_comp_dic[list(te_comp_dic)[-1]]['end']):
                                    ##if the lib location is covered
                                    if int(te_comp_dic[list(te_comp_dic)[-1]]['lib_ed']) >= int(te_dic[eachid]['lib_left']) and \
                                            int(te_comp_dic[list(te_comp_dic)[-1]]['lib_left']) <= int(te_dic[eachid]['lib_ed']):
                                        transfer(eachid, te_comp_dic, te_dic, 'Fusion')
                                    ##if the lib location is not covered
                                    else:
                                        transfer(eachid, te_comp_dic, te_dic, 'Unclear')
                                ##if the tes are not covered
                                else:
                                    ##if the lib location is covered
                                    if int(te_comp_dic[list(te_comp_dic)[-1]]['lib_ed']) >= int(te_dic[eachid]['lib_left']) and \
                                            int(te_comp_dic[list(te_comp_dic)[-1]]['lib_left']) <= int(te_dic[eachid]['lib_ed']):
                                        transfer(eachid, te_comp_dic, te_dic, 'Single')
                                    ##if the lib location is not covered
                                    else:
                                        transfer(eachid, te_comp_dic, te_dic, 'Unclear')

                        ##if they do not have the same direction update1.10
                        else:
                            ##if the tes are partially covered
                            if int(te_dic[eachid]['begin']) <= int(te_comp_dic[list(te_comp_dic)[-1]]['end']):
                                transfer(eachid, te_comp_dic, te_dic, 'Fusion')
                            ##if the tes are not covered
                            else:
                                ##the current will be stored marked as Single
                                transfer(eachid, te_comp_dic, te_dic, 'Single')

                    ##if they are not the same name not matter the direction
                    else:
                        ##if covered
                        # print('the situation is here')
                        if int(te_dic[eachid]['begin']) <= int(te_comp_dic[list(te_comp_dic)[-1]]['end']) and \
                                int(te_dic[eachid]['end']) >= int(te_comp_dic[list(te_comp_dic)[-1]]['end']):
                            transfer(eachid, te_comp_dic, te_dic, 'Fusion')
                        ##if not covered
                        else:
                            transfer(eachid, te_comp_dic, te_dic, 'Single')


        ##store the infomration from the te_comp_dic to the dic_line
        for key in te_comp_dic:
            ##generate the all te information line
            tar_te_line = te_comp_dic[key]['chr'] + '\t' + te_comp_dic[key]['begin'] + '\t' + te_comp_dic[key]['end'] + '\t' + \
                          te_comp_dic[key]['name'] + '\t' + te_comp_dic[key]['direct'] + '\t' + te_comp_dic[key]['lib_bg'] + '\t' + \
                          te_comp_dic[key]['lib_ed'] + '\t' + te_comp_dic[key]['lib_left'] + '\t' + te_comp_dic[key]['type']
            te_line_dic[tar_te_line] = 1
            #print(tar_te_line)


############################################################################
##step 3: conduct analysis for each dic in the list to generate the bed file
############################################################################

##############################################################
##initiate a function to conduct the iden_te_type for each dic
def iden_type (dic_list):
    te_line_dic = {}
    count = 1
    for eachdic in dic_list:
        count = count + 1
        print('the analyzed te_dic count is ' + str(count))
        iden_te_type(eachdic,te_line_dic)
    return (te_line_dic)


###################################################################
##initiate a function to generate bed file for the TEs with domains
def gen_te_dm_te_final_bed (te_line_dic):
    ##initiate a dic to store the domain te bed domain
    bed_line_dm_te_dic = {}
    line_all_te_dic = {}

    ##initiate a te_nm_dic to rename te name
    te_nm_dic = {}
    for eachline in te_line_dic:
        ##rename each te which let them to show differnet name for the same te
        col = eachline.strip().split()
        tenm = col[3]
        new_te_nm = ''
        if tenm in te_nm_dic.keys():
            te_nm_dic[tenm] = te_nm_dic[tenm] + 1
            new_te_nm = tenm + '_' + str(te_nm_dic[tenm])
        else:
            te_nm_dic[tenm] = 1
            new_te_nm = tenm + '_1'


        ##update 1.7.2019 to select the LTR, LINE, and DNA_nMITE should be kept.
        if re.match('.*LTR.*', new_te_nm) or re.match('.*LINE.*', new_te_nm) or \
                re.match('.*DNA_nMITE.*', new_te_nm):

            bed_line = col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + new_te_nm + '\t1\t' + col[4]
            bed_line_dm_te_dic[bed_line] = 1

        final_line = col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + new_te_nm + '\t' + col[4] + '\t' + col[5] + '\t' + col[6] + '\t' + col[7] + '\t' + col[8]
        line_all_te_dic[final_line] = 1

    return(bed_line_dm_te_dic,line_all_te_dic)

#################################################################
##initiate a function to get the match rate from the final te dic
def match_rate (line_all_te_dic):

    line_all_te_mtrt_dic = {}

    for eachline in line_all_te_dic:
        col = eachline.split()
        tenm = col[3]
        if not re.match('.*combined.*',tenm):  ##do not contain ltrfinder and ltrharvest results
            ##because the length of te is different from te length calculated by lib comparison, we use length of te of the location in the genome
            te_len = int(col[2]) - int(col[1])

            #print('the te line is ' + teline)
            if re.match('\+',col[4]):
                lib_te_len = int(col[6]) + int(col[7])
                rmk_te_len = int(col[6]) - int(col[5]) + 1

                mt_rate = "%.4f" % (rmk_te_len / lib_te_len)

                all_line_mt_rate = eachline + '\t' + str(te_len) + '\t' + str(lib_te_len) + '\t' + str(mt_rate)
                line_all_te_mtrt_dic[all_line_mt_rate] = 1

            if re.match('C',col[4]):
                lib_te_len = int(col[6]) + int(col[5])
                rmk_te_len = int(col[6]) - int(col[7]) + 1

                mt_rate = "%.4f" % (rmk_te_len / lib_te_len)
                all_line_mt_rate = eachline + '\t' + str(te_len) + '\t' + str(lib_te_len) + '\t' + str(mt_rate)
                line_all_te_mtrt_dic[all_line_mt_rate] = 1

        else: ##if the te is ltrfinder and ltrharvest result
            all_line_mt_rate = eachline + '\t' + 'na' + '\t' + 'na' + '\t' + 'na'
            line_all_te_mtrt_dic[all_line_mt_rate] = 1

    return (line_all_te_mtrt_dic)





