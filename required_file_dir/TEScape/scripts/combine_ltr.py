#!/usr/bin/env python

##update 1.11 generate the ltr information in three situations:
##combination between ltrfinder and ltrharvest
##ltrfinder
##ltrharvest

##This script is to generate the ltr filtered table

import re
from Bio import SeqIO
import pandas as pd


##Step 1: generate the final LTRs file by combining ltrfinder and ltrharvest

######################################################
##Step 1: combination between ltrfinder and ltrharvest
######################################################

###################################################
##1.1: generate the te line and connect all the tes
##inital a function to store te information from the ltrfinder and ltrdigest and ltrharvest
def store_te (ltrfinder_file,ltrharvest_file,ltrdigest_file,genome_file):
    ##ltrharvest could provide information about

    te_combined_dic = {}

    ##for the ltrfinder file
    ltrfinder_count = 0
    with open (ltrfinder_file,'r') as ltrfinder:
        for eachline in ltrfinder:
            eachline = eachline.strip()
            if re.match('^\[.+',eachline):
                ltrfinder_count = ltrfinder_count + 1
                ltrnm = 'ltrfdr_' + str(ltrfinder_count)
                col = eachline.split('\t')
                ##get the chromosome information
                chr = col[1]
                ##get the te begin and end information
                location = col[2]
                mt = re.match('(.+)\-(.+)',location)
                te_begin = mt.group(1)
                te_end = mt.group(2)
                ##get the direct information
                direct = col[12]
                ##get the similarity information
                #sim = col[15]
                ##get the ltr information
                ltr_d_len = col[3]
                mt_ltr = re.match('(.+)\,(.+)',ltr_d_len)
                lf_ltr_len = mt_ltr.group(1)
                rt_ltr_len = mt_ltr.group(2)

                lf_ltr_bg = te_begin
                lf_ltr_ed = int(te_begin) + int(lf_ltr_len) - 1

                rt_ltr_bg = int(te_end) - int(rt_ltr_len) + 1
                rt_ltr_ed = te_end

                ##generate the ltrfinder line
                te_line = chr + '\t' + ltrnm + '\t' + te_begin + '\t' + te_end + '\t' + direct + '\t' + lf_ltr_bg + '\t' + str(lf_ltr_ed) + '\t' + str(rt_ltr_bg) + '\t' + rt_ltr_ed
                te_combined_dic[te_line] = 1



    ##get the chr name in the genome file
    chr_nm_nb_dic = {}
    chr_count = 0
    for seq_record in SeqIO.parse(genome_file, "fasta"):
        chr_nm_nb_dic[str(chr_count)] = seq_record.id
        chr_count = chr_count + 1


    ##for the ltrdigest store the strand direction information, and use the ID name as key and chr, up and down and reference key
    #ltrdigest_dic = {}
    ltrdigest_chr_bg_dic = {}
    with open (ltrdigest_file,'r') as ltrdigest:
        for eachline in ltrdigest:
            eachline = eachline.strip()
            if not re.match('^#.+',eachline):
                col = eachline.split()
                if re.match('ID=LTR_retrotransposon.+',col[8]):
                    mt = re.match('ID=(.+?)\;',col[8])
                    id = mt.group(1)
                    mt2 = re.match('seq(.+)',col[0])
                    chr_num = str(mt2.group(1))
                    #ltrdigest_dic[id] = {'chr':chr_num,'begin':col[3],'end':col[4],'direct':col[6]}
                    ##create a chromosome and begin string in order to search in the ltrharvest file
                    chr_bg_line = chr_num + '\t' + col[3]
                    ltrdigest_chr_bg_dic[chr_bg_line] = {'tenm':id,'chr':chr_num,'begin':col[3],'end':col[4],'direct':col[6]}

    ##for the ltrharvest file
    ltrharvest_count = 0
    with open (ltrharvest_file,'r') as ltrharvest:
        for eachline in ltrharvest:
            eachline = eachline.strip()
            if not re.match('^#.+',eachline):
                ##get the name of ltr
                ltrharvest_count = ltrharvest_count + 1

                #ltrnm = 'repeat_region' + str(ltrharvest_count)  ##the ltrnm should be changed to the name same with ltr_digest

                eachline = eachline.strip()
                #print(eachline)

                col = eachline.split()
                ##get the chromosome name
                chr_nb_har = col[10]
                ##get the chr and begin line
                chr_bg_line_har = chr_nb_har + '\t' + col[0]

                if chr_bg_line_har in ltrdigest_chr_bg_dic.keys():

                    if chr_nb_har in chr_nm_nb_dic.keys():
                        chr_name_har = chr_nm_nb_dic[str(chr_nb_har)]
                        ltr_name = ltrdigest_chr_bg_dic[chr_bg_line_har]['tenm']
                        te_begin = col[0]
                        te_end = col[1]
                        lf_ltr_bg = col[3]
                        lf_ltr_ed = col[4]
                        rt_ltr_bg = col[6]
                        rt_ltr_ed = col[7]
                        ltr_direct = ltrdigest_chr_bg_dic[chr_bg_line_har]['direct']

                        if ltrdigest_chr_bg_dic[chr_bg_line_har]['direct'] != '?':

                            ##generate the te line
                            te_line = chr_name_har + '\t' + ltr_name + '\t' + te_begin + '\t' + te_end + '\t' + ltr_direct + '\t' + lf_ltr_bg + '\t' + lf_ltr_ed + '\t' + rt_ltr_bg + '\t' + rt_ltr_ed

                            te_combined_dic[te_line] = 1



    return (te_combined_dic)


##inital a function to order the te table
def order_location (te_ltr_file):

    ltr_fl = pd.read_table(te_ltr_file, header=None)  ##delete delimiter=r"\s+"
    ##sort chr and start region
    ltr_fl.columns = ['Chr','TE','Begin','End','Direct','Lltr_bg','Lltr_ed','Rltr_bg','Rltr_ed']
    ltr_sort_fl = ltr_fl.sort_values(by=['Chr','Begin'])
    return(ltr_sort_fl)


##initial a function to store all the te information in a list containing multiple dictionary, each dictionary contain each te group.
##the te group could be one te or connect te
##for the connected tes, we will regard them as
##connect each te
def te_connect (te_ltr_tb_file):

    ##store the chromosome information
    chr_dic = {}
    with open (te_ltr_tb_file,'r') as ipt_tb:
        for eachline in ipt_tb:
            if not re.match('.*Chr.+TE.+',eachline):
                col = eachline.split('\t')
                #print(col[1])
                chr_dic[col[1]] = 1

    ##calculate te number in each chromosome
    chr_te_num_dic = {}
    with open (te_ltr_tb_file,'r') as ipt_tb:
        for eachline in ipt_tb:
            if not re.match('.*Chr.+TE.+',eachline):
                col = eachline.split('\t')
                chr_nm = col[1]
                if chr_nm in chr_te_num_dic.keys():
                    chr_te_num_dic[chr_nm] = chr_te_num_dic[chr_nm] + 1
                else:
                    chr_te_num_dic[chr_nm] = 1

    ##initial a dictionary to store all the information of te, and use the chromosome name as the key
    te_all_dic = {}

    ##initial a dictionary to store the id number
    chr_id_dic = {}

    for eachchr in chr_dic:
        ##if the chr change to the next one, it should to initial an empty te_dic and an empty dic_list
        te_dic = {}
        dic_list = []

        with open (te_ltr_tb_file,'r') as ipt_tb:
            for eachline in ipt_tb:
                if not re.match('.*Chr.+TE.+Begin.+End.+',eachline):
                    eachline = eachline.strip('\n')
                    ##each chromosome should have its own list of te information
                    ##each chromosome has its own te id
                        ##split each line
                    col = eachline.split('\t')
                    chr_nm = col[1]

                    if  chr_nm == eachchr:
                        if eachchr in chr_id_dic.keys():
                            chr_id_dic[eachchr] = chr_id_dic[eachchr] + 1
                        else:
                            chr_id_dic[eachchr] = 1

                        id = str(chr_id_dic[eachchr])
                        te_nm = col[2]
                        te_begin = col[3]
                        te_end = col[4]
                        direct = col[5]

                        if id == str(1):
                            #print('the first id is ' + str(id))
                            te_dic[id] = {'chr': chr_nm, 'name': te_nm,'begin': te_begin, 'end': te_end,
                                          'direct': direct,'lltr_bg':col[6],'lltr_ed':col[7],'rltr_bg':col[8],
                                          'rltr_ed':col[9]}  ##no similarity
                            #print('the end is ' + te_dic[id]['end'])

                        if id > str(1):
                            #print('te_begin is ' + te_begin)
                            if te_begin <= te_dic[str(int(id)-1)]['end']:
                                te_dic[id] = {'chr': chr_nm, 'name': te_nm, 'begin': te_begin, 'end': te_end,
                                              'direct': direct, 'lltr_bg': col[6], 'lltr_ed': col[7], 'rltr_bg': col[8],
                                              'rltr_ed': col[9]} ##no similarity
                            else:
                                dic_list.append(te_dic)
                                te_dic = {}
                                te_dic[id] = {'chr': chr_nm, 'name': te_nm, 'begin': te_begin, 'end': te_end,
                                              'direct': direct, 'lltr_bg': col[6], 'lltr_ed': col[7], 'rltr_bg': col[8],
                                              'rltr_ed': col[9]} ##no similarity

                            if id == str(chr_te_num_dic[eachchr]):
                                dic_list.append(te_dic)

        te_all_dic[eachchr] = dic_list  ##this dictionary contain the key which is the chr name and value which is the dic list

    return(te_all_dic)


################
##1.2: filter te
##initial a function to store the connection information
def transfer (id_order,te_comp_dic,te_dic):

    te_comp_dic[id_order] = {'chr': te_dic[id_order]['chr'], 'name': te_dic[id_order]['name'], 'begin': te_dic[id_order]['begin'],'end': te_dic[id_order]['end'],'direct': te_dic[id_order]['direct'],'lltr_bg': te_dic[id_order]['lltr_bg'],'lltr_ed': te_dic[id_order]['lltr_ed'],'rltr_bg': te_dic[id_order]['rltr_bg'],'rltr_ed': te_dic[id_order]['rltr_ed']}

##initial a function to filter te covered by the same te
##the argument is the te which store all the information
def filter_te (te_all_dic):

    ##initial a dic to store the filtered te
    te_line_dic = {}

    ##te_all_dic's key is the chromosome dictionary of each store a list
    ##each list sotre multiple dictionary. Some dic contains one key and some dic contains multiple keys
    for eachchr_dic in te_all_dic:
        dic_list = te_all_dic[eachchr_dic]


        for tar_dic in dic_list: ##each dic store one key or multiple keys
            #tar_dic = te_all_dic[eachchr_dic][eachdic]  ##name a dic to a other name to decrease the length of dic


            if len(tar_dic) == 1:
                for eachid in tar_dic:  ##call the each key in this tar_dic
                    tar_te_line = tar_dic[eachid]['chr'] + '\t' + tar_dic[eachid]['name'] + '\t' + tar_dic[eachid]['begin'] + '\t' + \
                                  tar_dic[eachid]['end'] + '\t' + tar_dic[eachid]['direct'] + '\t' + tar_dic[eachid]['lltr_bg'] + '\t' + \
                                  tar_dic[eachid]['lltr_ed'] + '\t' + tar_dic[eachid]['rltr_bg'] + '\t' + tar_dic[eachid]['rltr_ed']  ##no sim

                    te_line_dic[tar_te_line] = 1


            ##if the length of tar_dic is not 1
            if len(tar_dic) > 1:
                ##create a another te_dic to store the comparision results from this te_dic
                te_comp_dic = {}

                for eachid in tar_dic:

                    ######################################
                    ##Analyze the first two keys in te_dic
                    if eachid != list(tar_dic)[0]:
                        ##extract the length of eachid
                        ##previous length
                        #pre_len = int(tar_dic[str(int(eachid)-1)]['end']) - int(tar_dic[str(int(eachid)-1)]['begin'])
                        ##current length
                        cur_len = int(tar_dic[eachid]['end']) - int(tar_dic[eachid]['begin'])

                        if eachid == list(tar_dic)[1]:

                            ##previous length
                            pre_len = int(tar_dic[str(int(eachid) - 1)]['end']) - int(tar_dic[str(int(eachid) - 1)]['begin'])

                            #cur_len = int(tar_dic[eachid]['end']) - int(tar_dic[eachid]['begin'])

                            ##if the first is longer than the second one
                            if pre_len >= cur_len:
                                ##import the first id to the te_comp_dic
                                transfer(str(int(eachid)-1),te_comp_dic,tar_dic)
                            ##if the second is longer than the first one
                            else:
                                ##import the second id to the te_comp_dic
                                transfer(eachid, te_comp_dic, tar_dic)

                        #################################
                        ##Analyze the rest keys in te_dic
                        ##the current of dic should compare with the te in the compare dictionary
                        else:
                            com_len = int(te_comp_dic[list(te_comp_dic)[-1]]['end']) -  int(te_comp_dic[list(te_comp_dic)[-1]]['begin'])
                            ##if there is cover between the current one and the compare dic one
                            if tar_dic[eachid]['begin'] < te_comp_dic[list(te_comp_dic)[-1]]['end']:
                                ##if current one is longer than the compared one
                                if cur_len >= com_len:
                                    ##remove the compared one and store the current one
                                    te_comp_dic.pop(list(te_comp_dic)[-1])
                                    transfer(eachid, te_comp_dic, tar_dic)
                                else:
                                    continue ##continue the analysis to ignore the current id

                            ##if there is no cover between the current one and the compare dic one
                            else:
                                transfer(eachid, te_comp_dic, tar_dic)


                ##store the te from compared dic to the te_line_dic
                for eachid in te_comp_dic:
                    #print('the target id count is ' + str(count))
                    tar_te_line = te_comp_dic[eachid]['chr'] + '\t' + te_comp_dic[eachid]['name'] + '\t' + te_comp_dic[eachid]['begin'] + '\t' + \
                                  te_comp_dic[eachid]['end'] + '\t' + te_comp_dic[eachid]['direct'] + '\t' + te_comp_dic[eachid]['lltr_bg'] + '\t' + \
                                  te_comp_dic[eachid]['lltr_ed'] + '\t' + te_comp_dic[eachid]['rltr_bg'] + '\t' + te_comp_dic[eachid]['rltr_ed'] ##no sim
                    te_line_dic[tar_te_line] = 1

    return(te_line_dic)

##initiate a function to generate bed line
def generate_bed (te_line_dic):
    bed_line_dic = {}
    for eachline in te_line_dic:
        col = eachline.split()
        direct = ''
        if col[4] == '-':
            direct = 'C'
        else:
            direct = col[4]
        bed_line = col[0] + '\t' + col[2] + '\t' + col[3] + '\t' + col[1] + '\t1\t' + direct
        bed_line_dic[bed_line] = 1
    return (bed_line_dic)



#################################
##Step 2: only generate ltrfinder
#################################
def store_ltrfinder_te (ltrfinder_file):

    ##store the ltrfinder_dic
    ltrfinder_line_dic = {}

    ##for the ltrfinder file
    ltrfinder_count = 0
    with open (ltrfinder_file,'r') as ltrfinder:
        for eachline in ltrfinder:
            eachline = eachline.strip()
            if re.match('^\[.+',eachline):
                ltrfinder_count = ltrfinder_count + 1
                ltrnm = 'ltrfdr_' + str(ltrfinder_count)
                col = eachline.split('\t')
                ##get the chromosome information
                chr = col[1]
                ##get the te begin and end information
                location = col[2]
                mt = re.match('(.+)\-(.+)',location)
                te_begin = mt.group(1)
                te_end = mt.group(2)
                ##get the direct information
                direct = col[12]
                ##get the similarity information
                #sim = col[15]
                ##get the ltr information
                ltr_d_len = col[3]
                mt_ltr = re.match('(.+)\,(.+)',ltr_d_len)
                lf_ltr_len = mt_ltr.group(1)
                rt_ltr_len = mt_ltr.group(2)

                lf_ltr_bg = te_begin
                lf_ltr_ed = int(te_begin) + int(lf_ltr_len) - 1

                rt_ltr_bg = int(te_end) - int(rt_ltr_len) + 1
                rt_ltr_ed = te_end

                ##generate the ltrfinder line
                te_line = chr + '\t' + ltrnm + '\t' + te_begin + '\t' + te_end + '\t' + direct + '\t' + lf_ltr_bg + '\t' + str(lf_ltr_ed) + '\t' + str(rt_ltr_bg) + '\t' + rt_ltr_ed
                ltrfinder_line_dic[te_line] = 1

    return (ltrfinder_line_dic)




##################################
##Step 3: only generate ltrharvest
##################################
def store_ltrharvest_te (ltrharvest_file,ltrdigest_file,genome_file):

    ltrharvest_line_dic = {}

    ##get the chr name in the genome file
    chr_nm_nb_dic = {}
    chr_count = 0
    for seq_record in SeqIO.parse(genome_file, "fasta"):
        chr_nm_nb_dic[str(chr_count)] = seq_record.id
        chr_count = chr_count + 1

    ##for the ltrdigest store the strand direction information, and use the ID name as key and chr, up and down and reference key
    # ltrdigest_dic = {}
    ltrdigest_chr_bg_dic = {}
    with open(ltrdigest_file, 'r') as ltrdigest:
        for eachline in ltrdigest:
            eachline = eachline.strip()
            if not re.match('^#.+', eachline):
                col = eachline.split()
                if re.match('ID=LTR_retrotransposon.+', col[8]):
                    mt = re.match('ID=(.+?)\;', col[8])
                    id = mt.group(1)
                    mt2 = re.match('seq(.+)', col[0])
                    chr_num = str(mt2.group(1))
                    # ltrdigest_dic[id] = {'chr':chr_num,'begin':col[3],'end':col[4],'direct':col[6]}
                    ##create a chromosome and begin string in order to search in the ltrharvest file
                    chr_bg_line = chr_num + '\t' + col[3]
                    ltrdigest_chr_bg_dic[chr_bg_line] = {'tenm': id, 'chr': chr_num, 'begin': col[3], 'end': col[4],
                                                         'direct': col[6]}

    ##for the ltrharvest file
    ltrharvest_count = 0
    with open(ltrharvest_file, 'r') as ltrharvest:
        for eachline in ltrharvest:
            eachline = eachline.strip()
            if not re.match('^#.+', eachline):
                ##get the name of ltr
                ltrharvest_count = ltrharvest_count + 1

                # ltrnm = 'repeat_region' + str(ltrharvest_count)  ##the ltrnm should be changed to the name same with ltr_digest

                eachline = eachline.strip()
                # print(eachline)

                col = eachline.split()
                ##get the chromosome name
                chr_nb_har = col[10]
                ##get the chr and begin line
                chr_bg_line_har = chr_nb_har + '\t' + col[0]

                if chr_bg_line_har in ltrdigest_chr_bg_dic.keys():

                    if chr_nb_har in chr_nm_nb_dic.keys():
                        chr_name_har = chr_nm_nb_dic[str(chr_nb_har)]
                        ltr_name = ltrdigest_chr_bg_dic[chr_bg_line_har]['tenm']
                        te_begin = col[0]
                        te_end = col[1]
                        lf_ltr_bg = col[3]
                        lf_ltr_ed = col[4]
                        rt_ltr_bg = col[6]
                        rt_ltr_ed = col[7]
                        ltr_direct = ltrdigest_chr_bg_dic[chr_bg_line_har]['direct']

                        if ltrdigest_chr_bg_dic[chr_bg_line_har]['direct'] != '?':
                            ##generate the te line
                            te_line = chr_name_har + '\t' + ltr_name + '\t' + te_begin + '\t' + te_end + '\t' + ltr_direct + '\t' + lf_ltr_bg + '\t' + lf_ltr_ed + '\t' + rt_ltr_bg + '\t' + rt_ltr_ed

                            ltrharvest_line_dic[te_line] = 1

    return(ltrharvest_line_dic)