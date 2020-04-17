#!/usr/bin/env python3.7

##updation5.6 identify the name should be consistent with the DNA

##dna_mite evaluation will help to detect the accuracy of dna mite
##in this situation do not consider the name because mite_finder2 and mite_hunter do not give us much information for family of mite

##updation4.29 change the FP to the TP if there is a cover one this TE
##from the evolution, some TE with gap cannot be connected each other.

##updation4.25 new simualted sequences containing mutation insertion and deletion and fusion
##we will focus on these four types one by one
##and calculate SE PR SP AC


##this script only focus output of the repeatmodeler
##because repmd generate different fam name in the output
##import modules
import re
import glob
import os
import sys
from Bio import SeqIO
import subprocess
from statistics import median




##input the location of real TE
input_real_TE_loc_file = sys.argv[1]
##input the opt file from the TEScape
input_opt_te_file = sys.argv[2]
input_opt_dir = sys.argv[3]

##input the fasta file to calculate the sequence length that will be used in identificaiton of identify_TN_FN
input_fasta_file = sys.argv[4]

###############################################
##Step 1 store the real TE location information
###############################################
def store_real_TE (input_real_TE_loc_file):
    real_te_line_dic = {}
    with open (input_real_TE_loc_file,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            te_line = col[0] + '\t' + col[1] + '\t' + col[2] ##the col[3] indicates the sequence string
            real_te_line_dic[te_line] = 1
    return (real_te_line_dic)

##############################################
##Step 2 store the opt TE location information
##############################################
def store_final_te_opt (input_opt_te_file):
    final_te_opt_dic = {}
    with open (input_opt_te_file,'r') as ipt:
        count = 0
        for eachline in ipt:
            count += 1
            eachline = eachline.strip('\n')
            final_line = eachline + '\t' + str(count)
            final_te_opt_dic[final_line] = 1
    return (final_te_opt_dic)



#######################################################
##Step 3 store the opt and real TE location information
#######################################################
##define a function to store each TEs (Real TEs) and their related TEs (TEScape TEs)
def store_te_and_final_te (real_te_line_dic,final_te_opt_dic):

    ##import 1: real_te_line_dic (Real TEs)
    ##import 2: final_te_opt_dic (TEScape TEs)

    ##initiate a all dic to store information
    ##the key is the each TE name
    ##the value is the list containing the op TE dic
    all_dic = {}

    re_te_dic = {}

    #count = 0
    for eachline in real_te_line_dic:
        #count += 1
        #print ('the analzyed eachline num is ' + str(count))

        eachline = eachline.strip('\n')
        col = eachline.strip().split()

        #if col[2] == 'transposon_fragment':
        re_chr = 'final_insert_seq'  ##the final_insert_seq is the name
        re_start = col[1]
        re_end = col[2]
        re_dir = '+'
        #re_annot = col[8]
        #re_annot_col = re_annot.split(';')

        #re_nm_col = re_annot_col[0].split(':')
        #mt = re.match('(.+)=(.+)', re_nm_col[0])
        #re_id = mt.group(2)

        re_nm = col[0]

        re_te_dic[re_nm] = {'chr':re_chr,'bg':re_start,'ed':re_end,"dir":re_dir,'nm':re_nm}


        op_te_list = []
        ##for each te in the final_te_opt_dic
        for eachline in final_te_opt_dic:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            op_chr = col[0]
            op_start = col[1]
            op_end = col[2]
            op_nm = col[3]
            op_ori_nm = col[15]


            if col[4] == 'C':
                op_dir = '-'
            else:
                op_dir = col[4]

            op_te_dic = {'chr':op_chr,'bg':op_start,'ed':op_end,"dir":op_dir,'nm':op_nm,'ori_nm':op_ori_nm}

            ##compare with the re_te
            ##chr should be the same
            if re_chr == op_chr:
                ##tes should be covered each other
                if int(op_end) >= int(re_start) and int(op_start) <= int(re_end):
                    op_te_list.append(op_te_dic)

        #print ('the opt_dir is ' + str(op_dir))
        all_dic[re_nm] = op_te_list


    return (all_dic,re_te_dic)




#########################
##Step 4 cluster analysis
#########################
##define a function to store the infor of opt TEs that cover the re TEs: one opt TE covers multiple rel TEs
##this function will be used in classify_TEs module
def store_cover_TE (te_rating_type_dic,op_nm,type,type_element):
    if op_nm in te_rating_type_dic.keys():
        #print('the op_nm is in te_rating_type_dic')
        if type in te_rating_type_dic[op_nm].keys():
            te_rating_type_dic[op_nm][type] = te_rating_type_dic[op_nm][type] + ';' + type_element
        else:
            te_rating_type_dic[op_nm][type] = type_element

    else:
        ##if te_rating_type_dic has no key op_nm, create an empty dic
        #print('the op_nm is not in the te_rating_type_dic')
        te_rating_type_dic[op_nm] = {}
        te_rating_type_dic[op_nm][type] = type_element



##define a function to calculate the number of TEs that cover the re TEs
def calculate_cover_TE (te_cover_num_dic,op_nm):
    if op_nm in te_cover_num_dic.keys():
        te_cover_num_dic[op_nm] += 1
    else:
        te_cover_num_dic[op_nm] = 1


##updation3.16: remove this function (do not use this function)
##define a function to re collect opt TEs
def re_collect_opt_TE (te_rating_type_dic,eachopt_te,ipt_model_tp,opt_te_rating_dic):


    ##opt_te_rating_dic store the final opt
    ##if one opt TE covers multiple rel TEs, we will choose the best one using this function,
    ##else, we will use opt_te_rating_dic generated in the following function

    ##get the order of this rating type
    col = te_rating_type_dic[eachopt_te]['model_type'].split(';')
    order = col.index(ipt_model_tp)

    ##get the model type order
    model_type = ipt_model_tp
    ##extact other type
    col = te_rating_type_dic[eachopt_te]['class_type'].split(';')
    class_type = col[order]
    col = te_rating_type_dic[eachopt_te]['redundance'].split(';')
    redundance = col[order]
    col = te_rating_type_dic[eachopt_te]['re_te'].split(';')
    re_te = col[order]
    opt_te_rating_dic[eachopt_te] = {'class_type': class_type,
                                     'model_type': model_type,
                                     'redundance': redundance,
                                     're_te': re_te}



##define a function to classify TEs.
def classify_TEs (all_dic,re_te_dic):

    ##all_dic contains information about a real te nm and its covered op te list
    ##re_te_dic contains a real te and its relative information
    ##te_full_annot_dic contains real te and its family name


    ##initiate a dic to store rating and its opt name and store the TP and FP
    opt_te_rating_dic = {}

    ##initiate a dic to store no cover TE and store the FN
    no_cover_re_te_dic = {}

    ##initiate a dic to store te rating type
    ##the key is the op_nm and the value is the lists of each classification
    te_rating_type_dic = {}
    ##initiate a dic to calculate the number of opt TEs covering to the real TE
    te_cover_num_dic = {}



    real_te_count = 0
    ##for each real te in the all_dic
    for eachte in all_dic:

        ##if the all_dic is not empty
        ##if opt_te cover to the real te
        if all_dic[eachte] != []:

            real_te_count += 1
            print('the analzyed te is ' + str(real_te_count))

            ##define the re_te information before evaluation
            re_nm = re_te_dic[eachte]['nm']
            re_bg = re_te_dic[eachte]['bg']
            re_ed = re_te_dic[eachte]['ed']
            re_dir = re_te_dic[eachte]['dir']


            ##get the fam of the real te (one of the 20 fasta file)
            mt = re.match('(.+TE\d+)_.+',re_nm)
            #mt = re.match('(.+)_\d+$',re_nm)
            re_fam = mt.group(1)

            print('re_fam is ' + re_fam)

            #print('the re te is ' + eachte + ' the re te dir is ' + re_dir + '\n')

            ##do not divide the number of TEs in the all_dic
            for each_op_te_dic in all_dic[eachte]:

                ##initiate a dic to store rating type for the tes
                ##use the rating level as the key
                ##this dic is to solve the problem that some TEs cover multiple reference TEs
                ##So we need to decide a best TE rating type

                ##obtain the information for these opt TEs
                op_bg = each_op_te_dic['bg']
                op_ed = each_op_te_dic['ed']
                op_nm = each_op_te_dic['nm']
                op_dir = each_op_te_dic['dir']
                #op_ori_nm = each_op_te_dic['ori_nm']  ##because the original name is the last col of the final_opt_tb

                ##get the fam of the opt te (one of the 20 fasta file)
                ##because the repeatmodeler or ltr-retriever have different name
                mt = re.match('(.+)_TE\d+_\d+$', op_nm)
                op_fam = mt.group(1)

                ##if the direct is the same
                if op_dir == re_dir:

                    ##if the fam name is the same
                    if 'SINE' in op_fam:

                        ##update4.29: if there is a cover on two sequences, which will be regarded as TP
                        model_type = 'TP'

                        ##papare the parameters

                        ##select the larger groups and lower group for the ref and opt

                        ##for the union_right and inter_right
                        if int(re_ed) >= int(op_ed):
                            union_right = int(re_ed)
                            inter_right = int(op_ed)
                        else:
                            union_right = int(op_ed)
                            inter_right = int(re_ed)

                        ##for the union_left and inter_left
                        if int(re_bg) >= int(op_bg):
                            union_left = int(op_bg)
                            inter_left = int(re_bg)
                        else:
                            union_left = int(re_bg)
                            inter_left = int(op_bg)

                        union_len = union_right - union_left
                        inter_len = inter_right - inter_left + 1 ##forbid to be zero

                        BIAS = union_len/inter_len
                        Accuracy = 1/BIAS

                        #Accuracy = 1

                        ##detect the model type
                        ##if accuracy >= 0.5, the TE will be regarded as TP
                        ##if accuracy < 0.5, the TE will be regarded as FP
                        ##generate the model_type information that will be stored in the following analysis
                        #model_type = ''
                        #if float(Accuracy) >= 0.8:
                        #    model_type = 'TP'

                        #else:
                        #    model_type = 'FP'

                        ##compare with the real TEs
                        ##the following will remove the rating type
                        ##only detect class_type
                        ##redundance

                        ##if start and end region
                        if int(op_bg) == int(re_bg) and int(op_ed) == int(re_ed):

                            class_type = 'Start and End'
                            redundance = 'left:' + str(0) + ',right:' + str(0)

                            ##store the information
                            opt_te_rating_dic[op_nm] = {'class_type': class_type,
                                                        'redundance': redundance,
                                                        'model_type': model_type,
                                                        're_te':eachte,
                                                        'accuracy':str(Accuracy)}

                            calculate_cover_TE(te_cover_num_dic, op_nm)
                            store_cover_TE(te_rating_type_dic, op_nm, 'class_type', class_type)
                            store_cover_TE(te_rating_type_dic, op_nm, 'model_type', model_type)
                            store_cover_TE(te_rating_type_dic, op_nm, 'redundance', redundance)
                            store_cover_TE(te_rating_type_dic, op_nm, 're_te', eachte)

                            ##store the accuracy
                            store_cover_TE(te_rating_type_dic, op_nm, 'accuracy', str(Accuracy))


                        ##if not perfect
                        else:

                            ##if opt TEs are within real TEs
                            #if int(op_bg) >= int(re_bg) and int(op_ed) <= int(re_ed):

                            ##if op_bg is equal to re_bg
                            if int(op_bg) == int(re_bg) and int(op_ed) < int(re_ed):

                                class_type = 'Start'
                                left = 0
                                right = '-' + str(int(re_ed) - int(op_ed))
                                redundance = 'left:' + str(left) + ',right:' + str(right)
                                ##store the information
                                opt_te_rating_dic[op_nm] = {'class_type': class_type,
                                                            'model_type': model_type,
                                                            'redundance': redundance,
                                                            're_te': eachte,
                                                            'accuracy': str(Accuracy)}

                                calculate_cover_TE(te_cover_num_dic, op_nm)
                                store_cover_TE(te_rating_type_dic, op_nm, 'class_type', class_type)
                                store_cover_TE(te_rating_type_dic, op_nm, 'model_type', model_type)
                                store_cover_TE(te_rating_type_dic, op_nm, 'redundance', redundance)
                                store_cover_TE(te_rating_type_dic, op_nm, 're_te', eachte)

                                ##store the accuracy
                                store_cover_TE(te_rating_type_dic, op_nm, 'accuracy', str(Accuracy))

                            ##if op_ed is equal to re_ed
                            if int(op_ed) == int(re_ed) and int(op_bg) > int(re_bg):

                                class_type = 'End'
                                right = 0
                                left = '+' + str(int(op_bg) - int(re_bg))
                                redundance = 'left:' + str(left) + ',right:' + str(right)
                                ##store the information
                                opt_te_rating_dic[op_nm] = {'class_type': class_type,
                                                            'model_type': model_type,
                                                            'redundance': redundance,
                                                            're_te': eachte,
                                                            'accuracy': str(Accuracy)}

                                calculate_cover_TE(te_cover_num_dic, op_nm)
                                store_cover_TE(te_rating_type_dic, op_nm, 'class_type', class_type)
                                store_cover_TE(te_rating_type_dic, op_nm, 'model_type', model_type)
                                store_cover_TE(te_rating_type_dic, op_nm, 'redundance', redundance)
                                store_cover_TE(te_rating_type_dic, op_nm, 're_te', eachte)

                                ##store the accuracy
                                store_cover_TE(te_rating_type_dic, op_nm, 'accuracy', str(Accuracy))

                            ##if op_te is totally within
                            if int(op_bg) > int(re_bg) and int(op_ed) < int(re_ed):

                                class_type = 'Totally within'
                                left = '+' + str(int(op_bg) - int(re_bg))
                                right = '-' + str(int(re_ed) - int(op_ed))
                                redundance = 'left:' + str(left) + ',right:' + str(right)
                                ##store the information
                                opt_te_rating_dic[op_nm] = {'class_type': class_type,
                                                            'model_type': model_type,
                                                            'redundance': redundance,
                                                            're_te': eachte,
                                                            'accuracy': str(Accuracy)}

                                calculate_cover_TE(te_cover_num_dic, op_nm)
                                store_cover_TE(te_rating_type_dic, op_nm, 'class_type', class_type)
                                store_cover_TE(te_rating_type_dic, op_nm, 'model_type', model_type)
                                store_cover_TE(te_rating_type_dic, op_nm, 'redundance', redundance)
                                store_cover_TE(te_rating_type_dic, op_nm, 're_te', eachte)

                                ##store the accuracy
                                store_cover_TE(te_rating_type_dic, op_nm, 'accuracy', str(Accuracy))

                            ##if opt_te is overlaps start
                            if int(op_bg) < int(re_bg) and int(op_ed) <= int(re_ed):
                                class_type = 'Overlaps start'
                                left = '-' + str(int(re_bg) - int(op_bg))
                                right = ''
                                if (int(re_ed) - int(op_ed)) == 0:
                                    right = str(0)
                                else:
                                    right = '-' + str(int(re_ed) - int(op_ed))
                                redundance = 'left:' + str(left) + ',right:' + str(right)
                                ##store the information
                                opt_te_rating_dic[op_nm] = {'class_type': class_type,
                                                            'model_type': model_type,
                                                            'redundance': redundance,
                                                            're_te': eachte,
                                                            'accuracy': str(Accuracy)}

                                calculate_cover_TE(te_cover_num_dic, op_nm)
                                store_cover_TE(te_rating_type_dic, op_nm, 'class_type', class_type)
                                store_cover_TE(te_rating_type_dic, op_nm, 'model_type', model_type)
                                store_cover_TE(te_rating_type_dic, op_nm, 'redundance', redundance)
                                store_cover_TE(te_rating_type_dic, op_nm, 're_te', eachte)

                                ##store the accuracy
                                store_cover_TE(te_rating_type_dic, op_nm, 'accuracy', str(Accuracy))

                            ##if opt_te is the overlaps end
                            if int(op_bg) >= int(re_bg) and int(op_ed) > int(re_ed):

                                class_type = 'Overlaps end'
                                left = ''
                                if (int(op_bg) - int(re_bg)) == 0:
                                    left = str(0)
                                else:
                                    left = '+' + str(int(op_bg) - int(re_bg))
                                right = '+' + str(int(op_ed) - int(re_ed))
                                redundance = 'left:' + str(left) + ',right:' + str(right)
                                ##store the information
                                opt_te_rating_dic[op_nm] = {'class_type': class_type,
                                                            'model_type': model_type,
                                                            'redundance': redundance,
                                                            're_te': eachte,
                                                            'accuracy': str(Accuracy)}

                                calculate_cover_TE(te_cover_num_dic, op_nm)
                                store_cover_TE(te_rating_type_dic, op_nm, 'class_type', class_type)
                                store_cover_TE(te_rating_type_dic, op_nm, 'model_type', model_type)
                                store_cover_TE(te_rating_type_dic, op_nm, 'redundance', redundance)
                                store_cover_TE(te_rating_type_dic, op_nm, 're_te', eachte)

                                ##store the accuracy
                                store_cover_TE(te_rating_type_dic, op_nm, 'accuracy', str(Accuracy))



                            ##if opt_te is the totally within
                            if int(op_bg) < int(re_bg) and int(op_ed) > int(re_ed):

                                class_type = 'Totally over'
                                left = '-' + str(int(re_bg) - int(op_bg))
                                right = '+' + str(int(op_ed) - int(re_ed))
                                redundance = 'left:' + str(left) + ',right:' + str(right)
                                ##store the information
                                opt_te_rating_dic[op_nm] = {'class_type': class_type,
                                                            'model_type': model_type,
                                                            'redundance': redundance,
                                                            're_te': eachte,
                                                            'accuracy': str(Accuracy)}

                                calculate_cover_TE(te_cover_num_dic, op_nm)
                                store_cover_TE(te_rating_type_dic, op_nm, 'class_type', class_type)
                                store_cover_TE(te_rating_type_dic, op_nm, 'model_type', model_type)
                                store_cover_TE(te_rating_type_dic, op_nm, 'redundance', redundance)
                                store_cover_TE(te_rating_type_dic, op_nm, 're_te', eachte)

                                ##store the accuracy
                                store_cover_TE(te_rating_type_dic, op_nm, 'accuracy', str(Accuracy))



                    ##if fam name is not the same
                    else:
                        print ('the te ' + op_nm  + ' has no consistent name for the te ref')
                        class_type = 'NA'
                        model_type = 'FP'
                        left = 'NA'
                        right = 'NA'
                        redundance = 'left:' + str(left) + ',right:' + str(right)

                        Accuracy = float(0)
                        ##store the information
                        opt_te_rating_dic[op_nm] = {'class_type': class_type,
                                                    'model_type': model_type,
                                                    'redundance': redundance,
                                                    're_te': eachte,
                                                    'accuracy':str(Accuracy)}

                        calculate_cover_TE(te_cover_num_dic, op_nm)
                        store_cover_TE(te_rating_type_dic, op_nm, 'class_type', class_type)
                        store_cover_TE(te_rating_type_dic, op_nm, 'model_type', model_type)
                        store_cover_TE(te_rating_type_dic, op_nm, 'redundance', redundance)
                        store_cover_TE(te_rating_type_dic, op_nm, 're_te', eachte)

                        ##store the accuracy
                        store_cover_TE(te_rating_type_dic, op_nm, 'accuracy', str(Accuracy))


                ##if the direct is not the same
                else:

                    #print ('the te ' + op_nm + 'has no consistent direction for the te ref' + '\n')
                    class_type = 'NA'
                    model_type = 'FP'
                    left = 'NA'
                    right = 'NA'
                    redundance = 'left:' + str(left) + ',right:' + str(right)

                    Accuracy = float(0)
                    ##store the information
                    opt_te_rating_dic[op_nm] = {'class_type': class_type,
                                                'model_type': model_type,
                                                'redundance': redundance,
                                                're_te': eachte,
                                                'accuracy': str(Accuracy)}

                    calculate_cover_TE(te_cover_num_dic, op_nm)
                    store_cover_TE(te_rating_type_dic, op_nm, 'class_type', class_type)
                    store_cover_TE(te_rating_type_dic, op_nm, 'model_type', model_type)
                    store_cover_TE(te_rating_type_dic, op_nm, 'redundance', redundance)
                    store_cover_TE(te_rating_type_dic, op_nm, 're_te', eachte)

                    ##store the accuracy
                    store_cover_TE(te_rating_type_dic, op_nm, 'accuracy', str(Accuracy))


        ##it appears no use
        ##if no opt_te cover to the real TE, it will be regarded as FN
        else:
            class_type = 'No cover'
            model_type = 'FN'
            ##store the information
            no_cover_re_te_dic[eachte] = {'class_type': class_type,
                                          'model_type': model_type}


    ##comment the follow script since the accuracy does not work

    ##consider the multiple opt TEs covering one real TE and choose the best status for the opt TE
    ##filter the name
    #print (te_cover_num_dic)
    #print (te_rating_type_dic)
    #for eachopt_te in te_cover_num_dic:
    #    if int(te_cover_num_dic[eachopt_te]) > 1:

            ##for this situation, it will compare accuracy for each opt TE and choose the highest one,
            ##the highest one that relates to the rating type is the model type

    #        col_acc = te_rating_type_dic[eachopt_te]['accuracy'].split(';')

    #        longest = float(0)
    #        for eachnum in col_acc:
    #            if float(eachnum) >= longest:
    #                longest = float(eachnum)

            ##get the order of the longest accuracy
    #        order = col_acc.index(str(longest))

            ##get the all the types information
    #        col = te_rating_type_dic[eachopt_te]['model_type'].split(';')
    #        model_type = col[order]
    #        col = te_rating_type_dic[eachopt_te]['class_type'].split(';')
    #        class_type = col[order]
    #        col = te_rating_type_dic[eachopt_te]['redundance'].split(';')
    #        redundance = col[order]
    #        col = te_rating_type_dic[eachopt_te]['re_te'].split(';')
    #        re_te = col[order]

            ##store the information into the opt_te_rating_dic, this will overlap the te stored in the previous step
     #       opt_te_rating_dic[eachopt_te] = {'class_type': class_type,
     #                                        'model_type': model_type,
     #                                        'redundance': redundance,
     #                                        're_te': re_te,
     #                                        'accuracy': str(longest)}


    return (opt_te_rating_dic, no_cover_re_te_dic)


##############################################
##define a function to calculate the TN and FN
##############################################
##first, identify the regions that not covered by real TEs.
##then, detect if these regions contains full FP, if yes, this region will be FN, otherwise, this region will be TN.
def identify_TN_FN (real_te_line_dic,final_te_opt_dic):

    ############################################
    ##step 1 store the inter region in a bed dic
    ##identify all the N regions
    ##generate a bed file for the inter sequence
    bed_line_dic = {}

    ##generate a list to contain all the location information
    location_list = []

    ##get the line count for the file
    line_count = len(list(real_te_line_dic.keys()))

    ##get all the sequence number from the fasta
    seq_len = 0
    for seq_record in SeqIO.parse(input_fasta_file, 'fasta'):
        seq_len = len(str(seq_record.seq))

    count = 0
    for eachline in real_te_line_dic:
        count += 1

        eachline = eachline.strip('\n')
        col = eachline.strip().split()
        start = col[1]
        end = col[2]

        if count == 1:

            ##from the begin to the the first TE
            location_list.append('1')
            location_list.append(start)
            location_list.append(end)

        else:

            ##if meet the final line
            if count == line_count:

                location_list.append(start)
                location_list.append(end)
                location_list.append(str(seq_len))

            else:
                location_list.append(start)
                location_list.append(end)

    if len(location_list) % 2 == 0:
        #print('the file has even line count')
        pair_num = len(location_list) / 2
        for i in range(0, int((pair_num))):
            # print('the range max is ' + str(int(pair_num) + 2))
            if int(location_list[2 * i + 1]) > int(location_list[2 * i]):
                bed_line = 'seq' + '\t' + location_list[2 * i] + '\t' + location_list[2 * i + 1] + '\t' + 'seq' + str(
                    i) + '\t' + '1' + '\t' + '+'
                bed_line_dic[bed_line] = 1

    #############################
    ##step 2 detect the FN and TN
    ##initiate a dic to store the TN and FN
    FN = 0
    TN = 0
    for inter_eachline in bed_line_dic:

        inter_col = inter_eachline.split()
        inter_st = inter_col[1]
        inter_ed = inter_col[2]

        ##calculate the number of FP in the inter_eachline
        N = 0
        for opt_eachline in final_te_opt_dic:

            col = opt_eachline.split()
            op_start = col[1]
            op_end = col[2]

            ##if there is opt_te in the inter_eachline, and this line will be considered as TN
            ##tes should be covered each other

            if int(op_start) >= int(inter_st) and int(op_end) <= int(inter_ed):
                N += 1
                break

        ##If there is a FP in the inter region
        if N == 1:
            FN += 1
        else:
            TN += 1

    ##The FN in this part indicates there are some FP in this region
    N_dic = {'FN':str(FN),'TN':str(TN)}

    return (N_dic)

###############################
##Step 4: write out the summary
###############################
##updation3.16: generate the average of accuracy, and remove the previous rating system.
##updation3.16: generate a file storing the FP, that is not covered with rel TEs. define another to see if there is domain on it.
##if there are domains for the FP, so we can regard these FP is novel TE.
def construct_summary (final_te_opt_dic,all_dic,opt_te_rating_dic,N_dic):

    ##N_dic will help to calculate TN and FN

    ##initiate a final line
    final_line_dic = {}

    ##calculate the total op_te
    total_opt_te_num = len(list(final_te_opt_dic.keys()))

    ##calculate all total re_te
    total_re_te_num = len(list(all_dic.keys()))

    ##calculate the covered op_te
    cover_opt_te_num = len(list(opt_te_rating_dic.keys()))

    ##calculate not covered op_te that will be the FP
    not_cover_opt_te_num = int(total_opt_te_num) - int(cover_opt_te_num)


    ##########################################
    ##calculate the number of seven class type
    ##initiate a dic to store number of seven class type for opt te
    class_type_opt_te_num_dic = {}
    for each_opt_te in opt_te_rating_dic:
        if opt_te_rating_dic[each_opt_te]['class_type'] in class_type_opt_te_num_dic.keys():
            class_type_opt_te_num_dic[opt_te_rating_dic[each_opt_te]['class_type']] += 1
        else:
            class_type_opt_te_num_dic[opt_te_rating_dic[each_opt_te]['class_type']] = 1


    ################################
    ##calculate total TN, FN, TP, FP
    ##FOR TN
    TN_num = int(N_dic['TN'])
    not_TN_num = int(N_dic['FN'])

    ##FOR FN
    ##calculate the covered re_te
    ##calculate FN number for the re_te also is the no covered re_te
    FN_num = 0
    cover_re_te_num = 0
    for eachte in all_dic:
        if all_dic[eachte] != []:
            cover_re_te_num += 1
        else:
            FN_num += 1

    ##FOR TP
    TP_num = 0
    for each_opt_te in opt_te_rating_dic:
        if opt_te_rating_dic[each_opt_te]['model_type'] == 'TP':
            TP_num += 1

    ##FOR FP
    FP_num = 0
    for each_opt_te in opt_te_rating_dic:
        if opt_te_rating_dic[each_opt_te]['model_type'] == 'FP':
            FP_num += 1

    ##this FP should add the FP that will be not covered by ref
    FP_num = FP_num + not_cover_opt_te_num

    ##############################
    ##calculate the SE, PR, SP, AC
    SE = TP_num / (TP_num + FN_num)
    PR = TP_num / (TP_num + FP_num)
    SP = TN_num / (TN_num + FP_num)
    AC = (TP_num + TN_num) / (TP_num + TN_num + FP_num + FN_num)

    #######################################
    ##calculate the average accuracy for TP
    total_accuracy_TP = 0
    for each_opt_te in opt_te_rating_dic:
        if opt_te_rating_dic[each_opt_te]['model_type'] == 'TP':
            accuracy = opt_te_rating_dic[each_opt_te]['accuracy']
            total_accuracy_TP = total_accuracy_TP + float(accuracy)
    average_accuracy = total_accuracy_TP/TP_num

    ######################################
    ##calculate the median accuracy for TP
    accuracy_TP_list = []
    for each_opt_te in opt_te_rating_dic:
        if opt_te_rating_dic[each_opt_te]['model_type'] == 'TP':
            accuracy = opt_te_rating_dic[each_opt_te]['accuracy']
            accuracy_TP_list.append(float(accuracy))
    acc_median = median(accuracy_TP_list)

    #######################################
    ##calculate TP number with 80% accuracy
    TP_num_over_08 = 0
    for each_opt_te in opt_te_rating_dic:
        if opt_te_rating_dic[each_opt_te]['model_type'] == 'TP':
            accuracy = opt_te_rating_dic[each_opt_te]['accuracy']
            if float(accuracy) >= float(0.8):
                TP_num_over_08 += 1

    #######################################
    ##calculate TP number with 50% accuracy
    TP_num_over_05 = 0
    for each_opt_te in opt_te_rating_dic:
        if opt_te_rating_dic[each_opt_te]['model_type'] == 'TP':
            accuracy = opt_te_rating_dic[each_opt_te]['accuracy']
            if float(accuracy) >= float(0.5):
                TP_num_over_05 += 1


    ########################################################################
    ##detect if other opt_te that is not covered by rel TE has domain or not
    ##do not consider the novel_te_count = 0
    #novel_te_count = 0

    ##############################################
    ##Generate the final line for the total number
    final_line = 'Total reference TE number' + '\t' + str(total_re_te_num)
    final_line_dic[final_line] = 1
    final_line = 'Detected reference TE number' + '\t' + str(cover_re_te_num)
    final_line_dic[final_line] = 1
    final_line = 'Total detected TE number' + '\t' + str(total_opt_te_num)
    final_line_dic[final_line] = 1
    ##add detected TE number in reference = all TP number

    #final_line = 'Total ref TE number' + '\t' + str(total_re_te_num)
    #final_line_dic[final_line] = 1
    #final_line = 'Covered opt TE number' + '\t' + str(cover_opt_te_num)
    #final_line_dic[final_line] = 1
    #final_line = 'Covered ref TE number' + '\t' + str(cover_re_te_num)
    #final_line_dic[final_line] = 1

    ##Gnerate the te_TP and te_FP number
    final_line = 'TP_num' + '\t' + str(TP_num)
    final_line_dic[final_line] = 1

    final_line = 'TP_average_accuracy' + '\t' + str(average_accuracy)
    final_line_dic[final_line] = 1

    final_line = 'TP_median_accuracy' + '\t' + str(acc_median)
    final_line_dic[final_line] = 1

    final_line = 'TP_num_over_0.8_accuracy' + '\t' + str(TP_num_over_08)
    final_line_dic[final_line] = 1

    final_line = 'TP_num_over_0.5_accuracy' + '\t' + str(TP_num_over_05)
    final_line_dic[final_line] = 1

    final_line = 'FP_num' + '\t' + str(FP_num)
    final_line_dic[final_line] = 1
    #final_line = 'Novel_TE_num' + '\t' + str(novel_te_count)
    #final_line_dic[final_line] = 1
    final_line = 'TN_num' + '\t' + str(TN_num)
    final_line_dic[final_line] = 1
    final_line = 'FN_num' + '\t' + str(FN_num)
    final_line_dic[final_line] = 1
    #final_line = 'not_TN_num' + '\t' + str(not_TN_num)  ##this will be removed after testing
    #final_line_dic[final_line] = 1


    final_line = 'SE' + '\t' + str(SE)
    final_line_dic[final_line] = 1
    final_line = 'PR' + '\t' + str(PR)
    final_line_dic[final_line] = 1
    final_line = 'SP' + '\t' + str(SP)
    final_line_dic[final_line] = 1
    final_line = 'AC' + '\t' + str(AC)
    final_line_dic[final_line] = 1


    ##Generate the final line for the class
    for eachclass in class_type_opt_te_num_dic:
        final_line = eachclass + '\t' + str(class_type_opt_te_num_dic[eachclass])
        final_line_dic[final_line] = 1


    return(final_line_dic)

##write out results
real_te_line_dic = store_real_TE (input_real_TE_loc_file)
final_te_opt_dic = store_final_te_opt (input_opt_te_file)
all_dic,re_te_dic = store_te_and_final_te (real_te_line_dic,final_te_opt_dic)
opt_te_rating_dic, no_cover_re_te_dic = classify_TEs (all_dic,re_te_dic)
N_dic = identify_TN_FN (real_te_line_dic,final_te_opt_dic)
final_line_dic = construct_summary (final_te_opt_dic,all_dic,opt_te_rating_dic,N_dic)


with open (input_opt_dir + '/opt_tmpt_evaluate_tb.txt','w+') as opt:
    for eachline in final_line_dic:
        opt.write(eachline + '\n')
