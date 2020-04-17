#!/usr/bin/env python


##update 1.9 change the combined name of ltr predicted from ltr software
##update 1.28 update the store_infor_gff function, add another function to add the id for each te, because we found one r
##cannot recognize id

##This script is to combine the ltr results into the gff results generated from repeatmasker

import re

##the inputfile
##01: output from the repeatmasker
##02: all the ltr information from the complete ltr results which could be generated from opt_compete_time_cul_tb but without time


##########################################################
##Step 1: create a dic to store all the te ltr information
def store_ltr (input_ltr_file):
    ltr_dic = {}
    with open (input_ltr_file, 'r') as ipt_ltr:
        for eachline in ipt_ltr:
            col = eachline.strip().split()
            ltr_nm = col[1]
            chr = col[0]
            start = col[2]
            end = col[3]
            direct = col[4]
            type = col[10]
            ##store the information into the dic
            ltr_dic[ltr_nm] = {'chr':chr,'bg':start,'ed':end,'dir':direct,'type':type}

    return (ltr_dic)

#######################################################
##Step 2: store the te information from the te_gff file
##do not match the line starting with SW and score and do not match the empty line
def store_infor_gff (input_te_gff):
    dic_te = {}
    rmk_dic = {} ##store the final line from the repeatmasker

    id_count = 0  ##update 1.28 add the id for each te

    right_line_list = []  ##update 1.29
    with open(input_te_gff, 'r') as ipt_gff:
        right_line_list = ipt_gff.readlines()[3:]

    
    for eachelem in right_line_list:

        id_count += 1
        id = str(id_count)
        col = eachelem.split()
        lib_left = ''
        lib_bg = ''
        #id = col[14]
        if re.match('C', col[8]):
            lib_left = col[13]
            mt = re.match('\((.+)\)', col[11])
            lib_bg = mt.group(1)

        if re.match('\+', col[8]):
            lib_bg = col[11]
            mt = re.match('\((.+)\)', col[13])
            lib_left = mt.group(1)

        dic_te[id] = {'chr': col[4], 'begin': col[5], 'end': col[6], 'direct': col[8],
                      'name': col[9], 'lib_bg': lib_bg, 'lib_ed': col[12], 'lib_left': lib_left}

        rmk_line = col[4] + '\t' + col[5] + '\t' + col[6] + '\t' + col[8] + '\t' + col[9] + '\t' + lib_bg + '\t' + \
                   col[12] + '\t' +  lib_left

        rmk_dic[rmk_line] = 1

    return (dic_te,rmk_dic)


#####################################################
##Step 3: extract TEs that belong to the range of ltr
def extract_te_under_ltr (ltr_dic,dic_te):
    ##initial a dic to store the ltr and its subset information
    ltr_te_dic = {}
    for eachltr in ltr_dic:
        ##set eachltr as key to search the gff file
        ##initial a list to store the candidate id which is covered into teach ltr
        id_list = []

        for eachid in dic_te:

            ##if the chr is equal that will increase the speed of analysis
            if ltr_dic[eachltr]['chr'] == dic_te[eachid]['chr']:
                #print(dic_te[eachid]['chr'])
                ##if the direction is the same
                direct_gff = ''
                if dic_te[eachid]['direct'] == 'C':
                    direct_gff = '-'
                if dic_te[eachid]['direct'] == '+':
                    direct_gff = '+'

                if ltr_dic[eachltr]['dir'] == direct_gff:

                    ##some tes in the border region is also concluded in the ltr
                    ##set the range of 100 bp
                    ltr_bg = int(ltr_dic[eachltr]['bg']) - 100
                    ltr_ed = int(ltr_dic[eachltr]['ed']) + 100


                    ##if the id in the gff is under the range of each ltr in the ltr_dic
                    if int(dic_te[eachid]['begin']) >= ltr_bg and int(dic_te[eachid]['end']) <= ltr_ed:
                        ##store the potential id to the list
                        id_list.append(eachid)

        ##store the list to the ltr dictionary
        ltr_te_dic[eachltr] = id_list

    #print(ltr_te_dic)

    return (ltr_te_dic)

######################################################################
##Step 4: detect whether the id in each ltr has the same name and type
##initial a function that will be used for the following filteration
def store_filtered_ltr (type,type_count,ltr_te_dic,ltr_dic,filter_ltr_dic,eachte,removed_ids_list):
    ##example: if all the subset is the Gypsy so the type should be Gypsy
    if type_count == len(ltr_te_dic[eachte]):
        ##if the annotate ltr_te from the domain is gypsy
        if re.match('.*RLG.*', ltr_dic[eachte]['type']):
            ##because the type has the RLG_C and RLG_N, so it needs to use re.match
            ##store this ltr
            filter_ltr_dic[eachte] = 1
            ##store the removed subset id
            for eachid in ltr_te_dic[eachte]:
                removed_ids_list.append(eachid)

        ##if the annotate ltr_te from the domain is copia
        if re.match('.*RLC.*', ltr_dic[eachte]['type']):
            ##think the RLC is right, store this ltr
            filter_ltr_dic[eachte] = 1
            ##store the removed subset id
            for eachid in ltr_te_dic[eachte]:
                removed_ids_list.append(eachid)

        ##if the annotate ltr_te from the domain is RLX
        if re.match('.*RLX.*', ltr_dic[eachte]['type']):

            ##the te type in the ltr_dic should be changed to Gypsy
            if type == 'Gypsy':
                ltr_dic[eachte]['type'] = 'RLG'   ##it should check if the 'RLG' has been successfully changed
                ##store this ltr
                filter_ltr_dic[eachte] = 1
                ##store the removed subset id
                for eachid in ltr_te_dic[eachte]:
                    removed_ids_list.append(eachid)

            if type == 'Copia':
                ltr_dic[eachte]['type'] = 'RLC'
                ##store this ltr
                filter_ltr_dic[eachte] = 1
                ##store the removed subset id
                for eachid in ltr_te_dic[eachte]:
                    removed_ids_list.append(eachid)

            if type == 'Other':
                ltr_dic[eachte]['type'] = 'Other'
                ##store this ltr
                filter_ltr_dic[eachte] = 1
                ##store the removed subset id
                for eachid in ltr_te_dic[eachte]:
                    removed_ids_list.append(eachid)
    ##no return

def filter_te (ltr_te_dic,dic_te,ltr_dic):
    ##initiate a dic to contain all the filtered ltr that will be combined with gff file
    filter_ltr_dic = {}
    ##initiate a list that contains ids that will be removed out from the gff file
    removed_ids_list = []

    count = 0
    for eachte in ltr_te_dic:

        count += 1
        print('the analzyed te number is ' + str(count))

        ##if the id list is empty, suggesting this ltr is inserted in the gff
        if len(ltr_te_dic[eachte]) == 0:
            filter_ltr_dic[eachte] = 1

        #############################
        ##if the id list is not empty
        if len(ltr_te_dic[eachte]) != 0:
            count_ltr = 0  ##this count will be used to count id that is LTR
            for eachid in ltr_te_dic[eachte]:
                ##if all the covered te is LTR
                if re.match('.*LTR.*',dic_te[eachid]['name']):
                    count_ltr += 1

            ###########################
            ##if all the subset is LTRs
            if count_ltr == len(ltr_te_dic[eachte]):
                count_gyp = 0  ##this count will be used to count Gypsy count
                count_cop = 0  ##this count will be used to count Copia count
                count_other = 0  ##this count will be used to count no Gypsy and no Copia
                for eachid in ltr_te_dic[eachte]:
                    ##there will be different situation
                    ##if these subset LTR are all Gypsy and the type for this ltr is Gypsy
                    if re.match('.*Gypsy.*',dic_te[eachid]['name'],flags=re.I, ):
                        count_gyp += 1
                    ##if these subset LTR are all Gypsy and the type for this ltr is Gypsy
                    if re.match('.*Copia.*',dic_te[eachid]['name'],flags=re.I, ):
                        count_cop += 1
                    ##if these subset LTR are LTR others
                    if not re.match('.*Gypsy.*',dic_te[eachid]['name'],flags=re.I, ) and not re.match('.*Copia.*',dic_te[eachid]['name'],flags=re.I, ):
                        count_other += 1

                ################################
                ##if all the subset is the Gypsy
                store_filtered_ltr('Gypsy', count_gyp, ltr_te_dic, ltr_dic, filter_ltr_dic, eachte, removed_ids_list)
                ##if all the subset is the Copia
                store_filtered_ltr('Copia', count_cop, ltr_te_dic, ltr_dic, filter_ltr_dic, eachte,removed_ids_list)
                ##if all the subset is others
                store_filtered_ltr('Other', count_other, ltr_te_dic, ltr_dic, filter_ltr_dic, eachte, removed_ids_list)

            ##if not all the subset is LTRs
            ##for the tes in this subset will not consider
            ##So do not need to store the ltr and store the removed ids list
    #print('the removed id list is following')
    #print(removed_ids_list)
    return (filter_ltr_dic,removed_ids_list)

########################################################
##Step 5: remove the ids and add the LTR in the gff file
##store all the information
def store_all_info (dic_te,removed_ids_list,filter_ltr_dic,ltr_dic):
    ##import dic_te contains all the te information
    ##import filter_ltr_dic contains the ltr that should be inserted into the gff
    ##import removed_ids_list contains the ids that should be removed from the gff file
    ##import ltr_dic contains all the ltr information

    ##initiate a dic to store the final opt
    final_line_dic = {}

    ##remove the id from the dic_te
    for eachid in dic_te:

        if not eachid in removed_ids_list:

            ##generate a final line
            final_line = dic_te[eachid]['chr'] + '\t' + dic_te[eachid]['begin'] + '\t' + dic_te[eachid]['end'] + '\t' + \
                         dic_te[eachid]['direct'] + '\t' + dic_te[eachid]['name'] + '\t' + dic_te[eachid]['lib_bg'] + '\t' + \
                         dic_te[eachid]['lib_ed'] + '\t' + dic_te[eachid]['lib_left']

            final_line_dic[final_line] = 1

    ##update1.11
    ##add the ltr te in the gff
    for eachte in filter_ltr_dic:

        ##divide different situation to update the name of LTRs
        if 'RLG' in ltr_dic[eachte]['type']:

            final_line = ltr_dic[eachte]['chr'] + '\t' + ltr_dic[eachte]['bg'] + '\t' + ltr_dic[eachte]['ed'] + '\t' + \
                         ltr_dic[eachte]['dir'] + '\t' + 'ClassI//LTR/Gypsy_TE_' + eachte + '_combined' + '\t' + \
                         'na' + '\t' + 'na' + '\t' + 'na'

            final_line_dic[final_line] = 1

        if 'RLC' in ltr_dic[eachte]['type']:

            final_line = ltr_dic[eachte]['chr'] + '\t' + ltr_dic[eachte]['bg'] + '\t' + ltr_dic[eachte]['ed'] + '\t' + \
                         ltr_dic[eachte]['dir'] + '\t' + 'ClassI//LTR/Copia_TE_' + eachte + '_combined' + '\t' + \
                         'na' + '\t' + 'na' + '\t' + 'na'

            final_line_dic[final_line] = 1

        if 'FU' in ltr_dic[eachte]['type'] or 'RLX' in ltr_dic[eachte]['type']:

            final_line = ltr_dic[eachte]['chr'] + '\t' + ltr_dic[eachte]['bg'] + '\t' + ltr_dic[eachte]['ed'] + '\t' + \
                         ltr_dic[eachte]['dir'] + '\t' + 'ClassI//LTR/ltr_Others_TE_' + eachte + '_combined' + '\t' + \
                         'na' + '\t' + 'na' + '\t' + 'na'

            final_line_dic[final_line] = 1


    return(final_line_dic)



