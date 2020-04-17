#!/usr/bin/env python

##update3.12 to exchange APE to EN, since other ENs will be added
##update3.11 to exchange EN to APE
##update3.4 to add the AP domain to the evalution
##update1.7 to change the EnSpm and MuDR to the no mite DNA and add the final original name for each te

##This script is to combine table with connection, match rate, and score information to get the final table
##And calculate the score of each te

import re
import statistics as s

##input file 1: te connection, match rate information from the 04_01_conect_tp_mt
##input file 2: te domain information from 04_02_domain

##############################
##step 1: store te information
##############################
##initiate a function to store mrt_conect infor
def store_mrt_conect (mrt_conect_file):
    mrt_conect_dic = {}
    with open (mrt_conect_file,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            te_nm = col[3]
            mrt = col[-1]
            mrt_conect_dic[te_nm] = {'all':eachline,'mrt':mrt}
    return (mrt_conect_dic)

##initiate a function to store domain infor
def store_domain (domain_file):
    domain_dic = {}
    with open (domain_file,'r') as ipt:
        for eachline in ipt:
            col = eachline.strip().split()
            te_nm = col[0]
            te_pat = col[1]
            te_pro = col[2]
            domain_dic[te_nm] = {'pat': te_pat, 'pro': te_pro}
    return (domain_dic)


#######################################
##step 2: combination of domain and mrt
#######################################
def combine_dm_mrt (mrt_conect_dic,domain_dic):

    ##initiate a dic store the combined line
    combine_dic = {}

    for eachte in mrt_conect_dic:
        combine_line = ''
        if eachte in domain_dic.keys():
            combine_line = mrt_conect_dic[eachte]['all'] + '\t' + domain_dic[eachte]['pat'] + '\t' + \
                           domain_dic[eachte]['pro']
        else:
            combine_line = mrt_conect_dic[eachte]['all'] + '\t' + 'na' + '\t' + 'na'
        combine_dic[combine_line] = 1

    return (combine_dic)


#############################
##step 3: calculate the score
#############################
##initiate a function to convert str to float
def convert(list):
    new_list = []
    for i in list:
        new_list.append(float(i))
    return new_list

##initate a function to calculate the score
def calculate_score (combine_dic):

    ##initate a dic to store the final line
    final_dic = {}

    count = 0
    for eachline in combine_dic:
        count += 1
        print ('the analyzed line is ' + str(count))

        col = eachline.strip().split()
        matchrt = ''
        if col[11] != 'na':
            matchrt = float(col[11])
        tenm = col[3]
        dm_compt = col[-1]
        score = ''

        #########################################################################
        ##it should consider use the combined situation: ltrfinder and ltrharvest
        if re.match('.*combined.*', tenm, flags=re.I):

            if dm_compt == 'na':
                score = 'na'  ##because there is match rate and this te lack all the domains
            ##if the domain competeness is not na
            else:
                dm_col = dm_compt.strip().split(',')
                ##tes may have multiple domains, and we need to calculate the average comp of each type of domain
                domain_dic = {}  ##initiate a dic to store domain name and its relative match rate
                for eachdm in dm_col:
                    mt = re.match('(.+)\:(.+)', eachdm)
                    dm_nm = mt.group(1)
                    mt_rt = mt.group(2)
                    if dm_nm in domain_dic.keys():
                        domain_dic[dm_nm] = domain_dic[
                                                dm_nm] + '\t' + mt_rt  ##use int to let the calculation more easier
                    else:
                        domain_dic[dm_nm] = mt_rt
                ##calculate the average for each type of domain
                m_GAG = 0
                m_RT = 0
                m_IN = 0
                m_RH = 0

                ##update 3.4
                m_AP = 0

                for dm in domain_dic:
                    if dm == 'GAG':
                        compt_list = domain_dic[dm].split()
                        new_list = convert(compt_list)
                        m_GAG = s.mean(new_list)
                    if dm == 'RT':
                        compt_list = domain_dic[dm].split()
                        new_list = convert(compt_list)
                        m_RT = s.mean(new_list)
                    if dm == 'IN':
                        compt_list = domain_dic[dm].split()
                        new_list = convert(compt_list)
                        m_IN = s.mean(new_list)
                    if dm == 'RH':
                        compt_list = domain_dic[dm].split()
                        new_list = convert(compt_list)
                        m_RH = s.mean(new_list)

                    ##update 3.4
                    if dm == 'AP':
                        compt_list = domain_dic[dm].split()
                        new_list = convert(compt_list)
                        m_AP = s.mean(new_list)

                ##update 3.4 change m_RT to 60
                score = (m_RT * 60 + m_IN * 10 + m_RH * 10 + m_AP * 10 + m_GAG * 10) / 100

        ##if not the ltr_combined
        else:
        ###########################
        ##if the te name is the LTR
            if re.match('.*LTR.*', tenm, flags=re.I):

                ##if the domain competeness na
                if dm_compt == 'na':
                    score = (matchrt * 10) / 100
                ##if the domain competeness is not na
                else:
                    dm_col = dm_compt.strip().split(',')
                    ##tes may have multiple domains, and we need to calculate the average comp of each type of domain
                    domain_dic = {}  ##initiate a dic to store domain name and its relative match rate
                    for eachdm in dm_col:
                        mt = re.match('(.+)\:(.+)', eachdm)
                        dm_nm = mt.group(1)
                        mt_rt = mt.group(2)
                        if dm_nm in domain_dic.keys():
                            domain_dic[dm_nm] = domain_dic[
                                                    dm_nm] + '\t' + mt_rt  ##use int to let the calculation more easier
                        else:
                            domain_dic[dm_nm] = mt_rt
                    ##calculate the average for each type of domain
                    m_GAG = 0
                    m_RT = 0
                    m_IN = 0
                    m_RH = 0

                    ##update 3.4
                    m_AP = 0


                    for dm in domain_dic:
                        if dm == 'GAG':
                            compt_list = domain_dic[dm].split()
                            new_list = convert(compt_list)
                            m_GAG = s.mean(new_list)
                        if dm == 'RT':
                            compt_list = domain_dic[dm].split()
                            new_list = convert(compt_list)
                            m_RT = s.mean(new_list)
                        if dm == 'IN':
                            compt_list = domain_dic[dm].split()
                            new_list = convert(compt_list)
                            m_IN = s.mean(new_list)
                        if dm == 'RH':
                            compt_list = domain_dic[dm].split()
                            new_list = convert(compt_list)
                            m_RH = s.mean(new_list)

                        ##update 3.4
                        if dm == 'AP':
                            compt_list = domain_dic[dm].split()
                            new_list = convert(compt_list)
                            m_AP = s.mean(new_list)

                    ##update 3.4 change m_RT to 50
                    score = (m_RT * 50 + m_IN * 10 + m_RH * 10 + m_AP * 10 + m_GAG * 10 + matchrt * 10) / 100

            ############################
            ##if the te name is the LINE
            if re.match('.*LINE.*', tenm, flags=re.I):
                ##if the domain competeness na
                if dm_compt == 'na':
                    score = (matchrt * 10) / 100

                ##if the domain competeness is not na
                else:
                    dm_col = dm_compt.strip().split(',')
                    ##tes may have multiple domains, and we need to calculate the average comp of each type of domain
                    domain_dic = {}  ##initiate a dic to store domain name and its relative match rate
                    for eachdm in dm_col:
                        mt = re.match('(.+)\:(.+)', eachdm)
                        dm_nm = mt.group(1)
                        mt_rt = mt.group(2)
                        if dm_nm in domain_dic.keys():
                            domain_dic[dm_nm] = domain_dic[
                                                    dm_nm] + '\t' + mt_rt  ##use int to let the calculation more easier
                        else:
                            domain_dic[dm_nm] = mt_rt
                    ##calculate the average for each type of domain
                    m_RT = 0
                    m_EN = 0

                    for dm in domain_dic:
                        if dm == 'RT':
                            compt_list = domain_dic[dm].split()
                            new_list = convert(compt_list)
                            m_RT = s.mean(new_list)

                        ##update3.12: exchange 'APE' to 'EN'
                        ##update3.11: exchange 'EN' to 'APE'
                        if dm == 'EN':
                            compt_list = domain_dic[dm].split()
                            new_list = convert(compt_list)
                            m_EN = s.mean(new_list)

                    score = (m_RT * 50 + m_EN * 30 + matchrt * 20) / 100

            #######################################
            ##if the te name is the DNA with domain

            ##update 1.7 change the EnSpm to the DNA_nMITE
            if re.match(".*DNA_nMITE.*", tenm, flags=re.I):
                ##if the domain competeness na
                if dm_compt == 'na':
                    score = (matchrt * 30) / 100

                ##if the domain competeness is not na
                else:
                    dm_col = dm_compt.strip().split(',')
                    ##tes may have multiple domains, and we need to calculate the average comp of each type of domain
                    domain_dic = {}  ##initiate a dic to store domain name and its relative match rate
                    for eachdm in dm_col:
                        mt = re.match('(.+)\:(.+)', eachdm)
                        dm_nm = mt.group(1)
                        mt_rt = mt.group(2)
                        if dm_nm in domain_dic.keys():
                            domain_dic[dm_nm] = domain_dic[
                                                    dm_nm] + '\t' + mt_rt  ##use int to let the calculation more easier
                        else:
                            domain_dic[dm_nm] = mt_rt
                    ##calculate the average for each type of domain
                    m_TR = 0

                    for dm in domain_dic:
                        if dm == 'TR':
                            compt_list = domain_dic[dm].split()
                            new_list = convert(compt_list)
                            m_TR = s.mean(new_list)

                    score = (m_TR * 70 + matchrt * 30) / 100

            #################
            ##if match others
            if not re.match('.*LTR.*', tenm, flags=re.I) and not re.match('.*LINE.*', tenm,flags=re.I) \
                    and not re.match(".*DNA_nMITE.*",tenm,flags=re.I):
                score = (matchrt * 100) / 100

        ######################################
        ##store the eachline to the final_line
        final_line = eachline + '\t' + str(score)
        final_dic[final_line] = 1

    return (final_dic)


###################################################################
##step 4: add the orginal name of TEs generated from the TE library
###################################################################
##define a function to add the original TE name
def add_org_nm (final_dic,name_file):

    ##initiate a dic to store the final dic
    final_orgnm_dic = {}

    ##initiate a dic to store the name from old new name file
    new_old_nm_dic = {}
    ##open the name_file to extract the old and new name
    with open (name_file, 'r') as ipt_nm:
        for eachline in ipt_nm:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            new_old_nm_dic[col[0]] = col[1]

    for eachline in final_dic:
        eachline = eachline.strip('\n')
        col = eachline.strip().split('\t')
        full_name = col[3]
        ##get the reduced name
        mt = re.match('(.+)_\d+$', full_name)
        re_nm = mt.group(1)
        if re_nm in new_old_nm_dic.keys():
            new_line = eachline + '\t' + new_old_nm_dic[re_nm]
            final_orgnm_dic[new_line] = 1
        else:
            new_line = eachline + '\t' + 'na'
            final_orgnm_dic[new_line] = 1

    return (final_orgnm_dic)











