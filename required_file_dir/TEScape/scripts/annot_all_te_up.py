#!/usr/bin/env python

##script updation: 3.12, remove the id.seq 1-6 and id.out 1-6 to make sure that the elder one will not disturb the current detection
##script updation: 3.12, change the APE to EN, because we will add other EN, and APE will also name EN
##do not update filter, since if there is a LTR domain in another family, which will be interesting, script updation: 3.12, add a new function to remove the domain name that will not be involved in the specific family
##script updation: 3.11, exchange EN to APE domain for the LINE
##script updation: 2.28, add AP domain for the LTR analysis
##script updation: 1.6 and not revise in the old version

##This script is to detect the domain pattern for each tes with domains
##This script is to detect the completeness of the pattern for each te in order to calculate the score

##################
##import functions
##################

##import function of translation
from Bio.SeqUtils import six_frame_translations
import re
from Bio import SeqIO
import subprocess
import pandas as pd
import os.path



##############################################################################
##Step 2: Translate each reading frame and identify the possible reading frame
##############################################################################

##update the translate for six reading frame for the no mite DNA
##define a reverse complement functions
def getRC(te_seq_rc):
    out = te_seq_rc.reverse_complement()
    return out

##define a function to translate all the six reading frame
def translate_no_mite (te_seq):
    ##initial a dictionary to store the protein sequence for the six reading frame
    pro_seq_dic = {}
    ##Runs a six-frame translation for each gene in the string ts2seq
    f6 = six_frame_translations(te_seq)
    ##Splits each gene in the above mentioned translation with a new line
    f6tmp = str(f6).split('\n')
    num_m = re.match('(.+)5$', str(len(f6tmp)))
    num_cnt = num_m.group(1)
    ####################################################################
    ##the following is the first three reading frame of this te sequence
    ##get the location information of sequences
    l1 = 5
    l2 = 6
    l3 = 7
    ##Initializes a hand full of lists we will later fill with the genomic data to analyze
    list1 = []
    list2 = []
    list3 = []
    ##Sets range that the genomic data we are looking at can be based on the 6-frame translation
    for i in range(0, (int(num_cnt))):
        l1 = 5 + i * 10
        list1.append(l1)
        l2 = 6 + i * 10
        list2.append(l2)
        l3 = 7 + i * 10
        list3.append(l3)
    ##extract the sequence information from the location informationa
    ##the following is the first three list
    allist_f = [list1, list2, list3]
    ##create a string to be the index for the key of the dictionary
    fr_count = 0
    for eachlist in allist_f:
        fr_count = fr_count + 1
        ##extrac the sequence information
        ##For each gene, adds the word 'list' for each additional gene after the start + end codon
        framenm = 'f' + str(fr_count)
        ##initial a list to store all the protein sequences
        f6aa = []
        # Setting the frame protein for each translation
        for i in eachlist:
            pseq = f6tmp[i]
            f6aa.append(pseq)
        ##combination of the sequence
        comptlist = ' '.join(f6aa)
        ##split the '' between each aa and then replace it with ,
        comptstr = comptlist.replace(' ', '')
        ##store sequence information on the dictionary
        pro_seq_dic[framenm] = comptstr

    ###############################################################
    ##the following is last three reading frame of this te sequence

    ##Inputs the data to get the reversed compliment
    seq_reverse = getRC(te_seq)
    f6 = six_frame_translations(seq_reverse)
    f6tmp = str(f6).split('\n')
    num_m = re.match('(.+)5$', str(len(f6tmp)))
    num_cnt = num_m.group(1)
    l4 = 5
    l5 = 6
    l6 = 7
    list4 = []
    list5 = []
    list6 = []
    for i in range(0, (int(num_cnt))):
        l4 = 5 + i * 10
        list4.append(l4)
        l5 = 6 + i * 10
        list5.append(l5)
        l6 = 7 + i * 10
        list6.append(l6)
    allist_f = [list4, list5, list6]
    re_count = 3
    for eachlist in allist_f:
        re_count = re_count + 1
        framenm = 'f' + str(re_count)
        f6aa = []
        for i in eachlist:
            pseq = f6tmp[i]
            f6aa.append(pseq)
        comptlist = ' '.join(f6aa)
        comptstr = comptlist.replace(' ', '')
        pro_seq_dic[framenm] = comptstr

    return (pro_seq_dic)


##intial a translation function
##only translate the first three reading frame because it should keep consensus with the PBS binding sites
def translate_ltr (te_seq):
    ##initial a dictionary to store the protein sequence for the six reading frame
    pro_seq_dic = {}
    ##Runs a six-frame translation for each gene in the string ts2seq
    f6 = six_frame_translations(te_seq)
    ##Splits each gene in the above mentioned translation with a new line
    f6tmp = str(f6).split('\n')
    num_m = re.match('(.+)5$', str(len(f6tmp)))
    num_cnt = num_m.group(1)
    ####################################################################
    ##the following is the first three reading frame of this te sequence
    ##get the location information of sequences
    l1 = 5
    l2 = 6
    l3 = 7
    ##Initializes a hand full of lists we will later fill with the genomic data to analyze
    list1 = []
    list2 = []
    list3 = []
    ##Sets range that the genomic data we are looking at can be based on the 6-frame translation
    for i in range(0, (int(num_cnt))):
        l1 = 5 + i * 10
        list1.append(l1)
        l2 = 6 + i * 10
        list2.append(l2)
        l3 = 7 + i * 10
        list3.append(l3)
    ##extract the sequence information from the location informationa
    ##the following is the first three list
    allist_f = [list1, list2, list3]
    ##create a string to be the index for the key of the dictionary
    fr_count = 0
    for eachlist in allist_f:
        fr_count = fr_count + 1
        ##extrac the sequence information
        ##For each gene, adds the word 'list' for each additional gene after the start + end codon
        framenm = 'f' + str(fr_count)
        ##initial a list to store all the protein sequences
        f6aa = []
        # Setting the frame protein for each translation
        for i in eachlist:
            pseq = f6tmp[i]
            f6aa.append(pseq)
        ##combination of the sequence
        comptlist = ' '.join(f6aa)
        ##split the '' between each aa and then replace it with ,
        comptstr = comptlist.replace(' ', '')
        ##store sequence information on the dictionary
        pro_seq_dic[framenm] = comptstr

    return (pro_seq_dic)

##initial a function to identify the domain for the six reading frame
##define a function to identify te domain and generate the domain table
##import the pro_seq_dic generated from the translate function
def iden_domain_six_rdframe (pro_seq_dic,seqid,hmmscan_exe,working_dir):

    ##initial a dictionary to store the results
    te_domain_six_rdframe_dic = {}

    ##initial a dictionary to store the name of six reading frame
    six_frame_nm_dic = {}
    for eachframe in pro_seq_dic:
        te_frame_nm = seqid + '_' + eachframe
        te_opt_nm = eachframe + '.seq'
        te_pro_seq = open(working_dir + '/' + te_opt_nm, 'w')
        ##this will generate each reading frame file
        te_pro_seq.write('>' + te_frame_nm + '\n' + str(pro_seq_dic[eachframe]))
        te_pro_seq.close()  ##it is important to close the te_pro_seq file
        six_frame_nm_dic[eachframe] = 1

    ##use the hmm to detect domain information for each frame
    for eachnm in six_frame_nm_dic:
        tefm_seq = eachnm + '.seq'
        opt_tefm = eachnm + '.out'
        cmd = hmmscan_exe + ' ' + working_dir + '/minifam ' + working_dir + '/' + tefm_seq + ' > ' + working_dir + '/' + opt_tefm
        subprocess.call(cmd, shell=True)

        ##open the output file from hmmscan
        domain_te_file = open(working_dir + '/' + opt_tefm, 'r')
        for line in domain_te_file:
            line.strip()
            if re.match('.+\d+\s+PF\d+_.+', line):
                mt = re.match('(.+PF\d+_[A-Z]+).+', line)
                new_line = mt.group(1)
                te_line = seqid + '_' + eachnm + '\t' + new_line
                te_domain_six_rdframe_dic[te_line] = 1

    return (te_domain_six_rdframe_dic)

#############################################################
##the following funcitons is to filter the best reading frame
##define a function to store the te name and its domain number information in the dictionary
def cal_dm_num (te_nm,domain_dtf): ##note: the domain_dtf is generated by the next function.
    ##creat a dictionary to store number information for each domain
    domain_num_dic = {}
    ##collect the target te
    ##filter the different domain name for the ltr
    te_ltr_rt_dt = ''
    te_ltr_in_dt = ''
    te_ltr_gag_dt = ''
    te_ltr_rh_dt = ''
    te_line_rt_dt = ''

    ##update:3.11 exc APE
    te_line_en_dt = ''
    te_dna_tr_dt = ''

    ##update: 2.28 add 'AP'
    te_ltr_ap_dt = ''


    if re.match('.*LTR.*', te_nm):
        ##because the domain_dtf cotains multiple te_nm and filter the specific name like only one readingframe
        te_dt = domain_dtf.loc[domain_dtf[domain_dtf.columns[0]] == te_nm]
        print('the following is the te_dt')
        print(te_dt)

        te_ltr_rt_dt = te_dt[te_dt[te_dt.columns[9]].str.contains('RT', flags=re.I, na=False)]
        te_ltr_in_dt = te_dt[te_dt[te_dt.columns[9]].str.contains('IN', flags=re.I, na=False)]
        te_ltr_gag_dt = te_dt[te_dt[te_dt.columns[9]].str.contains('GAG', flags=re.I, na=False)]
        te_ltr_rh_dt = te_dt[te_dt[te_dt.columns[9]].str.contains('RH', flags=re.I, na=False)]

        ##update: 2.28 add 'AP'
        te_ltr_ap_dt = te_dt[te_dt[te_dt.columns[9]].str.contains('AP', flags=re.I, na=False)]

    ##filter the LINE
    ##update: 3.12 exchange 'APE' to 'EN'
    ##update: 3.11 exchange 'EN' to 'APE'
    if re.match('.*LINE.*',te_nm):
        te_dt = domain_dtf.loc[domain_dtf[domain_dtf.columns[0]] == te_nm]
        te_line_rt_dt = te_dt[te_dt[te_dt.columns[9]].str.contains('RT', flags=re.I, na=False)]
        te_line_en_dt = te_dt[te_dt[te_dt.columns[9]].str.contains('EN', flags=re.I, na=False)]


    ##filter the DNA_nMITE
    if re.match(".*DNA_nMITE.*", te_nm):
        te_dt = domain_dtf.loc[domain_dtf[domain_dtf.columns[0]] == te_nm]
        te_dna_tr_dt = te_dt[te_dt[te_dt.columns[9]].str.contains('TR', flags=re.I, na=False)]

    ##update: 2.28 add 'AP'
    domain_num_dic[te_nm] = {'RT_ltr': len(te_ltr_rt_dt), 'IN_ltr': len(te_ltr_in_dt), 'GAG_ltr': len(te_ltr_gag_dt),
                             'RH_ltr': len(te_ltr_rh_dt), 'AP_ltr': len(te_ltr_ap_dt), 'RT_line': len(te_line_rt_dt),
                             'EN_line': len(te_line_en_dt),'TR_nomt': len(te_dna_tr_dt)}

    return domain_num_dic

##define a function to calculate the number of te_type
##this function will be used in the compre_rdframe
def cal_num (te_type):
    if len(str(te_type)) != 0:
        num = te_type
    else:
        num = 0
    return (num)

##intial a function to compare each reading frame and choose the best one
##the second output of the compare_rdframe function is to store the te domain number information
##import the pro_seq_dic from the translate function
def compare_rdframe (te_domain_six_rdframe_dic,working_dir):
    ##the import dic is the te_domain_six_rdframe_dic from the previous iden_domain_six_rdframe

    #########################################################################
    ##Step 1: write out because the next step will use the output of this dic
    with open(working_dir + '/opt_temp_six_rdframe_domain_hmm.tb', 'w+') as opt_te_ltr:
        for eachline in te_domain_six_rdframe_dic:
            opt_te_ltr.write(eachline + '\n')
    ##read the temp_domain_hmm table
    domain_df = pd.read_table(working_dir + '/opt_temp_six_rdframe_domain_hmm.tb', header=None, delimiter=r"\s+")
    te_nm = domain_df[domain_df.columns[0]].unique() ##contain information for six reading frame name
    ##the df dataframe only store te domain information for 0 or 1 that are used for filter the longest reading frame
    ##update: 2.28 add 'AP'
    df = pd.DataFrame(index=te_nm, columns=['RT_LTR', 'IN_LTR', 'GAG_LTR', 'RH_LTR','AP_LTR','RT_LINE','EN_LINE','TR_DNA'])
    ##transfer the domain_df to the domain_dtf
    domain_dtf = pd.DataFrame(domain_df)

    ##get the reduced name
    reduced_nm = ''
    ################################################
    ##Step 2: create the table of domain information
    for eachte in te_nm:
        ##call the domain dictionary from function cal_dm_num
        te_domain_num_dic = cal_dm_num(eachte,domain_dtf)
        # print(te_domain_num_dic[eachte]['RT_ltr'])
        ##save the number of domain to the table
        ##for ltr#######################
        if te_domain_num_dic[eachte]['RT_ltr'] != 0:
            df.loc[eachte, 'RT_LTR'] = 1
        else:
            df.loc[eachte, 'RT_LTR'] = 0
        if te_domain_num_dic[eachte]['IN_ltr'] != 0:
            df.loc[eachte, 'IN_LTR'] = 1
        else:
            df.loc[eachte, 'IN_LTR'] = 0
        if te_domain_num_dic[eachte]['GAG_ltr'] != 0:
            df.loc[eachte, 'GAG_LTR'] = 1
        else:
            df.loc[eachte, 'GAG_LTR'] = 0
        if te_domain_num_dic[eachte]['RH_ltr'] != 0:
            df.loc[eachte, 'RH_LTR'] = 1
        else:
            df.loc[eachte, 'RH_LTR'] = 0

        ##update: 2.28 add 'AP'
        if te_domain_num_dic[eachte]['AP_ltr'] != 0:
            df.loc[eachte, 'AP_LTR'] = 1
        else:
            df.loc[eachte, 'AP_LTR'] = 0


        ##for line######################
        if te_domain_num_dic[eachte]['RT_line'] != 0:
            df.loc[eachte, 'RT_LINE'] = 1
        else:
            df.loc[eachte, 'RT_LINE'] = 0
        if te_domain_num_dic[eachte]['EN_line'] != 0:
            df.loc[eachte, 'EN_LINE'] = 1
        else:
            df.loc[eachte, 'EN_LINE'] = 0
        ##for the nomit####################
        if te_domain_num_dic[eachte]['TR_nomt'] != 0:
            df.loc[eachte,'TR_DNA'] = 1
        else:
            df.loc[eachte, 'TR_DNA'] = 0

        ##get the reduced_nm
        mt = re.match('(.+)_f\d', eachte)
        reduced_nm = mt.group(1)



    ###############################################################################
    ##step 3: analyze each te with all reading frames and get the filtered te table
    ##create an empty dataframe to store the filtered information
    #df_ft = pd.DataFrame(index=set(te_nm_list), columns=['RT_LTR', 'IN_LTR', 'GAG_LTR', 'RH_LTR'])
    ##because df only contains six dataframe so do not need to create a new dt to store only dt for one dataframe
    max_te_nm = ''

    ltr_rt = ''
    ltr_in = ''
    ltr_rh = ''
    ltr_gag = ''
    line_rt = ''
    line_en = ''
    dna_tr = ''

    ##update: 2.28 add 'AP'
    ltr_ap = ''

    #print('the te_nm is '+ te_nm)
    ##for the ltr
    if re.match('.*LTR.*', reduced_nm):  ##te_nm is a list to contain reading frame information, and te_nm[1] means one of them,and contain type information

        ##update: 2.28 add 'AP'
        df['ltr_sum'] = df[['RT_LTR', 'IN_LTR', 'RH_LTR', 'AP_LTR']].sum(axis=1)
        ##only select the best te
        te_ltr_rt_dt = df[df['RT_LTR'] != 0]
        if (len(te_ltr_rt_dt) != 0):
            max_te_nm = te_ltr_rt_dt['ltr_sum'].idxmax()
            ltr_rt = df.loc[max_te_nm, 'RT_LTR']
            ltr_in = df.loc[max_te_nm, 'IN_LTR']
            ltr_rh = df.loc[max_te_nm, 'RH_LTR']
            ltr_ap = df.loc[max_te_nm, 'AP_LTR']
        else:
            max_te_nm = df['ltr_sum'].idxmax()
            ltr_rt = df.loc[max_te_nm, 'RT_LTR']
            ltr_in = df.loc[max_te_nm, 'IN_LTR']
            ltr_rh = df.loc[max_te_nm, 'RH_LTR']
            ltr_ap = df.loc[max_te_nm, 'AP_LTR']

        ##for the gag domain
        te_gag_nz_dt = df[df['GAG_LTR'] != 0]
        if (len(te_gag_nz_dt) != 0):
            ltr_gag = 1
        else:
            ltr_gag = 0

    ##for the line
    if re.match('.*LINE.*', reduced_nm):

        df['line_sum'] = df[['RT_LINE', 'EN_LINE']].sum(axis=1)
        ##only select the best te
        te_line_rt_dt = df[df['RT_LINE'] != 0]
        if (len(te_line_rt_dt) != 0):
            max_te_nm = te_line_rt_dt['line_sum'].idxmax()
            line_rt = df.loc[max_te_nm,'RT_LINE']
            line_en = df.loc[max_te_nm,'EN_LINE']
        else:
            max_te_nm = df['line_sum'].idxmax()
            line_rt = df.loc[max_te_nm,'RT_LINE']
            line_en = df.loc[max_te_nm,'EN_LINE']


    ##for the DNA_nMITE
    if re.match(".*DNA_nMITE.*", reduced_nm):

        df['no_mite_sum'] = df[['TR_DNA']].sum(axis=1)
        ##only select the best te
        te_dna_tr_dt = df[df['TR_DNA'] != 0]
        if (len(te_dna_tr_dt) != 0):
            max_te_nm = te_dna_tr_dt['no_mite_sum'].idxmax()
            dna_tr = df.loc[max_te_nm, 'TR_DNA']
        else:
            max_te_nm = df['no_mite_sum'].idxmax()
            dna_tr = df.loc[max_te_nm, 'TR_DNA']

    #########################################################
    ##step 4: store the number information for each max_te_nm
    ##get the reduced nm
    mt = re.match('(.+)_f\d',max_te_nm)
    domain_reduced_nm = mt.group(1)

    domain_num_dic = {} ##this dic store the number for the max_te_nm

    domain_num_dic[domain_reduced_nm] = {'RT_ltr': cal_num(ltr_rt), 'IN_ltr': cal_num(ltr_in),
                                         'GAG_ltr': cal_num(ltr_gag), 'RH_ltr': cal_num(ltr_rh),
                                         'AP_ltr':cal_num(ltr_ap),'RT_line': cal_num(line_rt),
                                         'EN_line': cal_num(line_en),'TR_nomt': cal_num(dna_tr)}

    return(max_te_nm,domain_num_dic)


#################################################
##Step 3: Identify the domain pattern for each te
#################################################
##the step 1 we obtain the best te, and this step we will go to each the best te seq to remove the covered te
##initial a function to obtain the seq name of max te
def iden_max_seq (max_te_nm):
    mt = re.match('.+_(f\d)',max_te_nm)
    rfnm = mt.group(1)
    seq_out = rfnm + '.out'
    return (seq_out)

##initial a function collect the output of the hmm
def collect_domain (seq_out_file):
    ##return the all information
    domain_location_dic = {}

    ##store the domain nm and its number information
    domain_nm_num_dic = {}

    with open(seq_out_file,'r') as ipt_seq:
        for eachline in ipt_seq:
            eachline.strip()
            if re.match('.+\d+\s+PF\d+_.+', eachline):
                col = eachline.strip().split()
                number = col[7]
                #print(number)
                domain_nm = col[8]
                #print(domain_nm)
                domain_nm_num_dic[domain_nm] = number

    ##store the lines to one dic
    seq_out_line_list = []
    with open(seq_out_file,'r') as ipt_seq:
        for eachline in ipt_seq:
            eachline = eachline.strip('\n')
            seq_out_line_list.append(eachline)

    ##initial a dic to store the domain location information
    #domain_location_dic = {}
    for eachline in seq_out_line_list:
        if '>>' in eachline:
            col = eachline.strip().split()
            te_domain_nm = col[1]
            te_domain_num = domain_nm_num_dic[te_domain_nm]
            ##get the line number
            line_th = seq_out_line_list.index(eachline)
            ##start loc number
            start = line_th + 3
            end = start + int(te_domain_num)

            ##extract the target information
            for i in range(start,end):
                tar_line = te_domain_nm + '\t' + seq_out_line_list[i]
                domain_location_dic[tar_line] = 1

    return (domain_location_dic)

##initial a function to remove the coverage of TE domain
def store_domain_info (domain_location_dic,working_dir):
    ##import file is the domain_nm_num_dic in the collect_domain function
    ##initial a domain location dic to store the domain information
    ##use the domain name as the last column in this dic
    domain_all_dic = {}

    new_domain_location_dic = {}

    #dm_nm_dic = {}

    for eachline in domain_location_dic:
        col = eachline.strip().split()
        ##filter out low evalue
        if float(col[5]) < 0.05 and float(col[6]) < 0.05:
            ##extract the domain name
            mt = re.match('(.+)_(.+)', col[0])
            dm_nm = mt.group(2)

            new_line = eachline + '\t' + dm_nm
            print(new_line)
            new_domain_location_dic[new_line] = 1

    ##write out the results as the temp file that will be as the input file for the next function which will use the pandas
    with open(working_dir + '/opt_temp_seq_domain_loc.txt', 'w+') as opt_dm_loc:
        for eachline in new_domain_location_dic:
            opt_dm_loc.write(eachline + '\n')

    with open(working_dir + '/opt_temp_seq_domain_loc.txt') as ipt_test:
        first = ipt_test.read(1)
        if not first:
            print('this file is empty')
        else:
            print('this file is not empty')
            ltr_fl = pd.read_table(working_dir + '/opt_temp_seq_domain_loc.txt', header=None, delimiter=r"\s+")
            ##sort chr and start region
            ltr_fl.columns = ['domain_name','1','2','3','4','5','6','7','8','9','start','11','12','13','14','15','16','type']
            ltr_sort_fl = ltr_fl.sort_values(by=['type','start'])
            ltr_sort_fl.to_csv(working_dir + '/opt_temp_seq_domain_loc_sort.csv', sep='\t')


            ##updation3.12: some temp seq will be empty after filter out the low evalue

            ##store the type information
            ##this step is to conduct a loop for analyzing each domain type
            domain_dic = {}
            with open(working_dir + '/opt_temp_seq_domain_loc_sort.csv', 'r') as ipt_tb:
                for eachline in ipt_tb:
                    if not re.match('.*domain_name.+', eachline):
                        eachline = eachline.strip('\n')
                        col = eachline.split('\t')
                        # print(col[1])
                        domain_dic[col[-1]] = 1

            ##this step is to detect the last domain for the further analysis
            ##calculate te number in each chromosome
            domain_te_num_dic = {}
            with open(working_dir + '/opt_temp_seq_domain_loc_sort.csv', 'r') as ipt_tb:
                for eachline in ipt_tb:
                    if not re.match('.*domain_name.+', eachline):
                        eachline = eachline.strip('\n')
                        col = eachline.split('\t')
                        type_nm = col[-1]
                        if type_nm in domain_te_num_dic.keys():
                            domain_te_num_dic[type_nm] = domain_te_num_dic[type_nm] + 1
                        else:
                            domain_te_num_dic[type_nm] = 1

            ##this step is to connect the domain information in the dictionary
            ##initial a dictionary to store all the information of domain, and use the domain type name as the key
            ##this domain_all_dic has been initialled in the begin of this function
            ##initial a dictionary to store the id number
            type_id_dic = {}
            for eachtype in domain_dic:
                ##if the chr change to the next one, it should to initial an empty te_dic and an empty dic_list
                te_dic = {}
                dic_list = []
                with open(working_dir + '/opt_temp_seq_domain_loc_sort.csv', 'r') as ipt_tb:
                    for eachline in ipt_tb:
                        if not re.match('.*domain_name.+', eachline):

                            eachline = eachline.strip('\n')

                            col = eachline.split('\t')
                            type_nm = col[-1]

                            if type_nm == eachtype:
                                if eachtype in type_id_dic.keys():
                                    type_id_dic[eachtype] = type_id_dic[eachtype] + 1
                                else:
                                    type_id_dic[eachtype] = 1

                                id = str(type_id_dic[eachtype])
                                type_nm = col[1]
                                type_begin = col[11]
                                type_end = col[12]
                                ##add the start amd emd information of hmm
                                ##this step is important because it can calculate the length of mapping proportion.
                                hmm_begin = col[8]
                                hmm_end = col[9]

                                if str(domain_te_num_dic[eachtype]) == str(1):
                                    te_dic[id] = {'Type_nm': col[1], 'Type_begin': type_begin, 'Type_end': type_end, 'Hmm_begin':hmm_begin, 'Hmm_end':hmm_end}
                                    dic_list.append(te_dic)
                                else:
                                    if id == str(1):
                                        te_dic[id] = {'Type_nm': col[1], 'Type_begin': type_begin, 'Type_end': type_end, 'Hmm_begin':hmm_begin, 'Hmm_end':hmm_end}
                                    if id > str(1):
                                        if type_begin <= te_dic[str(int(id) - 1)]['Type_end']:
                                            te_dic[id] = {'Type_nm': col[1], 'Type_begin': type_begin, 'Type_end': type_end, 'Hmm_begin':hmm_begin, 'Hmm_end':hmm_end}
                                        else:
                                            dic_list.append(te_dic)
                                            te_dic = {}
                                            te_dic[id] = {'Type_nm': col[1], 'Type_begin': type_begin, 'Type_end': type_end, 'Hmm_begin':hmm_begin, 'Hmm_end':hmm_end}

                                        if id == str(domain_te_num_dic[eachtype]):
                                            dic_list.append(te_dic)

                domain_all_dic[eachtype] = dic_list  ##this dictionary contain the key which is the chr name and value which is the dic list

    return(domain_all_dic)


##the following function is to filter te and then get the hmm length information for each domain
##initial a function to transfer data and store in another dictionary
def transfer (id_order,te_comp_dic,te_dic):

    te_comp_dic[id_order] = {'Type_nm': te_dic[id_order]['Type_nm'], 'Type_begin': te_dic[id_order]['Type_begin'], 'Type_end': te_dic[id_order]['Type_end'],'Hmm_begin': te_dic[id_order]['Hmm_begin'], 'Hmm_end': te_dic[id_order]['Hmm_end']}

##initial a function to filter te covered by the same te
##the argument is the te which store all the information
##the domain_all_dic is the dic from the previou function store_domain_info
##this filter_te will generate the output of proportion so the input will contain the full length information for each domain
##initial a function to store the full length information for each domain
def get_domain_len (input_domain_file):
    ##initial a dic to store length for each domain
    domain_len_dic = {}
    with open (input_domain_file, 'r') as ipt_dm:
        for eachline in ipt_dm:
            col = eachline.strip().split()
            domain_len_dic[col[0]] = col[1]
    return(domain_len_dic)

def filter_te (domain_all_dic,domain_len_dic,working_dir):

    ##initial a string to store the the domain pattern
    ##this function will return the domain string
    domain_str = ''

    ##initial a dic to store the filtered te
    ##there are some problems to store the te line because the order is wrong so we use the list to store
    #te_line_list = []
    te_line_dic = {}

    ##te_all_dic's key is the chromosome dictionary of each store a list
    ##each list sotre multiple dictionary. Some dic contains one key and some dic contains multiple keys
    for eachchr_dic in domain_all_dic:
        dic_list = domain_all_dic[eachchr_dic]

        for tar_dic in dic_list: ##each dic store one key or multiple keys
            #tar_dic = te_all_dic[eachchr_dic][eachdic]  ##name a dic to a other name to decrease the length of dic

            if len(tar_dic) == 1:
                for eachid in tar_dic:  ##call the each key in this tar_dic
                    tar_te_line = tar_dic[eachid]['Type_nm'] + '\t' + tar_dic[eachid]['Type_begin'] + '\t' + \
                                  tar_dic[eachid]['Type_end'] + '\t' + tar_dic[eachid]['Hmm_begin'] + '\t' + \
                                  tar_dic[eachid]['Hmm_end']

                    #te_line_list.append(tar_te_line)
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
                        cur_len = int(tar_dic[eachid]['Type_end']) - int(tar_dic[eachid]['Type_begin'])

                        if eachid == list(tar_dic)[1]:

                            pre_len = 0
                            if str(int(eachid) - 1) in tar_dic.keys():
                            ##previous length
                                pre_len = int(tar_dic[str(int(eachid) - 1)]['Type_end']) - int(tar_dic[str(int(eachid) - 1)]['Type_begin'])

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
                            com_len = int(te_comp_dic[list(te_comp_dic)[-1]]['Type_end']) -  int(te_comp_dic[list(te_comp_dic)[-1]]['Type_begin'])
                            ##if there is cover between the current one and the compare dic one
                            if tar_dic[eachid]['Type_begin'] < te_comp_dic[list(te_comp_dic)[-1]]['Type_end']:
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
                ##important: store the domain begin and domain end information
                for eachid in te_comp_dic:
                    #print('the target id count is ' + str(count))
                    tar_te_line = te_comp_dic[eachid]['Type_nm'] + '\t' + te_comp_dic[eachid]['Type_begin'] + '\t' + \
                                  te_comp_dic[eachid]['Type_end'] + '\t' + te_comp_dic[eachid]['Hmm_begin'] + '\t' + \
                                  te_comp_dic[eachid]['Hmm_end']

                    te_line_dic[tar_te_line] = 1


    ##transfer the te_line_dic to a dt because we will sort the begin of the te
    with open (working_dir + '/opt_temp_domain_order.tb','w+') as opt_dm_or_dt:
        for eachline in te_line_dic:
            opt_dm_or_dt.write(eachline + '\n')

    ##use pandas open the dt and order the dt according to the Type_begin
    temp_dm_df = pd.read_table(working_dir + '/opt_temp_domain_order.tb', header=None, delimiter=r"\s+")
    temp_dm_df.columns = ['domain_name', 'begin', 'end','hmm_begin','hmm_end']
    temp_dm_sort_fl = temp_dm_df.sort_values(by=['begin'])
    temp_dm_sort_fl.to_csv(working_dir + '/opt_temp_domain_order_sort.csv', sep='\t')

    ##open the sorted temp file
    ##generate proportion information
    domain_pro_list = []
    domain_pat_list = []
    domain_pat_gag_list = []
    with open (working_dir + '/opt_temp_domain_order_sort.csv') as ipt_dm_sort:
        for eachline in ipt_dm_sort:
            if not re.match('.*domain_name.*', eachline):
                col = eachline.strip().split()
                mt = re.match('(.+)_(.+)', col[1])
                dm_nm = mt.group(2)

                ##generate te proportion including gag example: RT:0.3,IN:0.5,IN:0.3
                dm_len = int(col[5])-int(col[4]) + 1
                dm_full_len = int(domain_len_dic[col[1]])
                dm_pro = dm_len/dm_full_len
                dm_nm_pro = dm_nm + ':' + str(dm_pro)
                domain_pro_list.append(dm_nm_pro)

                ##generate te pattern including gag example: (GAG)RTININ or (GAG)
                if dm_nm != 'GAG':
                    domain_pat_list.append(dm_nm)
                else:
                    domain_pat_gag_list.append(dm_nm)

    domain_pro_str = ','.join(domain_pro_list)

    domain_pat_str = ''
    if len(domain_pat_gag_list) != 0:
        domain_pat_str = '(GAG)' + ''.join(domain_pat_list)
    else:
        domain_pat_str = ''.join(domain_pat_list)


    return (domain_pro_str,domain_pat_str)


######################################################################
##Step 4: use the functions to generate the domain pattern for each te
######################################################################

##updation3.12: define a function to remove the temp file
def remove_temp_file (working_dir,file_nm):
    ##check the file exit or not
    PATH = working_dir + '/' + file_nm
    if os.path.isfile(PATH):
        print('file exists')
        cmd = 'rm ' + working_dir + '/' + file_nm
        subprocess.call(cmd, shell=True)
        print(cmd)
    else:
        print('file not exist')


def generate_domain_pattern (input_ltr_fas,hmmscan_exe,input_domain_length,working_dir):
    ##import the input_ltr_fas from the arg in the begin of the script

    ##initial a dic to store all the information about the te and domain
    te_domain_pattern_dic = {}

    ##initial a dic to store the number of domain for each te
    domain_num_line_dic = {}

    ##initial a te_count to calculate the number of analyzed te
    te_count = 0
    for seq_record in SeqIO.parse(input_ltr_fas, "fasta"):

        seqid = seq_record.id
        te_count = te_count + 1
        print('the number of analysis is ' + str(te_count))


        ##consider the no mite te that needs
        ##if the seq id is the LTR only transcript three reading frame
        if 'LTR' in seqid:
            ##function translate to get the pro seq for each te
            pro_seq_six_frame_dic = translate_ltr(seq_record.seq)
            ##function iden_domain_six_rdframe to identify domain for each reading frame
            te_domain_six_rdframe_dic = iden_domain_six_rdframe(pro_seq_six_frame_dic, seqid,hmmscan_exe,working_dir)

            if len(te_domain_six_rdframe_dic) != 0:
                ##function compare_rdframe will compare_rd_frame to get the ideal te_nm
                max_te_nm,domain_num_dic = compare_rdframe(te_domain_six_rdframe_dic,working_dir)
                ##function iden_max_seq could give us the target analyzed sequence.
                seq_out = iden_max_seq(max_te_nm)
                ##function collect_domain will collect information from the hmm opt
                domain_location_dic = collect_domain (working_dir + '/' + seq_out)
                ##function store_domain_info will store the domain infor in a dic
                domain_all_dic = store_domain_info(domain_location_dic,working_dir)


                ##updation3.12: if domain_all_dic is empty will not consider, because this step has filtered out low evalue
                if len(domain_all_dic) != 0:

                    ##function get_domain_len will get the domain len infor for the next step
                    domain_len_dic = get_domain_len (input_domain_length)
                    ##function filter_te will filter the domain infr and give us a domain string pro and pat
                    domain_pro_str,domain_pat_str = filter_te(domain_all_dic,domain_len_dic,working_dir)

                    tar_line = seqid + '\t' + domain_pat_str + '\t' + domain_pro_str
                    te_domain_pattern_dic[tar_line] = 1
                    print('the target line is ' + str(tar_line))

                    ##generate the line to store the number information
                    ##only generate that but do not write out this
                    domain_num_line = str(domain_num_dic[seqid]['RT_ltr']) + '\t' + str(domain_num_dic[seqid]['IN_ltr']) + '\t' + \
                                      str(domain_num_dic[seqid]['RH_ltr']) + '\t' + str(domain_num_dic[seqid]['GAG_ltr']) + '\t' + \
                                      str(domain_num_dic[seqid]['AP_ltr']) + '\t' + str(domain_num_dic[seqid]['RT_line']) + '\t' + \
                                      str(domain_num_dic[seqid]['EN_line']) + '\t' + str(domain_num_dic[seqid]['TR_nomt'])

                    domain_num_line_dic[seqid] = domain_num_line


        ##if the seq id is the LINE and nomite_DNA, it will use the six reading frame
        if 'LINE' in seqid or 'DNA_nMITE' in seqid:
            ##function translate to get the pro seq for each te
            pro_seq_six_frame_dic = translate_no_mite(seq_record.seq)
            ##function iden_domain_six_rdframe to identify domain for each reading frame
            te_domain_six_rdframe_dic = iden_domain_six_rdframe(pro_seq_six_frame_dic, seqid, hmmscan_exe, working_dir)

            if len(te_domain_six_rdframe_dic) != 0:
                ##function compare_rdframe will compare_rd_frame to get the ideal te_nm
                max_te_nm, domain_num_dic = compare_rdframe(te_domain_six_rdframe_dic, working_dir)
                ##function iden_max_seq could give us the target analyzed sequence.
                seq_out = iden_max_seq(max_te_nm)
                ##function collect_domain will collect information from the hmm opt
                domain_location_dic = collect_domain(working_dir + '/' + seq_out)
                ##function store_domain_info will store the domain infor in a dic
                domain_all_dic = store_domain_info(domain_location_dic, working_dir)

                ##updation3.12: if domain_all_dic is empty will not consider, because this step has filtered out low evalue
                if len(domain_all_dic) != 0:


                    ##function get_domain_len will get the domain len infor for the next step
                    domain_len_dic = get_domain_len(input_domain_length)
                    ##function filter_te will filter the domain infr and give us a domain string pro and pat
                    domain_pro_str, domain_pat_str = filter_te(domain_all_dic, domain_len_dic, working_dir)

                    tar_line = seqid + '\t' + domain_pat_str + '\t' + domain_pro_str
                    te_domain_pattern_dic[tar_line] = 1
                    print('the target line is ' + str(tar_line))

                    ##generate the line to store the number information
                    ##only generate that but do not write out this
                    domain_num_line = str(domain_num_dic[seqid]['RT_ltr']) + '\t' + str(domain_num_dic[seqid]['IN_ltr']) + '\t' + \
                                      str(domain_num_dic[seqid]['RH_ltr']) + '\t' + str(domain_num_dic[seqid]['GAG_ltr']) + '\t' + \
                                      str(domain_num_dic[seqid]['AP_ltr']) + '\t' + str(domain_num_dic[seqid]['RT_line']) + '\t' + \
                                      str(domain_num_dic[seqid]['EN_line']) + '\t' + str(domain_num_dic[seqid]['TR_nomt'])

                    domain_num_line_dic[seqid] = domain_num_line

        ##updation 3.12: remove all the temp file generated in this situation
        ##include seven temp files and all the seq out files
        ##removed temp:
        ##opt_temp_domain_order_sort.csv
        ##opt_temp_domain_order.tb
        ##opt_temp_seq_domain_loc_sort.csv
        ##opt_temp_seq_domain_loc.txt
        ##opt_temp_six_rdframe_domain_hmm.tb
        ##*.out
        ##*.seq
        cmd = 'rm ' + working_dir + '/*.out ' +  working_dir + '/*.seq '
        subprocess.call(cmd, shell=True)
        print(cmd)

        ##check the file exit or not
        remove_temp_file(working_dir, 'opt_temp_domain_order_sort.csv')
        remove_temp_file(working_dir, 'opt_temp_domain_order.tb')
        remove_temp_file(working_dir, 'opt_temp_seq_domain_loc_sort.csv')
        remove_temp_file(working_dir, 'opt_temp_seq_domain_loc.txt')
        remove_temp_file(working_dir, 'opt_temp_six_rdframe_domain_hmm.tb')


    return (te_domain_pattern_dic,domain_num_line_dic)


############################################################################
##Step 5: filter out the patterns that are not belong to the specific family
############################################################################
#def filter_pattern (te_domain_pattern_dic):

    ##initiate a dic to store filtered domain pattern
#    ft_dm_pat_line_dic = {}

#    for eachline in te_domain_pattern_dic:
#        eachline = eachline.strip('\n')
#        col = eachline.strip().split()
#        te_nm = col[0]
#        pat = col[1]
#        pro = col[2]

#        if '//LTR' in te_nm:






