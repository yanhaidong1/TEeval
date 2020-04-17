#!/usr/bin/env python

##update: EnSpm match includes En-Spm 2019.1.29

##The user will provide a homemade libraries that conform the formats of the reference TE name rule in this software
##This script will re-name TE name and generate a fasta file containing no-mite DNA to further detect if the DNA TEs is mite

import subprocess
from Bio import SeqIO
import re
##import the other script
from scripts import combine_lib_annot as mt_annot


##define a function to rename mite that will be used in the lib_combination function
def re_mite (mite_nm,te_nm_dic,te_id):
    new_nm = 'ClassII//DNA_MITE/'+mite_nm+'_TE'
    if new_nm in te_nm_dic.keys():
        te_nm_dic[new_nm] += 1
    else:
        te_nm_dic[new_nm] = 1
    final_nm = new_nm + str(te_nm_dic[new_nm]) + '\t' + te_id
    return(final_nm)

##define a function to rename non-mite that will be used in the lib_combination function
def re_no_mite (non_mite_nm,te_nm_dic,te_id):
    new_nm = 'ClassII//DNA_nMITE/'+non_mite_nm+'_TE'
    if new_nm in te_nm_dic.keys():
        te_nm_dic[new_nm] += 1
    else:
        te_nm_dic[new_nm] = 1
    final_nm = new_nm + str(te_nm_dic[new_nm]) + '\t' + te_id
    return(final_nm)

##define a function to rename LTR that will be used in the lib_combination function
def re_ltr (ltr_nm,te_nm_dic,te_id):
    new_nm = 'ClassI//LTR/' + ltr_nm + '_TE'
    if new_nm in te_nm_dic.keys():
        te_nm_dic[new_nm] += 1
    else:
        te_nm_dic[new_nm] = 1
    final_nm = new_nm + str(te_nm_dic[new_nm]) + '\t' + te_id
    return(final_nm)

##define a function to rename LINE and SINE that will be used in the lib_combination function
def re_no_ltr (non_ltr_nm,te_nm_dic,te_id):
    new_nm = 'ClassI//nLTR/' + non_ltr_nm + '_TE'
    if new_nm in te_nm_dic.keys():
        te_nm_dic[new_nm] += 1
    else:
        te_nm_dic[new_nm] = 1
    final_nm = new_nm + str(te_nm_dic[new_nm]) + '\t' + te_id
    return(final_nm)

##define a function to rename LINE and SINE that will be used in the lib_combination function
def te_rc (non_ltr_nm,te_nm_dic,te_id):
    new_nm = 'ClassIII//RC/' + non_ltr_nm + '_TE'
    if new_nm in te_nm_dic.keys():
        te_nm_dic[new_nm] += 1
    else:
        te_nm_dic[new_nm] = 1
    final_nm = new_nm + str(te_nm_dic[new_nm]) + '\t' + te_id
    return(final_nm)

##define a function combine all the libraries
def lib_combination (te_lib_file,hmmscan_exe,working_dir):

    ##initiate a dic to store the modified te name and te seq
    te_seq_dic = {}
    te_no_match_dic = {}
    te_nm_dic = {}

    ##initiate a dic to store the te with len over 800 and dm number is 0
    te_unclear_mite_dic = {}

    ##initiate a dic to store the old name and the new name
    te_new_old_name_dic = {}


    count = 0
    for seq_record in SeqIO.parse(te_lib_file, "fasta"):

        if 'X' not in seq_record.seq:


            final_nm = ''

            te_id = seq_record.id

            count += 1
            print('the analyzed lib te is ' + str(count))

            ###########
            ##match DNA
            if re.match('.*DNA.*', te_id, flags=re.I):

                ############
                ##match mite
                if re.match('.*MITE.*', te_id, flags=re.I):
                    if re.match('.*MuDR.*', te_id, flags=re.I):
                        final_nm = re_mite('MuDR', te_nm_dic,te_id)
                        te_seq_dic[final_nm] = seq_record.seq
                    if re.match('.*CACTA.*', te_id, flags=re.I) or re.match('.*EnSpm.*', te_id, flags=re.I) or re.match('.*En-Spm.*', te_id, flags=re.I):
                        final_nm = re_mite('CACTA', te_nm_dic,te_id)
                        te_seq_dic[final_nm] = seq_record.seq
                    if re.match('.*Tc.*', te_id):
                        final_nm = re_mite('Tc', te_nm_dic,te_id)
                        te_seq_dic[final_nm] = seq_record.seq
                    if re.match('.*Harbinger.*', te_id, flags=re.I):
                        final_nm = re_mite('Harbinger', te_nm_dic,te_id)
                        te_seq_dic[final_nm] = seq_record.seq
                    if re.match('.*hAT.*', te_id, flags=re.I):
                        final_nm = re_mite('hAT', te_nm_dic,te_id)
                        te_seq_dic[final_nm] = seq_record.seq
                    if re.match('.*Mutator.*', te_id, flags=re.I):
                        final_nm = re_mite('Mutator', te_nm_dic,te_id)
                        te_seq_dic[final_nm] = seq_record.seq
                    if re.match('.*Novosib.*', te_id, flags=re.I):
                        final_nm = re_mite('Novosib', te_nm_dic,te_id)
                        te_seq_dic[final_nm] = seq_record.seq
                    ##match mite but not the followings
                    if not re.match('.*CACTA.*', te_id, flags=re.I) and not re.match('.*Tc.*', te_id) \
                        and not re.match('.*Harbinger.*', te_id, flags=re.I) and not re.match('.*hAT.*', te_id, flags=re.I) \
                        and not re.match('.*Mutator.*', te_id, flags=re.I) and not re.match('.*Novosib.*', te_id, flags=re.I) \
                        and not re.match('.*EnSpm.*', te_id, flags=re.I) and not re.match('.*MuDR.*', te_id, flags=re.I) \
                        and not re.match('.*En-Spm.*', te_id, flags=re.I):
                        final_nm = re_mite('mite_Others', te_nm_dic,te_id)
                        te_seq_dic[final_nm] = seq_record.seq


                ##define a function to detect classify mite te
                ##the rule is length is below 800, and the number of dna_tr is 0
                ##if not rematch mite and this part of DNA will be stored in the no mite dna dic files to further annotate these DNAs
                if not re.match('.*MITE.*', te_id, flags=re.I):
                    ##translate te
                    pro_seq_dic = mt_annot.translate(seq_record.seq)
                    ##identify domain of six rdframe
                    te_domain_six_rdframe_dic = mt_annot.iden_domain_six_rdframe(pro_seq_dic, te_id,hmmscan_exe,working_dir)

                    ##if there is at one domain was found
                    if len(te_domain_six_rdframe_dic) != 0:
                        ##get the best te and its dic
                        domain_num_dic = mt_annot.compare_rdframe(te_domain_six_rdframe_dic,working_dir)
                        te_dm_num = domain_num_dic[te_id]['TR']
                        ##filter the te
                        ##if te_dm is equal to 0 and len of seq is below 800 bp, the te is regarded as mite
                        if int(te_dm_num) == 0:

                            ##this step will filter out seq len is over 800 bp
                            if len(seq_record.seq) <= 800:
                                #print('the te_seq is ' + str(len(seq_record.seq)))
                                ##this te will be regarded as the te mite
                                if re.match('.*CACTA.*', te_id, flags=re.I) or re.match('.*EnSpm.*', te_id, flags=re.I) or re.match('.*En-Spm.*', te_id, flags=re.I):
                                    final_nm = re_mite('CACTA', te_nm_dic,te_id)
                                    te_seq_dic[final_nm] = seq_record.seq
                                if re.match('.*MuDR.*', te_id, flags=re.I):
                                    final_nm = re_mite('MuDR', te_nm_dic,te_id)
                                    te_seq_dic[final_nm] = seq_record.seq
                                if re.match('.*Tc.*', te_id):
                                    final_nm = re_mite('Tc', te_nm_dic,te_id)
                                    te_seq_dic[final_nm] = seq_record.seq
                                if re.match('.*Harbinger.*', te_id, flags=re.I):
                                    final_nm = re_mite('Harbinger', te_nm_dic,te_id)
                                    te_seq_dic[final_nm] = seq_record.seq
                                if re.match('.*hAT.*', te_id, flags=re.I):
                                    final_nm = re_mite('hAT', te_nm_dic,te_id)
                                    te_seq_dic[final_nm] = seq_record.seq
                                if re.match('.*Mutator.*', te_id, flags=re.I):
                                    final_nm = re_mite('Mutator', te_nm_dic,te_id)
                                    te_seq_dic[final_nm] = seq_record.seq
                                if re.match('.*Novosib.*', te_id, flags=re.I):
                                    final_nm = re_mite('Novosib', te_nm_dic,te_id)
                                    te_seq_dic[final_nm] = seq_record.seq

                                ##match mite but not the followings
                                if not re.match('.*CACTA.*', te_id, flags=re.I) and not re.match('.*Tc.*', te_id,flags=re.I) \
                                        and not re.match('.*Harbinger.*', te_id, flags=re.I) and not re.match('.*hAT.*',te_id,flags=re.I) \
                                        and not re.match('.*Mutator.*', te_id, flags=re.I) and not re.match('.*Novosib.*',te_id,flags=re.I) \
                                        and not re.match('.*EnSpm.*', te_id, flags=re.I) and not re.match('.*MuDR.*', te_id,flags=re.I) \
                                        and not re.match('.*En-Spm.*', te_id, flags=re.I):
                                    final_nm = re_mite('mite_Others', te_nm_dic,te_id)
                                    te_seq_dic[final_nm] = seq_record.seq

                            else:
                                te_unclear_mite_dic[te_id] = seq_record.seq

                        ##if the te has domain, so it will be regarded as no mite te
                        else:
                            if re.match('.*EnSpm.*', te_id, flags=re.I) or re.match('.*CACTA.*', te_id, flags=re.I) or re.match('.*En-Spm.*', te_id, flags=re.I):
                                final_nm = re_no_mite('EnSpm', te_nm_dic,te_id)
                                te_seq_dic[final_nm] = seq_record.seq
                            if re.match('.*MuDR.*', te_id, flags=re.I):
                                final_nm = re_no_mite('MuDR', te_nm_dic,te_id)
                                te_seq_dic[final_nm] = seq_record.seq
                            if re.match('.*Tc.*', te_id):
                                final_nm = re_no_mite('Tc', te_nm_dic,te_id)
                                te_seq_dic[final_nm] = seq_record.seq
                            if re.match('.*Harbinger.*', te_id, flags=re.I):
                                final_nm = re_no_mite('Harbinger', te_nm_dic,te_id)
                                te_seq_dic[final_nm] = seq_record.seq
                            if re.match('.*hAT.*', te_id, flags=re.I):
                                final_nm = re_no_mite('hAT', te_nm_dic,te_id)
                                te_seq_dic[final_nm] = seq_record.seq
                            if re.match('.*Mutator.*', te_id, flags=re.I):
                                final_nm = re_no_mite('Mutator', te_nm_dic,te_id)
                                te_seq_dic[final_nm] = seq_record.seq
                            if re.match('.*Novosib.*', te_id, flags=re.I):
                                final_nm = re_no_mite('Novosib', te_nm_dic,te_id)
                                te_seq_dic[final_nm] = seq_record.seq

                            ###########
                            ##other DNA
                            if not re.match('.*CACTA.*', te_id, flags=re.I) and not re.match('.*Tc.*', te_id) \
                                    and not re.match('.*Harbinger.*', te_id, flags=re.I) and not re.match('.*hAT.*', te_id,flags=re.I) \
                                    and not re.match('.*Mutator.*', te_id, flags=re.I) and not re.match('.*Novosib.*',te_id, flags=re.I) \
                                    and not re.match('.*EnSpm.*', te_id, flags=re.I) and not re.match('.*MuDR.*', te_id,flags=re.I) \
                                    and not re.match('.*En-Spm.*', te_id, flags=re.I):

                                new_nm = 'ClassII//DNA_nMITE/dna_Others_TE'
                                if new_nm in te_nm_dic.keys():
                                    te_nm_dic[new_nm] += 1
                                else:
                                    te_nm_dic[new_nm] = 1
                                final_nm = new_nm + str(te_nm_dic[new_nm]) + '\t' + te_id
                                te_seq_dic[final_nm] = seq_record.seq

                    ##if there is no domain found
                    else:
                        ##this step will filter out seq len is over 800 bp
                        if len(seq_record.seq) <= 800:
                            #print('the te_seq is ' + str(len(seq_record.seq)))
                            ##this te will be regarded as the te mite
                            if re.match('.*CACTA.*', te_id, flags=re.I) or re.match('.*EnSpm.*', te_id, flags=re.I) or re.match('.*En-Spm.*', te_id, flags=re.I):
                                final_nm = re_mite('CACTA', te_nm_dic,te_id)
                                te_seq_dic[final_nm] = seq_record.seq
                            if re.match('.*MuDR.*', te_id, flags=re.I):
                                final_nm = re_mite('MuDR', te_nm_dic,te_id)
                                te_seq_dic[final_nm] = seq_record.seq
                            if re.match('.*Tc.*', te_id):
                                final_nm = re_mite('Tc', te_nm_dic,te_id)
                                te_seq_dic[final_nm] = seq_record.seq
                            if re.match('.*Harbinger.*', te_id, flags=re.I):
                                final_nm = re_mite('Harbinger', te_nm_dic,te_id)
                                te_seq_dic[final_nm] = seq_record.seq
                            if re.match('.*hAT.*', te_id, flags=re.I):
                                final_nm = re_mite('hAT', te_nm_dic,te_id)
                                te_seq_dic[final_nm] = seq_record.seq
                            if re.match('.*Mutator.*', te_id, flags=re.I):
                                final_nm = re_mite('Mutator', te_nm_dic,te_id)
                                te_seq_dic[final_nm] = seq_record.seq
                            if re.match('.*Novosib.*', te_id, flags=re.I):
                                final_nm = re_mite('Novosib', te_nm_dic,te_id)
                                te_seq_dic[final_nm] = seq_record.seq

                            ##match mite but not the followings
                            if not re.match('.*CACTA.*', te_id, flags=re.I) and not re.match('.*Tc.*', te_id) \
                                    and not re.match('.*Harbinger.*', te_id, flags=re.I) and not re.match('.*hAT.*', te_id,flags=re.I) \
                                    and not re.match('.*Mutator.*', te_id, flags=re.I) and not re.match('.*Novosib.*',te_id, flags=re.I) \
                                    and not re.match('.*EnSpm.*', te_id, flags=re.I) and not re.match('.*MuDR.*', te_id,flags=re.I) \
                                    and not re.match('.*En-Spm.*', te_id, flags=re.I):
                                final_nm = re_mite('mite_Others', te_nm_dic,te_id)
                                te_seq_dic[final_nm] = seq_record.seq

                        else:
                            te_unclear_mite_dic[te_id] = seq_record.seq



            ##################
            ##do not match DNA
            ##################
            if not re.match('.*DNA.*', te_id, flags=re.I):

                #################################################################
                ##sometimes, the library does not match DNA, but they are DNA TEs

                ############
                ##match MITE
                if re.match('.*MITE.*', te_id, flags=re.I):
                    if re.match('.*CACTA.*', te_id, flags=re.I) or re.match('.*EnSpm.*', te_id, flags=re.I) or re.match('.*En-Spm.*', te_id, flags=re.I):
                        final_nm = re_mite('CACTA', te_nm_dic,te_id)
                        te_seq_dic[final_nm] = seq_record.seq
                    if re.match('.*MuDR.*', te_id, flags=re.I):
                        final_nm = re_mite('MuDR', te_nm_dic,te_id)
                        te_seq_dic[final_nm] = seq_record.seq
                    if re.match('.*Tc.*', te_id):
                        final_nm = re_mite('Tc', te_nm_dic,te_id)
                        te_seq_dic[final_nm] = seq_record.seq
                    if re.match('.*Harbinger.*', te_id, flags=re.I):
                        final_nm = re_mite('Harbinger', te_nm_dic,te_id)
                        te_seq_dic[final_nm] = seq_record.seq
                    if re.match('.*hAT.*', te_id, flags=re.I):
                        final_nm = re_mite('hAT', te_nm_dic,te_id)
                        te_seq_dic[final_nm] = seq_record.seq
                    if re.match('.*Mutator.*', te_id, flags=re.I):
                        final_nm = re_mite('Mutator', te_nm_dic,te_id)
                        te_seq_dic[final_nm] = seq_record.seq
                    if re.match('.*Novosib.*', te_id, flags=re.I):
                        final_nm = re_mite('Novosib', te_nm_dic,te_id)
                        te_seq_dic[final_nm] = seq_record.seq
                    ##match mite but not the followings
                    if not re.match('.*CACTA.*', te_id, flags=re.I) and not re.match('.*Tc.*', te_id) \
                            and not re.match('.*Harbinger.*', te_id, flags=re.I) and not re.match('.*hAT.*', te_id,flags=re.I) \
                            and not re.match('.*Mutator.*', te_id, flags=re.I) and not re.match('.*Novosib.*', te_id,flags=re.I) \
                            and not re.match('.*EnSpm.*', te_id, flags=re.I) and not re.match('.*MuDR.*', te_id,flags=re.I) \
                            and not re.match('.*En-Spm.*', te_id, flags=re.I):
                        final_nm = re_mite('mite_Others', te_nm_dic,te_id)
                        te_seq_dic[final_nm] = seq_record.seq


                ################
                ##not match MITE
                if not re.match('.*MITE.*', te_id, flags=re.I):

                    ##if re.match the followings, tes will be analyzed using hmm
                    if re.match('.*CACTA.*', te_id, flags=re.I) or re.match('.*EnSpm.*', te_id, flags=re.I) \
                        or re.match('.*Tc.*', te_id) or re.match('.*Harbinger.*', te_id, flags=re.I) \
                        or re.match('.*hAT.*', te_id, flags=re.I) or re.match('.*Mutator.*', te_id, flags=re.I) \
                        or re.match('.*Novosib.*', te_id, flags=re.I) or re.match('.*MuDR.*', te_id,flags=re.I) \
                        or re.match('.*En-Spm.*', te_id, flags=re.I):

                        ##translate te
                        pro_seq_dic = mt_annot.translate(seq_record.seq)
                        ##identify domain of six rdframe
                        te_domain_six_rdframe_dic = mt_annot.iden_domain_six_rdframe(pro_seq_dic, te_id)

                        ##if there are at least one domains
                        if len(te_domain_six_rdframe_dic) != 0:
                            ##get the best te and its dic
                            domain_num_dic = mt_annot.compare_rdframe(te_domain_six_rdframe_dic)
                            te_dm_num = domain_num_dic[te_id]['TR']
                            ##filter the te
                            ##if te_dm is equal to 0 and len of seq is below 800 bp, the te is regarded as mite
                            if int(te_dm_num) == 0:
                                ##this step will filter out seq len is over 800 bp
                                if len(seq_record.seq) <= 800:
                                    #print('the te_seq is ' + str(len(seq_record.seq)))
                                    ##this te will be regarded as the te mite
                                    if re.match('.*CACTA.*', te_id, flags=re.I) or re.match('.*EnSpm.*', te_id, flags=re.I) or re.match('.*En-Spm.*', te_id, flags=re.I):
                                        final_nm = re_mite('CACTA', te_nm_dic,te_id)
                                        te_seq_dic[final_nm] = seq_record.seq
                                    if re.match('.*MuDR.*', te_id, flags=re.I):
                                        final_nm = re_mite('MuDR', te_nm_dic,te_id)
                                        te_seq_dic[final_nm] = seq_record.seq
                                    if re.match('.*Tc.*', te_id):
                                        final_nm = re_mite('Tc', te_nm_dic,te_id)
                                        te_seq_dic[final_nm] = seq_record.seq
                                    if re.match('.*Harbinger.*', te_id, flags=re.I):
                                        final_nm = re_mite('Harbinger', te_nm_dic,te_id)
                                        te_seq_dic[final_nm] = seq_record.seq
                                    if re.match('.*hAT.*', te_id, flags=re.I):
                                        final_nm = re_mite('hAT', te_nm_dic,te_id)
                                        te_seq_dic[final_nm] = seq_record.seq
                                    if re.match('.*Mutator.*', te_id, flags=re.I):
                                        final_nm = re_mite('Mutator', te_nm_dic,te_id)
                                        te_seq_dic[final_nm] = seq_record.seq
                                    if re.match('.*Novosib.*', te_id, flags=re.I):
                                        final_nm = re_mite('Novosib', te_nm_dic,te_id)
                                        te_seq_dic[final_nm] = seq_record.seq

                                    ##match mite but not the followings
                                    if not re.match('.*CACTA.*', te_id, flags=re.I) and not re.match('.*Tc.*', te_id,flags=re.I) \
                                            and not re.match('.*Harbinger.*', te_id, flags=re.I) and not re.match('.*hAT.*',te_id,flags=re.I) \
                                            and not re.match('.*Mutator.*', te_id, flags=re.I) and not re.match('.*Novosib.*',te_id,flags=re.I) \
                                            and not re.match('.*EnSpm.*', te_id, flags=re.I) and not re.match('.*MuDR.*', te_id,flags=re.I) \
                                            and not re.match('.*En-Spm.*', te_id, flags=re.I):
                                        final_nm = re_mite('mite_Others', te_nm_dic,te_id)
                                        te_seq_dic[final_nm] = seq_record.seq

                                else:
                                    te_unclear_mite_dic[te_id] = seq_record.seq

                            ##if the te has domain, so it will be regarded as no mite te
                            else:
                                if re.match('.*EnSpm.*', te_id, flags=re.I) or re.match('.*CACTA.*', te_id, flags=re.I) or re.match('.*En-Spm.*', te_id, flags=re.I):
                                    final_nm = re_no_mite('EnSpm', te_nm_dic,te_id)
                                    te_seq_dic[final_nm] = seq_record.seq
                                if re.match('.*MuDR.*', te_id, flags=re.I):
                                    final_nm = re_no_mite('MuDR', te_nm_dic,te_id)
                                    te_seq_dic[final_nm] = seq_record.seq
                                if re.match('.*Tc.*', te_id):
                                    final_nm = re_no_mite('Tc', te_nm_dic,te_id)
                                    te_seq_dic[final_nm] = seq_record.seq
                                if re.match('.*Harbinger.*', te_id, flags=re.I):
                                    final_nm = re_no_mite('Harbinger', te_nm_dic,te_id)
                                    te_seq_dic[final_nm] = seq_record.seq
                                if re.match('.*hAT.*', te_id, flags=re.I):
                                    final_nm = re_no_mite('hAT', te_nm_dic,te_id)
                                    te_seq_dic[final_nm] = seq_record.seq
                                if re.match('.*Mutator.*', te_id, flags=re.I):
                                    final_nm = re_no_mite('Mutator', te_nm_dic,te_id)
                                    te_seq_dic[final_nm] = seq_record.seq
                                if re.match('.*Novosib.*', te_id, flags=re.I):
                                    final_nm = re_no_mite('Novosib', te_nm_dic,te_id)
                                    te_seq_dic[final_nm] = seq_record.seq

                        ##if the domain is zero
                        else:

                            ##this step will filter out seq len is over 800 bp
                            if len(seq_record.seq) <= 800:
                                #print('the te_seq is ' + str(len(seq_record.seq)))
                                ##this te will be regarded as the te mite
                                if re.match('.*CACTA.*', te_id, flags=re.I) or re.match('.*EnSpm.*', te_id, flags=re.I) or re.match('.*En-Spm.*', te_id, flags=re.I):
                                    final_nm = re_mite('CACTA', te_nm_dic,te_id)
                                    te_seq_dic[final_nm] = seq_record.seq
                                if re.match('.*MuDR.*', te_id, flags=re.I):
                                    final_nm = re_mite('MuDR', te_nm_dic,te_id)
                                    te_seq_dic[final_nm] = seq_record.seq
                                if re.match('.*Tc.*', te_id):
                                    final_nm = re_mite('Tc', te_nm_dic,te_id)
                                    te_seq_dic[final_nm] = seq_record.seq
                                if re.match('.*Harbinger.*', te_id, flags=re.I):
                                    final_nm = re_mite('Harbinger', te_nm_dic,te_id)
                                    te_seq_dic[final_nm] = seq_record.seq
                                if re.match('.*hAT.*', te_id, flags=re.I):
                                    final_nm = re_mite('hAT', te_nm_dic,te_id)
                                    te_seq_dic[final_nm] = seq_record.seq
                                if re.match('.*Mutator.*', te_id, flags=re.I):
                                    final_nm = re_mite('Mutator', te_nm_dic,te_id)
                                    te_seq_dic[final_nm] = seq_record.seq
                                if re.match('.*Novosib.*', te_id, flags=re.I):
                                    final_nm = re_mite('Novosib', te_nm_dic,te_id)
                                    te_seq_dic[final_nm] = seq_record.seq

                                ##match mite but not the followings
                                if not re.match('.*CACTA.*', te_id, flags=re.I) and not re.match('.*Tc.*', te_id,flags=re.I) \
                                        and not re.match('.*Harbinger.*', te_id, flags=re.I) and not re.match('.*hAT.*',te_id,flags=re.I) \
                                        and not re.match('.*Mutator.*', te_id, flags=re.I) and not re.match('.*Novosib.*',te_id,flags=re.I) \
                                        and not re.match('.*EnSpm.*', te_id, flags=re.I) and not re.match('.*MuDR.*', te_id,flags=re.I) \
                                        and not re.match('.*En-Spm.*', te_id, flags=re.I):
                                    final_nm = re_mite('mite_Others', te_nm_dic,te_id)
                                    te_seq_dic[final_nm] = seq_record.seq

                            else:
                                te_unclear_mite_dic[te_id] = seq_record.seq


                    ###########
                    ##match LTR
                    ##some LTR has LINE name, and we will filter out this situation
                    if re.match('.*LTR.*', te_id, flags=re.I) and not re.match('.*SINE.*', te_id, flags=re.I) \
                            and not re.match('.*LINE.*', te_id, flags=re.I):
                        if re.match('.*Gypsy.*', te_id, flags=re.I):
                            final_nm = re_ltr('Gypsy', te_nm_dic,te_id)
                            te_seq_dic[final_nm] = seq_record.seq
                        if re.match('.*Copia.*', te_id, flags=re.I):
                            final_nm = re_ltr('Copia', te_nm_dic,te_id)
                            te_seq_dic[final_nm] = seq_record.seq
                        if re.match('.*Caulimovirus.*', te_id, flags=re.I):
                            final_nm = re_ltr('Caulimovirus', te_nm_dic,te_id)
                            te_seq_dic[final_nm] = seq_record.seq
                        ##other LTRs
                        if not re.match('.*Gypsy.*', te_id, flags=re.I) and not re.match('.*Copia.*', te_id, flags=re.I) \
                                and not re.match('.*Caulimovirus.*', te_id, flags=re.I):

                            final_nm = re_ltr('ltr_Others', te_nm_dic,te_id)
                            te_seq_dic[final_nm] = seq_record.seq

                    ############
                    ##match SINE
                    if re.match('.*SINE.*', te_id, flags=re.I) and not re.match('.*LTR.*', te_id, flags=re.I):
                        final_nm = re_no_ltr ('SINE', te_nm_dic,te_id)
                        te_seq_dic[final_nm] = seq_record.seq

                    ############
                    ##match LINE
                    if re.match('.*LINE.*', te_id, flags=re.I) and not re.match('.*LTR.*', te_id, flags=re.I):
                        final_nm = re_no_ltr('LINE', te_nm_dic,te_id)
                        te_seq_dic[final_nm] = seq_record.seq

                    ##########
                    ##match L1
                    if re.match('.*L1.*', te_id, flags=re.I) and not re.match('.*LINE.*', te_id, flags=re.I) \
                            and not re.match('.*SINE.*', te_id, flags=re.I) and not re.match('.*LTR.*', te_id, flags=re.I):
                        final_nm = re_no_ltr('LINE', te_nm_dic,te_id)
                        te_seq_dic[final_nm] = seq_record.seq

                    ##########
                    ##match RC  for the plant only use helitron for the mam: maverick
                    if re.match('.*Helitron.*', te_id, flags=re.I) and not re.match('.*LTR.*', te_id, flags=re.I):
                        final_nm = te_rc('Helitron', te_nm_dic,te_id)
                        te_seq_dic[final_nm] = seq_record.seq


                    if not re.match('.*LTR.*', te_id, flags=re.I) and not re.match('.*LINE.*', te_id, flags=re.I) \
                            and not re.match('.*SINE.*', te_id, flags=re.I) and not re.match('.*MITE.*', te_id, flags=re.I) \
                            and not re.match('.*L1.*', te_id, flags=re.I) and not re.match('.*Helitron.*', te_id, flags=re.I) \
                            and not re.match('.*CACTA.*', te_id, flags=re.I) and not re.match('.*Tc.*', te_id) \
                            and not re.match('.*Harbinger.*', te_id, flags=re.I) and not re.match('.*hAT.*', te_id, flags=re.I) \
                            and not re.match('.*Mutator.*', te_id, flags=re.I) and not re.match('.*Novosib.*', te_id, flags=re.I) \
                            and not re.match('.*EnSpm.*', te_id, flags=re.I) and not re.match('.*MuDR.*', te_id, flags=re.I) \
                            and not re.match('.*En-Spm.*', te_id, flags=re.I):

                        te_no_match_dic[te_id] = seq_record.seq

            te_new_old_name_dic[final_nm] = te_id


    return(te_seq_dic,te_no_match_dic,te_unclear_mite_dic,te_new_old_name_dic)



##define a function to divide the te into mite and no mite
def mite_classify (te_seq_dic):
    te_mite_seq_dic = {}
    te_no_mite_seq_dic = {}
    for eachte in te_seq_dic:
        if 'DNA_MITE' in eachte:
            te_mite_seq_dic[eachte] = te_seq_dic[eachte]
        else:
            te_no_mite_seq_dic[eachte] = te_seq_dic[eachte]

    return (te_mite_seq_dic,te_no_mite_seq_dic)

##define a function to write out the mite te
def write_to_mite_fam_dir (t_fam_nm,te_mite_seq_dic,te_family_dir):
    with open(te_family_dir + '/MITE_' + t_fam_nm + '.lib', 'w+') as opt_te_lib:
        for eachte in te_mite_seq_dic:
            mt = re.match('.+\/\/.+\/(.+)_TE.+', eachte)
            fam_nm = mt.group(1)
            if t_fam_nm == fam_nm:
                opt_te_lib.write('>' + eachte + '\n' + str(te_mite_seq_dic[eachte]) + '\n')

##define a function to write out the te lib to the family directory
def write_to_no_mite_fam_dir (t_fam_nm,te_no_mite_seq_dic,te_family_dir):
    with open (te_family_dir + '/' + t_fam_nm + '.lib','w+') as opt_te_lib:
        for eachte in te_no_mite_seq_dic:
            mt = re.match('.+\/\/.+\/(.+)_TE.+', eachte)
            fam_nm = mt.group(1)
            if t_fam_nm == fam_nm:
                opt_te_lib.write('>' + eachte + '\n' + str(te_no_mite_seq_dic[eachte]) + '\n')

##define a function to classify each type of te and store each of them to a directory
def classify_tes (te_mite_seq_dic,te_no_mite_seq_dic,te_family_dir):

    ##te_mite
    write_to_mite_fam_dir('CACTA', te_mite_seq_dic, te_family_dir)
    write_to_mite_fam_dir('Tc', te_mite_seq_dic, te_family_dir)
    write_to_mite_fam_dir('Harbinger', te_mite_seq_dic, te_family_dir)
    write_to_mite_fam_dir('hAT', te_mite_seq_dic, te_family_dir)
    write_to_mite_fam_dir('Mutator', te_mite_seq_dic, te_family_dir)
    write_to_mite_fam_dir('Novosib', te_mite_seq_dic, te_family_dir)
    write_to_mite_fam_dir('mite_Others', te_mite_seq_dic, te_family_dir)
    write_to_mite_fam_dir('EnSpm', te_mite_seq_dic, te_family_dir)
    write_to_mite_fam_dir('MuDR', te_mite_seq_dic, te_family_dir)
    ##no te_mite
    write_to_no_mite_fam_dir('Tc', te_no_mite_seq_dic, te_family_dir)
    write_to_no_mite_fam_dir('Harbinger', te_no_mite_seq_dic, te_family_dir)
    write_to_no_mite_fam_dir('hAT', te_no_mite_seq_dic, te_family_dir)
    write_to_no_mite_fam_dir('Mutator', te_no_mite_seq_dic, te_family_dir)
    write_to_no_mite_fam_dir('Novosib', te_no_mite_seq_dic, te_family_dir)
    write_to_no_mite_fam_dir('EnSpm', te_no_mite_seq_dic, te_family_dir)
    write_to_no_mite_fam_dir('MuDR', te_no_mite_seq_dic, te_family_dir)
    write_to_no_mite_fam_dir('dna_Others', te_no_mite_seq_dic, te_family_dir)
    ##LTR
    write_to_no_mite_fam_dir('Gypsy', te_no_mite_seq_dic, te_family_dir)
    write_to_no_mite_fam_dir('Copia', te_no_mite_seq_dic, te_family_dir)
    write_to_no_mite_fam_dir('Caulimovirus', te_no_mite_seq_dic, te_family_dir)
    write_to_no_mite_fam_dir('ltr_Others', te_no_mite_seq_dic, te_family_dir)
    ##nLTR
    write_to_no_mite_fam_dir('SINE', te_no_mite_seq_dic, te_family_dir)
    write_to_no_mite_fam_dir('LINE', te_no_mite_seq_dic, te_family_dir)
    ##RC
    write_to_no_mite_fam_dir('Helitron', te_no_mite_seq_dic, te_family_dir)


