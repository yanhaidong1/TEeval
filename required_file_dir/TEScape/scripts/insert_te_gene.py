#!/usr/bin/env python

##This script will be contained in TE_summary wrapper script
##This script is to generate information about proportion of gene elements covered by TEs.
##This script also generates a table contains information about insertion detailed information for each gene
##The elements include exon, intron, cds, UTR, five_prime_UTR, three_prime_UTR, gene, intergenic regions.
##This script is to use the TE summary from the TEscape output
##The output format is as follows:
##chr mRNAnm genetp start end direct annot tenm te_traits

##The input file 1 is gene gff files
##The input file 2 is the te gff files: it could be gff files generated from the TEscape and also user provided gff file
##make sure the ID in the gff3 has unique name

import re


##define a function to store mRNA information
def store_mRNA_gene_infor (input_gene_gff_file):
    ##initiate a dic to store mRNA information
    mRNA_dic = {}
    gene_dic = {}
    with open(input_gene_gff_file, 'r') as ipt_gene:
        for eachline in ipt_gene:
            if not eachline.startswith('#'):
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                ##generate store information
                chr = col[0]
                type = col[2]
                bg = col[3]
                ed = col[4]
                dir = col[6]
                annot = col[8]

                ##store mRNA information
                if type == 'mRNA':

                    col_annot = annot.split(';')
                    mt = re.match('ID=(.+)', col_annot[0])
                    ID = mt.group(1)
                    #print ('the id is ' + ID)

                    mRNA_dic[ID] = {'chr': chr, 'type': type, 'bg': bg, 'ed': ed, 'dir': dir, 'annot': annot}

                ##store gene information
                if type == 'gene':

                    if ';' in annot:  ##some annotation file do not contain ;
                        col_annot = annot.split(';')
                        mt = re.match('ID=(.+)', col_annot[0])
                        ID = mt.group(1)
                        gene_dic[ID] = {'chr': chr, 'type': type, 'bg': bg, 'ed': ed, 'dir': dir, 'annot': annot}

                    else:
                        mt = re.match('ID=(.+)', annot)
                        ID = mt.group(1)
                        gene_dic[ID] = {'chr': chr, 'type': type, 'bg': bg, 'ed': ed, 'dir': dir, 'annot': annot}

    return(mRNA_dic,gene_dic)


##define a function to store element information
##this function is similar with the previous store information, but will generate intron information
def store_mRNA_elem_adintron_infor (input_gene_gff_file,mRNA_dic):
    ##initiate a dic to store all the information
    all_dic = {}

    for eachmRNA_ID in mRNA_dic:
        ##initiate a dic to store all the elem into a list
        elem_list = []

        ##initiate a dic to store the compared CDS to generate the intron information
        comp_list = []

        ##store all the CDS into a list
        cds_store_list = []

        with open(input_gene_gff_file, 'r') as ipt_gene:
            for eachline in ipt_gene:
                if not eachline.startswith('#'):
                    eachline = eachline.strip('\n')
                    col = eachline.strip().split()
                    ##generate store information
                    chr = col[0]
                    type = col[2]
                    bg = col[3]
                    ed = col[4]
                    dir = col[6]
                    annot = col[8]
                    col_annot = annot.split(';')
                    mt = re.match('ID=(.+)', col_annot[0])
                    ID = mt.group(1)

                    if type != 'mRNA' and type != 'gene':
                        #print('the id is ' + ID)
                        ##parent ID is the transcript id, and it is also used to name the intron
                        mt = re.match('Parent=(.+)', col_annot[1])
                        Parent_ID = mt.group(1)
                        #print ('the parent id is ' + Parent_ID)

                        if eachmRNA_ID in ID:

                            #print ('the matched ID is ' + ID)
                            ##store all the element information in the elem_list
                            ##this elem_list should also contain information of the intron dected in the later part
                            elem_dic = {'chr': chr, 'type': type, 'bg': bg, 'ed': ed, 'dir': dir, 'annot': annot}
                            elem_list.append(elem_dic)

                            ##initiate a dic to store detail information for each cds
                            if type == 'CDS':
                                ##store the CDS information
                                cds_dic = {'chr': chr, 'type': type,'elem_nm':ID, 'bg': bg, 'ed': ed,
                                           'dir': dir,'annot': annot}
                                cds_store_list.append(cds_dic)
                                #print (cds_store_list)

        ##analyze each cds in the cds_store_list
        ##store the first cds_dic into the comp_list
        comp_list.append(cds_store_list[0])
        ##remove the first element of the cds_store_list
        cds_store_list.pop(0)
        ##compare eachdic from cds_store_list to the last element of the comp_list
        ##after each comparing, move the compared dic to the comp_list to be as the last element in the comp_list

        count = 0
        for eachdic in cds_store_list:
            count += 1 ##this is used to name the count number of intron
            ##get the location information from the comp_list
            ##extract the last element of the comp_list
            comp_list_dic = comp_list[-1]
            comp_ed = comp_list_dic['ed']

            ##extract eachdic information from the cds_store_list
            store_bg = eachdic['bg']

            ##get the intron location
            intron_bg = int(comp_ed) + 1
            intron_ed = int(store_bg) - 1
            ##name the intron
            intron_nm = eachmRNA_ID + '.intron' + str(count)
            intron_annot = 'ID=' + intron_nm + ';Parent=' + eachmRNA_ID
            ##store intron in to elem_dic
            elem_dic = {'chr': eachdic['chr'], 'type': 'intron', 'bg': str(intron_bg), 'ed': str(intron_ed), 'dir': eachdic['dir'], 'annot': intron_annot}
            elem_list.append(elem_dic)

            ##after get the intron information, move the analyzed store dic to the comp_list_dic
            ##the next dic comparision will happen between the last one of comp_list_dic and the current one in the cds_store_list
            comp_list.append(eachdic)

        ##store the elem_list to the all_dic
        ##the all_dic use each mRNA ID as key, and elem_list belonging to this mRNA as value
        all_dic[eachmRNA_ID] = elem_list

    return (all_dic)


##define a function to store TE information
def store_te_infor (input_te_gff):
    te_infor_dic = {}
    with open (input_te_gff,'r') as ipt_te:
        for eachline in ipt_te:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            chr = col[0]
            bg = col[1]
            ed = col[2]
            te_nm = col[3]
            te_rest = col[4] + '\t' + col[8] + '\t' + col[11] + '\t' + col[12] + '\t' + col[13] + '\t' + col[14]

            te_infor_dic[te_nm] = {'chr':chr,'bg':bg,'ed':ed, 'rest':te_rest}

    return (te_infor_dic)


##define a function to get detailed TE insertion for each element
def gene_elem_pro (all_dic,mRNA_dic,te_infor_dic,flank_range):
    ##initate a dic to store the final gene insertion line in TE
    TE_gene_insertion_dic = {}

    count = 0
    for eachmRNA_ID in all_dic:
        count += 1
        #print ('the analyzed mRNA count is ' + str(count))

        all_te_dic = {}

        ##get the mRNA position information
        mRNA_bg = mRNA_dic[eachmRNA_ID]['bg']
        mRNA_ed = mRNA_dic[eachmRNA_ID]['ed']
        mRNA_chr = mRNA_dic[eachmRNA_ID]['chr']
        ##get the flanking regions
        mRNA_up = int(mRNA_dic[eachmRNA_ID]['bg']) - flank_range
        mRNA_down = int(mRNA_dic[eachmRNA_ID]['ed']) + flank_range

        for eachte in te_infor_dic:
            ##get the te position information
            te_bg = te_infor_dic[eachte]['bg']
            te_ed = te_infor_dic[eachte]['ed']
            te_chr = te_infor_dic[eachte]['chr']

            ##if te and mRNA are in the same chr
            if mRNA_chr == te_chr:
                ##if te cover mRNA
                if int(te_ed) >= int(mRNA_bg) and int(te_bg) <= int(mRNA_ed):
                    ##store the information in the all_te_dic
                    all_te_dic[eachte] = 1

                ##if te cover upstream of mRNA
                if int(te_ed) <= int(mRNA_bg) and int(te_ed) >= int(mRNA_up):
                    ##if direction is +
                    if mRNA_dic[eachmRNA_ID]['dir'] == '+':
                        final_line = mRNA_chr + '\t' + eachmRNA_ID + '\t' + 'upstream' + '\t' + \
                                     str(mRNA_up) + '\t' + str(mRNA_bg) + '\t' + mRNA_dic[eachmRNA_ID]['dir'] + '\t' + \
                                     mRNA_dic[eachmRNA_ID]['annot'] + '\t' + eachte + '\t' + te_infor_dic[eachte]['bg'] + '\t' + \
                                     te_infor_dic[eachte]['ed'] + '\t' + te_infor_dic[eachte]['rest']
                        TE_gene_insertion_dic[final_line] = 1
                    ##if direction is -
                    else:
                        final_line = mRNA_chr + '\t' + eachmRNA_ID + '\t' + 'downstream' + '\t' + \
                                     str(mRNA_up) + '\t' + str(mRNA_bg) + '\t' + mRNA_dic[eachmRNA_ID]['dir'] + '\t' + \
                                     mRNA_dic[eachmRNA_ID]['annot'] + '\t' + eachte + '\t' + te_infor_dic[eachte]['bg'] + '\t' + \
                                     te_infor_dic[eachte]['ed'] + '\t' + te_infor_dic[eachte]['rest']
                        TE_gene_insertion_dic[final_line] = 1


                ##if te cover downstream of mRNA
                if int(te_bg) >= int(mRNA_ed) and int(te_bg) <= int(mRNA_down):
                    if mRNA_dic[eachmRNA_ID]['dir'] == '+':
                        final_line = mRNA_chr + '\t' + eachmRNA_ID + '\t' + 'downstream' + '\t' + \
                                     str(mRNA_ed) + '\t' + str(mRNA_down) + '\t' + mRNA_dic[eachmRNA_ID]['dir'] + '\t' + \
                                     mRNA_dic[eachmRNA_ID]['annot'] + '\t' + eachte + '\t' + te_infor_dic[eachte]['bg'] + '\t' + \
                                     te_infor_dic[eachte]['ed'] + '\t' + te_infor_dic[eachte]['rest']
                        TE_gene_insertion_dic[final_line] = 1
                    else:
                        final_line = mRNA_chr + '\t' + eachmRNA_ID + '\t' + 'upstream' + '\t' + \
                                     str(mRNA_ed) + '\t' + str(mRNA_down) + '\t' + mRNA_dic[eachmRNA_ID]['dir'] + '\t' + \
                                     mRNA_dic[eachmRNA_ID]['annot'] + '\t' + eachte + '\t' + te_infor_dic[eachte]['bg'] + '\t' + \
                                     te_infor_dic[eachte]['ed'] + '\t' + te_infor_dic[eachte]['rest']
                        TE_gene_insertion_dic[final_line] = 1


        ##get the list element information
        elem_list = all_dic[eachmRNA_ID]

        for eachte in all_te_dic:
            #print (eachte)
            te_bg = te_infor_dic[eachte]['bg']
            te_ed = te_infor_dic[eachte]['ed']

            for eachelem_dic in elem_list:
                #print (eachelem_dic)

                ##get the elem position information
                elem_bg = eachelem_dic['bg']
                #print ('the elem bg is ' + str(elem_bg))
                elem_ed = eachelem_dic['ed']

                ##if te cover the elem
                if int(te_ed) >= int(elem_bg) and int(te_bg) <= int(elem_ed):
                    ##generate the final line
                    final_line = mRNA_chr + '\t' + eachmRNA_ID + '\t' + eachelem_dic['type'] + '\t' + \
                                 elem_bg + '\t' + elem_ed + '\t' + eachelem_dic['dir'] + '\t' + \
                                 eachelem_dic['annot'] + '\t' + eachte + '\t' + te_infor_dic[eachte]['bg'] + '\t' + \
                                 te_infor_dic[eachte]['ed'] + '\t' + te_infor_dic[eachte]['rest']
                    #print('the final line is '+ final_line)
                    TE_gene_insertion_dic[final_line] = 1

    return (TE_gene_insertion_dic)


##define a function to generate the table
def generate_pro_tb (TE_gene_insertion_dic,all_dic,mRNA_dic,gene_dic):

    ##initiate a dic to store pro and number information
    pro_type_line_dic = {}

    ##initiate a dic to store the number of type information
    all_type_dic = {}
    ##initiate a dic to store the number of mRNA information
    mRNA_type_dic = {}
    mRNA_te_type_dic = {}
    ##initiate a dic to store the number of gene information
    gene_type_dic = {}
    gene_te_type_dic = {}

    ##########################################################
    ##calculate the number of all element number in the genome
    ##get the unique number of ID from the TE_gene_insertion_dic
    for eachmRNA in all_dic:
        elem_list = all_dic[eachmRNA]
        for eachelem_dic in elem_list:
            type = eachelem_dic['type']
            ##get the total number of proporiton of type
            if type in all_type_dic:
                all_type_dic[type] += 1
            else:
                all_type_dic[type] = 1

    ##initiate a dic to store the ID information to get the unique one
    ID_dic = {}
    for eachline in TE_gene_insertion_dic:
        col = eachline.strip().split()
        ID = col[6]
        type = col[2]
        ID_dic[ID] = type

    ##initiate a dic to store the number of type covered by TE
    te_type_dic = {}
    for eachID in ID_dic:
        type = ID_dic[eachID]
        if type in te_type_dic:
            te_type_dic[type] += 1
        else:
            te_type_dic[type] = 1

    ##generate the proportion table
    for eachall_type in all_type_dic:
        eachall_type_num = all_type_dic[eachall_type]

        for eachte_type in te_type_dic:
            eachte_type_num = te_type_dic[eachte_type]

            if eachall_type == eachte_type:
                pro = int(eachte_type_num)/int(eachall_type_num)
                final_line = eachall_type + '\t' + str(pro)
                pro_type_line_dic[final_line] = 1

    ########################################
    ##calculate the mRNA and gene proportion
    ##store the number of all mRNA and genes
    ##for mRNA
    for eachmRNA in mRNA_dic:
        mRNA_type_dic[eachmRNA] = 1
    mRNA_all_num = len(list(mRNA_type_dic.keys()))

    for eachline in TE_gene_insertion_dic:
        col = eachline.strip().split()
        mRNA_nm = col[1]
        if col[2] != 'upstream' and col[2] != 'downstream':
            mRNA_te_type_dic[mRNA_nm] = 1
    mRNA_te_num = len(list(mRNA_te_type_dic.keys()))
    pro_mRNA = int(mRNA_te_num)/int(mRNA_all_num)
    final_line = 'mRNA' + '\t' + str(pro_mRNA)
    pro_type_line_dic[final_line] = 1

    ##for gene
    for eachgene in gene_dic:
        gene_type_dic[eachgene] = 1

        ##get the parent name of the mRNA
        for eachline in TE_gene_insertion_dic:
            col = eachline.strip().split()
            mRNA_nm = col[1]
            mRNA_annot = mRNA_dic[mRNA_nm]['annot']
            col_mRNA_annot = mRNA_annot.split(';')

            annot_dic = {}
            gene_nm = ''
            left_annot = ''

            if col[2] != 'upstream' and col[2] != 'downstream':

                for eachpart_annot in col_mRNA_annot:
                    mt = re.match('(.+)=(.+)', eachpart_annot)
                    left_annot = mt.group(1)
                    right_annot = mt.group(2)
                    annot_dic[left_annot] = right_annot

                if left_annot == 'Parent':
                    gene_nm = annot_dic[left_annot]

                if eachgene == gene_nm:
                    gene_te_type_dic[eachgene] = 1

    gene_all_num = len(list(gene_type_dic.keys()))
    gene_te_num = len(list(gene_te_type_dic.keys()))
    pro_gene = int(gene_te_num)/int(gene_all_num)
    final_line = 'gene' + '\t' + str(pro_gene)
    pro_type_line_dic[final_line] = 1

    ##for the gene genetic regions
    intergenetic_pro = 1 - float(pro_gene)
    final_line = 'inter_genetic_regions' + '\t' + str(intergenetic_pro)
    pro_type_line_dic[final_line] = 1



    ##add this to the insert_te_gene.py
    ##for the flanking region
    mRNA_te_up_dic = {}
    mRNA_te_down_dic = {}
    for eachline in TE_gene_insertion_dic:
        col = eachline.strip().split()
        if col[2] == 'upstream':
            mRNA_nm = col[1]
            mRNA_te_up_dic[mRNA_nm] = 1
        else:
            mRNA_nm = col[1]
            mRNA_te_down_dic[mRNA_nm] = 1

    mRNA_te_up_num = len(list(mRNA_te_up_dic.keys()))
    mRNA_te_down_num = len(list(mRNA_te_down_dic.keys()))
    pro_mRNA_up = int(mRNA_te_up_num) / int(mRNA_all_num)
    pro_mRNA_down = int(mRNA_te_down_num) / int(mRNA_all_num)
    final_line = 'upstream' + '\t' + str(pro_mRNA_up)
    pro_type_line_dic[final_line] = 1
    final_line = 'downstream' + '\t' + str(pro_mRNA_down)
    pro_type_line_dic[final_line] = 1




    return(pro_type_line_dic,ID_dic)




##define a function to analyze flanking regions, and users will provide flanking regions ranges, default is 3000 bp



















