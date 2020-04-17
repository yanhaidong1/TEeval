#!/usr/bin/env python3.7

##updation7.9 add TSD into each inserted
##the final.lib should contains the DNA mite from the PGSB dataset with lower than 800 bp
##should increase the insertion number of the DNA mite that could be better

##updation5.9, add a count selection for duplicating the inserted sequences
##updation4.27, change the threshold to generate another set of clean sequences


##This script is to generate gold reference TE in a clean sequences
##we will have 7500 copies

##mutation 500
##deletion 500
##insertion 500
##no mutation 500

##mutation + deletion 500
##mutation + insertion 500
##deletion + insertion 500
##mutation + insertion + deletion 500
##mutation + insertion + deletion + fusion 500


##import modules
import re
import sys
from Bio import SeqIO
import subprocess
import os
import random
import glob


##input the library that will provide the raw TE sequence  20 TE copies 25 times
input_library_file = sys.argv[1]

##we will insert the LTR in this fasta file
input_clean_fasta_file = sys.argv[2]

##input the name of TE that will be analyzed
TE_insert_name = sys.argv[3]

##location of the software script
location_software = sys.argv[4]

##store TE insert dir
#insert_TE_dir = sys.argv[5]

##generate a temp file to store the TE information
temp_bed_fasta_dir = sys.argv[5]

##denote a count that shows the replicate number
replicate_count = sys.argv[6]
##denote a count that shows the total number of the replicated TE
total_replicate_count = sys.argv[7]


##################################
##Step 1: generate TE copies (500)
##################################

#############################
##copy the LTRs into 25 times
##updation5.9 determine how many TEs we will duplicate
def generate_TE_copy (input_library_file,TE_insert_name,replicate_count,total_replicate_count):


    ##determine how many raw TE we will choose
    ##such as we will generate 500 copies and allows each gene replicate 25, so we need 20 TEs
    choose_raw_TE_count = int(total_replicate_count)/int(replicate_count)
    print('the chosen raw TE count is ' + str(choose_raw_TE_count))

    ##initiate a dic to store the TE fasta after copying
    TE_fasta_dic = {}
    raw_TE_count = 0
    for seq_record in SeqIO.parse(input_library_file,'fasta'):
        if TE_insert_name in seq_record.id:
            raw_TE_count += 1
            TE_fasta_dic[str(raw_TE_count)] = {'id':seq_record.id,'seq':str(seq_record.seq)}
            #TE_fasta_dic[seq_record.id] = seq_record.seq

    print('total TE fasta number is ' + str(raw_TE_count))
    ##randomly select the id that will be duplicated
    ##choose  raw_TE_count times ex. choose copied 25 time
    my_randoms = []
    for i in range(int(choose_raw_TE_count)):

        random_num = random.randrange(1, raw_TE_count + 1, 1)

        while random_num in my_randoms:
            random_num = random.randrange(1, raw_TE_count + 1, 1)

        my_randoms.append(random_num)

    print('the copy number is ' + str(len(my_randoms)))

    ##initiate a dic to store choose_raw_TE fasta
    raw_TE_fasta_dic = {}
    for each_count_id in my_randoms:
        raw_TE_fasta_dic[TE_fasta_dic[str(each_count_id)]['id']] = TE_fasta_dic[str(each_count_id)]['seq']


    ##select the first 20 in the library
    #id_count = 0
    #for eachid in TE_fasta_dic:
    #   id_count += 1
    #   if id_count <= 20:
    #       TE_20te_fasta_dic[eachid] = TE_fasta_dic[eachid]


    ##initiate a dic to store TE fasta after copying
    TE_copy_fasta_dic = {}
    ##copy 25 times
    for eachid in raw_TE_fasta_dic:
        for i in range(1,int(replicate_count)+1):
            ##get the new name of TE
            new_id = eachid + '_' + str(i)
            TE_copy_fasta_dic[new_id] = raw_TE_fasta_dic[eachid]

    print('the new TE copy number is ' + str(len(list(TE_copy_fasta_dic.keys()))))

    return (raw_TE_fasta_dic,TE_copy_fasta_dic)


##################################################################
##insert into a clean sequence that is ready to go to the software
def insert_TEs (TE_copy_fasta_dic,input_clean_fasta_file,inter_insert_len,start_insert_loc):


    ##initiate a dic to store the final genome fasta file
    final_genome_fasta_dic = {}

    ##initiate a dic to store the inserted seq location
    insert_seq_loc_dic = {}

    ##store the seq of input_genome as a string
    raw_seq_string = ''
    clean_seq_string = ''
    for seq_record in SeqIO.parse(input_clean_fasta_file,'fasta'):
        raw_seq_string = str(seq_record.seq)
        clean_seq_string = str(seq_record.seq)


    count = 0
    ##initiate a insert_seq as empty
    insert_seq = ''
    ##initiate a insert_location to indicate the start of insertion
    insert_location = int(start_insert_loc)

    for eachid in TE_copy_fasta_dic:

        count += 1
        if count == 1:
            insert_location = insert_location + len(insert_seq)
        else:
            ##calculate the insert location
            insert_location = insert_location + int(inter_insert_len) + len(insert_seq)


        insert_seq = str(TE_copy_fasta_dic[eachid])

        ##divide the target genome sequence as two parts
        seq_part1 = raw_seq_string[:insert_location]
        seq_part2 = raw_seq_string[insert_location:]

        raw_seq_string = seq_part1 + insert_seq + seq_part2

        ##get the location information from each inserted seq and store them in to a dic
        seq_insert_line = eachid + '\t' + str(insert_location) + '\t' + str((int(insert_location) + len(insert_seq))) + \
                          '\t' + insert_seq
        insert_seq_loc_dic[seq_insert_line] = 1


    ##store the final fasta file
    final_genome_fasta_dic['final_seq'] = raw_seq_string
    len_clean_seq_string = len(clean_seq_string)


    return (final_genome_fasta_dic, insert_seq_loc_dic,len_clean_seq_string)



##############################################################
##Step 2: run software to generate mutation,deletion,insertion
##############################################################

###########################################
##generate gff file from insert_Seq_loc_dic
def generate_gff_file (insert_seq_loc_dic,len_clean_seq_string):

    store_gff_line_dic = {}

    ##add a first line in the top of the gff file
    first_line = 'final_seq' + '\t' +  'Genbank' + '\t' + 'region' + '\t' + '1' + '\t' + str(len_clean_seq_string) + '\t' + '.' + '\t' + \
        '+' + '\t' + '.' + '\t' + 'ID=gene0;Name=ClassI//LTR/ltr_Others_TE1_0'   ##for the clean sequence length

    store_gff_line_dic[first_line] = 1

    count = 0
    for eachline in insert_seq_loc_dic:

        count += 1

        col = eachline.strip().split()
        TE_nm = col[0]
        TE_start = col[1]
        TE_end = col[2]
        direct = '+'
        chr = 'final_seq'
        region_name = 'gene' ##should be consist with the software

        ##for the annotation part
        ID = 'gene' + str(count)

        gff_line = chr + '\t' + 'Genbank' + '\t' + region_name + '\t' + TE_start + '\t' + TE_end + '\t' + '.' + '\t' + \
                    direct + '\t' + '.' + '\t' + 'ID=' + ID + ';Name=' + TE_nm

        store_gff_line_dic[gff_line] = 1

    ##write ### in the end of the gff file
    last_gff_line = '###'
    store_gff_line_dic[last_gff_line] = 1

    return (store_gff_line_dic)

#################################################
##run the software simulome to generate variation
def run_software (opt_temp_final_fasta_file,opt_temp_gff_file,location_software,TE_copy_fasta_dic):

    ##opt_temp_final_fasta_file is opt_temp_final.fasta
    ##opt_temp_gff_file will be generated from store_gff_line_dic

    ##initiate a dic to store all the TEs
    all_dic = {}

    ##########################################
    ##calculate the average length of all LTRs
    ##target the opt_temp_gff_file
    total_length = 0
    count = 0
    with open (opt_temp_gff_file,'r') as ipt:
        for eachline in ipt:

            if not eachline.startswith('#'):

                eachline = eachline.strip('\n')
                col = eachline.strip().split()

                ##calculate the total length
                name = col[2]

                if name != 'region':
                    start = col[3]
                    end = col[4]
                    length_te = int(end) - int(start)

                    count += 1
                    total_length = total_length + length_te

    average_length = int(total_length/count)
    print('the average length is ' + str(average_length))

    ##################
    ##for the mutation
    ##################
    ##make mutation_dir
    #mutation_dir = 'mutation'
    #if not os.path.exists(mutation_dir):
    #    os.makedirs(mutation_dir)

    ##for the mutation it will generate 10% 20% and 30%
    ##for the --num_snp
    ##updation4.27 change to 0.05, 0.1, 0.2
    list_pro = [0.01,0.05,0.10]  ###this can be changed

    num_snp_std_dev = int(average_length*0.005)   ###this can be changed


    #pro_number = 0
    for eachnum in list_pro:

        #pro_number += 1
        ##make specific dir for in the mutation dir
        pro_dir = str(eachnum) + '_num_mutation'
        if not os.path.exists(pro_dir):
            os.makedirs(pro_dir)

        num_snp = int(average_length*eachnum)

        print('run mutation')
        cmd = 'python2.7 ' + location_software + \
              ' --genome=' + opt_temp_final_fasta_file + \
              ' --anno=' + opt_temp_gff_file + \
              ' --output=' + pro_dir + \
              ' --whole_genome=TRUE' + \
              ' --snp=TRUE' + \
              ' --snp_distrib=TRUE' + \
              ' --num_snp=' + str(num_snp) + \
              ' --snp_std_dev=' + str(num_snp_std_dev)
        print(cmd)
        subprocess.call(cmd, shell=True)


    ##################
    ##for the deletion
    ##################
    ##make deletion_dir
    #deletion_dir =   'deletion'
    #if not os.path.exists(deletion_dir):
    #    os.makedirs(deletion_dir)

    ##for the mutation it will generate 10% 20% and 30%
    ##for the --num_snp
    list_pro = [0.01,0.05,0.10]    ###this can be changed

    num_del_std_dev = int(average_length*0.005)   ###this can be changed

    #pro_number = 0
    for eachnum in list_pro:

        #pro_number += 1
        ##make specific dir for in the mutation dir
        pro_dir = str(eachnum) + '_num_deletion'
        if not os.path.exists(pro_dir):
            os.makedirs(pro_dir)

        del_len = int(average_length*eachnum)

        print('run deletion')
        cmd = 'python2.7 ' + location_software + \
              ' --genome=' + opt_temp_final_fasta_file + \
              ' --anno=' + opt_temp_gff_file + \
              ' --output=' + pro_dir + \
              ' --whole_genome=TRUE' + \
              ' --indel=2' + \
              ' --del_len=' + str(del_len) + \
              ' --del_distrib=TRUE' + \
              ' --num_del=2 ' + \
              ' --del_std_dev=' + str(num_del_std_dev)
        print(cmd)
        subprocess.call(cmd, shell=True)


    ###################
    ##for the insertion
    ###################
    ##make insertion_dir
    #insertion_dir = 'insertion'
    #if not os.path.exists(insertion_dir):
    #    os.makedirs(insertion_dir)

    ##for the mutation it will generate 10% 20% and 30%
    ##for the --num_snp
    list_pro = [0.01,0.05,0.10]    ###this can be changed

    num_ins_std_dev = int(average_length*0.005)   ###this can be changed

    #pro_number = 0
    for eachnum in list_pro:

        #pro_number += 1

        ##make specific dir for in the mutation dir
        pro_dir = str(eachnum) + '_num_insertion'
        if not os.path.exists(pro_dir):
            os.makedirs(pro_dir)

        ins_len = int(average_length*eachnum)

        print('run insertion')
        cmd = 'python2.7 ' + location_software + \
              ' --genome=' + opt_temp_final_fasta_file + \
              ' --anno=' + opt_temp_gff_file + \
              ' --output=' + pro_dir + \
              ' --whole_genome=TRUE' + \
              ' --indel=1' + \
              ' --ins_len=' + str(ins_len) + \
              ' --ins_distrib=TRUE' + \
              ' --num_ins=2 ' + \
              ' --ins_std_dev=' + str(num_ins_std_dev)
        print(cmd)
        subprocess.call(cmd, shell=True)

    ##############################
    ##for the mutation + insertion
    ##############################
    ##make mutation_insertion_dir
    mutation_insertion_dir = 'mutation_insertion'
    if not os.path.exists(mutation_insertion_dir):
        os.makedirs(mutation_insertion_dir)

    ##no loop for the mutation and insertion
    num_snp = int(average_length*0.01)
    ins_len = int(average_length*0.01)

    print('run mutation + insertion')
    cmd = 'python2.7 ' + location_software + \
          ' --genome=' + opt_temp_final_fasta_file + \
          ' --anno=' + opt_temp_gff_file + \
          ' --output=' + mutation_insertion_dir + \
          ' --whole_genome=TRUE' + \
          ' --snp=TRUE' + \
          ' --snp_distrib=TRUE' + \
          ' --num_snp=' + str(num_snp) + \
          ' --snp_std_dev=' + str(num_snp_std_dev) + \
          ' --indel=1' + \
          ' --ins_len=' + str(ins_len) + \
          ' --ins_distrib=TRUE' + \
          ' --num_ins=2 ' + \
          ' --ins_std_dev=' + str(num_ins_std_dev)
    print(cmd)
    subprocess.call(cmd, shell=True)


    #############################
    ##for the mutation + deletion
    #############################
    ##make mutation_deletion_dir
    mutation_deletion_dir = 'mutation_deletion'
    if not os.path.exists(mutation_deletion_dir):
        os.makedirs(mutation_deletion_dir)

    ##no loop for the mutation and deletion
    num_snp = int(average_length * 0.01)
    del_len = int(average_length * 0.01)

    print('run mutation + deletion')
    cmd = 'python2.7 ' + location_software + \
          ' --genome=' + opt_temp_final_fasta_file + \
          ' --anno=' + opt_temp_gff_file + \
          ' --output=' + mutation_deletion_dir + \
          ' --whole_genome=TRUE' + \
          ' --snp=TRUE' + \
          ' --snp_distrib=TRUE' + \
          ' --num_snp=' + str(num_snp) + \
          ' --snp_std_dev=' + str(num_snp_std_dev) + \
          ' --indel=2' + \
          ' --del_len=' + str(del_len) + \
          ' --del_distrib=TRUE' + \
          ' --num_del=2 ' + \
          ' --del_std_dev=' + str(num_del_std_dev)
    print(cmd)
    subprocess.call(cmd, shell=True)


    ##############################
    ##for the deletion + insertion
    ##############################
    ##make mutation_deletion_dir
    insertion_deletion_dir = 'insertion_deletion'
    if not os.path.exists(insertion_deletion_dir):
        os.makedirs(insertion_deletion_dir)

    ##no loop for the mutation and deletion
    del_len = int(average_length * 0.01)
    ins_len = int(average_length * 0.01)

    print('run insertion + deletion')
    cmd = 'python2.7 ' + location_software + \
          ' --genome=' + opt_temp_final_fasta_file + \
          ' --anno=' + opt_temp_gff_file + \
          ' --output=' + insertion_deletion_dir + \
          ' --whole_genome=TRUE' + \
          ' --indel=3' + \
          ' --del_len=' + str(del_len) + \
          ' --del_distrib=TRUE' + \
          ' --num_del=2 ' + \
          ' --del_std_dev=' + str(num_del_std_dev) + \
          ' --ins_len=' + str(ins_len) + \
          ' --ins_distrib=TRUE' + \
          ' --num_ins=2 ' + \
          ' --ins_std_dev=' + str(num_ins_std_dev)
    print(cmd)
    subprocess.call(cmd, shell=True)

    #########################################
    ##for the mutation + deletion + insertion
    #########################################
    mutation_deletion_insertion_dir = 'mutation_deletion_insertion'
    if not os.path.exists(mutation_deletion_insertion_dir):
        os.makedirs(mutation_deletion_insertion_dir)

    ##no loop for the mutation and deletion
    num_snp = int(average_length * 0.01)
    del_len = int(average_length * 0.01)
    ins_len = int(average_length * 0.01)

    print('run mutation + insertion + deletion')
    cmd = 'python2.7 ' + location_software + \
          ' --genome=' + opt_temp_final_fasta_file + \
          ' --anno=' + opt_temp_gff_file + \
          ' --output=' + mutation_deletion_insertion_dir + \
          ' --whole_genome=TRUE' + \
          ' --snp=TRUE' + \
          ' --snp_distrib=TRUE' + \
          ' --num_snp=' + str(num_snp) + \
          ' --snp_std_dev=' + str(num_snp_std_dev) + \
          ' --indel=3' + \
          ' --del_len=' + str(del_len) + \
          ' --del_distrib=TRUE' + \
          ' --num_del=2 ' + \
          ' --del_std_dev=' + str(num_del_std_dev) + \
          ' --ins_len=' + str(ins_len) + \
          ' --ins_distrib=TRUE' + \
          ' --num_ins=2 ' + \
          ' --ins_std_dev=' + str(num_ins_std_dev)
    print(cmd)
    subprocess.call(cmd, shell=True)


    #####################################
    ##extract information from output dir
    #####################################
    insert_TE_dir_list = glob.glob('./*tion')

    ##select the eachfl that will not be the 01 02 03
    ##for the 010203 TE_dir
    for eachdir in insert_TE_dir_list:

        ##get the name of the eachdir
        mt = re.match('.+/(.+)',eachdir)
        dir_nm = mt.group(1)


        ##initaite a sub_dic {'1' : {'te' : xx, 'seq' : sequence},'2' : {'te' : xx, 'seq' : sequence}
        sub_dic = {}

        ##consider combinatnion situations
        #if dir_nm != 'deletion' and dir_nm != 'mutation' and dir_nm != 'insertion':

        file_list = glob.glob(eachdir + '/*mutated_simulation*')

        for eachfl in file_list:
            ##generate the bedfile from gff and extract the sequence from the bed file
            if 'mutated_simulation.gff' in eachfl:

                print('the analzyed file is ' + eachfl)
                ##generate a bed dic to store the bed line
                bed_line_dic = {}
                with open (eachfl,'r') as ipt:
                    for eachline in ipt:

                        eachline = eachline.strip('\n')

                        if not eachline.startswith('#'):

                            col = eachline.strip().split()

                            chr = col[0]
                            start = col[3]
                            end = col[4]
                            direct = col[6]
                            type = col[2]

                            if type == 'gene':

                                annot_col = col[8].split(';')
                                mt = re.match('(.+)=(.+)',annot_col[1])
                                te_nm = mt.group(2)

                                bed_line = chr + '\t' + start + '\t' + end + '\t' + te_nm + '\t' + '1' + '\t' + direct
                                bed_line_dic[bed_line] = 1

                ##wirte a file in the temp dir
                with open (temp_bed_fasta_dir + '/temp.bed','w+') as opt:
                    for eachline in bed_line_dic:
                        opt.write(eachline + '\n')

        for eachfl in file_list:

            if 'mutated_simulation.fasta' in eachfl and 'mutated_simulation.fasta.fai' not in eachfl:

                ##generate the fasta file
                cmd = 'bedtools getfasta -fi ' + eachfl + ' -bed ' + temp_bed_fasta_dir + '/temp.bed -s -name -fo ' + \
                      temp_bed_fasta_dir + '/opt_te.fasta'
                print (cmd)
                subprocess.call(cmd, shell=True)


                ##next step is to extract the copies information including name and seq into a dic
                ##that contains all the information for the tes
                ##all_dic = {'deletion':{'1' : {'te' : xx, 'seq' : sequence},'2' : {'te' : xx, 'seq' : sequence}....},
                ##'mutation':{'1' : {'te' : xx, 'seq' : sequence},.....}

                ##extract the sequence information
                ##caculate the te_count
                te_count = 0
                for seq_record in SeqIO.parse(temp_bed_fasta_dir + '/opt_te.fasta','fasta'):

                    te_count += 1

                    id_code = str(te_count)
                    id_te = seq_record.id
                    seq_te = seq_record.seq

                    ##store them into the sub_dic {'1' : {'te' : xx, 'seq' : sequence},'2' : {'te' : xx, 'seq' : sequence}
                    sub_dic[id_code] = {'te':id_te,'seq':str(seq_te)}

        all_dic[dir_nm] = sub_dic


    ##store the no variation te in the all_dic that will be analyzed in the fusion
    no_var_te_count = 0
    no_var_sub_dic = {}
    for eachte in TE_copy_fasta_dic:
        no_var_te_count += 1

        id_code = str(no_var_te_count)
        id_te = eachte
        seq_te = str(TE_copy_fasta_dic[eachte])

        ##store them into the sub_dic {'1' : {'te' : xx, 'seq' : sequence},'2' : {'te' : xx, 'seq' : sequence}
        no_var_sub_dic[id_code] = {'te': id_te, 'seq': str(seq_te)}

    all_dic['novariation'] = no_var_sub_dic

    return (all_dic)



##write a script fuse two TEs selected from
##for the fusion, we will generate 10% of each big catogories as fusion status
##for the
# mutation,
# deletion,
# insertion,
# mutation_deletion,
# mutation_insertion,
# mutation_insertion_deletion
##for these six classes

##define a function to generate three types of fusions
##type1: delete one end and delete start of another, and combine
##type2: delete the second start and combine
##type3: delete the first end and combine
def fusion_type (two_fusion_id_list,eachtype,all_dic,type_number,deletion_rate):

    ##deletion rate = 0.5,  means that the random loc will start from the 0.5*first_te_len to the end
    ##deletion rate = 0.9, the random_loc = (len*0.9,first_te_len) for the left seq
    ##the random_loc = (0,len*0.1) for the right seq


    ##initiate a string to store the combine_seq
    combine_seq_dic = {}

    ##get the length of the first te
    first_seq = all_dic[eachtype][two_fusion_id_list[0]]['seq']
    first_te_len = len(all_dic[eachtype][two_fusion_id_list[0]]['seq'])
    ##get the random deletion for the end
    random_loc = random.randint(int(first_te_len*deletion_rate), first_te_len)  ##from 0 to 499 for the list
    ##get the left sequence after cut the down the random_loc to the end
    left_seq = all_dic[eachtype][two_fusion_id_list[0]]['seq'][0:random_loc]

    ##get the length of the second te
    second_seq = all_dic[eachtype][two_fusion_id_list[1]]['seq']
    second_te_len = len(all_dic[eachtype][two_fusion_id_list[1]]['seq'])
    ##get the random deletion for the start
    random_loc = random.randint(0,int(second_te_len * (1-deletion_rate)))  ##from 0 to 499 for the list
    right_seq = all_dic[eachtype][two_fusion_id_list[1]]['seq'][random_loc:second_te_len]

    ##splite the fusion into two copies instead of one copy
    if str(type_number) == '1':
        combine_seq_dic = {'lnm':all_dic[eachtype][two_fusion_id_list[0]]['te'],'lseq':left_seq,
                           'rnm':all_dic[eachtype][two_fusion_id_list[1]]['te'],'rseq':right_seq}

        #combine_seq = left_seq + right_seq

    if str(type_number) == '2':
        combine_seq_dic = {'lnm':all_dic[eachtype][two_fusion_id_list[0]]['te'],'lseq':first_seq,
                           'rnm':all_dic[eachtype][two_fusion_id_list[1]]['te'],'rseq':right_seq}

        #combine_seq = first_seq + right_seq

    if str(type_number) == '3':
        combine_seq_dic = {'lnm':all_dic[eachtype][two_fusion_id_list[0]]['te'],'lseq':left_seq,
                           'rnm':all_dic[eachtype][two_fusion_id_list[1]]['te'],'rseq':second_seq}


        #combine_seq = left_seq + second_seq

    return (combine_seq_dic)

def generate_fusion_seq (all_dic,total_replicate_count):

    ##initiate a dic to store the fusion dic
    fusion_dic = {}
    ##key is name of fusion te
    ##fusion_dic = {'deletion':{'fus_te_1' : seq}}

    ##10% date for each category

    ##all_dic contains information storing all the seq information
    for eachtype in all_dic:

        ##generate a list to contain all the ids
        ##for the part, we will randomly extract two id from the id list
        #id_list = list(all_dic[eachtype].keys())

        ##randomly select 10% from the id_list
        reduce_number = int(int(total_replicate_count)/10)
        fusion_number = reduce_number ##set 50 because 10% of the 500 copies as fusion number

        fusion_id_dic = {} ##this dic contains information for the fusion ids that will generate paired fusion te

        for x in range(fusion_number):

            ##change location mutation_len times:
            ##this is for eachtime
            random_loc = random.randint(1, int(total_replicate_count))  ##from 1 to 500 for the list ##id in all_dic[eachtype] has no zero
            ##check if random_loc in the occur_number_dic
            ##if in, run again,
            ##if not, the random_loc is the right one
            while str(random_loc) in fusion_id_dic.keys():

                ##from the 0 location to the seq_len-1 for the seq string
                random_loc = random.randint(1,int(total_replicate_count))

            fusion_id_dic[str(random_loc)] = 1


        ##extract two ids that will generate fusion
        fusion_id_list = list(fusion_id_dic.keys())

        ##set a loop to analyze fusion list
        ##set a fusion name
        fusion_count = 0
        for x in range(int(reduce_number/2)):

            fusion_count += 1

            two_fusion_id_list = random.sample(fusion_id_list, 2) ##this two samples are different

            fusion_id_list.remove(two_fusion_id_list[0])
            fusion_id_list.remove(two_fusion_id_list[1])

            ##now the two fusion id has got
            ##randomly select the fusion type in the fusion_type function
            random_type_number = random.randint(1, 3)
            combine_seq_dic = fusion_type (two_fusion_id_list,eachtype,all_dic,random_type_number,0.9) ##set the 0.9 for the deletion

            ##store the combine_seq to the store_dic
            ##fusion_nm should be the whole name for the two fusion te
            ##store left name
            fusion_dic[combine_seq_dic['lnm'] + '_fusion_' + eachtype] = combine_seq_dic['lseq']
            fusion_dic[combine_seq_dic['rnm'] + '_fusion_' + eachtype] = combine_seq_dic['rseq']

            ##store right name
            #fusion_nm = all_dic[eachtype][str(two_fusion_id_list[0])]['te'] + ';' + all_dic[eachtype][str(two_fusion_id_list[1])]['te'] + '_' + eachtype
            #fusion_dic[fusion_nm] = combine_seq


    ##becareful to store the eachtype to the fusion_dic

    return (fusion_dic)

###########################################
##change the store format in the all_te_dic
##this dic should be as the input for the randomly insertion
##A dic contains each sequences
##A = {'fullname_mutation':'ATCCG',te2_mutation:'GCCCG',te3_mutation:'GGGCG'}
##B list contains nul for the clean sequences
##B = ['A','T','G','G','C','C','C']

def change_all_dic_format (all_dic,fusion_dic,TE_copy_fasta_dic):

    ##initiate a dic to store all the information
    all_te_dic = {}

    ##store te information from the all_dic
    for eachtype in all_dic:

        ##all_dic = {'deletion':{'1' : {'te' : xx, 'seq' : sequence},'2' : {'te' : xx, 'seq' : sequence}....},
        for eachid in all_dic[eachtype]:

            ##generate full name that will be the key in the all_te_dic
            full_name = all_dic[eachtype][eachid]['te'] + '_' + eachtype
            te_seq = all_dic[eachtype][eachid]['seq']

            all_te_dic[full_name] = te_seq

    ##store the te information from the fusion_dic
    for eachnm in fusion_dic:

        all_te_dic[eachnm] = fusion_dic[eachnm]


    ##store the te without variation from the TE_copy_fasta_dic
    for eachid in TE_copy_fasta_dic:

        full_name = eachid + '_novariation'
        all_te_dic[full_name] = str(TE_copy_fasta_dic[eachid])


    return (all_te_dic)


##################################################
##Step 3: randomly insert into the clean sequences
##################################################
##generate the A dic and B list
##A dic contains each sequences
##A = {'te1':'ATCCG',te2:'GCCCG',te3:'GGGCG'}
##B list contains nul for the clean sequences
##B = ['A','T','G','G','C','C','C']

##updation 7.9 add TSD for the end of each insertion MITE
##updation 7.14 add TSD for the end of each insertion LTR
##generate a function to randomly create the TSD for the the LTR
def randomString(stringLength):
    """Generate a random string with the combination of lowercase and uppercase letters """
    letters = 'ATGC'
    return ''.join(random.choice(letters) for i in range(stringLength))


def inserte_generate_te (all_te_dic,input_clean_fasta_file):

    ##initate a dic to store the inserted te_seq loc information
    te_location_infor_dic = {}

    ##generate the dic to store the final seq
    te_final_seq_list = []
    te_final_seq_dic = {}
    clean_seq_list = []

    clean_seq = ''
    for seq_record in SeqIO.parse(input_clean_fasta_file,'fasta'):
        clean_seq = str(seq_record.seq)

    for i in range(0, len(clean_seq)):
        clean_seq_list.append(clean_seq[i])
        i += 1

    id_dic = {}
    for x in range(len(all_te_dic.keys())):

        te_key_list = list(all_te_dic.keys())

        te_key_in = te_key_list[random.randint(0, len(all_te_dic.keys()) - 1)]

        while te_key_in in id_dic.keys():
            te_key_in = te_key_list[random.randint(0, len(all_te_dic.keys()) - 1)]

        id_dic[te_key_in] = 1

        ##insert the sequence with name
        ##updation7.14 will randomly generate 4-6 bp TSD
        ##updation7.9 detects the mite and nmite
        ##updation7.9: add TSD; add TA since TA will be considered in the MITEfinder2 and MITEhunter

        ##randomly generate a TSD length
        tsd_len = random.choice('456')
        tsd_seq = randomString(int(tsd_len))

        #add_TSD_seq = 'TA' + str(all_te_dic[te_key_in]) + 'TA'
        add_TSD_seq = tsd_seq + str(all_te_dic[te_key_in]) + tsd_seq

        #insert_name_seq = te_key_in + ':' + all_te_dic[te_key_in]
        insert_name_seq = te_key_in + ':' + add_TSD_seq

        clean_seq_list.insert(random.randint(0, len(clean_seq_list)), insert_name_seq)

    ##generate the location for the inserted TEs
    ##generate a dic to store the te that have been detected
    detect_len = 0

    count_single = 0  ##for the each nul in the clean seq

    for eachitem in clean_seq_list:

        if len(eachitem) == 1: ##for each single nul in the list
            count_single += 1
            te_final_seq_list.append(eachitem)

        if len(eachitem) > 1: ##for the inserted seq
            #print(eachitem)
            if ':' in eachitem:
                mt = re.match('(.+):(.+)', eachitem)
                te_name = mt.group(1)
                #print(te_name)
                te_seq = mt.group(2)  ##the te_seq has already contains TA and TA

                te_final_seq_list.append(te_seq)

                te_loc_start = count_single + 1 + detect_len
                te_loc_end = te_loc_start + len(te_seq) - 1

                detect_len = detect_len + len(te_seq)

                ##generate the loc_line

                te_loc_line = te_name + '\t' + str(te_loc_start) + '\t' + str(te_loc_end)
                te_location_infor_dic[te_loc_line] = 1


    te_final_seq = ''.join(te_final_seq_list)

    ##store them into a dic
    te_final_seq_dic['final_insert_seq'] = te_final_seq

    return (te_location_infor_dic,te_final_seq_dic)


###############################
##Step 4: write out all results
###############################

##################################
##Step 1: generate TE copies (500)
print('run the first step: generate TE copies')

raw_TE_fasta_dic,TE_copy_fasta_dic = generate_TE_copy (input_library_file,TE_insert_name,replicate_count,total_replicate_count)
final_genome_fasta_dic,insert_seq_loc_dic,len_clean_seq_string = insert_TEs (TE_copy_fasta_dic,input_clean_fasta_file,5000,5000)

with open ('opt_temp_final.fasta','w+') as opt:
    for eachkey in final_genome_fasta_dic:
        opt.write('>' + eachkey + '\n' + final_genome_fasta_dic[eachkey])

with open ('opt_temp_insert_te_location.txt','w+') as opt:
    for eachline in insert_seq_loc_dic:
        opt.write(eachline + '\n')

with open ('opt_temp_20.lib','w+') as opt:
    for eachid in raw_TE_fasta_dic:
        opt.write('>' + eachid + '\n' + str(raw_TE_fasta_dic[eachid]) + '\n')

##############################################################
##Step 2: run software to generate mutation,deletion,insertion
print('run the second step: software to generate mutation,deletion,insertion')


store_gff_line_dic = generate_gff_file (insert_seq_loc_dic,len_clean_seq_string)

with open ('opt_temp_te.gff','w+') as opt:
    for eachline in store_gff_line_dic:
        opt.write(eachline + '\n')

all_dic = run_software ('opt_temp_final.fasta','opt_temp_te.gff',location_software,TE_copy_fasta_dic)

##store the all_dic
##all_dic = {'deletion':{'1' : {'te' : xx, 'seq' : sequence},'2' : {'te' : xx, 'seq' : sequence}....},
#with open('opt_test_all_dic.txt', 'w+') as opt:
#    for eachtype in all_dic:

#        for eachid in all_dic[eachtype]:

#            opt.write(eachtype + '\t' + )



fusion_dic = generate_fusion_seq (all_dic,total_replicate_count)
all_te_dic = change_all_dic_format (all_dic,fusion_dic,TE_copy_fasta_dic)

with open ('opt_temp_all_insert_te.txt','w+') as opt:
    for eachname in all_te_dic:
        opt.write(eachname + '\t' + all_te_dic[eachname] + '\n')


with open ('opt_temp_all_insert_te_name.txt','w+') as opt:
    for eachname in all_te_dic:
        opt.write(eachname + '\n')

##################################################
##Step 3: randomly insert into the clean sequences
print('run the third step: randomly insert into the clean sequences')

te_location_infor_dic,te_final_seq_dic = inserte_generate_te (all_te_dic,input_clean_fasta_file)

with open ('opt_te_location_infor.txt','w+') as opt:
    for eachline in te_location_infor_dic:
        opt.write(eachline + '\n')

with open ('opt_final_te.fasta','w+') as opt:
    for eachline in te_final_seq_dic:
        opt.write('>' + eachline + '\n' + te_final_seq_dic[eachline] + '\n')



