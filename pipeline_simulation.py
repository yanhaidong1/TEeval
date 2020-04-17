#!/usr/bin/env python3.7

##this script will generate a pipleline to simulate multiple times
##import modules
import sys
import subprocess
import os

##generate the simulation information

input_simulate_time = sys.argv[1]
replicate_count = sys.argv[2]
total_replicate_count = sys.argv[3]
store_simulate_results_dir = sys.argv[4]

##required files contains:
##mutation_table.dat  okay
##clean_noTE.fasta   okay
##opt_Gypsy_Copia.lib for LTR contains 500 Gypsy and 500 Copia okay
##final.lib, for DNA use  ClassII//DNA_MITE, for SINE use ClassI//nLTR/SINE  for LINE use ClassI//nLTR/LINE okay
##pipeline_generate_gold_reference_TE okay
##batch_running_evaluation.py  batch_running_landscape.py batch_combine_result.py okay
##raw_evalution_dir okay
##raw_landscape_dir okay
##TEScape okay
raw_simulate_required_fl_dir = sys.argv[5]
software_simulome_path = sys.argv[6]
all_simulate_output_dir = sys.argv[7]

##updation 032720 change the  family file name
##LTR for the ltrfinder and LTR harvest two lib
##and LTR is combining two tools
analyze_family_file = sys.argv[8] ##that should be LTR\tDNA

##updation add the clean_noTE.fasta as a argument
clean_noTE_fl = sys.argv[9]


def mul_simulate (replicate_count,total_replicate_count,input_simulate_time,store_simulate_results_dir,raw_simulate_required_fl_dir,software_simulome_path,all_simulate_output_dir,analyze_family_file):


    TE_family = []
    with open (analyze_family_file,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            TE_family.append(col[0])

    print(TE_family)

    #TE_family = ['LTR','DNA','LINE','SINE']
    ##construct four dir to contains 10 times simulation for each family
    LTR_output_dir = all_simulate_output_dir + '/LTR'
    if not os.path.exists(LTR_output_dir):
        os.makedirs(LTR_output_dir)
    LINE_output_dir = all_simulate_output_dir + '/LINE'
    if not os.path.exists(LINE_output_dir):
        os.makedirs(LINE_output_dir)
    SINE_output_dir = all_simulate_output_dir + '/SINE'
    if not os.path.exists(SINE_output_dir):
        os.makedirs(SINE_output_dir)
    DNA_output_dir = all_simulate_output_dir + '/DNA'
    if not os.path.exists(DNA_output_dir):
        os.makedirs(DNA_output_dir)


    ##create multiple dir to store the simulated data
    for i in range(int(input_simulate_time)):

        ##create the different times of simulation file containing simulate_dir and evaluation_dir
        simulate_dir = store_simulate_results_dir + '/simulate_' + str(i)
        if not os.path.exists(simulate_dir):
            os.makedirs(simulate_dir)


        ##create simulation_seq dir
        simulate_seq_dir = simulate_dir + '/simulate_seq_dir'
        if not os.path.exists(simulate_seq_dir):
            os.makedirs(simulate_seq_dir)

        ##generate the raw_data file in the simulation_seq dir
        raw_data_dir = simulate_seq_dir + '/raw_data_dir'
        if not os.path.exists(raw_data_dir):
            os.makedirs(raw_data_dir)

        ##cp required files to the raw_data_dir
        cmd = 'cp ' + clean_noTE_fl  + ' ' + \
              raw_simulate_required_fl_dir + '/final.lib ' + \
              raw_simulate_required_fl_dir + '/opt_Gypsy_Copia.lib ' + \
              raw_simulate_required_fl_dir + '/opt_MITE.lib ' + raw_data_dir + '/'
        print(cmd)
        subprocess.call(cmd, shell=True)

        ##create the evaluation dir
        evaluation_dir = simulate_dir + '/evaluation_dir'
        if not os.path.exists(evaluation_dir):
            os.makedirs(evaluation_dir)

        ##generate the different TE family dir
        for eachfm in TE_family:

            ############################################
            ##step 1: analyze for simulation of sequence
            ##generate eachfm dir
            eachfm_dir = simulate_seq_dir + '/' + eachfm
            if not os.path.exists(eachfm_dir):
                os.makedirs(eachfm_dir)

            ##cp the required files to eachfm dir
            cmd = 'cp ' + raw_simulate_required_fl_dir + '/pipeline_generate_gold_reference_TE.py ' + \
                  raw_simulate_required_fl_dir + '/mutation_table.dat ' + \
                  raw_simulate_required_fl_dir + '/opt_Gypsy_Copia.lib ' + \
                  raw_simulate_required_fl_dir + '/final.lib ' + \
                  raw_simulate_required_fl_dir + '/opt_MITE.lib ' + eachfm_dir + '/'
            print(cmd)
            subprocess.call(cmd,shell=True)


            ##decide which lib we will use

            ##decide which name we should look for
            search_name = ''
            if eachfm == 'LTR':
                search_name = 'LTR'
            if eachfm == 'LINE':
                search_name = 'ClassI//nLTR/LINE'
            if eachfm == 'SINE':
                search_name = 'ClassI//nLTR/SINE'
            if eachfm == 'DNA':
                search_name = 'ClassII//DNA_MITE'

            ##make a temp dir to store the temp_bed_fasta_dir
            temp_bed_fasta_dir = eachfm_dir + '/temp_bed_fasta_dir'
            if not os.path.exists(temp_bed_fasta_dir):
                os.makedirs(temp_bed_fasta_dir)

            #lib_path = ''  since we will change the current dir to eachfm_dir
            lib_path = ''
            if eachfm == 'LTR':
                lib_path = './opt_Gypsy_Copia.lib'
            if eachfm == 'DNA':
                lib_path = './opt_MITE.lib'

            #else:
            #    lib_path = './final.lib'

            ##running pipeline

            ##go to the target dir
            ##get the current dir
            cur_dirpath = os.getcwd()

            os.chdir(eachfm_dir)

            cmd = 'python3.7 ./pipeline_generate_gold_reference_TE.py ' + \
                  lib_path + ' ' + raw_data_dir + '/' + clean_noTE_fl + ' ' + search_name + \
                  ' ' + software_simulome_path + ' ' + './temp_bed_fasta_dir' + ' ' + str(replicate_count) + ' ' + \
                    str(total_replicate_count)
            print(cmd)
            subprocess.call(cmd,shell=True)

            ##go to the originial path
            os.chdir(cur_dirpath)


            ####################
            ##step 2: evaluation
            ##generate eachfm dir
            eval_eachfm_dir = evaluation_dir + '/' + eachfm
            if not os.path.exists(eval_eachfm_dir):
                os.makedirs(eval_eachfm_dir)

            ##cp the required files to eachfm dir
            cmd = 'cp ' + raw_simulate_required_fl_dir + '/batch_running_evaluation.py ' + \
                  raw_simulate_required_fl_dir + '/batch_running_landscape.py ' + \
                  raw_simulate_required_fl_dir + '/batch_combine_result.py ' + eval_eachfm_dir + '/'
            print(cmd)
            subprocess.call(cmd, shell=True)

            cmd = 'cp -r ' + raw_simulate_required_fl_dir + '/TEScape ' + eval_eachfm_dir + '/'
            print(cmd)
            subprocess.call(cmd, shell=True)

            ##cp the output from the simulate_seq to the current file
            cmd = 'cp ' + eachfm_dir + '/opt_final_te.fasta ' + eval_eachfm_dir + '/opt_final.fasta'
            print(cmd)
            subprocess.call(cmd, shell=True)
            cmd = 'cp ' + eachfm_dir + '/opt_te_location_infor.txt ' + eval_eachfm_dir + '/opt_insert_te_location.txt'
            print(cmd)
            subprocess.call(cmd, shell=True)

            ##make required dir
            accuracy_dir = eval_eachfm_dir + '/accuracy'
            if not os.path.exists(accuracy_dir):
                os.makedirs(accuracy_dir)

            landscape_dir = eval_eachfm_dir + '/landscape'
            if not os.path.exists(landscape_dir):
                os.makedirs(landscape_dir)

            lib_construct_dir = eval_eachfm_dir + '/lib_construction'
            if not os.path.exists(lib_construct_dir):
                os.makedirs(lib_construct_dir)

            raw_evaluation_dir = eval_eachfm_dir + '/raw_evaluation_dir'
            if not os.path.exists(raw_evaluation_dir):
                os.makedirs(raw_evaluation_dir)

            raw_input_lib_dir = eval_eachfm_dir + '/raw_input_lib_dir'
            if not os.path.exists(raw_input_lib_dir):
                os.makedirs(raw_input_lib_dir)

            raw_landscape_dir = eval_eachfm_dir + '/raw_landscape_dir'
            if not os.path.exists(raw_landscape_dir):
                os.makedirs(raw_landscape_dir)

            ##a) conduct lib construct analysis
            ##construct output and working_dir
            lib_working_ltf_dir = ''
            lib_output_ltf_dir = ''
            lib_working_har_dir = ''
            lib_output_har_dir = ''
            if eachfm == 'LTR':
                lib_working_dir = lib_construct_dir + '/working_dir'
                if not os.path.exists(lib_working_dir):
                    os.makedirs(lib_working_dir)

                lib_output_dir = lib_construct_dir + '/output_dir'
                if not os.path.exists(lib_output_dir):
                    os.makedirs(lib_output_dir)

                lib_working_ltf_dir = lib_construct_dir + '/working_ltf_dir'
                if not os.path.exists(lib_working_ltf_dir):
                    os.makedirs(lib_working_ltf_dir)

                lib_output_ltf_dir = lib_construct_dir + '/output_ltf_dir'
                if not os.path.exists(lib_output_ltf_dir):
                    os.makedirs(lib_output_ltf_dir)

                lib_working_har_dir = lib_construct_dir + '/working_har_dir'
                if not os.path.exists(lib_working_har_dir):
                    os.makedirs(lib_working_har_dir)

                lib_output_har_dir = lib_construct_dir + '/output_har_dir'
                if not os.path.exists(lib_output_har_dir):
                    os.makedirs(lib_output_har_dir)

            else:
                lib_working_dir = lib_construct_dir + '/working_dir'
                if not os.path.exists(lib_working_dir):
                    os.makedirs(lib_working_dir)

                lib_output_dir = lib_construct_dir + '/output_dir'
                if not os.path.exists(lib_output_dir):
                    os.makedirs(lib_output_dir)

            ##generate a software path
            ##different family need different software path
            lib_software_path = ''
            ##updation 032720
            lib_software_path_ltf = ''
            lib_software_path_har = ''

            if eachfm == 'LTR':
                lib_software_path = '--ltr_retriever ~/software/LTR_retriever/LTR_retriever --ltrfinder ' \
                                    '/data/haidongy_projects/software/LTR_Finder/source/ltr_finder ' \
                                    '--gt /data/haidongy_projects/software/bin/gt/bin/gt ' \
                                    '-para ' + eval_eachfm_dir + '/TEScape/test_input_dir/ltr_parameter.txt ' \
                                                                 '--repeatmodeler_dir /data/software/RepeatModeler --repmd_pro 10 --lrv_pro 15'

                ##updation 032720 only provide the ltr_finder
                ##only call ltrfinder
                lib_software_path_ltf = '--ltr_retriever ~/software/LTR_retriever/LTR_retriever --ltrfinder ' \
                                    '/data/haidongy_projects/software/LTR_Finder/source/ltr_finder ' \
                                    '-para ' + eval_eachfm_dir + '/TEScape/test_input_dir/ltr_parameter.txt --lrv_pro 15'

                ##only call ltrharvest
                lib_software_path_har = '--ltr_retriever ~/software/LTR_retriever/LTR_retriever ' \
                                    '--gt /data/haidongy_projects/software/bin/gt/bin/gt ' \
                                    '-para ' + eval_eachfm_dir + '/TEScape/test_input_dir/ltr_parameter.txt --lrv_pro 15'


            if eachfm == 'DNA':
                lib_software_path = '--repeatmodeler_dir /data/software/RepeatModeler --repmd_pro 10 --mite_finder2 ~/software/miteFinder/bin/miteFinder ' \
                                    '--mite_finder2_profile ~/software/miteFinder/profile/pattern_scoring.txt ' \
                                    '--mf2_thr 0.5 ' \
                                    '--mite_hunter_dir /data/haidongy_projects/software/MITE_Hunter/' ##only provide MITE_Hunter dir

            if eachfm == 'LINE' or eachfm_dir == 'SINE':
                lib_software_path = '--repeatmodeler_dir /data/software/RepeatModeler --repmd_pro 10'

            ##running lib_construction
            ##updation 032720
            if eachfm == 'LTR':
                cmd = 'python3.7 ' + eval_eachfm_dir + '/TEScape/TEScape_lib.py -fas ' + eval_eachfm_dir + '/opt_final.fasta ' + \
                      '-d ' + lib_working_dir + ' -o ' + lib_output_dir + ' ' + lib_software_path
                print(cmd)
                subprocess.call(cmd, shell=True)

                cmd = 'python3.7 ' + eval_eachfm_dir + '/TEScape/TEScape_lib.py -fas ' + eval_eachfm_dir + '/opt_final.fasta ' + \
                      '-d ' + lib_working_ltf_dir + ' -o ' + lib_output_ltf_dir + ' ' + lib_software_path_ltf
                print(cmd)
                subprocess.call(cmd, shell=True)

                cmd = 'python3.7 ' + eval_eachfm_dir + '/TEScape/TEScape_lib.py -fas ' + eval_eachfm_dir + '/opt_final.fasta ' + \
                      '-d ' + lib_working_har_dir + ' -o ' + lib_output_har_dir + ' ' + lib_software_path_har
                print(cmd)
                subprocess.call(cmd, shell=True)

            else:
                cmd = 'python3.7 ' + eval_eachfm_dir + '/TEScape/TEScape_lib.py -fas ' + eval_eachfm_dir + '/opt_final.fasta ' + \
                      '-d ' + lib_working_dir + ' -o ' + lib_output_dir + ' ' + lib_software_path
                print(cmd)
                subprocess.call(cmd, shell=True)



            ##updation 032720

            if eachfm == 'LTR':
                ##updation 032720
                ##add LTR information
                cmd = 'cp ' + lib_working_ltf_dir + '/combine/ltr_retriever.lib ' + raw_input_lib_dir + '/ltr_finder.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

                cmd = 'cp ' + lib_working_har_dir + '/combine/ltr_retriever.lib ' + raw_input_lib_dir + '/ltr_harvest.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

                ##cp the output from the lib to the raw_input_lib_dir
                cmd = 'cp ' + lib_working_dir + '/combine/repeatmodeler.lib ' + raw_input_lib_dir + '/rpmd.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

                cmd = 'cp ' + lib_working_dir + '/combine/ltr_retriever.lib ' + raw_input_lib_dir + '/ltr_retriever.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

                cmd = 'cp ' + eachfm_dir + '/opt_temp_20.lib ' + raw_input_lib_dir + '/ref.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

                ##generate combine information
                cmd = 'cat ' + raw_input_lib_dir + '/ref.lib ' + raw_input_lib_dir + '/rpmd.lib ' + \
                      raw_input_lib_dir + '/ltr_retriever.lib > ' + raw_input_lib_dir + '/ref_rpmd_ltr_retriever.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

                cmd = 'cat ' + raw_input_lib_dir + '/ref.lib ' + raw_input_lib_dir + '/rpmd.lib > ' + \
                      raw_input_lib_dir + '/ref_rpmd.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

                cmd = 'cat ' + raw_input_lib_dir + '/ref.lib ' + raw_input_lib_dir + '/ltr_retriever.lib > ' + \
                      raw_input_lib_dir + '/ref_ltr_retriever.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

                ##update.6.7 add the other combination
                cmd = 'cat ' + raw_input_lib_dir + '/rpmd.lib ' + raw_input_lib_dir + '/ltr_retriever.lib > ' + \
                      raw_input_lib_dir + '/rpmd_ltr_retriever.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

                ##updation 032720 add ltrfinder and ltrharvest to the pipeline
                cmd = 'cat ' + raw_input_lib_dir + '/ref.lib ' + raw_input_lib_dir + '/ltr_finder.lib > ' + \
                      raw_input_lib_dir + '/ref_ltr_finder.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

                cmd = 'cat ' + raw_input_lib_dir + '/ref.lib ' + raw_input_lib_dir + '/ltr_harvest.lib > ' + \
                      raw_input_lib_dir + '/ref_ltr_harvest.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

                cmd = 'cat ' + raw_input_lib_dir + '/rpmd.lib ' + raw_input_lib_dir + '/ltr_finder.lib > ' + \
                      raw_input_lib_dir + '/rpmd_ltr_finder.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

                cmd = 'cat ' + raw_input_lib_dir + '/rpmd.lib ' + raw_input_lib_dir + '/ltr_harvest.lib > ' + \
                      raw_input_lib_dir + '/rpmd_ltr_harvest.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

                cmd = 'cat ' + raw_input_lib_dir + '/ref.lib ' + raw_input_lib_dir + '/rpmd.lib ' + \
                      raw_input_lib_dir + '/ltr_finder.lib > ' + raw_input_lib_dir + '/ref_rpmd_ltr_finder.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

                cmd = 'cat ' + raw_input_lib_dir + '/ref.lib ' + raw_input_lib_dir + '/rpmd.lib ' + \
                      raw_input_lib_dir + '/ltr_harvest.lib > ' + raw_input_lib_dir + '/ref_rpmd_ltr_harvest.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)


            if eachfm == 'LINE' or eachfm == 'SINE':
                cmd = 'cp ' + eachfm_dir + '/opt_temp_20.lib ' + raw_input_lib_dir + '/ref.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

                cmd = 'cp ' + lib_working_dir + '/combine/repeatmodeler.lib ' + raw_input_lib_dir + '/rpmd.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

                ##generate combine information
                cmd = 'cat ' + raw_input_lib_dir + '/ref.lib ' + raw_input_lib_dir + '/rpmd.lib > ' + \
                      raw_input_lib_dir + '/ref_rpmd.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

            if eachfm == 'DNA':
                cmd = 'cp ' + eachfm_dir + '/opt_temp_20.lib ' + raw_input_lib_dir + '/ref.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

                cmd = 'cp ' + lib_working_dir + '/combine/repeatmodeler.lib ' + raw_input_lib_dir + '/rpmd.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

                ##add other siutations
                cmd = 'cp ' + lib_working_dir + '/combine/mite_finder2.lib ' + raw_input_lib_dir + '/mite_finder2.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

                cmd = 'cp ' + lib_working_dir + '/combine/mite_hunter.lib ' + raw_input_lib_dir + '/mite_hunter.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

                ##generate combine information
                cmd = 'cat ' + raw_input_lib_dir + '/ref.lib ' + raw_input_lib_dir + '/rpmd.lib > ' + \
                      raw_input_lib_dir + '/ref_rpmd.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

                cmd = 'cat ' + raw_input_lib_dir + '/ref.lib ' + raw_input_lib_dir + '/mite_finder2.lib > ' + \
                      raw_input_lib_dir + '/ref_mite_finder2.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

                cmd = 'cat ' + raw_input_lib_dir + '/ref.lib ' + raw_input_lib_dir + '/mite_hunter.lib > ' + \
                      raw_input_lib_dir + '/ref_mite_hunter.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

                cmd = 'cat ' + raw_input_lib_dir + '/ref.lib ' + raw_input_lib_dir + '/mite_hunter.lib ' + \
                      raw_input_lib_dir + '/mite_finder2.lib ' + raw_input_lib_dir + '/rpmd.lib > ' + \
                      raw_input_lib_dir + '/ref_mite_hunter_mite_finder2_rpmd.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

                ##update6.7 add other combinations
                cmd = 'cat ' + raw_input_lib_dir + '/rpmd.lib ' + raw_input_lib_dir + '/mite_hunter.lib > ' + \
                      raw_input_lib_dir + '/rpmd_mite_hunter.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

                cmd = 'cat ' + raw_input_lib_dir + '/rpmd.lib ' + raw_input_lib_dir + '/mite_finder2.lib > ' + \
                      raw_input_lib_dir + '/rpmd_mite_finder2.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)

                cmd = 'cat ' + raw_input_lib_dir + '/rpmd.lib ' + raw_input_lib_dir + '/mite_hunter.lib ' + \
                      raw_input_lib_dir + '/mite_finder2.lib > ' + raw_input_lib_dir + '/rpmd_mite_finder2_mite_hunter.lib'
                print(cmd)
                subprocess.call(cmd, shell=True)


            ##b) conduct landscape analysis
            ##cp files to the raw_landscape_dir
            cmd = 'cp -r ' + eval_eachfm_dir + '/opt_final.fasta ' + eval_eachfm_dir + '/TEScape ' + raw_landscape_dir + '/'
            print(cmd)
            subprocess.call(cmd, shell=True)

            cmd = 'python3.7 ' + eval_eachfm_dir + '/batch_running_landscape.py ' + raw_landscape_dir + ' ' + landscape_dir + ' ' + raw_input_lib_dir + \
                  ' ' + raw_landscape_dir + '/opt_final.fasta'
            print(cmd)
            subprocess.call(cmd, shell=True)

            ##c) conduct evaluation analysis
            ##cp files to the raw_evaluation_dir
            evalution_python_path = ''
            if eachfm == 'LTR':
                evalution_python_path = raw_simulate_required_fl_dir + '/evaluation_LTR.py'
            if eachfm == 'DNA':
                evalution_python_path = raw_simulate_required_fl_dir + '/evaluation_dna_mite.py'
            if eachfm == 'SINE':
                evalution_python_path = raw_simulate_required_fl_dir + '/evaluation_SINE.py'
            if eachfm == 'LINE':
                evalution_python_path = raw_simulate_required_fl_dir + '/evaluation_LINE.py'

            cmd = 'cp ' + eval_eachfm_dir + '/opt_final.fasta ' + eval_eachfm_dir + '/opt_insert_te_location.txt ' + \
                  raw_evaluation_dir + '/'
            print(cmd)
            subprocess.call(cmd, shell=True)


            ##run the batch evaluation
            cmd = 'python3.7 ' + eval_eachfm_dir + '/batch_running_evaluation.py ' + landscape_dir + ' ' + raw_evaluation_dir + ' ' + \
                  accuracy_dir + ' ' + evalution_python_path
            print(cmd)
            subprocess.call(cmd, shell=True)


            ##d) combine evaluation result
            cmd = 'python3.7 ' + eval_eachfm_dir + '/batch_combine_result.py ' + accuracy_dir + ' ' + raw_simulate_required_fl_dir + '/type_list.txt ' + \
                  eval_eachfm_dir
            print(cmd)
            subprocess.call(cmd, shell=True)

            ##########################
            ##store results into output
            output_path = ''
            if eachfm == 'LTR':
                output_path = LTR_output_dir + '/opt_LTR_' + str(i) + '.txt'
            if eachfm == 'DNA':
                output_path = DNA_output_dir + '/opt_DNA_' + str(i) + '.txt'
            if eachfm == 'LINE':
                output_path = LINE_output_dir + '/opt_LINE_' + str(i) + '.txt'
            if eachfm == 'SINE':
                output_path = SINE_output_dir + '/opt_SINE_' + str(i) + '.txt'

            cmd = 'cp ' + eval_eachfm_dir + '/opt_final_combine_all_genus.txt ' + output_path
            print(cmd)
            subprocess.call(cmd, shell=True)


##run the simulation
mul_simulate (replicate_count,total_replicate_count,input_simulate_time,store_simulate_results_dir,raw_simulate_required_fl_dir,software_simulome_path,all_simulate_output_dir,analyze_family_file)
