#!/usr/bin/env python3.7

##updation032920 add the threads for the ltr_retriever

##update: 5/13/2019 change the current dir for the mitehunter
##update: 5/2/2019
##add the mitehunter and mitefinder2
##update: 22/1/2019
##content: set the configuration.rc to the ltr_parameter.txt and set it as the argument that users can revise that

##this wrapper will help to generate the TE libraries if users do not provide the libraries
##if the users want to change the parameters for each recommended software, they can directly revise the parameters in
##the ltr_parameter.txt file in the test_input_dir.

##BUILT-IN MODULES
import argparse
import sys
from distutils.spawn import find_executable

import os
import subprocess
import glob
import re

##SCRIPTS
def get_parsed_args():
    parser = argparse.ArgumentParser(description="Generate TE lib information")

    ##############################
    ##parse the required arguments
    parser.add_argument("-d", dest='working_dir', default="./", help="Working directory to store the library generated"
                                                                     " by the pipline recommanded software. "
                                                                     "Default: ./ ")

    parser.add_argument('-fas', dest='genome_fas', help="Reference genome sequence information (fasta format)")

    parser.add_argument("-o", dest='output_dir', default="./", help="output directory to store the combined library. "
                                                                    "Default: ./ ")


    ####################
    ##parse of softwares
    parser.add_argument("--repeatmodeler_dir",dest='repeatmodeler_dir', help="File path to repeatmodeler dir"
                                                                             "This dir should contains BuildDatabase and RepeatModeler")
    ##optional processer number
    parser.add_argument('--repmd_pro', dest='repmd_pro',help='processor number in the repeatmodeler')

    parser.add_argument('--lrv_pro', dest='lrv_pro',help='processor number in the ltr_retriever')




    ##change the repeatmasker dir to the queryRepeatDatabase.pl
    parser.add_argument("--repeatmasker_dir",dest='repeatmasker_dir', help="File path to repeatmasker dir"
                                                                           "This dir should contains util/queryRepeatDatabase.pl"
                                                                           "If set this argument, users should provide species_list"
                                                                           "The species_list contains information about the species name and genus name"
                                                                           "If the species/genus is not exist in the database, the error report will be"
                                                                           "stored in the output dir")

    parser.add_argument("--ltr_retriever", dest="ltr_retriever", help="File path to LTR_retriever executable"
                                                                      "If call this argument, executable from ltrfinder or gt" 
                                                                      "(GenomeTools) should be provided")

    parser.add_argument("--ltrfinder", dest="ltr_finder",help="File path to ltrfinder executable"                                                          
                                                                "Parameters can be found in ltr_parameter.txt in test_input_dir")

    parser.add_argument("--gt", dest="gt",help="File path to ltrfinder executable "
                                                "ltr_harvest that will be used is wrapped into genometools"                                                      
                                                "Parameters can be found in ltr_parameter.txt in test_input_dir")

    ##updation5.2 add mitehunter and mitefinder2
    ##add mitefinder2
    parser.add_argument('--mite_finder2',dest='mite_finder2',help='File path to mitefinder2 executable')
    parser.add_argument('--mite_finder2_profile',dest='mite_finder2_profile',help='File path to mitefinder2 profile executable,(pattern_scoring.txt)')
    ##optional threshold
    parser.add_argument('--mf2_thr',dest='mf2_thr',help='threshold in the mitefinder2 parameter')

    ##add mitehunter
    parser.add_argument('--mite_hunter_dir',dest='mite_hunter_dir',help='File path to mitehunter dir')
    ##if provided mitehunter
    ##please make sure the
    ##formatdb
    ##blastall
    ##mdust
    ##muscle
    ##are environmental variable
    ##add install indication
    parser.add_argument('--mite_hunter_install',dest='install_yes',help='If specify the install, mitehunter will install')



    ##############################
    ##parse of optional te library
    parser.add_argument("-c", dest='other_lib',help="Provide other te libraries that can be combined into te libraries "
                                                    "generated from the recommended software")

    parser.add_argument("-p", dest='species_list', help="Provide the list of species from repbase as reference libraries"
                                                        "the species list contains species/genus name")

    parser.add_argument("-para", dest='ltr_parameter', help="Provide the ltr parameters when running ltr_retriever"
                                                            "the example of the parameter file is in the test_input_dir"
                                                            "Users can directly revise in this file")




    ##parse of parameters
    args = parser.parse_args()
    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv
    args = get_parsed_args()


    ##check for the input files (required)
    if args.genome_fas is None:
        print ('Cannot find input genome fasta file, please provide that')
        return
    else:
        try:
            file = open(args.genome_fas, 'r')
        except IOError:
            print('There was an error opening the genome file!')
            return


    ##check for other lib
    if args.other_lib is not None:
        try:
            file = open(args.other_lib, 'r')
        except IOError:
            print('There was an error opening the file!')
            return


    ##check for repeatmodeler
    if args.repeatmodeler_dir is not None:

        ##check the builddatabase
        repeatmodeler_dir = args.repeatmodeler_dir
        if find_executable(repeatmodeler_dir + '/BuildDatabase') is not None:
            print('repeatmodeler BuildDatabase can be found')
        else:
            print("Cannot find BuildDatabase executable, please check if it has been included in repeatmodeler.")
            return

        ##check the repeatmodeler
        if find_executable(repeatmodeler_dir + '/RepeatModeler') is not None:
            print('repeatmodeler can be found')
        else:
            print("Cannot find RepeatModeler executable, please check if it has been included in repeatmodeler.")
            return


    ##check for repeatmasker
    if args.repeatmasker_dir is not None:

        ##check the builddatabase
        repeatmasker_dir = args.repeatmasker_dir

        ##check the util/queryRepeatDatabase.pl
        if find_executable(repeatmasker_dir + '/util/queryRepeatDatabase.pl') is not None:
            print('repeatmasker queryRepeatDatabase can be found')
        else:
            print("Cannot find queryRepeatDatabase executable, please check if it has been included in repeatmasker.")
            return

        ##check the species list
        if args.species_list is None:
            print ('users should provide species_list after calling the repeatmasker argument')
            return
        else:
            try:
                file = open(args.species_list, 'r')
            except IOError:
                print('There was an error opening the species file!')
                return


    ##check for ltr_retriever
    if args.ltr_retriever is not None:


        ##check for the ltr_retriever
        if find_executable(args.ltr_retriever) is not None:
            print('ltr_retriever executable can be found')
        else:
            print("Cannot find ltr_retriever executable, please check if it has been installed.")
            return

        if args.ltr_finder is None and args.gt is None:
            print ('ltrfinder and ltrharvest should be provided for single each or both')
            return

        ##check for the ltrfinder
        if args.ltr_finder is not None:
            if find_executable(args.ltr_finder) is not None:
                print('ltrfinder executable can be found')
            else:
                print("Cannot find ltrfinder executable, please check if it has been installed.")
                return

        ##check for the gt
        if args.gt is not None:
            if find_executable(args.gt) is not None:
                print('gt executable can be found')
            else:
                print("Cannot find gt executable, please check if it has been installed.")
                return

        ##check for the parameter files
        if args.ltr_parameter is None:
            print ('users should provide ltr parameters after calling the ltr_retriever argument')
            return
        else:
            try:
                file = open(args.ltr_parameter, 'r')
            except IOError:
                print('There was an error opening the file!')
                return


    ##check for the ltrfinder
    if args.ltr_finder is not None:
        if find_executable(args.ltr_finder) is not None:
            print('ltrfinder executable can be found')
        else:
            print("Cannot find ltrfinder executable, please check if it has been installed.")
            return

        if args.ltr_retriever is None:
            print('ltr_retriever software need to be provided when ltrfinder is called')
            return


    ##check for the gt
    if args.gt is not None:
        if find_executable(args.gt) is not None:
            print('gt executable can be found')
        else:
            print("Cannot find gt executable, please check if it has been installed.")
            return

        if args.ltr_retriever is None:
            print('ltr_retriever software need to be provided when gt is called')
            return

    ##check for the mitefinder2
    if args.mite_finder2 is not None:

        if find_executable(args.mite_finder2) is not None:
            print('mite_finder2 executable can be found')
        else:
            print("Cannot find mite_finder2 executable, please check if it has been installed.")
            return

        ##check the profile file in the mitefinder2
        if args.mite_finder2_profile is None:
            print('users should provide mite_finder2_profile after calling the mite_finder2 argument')
            return
        else:
            try:
                file = open(args.mite_finder2_profile, 'r')
            except IOError:
                print('There was an error opening the mite_finder2_profile file!')
                return

    ##check for the mitehunter
    if args.mite_hunter_dir is not None:

        ##check the MITE_Hunter_manager
        mite_hunter_dir = args.mite_hunter_dir
        if find_executable(mite_hunter_dir + '/MITE_Hunter_manager.pl') is not None:
            print('MITE_Hunter_manager.pl can be found')
        else:
            print("Cannot find MITE_Hunter_manager.pl executable, please check if it has been included in mite_hunter dir.")
            return

        ##check the MITE_Hunter_Installer
        if find_executable(mite_hunter_dir + '/MITE_Hunter_Installer.pl') is not None:
            print('MITE_Hunter_Installer.pl can be found')
        else:
            print("Cannot find MITE_Hunter_Installer.pl executable, please check if it has been included in mite_hunter dir.")
            return



    ###########################################
    ##create the working and output directories
    working_dir = args.working_dir
    if not working_dir.endswith('/'):
        working_dir = working_dir + '/'
    else:
        working_dir = working_dir

    output_dir = args.output_dir
    if not output_dir.endswith('/'):
        output_dir = output_dir + '/'
    else:
        output_dir = output_dir

    ##########################################
    ##cp the genome into the working dir
    genome_file = args.genome_fas
    cmd = 'cp ' + genome_file + ' ' + working_dir + '/genome.fas'
    subprocess.call(cmd, shell=True)

    ##create a combine dir to store all the libraries
    combine_lib_dir = working_dir + '/combine/'
    if not os.path.exists(combine_lib_dir):  ##if dir not exit, it will create the dir
        os.makedirs(combine_lib_dir)



    ##############################################
    ##if the users type the repeatmodeler argument
    ##############################################
    if args.repeatmodeler_dir is not None:
        ##create a dir for the repeatmodeler
        repeatmodeler_dir = args.repeatmodeler_dir
        if not repeatmodeler_dir.endswith('/'):
            repeatmodeler_dir = repeatmodeler_dir + '/'
        else:
            repeatmodeler_dir = repeatmodeler_dir

        ##create the dir for the repeatmodeler
        repeatmodeler_w_dir = working_dir + '/repeatmodeler/'
        if not os.path.exists(repeatmodeler_w_dir):  ##if dir not exit, it will create the dir
            os.makedirs(repeatmodeler_w_dir)

        cmd = repeatmodeler_dir + '/BuildDatabase -name rpm_genome ' + working_dir + '/genome.fas '
        subprocess.call(cmd, shell=True)

        ##call the repeatmodeler
        ##decide the processor number
        if args.repmd_pro is not None:
            repmd_pro_num = args.repmd_pro
        else:
            repmd_pro_num = '1'

        cmd = repeatmodeler_dir + '/RepeatModeler' + ' -pa ' + repmd_pro_num + ' -database rpm_genome'
        print(cmd)
        subprocess.call(cmd, shell=True)

        ##cp the output of the repeatmodeler to the combine dir
        cmd = 'cp ./RM*/consensi.fa.classified ' + combine_lib_dir + '/repeatmodeler.lib'
        subprocess.call(cmd, shell=True)

        ##move the temp file to the repeatmodeler working dir
        cmd = 'mv ./RM* rpm* tmpConsensi.fa.cat.all unaligned.fa ' + repeatmodeler_w_dir
        subprocess.call(cmd, shell=True)


    #############################################
    ##if the users type the repeatmasker argument
    #############################################
    if args.repeatmasker_dir is not None:
        repeatmasker_dir = args.repeatmasker_dir
        if not repeatmasker_dir.endswith('/'):
            repeatmasker_dir = repeatmasker_dir + '/'
        else:
            repeatmasker_dir = repeatmasker_dir

        ##create the dir for the repeatmasker
        repeatmasker_w_dir = working_dir + 'repeatmasker/'
        if not os.path.exists(repeatmasker_w_dir):  ##if dir not exit, it will create the dir
            os.makedirs(repeatmasker_w_dir)

        ##store the name information from the species list
        species_file = args.species_list
        species_dic = {}
        with open (species_file, 'r') as ipt_sp:
            for eachline in ipt_sp:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                species_dic[col[0]] = 1


        ##create string to store TE
        TE_tp = 'LTR;LINE;SINE;DNA'
        col_TE_tp = TE_tp.split(';')
        ##create two dirs including right and error dirs
        lib_right_dir = repeatmasker_w_dir + '/lib_right_dir/'
        if not os.path.exists(lib_right_dir):  ##if dir not exit, it will create the dir
            os.makedirs(lib_right_dir)

        ##the following libs only appear lib_right_TEs_dir, lib_right_all_dir, lib_error_dir
        ##for the right dir, it classify the all lib including non TEs and classified dir including TEs solely
        lib_right_TEs_dir = lib_right_dir + '/lib_right_TEs_dir/'
        if not os.path.exists(lib_right_TEs_dir):  ##if dir not exit, it will create the dir
            os.makedirs(lib_right_TEs_dir)
        lib_right_all_dir = lib_right_dir + '/lib_right_all_dir/'
        if not os.path.exists(lib_right_all_dir):  ##if dir not exit, it will create the dir
            os.makedirs(lib_right_all_dir)
        ##for the error dir
        #lib_error_dir = repeatmasker_w_dir + '/lib_error_dir/'
        #if not os.path.exists(lib_error_dir):  ##if dir not exit, it will create the dir
        #    os.makedirs(lib_error_dir)

        ##initiate a dic to store the right and error dir
        right_dic = {}

        ##call the repbase for each name
        for eachnm in species_dic:
            print ('the analyzed name is ' + eachnm)
            ##run the query repbase function
            #cmd = repeatmasker_dir + 'util/queryRepeatDatabase.pl -species ' + eachnm + \
            #      ' &> ' + repeatmasker_w_dir + new_lib_nm
            #print ('the cmd for repeatmasker is ' + cmd)

            cmd = repeatmasker_dir + 'util/queryRepeatDatabase.pl -species ' + eachnm
            subprocess.call(cmd, stdout=open(eachnm + '.out','w'),stderr=open(eachnm + '.err','w'),shell=True)

            ##mv the out and err to the working dir
            cmd = 'mv ' + eachnm + '.out ' +  eachnm + '.err ' + repeatmasker_w_dir
            print (cmd)
            subprocess.call(cmd, shell=True)


            err_count = 0
            ##check the out lib and err lib to decide whether the species name is right or not
            with open (repeatmasker_w_dir + eachnm + '.err','r') as ipt_err:
                for eachline in ipt_err:
                    eachline = eachline.strip('\n')
                    if 'not in the database' in eachline:
                        err_count += 1

            ##if the typing name is wrong name
            if err_count == 1:
                ##move the lib to the error dir
                cmd = 'mv ' + repeatmasker_w_dir + eachnm + '.err ' + output_dir
                print("the wrong detection is " + cmd)
                subprocess.call(cmd, shell=True)
                ##delete the out, because out is empty
                cmd = 'rm ' + repeatmasker_w_dir + eachnm + '.out'
                print("the wrong detection is " + cmd)
                subprocess.call(cmd, shell=True)
            ##if the typing name is right
            else:
                ##move the lib to the right dir, and change the out name to the lib name
                cmd = 'mv ' + repeatmasker_w_dir + eachnm + '.out ' + lib_right_all_dir + eachnm + '.lib'
                subprocess.call(cmd, shell=True)
                print("the right detection is " + cmd)
                ##rm the err version
                cmd = 'rm ' + repeatmasker_w_dir + eachnm + '.err'
                subprocess.call(cmd, shell=True)
                print("the right detection is " + cmd)
                ##store the right name into a dic
                right_dic[eachnm] = 1


        ##call the repbase for each name from the right dir
        for eachnm in right_dic:
            ##each nm in right_dic will filter the TE: LTR/LINE/SINE/DNA
            ##for each type of TE
            for eachtp in col_TE_tp:
                ##generate the new name
                new_lib_nm = eachnm + '_' + eachtp + '.lib'
                ##run the repbase query
                cmd = repeatmasker_dir + 'util/queryRepeatDatabase.pl -species ' + eachnm + ' -class ' + eachtp + \
                      ' > ' + lib_right_TEs_dir + new_lib_nm
                subprocess.call(cmd, shell=True)

        ##combine each lib in the lib_right_TEs_dir to the combine TE
        te_lib_list = ''
        te_lib_file_list = glob.glob(lib_right_TEs_dir + '*.lib')
        for eachlib in te_lib_file_list:
            te_lib_list = eachlib + ' ' + te_lib_list
        ##combination
        cmd = 'cat ' + te_lib_list + ' > ' + combine_lib_dir + '/repeatmasker.lib'
        subprocess.call(cmd, shell=True)



    #####################################
    ##if the users type the ltr_retriever
    #####################################
    if args.ltr_retriever is not None:

        ##store the parameter information on a dic
        para_dic = {}

        para_file = args.ltr_parameter

        with open(para_file, 'r') as ipt_rc:
            for eachline in ipt_rc:
                eachline = eachline.strip('\n')
                # print ('eachline is ' + eachline)
                mt = re.match('(.+):(.+)', eachline)
                softwarenm = mt.group(1)
                para = mt.group(2)
                para_dic[softwarenm] = para

        ##get retriever executable file
        ltr_retriever_exe = args.ltr_retriever
        ##create a dir for the ltr retriever
        ltr_retriever_w_dir = working_dir + '/ltr_retriever/'
        if not os.path.exists(ltr_retriever_w_dir):
            os.makedirs(ltr_retriever_w_dir)

        ##set the threads
        if args.lrv_pro is not None:
            lrv_pro_num = args.lrv_pro
        else:
            lrv_pro_num = '1'


        ##if ltr_finder is provided
        if args.ltr_finder is not None:
            ltrfinder_exe = args.ltr_finder
            ##create a dir for the ltr finder
            ltr_finder_w_dir = working_dir + '/ltr_finder/'
            if not os.path.exists(ltr_finder_w_dir):
                os.makedirs(ltr_finder_w_dir)
            ##run the ltrfinder
            ##extract the parameters from the para_file file
            para = ''
            for eachkey in para_dic:
                if eachkey == 'ltrfinder':
                    para = para_dic[eachkey]
            ##call ltrfinder
            cmd = ltrfinder_exe + ' ' + para + ' ' + working_dir + '/genome.fas > ' + ltr_finder_w_dir + '/genome.finder.scn'
            subprocess.call(cmd, shell=True)
            ##cp the results to the ltr retriever
            cmd = 'cp ' + ltr_finder_w_dir + '/genome.finder.scn ' + ltr_retriever_w_dir
            subprocess.call(cmd, shell=True)


        ##if gt is provided
        if args.gt is not None:
            gt_exe = args.gt
            ##create a dir for the ltr harvest
            ltr_harvest_w_dir = working_dir + '/ltr_harvest/'
            if not os.path.exists(ltr_harvest_w_dir):
                os.makedirs(ltr_harvest_w_dir)
            ##run the ltrharvest
            ##create the index for the genome
            cmd = gt_exe + ' suffixerator -db ' + working_dir + '/genome.fas -indexname ' + working_dir + '/gt_genome' + ' -tis -suf -lcp -des -ssp -sds -dna'
            subprocess.call(cmd, shell=True)
            ##extract the parameters from the para_file
            para = ''
            for eachkey in para_dic:
                if eachkey == 'ltrharvest':
                    para = para_dic[eachkey]
            ##call ltrharvest
            cmd = gt_exe + ' ltrharvest -index ' + working_dir + '/gt_genome ' + para + ' > ' + ltr_harvest_w_dir + '/genome.harvest.scn'
            subprocess.call(cmd, shell=True)
            ##cp the results to the ltr retriever
            cmd = 'cp ' + ltr_harvest_w_dir + '/genome.harvest.scn ' + ltr_retriever_w_dir
            subprocess.call(cmd, shell=True)


        ##The output of ltr_retriever is in the dir same as the genome.fas
        ##So the output of the ltr_retriever is in the working dir name *LTRlib.fa
        ##call the ltr_retriever
        ##updation 032920 add threads
        if args.ltr_finder is not None and args.gt is not None:
            cmd = ltr_retriever_exe + ' -genome ' + working_dir + '/genome.fas -inharvest ' + ltr_retriever_w_dir + \
                  '/genome.harvest.scn -infinder ' +  ltr_retriever_w_dir + '/genome.finder.scn ' + '-threads ' + lrv_pro_num
            subprocess.call(cmd, shell=True)

        if args.ltr_finder is not None and args.gt is None:
            cmd = ltr_retriever_exe + ' -genome ' + working_dir + '/genome.fas -infinder ' + ltr_retriever_w_dir + \
                  '/genome.finder.scn ' + '-threads ' + lrv_pro_num
            subprocess.call(cmd, shell=True)

        if args.gt is not None and args.ltr_finder is None:
            cmd = ltr_retriever_exe + ' -genome ' + working_dir + '/genome.fas -inharvest ' + ltr_retriever_w_dir + \
                  '/genome.harvest.scn ' + '-threads ' + lrv_pro_num
            subprocess.call(cmd, shell=True)

        ##cp the output from the ltr_retriever to the combine_lib_dir
        all_work_file_list = glob.glob(working_dir + '/*.LTRlib.fa')
        for eachfl in all_work_file_list:
            if 'nmtf' not in eachfl:
                cmd = 'cp ' + eachfl + ' ' + combine_lib_dir + '/ltr_retriever.lib'
                subprocess.call(cmd, shell=True)


    ###################################
    ##if the users type the mitefinder2
    ###################################
    if args.mite_finder2 is not None:

        mite_finder2_exe = args.mite_finder2
        mite_finder2_profile = args.mite_finder2_profile

        ##create the dir for the mite_finder2
        mite_finder2_w_dir = working_dir + 'mite_finder2/'
        if not os.path.exists(mite_finder2_w_dir):  ##if dir not exit, it will create the dir
            os.makedirs(mite_finder2_w_dir)

        ##the default is 0.5
        if args.mf2_thr is not None:
            mf2_threshold = args.mf2_thr
        else:
            mf2_threshold = '0.5'

        ##run the mitefinder2
        cmd = mite_finder2_exe + ' -input ' + working_dir + '/genome.fas -output ' +  mite_finder2_w_dir + \
              '/mite_finder2.lib' + ' -pattern_scoring ' + mite_finder2_profile + ' -threshold ' + mf2_threshold
        print(cmd)
        subprocess.call(cmd, shell=True)

        ##cp the mitefinder2 lib to the combine_lib_dir
        cmd = 'cp ' + mite_finder2_w_dir + '/mite_finder2.lib' + ' ' + combine_lib_dir + '/mite_finder2.lib'
        print(cmd)
        subprocess.call(cmd, shell=True)


    ##################################
    ##if the users type the mitehunter
    ##################################
    if args.mite_hunter_dir is not None:

        mite_hunter_dir = args.mite_hunter_dir


        ##create the dir for the mite_finder2
        mite_hunter_w_dir = working_dir + 'mite_hunter/'
        if not os.path.exists(mite_hunter_w_dir):  ##if dir not exit, it will create the dir
            os.makedirs(mite_hunter_w_dir)


        ##install the MITE hunter
        if args.install_yes is not None:
            cmd = mite_hunter_dir + '/MITE_Hunter_Installer.pl' + ' -d ' + mite_hunter_dir + \
                  ' -f formatdb -b blastall -m mdust -M muscel'
            print(cmd)
            subprocess.call(cmd, shell=True)

        ##run the mite_hunter
        ##get the current dir
        current_dirpath = os.getcwd()
        ##change to mite_hunter_working_dir
        os.chdir(mite_hunter_w_dir)

        ##because the current dir changes, so we need to change the genome file dir
        cmd = mite_hunter_dir + '/MITE_Hunter_manager.pl' + ' -i ../genome.fas -g mite_hunter -S 12345678'
        print(cmd)
        subprocess.call(cmd, shell=True)

        ##change the dir to the previous current dir
        os.chdir(current_dirpath)


        ##cp the  mitehunter lib ot the combine_lib_dir
        cmd = 'cp ' + mite_hunter_w_dir + '/mite_hunter_Step8_singlet.fa ' + combine_lib_dir + '/mite_hunter.lib'
        subprocess.call(cmd, shell=True)

        ##redirect to the previous dirpath
        cmd = 'cd ' + current_dirpath
        subprocess.call(cmd, shell=True)

    #########################################
    ##cp the other lib to the combine lib dir
    #########################################
    if args.other_lib is not None:
        other_lib = args.other_lib
        ##cp the other_lib to the combine lib
        cmd = 'cp ' + other_lib + ' ' + combine_lib_dir + '/other.lib'
        subprocess.call(cmd, shell=True)

    #####################
    ##combine all the lib
    #####################
    lib_list = ''
    lib_file_list = glob.glob(combine_lib_dir + '/*.lib')
    for eachlib in lib_file_list:
        lib_list = eachlib + ' ' + lib_list

    cmd =  'cat ' + lib_list + ' > ' + output_dir + '/combine.lib'
    subprocess.call(cmd, shell=True)

    ##change the lowercase lib to the uppercase lib
    cmd = 'awk \'/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}\' ' + output_dir + '/combine.lib' + ' > '  + output_dir + '/combine_uppercase.lib'
    subprocess.call(cmd, shell=True)



if __name__ == "__main__":
    main()




