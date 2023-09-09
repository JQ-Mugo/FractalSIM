#!/usr/bin/python

#####################################################
# PARAMETERS:
#
recombination_rate = 1e-8
effective_population_size = 10000
#
#
#####################################################


import numpy as np
import os, sys, fileinput, exceptions, types, time,datetime,string,re
import logging, pickle
from random import random, choice, shuffle,randint,randrange,uniform,sample
import subprocess
import warnings
import itertools
from operator import mul
from os import system, popen
from math import sqrt,log,sqrt,tanh,exp
logger = logging.getLogger()
logger.setLevel(logging.INFO)
console = logging.StreamHandler()
console.setLevel(logging.INFO)
logger.addHandler(console)
warnings.filterwarnings("ignore")

class fracUtils:
    '''     Useful functions
    '''
    def __init__(self):
        return None

    def logical(self, Valeur):
        '''Convert a string to logical variable 
	Valeur: a string, NO or YES
        Return a logical variable
        '''
        if Valeur == "YES":
            return True
        else:
            return False

    def safe_open(self,filename):
        '''Opens a file in a safe manner or raise an exception if the file doesn't exist filename: a string, name of a file
        Return: file object'''

        if os.path.exists(filename):
            return open(filename)
        else:
            self.printHead()
            sys.stderr.write("ERROR ==> No such file or directory\"%s\"\n"%format(filename))
            return False

    def check_files_dicts(self,params):
        '''Check necessary parameters in a dict or a list params: a dictionary, of option:value '''
        dict = {} ; list=[]
        if type(params) == type(dict):
            for key in params:
                if params[key][0] in ['Fsnp','Flegend']:
                    if os.path.exists(params[key][1]) == False:
                       sys.stderr.write("ERROR ==> No such file or directory \"%s\"\n\n" % params[key][-1])
                       sys.exit(1)
        else:
            for i in params:
                if os.path.exists(i) == False:
                    sys.stderr.write("ERROR ==> No such file or directory \"%s\"\n\n" %i)
                    sys.exit(1)

    def check_other_params(self):
	'''Function checks the format of the input parameters in details'''
	
	#check the content of input file
	if self.homogeneous:
		if len(self.reference_haps) != len(self.positions_list):
        		sys.stderr.write('ERROR ==> Parameter entry: No of SNPs in the haplotype file differs with SNPs File\n')
                	self.terminate()
		no_of_haps = len(self.reference_haps[0])
		if no_of_haps%2 != 0:
			sys.stderr.write('ERROR ==> Parameter entry: No of haplotype should be even in the haplotype file\n')
                        self.terminate()


	#control evolution approach
        if self.homogeneous == self.admixed:
        	sys.stderr.write('ERROR ==> Parameter entry: You can either simulate a homogenous or an admixture population\n')
                self.terminate()
	
	# Disease params check
	if self.disease:
		if self.homogeneous:
			if self.null == self.causal:
                		sys.stderr.write("ERROR ==> You can either simulate a disease null model or a causal model\n")
                        	self.terminate()
		else:
			if self.post_admix_disease:
				if self.admix_null == self.admix_causal:
                                	sys.stderr.write("ERROR ==> You can either simulate a post admix_disease null model or a causal model\n")
                                	self.terminate()

	# Admixture params check
	if self.admixed:
		if len(self.anc_label) != len(self.parentalFiles):
			sys.stderr.write("ERROR ==> Parameter entry: Number of parental files is different from the ancestry labels\n")
                      	self.terminate()
		if self.pre_admix_sim:
			if len(self.pre_admix_option_dict) != len(self.parentalFiles):
        			sys.stderr.write("ERROR ==> Parameter entry: Number of parental files is different from the sections in pre_admixture simulation file\n")
                		self.terminate()
        
		if self.admix_single_point == self.admix_multi_point:
        		sys.stderr.write("ERROR ==> You can either simulate admixture single point or admixture multi-point\n")
                	self.terminate()
	   
		if self.selection:
			if self.post_admix_selection == False and len(self.pre_admix_select_pop) == 0:
				sys.stderr.write("ERROR ==> Either pre_admixture selection or post_admixture selection should be true to simulate selection under admixture simulation\n")
                		self.terminate()
			if self.post_admix_selection == True and len(self.pre_admix_select_pop) > 0:
				sys.stderr.write("ERROR ==> Either pre_admixture selection or post_admixture selection should be true and not both \n")
                                self.terminate()
			if len(self.pre_admix_select_pop) > 0:
				self.pre_admixed_selection = True
			else:
				self.pre_admixed_selection = False
		else:
			if self.post_admix_selection == True or len(self.pre_admix_select_pop) > 0:
				sys.stderr.write("WARNING ==> Select selection option to simulate either post-admixture selection or pre-adimxture selection\n")
				self.post_admix_selection = False

		if self.disease:
			if self.post_admix_disease == False and len(self.pre_admix_disease_pop) == 0:
				sys.stderr.write("ERROR ==> Either pre_admixture disease or post_admixture disease should be true to simulate disease under admixture simulation\n")
                       		self.terminate()
			if self.post_admix_disease == True and len(self.pre_admix_disease_pop) != len(self.anc_label):
				null_model_pop = ",".join([i for i in self.anc_label if i not in self.pre_admix_disease_pop.keys()])
				sys.stderr.write("WARNING ==> If admix disease is true then pre_admixture disease should be true for all parental populations \n")
				sys.stderr.write("%s to be simulated under a disease null model\n\n"%null_model_pop)
		else:
			if len(self.pre_admix_disease_pop) > 0:
				sys.stderr.write("WARNING ==> Disease simulation is false, all pre_admixture disease simulation set to false\n")
			if self.post_admix_disease == True:
				sys.stderr.write("WARNING ==> Disease simulation is false, post_admixture disease simulation set to false\n")
				self.post_admix_disease = False
		
		for ref_file in self.parentalFiles:
			F = open(ref_file);D = F.read(); F.close(); reference_haps = re.split ('\n', D)
			while "" in reference_haps: reference_haps.remove ('')
		

    def readPreAdmixOption(self):
        parafile = open(self.pre_admix_parafile,"rb")
	options =  ['parental_pop_label','simulation_size','pre_admixture_selection','fitness_para_file','pre_admixture_disease','null_model','causal_model','no_of_cases','no_of_controls','disease_SNPs_position','hete_homozygous_risk']	
	
	self.pre_admix_option_dict = {} ; return_dict = {} ; self.pre_admix_select_pop = [] ; self.pre_admix_disease_pop = {}
	no = 0
	return_dict[no] = ["\nPre_Admixture_Simulation_Options","\n"]
        for lineno,line in enumerate(parafile):
		if not line.startswith("#"):
			data = line.rstrip().split(":")
			if len(data) != 1:
				if data[0] not in options:
					sys.stderr.write("ERROR ==> Check pre_admixture parameter file line %s: unknown option\n\n"%str(lineno+1))
                                	self.terminate()
				if data[0] == "parental_pop_label" and data[1] not in self.anc_label:
					sys.stderr.write("ERROR ==> Check pre_admixture parameter file line %s parental label unknown\n\n"%str(lineno+1))
                                        self.terminate() 
				else:
					if data[0] == "parental_pop_label":
						options_dict = {}
						key = data[1]
						if lineno == 0:
							self.pre_admix_option_dict[key] = {}
				if data[0] in ['simulation_size','no_of_cases','no_of_controls']:
                                	try:
                                        	data[1] = int(data[1])
                                       	except (TypeError,NameError,ValueError,IndexError):
                                        	sys.stderr.write("ERROR ==> Check pre_admixture parameter file line %s \n\n"%str(lineno+1))
                                               	self.terminate()


				if data[0] in ['pre_admixture_selection','pre_admixture_disease','null_model','causal_model']:
					data[1] = self.logical(data[1])
				if self.disease == False and data[0] in ['pre_admixture_disease','null_model','causal_model']:
					data[1] = False

				if data[0] in ['null_model','causal_model','no_of_cases','no_of_controls','disease_SNPs_position','hete_homozygous_risk']:
					if options_dict['pre_admixture_disease'] == False:
						data[1] = False
					else:
						if data[0] in ['no_of_cases','no_of_controls']:
							try:
                                                		data[1] = int(data[1])

                                        		except (TypeError,NameError,ValueError,IndexError):
                                                		sys.stderr.write("ERROR ==> Check pre_admixture parameter file line %s \n\n"%str(lineno+1))
                                                		self.terminate()
						if data[0] == 'disease_SNPs_position':
                                        		try:
                                                		data[1] = [int(i) for i in data[1].split(",")]
                                        		except (TypeError,NameError,ValueError,IndexError):
                                                		sys.stderr.write("ERROR ==> Check pre_admixture parameter file line %s \n\n"%str(lineno+1))
                                                		self.terminate()
                                        		for pos in data[1]:
                                        			if pos not in self.positions_list:
                                                			sys.stderr.write("ERROR ==> Check pre_admixture parameter file line %s : position %s not in SNPs file \n\n"%(str(lineno+1),str(pos)))
                                                			self.terminate()

						if data[0] == 'hete_homozygous_risk':
                                        		try:
                                                		data[1] = [float(j) for j in data[1].split(",")]
                                        		except (TypeError,NameError,ValueError,IndexError):
                                                		sys.stderr.write("ERROR ==> Check pre_admixture parameter file line %s \n\n"%str(lineno+1))
                                                		self.terminate()


				if self.selection == False and data[0] in ['pre_admixture_selection','fitness_para_file']:
					data[1] = False
				else:
					if data[0] == 'fitness_para_file':
                                        	data[1] = self.infolder + data[1]
                                        	self.check_files_dicts([data[1]])
				if data[0] == 'null_model':
					if data[1] == False and options_dict['pre_admixture_disease'] == True:
						if self.admix_null == True:
							logger.info("ERROR => To simulate a null model for the post-admix disease simulation, all pre-admix disease simulations should be under a disease null model\n")
                        				self.terminate()
				if data[0] == 'pre_admixture_selection':
					if data[1] == True:
						self.pre_admix_select_pop.append(key)
				if data[0] == 'pre_admixture_disease':
					if data[1] == True:
						self.pre_admix_disease_pop[key] = []

				if data[0] == 'null_model' or data[0] == 'causal_model':
					if options_dict['pre_admixture_disease'] == True:
						if data[1] == True: 
							self.pre_admix_disease_pop[key] = data[0]

				options_dict[data[0]] = data[1]
				if key not in self.pre_admix_option_dict:
					self.pre_admix_option_dict[key] = options_dict
					no = no+1
					return_dict[no] = ["\n        "+key,""]
				if 'pre_admixture_selection' in options_dict:
					if data[0] == 'fitness_para_file' and options_dict['pre_admixture_selection'] == False:
						data[1] = False
				if 'pre_admixture_disease' in options_dict:
					if data[0] in ['null_model','causal_model','no_of_cases','no_of_controls','disease_SNPs_position','hete_homozygous_risk'] and options_dict['pre_admixture_disease'] == False:
						data[1] = False
				if data[1] != False and data[0] != "parental_pop_label":
					no = no+1
					if type(data[1]) == type([]):
						data[1] = ",".join([str(i) for i in data[1]])
                                        return_dict[no] = [data[0],data[1]]
	return return_dict			

    def terminate(self):
        '''Terminate the process'''
        log_rem = os.getcwd()
        try:

            pkl_files = [fil for fil in os.listdir(self.outfolder) if fil.endswith('pkl')]
            for fil in pkl_files:
              pass
            os.system("cp"+" "+log_rem+'/'+self.logFile+" "+self.outfolder)
        except:

            filenames = os.listdir(os.curdir)
            for filename in filenames:
                if os.path.isfile(filename) and filename.endswith('.log'):
                        os.system("rm"+" "+ log_rem+"/*.log")
        finally:
            log = log_rem+'/'+self.outfolder+self.logFile
            print'Log file generated in '+log+'\nHave a nice day!\n'

            filenames = os.listdir(os.curdir)
            for filename in filenames:
                if os.path.isfile(filename) and filename.endswith('.log'):
                        pass
            sys.exit(1)


class fractInit(fracUtils):

    try:
        logger.removeHandler(logger.handlers[1])
    except:
        pass
    log_rem = os.getcwd()
    logger.setLevel(logging.INFO)
    filenames = os.listdir(os.curdir)
    for filename in filenames:
        if os.path.isfile(filename) and filename.endswith('.log'):
                os.system("rm"+" "+log_rem+"/*.log")

    logFile = 'fractalsim-'+str(time.time()).split('.')[0]+'.log'
    fh = logging.FileHandler(logFile, mode='w')
    logger.addHandler(fh)

    logger.info("\n***********************************************************************************************************************  ")
    logger.info("               FractalSIM: A Multi-scenario Whole Genome-based Medical Population Genetics Simulation Framework    ")
    logger.info("                                        Computational Biology Group                                                 ") 
    logger.info("                                2016, University of Cape Town, South Africa                                          ")
    logger.info("                                            Verson 1.0 Beta                                                         ")
    logger.info(  "***********************************************************************************************************************\n")


    def __init__(self, argv, logFile=logFile):

        '''Initializing FractalSIM by reading the parameter file'''
        self.argv = [argv]
        self.logFile = logFile
        popnames = []
        if len(self.argv) == 0 or self.argv == ['']:
            logger.info('Command line usage: %s <parameter file>  ' % sys.argv[0])
            logger.info('eg: python fractalsim.py parancFractal.txt\n')
            sys.exit(1)
        elif len(self.argv) == 1:
            try:
                self.paramFile = self.argv[0]
                if os.path.exists(os.getcwd()+'/'+self.paramFile):
                    inFile = open(os.getcwd()+'/'+self.paramFile)
                elif os.path.exists(self.paramFile):
                    inFile = open(self.paramFile)
                else:
                    logger.info('\nERROR ==> Failed to process the input, check the parameters!\n')
                    self.terminate()
                    sys.exit()
                rows = {}
                self.opts = [];self.opts1 = {};myList1 = []
		myList = ['homogenous','admixture','disease','selection','null_model','causal_model','pre_admixture_simulation','post_admix_selection','post_admix_disease','admix_single_point','admix_multi_point']
                for line in inFile:
                    data = line.split()
                    if not line.startswith("#"):
                        if len(data) != 0:
                            rows[data[0].split(":")[0]] = data[0].split(":")[1]
                            data1 = data[0].split(":")[0]
                            data2 = data[0].split(":")[1]
                            if data1 in myList:
                                data2 = self.logical(data2)
                            self.opts.append(data2)
                            self.opts1[data1] = data2
                            myList1.append(data1)
                for file in rows:
                        if file in ["Fhap","Fsnp","Fgenetic_map","Fitness","pre_admix_option_file","multi_point_paraFile"]:
                                rows[file] = rows['infolder']+rows[file]
                if len(self.opts1) != 30:

                    sys.stderr.write('ERROR ==> Missing parameters!!! Failed to process the input, check the parameters!\n')
                    sys.exit()
                try:
                    self.outfolder = self.opts1['outfolder']
                    path = os.getcwd()
                    if os.path.exists(self.outfolder):
                        pass
                    else:
                        os.makedirs(self.outfolder)

                except IndexError:
                    sys.stderr.write('ERROR ==> Can not create directory. Please create a directory OUT in your working directory\n')
                    sys.exit(1)
		allList =['infolder','homogenous','admixture','disease','selection','Fsnp','Fgenetic_map','Ffitness','Fhap','samples_to_add','pre_admixture_simulation','pre_admix_option_file','admix_single_point','admix_multi_point','multi_point_paraFile','post_admix_selection','post_admix_disease','post_admix_sim_sample_size','parental_files','anc_label','anc_proportions','samples_to_simulate','generations_to_simulate','null_model','causal_model','no_of_cases','no_of_controls','disease_SNPs_position','hete_homozygote_risk','outfolder']
                List1 = []
                for i in range(len(myList1)):
                    if myList1[i] != allList[i]:
                        List1.append(myList1[i])
                if len(List1) != 0:
                    sys.stderr.write('\nInvalid option: '+','.join(List1))
                    sys.stderr.write('\nERROR ==> Failed to process the input, check the parameters!\n')
                    self.terminate()
                    sys.exit()

                for param in myList:
                    rows[param] = self.logical(rows[param])

                self.infolder = self.opts1['infolder']

                #simulation options
                self.homogeneous = self.opts1['homogenous']
                self.admixed = self.opts1['admixture']
                self.disease = self.opts1['disease']
		self.selection = self.opts1['selection']
                
		# compulsory input files
                self.snpFile = self.infolder + self.opts1['Fsnp']
                self.geneticFile = self.infolder + self.opts1['Fgenetic_map']
		self.check_files_dicts([self.snpFile,self.geneticFile])

		self.genetic_map_dict = {}
                for line in fileinput.input(self.geneticFile):
                	if fileinput.lineno() > 1:
                        	data = line.split()
                                self.genetic_map_dict[int(data[0])] = float(data[1])

		#create a dictionary for the positions per chrom and a list of the positions
                self.positions_dict = {} ; self.positions_list = []
                for line in fileinput.input(self.snpFile):
                        data = line.split()
                        if fileinput.lineno() == 1:
                                self.positions_dict[int(data[1])] =  [int(data[3])]
                        else:
                                if int(data[1]) in self.positions_dict:
                                        self.positions_dict[int(data[1])].append(int(data[3]))
                                else:
                                        self.positions_dict[int(data[1])] = [int(data[3])]
                        self.positions_list.append(int(data[3]))

		if self.selection:
                        self.fitness_file = self.infolder + self.opts1["Ffitness"]
                        self.check_files_dicts([self.fitness_file])
                else:
                        del rows["Ffitness"]


		if self.disease:
			if self.admixed:
				self.admix_cases = int(self.opts1["no_of_cases"])*2
                        	self.admix_controls = int(self.opts1["no_of_controls"])*2
				self.admix_null = self.opts1["null_model"]
                        	self.admix_causal = self.opts1["causal_model"]

			else:
				self.cases = int(self.opts1["no_of_cases"])*2
                        	self.controls = int(self.opts1["no_of_controls"])*2
				self.null = self.opts1["null_model"]
                        	self.causal = self.opts1["causal_model"]


                 	#Convert positions to integers and risks allele to float as expected
                     	self.hete_homo_risk = []; self.disease_snp_pos = [] ; self.risk_snp_no = []
                        for pos in self.opts1['disease_SNPs_position'].split(","):
                        	self.disease_snp_pos.append(int(pos))
                     	for allele in self.opts1['hete_homozygote_risk'].split(","):
                        	self.hete_homo_risk.append(float(allele))
                      	for i,position in enumerate(self.positions_list):
                        	if int(position) in self.disease_snp_pos:
                        		self.risk_snp_no.append(i)
		else:
                	for option in ['null_model','causal_model','no_of_cases','no_of_controls','disease_SNPs_position','hete_homozygote_risk']:
                        	del rows[option]


		#homogenous input file
                if self.homogeneous:
                        self.hapFile = self.infolder + self.opts1['Fhap']
                        self.samples_to_add = int(self.opts1['samples_to_add'])*2
                        self.check_files_dicts([self.hapFile])

                        #create a dictionary of reference population for homogenous simulation
                        F = open(self.hapFile);D = F.read(); F.close(); self.reference_haps = re.split ('\n', D)
                        while "" in self.reference_haps: self.reference_haps.remove ('')

                        if self.selection:
                                self.homog_selection = True
			else:
				self.homog_selection = False
                else:
                        for option in ['Fhap','samples_to_add']:
                                del rows[option]


                #Admixture parameters only read if admixture is true
                if self.admixed:
			parentalFiles = self.opts1['parental_files'].split(",")
			self.anc_label = self.opts1['anc_label'].split(",")
				
			for j in range(len(parentalFiles)):
                                parentalFiles[j] = self.infolder + parentalFiles[j]
                        self.check_files_dicts(parentalFiles)
                        self.parentalFiles = parentalFiles			

			self.pre_admix_sim = self.opts1["pre_admixture_simulation"]
			if self.pre_admix_sim:
				self.pre_admix_parafile = self.infolder + self.opts1["pre_admix_option_file"]
				self.check_files_dicts([self.pre_admix_parafile])
				pre_admix_options_dict = self.readPreAdmixOption()
			else:
				del rows["pre_admix_option_file"]

                        self.admix_single_point = self.opts1["admix_single_point"]
                        self.admix_multi_point = self.opts1["admix_multi_point"]
                        if self.admix_single_point:
                                self.generations = int(self.opts1['generations_to_simulate'])
                                self.anc_props = [float(i) for i in self.opts1['anc_proportions'].split(",")]
				if self.disease:
					self.no_to_simulate = self.admix_cases + self.admix_controls
				else:
					self.no_to_simulate = int(self.opts1['samples_to_simulate'])*2
				del rows['multi_point_paraFile']				
				
                        else:
                                self.multipoint_parafile = self.infolder + self.opts1['multi_point_paraFile']
				self.check_files_dicts([self.multipoint_parafile])
				for option in ['samples_to_simulate','generations_to_simulate','anc_proportions']:
					del rows[option]

			self.post_admix_selection =  self.opts1['post_admix_selection']
			self.post_admix_sim_size = int(self.opts1['post_admix_sim_sample_size'])*2
			if self.selection and self.post_admix_selection == False:
				del rows["Ffitness"] ; del rows["post_admix_sim_sample_size"]

			if self.selection == False and self.post_admix_selection == True:
				rows['post_admix_selection'] = False
			if self.post_admix_selection == False:
				del rows["post_admix_sim_sample_size"]
				
                        self.post_admix_disease =  self.opts1['post_admix_disease']
			if self.disease and self.post_admix_disease == False:
				self.post_admix_disease = True

			if self.disease == False and self.post_admix_disease == True:
				 del rows['post_admix_disease']
		else:
			for option in ['pre_admixture_simulation','post_admix_selection','post_admix_sim_sample_size','post_admix_disease','admix_single_point','admix_multi_point','parental_files','anc_label','pre_admix_option_file','samples_to_simulate','generations_to_simulate','anc_proportions','multi_point_paraFile']:
				del rows[option]
	
		i = 1
                self.Params = {}
                for param in myList1:
                    if param in rows and rows[param] != False:
                        self.Params[i] = [param,rows[param]]
                        i+=1
		if self.admixed:
			if self.pre_admix_sim:
				idxs = self.Params.keys() ; idxs.sort()
				init_idx = idxs[-1]+1
				op_keys =  pre_admix_options_dict.keys()
				op_keys.sort()
				for op in op_keys:
					self.Params[init_idx] = pre_admix_options_dict[op]
					init_idx = init_idx + 1
                self.check_files_dicts(self.Params)


            except IndexError, TypeError:
                sys.stderr.write('ERROR ==> Failed to process the input, check the parameters!\n\n')
                self.terminate()
                sys.exit(1)
        else:
            logger.info('Command line usage: %s <parameter file>  ' % sys.argv[0])
            logger.info('eg. python fractalSIM.py parameters_file.txt\n')
            self.terminate()
            sys.exit(1)

class OtherUtlis:
    '''extra pvalues functions'''
    def __init__(self):
    	return None

    def transpose_file(self,Infile,Outfile,get_control_risk_snp):
	'''inverts the haplotype data set so that SNP becomes the columns and individuals are in the rows'''

        transpose_starttime = datetime.datetime.today()
        logger.info("Transposing file %s ...\n\n"%str(Infile))
		
	rows={}
       	transp_file = self.outfolder+Outfile
       	transp = open(transp_file,"wt")
	infile = open(Infile,"r")
		
	conc = False ; not_conc = False		
        for line in infile:
	        data = line.split()
                if len(data) == 1:
			conc = True
			ref_length = len(data[0])
		else:
			not_conc = True
			ref_length = len(data)
                break
	if conc:
		for line in fileinput.input(Infile):
               		data=line.split()
              		rows[fileinput.lineno()]=data[0]
               		sample = len(data[0])
			if sample == ref_length:
				pass
			else:
				logger.info("ERROR => line %s in file %s inconsistent \n"%(str(fileinput.lineno()),str(Infile)))
                		self.terminate()
	elif not_conc:
		for line in fileinput.input(Infile):
               		data=line.split()
               		rows[fileinput.lineno()]=data
               		sample=len(data)
			if sample == ref_length:
                                pass
                        else:
                                logger.info("ERROR => line %s in file %s inconsistent \n"%(str(fileinput.lineno()),str(Infile)))
                                self.terminate()
	else:
		logger.info("ERROR ==> Haplotype file format unknown")
		self.terminate()
	
	if get_control_risk_snp:
		self.control_risk_snps = {}
		for snp in self.risk_snp_no:
			self.control_risk_snps[snp] = []
		
	for i in range(sample):
               	for des in rows:
                       	transp.writelines(rows[des][i])
			if get_control_risk_snp and i in self.control_risk_snps:
				self.control_risk_snps[i].append(rows[des][i])
               	transp.write("\n")

       	transp.close()
        infile.close()
		

    def invert_dict(self, dict,assign_positions):
	'''inverts a dictionary so that  the dictinary length becomes the length of the key elements'''
       	column_no = len(dict[dict.keys()[0]])
       	inverted_dict = {}
        if assign_positions:
                for i,pos in enumerate(self.position_list):
                        inverted_dict[pos]=[]
                        new_hap = []
                        for j in dict.keys():
                                new_hap.append(dict[j][i])
                        inverted_dict[pos] = ''.join(new_hap)
        else:
                for k in range(column_no):
                        inverted_dict[k+1] = []
                        new_hap = []
                        for j in dict.keys():
                                new_hap.append(dict[j][k])
                        inverted_dict[k+1] = ''.join(new_hap)
        return inverted_dict

class readdata():
    '''Mapping module, which maps snps to genes and computes the p-value at gene level.'''
    def __init__(self):
        return None


class statUtils():
    '''Estimates the egde weight using the LD between snps in genes and constructs the LD-weighted PPI network.'''
    def __init__(self):
        return None

    def binomial (self, n, k, p):
        return np.exp(self.log_choose(n,k) + k*log(p) + (n-k)*log(1-p))


    def log_choose (self, a, b):
        if (a == 0 or b == 0): return 0
        x = self.log_factorial (a)
        y = self.log_factorial (b)
        z = self.log_factorial (a-b)
        return x - y - z


    def log_factorial (self, x):
        if ((x)==0): return 0
        elif ((x)==1): return 0
        elif ((x)==2): return 0.693147180559945
        elif ((x)==3): return 1.79175946922805
        elif ((x)==4): return 3.17805383034795
        elif ((x)==5): return 4.78749174278205
        elif ((x)==6): return 6.57925121201010
        elif ((x)==7): return 8.52516136106541
        elif ((x)==8): return 10.6046029027453
        elif ((x)==9): return 12.8018274800815
        elif ((x)==10): return 15.1044125730755
        elif ((x)==11): return 17.5023078458739
        elif ((x)==12): return 19.9872144956619
        else: return .5 * log (2*np.pi) + ((x)+.5)*log((x)) - (x) + 1/(12*(x))

class  Homogeneous_Simulator(fractInit,OtherUtlis,statUtils):
    '''This module computes the recombination, selection and mutation points, it inverses the given haplotype file and amplifies the data to the 
    number stated in the parameters file.'''

    def __init__(self):
	return None

    def check_simulation_params(self):
        if self.selection:
		self.fitness_dict = {}
		self.selection_positions = []		
                for line in fileinput.input(self.fitness_file):
                        if fileinput.lineno() > 1:
                                data = line.split()
                                if data[0] not in [str(i) for i in self.positions_list]:
                                        sys.stderr.write('ERROR ==> Fitness file %s : position %s not in list of SNPs \n'%(str(self.fitness_file),str(data[0])))
                                        self.terminate()
                                if float(data[1]) != 1.0:
                                        sys.stderr.write('ERROR ==> Fitness file: s0 in line %i should be 1.0\n'%fileinput.lineno())
                                        self.terminate()
                                if float(data[2]) < 0.0 or float(data[3]) < 0.0:
                                        sys.stderr.write('ERROR ==> Fitness file: s1 and s2 in line %i should be greater than 0\n'%fileinput.lineno())
                                        self.terminate()
				self.selection_positions.append(int(data[0]))
				self.fitness_dict[self.positions_list.index(int(data[0]))] = [float(data[1]),float(data[2]),float(data[3])]

    def rand_exponential (self, rate):
        u = random()
        x = -(1/rate)*(log(1-u))
        return x
   
    def recomb_points(self, m):
        rec_points = []
        chrom_initial_pos = 0
        for chrom in self.positions_dict:
                chrom_positions = self.positions_dict[chrom]
                no_of_chrom_positions = len(chrom_positions)
                number_of_base_pairs = chrom_positions[-1] - chrom_positions[0]
                chrom_rec_points = []
                p = 0
                rec_pos = 0 ; gen_dist_list = []
                for i,pos in enumerate(chrom_positions):
                        if pos in self.genetic_map_dict:
				physical_dist =  (chrom_positions[i] - chrom_positions[i-1]) *1e-6
				recomb_rate = self.genetic_map_dict[pos]
                        else:
				recomb_rate = recombination_rate
				physical_dist =  (chrom_positions[i] - chrom_positions[i-1]) 

                        if i == 0:
                        	gen_dist = 0.0
                        else:
                        	gen_dist = physical_dist * recomb_rate
                                gen_dist_list.append(gen_dist)

			delta = (4*effective_population_size*gen_dist)/float(m)
			prob1 = exp(-delta) + (1 - exp(-delta))/float(m)
			prob2 = (1 - exp(-delta))/float(m)
			
			if prob2 >= prob1:
                                rec_pos = i
                                pos_append = i+chrom_initial_pos
                                chrom_rec_points.append (pos_append)

                chrom_initial_pos = chrom_initial_pos + no_of_chrom_positions
                chrom_rec_points.append(chrom_initial_pos)
                rec_points.extend(chrom_rec_points)

        return rec_points

    def mut_points (self, mutation_rate,no_of_snps):
        self.mut_list = []
        max_no_of_mutations = min (int (no_of_snps * mutation_rate * 10 + 10), no_of_snps)
        mutation_rate_vector = [float] * (max_no_of_mutations+1)

        for i in xrange (0, max_no_of_mutations+1):
                mutation_rate_vector[i] = self.binomial (no_of_snps, i, mutation_rate)
        for i in xrange (0, max_no_of_mutations+1):
                mutation_rate_vector[i] += mutation_rate_vector[i-1]

        r = random()
        for no_of_mutations in xrange (0, max_no_of_mutations+1):
                if r < mutation_rate_vector[no_of_mutations]:
                       break
        for i in xrange (no_of_mutations):
                self.mut_list.append (randint (0, no_of_snps-1))

    def apply_mut (self, hap):
        hap = list (hap)
        for p in self.mut_list:
                a = hap[p]
                if a == '0': b = '1'
                elif a == '1': b = '0'
                elif a == '2': b = '3'
                elif a == '3': b = '2'
                elif a == '4': b = '5'
                elif a == '5': b = '4'
                elif a == '6': b = '7'
                elif a == '7': b = '6'
                hap[p] = b
        return ''.join (hap)

    def simulation_option(self,no_of_haps):
	'''Determine the samples to be simulated depending on the option of simulation given'''
	if not self.disease:
		simu_hapsize = [int(no_of_haps + self.samples_to_add)]
	else:
		simu_hapsize = [int(no_of_haps + self.controls + self.cases)]

	return simu_hapsize


    def homogeneous_simulation(self, input_file):
	'''Main simulation function, simulates homogeneous population with account to selection, recombination, mutation.'''
	self.check_simulation_params()
        amplify_starttime = datetime.datetime.today()
        logger.info("Start amplification of %s at time:%s\n\n"%(str(input_file),str(amplify_starttime)))

        #Transpose the input haplotype file
        self.transpose_file(input_file,"Ref.inv.haps",False)
        transp_file = self.outfolder+"Ref.inv.haps"
        transposed = open(transp_file);D = transposed.read(); transposed.close(); transp_haps = re.split ('\n', D)
        while "" in transp_haps: transp_haps.remove ('')
	
	no_of_snps_hapFile = len(list(transp_haps[0]))
        no_of_haps = len(transp_haps) 
	no_of_snps = len (self.positions_list) 

	# check the SNPs file and the haplotype file compatibility
	if no_of_haps %2 != 0:
		sys.stderr.write("ERROR ==> Number of haplotypes should be even\n")
               	self.terminate()
	elif no_of_snps != no_of_snps_hapFile:
		sys.stderr.write("ERROR ==> Number of SNPs in the haplotypes file %s differs with SNPs file\n"%str(input_file))
                self.terminate()
	else:
		pass

	# Specifications for output file
        if self.disease:
		if self.causal:
			outfile = "".join([str(i) for i in input_file.split("/")[-1].split(".")[0]])+".controls.haps"
		else:
			outfile = "".join([str(i) for i in input_file.split("/")[-1].split(".")[0]])+".amp"
	else: 
                outfile = "".join([str(i) for i in input_file.split("/")[-1].split(".")[0]])+".amp"
	
	self.ampOutFile = self.outfolder+outfile
       	ampFile =  open (self.ampOutFile, 'wt')

	# obtain the simulation sample size depending on the simulation option
	simu_hapsize = self.simulation_option(no_of_haps)

	simulated_haps = {}
        for count,size in enumerate(simu_hapsize):
		simulated_haps.clear()
		
		#Calculate theta; used to calculate the mutation rate
                theta_deno = 0
		samples_to_add = size - no_of_haps
                for j in range(1,samples_to_add):
                       	theta_deno = theta_deno + 1/float(j)
                theta_est = 1/float(theta_deno)

		# write the reference haplotypes in a dict
                for line in fileinput.input(transp_file):
                	data = line.rstrip()
                        simulated_haps[fileinput.lineno()] = data

		# partition reference haplotypes depending on disease alleles and allele under selection
		if self.disease or self.homog_selection:
                        partitioned_haplotype_dict = self.getPartitionedHapsDict(input_file)
			
		# Obtain the risk alleles for the disease SNPs, from disease model
		if self.disease:
			if self.causal:
				simulated_SNPs = self.getDsl(admix=False)
			else: # self.null
				if self.homog_selection:
					simulated_SNPs = self.getFitSNPs(size)
		
		# Obtain the fit alleles for the SNPs under selection, from the selection model
		if self.homog_selection and not self.disease:
			simulated_SNPs = self.getFitSNPs(size)
		

		next_h_dict = {}						
                for m in xrange (no_of_haps+1,size+1):
			ind_no_list = range(1,m)
			hap_no = m - (no_of_haps+1)
                        
			# Obtain the mutation and recombination points
			mutation_rate = theta_est/float(m+theta_est)
			self.mut_points (mutation_rate,no_of_snps)
			cut_list = self.recomb_points (m)

			# haplotype simulation
                     	st = 0
                        new_hap = []
                        for cut_idx,cut in enumerate(cut_list):
				condition = "choose_h"
				# choose a case_ind segment that contains a risk SNP
				if self.disease: 
					if self.causal:
						dsl_list = []
						for dsl_no in self.risk_snp_no:
							if dsl_no in range(st,cut):
								dsl_list.append(dsl_no)

						if len(dsl_list) != 0:
							condition = "h_choosen"
							if self.homog_selection:
								for sel_no in self.fitness_dict:
                                                                	if sel_no in range(st,cut) and sel_no not in dsl_list:
                                                                        	dsl_list.append(sel_no)

							prev_h,hap_seg = self.getSeg(st,cut,hap_no,simulated_haps,partitioned_haplotype_dict,simulated_SNPs,dsl_list)

				if condition == "choose_h":
					# choose a segment that contains a SNP under selection
					if self.homog_selection:
						condition2 = "choose_h"
						fit_list = []
						for sel_no in self.fitness_dict:
							if sel_no in range(st,cut):
								fit_list.append(sel_no)
					
						if len(fit_list) != 0:
							condition2 = "h_choosen"
							prev_h,hap_seg = self.getSeg(st,cut,hap_no,simulated_haps,partitioned_haplotype_dict,simulated_SNPs,fit_list)
					
					# choose a segment that contains neither a  risk SNP nor a SNP under selection
					if self.homog_selection == False  or condition2 == "choose_h":
						#if cut_idx == 0:
						#	h = randint (1, m-1)
						#	prev_h = h
						#else:
						#	h = randint (1, m-1)
						#	while h == prev_h:
						#		h = randint (1, m-1)
						#	prev_h = h
						h = choice(ind_no_list)
						prev_h = h
						ind_no_list.pop(ind_no_list.index(h))
						hap_seg = list(simulated_haps[h])[st:cut]

				new_hap.extend(hap_seg)
				st = cut
				# Refill the individuals list once its empty
				if len(ind_no_list) < 1:
					ind_no_list = range(1,m)

			# apply mutation to the simulated haplotype
			new_hap = self.apply_mut (''.join (new_hap))

                     	simulated_haps[m] = new_hap
			
			# Update the haplotype partitioned dict
			if self.disease or self.homog_selection:
				partitioned_haplotype_dict = self.updatePartitionList(m,new_hap,partitioned_haplotype_dict)

			# write new hap in a file immediately
                  	if count == 0:
                        	print >> ampFile, new_hap
                       	 	ampFile.flush()

	# Remove the initial haplotypes
        for key in range(1,no_of_haps+1):
                del simulated_haps[key]
	
	#Transpose the simulated data if not disease
	if self.homogeneous and not self.disease:
                controls_haps_file = "".join([str(i) for i in input_file.split("/")[-1].split(".")[0]])+".amplified.haps"
                self.transpose_file(self.outfolder+outfile, controls_haps_file,False)
                system("rm"+" "+self.outfolder+outfile)

        logger.info('Amplification of %s complete!'%input_file)
        logger.info('Current population size: %i haplotypes'% int(samples_to_add))
        amplify_finishtime = datetime.datetime.today()
        #logger.info("Finish amplification at time:%s\n\n"%str(amplify_finishtime))
        #logger.info("Amplification took time:%s\n\n"%str(amplify_finishtime - amplify_starttime))
        
	if self.disease:
		return simulated_haps

class case_controlSimulator(fractInit,fracUtils,statUtils,Homogeneous_Simulator):
    '''This module simulates the case_control data, according to the given parameters in the parameters file'''
 
    def __init__(self):
        return none

    def check_disease_parameters(self):
        if len(self.hete_homo_risk)%2 != 0 or len(self.hete_homo_risk) != 2*len(self.disease_snp_pos):
                logger.info("\nERROR ==>: Please check the risk parameters\n")
                self.terminate()
        for i in self.disease_snp_pos:
                if i not in self.positions_list:
                        logger.info("\nERROR ==>: Position %i of the disease SNP not in SNPs File\n"%i)
                        self.terminate()

    def determinePop(self,hapfile):
	'''Assigns the population to be simulated and calls homogenous population simulation function'''
        amp_pops = {}
	if self.homogeneous:
		prefix_filename = "".join([str(i) for i in self.hapFile.split("/")[-1].split(".")[0]])
		#if self.null:
		#	samples_to_add = [self.controls + self.cases]
		#else:
		#	samples_to_add = [self.controls + self.cases]
		
		amp_pops = self.homogeneous_simulation(self.hapFile)

        else: #self.admixed true
		prefix_filename = "".join([str(i) for i in hapfile.split("/")[-1].split(".")[0]])
                amp_pops = self.homogeneous_simulation(hapfile)

        return prefix_filename, amp_pops

    def getPartitionedHapsDict(self,reference_haps_file):
	'''This function splits the specified snps (disease snps or under selection) in the reference haplotypes according to the allele in 
	each haplotype, and writes into a dictionary, with the disease snp no as the key, and the values as the haplotypes in a list'''

	if self.homog_selection and self.disease:
		pos_list = self.fitness_dict.keys()
		if self.causal:
			for pos in self.risk_snp_no:
				if pos not in pos_list:
					pos_list.append(pos)
	else:
		if self.disease:
			if self.causal:
				pos_list = self.risk_snp_no
			else:
				pos_list = []
		else:
			pos_list = self.fitness_dict.keys()			

	partitioned_haps_dict = {}
	if len(pos_list) != 0: 
		F = open(reference_haps_file);D = F.read(); F.close(); reference_haps = re.split ('\n', D)
        	while "" in reference_haps: reference_haps.remove ('')
			
		for snp_no in range(len(reference_haps)):
			if snp_no in pos_list:
				snp = reference_haps[snp_no].rstrip().split(" ")
				if len(snp) == 1:
					snp = list(reference_haps[snp_no])

				list0 = [] ; list1 = []
				for hap_no in range(len(snp)):
					if snp[hap_no] == str(0):
						list0.append(hap_no+1)
					else:
						list1.append(hap_no+1)
				partitioned_haps_dict[snp_no] = [list0,list1]
			else:
				pass
	return partitioned_haps_dict

    def updatePartitionList(self,m,new_hap,partitioned_haplotype_dict):
	'''The function update the partitioned_haplotype_dict obtained in function getPartitionedHapsDict, once a new haplotype is simulated 
	and joins the sample haps'''
	
	if self.homog_selection and self.disease:
                pos_list = self.fitness_dict.keys()
                if self.causal:
                        for pos in self.risk_snp_no:
                                if pos not in pos_list:
                                        pos_list.append(pos)
        else:
                if self.disease:
                        if self.causal:
                                pos_list = self.risk_snp_no
			else:
				pos_list = []
                else:
                        pos_list = self.fitness_dict.keys()
	
	if len(pos_list) != 0:
		for snp in pos_list:
			hap_allele = new_hap[snp]
			if hap_allele == str(0):
				partitioned_haplotype_dict[snp][0].append(m)
			else:
				partitioned_haplotype_dict[snp][1].append(m)
	return partitioned_haplotype_dict

    def getFitAlleles(self):
        '''The function implements the fitness model to generate the genotypes of the SNPs under selection pressure'''

        total_selection_snps = len(self.fitness_dict)
        fitness_perSNP_list = [self.fitness_dict[fit] for fit in self.fitness_dict]

	overall_fitness_dict = {} ; overall_fitness_list = []
	for fit_list in itertools.product(*fitness_perSNP_list):
		F = reduce(mul, list(fit_list))
		overall_fitness_list.append(F)
		if F not in overall_fitness_dict:
			overall_fitness_dict[F] = [list(fit_list)]
		else:
			overall_fitness_dict[F].append(list(fit_list))

	cut_off = (overall_fitness_list[0] + overall_fitness_list[-1])/2.0
	overall_fitness_list.sort()			
	cut_off_list = [overall_fitness_list[-2],overall_fitness_list[-1]]
	
	fit_snps = self.fitness_dict.keys()
        picked  = False
        while not picked:
                genotype_list = []
                fit_choice = [] 
                for snp_no in self.fitness_dict:
                        geno = randint(0,2)
                        genotype_list.append(geno)

		F_sim = 1.0
                for i,g in enumerate(genotype_list):
                     	F_sim = F_sim * self.fitness_dict[fit_snps[i]][g]   

		if F_sim in cut_off_list:
                	picked = True

        del overall_fitness_dict
        return genotype_list

    def getFitSNPs(self,inds):
        '''The function splits the simulated genotypes for the SNPs under selection into haplotype format'''

        fitness = {}
	for ind in range(inds):
                ind_genotype = self.getFitAlleles()
                for i,snp in enumerate(self.fitness_dict):
                        if ind_genotype[i] == 0:
                                if snp in fitness:
                                        fitness[snp].extend(['0','0'])
                                else:
                                        fitness[snp] = []
                                        fitness[snp].extend(["0","0"])
                        elif ind_genotype[i] == 1:
                                if snp in fitness:
                                        fitness[snp].extend(["0","1"])
                                else:
                                        fitness[snp] = []
                                        fitness[snp].extend(["0","1"])
			else:
                                if snp in fitness:
                                        fitness[snp].extend(["1","1"])
                                else:
                                        fitness[snp] = []
                                        fitness[snp].extend(["1","1"])

        #for h in fitness:
        #        print h,fitness[h].count("0")/float(len(fitness[h])), fitness[h].count("1")/float(len(fitness[h]))
        return fitness

    def getRisk(self,pop):
        ''' The function implements the regression model to generate the genotypes of the disease risk SNPs'''

        total_disease_snps = len(self.risk_snp_no)
        risk_dict = {} #dictionary of the risks given 

        l=0
        for k in range(1,total_disease_snps+1):
                 for j in range(1,3):
                        risk_dict[k,j] = self.hete_homo_risk[l+j-1]
                 l += 2
        t = 0
        highest_risk = 1 ;  lowest_risk = 1
        for dsl_no in self.risk_snp_no:
                risk1 = 1
                risk2 = self.hete_homo_risk[t]
                risk3 = self.hete_homo_risk[t+1]
                high_risk = max(risk1,risk2,risk3)
                low_risk = min(risk1,risk2,risk3)
                highest_risk = highest_risk * high_risk
                lowest_risk = lowest_risk * low_risk
                t += 2

        picked  = False
        while not picked:
                genotype_list = []
                risks_choice = [] ; base_list = []
                for snp_no in self.risk_snp_no:
                        geno = randint(0,2)
                        genotype_list.append(geno)
                for m,i in enumerate(genotype_list):
                        m += 1
                        if i == 0:
                                risks_choice.append(log(1))
                        else:
                                risks_choice.append(log(risk_dict[m,i]))

                risk_total = sum(risks_choice) + 0.0001
		base_total = log((highest_risk + lowest_risk)/2.0) + 0.0001		

                prob_risk = 1/float(1 + exp(-risk_total))
                base_prob = 1/float(1 + exp(-base_total))
		
		if pop == "case":
                	if prob_risk >= base_prob:
                        	picked = True
		else:
			if prob_risk < base_prob:
                                picked = True	
        del risk_dict
        return genotype_list
	
    def getDsl(self,admix):
	'''The function splits the simulated genotypes for the disease SNPs into haplotype format'''

	dsl = {} 
	if admix == True:
		size_list = [self.admix_controls/2,self.admix_cases/2]
	else:
		size_list = [self.controls/2, self.cases/2]

	for s,inds in enumerate(size_list):
        	if s == 0:
			pop = "control"
		else:
			pop = "case"
		for ind in range(inds):
			ind_genotype = self.getRisk(pop)
			for i,snp in enumerate(self.risk_snp_no):
				if ind_genotype[i] == 0:
					if snp in dsl:
						dsl[snp].extend(['0','0'])
					else:
						dsl[snp] = []
						dsl[snp].extend(["0","0"])
				elif ind_genotype[i] == 1:
					if snp in dsl:
                                        	dsl[snp].extend(["0","1"])
                                	else:
                                        	dsl[snp] = []
                                        	dsl[snp].extend(["0","1"])
				else:
					if snp in dsl:
                                        	dsl[snp].extend(["1","1"])
                                	else:
                                        	dsl[snp] = []
                                        	dsl[snp].extend(["1","1"])
	
	if admix != True:
		if self.homog_selection:
                	fitness_SNPs = self.getFitSNPs(sum(size_list))
                	for no in self.risk_snp_no: 
				if no not in self.fitness_dict:
					fitness_SNPs[no] = dsl[no]
			dsl = fitness_SNPs
	for h in dsl:
		if admix != True:
			control_list = dsl[h][:self.controls] ; case_list = dsl[h][self.controls:]
		else:
			control_list = dsl[h][:self.admix_controls] ; case_list = dsl[h][self.admix_controls:]

		#print "controls",h,":",control_list.count("0")/float(len(control_list)), control_list.count("1")/float(len(control_list))
                #print "cases",h,":",case_list.count("0")/float(len(case_list)), case_list.count("1")/float(len(case_list)),"\n"	
	return dsl

	
    def getSeg(self,st,cut,hap_no,simulated_haps,partitioned_haps,simulated_snps,snp_no_list):
	'''The function selects a segment from the sample hap, if the segment has a disease risk SNP in absence of selection'''
	hap_seg = []
	
	# When only 1 snp, choose a h that satisfy the simulated snp if found, else "mutate" a random h
	if len(snp_no_list) == 1:
		snp_no = snp_no_list[0]
                if simulated_snps[snp_no][hap_no] == str(0):
                	if len(partitioned_haps[snp_no][0]) == 0:
                        	h = choice(partitioned_haps[snp_no][1])
				#h = np.random.choice(partitioned_haps[snp_no][1])
                                choosen_hap = list(simulated_haps[h])
                                choosen_hap[snp_no] = "0"
                                hap_seg = choosen_hap[st:cut]
			else:
                        	h = choice(partitioned_haps[snp_no][0])
				#h = np.random.choice(partitioned_haps[snp_no][0])
                                hap_seg = list(simulated_haps[h])[st:cut]
		else:
                	if len(partitioned_haps[snp_no][1]) == 0:
                        	h = choice(partitioned_haps[snp_no][0])
				#h = np.random.choice(partitioned_haps[snp_no][0])
                                choosen_hap = list(simulated_haps[h])
                                choosen_hap[snp_no] = "1"
                                hap_seg = choosen_hap[st:cut]
              		else:
                        	h = choice(partitioned_haps[snp_no][1])
				#h = np.random.choice(partitioned_haps[snp_no][1])
                                hap_seg = list(simulated_haps[h])[st:cut]
		prev_h = h
	else:
		seg = []
		while len(snp_no_list) >= 2:
			# Option1: find a h that satisfies all the simulated snps
			common_haps = {}
			for snp_counter,snp_no in enumerate(snp_no_list):
				if snp_counter == 0:
					if simulated_snps[snp_no][hap_no] == str(0): 
						prev_list = partitioned_haps[snp_no][0]
					else:
						prev_list = partitioned_haps[snp_no][1]

				else:
					if simulated_snps[snp_no][hap_no] == str(0):
						prev_list = list( set(prev_list) & set(partitioned_haps[snp_no][0]) )
					else:
						prev_list = list( set(prev_list) & set(partitioned_haps[snp_no][1]) )
					if len(prev_list) > 0:
						common_haps[snp_counter] = prev_list

			if len(prev_list) > 0:
				h = choice(prev_list)
				#h = np.random.choice(prev_list)
				seg.extend(list(simulated_haps[h])[st:cut])
				del snp_no_list[:]
				prev_h = h
				break
			else:
				# Option 2: Choose a h that satisfies atleast two of the simulated snps (in order)
				# then split the segment into two. 
				if len(common_haps) != 0:
					large_seg_dsl = snp_no_list[max(common_haps.keys())]
					h1 = choice(common_haps[max(common_haps.keys())])
					#h1 = np.random.choice(common_haps[max(common_haps.keys())])
					seg.extend(list(simulated_haps[h1])[st:large_seg_dsl])
					st = large_seg_dsl
					discard_list = [snp_no_list[d] for d in range(max(common_haps.keys())+1)]
					[snp_no_list.remove(i) for i in discard_list]
					prev_h = h1
				else:
					break

		# Option 3: If more than 2 snps still remain, split them and obtain a h for each segments
		if len(snp_no_list) >= 2:
			for j in range(1,len(snp_no_list)):
				diff = snp_no_list[j] - snp_no_list[j-1]
				if diff >= 2:
					separator = snp_no_list[j-1]+int(diff/2.0)
				else:
					separator = snp_no_list[j-1]
				if simulated_snps[snp_no_list[j-1]][hap_no] == str(0) and len(partitioned_haps[snp_no_list[j-1]][0]) != 0:
					h2 = choice(partitioned_haps[snp_no_list[j-1]][0])
					#h2 = np.random.choice(partitioned_haps[snp_no_list[j-1]][0])
				else:
					if len(partitioned_haps[snp_no_list[j-1]][1]) > 0:
						h2 = choice(partitioned_haps[snp_no_list[j-1]][1])
						#h2 = np.random.choice(partitioned_haps[snp_no_list[j-1]][1])
					else:
						h2 = choice(partitioned_haps[snp_no_list[j-1]][0])
						#h2 = np.random.choice(partitioned_haps[snp_no_list[j-1]][0])
				seg.extend(list(simulated_haps[h2])[st:separator])
				st = separator
				prev_h = h2

		# Option 4, if only 1 snp is remaining, obtain a h that satisfies the simulated snp 
		if len(snp_no_list) != 0:
			last_dsl = snp_no_list[-1]
			if simulated_snps[last_dsl][hap_no] == str(0) and len(partitioned_haps[last_dsl][0]) != 0:
                		h3 = choice(partitioned_haps[last_dsl][0])
				#h3 = np.random.choice(partitioned_haps[last_dsl][0])
                	else:
         			if len(partitioned_haps[last_dsl][1]) > 0:
                        		h3 = choice(partitioned_haps[last_dsl][1])
					#h3 = np.random.choice(partitioned_haps[last_dsl][1])
                        	else:
                        		h3 = choice(partitioned_haps[last_dsl][0])
					#h3 = np.random.choice(partitioned_haps[last_dsl][0])
			seg.extend(list(simulated_haps[h3])[st:cut])
			prev_h = h3
		hap_seg = seg	

	return prev_h,hap_seg
		
    def split_case_control(self,hapfile):
	'''Splits the population received from the determinePop function into cases and controls'''
        if self.homogeneous:
		pop_name,received_dict = self.determinePop(False)
	else:
		pop_name,received_dict = self.determinePop(hapfile)

        amp_pops = {}
        for i,key in enumerate(received_dict):
                amp_pops[i+1] = received_dict[key]
	cases_haps = {}
        controls_haps = {}
	if self.null:
        	rand_list=[]
        	curr_hapsize = len(amp_pops)
        	for ind in range(self.cases):
                	rand_choose= randrange(1,curr_hapsize)
                	while rand_choose in rand_list:
                        	rand_choose = randrange(1,curr_hapsize)
                	cases_haps[ind+1] = amp_pops[rand_choose]
                	rand_list.append(rand_choose)
        	counter = 0
        	for ind2 in range(1,curr_hapsize+1):
                	if ind2 not in rand_list:
                        	controls_haps[counter+1] = amp_pops[ind2]
                        	counter += 1
                	if counter == self.controls:
                        	break
	else:
		for ind in range(1,self.controls+1):
			controls_haps[ind] = amp_pops[ind]
			del amp_pops[ind]

		for i,ind2 in enumerate(amp_pops.keys()):
			cases_haps[i+1] = amp_pops[ind2]

        del amp_pops
        return [pop_name,cases_haps,controls_haps]
    

    def case_control_output(self,hapfile):
	'''Write the simulated cases and controsl in the respective files'''
        split_haps = self.split_case_control(hapfile)
        prefix_filename = split_haps[0]; split_cases_haps = split_haps[1] ; split_controls_haps = split_haps[2]

	if self.homogeneous:
        	cases_snps = self.invert_dict(split_cases_haps,assign_positions=False)
        	controls_snps = self.invert_dict(split_controls_haps,assign_positions=False)
	else:
		cases_snps = split_cases_haps ; controls_snps = split_controls_haps

        number_of_alleles = self.cases + self.controls
	
	if self.homogeneous:
		if self.null:
        		controls_fout = open(self.outfolder+"controls_homogeneous_nullmodel.haps","wt")
        		cases_fout = open(self.outfolder+"cases_homogeneous_nullmodel.haps","wt")
		else:
			controls_fout = open(self.outfolder+"controls_homogeneous_causalmodel.haps","wt")
                	cases_fout = open(self.outfolder+"cases_homogeneous_causalmodel.haps","wt")
	else:
		controls_fout = open(self.outfolder +prefix_filename+".controls.haps",'wt')
		cases_fout = open(self.outfolder +prefix_filename+".cases.haps",'wt')

	for key in cases_snps.keys():
		data1 = list(cases_snps[key])
		cases_fout.write("".join([str(i) for i in data1])+"\n")

	for key2 in controls_snps.keys():
		data2 = list(controls_snps[key2])
                controls_fout.write("".join([str(i) for i in data2])+"\n")

	cases_fout.close()
        controls_fout.close()
	del cases_snps, controls_snps, split_haps, split_cases_haps, split_controls_haps 

    def case_control_homogeneous_simulation(self,hapfile):
	'''Check simulation parameters and simulate cases and controls'''
        self.check_disease_parameters()
	self.case_control_output(hapfile)
        logger.info("\ncase-control homogeneous simulation complete!\n\n")

	
class Admixture(fracUtils,statUtils,Homogeneous_Simulator,case_controlSimulator):
    '''Simulates Admixture both multi and single point, case control admixed population'''
    def __init__(self):
        return None

    def check_admixture_parameter(self):
	# check paremeters specific to multi-point admixture simulation
        if self.admix_multi_point:
                for line in fileinput.input(self.multipoint_parafile):
                        data = line.split()
                        if fileinput.lineno() == 1:
                                self.no_to_simulate = int(data[0])*2
                                anc_labels = data[1:]
                                if len(anc_labels) != len(self.parentalFiles)+1:
                                        sys.stderr.write("ERROR ==> The number of ancestry labels in the multipoint admixture parameter file differs parental files\n")
                                        self.terminate()
			elif fileinput.lineno() == 2:
				prev_gen = int(data[0])
				if float(data[1]) != 0.0:
                                        sys.stderr.write("ERROR ==> Admixture proportion in line %i should always be 0.0\n"%fileinput.lineno())
                                        self.terminate()
				anc_props = [float(i) for i in data[1:]]
                                if sum(anc_props) != 1.0:
                                	sys.stderr.write("ERROR ==> Please check the ancestry proportions in line %i of multipoint parameter file, they should add to 1\n"%fileinput.lineno())
                                        self.terminate()
				all_anc_labels = []
				for prp in anc_props:
					if float(prp) > 0.0:
						all_anc_labels.append(anc_labels[anc_props.index(prp)])
				if len(all_anc_labels) < 2:
					sys.stderr.write("ERROR ==> You can only simulate 2-way admixed and above, check line %i in multi-point parameter file\n"%fileinput.lineno())
                                        self.terminate()
                        else:
				curr_gen = int(data[0])
                                if len(data[1:]) != len(self.parentalFiles)+1:
                                        sys.stderr.write("ERROR ==> Please check the ancestry proportions in line %i, they do not match the parental files\n"%fileinput.lineno())
                                        self.terminate()
				if float(data[1]) == 0.0 or float(data[1]) == 1.0:
					sys.stderr.write("ERROR ==> Admixture proportion in line %i should not be 0.0 or 1.0\n"%fileinput.lineno())
                                        self.terminate()
                                anc_props = [float(i) for i in data[1:]]
                                if sum(anc_props) != 1.0:
                                        sys.stderr.write("ERROR ==> Please check the ancestry proportions in line %i of multipoint parameter file, they should add to 1\n"%fileinput.lineno())
                                        self.terminate()
				for prp2 in anc_props:
                                        if float(prp2) > 0.0 and anc_labels[anc_props.index(prp2)] in all_anc_labels:
						sys.stderr.write("ERROR ==> population %s was included in  generation %i, cannot be included again in generation %i \n"%(anc_labels[anc_props.index(prp2)],int(prev_gen),int(curr_gen)))
						self.terminate()
					elif float(prp2) > 0.0 and anc_labels[anc_props.index(prp2)] not in all_anc_labels and anc_labels[anc_props.index(prp2)] != "admixed":
                                                all_anc_labels.append(anc_labels[anc_props.index(prp2)])
					else:
						pass
				if curr_gen <= prev_gen:
					sys.stderr.write("ERROR ==> The generation in line %i should be higher than previous generation in the multi-point parameter file\n"%fileinput.lineno())
					self.terminate()
				prev_gen = curr_gen

        else:
		# checks parameters specific to single point admixture simulation
                if len(self.parentalFiles) != len(self.anc_label):
                        sys.stderr.write('ERROR ==> Parameter entry: The number of parental files is different from ancestry labels \n')
                        self.terminate()

                if len(self.anc_props) != len(self.anc_label):
                        sys.stderr.write('ERROR ==> Parameter entry: The number of populations label is different from the number of population proportions \n')
                        self.terminate()
                if round(sum(self.anc_props)) != 1.0:
                        sys.stderr.write('ERROR ==> Parameter entry: The first population is Admixed with proportion 1 \n')
                        self.terminate()
	
	print "\nparameter loading for admixture done!!!\n"

    def checkPreAdmixParam(self,pop):
	if self.pre_admix_disease:
		if self.cases < 2*(self.admix_cases):
			logger.info("WARNING ==> No. of cases in pre-admixture disease simulation should be twice the admixture cases to simulate")
                	logger.info("cases sample size for %s changed to %s\n\n"%(pop,str(self.admix_cases)))
			self.cases = 2*(self.admix_cases)

		if self.controls < 2*(self.admix_controls):
                        logger.info("WARNING ==> No. of controls in pre-admixture disease simulation should be twice the admixture controls to simulate")
                        logger.info("controls sample size for %s changed to %s\n\n"%(pop,str(self.admix_controls)))
			self.controls = 2*(self.admix_controls)

	else:
		if self.post_admix_disease == False:
			if self.samples_to_add < 2*self.no_to_simulate:
				logger.info("WARNING ==> simulation size in pre-admixture disease simulation should be twice the admixture sample to simulate")
                        	logger.info("simulation size for %s changed to %s\n"%(pop,str(self.no_to_simulate)))
				self.samples_to_add = 2*self.no_to_simulate
 
    def getPopParams(self,pop):
	params = self.pre_admix_option_dict[pop]
	self.pre_admix_selection = params['pre_admixture_selection']
	self.pre_admix_disease = params['pre_admixture_disease']
	if self.pre_admix_disease:
		self.null = params['null_model'] 
		self.causal = params['causal_model']
		self.cases = int(params['no_of_cases'])*2
                self.controls = int(params['no_of_controls'])*2
		if self.causal:
			self.hete_homo_risk = params['hete_homozygous_risk']
			self.disease_snp_pos = params['disease_SNPs_position']
			self.risk_snp_no =  []
			for i,position in enumerate(self.positions_list):
                                if int(position) in self.disease_snp_pos:
                                        self.risk_snp_no.append(i)

	elif self.post_admix_disease  and self.pre_admix_disease == False:
		self.null = True
		self.causal = False
		self.cases = 2*(self.admix_cases)
		self.controls = 2*(self.admix_controls)
		self.pre_admix_disease_pop[pop] = "null_model"
	else: 
		self.samples_to_add = int(params['simulation_size'])*2
	if self.pre_admix_selection:
		self.fitness_file = params['fitness_para_file']
		self.homog_selection = True
	else:
		self.homog_selection = False

	self.checkPreAdmixParam(pop) 

    def shuffle_list(self,list_objects):
        ''' receive a list of values and a integer number of sampling. It returns a shuffling list of values'''
        if self.admix_single_point:
                gen = int(self.generations)
        else:
                gen = int(self.generations[-1])
        GEN = gen*2+1
        for i in range(GEN):
                shuffle(list_objects)
        return list_objects

    def choosePop(self,state):
        anc_idx_list = [idx+1 for idx in range(len(self.anc_props))]
	#pick_anc = np.random.choice(anc_idx_list,p=self.anc_props)
        pick_anc = choice(anc_idx_list)
        return pick_anc


    def break_points(self,state,positions,generation):
	'''Returns the individual segments for a specific chromosome and the corresponding ancestry '''
	break_points_list = []; pop_choice_list = [] #; prev_anc_len = 0
	for i,pos in enumerate(positions):
		if i == 0:
			pass
		else:
			if int(pos) in self.genetic_map_dict and self.genetic_map_dict[pos] > 0:
				recomb_rate = self.genetic_map_dict[pos]
			else:
				recomb_rate = 4 * effective_population_size * recombination_rate

			physical_dist = (positions[i] - positions[i-1])*1e-6
			genetic_dist =  physical_dist * recomb_rate
		
			break_prob1 =  1 - exp(-(genetic_dist)*generation)
			break_prob2 = exp(-(genetic_dist)*generation) + (1 - exp(-(genetic_dist)*generation))

			if break_prob1 == break_prob2:
				if len(self.anc_props) == 1.0:
					pop_choice = 1
				else:
					pop_choice = self.choosePop(state)
				break_points_list.append(pos)
				pop_choice_list.append(pop_choice)
				state = pop_choice 
			else:
				continue
			
	if len(break_points_list) == 0 or break_points_list[-1] != positions[-1]:
		break_points_list.append(positions[-1])
		if len(self.anc_props) == 1:
			pop_choice = 1
		else:
			pop_choice = self.choosePop(state)
		pop_choice_list.append(pop_choice)
		state = pop_choice 
	return [state, break_points_list,pop_choice_list]			

    def check_generation(self,segment):
        ''' Receive a list of ancestry labels in order their segment contribution along a chromosome. It returns the number of breakpoints = Generation = ancestry switch '''
        gen = 0
        for seg in range(len(segment)):
                if seg ==0 :
                        anc = segment[seg]

                else:
                        if anc == segment[seg]:
                                pass
                        else:
                                gen =gen +1
                                anc = segment[seg]
        return gen


    def check_seg_generation(self, SEG, ANC):
        ''' This function receives a first round of ancestry segment and related label, Check the breakpoint regarding the given generations. 
        The function aims to break down long segments into chuck following based on the generation. It returns the individual final ancestry 
        segments and generation.'''

        get_all_seg = []
        GEN = 0

        for chrm_seg in range(len(SEG)):
                GEN = GEN + self.check_generation(ANC[chrm_seg]) ### NEED to be EXpanded @@@@@@@@@@@
                for point in range(len(SEG[chrm_seg])):
                        get_all_seg.append(str(ANC[chrm_seg][point])+str(':')+str(SEG[chrm_seg][point]))
        return get_all_seg, GEN

    def genome_recordSegment(self,GEN):
        '''This function aims to write the simulated admixed ancestry segments. It calls two functions [break_points() and
         check_seg_generati_generation ()]'''

	segments_starttime = datetime.datetime.today()
        if self.admix_single_point:
                out = self.outfolder+"ancestry_segments_haps.out"
        else:
                out = self.outfolder+"multipoint_ancestry_segment_haps.out"+str(GEN)
        fin = open(out,"wt")
        fin.write("\t".join(self.anc_label)+"\n")
	
        genome_av_gen = []
	state = 0
	segment_dict = {}
	simulated_samples = self.no_to_simulate * 5
	chrom_list = self.getMapPos.keys()
	chrom_list.sort()

        for ind in range(1,simulated_samples): # Loop for over the given number of the individual to be simulated
		gen = randint(1,(2*GEN)-1)
                rec_genome_seg = []; rec_genome_anc = []#;state = [0]
                for chrm in chrom_list:
                        positions = self.getMapPos[chrm]; positions.sort()
			chrm_seg_result =  self.break_points(state,positions,gen)
			rec_genome_seg.append(chrm_seg_result[1]); rec_genome_anc.append(chrm_seg_result[-1]); state = chrm_seg_result[0]
                get_all, gen = self.check_seg_generation(rec_genome_seg, rec_genome_anc) ### CHECK
		segment_dict[ind] = get_all 
                genome_av_gen.append(gen)

	ind_list = sample(range(1,len(segment_dict)),self.no_to_simulate)
	ind_list.sort()
	for ind2 in ind_list:
                fin.write("\t".join(segment_dict[ind2])+"\n")

	del segment_dict
        fin.close()
	
	segments_endtime = datetime.datetime.today()
        logger.info("\n\nSegments generation took %s time..\n\n"%(str(segments_endtime - segments_starttime)))
        #logger.info('Theorethical %s and Simulated %s generation since admixture happened ' % (str(self.generations),str(np.mean(genome_av_gen))))
        logger.info('The true ancestry Segments written at %s'% out)

    def convert_hap2geno(self, hap_file, geno_file):
	'''Converts a haplotype file into  genotype format'''
        rows={}
        rows1={}
        fix=open(geno_file,"wt")
        print"\nConverting haplotypes to genotypes data format ..."
        for line in fileinput.input(hap_file):
                data=line.split()
                if fileinput.lineno()%2==0:
                        rows[header]=data[0]
                elif fileinput.lineno()%2!=0:
                        header=fileinput.lineno()
                        rows1[fileinput.lineno()]=data[0]
        if len(rows)==len(rows1):

		for idx in rows:
			for i in range(len(rows[rows.keys()[0]])):
                                if rows[idx][i]!=rows1[idx][i]:
                                        if rows[idx][i]=="1" and rows1[idx][i]=="0":
                                                fix.writelines("1")
                                        elif rows[idx][i]=="0" and rows1[idx][i]=="1":
                                                fix.writelines("1")
                                        elif rows[idx][i]=="0" and rows1[idx][i]=="?":
                                                fix.writelines("1")
                                        elif rows[idx][i]=="?" and rows1[idx][i]=="0":
                                                fix.writelines("1")
                                        elif rows[idx][i]=="?" and rows1[idx][i]=="1":
                                                fix.writelines("0")
                                        elif rows[idx][i]=="1" and rows1[idx][i]=="?":
                                                fix.writelines("0")
                                else:
                                        if rows[idx][i]=="0" and rows1[idx][i]=="0":
                                                fix.writelines("2")
                                        elif rows[idx][i]=="1" and rows1[idx][i]=="1":
                                                fix.writelines("0")
                                        else:
                                                fix.writelines("?")
                                                fix.writelines("?")
                        fix.write("\n")
                fix.close()
                print "Done with conversion \n\n"
        else:
                print"\nThere is inconsistancy in the genotype length!"
                sys.exit()

    def evolve_parental_pops(self):
	'''Homogeneous simulation of the parental populations to be used in the admixture simulation'''
        haplotype_amp_files = [] 
        for i,file in enumerate(self.parentalFiles):
		self.getPopParams(self.anc_label[i])
		if self.disease:
			#simulate the homogeneous under disease, call the homogeneous simulation
			self.case_control_homogeneous_simulation(file)

			prefix_filename = "".join([str(j) for j in file.split("/")[-1].split(".")[0]])
                        cases_outfile = self.outfolder+prefix_filename+".cases.haps"
                        controls_outfile = self.outfolder+prefix_filename+".controls.haps" 

                        haplotype_file = [controls_outfile,cases_outfile]

		else:
			self.homogeneous_simulation(file)
			amp_file = self.outfolder+"".join([str(m) for m in file.split("/")[-1].split(".")[0]])+".amp"
			hap_file = self.outfolder+"".join([str(n) for n in file.split("/")[-1].split(".")[0]]) +".haps"
			system("mv" + " "+ amp_file + " "+ hap_file)
			haplotype_file = [hap_file]

  		haplotype_amp_files.append(haplotype_file)
        self.parentalFiles = haplotype_amp_files

    def split_segs(self):
	'''Splits the simulated admixture segments into cases and controls'''
        sim_segs = {}
        discard_list=[]
        cases_outfile = open(self.outfolder+"admixture_cases_null_segments_haps.out","wt")
        if self.admix_null:
                controls_outfile = open(self.outfolder+"admixture_controls_null_segments_haps.out","wt")
        else:
                controls_outfile = open(self.outfolder+"admixture_controls_causal_segments_haps.out","wt")

        curr_popsize = int(self.no_to_simulate)

        if self.admix_single_point:
                seg_file = self.outfolder+"ancestry_segments_haps.out"
        else:
                seg_file = self.outfolder+"multipoint_ancestry_segments_haps.out"+str(self.generations[-1])

        cases_outfile.write('\t'.join([i for i in self.anc_label])+"\n")
        controls_outfile.write('\t'.join([i for i in self.anc_label])+"\n")
        for line in fileinput.input(seg_file):
                if fileinput.lineno() > 1:
                        line = line.rstrip()
                        sim_segs[fileinput.lineno()-1] = line
        rand_choose= randint(1,curr_popsize)
	
        for ind in range(self.admix_cases):
                rand_choose= randint(1,curr_popsize)
                while rand_choose in discard_list:
                        rand_choose = randint(1,curr_popsize)
                cases_outfile.write(''.join([sim_segs[rand_choose]])+"\n")
                discard_list.append(rand_choose)
        
	counter = 0
        for ind2 in range(1,curr_popsize+1):
                if ind2 not in discard_list:
                        controls_outfile.write(''.join([sim_segs[ind2]])+"\n")
                        counter += 1
                if counter == self.admix_controls:
                        break
	del sim_segs
        cases_outfile.close()
        controls_outfile.close()

    def get_case_segments(self):
	'''Alters the segments with the disease SNPs in the case segments under a disease causal model'''
        position_dict = self.getMapPos
        chrom_list = position_dict.keys()
	chrom_list.sort()
        fout = open(self.outfolder+"admixture_cases_causal_segments_haps.out","wt")

        fout.write('\t'.join([i for i in self.anc_label])+"\n")
        for line in fileinput.input(self.outfolder+"admixture_cases_null_segments_haps.out"):
                if fileinput.lineno() > 1:
			chrom = chrom_list[0]
			counter = 0
        		positions = position_dict[chrom]
                        data1 = line.split()
			chrom_segmnts = [] ; chrom_segments_sublist = [] ; chrom_counter = 0
			for seg1 in data1:  # break the segments for a given ind into different chroms
				pop_s,pos_s = seg1.split(':')
				positions_2 = position_dict[chrom_list[chrom_counter]]
				if int(pos_s) == positions_2[-1]:
					chrom_segments_sublist.append(seg1)
					chrom_segmnts.append(chrom_segments_sublist) 
					chrom_segments_sublist = [] ; chrom_counter += 1
				else: 
					chrom_segments_sublist.append(seg1)
			for chrom_segs in chrom_segmnts:
				ind_segment_list = []
	                	positions_list = [] ; anc_list = [] ; endMarker_list = []
                        	for i in range(len(chrom_segs)):
                                	pop,pos = chrom_segs[i].split(':')
                                	anc_list.append(pop)
                                	endMarker_list.append(pos)
                        	segment_taken = False
                        	for seg in range(len(chrom_segs)):
                                	if segment_taken == False:
                                        	anc = anc_list[seg] ; endMarker = endMarker_list[seg]
                                        	endMarker_idx = positions.index(int(endMarker))
                                        	if seg == 0 or len(positions_list) == 0:
                                                	seg_position_list = positions[0:endMarker_idx+1]
                                        	else:
                                                	startMarker_idx = positions.index(positions_list[-1])+1
                                                	seg_position_list = positions[startMarker_idx:endMarker_idx+1]
                                        	disease_seg = False
                                        	for i in self.disease_snp_pos:
                                                	if i in seg_position_list:
                                                        	disease_seg = True
                                                        	break
                                        	if disease_seg:
                                                	distance = int(seg_position_list[-1]) - int(seg_position_list[0])
                                                	if distance < 1000000: #merge two segment and skip next
                                                        	if seg != len(chrom_segs)-1: # not last segment
                                                                	pick_anc_idx = anc_list.index(anc) + 1
                                                                	pick_anc = anc_list[pick_anc_idx]
                                                                	ind_segment_list.append(str(anc)+":"+str(endMarker_list[seg+1]))
                                                                	segment_taken = True
                                                        	else: #last segment
                                                                	if len(ind_segment_list) == 0.0:
                                                                        	ind_segment_list.append(str(anc)+":"+str(endMarker))
                                                                	else:
                                                                        	ind_segment_list.remove(ind_segment_list[-1])
                                                                        	ind_segment_list.append(str(anc)+":"+str(endMarker))
                                                	else: # distance > 1000000
                                                        	seg_len = len(seg_position_list)
                                                        	new_endMarker_idx = int(round(1/4.0 * seg_len))
                                                        	new_endMarker = seg_position_list[new_endMarker_idx]
                                                        	ind_segment_list.append(str(anc)+":"+str(new_endMarker))
                                                        	pick_anc = self.anc_props.index(max(self.anc_props))+1
                                                        	if pick_anc == anc:
                                                                	anc_props = self.anc_props
                                                                	pick_anc = max(anc_props.remove(pick_anc))+1
                                                                	if pick_anc == anc_list[seg+1]:
                                                                        	anc_props = self.anc_props
                                                                        	pick_anc = max(anc_props.remove(pick_anc))+1
                                                        	ind_segment_list.append(str(pick_anc)+":"+str(endMarker))
                                        	else:
                                                	ind_segment_list.append(str(anc)+":"+str(endMarker))
                                	else: #Segment taken is True
                                        	endMarker = endMarker_list[seg]
                                        	endMarker_idx = positions.index(int(endMarker))+1
                                        	startMarker_idx = positions.index(positions_list[-1])+1
                                        	seg_position_list = positions[startMarker_idx:endMarker_idx]
                                        	segment_taken = False
                                	[positions_list.append(i) for i in seg_position_list]
                          		if positions_list[-1] == self.getMapPos[chrom][-1]:
                                        	if len(chrom_list) > 1 and chrom != chrom_list[-1]:
                                                	counter += 1
                                               		chrom = chrom_list[counter]
                                                	positions = position_dict[chrom]
                                        	else:
                                                	pass
				fout.writelines("\t".join([str(i) for i in ind_segment_list])+"\t")
                        fout.writelines("\n")
        fout.close()

    def multipoint_seg_breakpoints(self,seg_position_list,anc_props,anc_idx,anc,prev_anc):
        len_seg_pos_list = len(seg_position_list)
        split_pos = choice(seg_position_list)
	#split_pos = np.random.choice(seg_position_list)
        split_pos_idx = seg_position_list.index(split_pos)
        split_segs = [seg_position_list[0:split_pos_idx+1],seg_position_list[split_pos_idx+1:len_seg_pos_list]]
        len_seg1 = len(split_segs[0]); len_seg2 = len(split_segs[1])
        if len(anc_props) == 1:
                if len_seg1 > len_seg2:
                        return [anc_idx[0],split_pos,anc,seg_position_list[-1]]
                else: #len_seg1 < len_seg2
                        if prev_anc == anc:
                                return [anc_idx[0],split_pos,anc,seg_position_list[-1]]
                        else:
                                return [anc,split_pos,anc_idx[0],seg_position_list[-1]]
        else: #len(anc_props) > 1
                randVal = uniform(0,1)
                sample_pop = []
                for i in range(len(anc_props)):
                        x = int(anc_props[i]*100)
                        sample_pop = sample_pop + [anc_idx[i]]*x
                sample_pop = self.shuffle_list(sample_pop)
                pick_anc = choice(sample_pop)
		#pick_anc = np.random.choice(sample_pop)
                if len_seg1 > len_seg2:
                        return [pick_anc,split_pos,anc,seg_position_list[-1]]
                else:
                        if prev_anc == anc:
                                return [pick_anc,split_pos,anc,seg_position_list[-1]]
                        else:
                                return [anc,split_pos,pick_anc,seg_position_list[-1]]
    

    def multipoint_ind_breakpoints(self,data,prev_gen_props,all_anc_labels):
        '''Alters the ancestry segments to account to multi point approach'''
        position_dict = self.getMapPos
        chrom_list = position_dict.keys()
        state = -1
        prev_gen_seg_dict = {}
        fin = open(self.outfolder+"multipoint_ancestry_segments_haps.out"+str(self.generations[-1]),"rb")
        order_anc_labels = []
        for line_no,line in enumerate(fin):
                if line_no == 0:
                        order_anc_labels.extend(line.rstrip().split("\t"))
                else:
                        prev_gen_seg_dict[line_no] = line.rstrip().split("\t")
        fin.close()

        self.generations.append(int(data[0])); true_anc_idx = [] ; all_anc_idx = []
        for i,prop in enumerate(data[1:]):
                if float(prop) > 0:
                        if all_anc_labels[i] != "admixed":
                                order_anc_labels.append(all_anc_labels[i])
                                self.anc_label.append(all_anc_labels[i])
                                true_anc_idx.append(self.anc_label.index(all_anc_labels[i]))
                self.anc_props.append(float(prop))
 
	fout = open(self.outfolder+"multipoint_ancestry_segments_haps.out"+str(self.generations[-1]),"wt")
        fout.writelines("\t".join([i for i in all_anc_labels[1:]])+"\n")
        for ind in range(1,self.no_to_simulate+1):
                gen = randint(1,2*(self.generations[-1])-1)
                for chrm in chrom_list:
                        positions = self.getMapPos[chrm]; positions.sort()
                        last_pos = int(positions[-1])
                        first_pos = int(positions[0])
                        state, break_points_list, pop_choice_list = self.break_points(state,positions,gen)
			anc_list = [int(i)-1 for i in pop_choice_list] 
                        point = int(break_points_list[0])
                        point_counter = 0
                        for seg in prev_gen_seg_dict[ind]:
                                anc,endMarker = seg.split(":")
                                seg_cur = []
                                if int(endMarker) not in positions:
                                        continue
                                else:
                                        if int(endMarker) > int(point):
                                                while int(endMarker) > int(point):
							if int(anc_list[point_counter]) == 0: #admixed
                                                        	#new_anc = np.random.choice(range(1,len(prev_gen_props)+1),p=prev_gen_props)
                                                        	new_anc = choice(range(1,len(prev_gen_props)+1))
                                                	else:
                                                        	new_anc = anc_list[point_counter]
                                                        seg_cur.append(str(new_anc)+":"+str(point))
                                                        if point_counter >= len(break_points_list)-1: break
                                                        point_counter = point_counter + 1
                                                        point = break_points_list[point_counter]
                                                
						if int(endMarker) not in break_points_list: 
                                                        seg_cur.append(str(anc)+":"+str(endMarker))
						if int(endMarker) == positions[-1]:
							if int(anc_list[point_counter]) == 0: 
								new_anc = choice(range(1,len(prev_gen_props)+1))
								#new_anc = np.random.choice(range(1,len(prev_gen_props)+1))
							else:
								new_anc = anc_list[point_counter]
							seg_cur.append(str(new_anc)+":"+str(endMarker))
                                        else:
                                                if int(endMarker) == int(point):
                                                        if point_counter >= len(break_points_list)-1:
                                                                pass
                                                        else:
                                                                point_counter = point_counter + 1
                                                                point = break_points_list[point_counter]
                                                seg_cur.append(seg)
                                        fout.writelines("\t".join(seg_cur)+"\t")
					if int(endMarker) == positions[-1]:
                                        	break
        	fout.writelines("\n")
	fout.close()

	
    def calling(self):
	for ind in range(10):
        	prev_segs =  prev_gen_seg_dict[ind]
		print prev_segs
		prev_chrom_endpos = 0
                for chrom in chrom_list:
                        positions = position_dict[chrom]
                        state, break_points_list, pop_choice_list = self.break_points(state,positions,self.generations[-1])
			print break_points_list, pop_choice_list
			anc_list = [int(i)-1 for i in pop_choice_list]
			point = break_points_list[0]
			point_counter = 0
			for seg in prev_segs:
				anc,endMarker = seg.split(":")
				if int(endMarker) <= prev_chrom_endpos:
					continue 
				else:		
					if int(endMarker) > int(point):
						if int(anc_list[point_counter]) == 0: #admixed
							#new_anc = np.random.choice(range(1,len(prev_gen_props)+1),p=prev_gen_props)
							new_anc = choice(range(1,len(prev_gen_props)+1))
						else:
							new_anc = anc_list[point_counter]
                                        	seg_cur = str(anc)+":"+str(point) +"\t"+str(new_anc)+":"+str(endMarker)
						if int(endMarker) in break_points_list: break_points_list.remove(int(endMarker))
						if point_counter >= len(break_points_list)-1:
							pass
						else:
							point_counter = point_counter + 1
							point = break_points_list[point_counter]
					else:
						if int(endMarker) == int(point):
							if point_counter >= len(break_points_list)-1:
                                                        	pass
                                                	else:
                                                        	point_counter = point_counter + 1
                                                        	point = break_points_list[point_counter]
						seg_cur = seg
					print seg_cur
                                	fout.writelines(seg_cur + "\t")
					if int(endMarker) == positions[-1]:
						prev_chrom_endpos = int(endMarker)
						break
		raw_input()
                fout.writelines("\n")
        fout.close()
        logger.info('The true ancestry Segment file is written at %s'% self.outfolder+"multipoint_ancestry_segments_haps.out"+str(self.generations[-1]))

    def simulate_multipoint_admixture_segs(self):
        '''Simualtes the multi point approach'''
        for line in fileinput.input(self.multipoint_parafile):
                data = line.split()
                if fileinput.lineno() == 1:
                        all_anc_labels = [i for i in data[1:]]
                elif fileinput.lineno() == 2:
                        self.anc_label = [] ; self.anc_props = []; self.generations = [] ; prev_gen_props  = []
                        self.generations.append(int(data[0]))
                        for i,prop in enumerate(data[1:]):
                                if float(prop) > 0:
                                        self.anc_label.append(all_anc_labels[i])
				if all_anc_labels[i] != "admixed":
                        		self.anc_props.append(float(prop))
                        prev_gen_props.extend([i for i in self.anc_props])
                        self.genome_recordSegment(self.generations[0])
                        del self.anc_props
                else:
                        self.last_gen_labels = [] ;  self.anc_props = []
                        self.multipoint_ind_breakpoints(data,prev_gen_props,all_anc_labels)
                        last_gen = int(data[0])
                        for j,p in enumerate(data[1:]):
                                if float(p) > 0 and all_anc_labels[j] != "admixed":
                                        prev_gen_props.append(float(p))
                                        self.last_gen_labels.append(all_anc_labels[j])
        return last_gen
    
    def sim_selection_segments(self, seg_file):
        '''Alter the segments that contain the SNP under selection by assigning a high proportion population to the segment'''
     
        fout = open(self.outfolder +"ancestry_selection_segments_haps.out","wt")
        fitness_index_list = self.fitness_dict.keys()
        if self.admix_multi_point:
                self.anc_props.pop(0) #remove prop corresponding to admixed pop

        if len(self.anc_props) == 2 or len(self.anc_props) == 1:#two-way or one-way 
                choicePop_idx = [self.anc_props.index(max(self.anc_props))]
        else:
                new_anc_props = []
                [new_anc_props.append(label) for label in self.anc_props]
                choicePop_idx = [self.anc_props.index(max(self.anc_props))] # choose the highest
                new_anc_props.remove(self.anc_props[choicePop_idx[0]])
                choicePop_idx.append(self.anc_props.index(max(new_anc_props))) # choose the second highest
 
        if self.admix_multi_point:
                choicePop = []
                for idx in choicePop_idx:
                        choicePop.append(self.anc_label.index(self.last_gen_labels[idx])+1)
        else:
                choicePop = [int(i)+1 for i in choicePop_idx]

        for line in fileinput.input(seg_file):
                if fileinput.lineno() == 1:
                        fout.write(line)
                else:
                        data = line.rstrip().split("\t")
                        for i,seg in enumerate(data):
                                anc,endMarker = seg.split(":")
                                ind_seg_list = []
                                if i == 0:
                                        seg_range = [0,self.positions_list.index(int(endMarker))+1]
                                        prev_endMarker = endMarker
                                else:
                                        prev_endMarker_idx = self.positions_list.index(int(prev_endMarker))
                                        seg_range = [prev_endMarker_idx+1, self.positions_list.index(int(endMarker))+1]
                                        prev_endMarker = endMarker
                                for fit in fitness_index_list:
                                        if fit in range(seg_range[0],seg_range[1]):
                                                new_pop = choice(choicePop)
						# new_pop = np.random.choice(choicePop)
                                                new_seg = str(new_pop)+":"+str(endMarker)
                                                if new_seg not in ind_seg_list:
                                                        ind_seg_list.append(new_seg)
                                                else:
                                                        pass
                                if len(ind_seg_list) == 0:
                                        fout.write("".join(seg)+"\t")
                                else:
                                        fout.write("".join(new_seg)+"\t")
                        fout.write("\n")

    def file_to_dict(self):
	'''Reads a given file to a dictionary, with chromosome and position as the keys to the main dict and sub_dict resp.'''
       	self.return_dict = {}
	for idx,file in enumerate(self.parentalFiles):
		if len(file) == 1:
			dict_name = {}
			for line in fileinput.input(file):
				data = line.split()
				dict_name[fileinput.lineno()] = data[0]
			self.return_dict[idx+1] = [dict_name]
		else: # cases and controls
			sub_file_list = []
			for sub_file in file:
				dict_name = {}
				for line in fileinput.input(sub_file):
                                	data = line.split()
                                	dict_name[fileinput.lineno()] = data[0]
				sub_file_list.append(dict_name) 
                        self.return_dict[idx+1] = sub_file_list

    def chooseCaseInd(self,snp_no,chrom,hap_no,anc):
	chrm_risk_snp_no = [] ; no_dict = {}
	for chr in self.risk_snp_dict:
		chrm_risk_snp_no.extend(self.risk_snp_dict[chr])
	for i,no in enumerate(self.risk_snp_no):
		no_dict[no] = chrm_risk_snp_no[i]
	
	idx = [n for n in no_dict if no_dict[n] == snp_no][0]

	if int(self.simulated_disease_alleles2[chrom][snp_no][hap_no]) == 1:
		if len(self.separated_haps_dict[anc][idx][1]) > 0:
			ind = choice(self.separated_haps_dict[anc][idx][1])
		else:
			ind = choice(self.separated_haps_dict[anc][idx][0])
	else:
		if len(self.separated_haps_dict[anc][idx][0]) > 0:
                        ind = choice(self.separated_haps_dict[anc][idx][0])
                else:
                        ind = choice(self.separated_haps_dict[anc][idx][1])
	return ind

    def getFitIndDict(self):
	admix_fitness_dict = {}
	for anc in self.return_dict:
		simulated_haps = self.return_dict[anc][0]
		fitness_dict = self.get_inds_fitness(simulated_haps, False)
		max_fit = max(fitness_dict.keys())
		admix_fitness_dict[anc] = fitness_dict[max_fit]
	return admix_fitness_dict

    def getAdmixSimulatedDsl(self,pop):
	if pop == 'control':
		dsl = self.getDsl(admix=True)
		self.control_risk_snp = {} ; self.cases_risk_snp = {}
		for h in dsl:
			self.control_risk_snp[h] = dsl[h][:self.admix_controls] ; self.cases_risk_snp[h] = dsl[h][self.admix_controls:]
		simulated_disease_alleles = self.control_risk_snp
	else:
		simulated_disease_alleles = self.cases_risk_snp

        actual_risk_no = simulated_disease_alleles.keys()
        counter = 0 ; simulated_disease_alleles2 = {}
        for chrm in self.risk_snp_dict:
        	simulated_disease_alleles2[chrm] = {}
       		for ps in self.risk_snp_dict[chrm]:
                	simulated_disease_alleles2[chrm][ps] = simulated_disease_alleles[actual_risk_no[counter]]
               		counter = counter + 1
	return simulated_disease_alleles2

    def getAdmixFitnessDict(self):
	'''The function arranges the positions under selection. origin chromosome:actual index in the chromosome'''
	counter = 0
        pos_no_list = self.fitness_dict.keys()
      	self.selection_dict = {}
        for chrom in self.getMapPos:
        	self.selection_dict[chrom] = {}
        	for pos in self.selection_positions:
                	if pos in self.getMapPos[chrom]:
                        	self.selection_dict[chrom][self.getMapPos[chrom].index(pos)] = self.fitness_dict[pos_no_list[counter]]
                        	counter = counter + 1	

    def getDiseaseAncs(self):
	'''Obtain the parental population simulated under a disease causal model and the respective ancestry proportion'''
	
	diseased_anc_causal = [] ; diseased_anc_props_causal = [] 
	diseased_anc_null = [] ; diseased_anc_props_null = []
	for anc in self.pre_admix_disease_pop:
		anc_idx = self.anc_label.index(anc)

		if self.pre_admix_disease_pop[anc] == 'causal_model':
			diseased_anc_causal.append(anc_idx+1) 
			diseased_anc_props_causal.append(self.anc_props[anc_idx])
		else:
			diseased_anc_null.append(anc_idx+1)
                        diseased_anc_props_null.append(self.anc_props[anc_idx])
	
	if len(diseased_anc_causal) == 0:
		return diseased_anc_null,diseased_anc_props_null
	else:
		return diseased_anc_causal,diseased_anc_props_causal

    def getAncsUnderSelection(self):
	'''Obtain the ancestry of the parental population simulated under selection and the coresponding ancestry proportion'''
	
	ancs_under_selection = [] ; selection_ancs_props = []
	for anc in self.pre_admix_select_pop:
        	anc_idx = self.anc_label.index(anc)
                ancs_under_selection.append(anc_idx+1)
                selection_ancs_props.append(self.anc_props[anc_idx])

	return ancs_under_selection,selection_ancs_props

    def segment_genotype_sampling(self,seg_file,genotype_outfile,file_idx,post_admix):
	'''The function reads the seg_file and copy the corresponding genotype from the parental population genotype file'''

	logger.info("\nWriting genotypes of file %s to %s ... \n\n"%(str(seg_file),str(genotype_outfile)))

	fout = open(genotype_outfile,"wt")
	
	# create  a dictionary to store discarded individuals while writing
	discard_dict = {}
	for idx in range(len(self.parentalFiles)):
		discard_dict[idx+1] = []
	
	# determine whether the segments are controls segments or case segments and 
	# obtain the alleles for the risk SNPs, parental population simulated with disease and corresponding ancestry 
	if self.post_admix_disease:
		if self.admix_causal:
			if file_idx == 0:
				pop = "control"
			else:
				pop = "case"
			self.simulated_disease_alleles2 = self.getAdmixSimulatedDsl(pop)
			diseased_anc,diseased_anc_props = self.getDiseaseAncs()

        # obtain the alleles of SNPs simulated under selection, and the parental population simulated with selection, and 
	# the corresponding ancestry. Done for pre-admixture selection
	if self.selection:
		self.getAdmixFitnessDict()
		if post_admix:
			admix_fitness_dict = self.getFitIndDict()
		ancs_under_selection,selection_ancs_props = self.getAncsUnderSelection()

	# Read ancestry segment file, per haplotype, per chromsome and write corresponding genotype information.
	for line in fileinput.input(seg_file):
		
		# Initial conditions
		chrom_list = self.getMapPos.keys()
		chrom_list.sort()
		no_of_chrom = len(chrom_list)
		
        	chrom = chrom_list[0]
        	chrom_positions =  self.getMapPos[chrom]
		len_chrom = int(len(chrom_positions))
		len_prev_chrom  = 0
		curr_marker_chrom =  -1
		length_seg = 0
		
		# Readind start from second line, first line is header
		if fileinput.lineno() > 1:
			data = line.split()
			counter = 0
			if fileinput.lineno()%2 == 0: fit_h_dict = {}
			for seg in data:
				# set the segment length(len_seg), the index of the startMarker and endMarker
				prev_marker_chrom = curr_marker_chrom	
				curr_marker_chrom = int(chrom)
				anc,endMarker = seg.split(":") ; anc_idx = int(anc)
				if prev_marker_chrom != curr_marker_chrom:
					startMarker_idx = 0
				else:
					startMarker_idx = chrom_positions.index(int(prev_endMarker))+1
				endMarker_idx = chrom_positions.index(int(endMarker))+1

				len_seg = len(chrom_positions[startMarker_idx:endMarker_idx])
				
				# Empty the discard list if all individuals are used.
				if len(discard_dict[anc_idx]) == len(self.return_dict[anc_idx][file_idx]):
					discard_dict[anc_idx] = []
				
				# Randomly choose a new individual not in the discard dict from respective ancestry
				rand_ind = choice(self.return_dict[anc_idx][file_idx].keys())
				while rand_ind in discard_dict[anc_idx]:
					rand_ind = choice(self.return_dict[anc_idx][file_idx].keys())
				discard_dict[anc_idx].append(rand_ind)
 	
				# Obtain genotype for segment with a SNP under selection (post-admix selection)
				condition = 'choosing_seg'


				# obtain the copying start position and ancestry, for a segment with a SNP under selection
				# (pre-admixture selection)
				if self.selection and condition == 'choosing_seg':
					if self.pre_admixed_selection:
                                               	for pos_idx in self.selection_dict[chrom]:
                                                       	if pos_idx in range(startMarker_idx,endMarker_idx):
                                                               	startpos = startMarker_idx+len_prev_chrom
								anc_idx = choice(ancs_under_selection)
                                                               	condition = 'seg_chosen'
                                                               	break

				# obtain the copying start position and ancestry, for a segment with a risk SNP under a causal model
				if self.disease and condition == 'choosing_seg':
                                	if self.admix_causal:
                                        	for snp_no in self.risk_snp_dict[chrom]:
                                                	if snp_no in range(startMarker_idx,endMarker_idx):
                                                        	startpos = startMarker_idx + len_prev_chrom
                                                                anc_idx = choice(diseased_anc)
                                                                condition = 'seg_chosen'

				# obtain the copying start position if the segment does not contains a SNP either under selection 
				# or risk SNP. 
				if condition == 'choosing_seg':
					limit = (len_prev_chrom+len_chrom) - len_seg
					startpos = randint(len_prev_chrom,limit)
				
				endpos = startpos + len_seg
				geno = list(self.return_dict[anc_idx][file_idx][rand_ind])[startpos:endpos]
				fout.writelines("".join(geno))
				
				if int(endMarker) == int(self.getMapPos[chrom][-1]):
                                	if len(chrom_list) > 1 and chrom != chrom_list[-1]:
                                        	len_prev_chrom = len_prev_chrom + len_chrom
                                                counter += 1
                                                chrom = chrom_list[counter]
                                                chrom_positions = self.getMapPos[chrom]
                                                len_chrom = int(len(chrom_positions))
                                        else:
                                                pass
				prev_endMarker = endMarker
                        	length_seg = length_seg+len_seg
			fout.writelines("\n")
	del discard_dict			   
	fout.close()

    def separateHaps(self,partitioned_haps_dict,update):
	'''partitions the haplotypes depending on the allele at the positions under selection, to allow easy choice of the haplotype to be 
	used at a segment with risk SNPs'''

	if not update:
		inversed_haps_dict = self.invert_dict(self.migrated_haps_dict,False) 
		for snp_no in inversed_haps_dict:
        		if (snp_no-1) in self.fitness_dict:
                		snp = inversed_haps_dict[snp_no]
                		if len(snp) == 1:
                        		snp = list(inversed_haps_dict[snp_no])

                        	list0 = [] ; list1 = []
				for hap_no in range(len(snp)):
                        		if snp[hap_no] == str(0):
                                		list0.append(hap_no+1)
                             		else:
                                		list1.append(hap_no+1)
             			partitioned_haps_dict[snp_no-1] = [list0,list1]
			else:
                		pass
	else:
		haps_no = self.migrated_haps_dict.keys() ; haps_no.sort()
		last_hap_no = haps_no[-1]
		added_hap = self.migrated_haps_dict[last_hap_no]
		for snp_no in self.fitness_dict:
			if added_hap[snp_no] == str(0):
				partitioned_haps_dict[snp_no][0].append(last_hap_no)
			else:
				partitioned_haps_dict[snp_no][1].append(last_hap_no)
					
        return partitioned_haps_dict

    def getFitSeg(self,st,cut,hap_no,simulated_haps,partitioned_haps,simulated_snps,snp_no_list):
        '''The function selects a segment from the sample hap, if the segment has a disease risk SNP in absence of selection'''
        hap_seg = []

        # When only 1 snp, choose a h that satisfy the simulated snp if found, else "mutate" a random h
        if len(snp_no_list) == 1:
                snp_no = snp_no_list[0]
                #print snp_no,hap_no, len(simulated_snps[snp_no])
                if simulated_snps[snp_no][hap_no] == str(0):
                        if len(partitioned_haps[snp_no][0]) == 0:
                                h = choice(partitioned_haps[snp_no][1])
                                #h = np.random.choice(partitioned_haps[snp_no][1])
                                choosen_hap = list(simulated_haps[h])
                                choosen_hap[snp_no] = "0"
                                hap_seg = choosen_hap[st:cut]
                        else:
                                h = choice(partitioned_haps[snp_no][0])
				 #h = np.random.choice(partitioned_haps[snp_no][0])
                                hap_seg = list(simulated_haps[h])[st:cut]
                else:
                        if len(partitioned_haps[snp_no][1]) == 0:
                                h = choice(partitioned_haps[snp_no][0])
                                #h = np.random.choice(partitioned_haps[snp_no][0])
                                choosen_hap = list(simulated_haps[h])
                                choosen_hap[snp_no] = "1"
                                hap_seg = choosen_hap[st:cut]
                        else:
                                h = choice(partitioned_haps[snp_no][1])
                                #h = np.random.choice(partitioned_haps[snp_no][1])
                                hap_seg = list(simulated_haps[h])[st:cut]
                prev_h = h
        else:
                # Option1: find a h that satisfies all the simulated snps
 		common_haps = {}
                for snp_counter,snp_no in enumerate(snp_no_list):
                	if snp_counter == 0:
                        	if simulated_snps[snp_no][hap_no] == str(0):
                                	prev_list = partitioned_haps[snp_no][0]
                              	else:
                                	prev_list = partitioned_haps[snp_no][1]

                        else:
                        	if simulated_snps[snp_no][hap_no] == str(0):
                                	prev_list = list( set(prev_list) & set(partitioned_haps[snp_no][0]) )
                             	else:
                                	prev_list = list( set(prev_list) & set(partitioned_haps[snp_no][1]) )
                                if len(prev_list) > 0:
                                	common_haps[snp_counter] = prev_list

		if len(prev_list) > 0:
                	h = choice(prev_list)		
		        #h = np.random.choice(prev_list)
                	hap_seg.extend(list(simulated_haps[h])[st:cut])
                        prev_h = h
		else:
                	# Option 2: Choose a h that satisfies atleast two of the simulated snps (in order)
                        if len(common_haps) != 0:
                              	h1 = choice(common_haps[max(common_haps.keys())])
                               	#h1 = np.random.choice(common_haps[max(common_haps.keys())])
                                hap_seg.extend(list(simulated_haps[h1]))
                                prev_h = h1
			else:
				# Option 3: Choose a random h and mutate the required positions
				h1 = choice(simulated_haps)
				prev_h = h1
				choosen_hap = list(simulated_haps[h1])
				for s in snp_no_list:
					if simulated_snps[snp_no][hap_no] == str(0):
                                		choosen_hap[snp_no] = "0"
					else:
						choosen_hap[snp_no] = "1"
				hap_seg.extend(choosen_hap)
	return prev_h,hap_seg
	
	
    def getMultiPointParams(self):
	'''Adjust the self.parentalFiles order and self.anc_props to match the continuous order of the multipoint parameter file'''
	file_anc_labels = [] ; param_dict = {}
	for line in fileinput.input(self.multipoint_parafile):
		data = line.split()
		if fileinput.lineno() == 1:
			file_anc_labels.extend(data[1:])

		elif fileinput.lineno() == 2:
			param_dict[fileinput.lineno()-1] = data[1:]
		else:
			new_props = []
			for i,anc in enumerate(file_anc_labels):
				if float(param_dict[fileinput.lineno()-2][i]) > 0:
					prop = (float(param_dict[fileinput.lineno()-2][i])/1.0)*float(data[1])
		 			new_props.append(prop)
				else:
					new_props.append(float(data[i+1]))
			param_dict[fileinput.lineno()-2] = new_props
		last_key = fileinput.lineno()-2
	self.anc_props =  param_dict[last_key]
	ordered_parentalFiles = list(np.zeros(len(self.parentalFiles)))
	for i,file in enumerate(self.parentalFiles):
		key = self.parentalFiles[i][0].split("/")[-1].split(".")[0]
		file_idx = self.anc_label.index(key)
		ordered_parentalFiles[file_idx] = file
	self.parentalFiles = ordered_parentalFiles

    def simulatePostAdmixSelection(self,segs_file,haps_file,segs_outfile,haps_outfile):
	'''Simulate the segments for post_admixture selection, this are created through recombination of the 'migrated' set'''
	
	segs_out = open(segs_outfile,"wt")
	haps_out = open(haps_outfile,"wt")

	#Randomly choose 1/3 of the population 'migrated set' (odd index of the haplotypes)
	initial_sample_size = ((self.no_to_simulate/2)/3)
	random_sample = sample(range(1,self.no_to_simulate,2),initial_sample_size)
	
	#Extract corresponding haplotype segments, in pairs
	segment_dict = {} ; enum = 1 ; haps_list = []
	for line in fileinput.input(segs_file):
		if fileinput.lineno() == 1:
			segs_out.writelines(line)
		if fileinput.lineno() > 1:
			idx = fileinput.lineno()-1
			if idx in random_sample or idx-1 in random_sample :
				haps_list.append(idx)
				data = line.split()
				segment_dict[enum] = data
				enum = enum + 1
	fileinput.close()
	haps_sample_size = len(segment_dict)

	#Extract corresponding haplotypes from the initial admixture population
	self.migrated_haps_dict = {} ; enum2 = 1
	for line1 in fileinput.input(haps_file):
		data1 = list(line1.rstrip())
		if fileinput.lineno() in haps_list:
			self.migrated_haps_dict[enum2] = data1
			enum2 = enum2 + 1
	fileinput.close()

	#separate the generated haplotype set depending on the allele at the SNPs under selection
	partitioned_haps = {}
	partitioned_haps = self.separateHaps(partitioned_haps,update = False)
	
	#Obtain the simulated alleles per SNP under selection
	sim_fit_allele = self.getFitSNPs(self.post_admix_sim_size/2)
	
        
	chrom_list = self.getMapPos.keys()
        chrom_list.sort()
        for m in range(haps_sample_size+1, haps_sample_size+self.post_admix_sim_size+1):
		hap_idx = m - (haps_sample_size+1) 
		counter = 0
        	chrom = chrom_list[counter]
		break_points_list = self.recomb_points(m)
		startpoint_idx = 0
		segs_list = []
		chrom_positions = self.getMapPos[chrom]; chrom_positions.sort()
		last_chrom_pos = chrom_positions[-1]
		prev_point = 0
		haps = []
		new_ind_segs = []
		haps_list = self.migrated_haps_dict.keys()
		for point_idx in break_points_list:
			endpoint_idx = int(point_idx)
			point = self.positions_list[point_idx-1]
			if int(point) not in chrom_positions:
				pass
			else:
				fit_list = []
				for fit_pos in self.fitness_dict:
					if fit_pos in range(startpoint_idx,endpoint_idx):
						fit_list.append(fit_pos)
				if len(fit_list) != 0:
					rand_hap,geno_seg = self.getFitSeg(startpoint_idx,endpoint_idx,hap_idx,self.migrated_haps_dict,partitioned_haps,sim_fit_allele,fit_list)
				
				else:	
					rand_hap = choice(haps_list)
					geno_seg = self.migrated_haps_dict[rand_hap][startpoint_idx:endpoint_idx]

				haps.extend(geno_seg)
				haps_out.writelines(''.join(geno_seg))
				for seg in segment_dict[rand_hap]:
					anc,endMarker = seg.split(":")
					if int(endMarker) not in chrom_positions:
						pass
					else:
						if int(endMarker) >= int(point):
							sim_seg = str(anc)+":"+str(point)
							new_ind_segs.append(sim_seg)
							segs_out.writelines("\t"+sim_seg)
							segs_list.append(sim_seg)
							break
						else:
							if int(endMarker) > prev_point:
								sim_seg = str(anc)+":"+str(endMarker)
								new_ind_segs.append(sim_seg)
								segs_out.writelines("\t"+sim_seg)
								segs_list.append(sim_seg)
			prev_point = point
			startpoint_idx = int(point_idx)
			if int(point) == chrom_positions[-1] and point_idx != break_points_list[-1]:
				counter = counter + 1
				chrom = chrom_list[counter]
				chrom_positions = self.getMapPos[chrom]; chrom_positions.sort()
				prev_point = 0

		self.migrated_haps_dict[m] = haps
		segment_dict[m] = new_ind_segs
		partitioned_haps = self.separateHaps(partitioned_haps,update = True)
		segs_out.writelines("\n")
		haps_out.writelines("\n")
        segs_out.close()
	haps_out.close()
	del self.migrated_haps_dict
	del segment_dict
        logger.info("Post admixture segment simulation complete!\n\n")

    def merge_haps_segment(self,haps_segfile,geno_segfile):
	'''The function merges the haplotypes segments (haps_segfile) into individual segments(geno_segfile)'''

	outfile = open(geno_segfile,"wt")
	haps_segs = open(haps_segfile); read_input = haps_segs.read(); haps_segs.close(); input_segs = re.split ('\n', read_input)
	while "" in input_segs: input_segs.remove ('')

	header = input_segs[0]
	input_segs = input_segs[1:]

	no_of_haps = len(input_segs)

	chrom_list = self.getMapPos.keys()
	chrom_list.sort()

	outfile.writelines(header +"\n")
	for i in range(0,no_of_haps,2):
        	counter = 0
        	data1 = input_segs[i].split('\t')
		data2 = input_segs[i+1].split('\t')
        	while '' in data1: data1.remove('')
        	while '' in data2: data2.remove('')
        	curr_chrom = -1
        	for chrom in chrom_list:
                	out = [] ; position_list = [] ; anc_hap1 = [] ; anc_hap2 = [] ; pos_hap1 = [] ; pos_hap2 = []
                	prev_chrom = curr_chrom
                	curr_chrom = chrom
                	chrom_positions = self.getMapPos[chrom]
                	for seg1 in data1:
                        	anc,endMarker = seg1.split(':')
                        	if int(endMarker) not in chrom_positions:
                                	pass
                        	else:
					position_list.append(int(endMarker))
                                	endpoint = chrom_positions.index(int(endMarker))+1
                                	if prev_chrom != curr_chrom:
                                        	startpoint = 0
                                	else:
                                        	startpoint = chrom_positions.index(int(prev_endMarker))+1
                                	anc_hap1.append(anc) ; pos_hap1.append([startpoint,endpoint])
                                	prev_endMarker = endMarker
                	for seg2 in data2:
                        	anc,endMarker = seg2.split(':')
                        	if int(endMarker) not in chrom_positions:
                                	pass
                        	else:
                                	if int(endMarker) not in  position_list:
                                        	position_list.append(int(endMarker))
                                	endpoint = chrom_positions.index(int(endMarker))+1
                                	if prev_chrom != curr_chrom:
						startpoint = 0
                                	else:
                                        	startpoint = chrom_positions.index(int(prev_endMarker2))+1
                                	anc_hap2.append(anc) ; pos_hap2.append([startpoint,endpoint])
                                	prev_endMarker2 = endMarker
                	position_list.sort()
                	for pos in position_list:
                        	if int(pos) not in chrom_positions:
                                	pass
                        	else:
                                	for j,seg in enumerate(pos_hap1):
                                        	if int(pos) in chrom_positions[seg[0]:seg[1]]:
                                                	anc1 = anc_hap1[j]
                                                	break
                                	for k,seg2 in enumerate(pos_hap2):
                                        	if int(pos) in chrom_positions[seg2[0]:seg2[1]]:
                                                	anc2 = anc_hap2[k]
    							break
                                	if anc2 == anc1:
                                        	fin_anc = anc1
                                	else:
                                        	fin_anc = [anc1,anc2] ; fin_anc.sort()
                                        	fin_anc = ''.join(fin_anc)

                                	if pos == position_list[0]:
                                        	prev_anc = fin_anc
                                        	out.append([fin_anc,self.positions_list.index(pos)])
                                	else:
                                        	if int(prev_anc) == int(fin_anc):
                                                	out.pop(-1)
                                        	out.append([fin_anc,self.positions_list.index(pos)])
                                        	prev_anc = fin_anc
			outfile.writelines(''.join(['\t'+str(i[0])+':'+str(i[1]) for i in out]))

        	outfile.writelines("\n")
	outfile.close()

    def admixture_sampling(self):
	'''The function determines which segment file is required, and the opens the output files to write the sampled genotype, 
	then calls the necessary functions for genotype sampling'''
	sampling_starttime = datetime.datetime.today()
	if not self.disease:
		if self.admix_multi_point:
			seg_file = self.outfolder+"multipoint_ancestry_segments_haps.out"+str(self.generations[-1])
		else:
			seg_file = self.outfolder + "ancestry_segments_haps.out"
		if self.selection:
			if self.pre_admixed_selection:
				genotype_outfile = self.outfolder + "pre_admixture_selection.haps"
			else:
				genotype_outfile = self.outfolder + "admixture.haps"
		else:
			genotype_outfile = self.outfolder + "admixture.haps"
		self.segment_genotype_sampling(seg_file,genotype_outfile,0,False)
		geno_segfile = seg_file.split(".") ; geno_segfile.pop(-1)
		geno_segfile = ''.join(geno_segfile)+"_inds.out"
		self.merge_haps_segment(seg_file,geno_segfile)
		
		if self.selection:
			if self.post_admix_selection:
				post_seg_file = self.outfolder+"post_admixture_selection_segment_haps.out"
				post_genotype_outfile = self.outfolder + "post_admixture_selection.haps"
				self.simulatePostAdmixSelection(seg_file,genotype_outfile,post_seg_file,post_genotype_outfile)
				geno_segfile = post_seg_file.split(".") ; geno_segfile.pop(-1)
                                geno_segfile = ''.join(geno_segfile)+"_inds.out"
                                self.merge_haps_segment(post_seg_file,geno_segfile)

	else:
		if self.admix_null:
			seg_files = [self.outfolder +"admixture_controls_null_segments_haps.out",self.outfolder+"admixture_cases_null_segments_haps.out"]
			genotype_outfile = [self.outfolder+"admixture_controls_nullmodel.haps",self.outfolder+"admixture_cases_nullmodel.haps"]
		else:# self.admix_causal
			seg_files = [self.outfolder +"admixture_controls_causal_segments_haps.out",self.outfolder+"admixture_cases_causal_segments_haps.out"]
			genotype_outfile = [self.outfolder+"admixture_controls_causalmodel.haps",self.outfolder+"admixture_cases_causalmodel.haps"]
		for idx in range(len(seg_files)):
			self.segment_genotype_sampling(seg_files[idx],genotype_outfile[idx],idx,False)
			geno_segfile = seg_files[idx].split(".") ; geno_segfile.pop(-1)
                	geno_segfile = ''.join(geno_segfile)+"_inds.out"
                	self.merge_haps_segment(seg_files[idx],geno_segfile)
			if self.selection:
                        	if self.post_admix_selection:
					
					if idx == 0:
						post_seg_outfile = self.outfolder+"post_admixture_selection_controls_segments.haps"
						post_genotype_outfile = self.outfolder+"post_admixture_selection_controls.haps"
					else:
						post_seg_outfile = self.outfolder+"post_admixture_selection_cases_segments.haps"
						post_genotype_outfile = self.outfolder+"post_admixture_selection_cases.haps"
                                	self.simulatePostAdmixSelection(seg_files[idx],genotype_outfile[idx],post_seg_outfile,post_genotype_outfile)
					geno_segfile = post_seg_outfile.split(".") ; geno_segfile.pop(-1)
                        		geno_segfile = ''.join(geno_segfile)+"_inds.out"
                			self.merge_haps_segment(post_seg_outfile,geno_segfile)
					
	sampling_endtime = datetime.datetime.today()
	logger.info("Writing genotypes took %s time ...\n\n"%str(sampling_endtime - sampling_starttime))


    def admixture_simulation(self):
        ''' Main function. 
        1. This function load the reference genotype data and check consistancy with the provided SNP file. 
        2. Call the ancestry segment [self.genome_recordSegment()] function to process admixture scenario based on givenn number of generation$
        3. Call the function for drawing the genotypes of admixed popualtion [sampling_simulation_genotype()] based on ancestry segment file a$
        '''

        self.getMapPos = {}; getMarker = []; ref_file = self.parentalFiles ; self.risk_snp_dict = {}
        for line in fileinput.input(self.snpFile):
                data = line.split()
                getMarker.append(data[0])
                if fileinput.lineno() == 1:
                        self.getMapPos[int(data[1])] = [int(data[3])]
                else:
                        if int(data[1]) in self.getMapPos:
                                self.getMapPos[int(data[1])].append(int(data[3]))
                        else:
                                self.getMapPos[int(data[1])] = [int(data[3])]
	if self.disease:
		for chrm in self.getMapPos:
			self.risk_snp_dict[chrm] = [] 
			for pos in self.disease_snp_pos:
				if pos in self.getMapPos[chrm]:
					self.risk_snp_dict[chrm].append(self.getMapPos[chrm].index(pos))
								

       	# check that the number of SNPs in position file equals no of snps in the references file.
	for filename in ref_file:
                file_index = ref_file.index(filename);reference_data = {}
                for line in fileinput.input(filename):
                	data = line.split()
                        try:
                                reference_data[getMarker[fileinput.lineno()-1]] = data[0]
                        except (TypeError,NameError,ValueError,IndexError):
                                sys.stderr.write('ERROR ==> Failed to process the input, reference data do not match the position file\n')
                                self.terminate()
		if len(reference_data) != len(getMarker):
			sys.stderr.write('ERROR ==> Failed to process the input, reference data %s do not match the position file\n\n'%filename)
                	self.terminate()
                reference_data.clear()

        #Record ancestry segment, save them into a file and sample the simulation genotype data from the saved ancestry segment file
        self.check_admixture_parameter()
        self.evolve_parental_pops()
	if self.admix_single_point:
		if self.selection:
			self.genome_recordSegment(self.generations)
			single_point_seg_file = self.outfolder+"ancestry_segments_haps.out"
			self.sim_selection_segments(single_point_seg_file)
			system("mv"+" "+self.outfolder+"ancestry_selection_segments_haps.out"+" "+single_point_seg_file)
		else:
			self.genome_recordSegment(self.generations)

		self.file_to_dict()
                if self.disease:
			self.split_segs()
			if self.admix_causal:
				self.get_case_segments()
				system("mv"+" "+ self.outfolder+"admixture_cases_null_segments_haps.out"+" "+self.outfolder+"admixture_cases_causal_segments_haps.out")

		self.admixture_sampling()

        if self.admix_multi_point:
		self.simulate_multipoint_admixture_segs()
		self.file_to_dict()
		self.getMultiPointParams()
                if self.disease:
                        self.split_segs()
			if self.admix_causal:
                                self.get_case_segments()

		self.admixture_sampling()
        logger.info("simulation of admixed population complete!\n\n")

class fractalsim(fractInit,OtherUtlis,statUtils,Homogeneous_Simulator,case_controlSimulator,Admixture):
    '''Performs the overrall steps of FractalSIM as described in the ancGWAS method (Chimusa et al 2013)'''

    def fractal(self):
        '''        
                Running fractalsim
        '''
        starting_time = datetime.datetime.today()
        logger.info('Starting at time:%s' % str(starting_time))
        self.res = fractalsim(argv)
        logger.info("Loading parameters from %s ..."%os.path.abspath(argv))
        logger.info("Options in effect:")
        for param in sorted(self.res.Params):
            logger.info('             '+str(self.res.Params[param][0])+': '+str(self.res.Params[param][1]))
        print ''
	
	self.check_other_params()
	if self.homogeneous:
		if not self.disease:
			self.homogeneous_simulation(self.hapFile)
		else:
			self.case_control_homogeneous_simulation(False)
	if self.admixed:
		self.admixture_simulation()
	
        currenttime1 = time.asctime(time.localtime())
        finishing_time = datetime.datetime.today()
        logger.info("\n\nFinish at time:%s"%str(finishing_time))
        logger.info("Simulation took %s time to finish.\n\n\n"%str(finishing_time-starting_time))
        self.res.terminate()

if __name__ == '__main__':
    try:
        global gpv1
        argv = sys.argv[1]
    except IndexError:
        argv = ''
    finally:
    	run = fractalsim(argv)
        run.fractal()



