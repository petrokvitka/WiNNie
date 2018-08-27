
"""
WiNNie - Weighted motif eNrichmeNt analysis uses the weights received from TOBIAS in bigwig files to find enriched motifs
@author: Anastasiia Petrova
@contact: anastasiia.petrova(at)mpi-bn.mpg.de
"""

import argparse
import sys
import os
import re
import time
import multiprocessing
import logging
import subprocess
from Bio import SeqIO
import Bio.SeqIO.FastaIO as bio
import numpy as np
from collections import defaultdict
from scipy import stats
import pyBigWig
from statsmodels.sandbox.stats.multicomp import multipletests #for bonfferoni
import matplotlib.pyplot as plt
import random

import textwrap


import MOODS.scan
import MOODS.tools
import MOODS.parsers

logger = logging.getLogger('my_ame')
logger.setLevel(logging.INFO)

formatter = logging.Formatter("%(asctime)s : %(message)s", "%Y-%m-%d %H:%M")

fh = logging.FileHandler('winnie_log.txt')
fh.setLevel(logging.INFO)
fh.setFormatter(formatter)
logger.addHandler(fh)

ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)
logger.addHandler(ch)


#catch all the information about input and output files as well as information on the used tool (fimo or moods)
def parse_args():
	
	parser = argparse.ArgumentParser(prog = 'winnie', description = textwrap.dedent('''                           

		This script takes a list of motifs loaded from jaspar.genereg.net as a combined text file in .MEME or .PFM format, two bigWig-files containing scores for two different conditions, a genome file in FASTA format and optionaly a .bed file with regions of interest as input. The output is a file containing information about the mean and standard deviation of the data and enriched motifs sorted by p value. Please note, if you want to have all intermediate output files, enter --clean nothing.
		'''), epilog='That is what you need to make this script work for you. Enjoy it')
	
	required_arguments = parser.add_argument_group('required arguments')
	required_arguments.add_argument('-m', '--motifs', help='file in .MEME or .PFM format with mofits loaded from jaspar.genereg.net', required=True)
	required_arguments.add_argument('-g', '--genome', help='a whole genome file or regions of interest in FASTA format to be scanned with motifs', required=True)
	required_arguments.add_argument('-c1', '--condition1', help='a bigWig-file with scores for the first condition', required=True)
	required_arguments.add_argument('-c2', '--condition2', help='a bigWig-file with scores for the second condition', required=True)

	#all other arguments are optional
	parser.add_argument('-o', '--output_directory',  default='output', const='output', nargs='?', help='output directory, default ./output/')
	parser.add_argument('-b', '--bed_file',  nargs='?', help='a .bed file to be merged with the whole genome file to find regions of interest')
	parser.add_argument('--clean', nargs='*', choices=['nothing', 'all', 'cut_motifs', 'merge_output'], dest='cleans', help='choose the files you want to delete from the output directory, the default is deleting all the temporary files from the directory', default=['all'])
	parser.add_argument('--cores', type=int, help='number of cores allowed to use by this tool, by default the tool uses 2 cores', default=2)
	parser.add_argument('-p', '--p_value', type=float, help='enter the p value, the default p value is 1e-4', default=0.0001)
	parser.add_argument('--silent', action='store_true', help='while working with data write the information only into ./log.txt')
	parser.add_argument('--moods_bg', nargs='+', type=float, help='set the bg for moods, by default moods uses the bg is 0.25 0.25 0.25 0.25', default=[0.25, 0.25, 0.25, 0.25])
	parser.add_argument('--output_number', type=int, help='number of enriched motifs that should be written to output, by default 40', default=40)
	args = parser.parse_args()

	return args

def check_directory(directory):
	if not os.path.exists(directory):
		os.makedirs(directory)
		logger.info('a new directory ' + directory + ' was created')
					
#merge the whole genome with the regions mentioned in .bed file
def merge(genome, bed_file, output_directory):
	logger.info('the merging of files ' + genome + ' and ' + bed_file + ' will end soon, the result file is output_merge.fa')
	output_merge = os.path.join(output_directory, "output_merge.fa")
	if os.path.isfile(output_merge):
		output_merge = os.path.join(output_directory, "custom_output_merge.fa")
	os.system("bedtools getfasta -fi " + genome + " -bed " + bed_file + " -fo " + output_merge)
	return output_merge

#split the motifs each in other file
def split_motifs(motifs, output_directory):
	logger.info("the file with motifs " + motifs + " will be checked for motifs and if needed splitted in files each containing only one motif")
	
	first_line = subprocess.getoutput("head -1 " + motifs) #find the first line of the input file

	if first_line.startswith(">"):
		#the motif file probably has the .pfm format, try to read and split it
		splitted_motifs = read_pfm(motifs, output_directory)
	else: #maybe the file with motifs is in MEME format, so try to convert it
		logger.info("the input file has not the expected format, I will try to convert it to .pfm format")
		splitted_motifs = convert_meme_to_pfm(motifs, output_directory)

	return splitted_motifs

def read_pfm(motifs, output_directory):
	splitted_motifs = [] #to save the names of files after splitting
	motif = [] #to save the motif itself, which will be written to the file

	with open(motifs) as read_file:
		lines = 0
		for line in read_file:
			#as the motif has first line with the name and 4 lines with information, if the 5th line is something else than the name of the next motif, the exit will be forced
			if lines == 5 and not line.startswith(">"):
				logger.info('please make sure that the file with motifs has a right format and the number of lines is right in the motif file')
				sys.exit()
			else:
				if line.startswith(">"):
					if 'written_file' in locals():
						written_file.write(''.join(motif))
						motif = []
						lines = 0
						written_file.close()

					motif_alternate_name = check_name(re.split(' ', line)[1].rstrip())
					motif_id = re.split(' ', line[1:])[0] #[1:] meands do not use the first character
					motif_name = os.path.join(output_directory, motif_alternate_name + '_' + motif_id + '.pfm')
							
					splitted_motifs.append(motif_name)
					written_file = open(motif_name, 'w')
						
			if lines >= 1 and lines <= 4: #one motif has 5 lines, the first consists the name, the next 4 - the information we need to proceed the data within moods
				motif.append(line)
					
			lines = lines + 1
	written_file.write(''.join(motif))
	written_file.close()

	return splitted_motifs

def convert_meme_to_pfm(motifs, output_directory):
	#i can only convert the file to pfm if the motifs file is in MEME format

	splitted_motifs = [] #to save the names of files after splitting
	rows = [[] for row in range(4)]

	with open(motifs) as read_file:
		lines = 0
		for line in read_file:
			if lines == 0 and not line.startswith("MEME version"):
				logger.info('please make sure that the file with motifs has a right format and the number of lines is right in the motif file')
				sys.exit()
			else:
				#search for motifs and save each to another file
				if line.startswith("MOTIF"):
					
					if 'written_file' in locals():
						for row in rows:
							written_file.write('\t'.join(row) + '\n')

						rows = [[] for row in range(4)]

						written_file.close()
					
					#the alternate name will be checked for validity and the invalid chars will be replaced with '_'
					if len(re.split(' ', line.rstrip())) == 3: #in the input motif file the motif id and the alternate name are splitted using the tab, otherwise they are splitted using _, but we do not want to change it if so
						motif_alternate_name = check_name(re.split(' ', line)[2].rstrip())
						motif_id = re.split(' ', line)[1]
						motif_name = os.path.join(output_directory, motif_alternate_name + '_' + motif_id + '.pfm')

					else: 
						motif_alternate_name = check_name(re.split(' ', line)[1].rstrip())
						motif_name = os.path.join(output_directory, motif_alternate_name + '.pfm')
					
					#make a list with all the motif names to know which files to iterate when fimo is called
					splitted_motifs.append(motif_name)

					written_file = open(motif_name, 'w')
					
				elif line.startswith("letter-probability matrix"):
					columns = int(re.split(' ', re.split('w= ', line)[1])[0]) #find the number of columns from the line out of the motifs file
					nsites = int(re.split(' ', re.split('nsites= ', line)[1])[0]) #find the nsites to count the frequency count for .pfm file

				elif line.startswith(' '): #each line with information about frequency starts in MEME format with ' '
					for i in range(len(rows)):
						rows[i].append(str(round(float(re.findall(r'\S+', line)[i])*nsites))) #split the line, do not mention how much whitespaces are in between, multiply it with nsites and save it to the corresponding row

			lines = lines + 1
					
		#write the last motif
		for row in rows:
			written_file.write('\t'.join(row) + '\n')		

		written_file.close()
	read_file.close()

	return splitted_motifs

#if there are chars that are not allowed, they will be replaced with '_', to the possibly invalid names there will be added '_' at the beginning of the name
def check_name(name_to_test):
	badchars= re.compile(r'[^A-Za-z0-9_. ]+|^\.|\.$|^ | $|^$')
	badnames= re.compile(r'(aux|com[1-9]|con|lpt[1-9]|prn)(\.|$)')

	#replace all the chars that are not allowed with '_'
	name = badchars.sub('_', name_to_test)

	#check for the reserved by the os names
	if badnames.match(name):
		name = '_' + name
	return name

def call_moods(one_motif, genome, output_directory, p_value, moods_bg, fisher_dict_bg, fisher_dict, control, overexpression, control_dict, overexpression_dict, differences):

	#check if this is a bigwig file
	bw_control = pyBigWig.open(control)
	bw_overexpression = pyBigWig.open(overexpression)
	if not bw_control.isBigWig() or not bw_overexpression.isBigWig(): # <---------------- if the function to make random score stays, delete this if!!!
		logger.info("please provide the bigwig file!")
		sys.exit()
	else:
		# prepare everything for moods
		# setting standard parameters for moods
		pseudocount = 0.0001

		bg = tuple(moods_bg)

		#logger.info("moods will work with the p_value " + str(p_value) + " and the bg " + str(bg))

		motif_name = os.path.basename(one_motif)

		matrix_names = [os.path.basename(one_motif)]

		matrices = []
		matrices_rc = []

		valid, matrix = pfm_to_log_odds(one_motif, bg, pseudocount)

		key_for_bed_dict = ''

		if valid:

			#logger.info("please be patient, moods is working on the data")

			matrices.append(matrix)
			matrices_rc.append(MOODS.tools.reverse_complement(matrix,4))
			matrices_all = matrices + matrices_rc
			thresholds = [MOODS.tools.threshold_from_p(m, bg, p_value, 4) for m in matrices_all]

			scanner = MOODS.scan.Scanner(7)
			scanner.set_motifs(matrices_all, bg, thresholds)

			with open(genome) as handle:

				seq_iterator = bio.SimpleFastaParser(handle)

				for header, seq in seq_iterator:				

					header_splitted = re.split(r':', header)

					if len(header_splitted) == 1: #if there are no positions given
						header = header + ":0-" #set the first position as 0 and split it once more
						header_splitted = re.split(r':', header)
						logger.info("moods works with " + header)
					else: #the given genome file is a file with peaks, so use the header of the peak as a key to search in the bed dictionary for additional information later on
						key_for_bed_dict = header

						#------------------------delete?-------------------------------------
						if not key_for_bed_dict in fisher_dict_bg.keys():
							fisher_dict_bg[key_for_bed_dict] = 0
							fisher_dict[key_for_bed_dict] = 0
							#peak_length = int(re.split(r'-', header_splitted[1])[1]) - int(re.split(r'-', header_splitted[1])[0]) #find the length of the peak
							#fisher_dict[key_for_bed_dict] = np.zeros(peak_length) #initialize array of 0.0 with the length of the peak
						#----------------------------------------------------------------

					chromosom = header_splitted[0]
					positions = re.split(r'-', header_splitted[-1])

					#compute the background difference for this peak
					#bw_global_score_control = np.mean(np.nan_to_num(np.array(list(bw_control.values(chromosom, int(positions[0]), int(positions[1]))))))
					#bw_global_score_overexpression = np.mean(np.nan_to_num(np.array(list(bw_overexpression.values(chromosom, int(positions[0]), int(positions[1]))))))
					#
					#bw_global_difference = bw_global_score_control - bw_global_score_overexpression
					#print(bw_global_difference)

					results = scanner.scan(seq)

					fr = results[:len(matrix_names)] #forward strand
					rr = results[len(matrix_names):] #reverse strand

					results = [[(r.pos, r.score, '+', ()) for r in fr[i]] + 
						[(r.pos, r.score, '-', ()) for r in rr[i]] for i in range(len(matrix_names))] #use + and - to indicate strand

					for (matrix, matrix_name, result) in zip(matrices, matrix_names, results):

						motif_id = re.split(r'_', matrix_name)[-1].replace(".pfm", '') #find the id of the given morif
						motif_alternate_name = matrix_name.replace(motif_id, '')[:-1] #the alternate name of the motif is the name of the file without id and with cutted last character, that is _

						if len(matrix) == 4:
							l = len(matrix[0])
						if len(matrix) == 16:
							l = len(matrix[0] + 1)
						for r in sorted(result, key=lambda r: r[0]):
							strand = r[2]
							pos = r[0]
							hitseq = seq[pos:pos+l] #sequence
							score = format(r[1], '.15f') #round to 15 digits after floating point, already type str

							if key_for_bed_dict != '':
								start = pos + 1
								end = pos + len(hitseq)
								#chromosom = key_for_bed_dict #instead of only the name of chromosom write the key to search in the bed_file					
							else:
								start = int(positions[0]) + pos + 1
								end = start + len(hitseq) - 1

							"""# ---------------------------------------delete??-------------------------------
							if key_for_bed_dict not in scores_dict: #if this peak has no matches yet
								#fill the peak with 0
								peak_length = int(positions[1]) - int(positions[0]) + 1
								scores_dict[key_for_bed_dict] = np.full(peak_length, 0.0)
							"""#-------------------------------------------------------------------------------

							#find the real start and end positions on the chromosom		
							real_start = int(positions[0]) + int(start) #start of the peak + start of the motif within the peak, do not add 1, as bigwig is 0-based
							real_end = real_start + len(hitseq)

							#count this match to the fisher_dict_bg
							fisher_dict_bg[key_for_bed_dict] += 1

							#get the values from bw file 
							bw_scores_control = np.mean(np.nan_to_num(np.array(list(bw_control.values(chromosom, real_start, real_end)))))
							bw_scores_overexpression = np.mean(np.nan_to_num(np.array(list(bw_overexpression.values(chromosom, real_start, real_end)))))

							control_dict = save_bw_score(one_motif, key_for_bed_dict, control_dict, bw_scores_control, hitseq, score)
							overexpression_dict = save_bw_score(one_motif, key_for_bed_dict, overexpression_dict, bw_scores_overexpression, hitseq, score)

							bw_difference = abs(bw_scores_overexpression - bw_scores_control)
							
							if not np.isnan(bw_difference) and bw_difference != 0.0: #do not need to check for nan
								differences.append(bw_difference)

							"""#--------------------------------------------------delete????------------------------------
								if not key_for_bed_dict in matches_dict: 
								# and not start in matches_dict[key_for_bed_dict]: #and float(matches_dict[key_for_bed_dict][start]['difference']) < bw_difference: #if there is no information saved about this match yet and if so, the saved difference is smaller that the one we have just found
									matches_dict[key_for_bed_dict] = matches_dict.get(key_for_bed_dict, {})
									matches_dict[key_for_bed_dict][start] = {'moods_score': score, 'hitseq': hitseq, 'difference': bw_difference}
							"""

							"""
							for i in range(len(hitseq)):
								#check for the score, if there is already some score saved on this position, save the greatest (better) score
								if scores_dict[key_for_bed_dict][int(start) + i] != 0.0:
									#print(scores_dict[key_for_bed_dict][int(start) + i])
									scores_dict[key_for_bed_dict][int(start) + i] = max(bw_scores[i], scores_dict[key_for_bed_dict][int(start) + i])
									#print(scores_dict[key_for_bed_dict][int(start) + i])
								else:
									scores_dict[key_for_bed_dict][int(start) + i] = bw_scores[i]

							#if we already have background, use it to make the dictionary for fisher exact test
							if background_dict:
								#print(scores_dict[key_for_bed_dict][int(start):int(start) + len(hitseq)])
								overexpression_mean = np.mean(np.array(scores_dict[key_for_bed_dict][int(start):int(start) + len(hitseq)]))
								if not np.isnan(overexpression_mean) and overexpression_mean > background_dict[key_for_bed_dict]:
										motifs_dict[key_for_bed_dict] += 1
							else:
								motifs_dict[key_for_bed_dict] += 1
							"""#-----------------------------------------------------------------------------------------------------

			#one doesnt need to close file that was opened like so, as python does it on itself. file.closed says True			
			return fisher_dict_bg, fisher_dict, control_dict, overexpression_dict, differences
		else:
			logger.info("The input for moods was not validated by the MOODS.parsers.pfm. Please check if it has the right format (note that the MOODS accepts only the old version of .pfm files, that is one without the header containing the name and id of the motif)")
			sys.exit()

#check if the bw score is already saved, if so check if it is bigger than the new one <--------------------------- do i need to save motif name in this matrix???
def save_bw_score(motif, key_for_bed_dict, matches_dict, bw_score, hitseq, score):

	if np.isnan(bw_score): bw_score = 0.0

	#bw_score = float(bw_score) * float(score) #apply the moods score

	if key_for_bed_dict in matches_dict:

		""" #---------------------------delete or leave as a feature??-----------
		if float(matches_dict[key_for_bed_dict]['bw_score']) < bw_score: #save the match with the bigger score
			#print(matches_dict[key_for_bed_dict]['bw_score'])
			matches_dict[key_for_bed_dict] = {'bw_score': bw_score, 'hitseq': hitseq, 'moods_score': score}
			#print(matches_dict[key_for_bed_dict]['bw_score'])
		"""
		#save the mean of the boths scores
		#matches_dict[key_for_bed_dict] = np.mean([matches_dict[key_for_bed_dict], bw_score])
		#----------------------------------------------------------------------

		#save the biggest of scores
		if matches_dict[key_for_bed_dict] < bw_score:
			matches_dict[key_for_bed_dict] = bw_score
	else:
		matches_dict[key_for_bed_dict] = bw_score

	return matches_dict

#help function for the moods call, convert pfm to log odds
def pfm_to_log_odds(filename, bg, pseudocount):
	if pfm(filename):
		mat = MOODS.parsers.pfm_to_log_odds(filename, bg, pseudocount)
		if len(mat) != 4: #if something went wrong, the empty list will be returned
			return False, mat
		else:
			return True, mat
	else:
		logger.info('please make sure the motif file has a .pfm format needed for moods')
		sys.exit()

#help function for the moods call, check if the file is in a pfm format using moods
def pfm(filename):
	mat = MOODS.parsers.pfm(filename)
	if len(mat) != 4:
		return False
	else:
		return True

def remove_file(file):
	if os.path.isfile(file):
		os.remove(file)

def clean_directory(cleans, output_directory, motif, tool_output_file): #<-----------------------check if it is working

	moods_output_unsorted = os.path.join(tool_output_file.replace("moods_output", "moods_output_unsorted"))
	
	for clean in cleans:
		if clean == 'all':
			remove_file(motif)
			remove_file(tool_output_file)
			remove_file(moods_output_unsorted)
		elif clean == 'cut_motifs':
			remove_file(motif)
		elif clean == 'moods_output':
			if os.path.basename(tool_output_file).startswith("moods"):
				remove_file(tool_output_file)

def tool_make_output(motif, genome, output_directory, cleans, p_value, bed_dictionary, moods_bg, bw_overexpression, bw_control, global_mean, global_std):

	standard_moods_bg = 0.25, 0.25, 0.25, 0.25

	fisher_dict = defaultdict(int)	
	fisher_dict_bg = defaultdict(int)	

	control_dict = {}
	overexpression_dict = {}

	differences = []
	#--------------------do we need the  differences here?
	fisher_dict_bg, fisher_dict, control_dict, overexpression_dict, differences = call_moods(motif, genome, output_directory, p_value, standard_moods_bg, fisher_dict_bg, fisher_dict, bw_control, bw_overexpression, control_dict, overexpression_dict, differences) 
	
	#make arrays of the dictionaries
	control_array = []
	overexpression_array = []

	for key in control_dict:
		control_array.append(control_dict[key])
		overexpression_array.append(overexpression_dict[key])

	motif_name = motif.replace(output_directory, '')
	motif_name = motif_name.replace('/', '')

	#make the wilcoxon signed-rank test
	my_wilcoxon_pvalue, direction, differences, differences_normalized, motif_std = my_wilcoxon(bw_control, bw_overexpression, control_array, overexpression_array, global_mean, global_std, motif_name, correction = False)

	return my_wilcoxon_pvalue, fisher_dict, direction, differences, differences_normalized, motif_std

def make_normal_distribution_plot(input_array, motif_name):
	figure_name = motif_name.replace("pfm", "png")

	output_directory = "./plots3"

	mu = np.mean(input_array)
	std = np.std(input_array, ddof = 1)

	#axes = plt.gca() #gca = get the current axes
	#axes.set_ylim([0, 80])

	fig, ax = plt.subplots()

	plt.hist(input_array, bins = len(input_array), color = 'g', normed=True)
	xmin, xmax = plt.xlim()
	x = np.linspace(xmin, xmax, 100)
	p = stats.norm.pdf(x, mu, std)
	plt.plot(x, p, 'r', linewidth = 3)
	motif_name = motif_name.replace(".pfm", "")
	title = motif_name + " fit results: mu = %.4f, std = %.4f" % (mu, std)

	plt.title(title)
	plt.grid()

	plt.ylim(0, 50)

	fig.savefig(os.path.join(output_directory, figure_name))
	#plt.show()

#modified from scipy.wilcoxon
def my_wilcoxon(bw_control, bw_overexpression, x, y, global_mean, global_std, motif_name, correction = False):

	x, y = map(np.asarray, (x, y)) #apply np.asarray to both input arrays

	if len(x) != len(y):
		raise ValueError("The length of both arrays in Wilcoxon test should be the same. Aborting")

	#normalize both arrays considering the global mu and std for overexpression and control
	
	#make_normal_distribution_plot(x, "control_no_normalization")
	#make_normal_distribution_plot(y, "overexpression_no_normalization")

	#x = (x - mu_control) / std_control
	#y = (y - mu_overexpression) / std_overexpression

	#make_normal_distribution_plot(x, "control_with_normalization")
	#make_normal_distribution_plot(y, "overexpression_with_normalization")

	d = x - y #find the difference

	#make_normal_distribution_plot(d, "differences_with_normalization")

	#keep all non-zero differences
	d = np.compress(np.not_equal(d, 0), d) #in scipy axis = -1, in my case it does not matter, as i have flattened array

	#make_normal_distribution_plot(d, motif_name)

	#correct the differences according to the normal distribution
	#mu = np.mean(d)
	#std = np.std(d, ddof = 1) #ddof = 1 calculates corrected sample sd which is sqrt(N/(N-1)) times the population sd where N is the number of points, interpretes the data as samples, estimates true variance
	#https://stackoverflow.com/questions/15389768/standard-deviation-of-a-list

	#correct the differences according to the global mean and std
	d_normalized = (d - global_mean) / global_std

	#motif_name = motif_name.replace(".pfm", "_normalized.pfm")
	#make_normal_distribution_plot(d, motif_name)

	count = len(d_normalized)

	if count < 10:
		logger.info("The sampe size is too small for normal approximation")

	r = stats.rankdata(abs(d_normalized)) #assign ranks to data, dealing with ties appropriately

	r_plus = np.sum((d_normalized > 0) * r, axis = 0)
	r_minus = np.sum((d_normalized < 0) * r, axis = 0)

	"""
	#direction towards first array, plus
	first_array = max(d)
	#direction towards second array, minus
	second_array = abs(min(d))

	if first_array > second_array:
		#x was bigger
		direction = "first_array"
	else:
		#y was bigger
		direction = "second_array"

	#print(r_plus, r_minus)
	"""

	T = min(r_plus, r_minus)
	mn = count * (count + 1.) * 0.25
	se = count * (count + 1.) * (2. * count + 1.)

	replist, repnum = stats.find_repeats(r)

	if repnum.size != 0:
		#correction for repeated elements
		se -= 0.5 * (repnum * (repnum * repnum -1)).sum()

	se = np.sqrt(se / 24)
	correction = 0.5 * int(bool(correction)) * np.sign(T - mn)
	z = (T - mn - correction) / se
	#scale down
	prob = 2. * stats.norm.sf(abs(z), scale = 1) #before scale = 2.5 <-----------------------------------------------------

	motif_std = np.std(d_normalized, ddof = 1)
	motif_mu = np.mean(d_normalized)

	#direction = "control" #first_array
	direction = get_name_from_path(bw_control)
	if motif_mu < 0:
		#direction = "overexpression" #second_array
		direction = get_name_from_path(bw_overexpression)

	return prob, direction, d, d_normalized, motif_std

def make_name_from_path(full_path, output_directory, ending):
	return os.path.join(output_directory, get_name_from_path(full_path) + ending)
	
def get_name_from_path(full_path):
	return os.path.splitext(os.path.basename(full_path))[0]

def multiprocess(motifs, genome, output_directory, cleans, p_value, bed_dictionary, cpu_count, moods_bg, condition1, condition2, output_number):

	pool = multiprocessing.Pool(cpu_count) #by default cpu_count is 2

	#length = len(motifs) #the number of the motifs to find the percentage of the job that was done
	#step = 100/length #the percentage that should be added after the job with each single motif is done
	#tasks = [] #here the jobs done by processes are saved

	fisher_dict = dict()
	motifs_array = []
	p_values_array = []
	directions_array = []

	global_mean, global_std = compute_differences(bed_dictionary, condition1, condition2)

	all_differences = []
	all_differences_normalized = []
	motifs_p_values = open("test_motifs_p_values11.txt", 'w')

	all_stds = {}

	#first find motifs using standard background
	for motif in motifs:
		wilcoxon_p_value, fisher_dict, direction, differences, differences_normalized, motif_std = pool.apply_async(tool_make_output, args = (motif, genome, output_directory, cleans, p_value, bed_dictionary, moods_bg, condition1, condition2, global_mean, global_std, )).get()

		motif_name = motif.replace(output_directory, "") #find an elegant way to do it
		motif_name = motif_name.replace("/", "")
		motif_name = motif_name.replace(".pfm", "")

		motifs_array.append(motif_name)
		p_values_array.append(wilcoxon_p_value)
		directions_array.append(direction)

		#-------------
		all_stds[motif_name] = motif_std
		all_differences.extend(differences)
		all_differences_normalized.extend(differences_normalized)
		#-------------
		print(motif_name, wilcoxon_p_value, direction)
		motifs_p_values.write('\t'.join([str(motif_name), str(wilcoxon_p_value), direction]) + '\n')

	motifs_p_values.close()

	logger.info('Applying Bonferroni correction')
	p_values_adjusted = multipletests(p_values_array, method = 'bonferroni')[1]

	dict_motifs_p_values = {}

	for i in range(len(motifs_array)):
		motif = motifs_array[i]
		dict_motifs_p_values[motif] = dict_motifs_p_values.get(motif, {})
		dict_motifs_p_values[motif] = {'p_value': p_values_array[i], 'adjusted_p_value': p_values_adjusted[i], 'direction': directions_array[i]}

	sorted_dict = sorted(dict_motifs_p_values.items(), key = lambda x : (x[1]['p_value']), reverse = False)

	logger.info('The most enriched motifs are:')
	logger.info("{:30s} | {:30s} | {:30s} | {:30s}".format('Motif', 'p_value', 'adjusted_p_value', 'direction'))
	logger.info('-' * 120)
	#for i in range(output_number):
	for i in range(len(motifs_array)):
		#logger.info(sorted_dict[i])
		logger.info("{:30s} | {:30s} | {:30s} | {:30s}".format(str(sorted_dict[i][0]), str(sorted_dict[i][1]['p_value']), str(sorted_dict[i][1]['adjusted_p_value']), str(sorted_dict[i][1]['direction'])))

	logger.info('\n')

	for motif in sorted_dict:
		if str(motif[0]).startswith("DUX") or str(motif[0]).startswith("Dux"):
			logger.info(motif)

	logger.info('\n')

	max_pvalues_array = max(p_values_array)

	logger.info("the greatest pvalue is " + str(max_pvalues_array))
	logger.info("the smallest pvalue is " + str(min(p_values_array)))
	logger.info("the number of all motifs is " + str(len(p_values_array)))
	logger.info("the number of significant motifs is " + str(sum(i <= 1e-5 for i in p_values_array)))

	
	mu = np.mean(all_differences)
	std = np.std(all_differences, ddof = 1) #ddof = 1 calculates corrected sample sd which is sqrt(N/(N-1)) times the population sd where N is the number of points, interpretes the data as samples, estimates true variance
	#https://stackoverflow.com/questions/15389768/standard-deviation-of-a-list
	plt.hist(all_differences, bins = 1500, color = 'g', normed=True)
	xmin, xmax = plt.xlim()
	x = np.linspace(xmin, xmax, 25000)
	p = stats.norm.pdf(x, mu, std)
	#p = stats.norm.pdf(x)
	plt.plot(x, p, 'r', linewidth = 1)
	title = "All differences before normalization: mu = %.4f, std = %.4f" % (mu, std)
	plt.title(title)
	plt.grid()
	plt.show()

	"""
	mu = np.mean(all_differences)
	std = np.std(all_differences, ddof = 1) #ddof = 1 calculates corrected sample sd which is sqrt(N/(N-1)) times the population sd where N is the number of points, interpretes the data as samples, estimates true variance
	#https://stackoverflow.com/questions/15389768/standard-deviation-of-a-list
	plt.hist(all_differences, bins = 1500, color = 'g', normed=True)
	xmin, xmax = plt.xlim()
	x = np.linspace(xmin, xmax, 25000)
	#p = stats.norm.pdf(x, mu, std)
	p = stats.norm.pdf(x)
	plt.plot(x, p, 'r', linewidth = 1)
	title = "All differences before normalization: mu = %.4f, std = %.4f" % (mu, std)
	plt.title(title)
	plt.grid()
	plt.show()
	"""


	mu_normalized = np.mean(all_differences_normalized)
	std_normalized = np.std(all_differences_normalized, ddof = 1)
	plt.hist(all_differences_normalized, bins = 1500, color = 'b', normed=True)
	xmin, xmax = plt.xlim()
	x_normalized = np.linspace(xmin, xmax, 25000)
	p_normalized = stats.norm.pdf(x_normalized, mu_normalized, std_normalized)
	#p = stats.norm.pdf(x)
	plt.plot(x_normalized, p_normalized, 'r', linewidth = 1)
	title_normalized = "All differences after normalization: mu = %.4f, std = %.4f" % (mu_normalized, std_normalized)
	plt.title(title_normalized)
	plt.grid()
	plt.show()

	mu_normalized = np.mean(all_differences_normalized)
	std_normalized = np.std(all_differences_normalized, ddof = 1)
	plt.hist(all_differences_normalized, bins = 1500, color = 'b', normed=True)
	xmin, xmax = plt.xlim()
	x_normalized = np.linspace(xmin, xmax, 25000)
	#p_normalized = stats.norm.pdf(x_normalized, mu_normalized, std_normalized)
	p = stats.norm.pdf(x)
	plt.plot(x_normalized, p_normalized, 'r', linewidth = 1)
	title_normalized = "All differences after normalization: mu = %.4f, std = %.4f" % (mu_normalized, std_normalized)
	plt.title(title_normalized)
	plt.grid()
	plt.show()

	"""
	if mu != 0.0 and std != 1.0:
		all_differences_corrected = (all_differences - mu ) / std

		mu_corrected = np.mean(all_differences_corrected)
		std_corrected = np.std(all_differences_corrected, ddof = 1)
		plt.hist(all_differences_corrected, bins = 1500, color = 'b', normed=True)
		xmin, xmax = plt.xlim()
		x = np.linspace(xmin, xmax, 25000)
		p = stats.norm.pdf(x, mu, std)
		#p = stats.norm.pdf(x)
		plt.plot(x, p, 'r', linewidth = 1)
		title = "corrected results: mu = %.4f, std = %.4f" % (mu_corrected, std_corrected)
		plt.title(title)
		plt.grid()
		plt.show()
	"""



	"""
	tasks.append(pool.apply_async(tool_make_output, args = (usage, motif, genome, output_directory, cleans, p_value, bed_dictionary, fimo_data, resolve_overlaps, moods_bg, )))
	tasks_done = sum([task.ready() for task in tasks]) #the number of the processes that ended their job
	#check the number of the processes that are ready till the number of them reaches the number of processes started in the pool
	while tasks_done < len(tasks):
		#if the number of ready processes has changed, save the new number and print the percentage of the job done
		if sum([task.ready() for task in tasks]) != tasks_done:
			tasks_done = sum([task.ready() for task in tasks])
			sys.stdout.write("%-100s| %d%% \r" % (''*tasks_done, step*tasks_done))
			sys.stdout.flush()
			sys.stdout.write("\n")
		#update the number of ready processes each 0.05 seconds
		time.sleep(0.05)
	"""

	pool.close()
	pool.join() #make sure all the processes are done and exit
		
	#the processes should not delete the merged genome file
	#so make sure if this file is needed, otherwise delete it
	for clean in cleans:
		if clean == 'all' or clean == 'merge_output':
			for filename in os.listdir(output_directory):
				if filename == "output_merge.fa":
					remove_file(genome)

		if clean != 'nothing':
			logger.info('the directory ' + output_directory + ' was cleaned, only the required files were left')

def compute_differences(bed_dictionary, condition1, condition2):
	logger.info("the mean and standard deviation for the differences in peaks will be count now")

	bw_cond1 = pyBigWig.open(condition1)
	bw_cond2 = pyBigWig.open(condition2)

	global_differences = {} #dict
	differences_array = [] #to compute the mean at the end
	cond1_array = []
	cond2_array = []

	for header in bed_dictionary:
		header_splitted = re.split(r':', header)
		chromosom = header_splitted[0]
		positions = re.split(r'-', header_splitted[-1])

		#compute the background difference for this peak
		bw_global_score_cond1 = np.mean(np.nan_to_num(np.array(list(bw_cond1.values(chromosom, int(positions[0]), int(positions[1]))))))
		bw_global_score_cond2 = np.mean(np.nan_to_num(np.array(list(bw_cond2.values(chromosom, int(positions[0]), int(positions[1]))))))
		bw_global_difference = bw_global_score_cond2 - bw_global_score_cond1

		global_differences[header] = bw_global_difference

		cond1_array.append(bw_global_score_cond1)
		cond2_array.append(bw_global_score_cond2)

		differences_array.append(bw_global_difference)

	bw_cond1.close()
	bw_cond2.close()

	mu = np.mean(differences_array)
	std = np.std(differences_array, ddof = 1)

	mu_cond1 = np.mean(cond1_array)
	#std_cond1 = np.std(cond1_array, ddof = 1)

	mu_cond2 = np.mean(cond2_array)
	#std_cond2 = np.std(cond2_array, ddof = 1)

	#make_normal_distribution_plot(cond1_array, "global_cond1_peaks")
	#make_normal_distribution_plot(cond2_array, "global_cond2_peaks")

	cond1_name = get_name_from_path(condition1)
	cond2_name = get_name_from_path(condition2)

	#if (mu_cond1 != mu_cond2):
	#	logger.info("need to normalize! The mean of " + str(cond1_name) + " is " + str(mu_cond1) + ", the mean of " + str(cond2_name) " is " + str(mu_cond2))
	return mu, std

def is_fasta(check_fasta):
	if not os.path.isfile(check_fasta):
		logger.info('there is no file with genome, the exit is forced')
		sys.exit()
	else:
		# modified code from https://stackoverflow.com/questions/44293407/how-can-i-check-whether-a-given-file-is-fasta
	    with open(check_fasta, "r") as handle:
	        fasta = SeqIO.parse(handle, "fasta")
	        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file

def check_existing_input_files(args):

	"""
	cond1 = pyBigWig.open(args.condition1)
	cond2 = pyBigWig.open(args.condition2)

	elif not cond1.isBigWig() or not cond2.isBigWig():
		logger.info('please provide bigWig files for the both conditions you want to compare')
		sys.exit()

	cond1.close()
	cond2.close()
	"""
	if not is_fasta(args.genome):
		logger.info('please make sure the input genome file has a fasta format')
		sys.exit()
	if not os.path.isfile(args.condition1) or not os.path.isfile(args.condition2):
		logger.info('please make sure the both files with conditions to compare exist')
		sys.exit()
	if not args.condition1.endswith('.bw') or not args.condition2.endswith('.bw'):
		logger.info('please check if the both conditions files are in bigWig format')
		sys.exit()
	#check if the file with motifs exists
	if not os.path.isfile(args.motifs):
		logger.info('there is no file with motifs, the exit is forced')
		sys.exit()
	#if the bedfile was given as input, check if this file exists
	if args.bed_file:
		if not os.path.isfile(args.bed_file):
			logger.info('there is no such bed file ' + args.bed_file + ', the exit is forced')
			sys.exit()

def make_bed_dictionary(bed_file):

	bed_dictionary = {}

	with open(bed_file) as read_bed_file:
		for bed_line in read_bed_file:
			bed_line_array = re.split(r'\t', bed_line.rstrip('\n'))
			if bed_line_array[1].isdigit() and bed_line_array[2].isdigit() and int(bed_line_array[1]) <= int(bed_line_array[2]): #in the real bedfile the second column is a start position, and the third column is an end position, so we are checking if these are integers and if the start position is smaller than the end one
				key = bed_line_array[0] + ":" + bed_line_array[1] + "-" + bed_line_array[2]
				value = []
				for i in range(3, len(bed_line_array)):
					value.append(bed_line_array[i]) 

				bed_dictionary[key] = value
			else: #this is not a bed file, force exit
				logger.info('please make sure the input bed file has a right format, the problem occured on the line ' + bed_line)
				sys.exit()

	read_bed_file.close()

	return bed_dictionary

def print_small_logo():
	winnie_logo = """\
	 __      ___ _  _ _  _ _     
	 \ \    / (_) \| | \| (_)___ 
	  \ \/\/ /| | .` | .` | / -_)
	   \_/\_/ |_|_|\_|_|\_|_\___|
	                             
	"""
	logger.info(winnie_logo)

def print_big_logo():
	winnie_big_logo = """\


     
           ▄▒▀▀▀▒▄             ▄▄▀▀▒▄
        ▄▄▒▀       ▀▒▄▄     ▄▄▒▀       ▀▒▄
      ▐▒               ▐▒ ▒▀               ▒
      ▐▌               ▐▒ ▒                ▐▌
      ▐▌               ▐▒ ▒                ▐▌
      ▐▌               ▐▒ ▒          ▒█▌   ▐▌                      ▒█▌
      ▐▌        ▐█▄    ▐▒ ▒      █▀   ▀▀   ▐█▄   ▐█     █▄    █     ▀▀
        ▀▒▄▄      █▄▄▒▀ ▐█  ▀▒▄▄█      ▄▄▒▀▐██▄  ▐█    ▐██▌   █            ▄█▀▌
            ▀▒▄▄▄▒▀█▌▄▒▐█▀█▒▄▄██▀▄▄▄▄▐█    ▐█ ▀▌ ▐█    ▐█ ▀█  █     █     ▓█▄▄▐█▄
                 ▄▄▒▀█ █▌ █▌ █▀░▄▄   ▐█    ▐█  ▀█▄█    ▐█  ▀█▄█     █     ▀█▀▀▀▀
                ▐▌   ▀██   ██    ▐▌  ▐█    ▐█    ██    ▐█    ██     █       █▄
                ▐                ▐▌
                ▐                ▐▌
                ▐                ▐▌
                ▐▒▄             ▄▒▌
                   ▀▒▄▄     ▄▄▒▀
                       ▀▒▄▒▀
     
    

	"""
	logger.info(winnie_big_logo)

def main():

	start = time.time()
	print_big_logo()

	args = parse_args()

	#if user do not want to see the information about the status of jobs, remove the handler, that writes to the terminal
	if args.silent:
		logger.removeHandler(ch)

	#motifs = "../../PaperInPrep/TOBIAS/buenrostro_analysis/data/buenrostro_motifs.meme"
	#genome = "../analysis_my_tool/chipseq/hg19.fasta"
	#bed_file = "../analysis_my_tool/chipseq/hg19_peaks.bed"

	#bw_overexpression = "overexpression_footprints.bw"
	#bw_control = "control_footprints.bw"

	#bed_file = "../TOBIAS/tobias/TFBScan/peaks.bed"

	"""
	bw_overexpression = "../../mette.bentsen/to_anastasiia/duxbl_footprints.bw"
	bw_control = "../../mette.bentsen/to_anastasiia/gfp_footprints.bw"

	#to check with ame
	motifs = "../my_ame/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.meme"

	genome = "../my_ame/mm10.fa"
	#genome = "../my_ame/overexpression_peaks.fa"

	bed_file = "../../mette.bentsen/to_anastasiia/all_merged.bed"
	#bed_file = "../my_ame/overexpression_peaks.bed"

	"""

	logger.info("WiNNie was called using these parameters:")
	logger.info(vars(args))
	check_existing_input_files(args)

	#check if there is an existing directory that user gave as input, otherwise create this directory from the path provided from the user
	check_directory(args.output_directory)

	splitted_motifs = split_motifs(args.motifs, args.output_directory)
	
	#check if there is a .bed file to merge the genome file with. If so, merge them
	if args.bed_file:
		bed_dictionary = make_bed_dictionary(args.bed_file)
		args.genome = merge(args.genome, args.bed_file, args.output_directory)
	else:
		bed_dictionary = {}

	multiprocess(splitted_motifs, args.genome, args.output_directory, args.cleans, args.p_value, bed_dictionary, args.cores, args.moods_bg, args.condition1, args.condition2, args.output_number)

	logger.info("WiNNie needed %s seconds to make the output" % (time.time() - start))
	
	for handler in logger.handlers:
		handler.close()
		logger.removeFilter(handler)
	
if __name__ == "__main__":
	main()