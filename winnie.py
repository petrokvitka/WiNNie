
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

logger = logging.getLogger('winnie')
logger.setLevel(logging.INFO)

formatter = logging.Formatter("%(asctime)s : %(message)s", "%Y-%m-%d %H:%M")

#catch all the information about input and output files as well as information on the used tool (fimo or moods)
def parse_args():
	
	parser = argparse.ArgumentParser(prog = 'winnie', description = textwrap.dedent('''                           

		This script takes a list of motifs loaded from jaspar.genereg.net as a combined text file in .MEME or .PFM format, two bigWig-files containing scores for two different conditions, a genome file in FASTA format and a .bed file with regions of interest as input. The output is a file containing enriched motifs sorted by adjusted p value. Please note, if you want to have all intermediate output files, enter --clean nothing.
		'''), epilog='That is what you need to make this script work for you. Enjoy it')
	
	required_arguments = parser.add_argument_group('required arguments')
	required_arguments.add_argument('-m', '--motifs', help='file in .MEME or .PFM format with mofits loaded from jaspar.genereg.net', required=True)
	required_arguments.add_argument('-g', '--genome', help='a whole genome file or regions of interest in FASTA format to be scanned with motifs', required=True)
	required_arguments.add_argument('-b', '--bed_file',  nargs='?', help='a .bed file to be merged with the whole genome file to find regions of interest')
	required_arguments.add_argument('-c1', '--condition1', help='a bigWig-file with scores for the first condition', required=True)
	required_arguments.add_argument('-c2', '--condition2', help='a bigWig-file with scores for the second condition', required=True)

	#all other arguments are optional
	parser.add_argument('-o', '--output_directory',  default='output', const='output', nargs='?', help='output directory, default ./output/')
	parser.add_argument('--clean', nargs='*', choices=['nothing', 'all', 'cut_motifs', 'merge_output'], dest='cleans', help='choose the files you want to delete from the output directory, the default is deleting all the temporary files from the directory', default=['all'])
	parser.add_argument('--cores', type=int, help='number of cores allowed to use by this tool, by default the tool uses 2 cores', default=2)
	parser.add_argument('-p', '--p_value', type=float, help='enter the p value, the default p value is 1e-4', default=0.0001)
	parser.add_argument('--silent', action='store_true', help='while working with data write the information only into ./winnie_log.txt')
	parser.add_argument('--moods_bg', nargs='+', type=float, help='set the bg for moods, by default moods uses the bg is 0.25 0.25 0.25 0.25', default=[0.25, 0.25, 0.25, 0.25])
	parser.add_argument('--output_number', type=int, help='number of enriched motifs that should be written to terminal or log-file, by default 40', default=40)
	parser.add_argument('--score', choices=['mean', 'greater'], help='apply mean of the scores or the greatest score from the condition file', default='greater')
	args = parser.parse_args()

	return args

def check_directory(directory):
	if not os.path.exists(directory):
		os.makedirs(directory)
		#logger.info('a new directory ' + directory + ' was created')
		print('a new directory ' + directory + ' was created')
					
#merge the whole genome with the regions mentioned in .bed file
def merge(genome, bed_file, output_directory):	
	output_merge = os.path.join(output_directory, "output_merge.fa")
	logger.info('the merging of files ' + genome + ' and ' + bed_file + ' will end soon, the result file is ' + output_merge)
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
		logger.info("the input file has an unexpected format, i will try to convert it to .pfm format")
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

def call_moods(one_motif, genome, output_directory, p_value, moods_bg, condition2, condition1, control_dict, overexpression_dict, differences, which_score):

	#check if this is a bigwig file
	bw_condition2 = pyBigWig.open(condition2)
	bw_condition1 = pyBigWig.open(condition1)
	if not bw_condition2.isBigWig() or not bw_condition1.isBigWig(): 
		logger.info("please provide the bigwig file!")
		sys.exit()
	else:
		# prepare everything for moods
		# setting standard parameters for moods
		# this code was token and modified from gitHub MOODS page
		pseudocount = 0.0001

		bg = tuple(moods_bg)

		matrix_names = [os.path.basename(one_motif)]

		matrices = []
		matrices_rc = []

		valid, matrix = pfm_to_log_odds(one_motif, bg, pseudocount)

		key_for_bed_dict = ''

		if valid:

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

					chromosom = header_splitted[0]
					positions = re.split(r'-', header_splitted[-1])

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

							#find the real start and end positions on the chromosom		
							real_start = int(positions[0]) + int(start) #start of the peak + start of the motif within the peak, do not add 1, as bigwig is 0-based
							real_end = real_start + len(hitseq)

							#get the values from bw file 
							bw_scores_control = np.mean(np.nan_to_num(np.array(list(bw_condition2.values(chromosom, real_start, real_end)))))
							bw_scores_overexpression = np.mean(np.nan_to_num(np.array(list(bw_condition1.values(chromosom, real_start, real_end)))))

							control_dict = save_bw_score(key_for_bed_dict, control_dict, bw_scores_control, float(score), which_score)
							overexpression_dict = save_bw_score(key_for_bed_dict, overexpression_dict, bw_scores_overexpression, float(score), which_score)

							bw_difference = abs(bw_scores_overexpression - bw_scores_control)
							
							if not np.isnan(bw_difference) and bw_difference != 0.0: #do not need to check for nan
								differences.append(bw_difference)

			#one doesnt need to close file that was opened like so, as python does it on itself. file.closed says True			
			return control_dict, overexpression_dict, differences
		else:
			logger.info("The input for moods was not validated by the MOODS.parsers.pfm. Please check if it has the right format (note that the MOODS accepts only the old version of .pfm files, that is one without the header containing the name and id of the motif)")
			sys.exit()

#check if the bw score is already saved, if so check if it is bigger than the new one
def save_bw_score(key_for_bed_dict, matches_dict, bw_score, moods_score, which_score):

	if np.isnan(bw_score): bw_score = 0.0

	#bw_score = moods_score * bw_score #apply moods score as well 

	if key_for_bed_dict in matches_dict:

		if which_score == "mean":
			#save the mean of both scores
			matches_dict[key_for_bed_dict] = np.mean([matches_dict[key_for_bed_dict], bw_score])
		
		elif which_score == "greater":
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

def clean_directory(cleans, output_directory, motif):
	
	for clean in cleans:
		if clean == 'all' or clean == 'cut_motifs':
			remove_file(motif)
			#the output_merge.fa will be deleted after processing of the multiprocessing

def tool_make_output(motif, genome, output_directory, cleans, p_value, bed_dictionary, moods_bg, condition1, condition2, global_mean, global_std, which_score):

	standard_moods_bg = 0.25, 0.25, 0.25, 0.25

	control_dict = {}
	overexpression_dict = {}

	differences = []
	control_dict, overexpression_dict, differences = call_moods(motif, genome, output_directory, p_value, standard_moods_bg, condition2, condition1, control_dict, overexpression_dict, differences, which_score) 
	
	#make arrays of the dictionaries
	control_array = []
	overexpression_array = []

	for key in control_dict:
		control_array.append(control_dict[key])
		overexpression_array.append(overexpression_dict[key])

	motif_name = motif.replace(output_directory, '')
	motif_name = motif_name.replace('/', '')

	#make the wilcoxon signed-rank test
	my_wilcoxon_pvalue, direction, differences, differences_normalized, motif_std = my_wilcoxon(condition2, condition1, control_array, overexpression_array, global_mean, global_std, motif_name, correction = False)

	clean_directory(cleans, output_directory, motif)

	return my_wilcoxon_pvalue, direction, differences, differences_normalized, motif_std

def make_normal_distribution_plot(input_array, figure_name, output_directory, figure_color):

	mu = np.mean(input_array)
	std = np.std(input_array, ddof = 1) #ddof = 1 calculates corrected sample sd which is sqrt(N/(N-1)) times the population sd where N is the number of points, interpretes the data as samples, estimates true variance
	#https://stackoverflow.com/questions/15389768/standard-deviation-of-a-list

	fig, ax = plt.subplots()

	plt.hist(input_array, bins = 3000, color = figure_color, normed=True, label = 'differences')
	xmin, xmax = plt.xlim()
	x = np.linspace(xmin, xmax, 3000)
	p = stats.norm.pdf(x, mu, std)
	plt.plot(x, p, 'r', linewidth = 1, label = 'normal distribution')
	title = figure_name + ": mu = %.4f, std = %.4f" % (mu, std)

	plt.title(title)
	plt.grid()

	plt.legend()

	x_lim = np.percentile(input_array, 95)

	plt.xlim(-x_lim, x_lim)

	fig.savefig(os.path.join(output_directory, figure_name))

#modified from scipy.wilcoxon
def my_wilcoxon(condition2, condition1, x, y, global_mean, global_std, motif_name, correction = False):

	x, y = map(np.asarray, (x, y)) #apply np.asarray for both input arrays

	if len(x) != len(y):
		raise ValueError("The length of both arrays in Wilcoxon test should be the same. Aborting")

	d = x - y #find the difference

	#keep all non-zero differences
	d = np.compress(np.not_equal(d, 0), d) #in scipy axis = -1, in my case it does not matter, as i have a flattened array

	#correct the differences according to the global mean and std
	d_normalized = (d - global_mean) / global_std

	count = len(d_normalized)

	if count < 10:
		logger.info("The sampe size is too small for normal approximation")

	r = stats.rankdata(abs(d_normalized)) #assign ranks to data, dealing with ties appropriately

	r_plus = np.sum((d_normalized > 0) * r, axis = 0)
	r_minus = np.sum((d_normalized < 0) * r, axis = 0)

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
	prob = 2. * stats.norm.sf(abs(z), scale = 1) #do not scale

	motif_std = np.std(d_normalized, ddof = 1)
	motif_mu = np.mean(d_normalized)

	direction = get_name_from_path(condition2)
	if motif_mu < 0:
		direction = get_name_from_path(condition1)

	return prob, direction, d, d_normalized, motif_std

def make_name_from_path(full_path, output_directory, ending):
	return os.path.join(output_directory, get_name_from_path(full_path) + ending)
	
def get_name_from_path(full_path):
	return os.path.splitext(os.path.basename(full_path))[0]

def multiprocess(motifs, genome, output_directory, cleans, p_value, bed_dictionary, cpu_count, moods_bg, condition1, condition2, output_number, global_mean, global_std, which_score):

	pool = multiprocessing.Pool(cpu_count) #by default cpu_count is 2

	motifs_array = []
	p_values_array = []
	directions_array = []

	all_differences = []
	all_differences_normalized = []

	all_stds = {}

	count = 1 #a count for printing the number of motif the processes are working with

	for motif in motifs:
		motif_name = get_name_from_path(motif)
		logger.info("i am working with motif " + str(motif_name) + "\t" + str(count) + "/" + str(len(motifs)))
		wilcoxon_p_value, direction, differences, differences_normalized, motif_std = pool.apply_async(tool_make_output, args = (motif, genome, output_directory, cleans, p_value, bed_dictionary, moods_bg, condition1, condition2, global_mean, global_std, which_score, )).get()

		motifs_array.append(motif_name)
		p_values_array.append(wilcoxon_p_value)
		directions_array.append(direction)

		#-------------
		all_stds[motif_name] = motif_std
		all_differences.extend(differences)
		all_differences_normalized.extend(differences_normalized)
		#-------------

		count = count + 1

	logger.info('Applying Bonferroni correction')
	p_values_adjusted = multipletests(p_values_array, method = 'bonferroni')[1]

	dict_motifs_p_values = {}

	for i in range(len(motifs_array)):
		motif = motifs_array[i]
		dict_motifs_p_values[motif] = dict_motifs_p_values.get(motif, {})
		dict_motifs_p_values[motif] = {'p_value': p_values_array[i], 'adjusted_p_value': p_values_adjusted[i], 'direction': directions_array[i]}

	sorted_dict = sorted(dict_motifs_p_values.items(), key = lambda x : (x[1]['adjusted_p_value']), reverse = False) #sort after adjusted p value

	output_file = make_name_from_path("winnie_output", output_directory, ".txt")
	opened_output_file = open(output_file, 'w')

	logger.info('i will write the motifs to the file ' + output_file)
	
	opened_output_file.write('\t'.join(['Motif', 'p_value', 'adjusted_p_value', 'direction']) + '\n')
	for i in range(len(motifs_array)):
		opened_output_file.write('\t'.join([str(sorted_dict[i][0]), str(sorted_dict[i][1]['p_value']), str(sorted_dict[i][1]['adjusted_p_value']), str(sorted_dict[i][1]['direction'])]) + '\n')

	opened_output_file.close()

	logger.info('The ' + str(output_number) + ' of the most enriched motifs are:')
	logger.info("{:30s} | {:30s} | {:30s} | {:30s}".format('Motif', 'p_value', 'adjusted_p_value', 'direction'))
	logger.info('-' * 120)
	for i in range(output_number):
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

	make_normal_distribution_plot(all_differences, "before_normalisation", output_directory, 'blue')
	make_normal_distribution_plot(all_differences_normalized, "after_normalisation", output_directory, 'green')

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

	mu_cond2 = np.mean(cond2_array)

	cond1_name = get_name_from_path(condition1)
	cond2_name = get_name_from_path(condition2)

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

	if not is_fasta(args.genome):
		#logger.info('please make sure the input genome file has a fasta format')
		print('please make sure the input genome file has a fasta format')
		sys.exit()
	if not os.path.isfile(args.condition1) or not os.path.isfile(args.condition2):
		#logger.info('please make sure the both files with conditions to compare exist')
		print('please make sure the both files with conditions to compare exist')
		sys.exit()
	if not args.condition1.endswith('.bw') or not args.condition2.endswith('.bw'):
		#logger.info('please check if the both conditions files are in bigWig format')
		print('please check if the both conditions files are in bigWig format')
		sys.exit()
	#check if the file with motifs exists
	if not os.path.isfile(args.motifs):
		#logger.info('there is no file with motifs, the exit is forced')
		print('there is no file with motifs, the exit is forced')
		sys.exit()
	#check if the bed file exists
	if not os.path.isfile(args.bed_file):
		#logger.info('there is no such bed file ' + args.bed_file + ', the exit is forced')
		print('there is no such bed file ' + args.bed_file + ', the exit is forced')
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

	args = parse_args()

	check_existing_input_files(args)
	#check if there is an existing directory that user gave as input, otherwise create this directory from the path provided from the user
	check_directory(args.output_directory)

	fh = logging.FileHandler(os.path.join(args.output_directory, "winnie_log.txt"))
	fh.setLevel(logging.INFO)
	fh.setFormatter(formatter)
	logger.addHandler(fh)

	ch = logging.StreamHandler()
	ch.setLevel(logging.INFO)
	ch.setFormatter(formatter)
	logger.addHandler(ch)

	#if user do not want to see the information about the status of jobs, remove the handler, that writes to the terminal
	if args.silent:
		logger.removeHandler(ch)

	print_big_logo()
	logger.info("WiNNie was called using these parameters:")
	logger.info(vars(args))

	splitted_motifs = split_motifs(args.motifs, args.output_directory)
	
	bed_dictionary = make_bed_dictionary(args.bed_file)
	args.genome = merge(args.genome, args.bed_file, args.output_directory)

	global_mean, global_std = compute_differences(bed_dictionary, args.condition1, args.condition2)

	multiprocess(splitted_motifs, args.genome, args.output_directory, args.cleans, args.p_value, bed_dictionary, args.cores, args.moods_bg, args.condition1, args.condition2, args.output_number, global_mean, global_std, args.score)

	logger.info("WiNNie needed %s seconds to generate the output" % (time.time() - start))
	
	for handler in logger.handlers:
		handler.close()
		logger.removeFilter(handler)
	
if __name__ == "__main__":
	main()