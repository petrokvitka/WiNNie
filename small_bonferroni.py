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

def print_small_logo():
	winnie_logo = """\
	 __      ___ _  _ _  _ _     
	 \ \    / (_) \| | \| (_)___ 
	  \ \/\/ /| | .` | .` | / -_)
	   \_/\_/ |_|_|\_|_|\_|_\___|
	                             
	"""
	print(winnie_logo)

def main():

	read_file = "test_motifs_p_values9.txt"
	write_file = "winnie_output.txt"
	motifs_array = []
	p_values_array = []
	directions = []
	counts = []

	with open(read_file) as r_file:
		for r_line in r_file:
			r_line_array = re.split(r'\s', r_line.rstrip('\n'))
			#if float(r_line_array[1]) > 9.7e-90:
			motifs_array.append(r_line_array[0])
			p_values_array.append(float(r_line_array[1]))
			directions.append(r_line_array[2])
			counts.append(r_line_array[3])
			
	r_file.close()

	#apply bonferroni correction
	p_adjusted = multipletests(p_values_array, method = 'bonferroni')[1]

	dict_motifs_pvalues = {}

	for i in range(len(motifs_array)):
		motif = motifs_array[i]
		dict_motifs_pvalues[motif] = dict_motifs_pvalues.get(motif, {})
		#test_p_value = p_adjusted[i] * int(counts[i])
		dict_motifs_pvalues[motif] = {'p_value': p_values_array[i], 'adjusted_p_value': p_adjusted[i], 'direction': directions[i]} #, 'count': counts[i], 'test_p_value': test_p_value}
		#if p_values_array[i] == 0 or p_adjusted[i] == 0:
			#print(motif, "fisher_p_value", p_values_array[i], "adjusted_p_value", p_adjusted[i])

	#now sort the data after adjusted p_value
	sorted_dict = sorted(dict_motifs_pvalues.items(), key = lambda x : (x[1]['p_value']), reverse = False)

	for i in range(len(motifs_array)): 
		print(sorted_dict[i])

	print()

	for motif in sorted_dict:
		if str(motif[0]).startswith("DUX") or str(motif[0]).startswith("Dux"):
			print(motif)

	print()

	#find the number of significant motifs, so the ones, where the pvalue is smaller than treshold 1e-5
	print("the greatest pvalue is ", max(p_values_array))
	print("the smallest pvalue is ", min(p_values_array))
	print("the number of all motifs is ", len(p_values_array))
	#print("the number of significant motifs is ", sum(i < 1e-5 for i in p_values_array))
	#print("adjusted p values")
	print("the number of significant motifs is ", sum(i <= 1e-5 for i in p_values_array))

	
if __name__ == "__main__":
	main()