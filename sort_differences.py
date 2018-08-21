
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
from statsmodels.sandbox.stats.multicomp import multipletests
import matplotlib.pyplot as plt
import random

def make_normal_distribution_plot(input_array):
	print("entered make_normal_distribution_plot")
	#figure_name = motif_name.replace("pfm", "png")

	#output_directory = "./plots3"

	mu = np.mean(input_array)
	std = np.std(input_array, ddof = 1)

	print(mu, std)

	fig, ax = plt.subplots()

	plt.hist(input_array, bins = len(input_array), color = 'g') #, normed=True)
	xmin, xmax = plt.xlim()
	x = np.linspace(xmin, xmax, 100)
	p = stats.norm.pdf(x, mu, std)
	plt.plot(x, p, 'r', linewidth = 1)
	motif_name = motif_name.replace(".pfm", "")
	title = motif_name + " fit results: mu = %.4f, std = %.4f" % (mu, std)

	plt.title(title)
	plt.grid()

	plt.ylim(0, 50)
	#plt.xlim(-40, 40)

	#fig.savefig(os.path.join(output_directory, figure_name))
	print("show")
	plt.show()

def main():

	read_file = "test_differences8.txt"
	differences_array = []

	with open(read_file) as r_file:
		for r_line in r_file:
			r_line = re.sub(r'\s+', '', r_line)
			#r_line = r_line[:-1] #delete last character which is ,
			differences_array = re.split(r',', r_line)
			differences_array.pop() #delete last element which is empty because of ,
			
	r_file.close()

	#print(len(differences_array))
	#print(differences_array[0])


	"""
	differences_array.sort()

	for i in differences_array:
		if i.startswith("0.00000000"):
			print(i)
		
	for i in range(40):
		print(differences_array[i])

	print()

	print(max(differences_array))
	"""

	#to_print = int(len(differences_array) / 591)

	#for i in range(to_print):
	#	print(differences_array[i])

	#for i in range(len(differences_array)):
	#	print(differences_array[i])
	#	differences_array[i] = float(differences_array[i])

	#print(differences_array[0])
	#print(type(differences_array[0]))
	#print(float(differences_array[0]))

	differences_array = [float(i) for i in differences_array]
	#differences_array = np.array(differences_array)

	#differences_array.sort()

	print(np.percentile(differences_array, 95))
	print(np.percentile(differences_array, 50)) #median

	print()

	print(max(differences_array))
	print(min(differences_array))

	"""
	#make rectangular box plot
	fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize = (9, 4))
	bplot1 = axes[0].boxplot(differences_array,
							vert = True, #vertical box alignment
							patch_artist = True) #fill with color
	axes[0].set_title('Rectangular box plot')

	blot2 = axes[1].boxplot(differences_array,
							0, '')
							#notch = True,
							#vert = True,
							#patch_artist = True)
	axes[1].set_title('Without outlier points')

	plt.show()
	"""

	differences_array_small = []

	for i in range(len(differences_array)):
		if differences_array[i] <= 40 and differences_array[i] >= -40:
			differences_array_small.append(differences_array[i])

	print(len(differences_array))
	print(len(differences_array_small))

	mu = np.mean(differences_array_small)
	std = np.std(differences_array_small, ddof = 1)

	fig, ax = plt.subplots()

	plt.hist(differences_array_small, bins = len(differences_array_small), color = 'g') #, normed=True)
	xmin, xmax = plt.xlim()
	x = np.linspace(xmin, xmax, 100)
	p = stats.norm.pdf(x, mu, std)
	plt.plot(x, p, 'r', linewidth = 1)
	#motif_name = motif_name.replace(".pfm", "")
	#title = motif_name + " fit results: mu = %.4f, std = %.4f" % (mu, std)

	#plt.title(title)
	plt.grid()

	plt.ylim(0, 50)
	#plt.xlim(-40, 40)

	#fig.savefig(os.path.join(output_directory, figure_name))
	#print("show")
	plt.show()

	
if __name__ == "__main__":
	main()