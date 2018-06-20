import numpy as np
#import sys
import argparse
import os
import time
import ntpath
#import subprocess as sp 
#import fnmatch 
#import glob 
import s6fits
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.backends.backend_pdf import PdfPages
import pdb # for debug

def plot_multi_with_gnuplot(filelist):
#	return sort_meanpower[:q], sort_ifreq[:q], base_filename
	base_filename = ','.join([os.path.split(f)[-1] for f in filelist])
	data_filename = base_filename+".bb_data"
	plot_filename = base_filename+".bb.png"
	for filename in filelist:
		os.system('python load_file.py -f {} -t {}'.format(filename, data_filename))
	
	# Reload data from `data_filename`, sort adn write
	with open(data_filename,'r') as f:
		data = f.read()
		power, ifreq = [], []
	#	pdb.set_trace()
		for line in data.split('\n'):
			if line=='':continue
			pstr, fstr = line.split()
			power.append(float(pstr))
			ifreq.append(float(fstr))

	# Resort
	ifreq_sort_indx = np.argsort(ifreq)
	power = np.array(power)[ifreq_sort_indx]
	ifreq = np.array(ifreq)[ifreq_sort_indx]
	# Overwrite `data_filename`
	with open(data_filename,'w') as f:
		for p, i in zip(power, ifreq):
			print >> f, p, i	
 	
	cmd = 'set term png; set out "%s";                \
			set title "%s"; unset key;                \
		   	set xl "ifreq"; set yl "Mean Power"; \
			plot "%s" using 2:1 with line;'           \
			% (plot_filename, base_filename, data_filename)
	cmd_f=open("baseband_gnuplot_cmds", 'w')
	print >> cmd_f, cmd
	cmd_f.close()
	os.system('gnuplot baseband_gnuplot_cmds')
	os.remove("baseband_gnuplot_cmds")
#MARK	os.remove('gnuplot_batchplot.bb_data')



if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Plot multiple files with gnuplot")
	parser.add_argument('--files',metavar='N',type=str, nargs='+',help='filenames that will push into gnuplot. e.g. python s6_plot.py --files 1.fits 2.fits 3.fits')
	parser.add_argument('--folder',type=str,default='',help='if set, plot the files in the folder. This keyword will overwrite the filenames in "--files" keyword')
	args = parser.parse_args()
	filelist = args.files
	folder = args.folder
	if len(folder)>0:
		filelist = os.listdir(folder)
#	print(filelist)
	plot_multi_with_gnuplot(filelist)


