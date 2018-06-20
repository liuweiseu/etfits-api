import numpy as np
#import sys
import argparse
import os
import time
import ntpath
import subprocess as sp 
import fnmatch 
import glob 
import s6fits
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.backends.backend_pdf import PdfPages
import pdb # for debug

def load_file(filename,spec):
#	spec = s6fits.s6dataspec_t()
	spec.filename = filename
	spec.sortby_rfreq = 0
	spec.sortby_time = 0
	spec.sortby_bors = 0
	spec.sortby_ifreq = 0
	# Setting up this mode will remove all those hits which shows change of RF due to 
	# smart frequency switching experiment. It only reads those hits which occur maximum number of time. 
	spec.filterby_rf_center_mode = 1
	#Get hits from file
	start_time = time.time()
	s6fits.get_s6data(spec)
	end_time = time.time()
	print("read '%s' time was %g seconds" % (filename, end_time - start_time))

	hitnum = int(s6fits.get_hits_over_file(spec.filename))
	
	m_ifreq = [] 
	m_unix_time = [] 
	m_meanpower = []
	 
	base_filename = ntpath.basename(spec.filename)
	
	# debug
	#pdb.set_trace()
	
	start_time = spec.s6hits[0].unix_time
#MARK	for i in range(0,hitnum-1):
	for i in range(hitnum):
		if spec.s6hits[i].unix_time == start_time:
			m_ifreq.append(spec.s6hits[i].ifreq)
			m_unix_time.append(spec.s6hits[i].unix_time)
			m_meanpower.append(spec.s6hits[i].mean_power)
		else:	
			break
				
	sort_ifreq_index = np.argsort(m_ifreq)
	sort_ifreq = np.array(m_ifreq)[sort_ifreq_index]
	sort_meanpower = np.array(m_meanpower)[sort_ifreq_index]
	
	# Set break point
	#pdb.set_trace()	
	return sort_meanpower, sort_ifreq , base_filename

def plot_multi_with_gnuplot(filelist):
	spec = s6fits.s6dataspec_t()
	
#	return sort_meanpower[:q], sort_ifreq[:q], base_filename
	sort_meanpower = []
	sort_ifreq = []
	base_filename = ''
	for filename in filelist:
		sort_meanpower_, sort_ifreq_, base_filename_ = load_file(filename,spec)
		sort_meanpower.append(sort_meanpower_)
		sort_ifreq.append(sort_ifreq_)
		base_filename = base_filename + base_filename_ + ','
		# Set a break points
		#pdb.set_trace()

	# Set breakpoint
	#pdb.set_trace()

	meanpower = np.hstack(sort_meanpower)
	ifreq = np.hstack(sort_ifreq)
	base_filename = base_filename[:-1]
	assert len(ifreq.shape) == 1,'ifreq should be of shape (n,)'
	
	data_filename = base_filename+".bb_data"
	plot_filename = base_filename+".bb.png"
	 	
	data_f=open(data_filename, 'w')
	for x in range(len(ifreq)):
		print >> data_f, meanpower[x], ifreq[x]	
	data_f.close()
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


