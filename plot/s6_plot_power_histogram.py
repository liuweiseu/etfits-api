import numpy as np
import sys
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

def plot_data_with_gnuplot(spec):
	##hitrec = spec.s6hits
	##s6fits.print_hits_table(spec.s6hits)
	#spechdr = spec.s6hitsheaders

	#print spec.total_missedpk
	#print spec.rf_center_mode,spec.filter_rf_center_mode,spec.cfitsio_error
	os.system('gnuplot pow_gnuplot_cmds')

	flen = float(s6fits.get_time_over_file(spec.filename))
	hitnum = int(s6fits.get_hits_over_file(spec.filename))

	base_filename = ntpath.basename(spec.filename)
	data_filename = base_filename+".pow_data"
	plot_filename = base_filename+".pow.png"

	histo = [0] * 100
	data_f=open(data_filename, 'w')
	for i in range(1, hitnum-1):
		snr = spec.s6hits[i].detected_power/spec.s6hits[i].mean_power
		print >> data_f, snr
	
	# see http://psy.swansea.ac.uk/staff/carter/gnuplot/gnuplot_frequency.htm
	cmd = 'set term png; set out "%s";                						\
			set boxwidth 0.05 absolute ; set style fill solid 1.0 noborder; \
			set title "%s"; unset key;                						\
		   	set xl "SNR"; set yl "number of hits"; 							\
																			\
			bin_width  = 0.1; bin_number(x) = floor(x/bin_width); 			\
			rounded(x) = bin_width * ( bin_number(x) + 0.5 );				\
																			\
			plot "%s" using (rounded($1)):(1) smooth frequency with boxes;' \
			% (plot_filename, base_filename, data_filename)
	cmd_f=open("pow_gnuplot_cmds", 'w')
	print >> cmd_f, cmd
	cmd_f.close()

	os.system('gnuplot pow_gnuplot_cmds')

	#os.remove("pow_gnuplot_cmds")
	#os.remove(data_filename)

	print flen,hitnum

def plot_data_with_matplotlib(spec):

	start_time = time.time()
	hitnum = int(s6fits.get_hits_over_file(spec.filename))
	end_time = time.time()
	print("get hit num time was %g seconds" % (end_time - start_time))

	snr = spec.s6hits[0].detected_power/spec.s6hits[0].mean_power
	x = np.array([snr])
	start_time = time.time()
	for i in range(1, hitnum-1):
		snr = spec.s6hits[i].detected_power/spec.s6hits[i].mean_power
		x = np.append(x, snr)
		#print("SNR %f" % (snr))
		#print("det %s" % (spec.s6hits[i].detected_power))
	end_time = time.time()
	print("make snr array was %g seconds" % (end_time - start_time))

    #plt.plot(y,x,marker='.',)
	start_time = time.time()
	plt.hist(x, bins=1000)
	end_time = time.time()
	print("make plot was %g seconds" % (end_time - start_time))

	plt.xticks(np.arange(min(x), max(x)+1, 10, dtype=int))

	start_time = time.time()
	plt.show()
	end_time = time.time()
	print("show plot was %g seconds" % (end_time - start_time))


if __name__ == "__main__":
	spec = s6fits.s6dataspec_t()

	spec.filename = sys.argv[1]

	#Set various data spec 
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
	print("read time was %g seconds" % (end_time - start_time))
	#s6fits.get_s6hitsheaders(spec)

	#plot_data_with_matplotlib(spec)
	plot_data_with_gnuplot(spec)
