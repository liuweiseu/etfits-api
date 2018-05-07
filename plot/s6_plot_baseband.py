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
			
	#flen = float(s6fits.get_time_over_file(spec.filename))
	hitnum = int(s6fits.get_hits_over_file(spec.filename))
	
	base_filename = ntpath.basename(spec.filename)
	data_filename = base_filename+".bb_data"
	plot_filename = base_filename+".bb.png"	
	
	m_ifreq = [] 
	m_unix_time = [] 
	m_meanpower = [] 

	start_time = spec.s6hits[0].unix_time
	for i in range(0,hitnum-1):
		if spec.s6hits[i].unix_time == start_time:
			m_ifreq.append(spec.s6hits[i].ifreq)
			m_unix_time.append(spec.s6hits[i].unix_time)
			m_meanpower.append(spec.s6hits[i].mean_power)
		else:	
			break
				
	sort_ifreq_index = np.argsort(m_ifreq)
	sort_ifreq = np.array(m_ifreq)[sort_ifreq_index]
	sort_meanpower = np.array(m_meanpower)[sort_ifreq_index]
	
	data_f=open(data_filename, 'w')
	for x in range(0,len(m_ifreq)-1):
		print >> data_f, sort_meanpower[x], sort_ifreq[x]	
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
	os.remove(data_filename)

	#print flen,hitnum

# return value added: ifreq and unix_time
def plot_data_with_matplotlib(spec):

	start_time = time.time()
	hitnum = int(s6fits.get_hits_over_file(spec.filename))
	end_time = time.time()
	print("get hit num time was %g seconds" % (end_time - start_time))

	m_ifreq = [] # new
	m_unix_time = [] # new
	m_meanpower = [] # new
	
	start_time = time.time()
	starttime = spec.s6hits[0].unix_time
	for i in range(0,hitnum-1):
		if spec.s6hits[i].unix_time == starttime:
			m_ifreq.append(spec.s6hits[i].ifreq)
			m_meanpower.append(spec.s6hits[i].mean_power)
        		m_unix_time.append(spec.s6hits[i].unix_time)
		else:
			break
	sort_ifreq_index = np.argsort(m_ifreq)
	sort_ifreq = np.array(m_ifreq)[sort_ifreq_index]
	sort_meanpower = np.array(m_meanpower)[sort_ifreq_index]

	#unique_unix_time, unique_index = np.unique(m_unix_time, return_index=True)
	#unique_ifreq = np.array(m_ifreq)[unique_index]
  	#unique_meanpower = np.array(m_meanpower)[unique_index]
	#unique_ifreq_index = np.argsort(unique_ifreq)
	#sort_unique_ifreq = np.array(unique_ifreq)[unique_ifreq_index]	
	#sort_unique_meanpower = np.array(unique_meanpower)[unique_ifreq_index]
	end_time = time.time()
	
	print("make snr array was %g seconds" % (end_time - start_time))

    #plt.plot(unix_time, ifreq,marker='.-',)
	start_time = time.time()
	#df.plot(kind='line',x="unix_time", y="ifreq")
	#plt.plot(unique_unix_time, unique_ifreq, '.-')
	plt.plot(sort_ifreq, sort_meanpower, '.-')
	plt.xlabel('ifreq')
	plt.ylabel('mean power')
	end_time = time.time()
	print("make plot was %g seconds" % (end_time - start_time))

#	plt.xticks(np.arange(min(x), max(x)+1, 10, dtype=int))

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
	
	# Uncomment this line below.
	#plot_data_with_matplotlib(spec)
	plot_data_with_gnuplot(spec)

#	print('-'*10 + '\n' + 'Information of `spec`:\n')
#	print("spec.s6hits: length={}, dtype={}\n".format(len(spec.s6hits), type(spec.s6hits)))
#	print("spec.s6hits: attributes = {}\n".format(dir(spec.s6hits[0])))
#	print("spec.s6hits:\n{}\n".format(spec.s6hits[0]))

