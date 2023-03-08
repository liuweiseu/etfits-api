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

    #print spec.s6hits[0].unix_time, spec.s6hits[1].unix_time

    flen = float(s6fits.get_time_over_file(spec.filename))
    hitnum = int(s6fits.get_hits_over_file(spec.filename))

    base_filename = ntpath.basename(spec.filename)
    data_filename = base_filename+".wf_data"
    plot_filename = base_filename+".wf.png"

    data_f=open(data_filename, 'w')
    start_time = spec.s6hits[0].unix_time
    start_jtime = spec.s6hits[0].julian_date
    for x in range(0, hitnum-1):
        #print >> data_f, spec.s6hits[x].unix_time - start_time, spec.s6hits[x].ifreq
        #print >> data_f, int((spec.s6hits[x].julian_date-start_jtime)*86400), spec.s6hits[x].rfreq
        print(int((spec.s6hits[x].julian_date-start_jtime)*86400), spec.s6hits[x].rfreq, file=data_f)
    data_f.close()

    cmd = 'set term png; set out "%s";                \
            set title "%s" noenhanced; unset key;                \
               set xl "Freq (MHz)"; set yl "Time (Sec)"; \
            plot "%s" using 2:1 with dots;'           \
            % (plot_filename, base_filename, data_filename)
    cmd_f=open("wf_gnuplot_cmds", 'w')
    #print >> cmd_f, cmd
    print(cmd, file=cmd_f)
    cmd_f.close()

    os.system('gnuplot wf_gnuplot_cmds')

    os.remove("wf_gnuplot_cmds")
    os.remove(data_filename)
    

    #print flen,hitnum
    print(flen,hitnum)

def plot_data_with_matplotlib(spec):

    start_time = time.time()
    hitnum = int(s6fits.get_hits_over_file(spec.filename))
    end_time = time.time()
    print("get hit num time was %g seconds" % (end_time - start_time))

    start_time = time.time()
    z = np.dtype(spec.s6hits[0])
    z.names
    end_time = time.time()
    print("make freq array one line was %g seconds" % (end_time - start_time))

    hit_start_time = spec.s6hits[0].unix_time
    #print("start time %s" % (hit_start_time))
    x = np.array([spec.s6hits[0].unix_time-hit_start_time])
    start_time = time.time()
    for i in range(1, hitnum-1):
        x = np.append(x, [spec.s6hits[i].unix_time-hit_start_time])
        #print("time %s x %s" % (spec.s6hits[i].unix_time, x))
    end_time = time.time()
    print("make time array was %g seconds" % (end_time - start_time))

    start_time = time.time()
    y = np.array([spec.s6hits[0].ifreq])
    for i in range(1, hitnum-1):
        y = np.append(y, [spec.s6hits[i].ifreq])
        #print("ifreq %s" % (spec.s6hits[i].ifreq))
    end_time = time.time()
    print("make freq array was %g seconds" % (end_time - start_time))

    #plt.plot(y,x,marker='.',)
    start_time = time.time()
    plt.scatter(y,x,s=.1)
    end_time = time.time()
    print("make plot was %g seconds" % (end_time - start_time))

    start_time = time.time()
    plt.show()
    end_time = time.time()
    print("show plot was %g seconds" % (end_time - start_time))

    #print flen,hitnum
    print(flen,hitnum)

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
    #print spec.s6hits[0].unix_time, spec.s6hits[0].fine_channel_bin, spec.s6hits[0].julian_date
    plot_data_with_gnuplot(spec)
