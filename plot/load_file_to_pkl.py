import numpy as np
import argparse
import time
import ntpath
#import subprocess as sp 
#import fnmatch 
#import glob 
import pickle
import s6fits

def load_file(filename):
    spec = s6fits.s6dataspec_t()
    spec.filename = filename
    spec.sortby_rfreq = 0
    spec.sortby_time = 0
    spec.sortby_bors = 0
    spec.sortby_ifreq = 0
    # Setting up this mode will remove all those hits which shows change of RF due to 
    # smart frequency switching experiment. It only reads those hits which occur maximum number of time. 
    spec.filterby_rf_center_mode = 1
    #Get hits from file
    s6fits.get_s6data(spec)

    hitnum = int(s6fits.get_hits_over_file(spec.filename))

    m_ifreq = [] 
    m_unix_time = [] 
    m_meanpower = []
 
    base_filename = ntpath.basename(spec.filename)
    
    start_time = spec.s6hits[0].unix_time
    #MARK for i in range(0,hitnum-1):
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
    
    return sort_meanpower, sort_ifreq , base_filename

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="load partial data from source to target")
    parser.add_argument('-f',type=str, help='filename that will push into gnuplot. e.g. python sort_print.py -f 1.fits')
    parser.add_argument('-t',type=str,default='sort_print.pkl',help='print the value to the target file. e.g. python sort_print.py -t sort_print.pkl')
    args = parser.parse_args()
    filein = args.f
    fileto = args.t

    mean_power, ifreq, _ = load_file(filein)
    
    data = {"mean_power":mean_power, "ifreq":ifreq}
    file_io = open(fileto,'wb')
    pickle.dump(data, file_io)
    file_io.close()
