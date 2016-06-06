import s6fits
import pyfits
import os
import FitsData
import time
from timeit import Timer

def whole_folder(path):
  file_list = []
  for (dirpath, dirnames, filenames) in os.walk(path):
    for filename in filenames:
      if filename.endswith('fits'):
        file_c = os.sep.join([dirpath, filename])
        file_list = get_size_and_hits(file_c, file_list)
  return file_list

def get_size_and_hits(filename, file_list):
  hits_size_file = []
  spec = s6fits.s6dataspec_t()
  spec.filename = filename  
  s6fits.get_s6data(spec)
  total_hits = 0
  for hit in spec.s6hits:
    total_hits += 1 
  hits_size_file = [filename, 
                    (os.stat(filename).st_size/(1024*1024)),
                    total_hits]
  file_list.append(hits_size_file)  
  return file_list

def s6fits_test(filename):
  spec = s6fits.s6dataspec_t()
  spec.filename = filename
  s6fits.get_s6data(spec)

def pure_python_test(filename):
  data = FitsData.FitsData(filename) 

if __name__ == '__main__':
  path = "/disks/bucket/b/users/seti/serendip6_data/ao/s6c0/2016"
  filename = "serendip6_eth2_AO_ALFA_1006_20160127_000604.fits"
  #run_count = 10
  #s6time = Timer(lambda: s6fits_test(filename))
  #print "s6 time (average of 10): " + str(s6time.timeit(run_count)/run_count)
  #pytime = Timer(lambda: pure_python_test(filename))
  #print "python time (average of 10): " + str(pytime.timeit(number=run_count)/run_count)
  #get_size_and_hits(filename)
  file_list = whole_folder(path)
  file_list.sort(key=lambda x: x[1])
  for i in file_list:
    print i[0]
    print "file size: " + str(i[1])
    print "number of hits: " + str(i[2])
