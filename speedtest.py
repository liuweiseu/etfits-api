import s6fits
import pyfits
import os
from timeit import Timer

#def whole_folder(path):
#  for (dirpath, dirnames, filenames) in os.walk(path):
#    for filename in filenames:
#      if filename.endswith('fits'):
#        file_c = os.sep.join([dirpath, filename])
#        s6time = Timer(lambda: s6fits_test(file_c))
#        print s6time.timeit(number=1)

def s6fits_test(filename):
  time_list = []
  spec = s6fits.s6dataspec_t()
  spec.filename = filename
  s6fits.get_s6data(spec)
  vec = spec.s6hits
  
  for i in vec:
    time_list.append(i.time)

def pyfits_test(filename):
  time_list = []
  hdulist = pyfits.open(filename, memmap=True)
  for hdu in hdulist:
    try:
      time = hdu.header['TIME']
      time_list.append(time)
    except KeyError:
      pass 
  hdulist.close()

if __name__ == '__main__':
  path = "/disks/bucket/b/users/seti/serendip6_data/ao/s6c0/2016"
  filename = "serendip6_eth2_AO_ALFA_1006_20160127_000604.fits"
  s6time = Timer(lambda: s6fits_test(filename))
  print s6time.timeit(number=10)
  pytime = Timer(lambda: pyfits_test(filename))
  print pytime.timeit(number=10)
