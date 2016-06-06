import s6fits
import sys

def run_s6(filename):
  spec = s6fits.s6dataspec_t()
  spec.filename = filename  
  s6fits.get_s6data(spec)

if __name__ == '__main__':
  #filename = "serendip6_eth2_AO_ALFA_1006_20160127_000604.fits"
  filename = sys.argv[1]
  run_s6(filename)
  
