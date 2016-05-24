import s6fits

a = s6fits.s6dataspec_t()
a.filename = "serendip6_eth2_AO_ALFA_2820_20150609_130436.fits"
s6fits.get_s6data(a)
vec = a.s6hits

for i in vec:
  print (i.ra)
