#include <stdio.h>
#include <time.h>
#include <vector>
#include "s6fits.h"

int main(int argc, char *argv[]) 
{
  s6dataspec_t s6dataspec; 
  //s6dataspec.filename = const_cast<char*>("/disks/bucket/b/users/seti/serendip6_data/ao/s6c0/2016/serendip6_eth2_AO_ALFA_1006_20160127_000604.fits");
  s6dataspec.filename = const_cast<char*>("serendip6_eth2_AO_ALFA_2820_20150609_130436.fits");  
  std::vector<s6hits_t> hits;
  s6dataspec.s6hits = hits;
  s6dataspec.sortby_freq = 0;
  s6dataspec.sortby_time = 2;
  std::vector<int> bors;
  bors.push_back(3);
  bors.push_back(6);
  s6dataspec.bors = bors;
  s6dataspec.sortby_bors = 1;
  s6dataspec.threshold = 0; 
  get_s6data(&s6dataspec);
  print_hits_table(s6dataspec.s6hits);
}
