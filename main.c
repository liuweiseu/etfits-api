#include <stdio.h>
#include <time.h>
#include <vector>
#include "s6fits.h"


int main(int argc, char *argv[]) 
{
  s6dataspec_t s6dataspec; 
  char * filename = argv[1];
  //std::vector<s6hits_t> hits;
  s6dataspec.filename = filename;
  //s6dataspec.s6hits = hits;
  s6dataspec.sortby_time = 0;
  //std::vector<int> bors;
  //bors.push_back(3);
  //bors.push_back(6);
  //s6dataspec.bors = bors;
  s6dataspec.sortby_bors = 2;
  //we need to be careful here since this will probably break hardcore if
  //someone accidentally sorts by both ifreq and rfreq, or if they try to sort
  //by any frequency when trying to get only headers
  s6dataspec.sortby_ifreq = 1; 
  s6dataspec.sortby_rfreq = 0;
  s6dataspec.threshold = 20.0; 
  get_s6data(&s6dataspec);
  //it's not sorting ;_;
  print_hits_table(s6dataspec.s6hits);
  return 0;
}
