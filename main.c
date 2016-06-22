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
  //s6dataspec.sortby_freq = 0;
  s6dataspec.sortby_time = 0;
  //std::vector<int> bors;
  //bors.push_back(3);
  //bors.push_back(6);
  //s6dataspec.bors = bors;
  s6dataspec.sortby_bors = 0;
  s6dataspec.threshold = 20.0; 
  get_s6hitsheaders(&s6dataspec);
  print_hits_header_table(s6dataspec.s6hitsheaders);
  //print_hits_table(s6dataspec.s6hits);
  //if (argc == 3) get_num_hits(s6dataspec.s6hits);
  //time_t time = get_time_over_file (filename);
  //int num_hits = get_hits_over_file (filename);
  return 0;
}
