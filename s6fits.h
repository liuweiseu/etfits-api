#ifndef S6FITS_H_INCLUDED
#define S6FITS_H_INCLUDED

#include <stdio.h>
#include <time.h>
#include <vector>
#include <stdint.h>

typedef struct 
{
  double julian_date;
  time_t unix_time;
  double ra;
  int bors;
  double dec;
  int missedpk;
  float detected_power;
  float mean_power;
  int32_t fine_channel_bin;
  unsigned short coarse_channel_bin;
  double ifreq;
  double rfreq;
} s6hits_t;

typedef struct
{
  double julian_date;
  time_t unix_time;
  double ra;
  int bors;
  double dec;
  int nhits;
  int missedpk;  
} s6hitsheader_t;

typedef struct s6dataspec_t
{
  char * filename;
  //vector where hits requested will be added onto.
  std::vector<s6hits_t> s6hits;
  std::vector<s6hitsheader_t> s6hitsheaders;
  //these can't both be > 0
  int sortby_ifreq;
  int sortby_rfreq;
  int sortby_time;
  int sortby_bors;
  float threshold;
  std::vector<int> bors;
  //same as above but with channels
  std::vector<int> channels;
  //status in cfitsio
  int errorcode;
  s6dataspec_t() : sortby_ifreq(0),
                   sortby_rfreq(0),
                   sortby_time(0),
                   sortby_bors(0),
                   threshold(20.0) {}
} s6dataspec_t;



int get_s6data (s6dataspec_t * s6dataspec);

int get_s6hitsheaders (s6dataspec_t * s6dataspec);

time_t get_time_over_file (char * filename);

int get_hits_over_file (char * filename);

void print_hits_structure (std::vector<s6hits_t> s6hits);

void print_hits_table (std::vector<s6hits_t> s6hits);

void print_hits_header_table (std::vector<s6hitsheader_t> s6hitsheaders);

#endif
