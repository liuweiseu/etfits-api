#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include "fitsio.h"
#include "s6fits.h"
#include <math.h>
#include <vector>
#include <algorithm>

/*function declarations*/
void get_telescope (fitsfile * fptr, int * status, char * telescope);

time_t get_time (fitsfile * fptr, int * status);

double get_RA (fitsfile * fptr, int * status);

int get_bors (fitsfile * fptr, char * obs, int * status);

double get_dec (fitsfile * fptr, int * status);

int get_nhits (fitsfile * fptr, int * status);

int get_missedpk (fitsfile * fptr, int * status);

int get_coarchid (fitsfile * fptr, int * status);
        
double get_clock_freq (fitsfile * fptr, int * status);

int32_t get_signed_fc (int32_t fc, char * obs);

double get_if1synhz (fitsfile * fptr, int * status);

double get_if2synhz (fitsfile * fptr, int * status);

double calc_ifreq(double clock_freq, int coarchid, char * obs,  
                  int32_t signed_fc, unsigned short cc);

double calc_rfreq(char * telescope, double if1synhz, 
                  double if2synhz, double ifreq);

int is_aoscram(fitsfile * fptr, int * status, 
               char * obs, int * obs_flag);

int is_gbtstatus(fitsfile * fptr, int * status, 
                 char * obs, int * obs_flag);

int is_ethits (fitsfile * fptr, int * status);

int is_desired_bors (int bors, std::vector<int> desired_bors);

int is_desired_coarchan (int cc, std::vector<int> desired_coarchan);

double get_julian_from_unix (int unix_time);

int find (int val, std::vector<int> vec);

void sort_bors(std::vector<s6hits_t> s6hits);

void sort_time(std::vector<s6hits_t> s6hits);

bool cmpbors(const s6hits_t &lhs, const s6hits_t &rhs);

bool cmptime(const s6hits_t &lhs, const s6hits_t &rhs);

void sort(s6dataspec_t * s6dataspec);

/*get_s6data will populate a s6hits_t data type according to the data specs
 * given by the user in the s6dataspec argument. For more details please see
 * s6hits.h, which contains the two data types used here.*/

int get_s6data(s6dataspec_t * s6dataspec) 
{
  fitsfile *fptr;
  /*status must be initialized to zero!*/
  int status = 0;
  int hdupos = 0, nkeys;
  int flag = 0;
  int obs_flag = 0;
  char * filename = s6dataspec->filename;
  int coarchid, clock_freq;
  /*observatory, if another gets added with length > 3 change this*/
  char telescope[FLEN_VALUE];
  char obs[10];
  double if2synhz, if1synhz;
  if (!fits_open_file(&fptr, filename, READONLY, &status))
  {
    //fits_get_hdu_num(fptr, &hdupos); 
    for (; !status; hdupos++) 
    {
      //primary header
      if (hdupos == 0)
      {
        get_telescope(fptr, &status, telescope);
      }
      //checks if the header we're looking at correspond to aoscram or gbt, and
      //if they do, get the coarchid and clock frequency once. Also sets the
      //observatory.
      if (is_aoscram(fptr, &status, obs, &obs_flag) || 
          is_gbtstatus(fptr, &status, obs, &obs_flag))
      {
        if (!flag) 
        {
          coarchid = get_coarchid(fptr, &status);
          clock_freq = get_clock_freq(fptr, &status);
          flag = 1;
        }
        if (strcmp(telescope, "AO_327MHz") == 0)
        {
          if2synhz = get_if2synhz(fptr, &status);
          printf("AO_327MHz");
        }
        if1synhz = get_if1synhz(fptr, &status); 
      }
      /* calculated and saved before we make the hit to save us the hassle of
 * grabbing it for each new hit in the binary table*/
      time_t time = get_time(fptr, &status);
      double ra = get_RA(fptr, &status);
      int bors = get_bors(fptr, obs, &status);
      double dec = get_dec(fptr, &status);
      int missedpk = get_missedpk(fptr, &status); 
      int nhits = get_nhits(fptr, &status);
      if (is_ethits(fptr, &status) && nhits > 0 
          && is_desired_bors(bors, s6dataspec->bors))
      {
        int ncols, anynul;
        long nrows; 

        fits_get_num_rows(fptr, &nrows, &status);
        fits_get_num_cols(fptr, &ncols, &status); 
        
        float detpow_array[nrows];
        fits_read_col (fptr, TFLOAT, 1, 1, 1, nrows, NULL,
            &detpow_array, &anynul, &status); 

        float meanpow_array[nrows];
        fits_read_col (fptr, TFLOAT, 2, 1, 1, nrows, NULL,
            &meanpow_array, &anynul, &status); 

        unsigned short coarch_array[nrows];
        fits_read_col (fptr, TUSHORT, 3, 1, 1, nrows, NULL,
            &coarch_array, &anynul, &status); 
        
        unsigned long finechan_array[nrows];
        fits_read_col (fptr, TULONG, 4, 1, 1, nrows, NULL,
            &finechan_array, &anynul, &status); 

        for (int i = 0; i < nhits; i++) 
        {
          if (is_desired_coarchan(coarch_array[i], s6dataspec->channels))
          {
            s6hits_t hit;
            hit.unix_time = time;
            hit.julian_date = get_julian_from_unix((int) time);
            hit.ra = ra;
            hit.bors = bors;
            hit.dec = dec;
            hit.nhits = nhits;
            hit.missedpk = missedpk;
            hit.detected_power = detpow_array[i];
            hit.mean_power = meanpow_array[i];
            int32_t signed_fc = get_signed_fc(finechan_array[i], obs);
            hit.fine_channel_bin = signed_fc;
            //hit.fine_channel_bin = finechan_array[i];
            hit.coarse_channel_bin = coarch_array[i];
            hit.ifreq = calc_ifreq(clock_freq, coarchid, obs, 
                                   signed_fc, coarch_array[i]);
            s6dataspec->s6hits.push_back(hit);
          }
        } 
      }
      
      /*moves to next HDU*/
      fits_get_hdrspace(fptr, &nkeys, NULL, &status); 
      fits_movrel_hdu(fptr, 1, NULL, &status);
    }  
  }    
  status = 0;
  fits_close_file(fptr, &status);
 
  if (status) fits_report_error(stderr, status);
  
  s6dataspec->errorcode = status; 

  sort(s6dataspec);
  
  return (status);
}

/*gets telescope found in primary header*/
void get_telescope (fitsfile * fptr, int * status, char * telescope)
{
  fits_read_key(fptr, TSTRING, "TELESCOP", telescope, NULL, status);
  //telescope is in the primary header, this should only run once. If there's an
  //error, something is very wrong.
  if (*status) fits_report_error(stderr, *status);
}

/*gets time found in the binary table header block*/
time_t get_time (fitsfile * fptr, int * status)
{
  time_t time;
  fits_read_key(fptr, TINT, "TIME", &time, NULL, status);
  if (*status == KEY_NO_EXIST) time = -1, *status=0;
  return ((time_t) time);
}

/*gets RA found in the binary table header block*/
double get_RA (fitsfile * fptr, int * status) 
{
  double ra;
  fits_read_key(fptr, TDOUBLE, "RA", &ra, NULL, status);
  if (*status == KEY_NO_EXIST) ra = -1, *status=0;
  return (ra);
}

/*gets beampol found in the binary table header block*/
int get_bors (fitsfile * fptr, char * obs, int * status)
{
  int bors;
  fits_read_key(fptr, TINT, "BORSPOL", &bors, NULL, status);
  if (*status == KEY_NO_EXIST) 
  {
    *status = 0;
    if (strcmp(obs, "AO") == 0)
    {
      fits_read_key(fptr, TINT, "BEAMPOL", &bors, NULL, status);
    }
    else if (strcmp(obs, "GBT") == 0)
    {
      /*does this word come up???*/
      fits_read_key(fptr, TINT, "SPECTRA", &bors, NULL, status);
    }
    if (*status == KEY_NO_EXIST) bors = -1, *status=0;
  }
  return (bors);
}

/*gets declination found in the binary table header block*/
double get_dec (fitsfile * fptr, int * status)
{
  double dec;
  fits_read_key(fptr, TDOUBLE, "DEC", &dec, NULL, status);
  if (*status == KEY_NO_EXIST) dec = -1, *status=0;
  return (dec);
}

/*gets nhits found in the binary table header block. Used for testing to see if
 * we should look in the binary table or not and finding total hits in a file*/
int get_nhits (fitsfile * fptr, int * status)
{
  int nhits;
  fits_read_key(fptr, TINT, "NHITS", &nhits, NULL, status);
  if (*status == KEY_NO_EXIST) nhits = -1, *status=0;
  return (nhits);
}

/*gets missed packets found in the binary table header block*/
int get_missedpk (fitsfile * fptr, int * status)
{
  int missedpk;
  fits_read_key(fptr, TINT, "MISSEDPK", &missedpk, NULL, status);
  if (*status == KEY_NO_EXIST) missedpk = -1, *status=0;
  return (missedpk);
}

/*coarchid and clock_freq are found in AOSCRAM or GBTSTATUS. They are both
 * used in calculating intermediate frequency.*/
int get_coarchid (fitsfile * fptr, int * status) 
{
  int coarchid;
  fits_read_key(fptr, TINT, "COARCHID", &coarchid, NULL, status);
  if (*status == KEY_NO_EXIST) coarchid = -1, *status=0;
  return (coarchid);
}
        
double get_clock_freq (fitsfile * fptr, int * status)
{
  int clock_freq;
  fits_read_key(fptr, TINT, "CLOCKFRQ", &clock_freq, NULL, status);
  if (*status == KEY_NO_EXIST) clock_freq = -1, *status=0;
  return (clock_freq);
}

/*if1synhz and if2synhz are used to calculate sky frequency*/
double get_if1synhz (fitsfile * fptr, int * status)
{
  double if1synhz;
  fits_read_key(fptr, TDOUBLE, "IF1SYNHZ", &if1synhz, NULL, status);
  if (*status == KEY_NO_EXIST) if1synhz = -1, *status=0;
  return (if1synhz);
}

double get_if2synhz (fitsfile * fptr, int * status)
{
  double if2synhz;
  fits_read_key(fptr, TDOUBLE, "IF2SYNHZ", &if2synhz, NULL, status);
  if (*status == KEY_NO_EXIST) if2synhz = -1, *status=0;
  return (if2synhz);
}

int32_t get_signed_fc (int32_t fc, char * obs) 
{
  int32_t signed_fc;
  if (strcmp(obs, "AO") == 0)
    signed_fc = (fc & 0x0001FFFF) | ((fc & 0x00010000) ? 0xFFFE0000 : 0);	
  else if (strcmp(obs, "GBT") == 0)
    signed_fc = (fc & 0x0007FFFF) | ((fc & 0x00040000) ? 0xFFF80000 : 0);	
  else signed_fc = -1;
  return signed_fc;
}

/*now works with both ao and gbt. Some finechannel per coarsechannel is
 * hardcoded in depending on which observatory is used. Should another
 * observatory be added, this function will need to be updated*/
double calc_ifreq(double clock_freq, int coarchid, char * obs, 
                  int32_t signed_fc, unsigned short cc)
{
  //differing constant
  int32_t fc_per_cc;
  if (strcmp(obs, "GBT") == 0) fc_per_cc = pow(2.0, 19); 
  else if (strcmp(obs, "AO") == 0) fc_per_cc = pow(2.0, 17);
  int cc_per_sys = 4096;
  double band_width = clock_freq/2;
  double fc_bin_width = band_width/(cc_per_sys * fc_per_cc);
  double resolution = fc_bin_width * 1000000;    
  double sys_cc = coarchid + cc;
  double sys_fc = fc_per_cc * sys_cc + signed_fc;
  double ifreq = ((sys_fc + fc_per_cc) * resolution) / 1000000;   
  return ifreq;
}

/*calculates sky frequency, only works for AO*/
double calc_rfreq(char * telescope, double if1synhz, 
                  double if2synhz, double ifreq)
{
  double rf;
  if (strcmp(telescope, "AO_ALFA") == 0)
    rf = (if1synhz / 1000000) - ifreq;
  else if (strcmp(telescope, "AO_327MHz") == 0)
    rf = (if1synhz / 1000000) - ((if2synhz / 1000000) - ifreq);
  else rf = -1;
  return rf;
}

int is_aoscram (fitsfile * fptr, int * status, char * obs, int * obs_flag)
{
  char extname[FLEN_VALUE];
  fits_read_key(fptr, TSTRING, "EXTNAME", &extname, NULL, status);
  if (strcmp(extname, "AOSCRAM") == 0)
  {
    if (!*obs_flag) 
    {
      strcpy(obs, "AO");
      *obs_flag = 1;
    }
    return 1;
  }
  if (*status == KEY_NO_EXIST) *status=0;
  return 0;
}

int is_gbtstatus (fitsfile * fptr, int * status, char * obs, int * obs_flag)
{
  char extname[FLEN_VALUE];
  fits_read_key(fptr, TSTRING, "EXTNAME", &extname, NULL, status);
  if (strcmp(extname, "GBTSTATUS") == 0)
  {
    if (!*obs_flag) 
    {
      strcpy(obs, "GBT");
      *obs_flag = 1;
    }
    return 1;
  }
  if (*status == KEY_NO_EXIST) *status=0;
  return 0; 
}

/*determine if the hdu we're looking at will include the hits or not*/
int is_ethits (fitsfile * fptr, int * status)
{
  char extname[16];
  fits_read_key(fptr, TSTRING, "EXTNAME", &extname, NULL, status);  
  if (*status == KEY_NO_EXIST) *status=0;
  if (strcmp(extname, "ETHITS") == 0) return 1;
  else return 0; 
}

/*convert from unix time to julian date*/
double get_julian_from_unix (int unix_time)
{
  return (double(unix_time) / 86400.0) + 2440587.5;
}

/*determines if we want to grab hits depending on the beampol/spectra provided
 * by the user. If no specifications are given (empty bors vector), then we'll
 * return 1 for everything*/
int is_desired_bors (int bors, std::vector<int> desired_bors) 
{
  if (desired_bors.empty()) return 1;
  if (find(bors, desired_bors)) return 1;
  else return 0;
}

/*same thing as above but for coarse channels*/
int is_desired_coarchan (int cc, std::vector<int> desired_coarchan) 
{
  if (desired_coarchan.empty()) return 1;
  if (find(cc, desired_coarchan)) return 1;
  else return 0;
}

/*might replace this with std::find since I have to include algorithm for sort
 * anyway*/
int find (int val, std::vector<int> vec)
{
  for (std::vector<int>::iterator it = vec.begin(); it != vec.end(); ++it)
  {
    if (val == * it) return 1;
  }
  return 0;
}

void sort(s6dataspec_t * s6dataspec)
{
  char const* sort_order[3] = {"", "", ""};
  if (s6dataspec->sortby_bors > 0) {
    sort_order[s6dataspec->sortby_bors] = "bors";
  }
  if (s6dataspec->sortby_time > 0) {
    sort_order[s6dataspec->sortby_time] = "time";
  } 
  //if (&s6dataspec->sortby_freq > 0) sort_order[&s6dataspec->sortby_freq] = "freq";
  for (int i=2; i >= 0; i--) 
  {
    if (strcmp(sort_order[i], "bors") == 0)
      std::stable_sort(s6dataspec->s6hits.begin(), 
                       s6dataspec->s6hits.end(), 
                       cmpbors);
    else if (strcmp(sort_order[i], "time") == 0) sort_time(s6dataspec->s6hits);
    //else if (sort_order[i] == "freq") sort_freq(s6dataspec->s6hits);
    else {;}
  } 
}

bool cmpbors(const s6hits_t &lhs, const s6hits_t &rhs)
{
  return lhs.bors < rhs.bors;
}

void sort_time(std::vector<s6hits_t> s6hits)
{
  std::stable_sort(s6hits.begin(), s6hits.end(), cmptime);
}

bool cmptime(const s6hits_t &lhs, const s6hits_t &rhs)
{
  return lhs.unix_time < rhs.unix_time;
}

time_t get_time_over_file (char * filename) 
{
  fitsfile *fptr;
  int status = 0;
  int hdupos = 0, nkeys;
  time_t initial = -1;
  time_t last;
  if (!fits_open_file(&fptr, filename, READONLY, &status))
  { 
    for (; !status; hdupos++) 
    {
      time_t temp_time = get_time(fptr, &status);
      if (temp_time != -1)
      {
        if (initial == -1) initial = temp_time; 
        last = temp_time;
      }
      fits_get_hdrspace(fptr, &nkeys, NULL, &status); 
      fits_movrel_hdu(fptr, 1, NULL, &status);
    }
  }
  status = 0;
  fits_close_file(fptr, &status); 
  if (status) fits_report_error(stderr, status);
  return (last - initial);
}

int get_hits_over_file (char * filename) 
{
  fitsfile *fptr;
  int status = 0;
  int hdupos = 0, nkeys;
  int total_hits = 0;
  if (!fits_open_file(&fptr, filename, READONLY, &status))
  { 
    for (; !status; hdupos++) 
    {
      int temp_hits = get_nhits(fptr, &status);
      if (temp_hits != -1) total_hits += temp_hits;
      fits_get_hdrspace(fptr, &nkeys, NULL, &status); 
      fits_movrel_hdu(fptr, 1, NULL, &status);
    }
  }
  status = 0;
  fits_close_file(fptr, &status); 
  if (status) fits_report_error(stderr, status);
  return total_hits;
}

void print_hits_structure (std::vector<s6hits_t> s6hits) 
{
    for (std::vector<s6hits_t>::iterator hit = s6hits.begin() ; 
         hit != s6hits.end(); ++hit)
    {
      /*header metadata*/
      printf("unix time: %d\n", (int) hit->unix_time);
      printf("julian date: %g\n", (float) hit-> julian_date);
      printf("RA: %g\n", hit->ra);
      printf("beampol: %d\n", hit->bors);
      printf("dec: %g\n", hit->dec);
      printf("missedpk: %d\n", hit->missedpk);
      /*binary table data*/
      printf("detected power: %f\n", hit->detected_power);
      printf("mean power: %f\n", hit->mean_power);
      printf("fine channel bin: %zu\n", (size_t) hit->fine_channel_bin);
      printf("coarse channel bin: %hu\n", hit->coarse_channel_bin);
      /*calculated frequencies*/ 
      printf("ifreq: %g\n", hit->ifreq);
      /*
      printf("rfreq: %g\n", hit->rfreq); 
      */
      printf("\n");
    }
}

void print_hits_table (std::vector<s6hits_t> s6hits)
{
  printf("%15s", "Unix Time");
  printf("%25s", "Julian Date");
  printf("%10s", "RA");
  printf("%10s", "BorS");
  printf("%15s", "DEC");
  printf("%10s", "MSDPK");
  printf("%25s", "DETPOW");
  printf("%20s", "MEANPOW");
  printf("%20s", "FINECHAN");
  printf("%15s", "COARCHID"); 
  printf("%15s\n", "ifreq");
  for (std::vector<s6hits_t>::iterator hit = s6hits.begin() ; 
         hit != s6hits.end(); ++hit)
  {
    printf("%15d", (int) hit->unix_time);
    printf("%25f", hit->julian_date);
    printf("%10g", hit->ra);
    printf("%10d", hit->bors);
    printf("%15g", hit->dec);
    printf("%10d", hit->missedpk);
    printf("%25f", hit->detected_power);
    printf("%20f", hit->mean_power);
    printf("%20d", (int)  hit->fine_channel_bin);
    printf("%15hu", hit->coarse_channel_bin);
    printf("%15f\n", hit->ifreq); 
  }
}

