#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include "fitsio.h"
#include "s6fits.h"
#include <math.h>
#include <vector>
#include <algorithm>

/*function declarations*/
int get_nsubband (fitsfile * fptr, int * status);
void get_telescope (fitsfile * fptr, int * status, char * telescope);
float get_threshold (fitsfile * fptr, int * status);
time_t get_time (fitsfile * fptr, int * status);
double get_RA (fitsfile * fptr, int * status);
int get_bors (fitsfile * fptr, char * obs, int * status);
double get_dec (fitsfile * fptr, int * status);
int get_nhits (fitsfile * fptr, int * status);
int get_missedpk (fitsfile * fptr, int * status);
int get_coarchid (fitsfile * fptr, int * status);
double get_clock_freq (fitsfile * fptr, int * status);
int32_t get_signed_fc (int32_t fc, char * obs);
double calc_ifreq(double clock_freq, int coarchid, char * obs,  
                  int32_t signed_fc, unsigned short cc);
int is_good_data_gb(double ifv1iffq, double rf_reference, double bammpwr1, double bammpwr2);
int is_ethits (fitsfile * fptr, int * status);
int is_ccpowers (fitsfile * fptr, int * status);
int is_desired_bors (int bors, std::vector<int> desired_bors);
int is_desired_coarchan (int cc, std::vector<int> desired_coarchan);
double get_julian_from_unix (int unix_time);
int find (int val, std::vector<int> vec);
void filter(s6dataspec_t * s6dataspec, std::vector<double> &rf_reference_vec);
void sort_bors(std::vector<s6hits_t> s6hits);
void sort_time(std::vector<s6hits_t> s6hits);
void sort_ifreq(std::vector<s6hits_t> s6hits);
void sort_rfreq(std::vector<s6hits_t> s6hits);
bool cmpbors(const s6hits_t &lhs, const s6hits_t &rhs);
bool cmptime(const s6hits_t &lhs, const s6hits_t &rhs);
bool cmpifreq(const s6hits_t &lhs, const s6hits_t &rhs);
bool cmprfreq(const s6hits_t &lhs, const s6hits_t &rhs);
void sort(s6dataspec_t * s6dataspec);

// AO specific
int is_aoscram(fitsfile * fptr, int * status, 
               char * obs, int * obs_flag);
double get_if1synhz (fitsfile * fptr, int * status);
double get_if2synhz (fitsfile * fptr, int * status);
void get_ifv1ssb (fitsfile * fptr, int * status, char * ifv1ssb);
double calc_ao_rfreq(double ifreq, double rf_reference, double if2synhz, char * telescope);

// GBT specific
int is_gbtstatus(fitsfile * fptr, int * status, 
                 char * obs, int * obs_flag);
double get_ifv1csfq (fitsfile * fptr, int * status);
double get_ifv1iffq (fitsfile * fptr, int * status);
double get_bammpwr1 (fitsfile * fptr, int * status);
double get_bammpwr2 (fitsfile * fptr, int * status);
double calc_gbt_rfreq(double ifreq, double rf_reference, double ifv1iffq, char * ifv1ssb);

// FAST & MRO specific
int is_faststatus(fitsfile * fptr, int * status, 
                 char * obs, int * obs_flag);
int is_mrostatus(fitsfile * fptr, int * status, 
                 char * obs, int * obs_flag);
double calc_fast_rfreq(double ifreq, double rf_reference);

// represents an array of booleans
static char _bits[8] = {1, 2, 4, 8, 16, 32, 64, (char)128};
struct BITMAP {
    char* data;
    static inline long get_size(long n) {
        return (n+7)/8;
    }
    void init(char* d) {
        data = d;
    }
    int write_file(const char* filename, long nbits);
    inline void set(long i) {
        data[i/8] |= _bits[i%8];
    }
    inline int get(long i) {
        return data[i/8]&_bits[i%8];
    }
};

//------------------------------------------------------------------------------
int process_primary(fitsfile *fptr, s6dataspec_t * s6dataspec, int &testmode) {
//------------------------------------------------------------------------------

  int status = 0;
  char telescope[FLEN_VALUE];

  get_telescope(fptr, &status, telescope);
  if(strcmp(telescope, "S6TEST") == 0) testmode = 1;
  strncpy(s6dataspec->telescope, telescope, FLEN_VALUE);
  s6dataspec->threshold = get_threshold(fptr, &status);

  return(status);
}

//------------------------------------------------------------------------------
int process_meta_data(int &coarchid, 
					  int &clock_freq, 
					  char * obs, 
  					  fitsfile *fptr,
					  s6dataspec_t * s6dataspec, 
  					  double &rf_reference,
  					  std::vector<double> &rf_reference_vec,
					  int &good_data) {
//------------------------------------------------------------------------------

  double bammpwr1, bammpwr2;
  double if2synhz, if1synhz;
  double ifv1iffq;
  char ifv1ssb[FLEN_VALUE];
  int status = 0;
  int testmode = 0;		// set to 1 for testing, 0 for production
  // grab the meta data that is present for all observatories
  coarchid = get_coarchid(fptr, &status);
  clock_freq = get_clock_freq(fptr, &status);

  //get LO settings for RF calculation and determine if data are good

  if (strcmp(obs, "AO") == 0) {
    if1synhz = get_if1synhz(fptr, &status);	// need for both ALFA and 327MHz 
    if2synhz = get_if2synhz(fptr, &status);	// just need for 327MHz
    rf_reference = if1synhz / 1000000;
    //if (strcmp(telescope, "AO_327MHz") == 0) {
    //		rf_reference = (if1synhz / 1000000) - (if2synhz / 1000000);
    //     } else {		// assume ALFA
	//	  }
	 good_data = 1;	// TODO is assuming good data a good idea for AO?
  }  // end AO 

  else if (strcmp(obs, "GBT") == 0) {
    if(testmode) {
		good_data = 1;		// data are always good in test mode!
    } else {
      	get_ifv1ssb(fptr, &status, ifv1ssb);
       	ifv1iffq = get_ifv1iffq(fptr, &status);
       	rf_reference = get_ifv1csfq(fptr, &status); 
	  	bammpwr1 = get_bammpwr1(fptr, &status);
		bammpwr2 = get_bammpwr2(fptr, &status);
	  	good_data = is_good_data_gb(ifv1iffq, rf_reference, bammpwr1, bammpwr2);
	}
  }  // end GBT 

  else if (strcmp(obs, "FAST") == 0) {
		// no meta data needed for FAST rf_reference, so ne need to test for testmode
		good_data = 1;		//  TODO is assuming good data a good idea for FAST?
		int nyquist_zone = 3;
       	rf_reference = (nyquist_zone-1) * (double)clock_freq/2.0;	
								// For FAST, we sample sky frequency at the third nyquist zone.
								// As an example, with a clock freq of 1000MHz, the third 
								// nyquist zone 1000MHz to 1500MHz converts down to 0 to 500MHz. 
								// It aliases down all by itself and it's not flipped.   
								// TODO - this value should be in the FITS file.
  }  // end FAST
  else if (strcmp(obs, "MRO") == 0) {
		// no meta data needed for FAST rf_reference, so ne need to test for testmode
		good_data = 1;		//  TODO is assuming good data a good idea for FAST?
		int nyquist_zone = 3;
        clock_freq = 1000;
       	//rf_reference = (nyquist_zone-1) * (double)clock_freq/2.0;	 
        fits_read_key(fptr, TDOUBLE, "LO_FREQ", &rf_reference, NULL, &status);
        //status = 0;
								// For FAST, we sample sky frequency at the third nyquist zone.
								// As an example, with a clock freq of 1000MHz, the third 
								// nyquist zone 1000MHz to 1500MHz converts down to 0 to 500MHz. 
								// It aliases down all by itself and it's not flipped.   
								// TODO - this value should be in the FITS file.
  }  // end MRO

  s6dataspec->filterby_rf_reference_mode = 1;		// for testing
  //if(s6dataspec->filterby_rf_reference_mode && rf_reference >= 0) 
  if(s6dataspec->filterby_rf_reference_mode && rf_reference != 0) {
  	rf_reference_vec.push_back(rf_reference);
  }

  return(status);
}

//------------------------------------------------------------------------------
int process_hits(int &coarchid, 
					  int &clock_freq, 
					  char * obs, 
					  char * telescope, 
  					  fitsfile *fptr,
					  s6dataspec_t * s6dataspec, 
  					  double rf_reference,
  					  std::vector<double> &rf_reference_vec,
					  int &good_data) {
//------------------------------------------------------------------------------
  int status = 0;
  int total_missedpk;
  double if2synhz;
  double ifv1iffq;
  char ifv1ssb[FLEN_VALUE];

  time_t time = get_time(fptr, &status);
  double ra = get_RA(fptr, &status);
  int bors = get_bors(fptr, obs, &status);
  double dec = get_dec(fptr, &status);
  int missedpk = get_missedpk(fptr, &status); 
  if (missedpk > 0) {
    total_missedpk += missedpk;
  }
  int nhits = get_nhits(fptr, &status);

//fprintf(stderr, "process_hits telescope %s\n", s6dataspec->telescope); 

  if (nhits > 0 && is_desired_bors(bors, s6dataspec->bors)) {
    int ncols, anynul;
    long nrows; 

    fits_get_num_rows(fptr, &nrows, &status);
    fits_get_num_cols(fptr, &ncols, &status); 
        
    double detpow_array[nrows];
    fits_read_col (fptr, TDOUBLE, 1, 1, 1, nrows, NULL, &detpow_array, &anynul, &status); 

    double meanpow_array[nrows];
    fits_read_col (fptr, TDOUBLE, 2, 1, 1, nrows, NULL, &meanpow_array, &anynul, &status); 

    unsigned short coarch_array[nrows];
    fits_read_col (fptr, TUSHORT, 3, 1, 1, nrows, NULL, &coarch_array, &anynul, &status); 
        
    unsigned long finechan_array[nrows];
    fits_read_col (fptr, TULONG, 4, 1, 1, nrows, NULL, &finechan_array, &anynul, &status); 

    for (int i = 0; i < nhits; i++) {
      if (is_desired_coarchan(coarch_array[i], s6dataspec->channels)) {
    		s6hits_t hit;
        	hit.unix_time = time;
        	hit.julian_date = get_julian_from_unix((int) time);
        	if (strcmp(obs, "GBT") == 0)  
          		hit.ra = ra/15;
        	else  
          		hit.ra = ra;
        	hit.bors = bors;
        	hit.dec = dec;
        	hit.missedpk = missedpk;
        	hit.user_flag = -1;
        	hit.detected_power = detpow_array[i];
        	hit.mean_power = meanpow_array[i];
        	int32_t signed_fc = get_signed_fc(finechan_array[i], obs);
        	hit.fine_channel_bin = signed_fc;
        	hit.coarse_channel_bin = coarch_array[i];
        	hit.ifreq = calc_ifreq(clock_freq, coarchid, obs, signed_fc, coarch_array[i]);
  			//hit.rfreq = rf_reference >= 0 ? rf_reference - hit.ifreq : -1;	// an rf_reference of -1 indicates a header error
			hit.rf_reference = rf_reference;
        	if (strcmp(obs, "AO") == 0) 
{
//fprintf(stderr, "calling telescope %s\n", s6dataspec->telescope);
       			hit.rfreq = calc_ao_rfreq(hit.ifreq, rf_reference, if2synhz, s6dataspec->telescope);
}
        	else if (strcmp(obs, "GBT") == 0) 
              	hit.rfreq = calc_gbt_rfreq(hit.ifreq, rf_reference, ifv1iffq, ifv1ssb);
        	else if (strcmp(obs, "FAST") == 0) 
              	hit.rfreq = calc_fast_rfreq(hit.ifreq, rf_reference);
          else if (strcmp(obs, "MRO") == 0) 
              	hit.rfreq = calc_fast_rfreq(hit.ifreq, rf_reference);
        	else 
				hit.rfreq = -1;
        	s6dataspec->s6hits.push_back(hit);
	  }  // end if(is_desired_coarchan)
    }  // end for(i < nhits)
  }  // end if(nhits > 0 && is_desired_bors()) 

  //fprintf(stderr, "nhits = %ld", nhits);
  return(status);
}

/*get_s6data will populate a s6hits_t data type according to the data specs
 * given by the user in the s6dataspec argument. For more details please see
 * s6hits.h, which contains the two data types used here.*/

//------------------------------------------------------------------------------
int get_s6data(s6dataspec_t * s6dataspec) {
//------------------------------------------------------------------------------
  fitsfile *fptr;
  int status = 0;		// status must be initialized to zero!
  int hdupos = 0, nkeys;
  int testmode = 0;		// set to 1 for testing, 0 for production
  int obs_flag = 0;
  int good_data = 0;
  int total_missedpk;
  char * filename = s6dataspec->filename;
  int coarchid, clock_freq;
  /*observatory, if another gets added with length > 3 change this*/
  char telescope[FLEN_VALUE];
  char obs[10];
  double rf_reference;
  std::vector<double> rf_reference_vec;
  if (!fits_open_file(&fptr, filename, READONLY, &status)) {
//fprintf(stderr, "pre main loop %s\n", telescope);
    //fits_get_hdu_num(fptr, &hdupos); 
    for (; !status; hdupos++) {

      if (hdupos == 0) {
		status = process_primary(fptr, s6dataspec, testmode);
      }

      else if (is_aoscram    (fptr, &status, obs, &obs_flag) || 
          	   is_gbtstatus  (fptr, &status, obs, &obs_flag) ||
          	   is_faststatus (fptr, &status, obs, &obs_flag) ||
               is_mrostatus  (fptr, &status, obs, &obs_flag)) {
		status = process_meta_data(	coarchid, 
						  		   	clock_freq, 
						  			obs, 
  						  			fptr,
						  			s6dataspec, 
  					  			   	rf_reference,
  					  			   	rf_reference_vec,
					  			   	good_data); 

	  }
      else if (is_ethits(fptr, &status) && good_data) {
		status = process_hits(	coarchid, 
						  		   	clock_freq, 
						  			obs, 
									telescope,
  						  			fptr,
						  			s6dataspec, 
  					  			   	rf_reference,
  					  			   	rf_reference_vec,
					  			   	good_data); 
      }  // end if(is_ethits() && good_data)
      
      /*moves to next HDU*/
      fits_get_hdrspace(fptr, &nkeys, NULL, &status); 
      fits_movrel_hdu(fptr, 1, NULL, &status);
    }  // end for(; !status; hdupos++)
  }  // end if(!fits_open_file(&fptr, filename, READONLY, &status))

  status = 0;
  fits_close_file(fptr, &status);
 
  if (status) fits_report_error(stderr, status);
   
  s6dataspec->total_missedpk = total_missedpk; 
  s6dataspec->ciftsio_error = status; 

  if (status > 0) s6dataspec->errorcode += 1;

  filter(s6dataspec, rf_reference_vec);
  sort(s6dataspec);
  
  return (status);
}  // end int get_s6data(s6dataspec_t * s6dataspec)


//------------------------------------------------------------------------------
int get_s6ccpowers(s6dataspec_t * s6dataspec) 
//------------------------------------------------------------------------------
{
  fitsfile *fptr;
  /*status must be initialized to zero!*/
  int status = 0;
  int hdupos = 0, nkeys;
  int header_flag = 0;
  int obs_flag = 0;
  char * filename = s6dataspec->filename;
  int coarchid, clock_freq;
  /*observatory, if another gets added with length > 3 change this*/
  char telescope[FLEN_VALUE];
  char obs[10];
  int nsubband;
  //double if2synhz, if1synhz;
  //double ifv1iffq, ifv1csfq;
  //char ifv1ssb[FLEN_VALUE];
  if (!fits_open_file(&fptr, filename, READONLY, &status))
  {
    //fits_get_hdu_num(fptr, &hdupos); 
    for (; !status; hdupos++) 
    {
      //primary header
      if (hdupos == 0)
      {
        get_telescope(fptr, &status, telescope);
        s6dataspec->threshold = get_threshold(fptr, &status);
        nsubband = get_nsubband (fptr, &status);
      }
      //checks if the header we're looking at correspond to aoscram or gbt, and
      //if they do, get the coarchid and clock frequency once. Also sets the
      //observatory.
      if (is_aoscram(fptr, &status, obs, &obs_flag) || 
          is_gbtstatus(fptr, &status, obs, &obs_flag))
      {
        if (!header_flag) 
        {
          coarchid = get_coarchid(fptr, &status);
          clock_freq = get_clock_freq(fptr, &status);
          header_flag = 1;
        }
        /*if (strcmp(obs, "AO") == 0) 
        {
          if (strcmp(telescope, "AO_327MHz") == 0)
          {
            if2synhz = get_if2synhz(fptr, &status);
          }
          if1synhz = get_if1synhz(fptr, &status); 
        }
        else if (strcmp(obs, "GBT") == 0) 
        {
          get_ifv1ssb(fptr, &status, ifv1ssb);
          ifv1iffq = get_ifv1iffq(fptr, &status);
          ifv1csfq = get_ifv1csfq(fptr, &status); 
        }*/
      }
      /* calculated and saved before we make the hit to save us the hassle of
 * grabbing it for each new hit in the binary table*/
      time_t time = get_time(fptr, &status);
      double ra = get_RA(fptr, &status);
      double dec = get_dec(fptr, &status);
      if (is_ccpowers(fptr, &status))
      {
        int ncols, anynul;
        long nrows; 

        fits_get_num_rows(fptr, &nrows, &status);
        fits_get_num_cols(fptr, &ncols, &status);        

        double polx_array[nrows];
        fits_read_col (fptr, TDOUBLE, 1, 1, 1, nrows, NULL,
            &polx_array, &anynul, &status); 

        double poly_array[nrows];
        fits_read_col (fptr, TDOUBLE, 2, 1, 1, nrows, NULL,
            &poly_array, &anynul, &status); 

        for (int i = 0; i < nrows; i++) 
        {
          s6ccpowers_t ccpowers;
          ccpowers.unix_time = time;
          ccpowers.julian_date = get_julian_from_unix((int) time);
          ccpowers.ra = ra;
          ccpowers.dec = dec;
          if (strcmp(obs, "AO") == 0) 
          {
            ccpowers.beam = i/(nsubband-1);
            ccpowers.coarse_channel_bin = i - (nsubband * ccpowers.beam);
          }
          else
          { 
            ccpowers.beam = 1;
            ccpowers.coarse_channel_bin = i;
          }
          ccpowers.power_x = polx_array[i];
          ccpowers.power_y = poly_array[i];
          /*hit.ifreq = calc_ifreq(clock_freq, coarchid, obs, 
                                 signed_fc, coarch_array[i]);
          if (strcmp(obs, "AO") == 0) 
            hit.rfreq = calc_ao_rfreq(telescope, if1synhz, 
                                      if2synhz, hit.ifreq);
          else if (strcmp(obs, "GBT") == 0) 
            hit.rfreq = calc_gbt_rfreq(hit.ifreq, ifv1csfq, 
                                       ifv1iffq, ifv1ssb);
          else hit.rfreq = -1;*/
          s6dataspec->s6ccpowers.push_back(ccpowers);
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
  
  s6dataspec->ciftsio_error = status; 
  if (status > 0) s6dataspec->errorcode += 1;

  sort(s6dataspec);
  
  return (status);
}

//------------------------------------------------------------------------------
int get_s6hitsheaders(s6dataspec_t * s6dataspec) 
//------------------------------------------------------------------------------
{
  fitsfile *fptr;
  /*status must be initialized to zero!*/
  int status = 0;
  int hdupos = 0, nkeys;
  char * filename = s6dataspec->filename;
  int obs_flag = 0;
  char obs[10];
  if (!fits_open_file(&fptr, filename, READONLY, &status))
  {
    for (; !status; hdupos++) 
    {
      /* these below functions are necessary just to grab the observatory for
       * bors */
      is_aoscram(fptr, &status, obs, &obs_flag);
      is_gbtstatus(fptr, &status, obs, &obs_flag);
      int bors = get_bors(fptr, obs, &status);
      if (is_ethits(fptr, &status) && is_desired_bors(bors, s6dataspec->bors))
      {
        s6hitsheader_t header;
        header.unix_time = get_time(fptr, &status);
        header.julian_date = get_julian_from_unix((int) header.unix_time);
        header.ra = get_RA(fptr, &status);
        header.bors = bors;
        header.dec = get_dec(fptr, &status);
        header.missedpk = get_missedpk(fptr, &status); 
        header.nhits = get_nhits(fptr, &status);
        
        s6dataspec->s6hitsheaders.push_back(header);
      }
      /* moves to next hdu */
      fits_get_hdrspace(fptr, &nkeys, NULL, &status); 
      fits_movrel_hdu(fptr, 1, NULL, &status);
    }  
  }    
  status = 0;
  fits_close_file(fptr, &status);
 
  if (status) fits_report_error(stderr, status);
  
  s6dataspec->ciftsio_error = status; 
  if (status > 0) s6dataspec->errorcode += 1;

  sort(s6dataspec);
  
  return (status);
}

//delete later
//------------------------------------------------------------------------------
int get_nsubband (fitsfile * fptr, int * status) 
//------------------------------------------------------------------------------
{
  int nsubs;
  fits_read_key(fptr, TINT, "NSUBBAND", &nsubs, NULL, status);
  if (*status == KEY_NO_EXIST) nsubs = -1, *status=0;
  return (nsubs);
}

//------------------------------------------------------------------------------
float get_threshold (fitsfile * fptr, int * status) 
//------------------------------------------------------------------------------
{
  float threshold;
  fits_read_key(fptr, TFLOAT, "THRSHOLD", &threshold, NULL, status);
  if (*status == KEY_NO_EXIST) threshold = 20.0, *status=0;
  return (threshold);
}

/*gets telescope found in primary header*/
//------------------------------------------------------------------------------
void get_telescope (fitsfile * fptr, int * status, char * telescope)
//------------------------------------------------------------------------------
{
  fits_read_key(fptr, TSTRING, "TELESCOP", telescope, NULL, status);
//fprintf(stderr, "TELESCOP %s\n", telescope);
  //telescope is in the primary header, this should only run once. If there's an
  //error, something is very wrong.
  if (*status) fits_report_error(stderr, *status);
}

/*gets time found in the binary table header block*/
//------------------------------------------------------------------------------
time_t get_time (fitsfile * fptr, int * status)
//------------------------------------------------------------------------------
{
  time_t time;
  fits_read_key(fptr, TINT, "TIME", &time, NULL, status);
  if (*status == KEY_NO_EXIST) time = -1, *status=0;
  return ((time_t) time);
}

/*gets RA found in the binary table header block*/
//------------------------------------------------------------------------------
double get_RA (fitsfile * fptr, int * status) 
//------------------------------------------------------------------------------
{
  double ra;
  fits_read_key(fptr, TDOUBLE, "RA", &ra, NULL, status);
  if (*status == KEY_NO_EXIST) ra = -1, *status=0;
  return (ra);
}

/*gets beampol found in the binary table header block*/
//------------------------------------------------------------------------------
int get_bors (fitsfile * fptr, char * obs, int * status)
//------------------------------------------------------------------------------
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
//------------------------------------------------------------------------------
double get_dec (fitsfile * fptr, int * status)
//------------------------------------------------------------------------------
{
  double dec;
  fits_read_key(fptr, TDOUBLE, "DEC", &dec, NULL, status);
  if (*status == KEY_NO_EXIST) dec = -1, *status=0;
  return (dec);
}

/*gets nhits found in the binary table header block. Used for testing to see if
 * we should look in the binary table or not and finding total hits in a file*/
//------------------------------------------------------------------------------
int get_nhits (fitsfile * fptr, int * status)
//------------------------------------------------------------------------------
{
  int nhits;
  fits_read_key(fptr, TINT, "NHITS", &nhits, NULL, status);
  if (*status == KEY_NO_EXIST) nhits = -1, *status=0;
  return (nhits);
}

/*gets missed packets found in the binary table header block*/
//------------------------------------------------------------------------------
int get_missedpk (fitsfile * fptr, int * status)
//------------------------------------------------------------------------------
{
  int missedpk;
  fits_read_key(fptr, TINT, "MISSEDPK", &missedpk, NULL, status);
  if (*status == KEY_NO_EXIST) missedpk = -1, *status=0;
  return (missedpk);
}

/*coarchid and clock_freq are found in AOSCRAM or GBTSTATUS. They are both
 * used in calculating intermediate frequency.*/
//------------------------------------------------------------------------------
int get_coarchid (fitsfile * fptr, int * status) 
//------------------------------------------------------------------------------
{
  int coarchid;
  fits_read_key(fptr, TINT, "COARCHID", &coarchid, NULL, status);
  if (*status == KEY_NO_EXIST) coarchid = -1, *status=0;
  return (coarchid);
}
        
//------------------------------------------------------------------------------
double get_clock_freq (fitsfile * fptr, int * status)
//------------------------------------------------------------------------------
{
  int clock_freq;
  fits_read_key(fptr, TINT, "CLOCKFRQ", &clock_freq, NULL, status);
  if (*status == KEY_NO_EXIST) clock_freq = -1, *status=0;
  return (clock_freq);
}

/*if1synhz and if2synhz are used to calculate sky frequency*/
//------------------------------------------------------------------------------
double get_if1synhz (fitsfile * fptr, int * status)
//------------------------------------------------------------------------------
{
  double if1synhz;
  fits_read_key(fptr, TDOUBLE, "IF1SYNHZ", &if1synhz, NULL, status);
  if (*status == KEY_NO_EXIST) if1synhz = -1, *status=0;
  return (if1synhz);
}

//------------------------------------------------------------------------------
double get_if2synhz (fitsfile * fptr, int * status)
//------------------------------------------------------------------------------
{
  double if2synhz;
  fits_read_key(fptr, TDOUBLE, "IF2SYNHZ", &if2synhz, NULL, status);
  if (*status == KEY_NO_EXIST) if2synhz = -1, *status=0;
  return (if2synhz);
}

//------------------------------------------------------------------------------
int32_t get_signed_fc (int32_t fc, char * obs) 
//------------------------------------------------------------------------------
{
  // signed fine channels should go as :
  // 	max neg..........0..........max pos
  // unsigned fine channels should go as :
  //    0---------------------------max pos 	
  //
 
  int32_t signed_fc;

  // 		     select desired bits                     make int  leave int
  //                  2**N-1        "our" sign bit neg ? negative  positive
  if (strcmp(obs, "AO") == 0)
    signed_fc = (fc & 0x0001FFFF) | ((fc & 0x00010000) ? 0xFFFE0000 : 0);	// 2**17 signed fc's 	
  else if (strcmp(obs, "GBT") == 0)
    signed_fc = (fc & 0x0007FFFF) | ((fc & 0x00040000) ? 0xFFF80000 : 0);	// 2**19 signed fc's	
  else if (strcmp(obs, "FAST") == 0)
    signed_fc = fc;															// fc's are unsigned
    //signed_fc = (fc & 0x0FFFFFFF) | ((fc & 0x08000000) ? 0xF0000000 : 0);	// 2**28 signed fc's	
  else if (strcmp(obs, "MRO") == 0)
    signed_fc = fc;		
  else signed_fc = -1;
  return signed_fc;
}

//------------------------------------------------------------------------------
double get_ifv1csfq (fitsfile * fptr, int * status)
//------------------------------------------------------------------------------
{
  double ifv1csfq;
  fits_read_key(fptr, TDOUBLE, "IFV1CSFQ", &ifv1csfq, NULL, status);
  if (*status == KEY_NO_EXIST) ifv1csfq = -1, *status=0;
  return (ifv1csfq);
}

//------------------------------------------------------------------------------
double get_ifv1iffq (fitsfile * fptr, int * status)
//------------------------------------------------------------------------------
{
  double ifv1iffq;
  fits_read_key(fptr, TDOUBLE, "IFV1IFFQ", &ifv1iffq, NULL, status);
  if (*status == KEY_NO_EXIST) ifv1iffq = -1, *status=0;
  return (ifv1iffq);
}

//------------------------------------------------------------------------------
void get_ifv1ssb (fitsfile * fptr, int * status, char * ifv1ssb)
//------------------------------------------------------------------------------
{
  fits_read_key(fptr, TSTRING, "IFV1SSB", ifv1ssb, NULL, status);
  if (*status) fits_report_error(stderr, *status);
}

//------------------------------------------------------------------------------
double get_bammpwr1 (fitsfile * fptr, int * status)
//------------------------------------------------------------------------------
{
  double bammpwr1;
  fits_read_key(fptr, TDOUBLE, "BAMMPWR1", &bammpwr1, NULL, status);
  if (*status) fits_report_error(stderr, *status);
  return(bammpwr1);
}

//------------------------------------------------------------------------------
double get_bammpwr2 (fitsfile * fptr, int * status)
//------------------------------------------------------------------------------
{
  double bammpwr2;
  fits_read_key(fptr, TDOUBLE, "BAMMPWR2", &bammpwr2, NULL, status);
  if (*status) fits_report_error(stderr, *status);
  return(bammpwr2);
}


/*now works with both ao and gbt. Some finechannel per coarsechannel is
 * hardcoded in depending on which observatory is used. Should another
 * observatory be added, this function will need to be updated*/
//------------------------------------------------------------------------------
double calc_ifreq(double clock_freq, int coarchid, char * obs, 
                  int32_t signed_fc, unsigned short cc)
//------------------------------------------------------------------------------
{
  //differing constant
  int32_t fc_per_cc;
  double cc_per_sys;
//fprintf(stderr, "obs %s clock %lf coarchid %d signed_fc %ld cc %d\n", obs, clock_freq, coarchid, signed_fc, cc);
  if (strcmp(obs, "GBT") == 0) 
  {
    cc_per_sys = 4096;
    fc_per_cc = pow(2.0, 19); 
    clock_freq = clock_freq * 2;
  }
  else if (strcmp(obs, "AO") == 0) 
  { 
    cc_per_sys = 4096;
    fc_per_cc = pow(2.0, 17);
  }
  else if (strcmp(obs, "FAST") == 0) 
  { 
    cc_per_sys = 1;
    fc_per_cc = pow(2.0, 27);		// 512M real channels becomes 256M complex channels
  }
  else if  (strcmp(obs, "MRO") == 0)
  {
    cc_per_sys = 1;
    fc_per_cc = pow(2.0, 29);		// 512M real channels becomes 256M complex channels
  }
  double band_width 	= clock_freq/2;							// MHz
  double fc_bin_width 	= band_width/(cc_per_sys * fc_per_cc);	// MHz
  double resolution 	= fc_bin_width * 1000000;   			// Hz 
  long sys_cc 			= coarchid + cc;
  long sys_fc 			= fc_per_cc * sys_cc + signed_fc;
  double ifreq 			= (sys_fc * resolution) / 1000000;   	// MHz

//fprintf(stderr, "data from %s total BW %lf FC BW %lf MHz resolution %lf Hz this CC %ld this FC %ld IF %lf\n", obs, band_width,fc_bin_width*1e6,resolution,sys_cc,sys_fc,ifreq);

  return ifreq;
}

/*calculates sky frequency for AO*/
//------------------------------------------------------------------------------
double calc_ao_rfreq(double ifreq, double rf_reference, double if2synhz, char * telescope)
//------------------------------------------------------------------------------
//double calc_ao_rfreq(char * telescope, double if1synhz, 
//                  double if2synhz, double ifreq)
{
  double rf;
  if (strcmp(telescope, "AO_ALFA") == 0)
  	rf = rf_reference >= 0 ? rf_reference - ifreq : -1;
  else if (strcmp(telescope, "AO_327MHz") == 0)
	rf = rf_reference >= 0 ? rf_reference - ((if2synhz / 1000000) - ifreq) : -1;
//fprintf(stderr, "telescope %s rf_reference %lf ifreq %lf rf %lf\n", telescope, rf_reference, ifreq, rf);
  return rf;   // -1 indicates an error
}

//------------------------------------------------------------------------------
double calc_gbt_rfreq(double ifreq, double rf_reference, 
                      double ifv1iffq, char * ifv1ssb) 
//------------------------------------------------------------------------------
{
  double rf;
  if (strcmp(ifv1ssb, "upper") == 0)
    rf = rf_reference + (ifreq - ifv1iffq);
  else if (strcmp(ifv1ssb, "lower") == 0)
    rf = rf_reference - (ifreq - ifv1iffq); 
  else 
    rf = -1;
  return rf;
}

//------------------------------------------------------------------------------
double calc_fast_rfreq(double ifreq, double rf_reference)
//------------------------------------------------------------------------------
{
  return(ifreq + rf_reference);
}

//------------------------------------------------------------------------------
int is_good_data_gb(double ifv1iffq, double rf_reference, double bammpwr1, double bammpwr2) {
//------------------------------------------------------------------------------

	int rv = 1;		// assume good until proven bad

	// if either of these is zero, there was not a usable 
	// signal at acquisition time
	if(ifv1iffq == 0) 		rv = 0;
	if(rf_reference == 0) 	rv = 0;

	// Make sure input power levels are in range
	//if(bammpwr1 < -30.0 || bammpwr1 > -10.0) rv = 0;		// TODO why commented out?
	if(bammpwr2 < -30.0 || bammpwr2 > -10.0) rv = 0;

	return(rv);
}

//------------------------------------------------------------------------------
int is_aoscram (fitsfile * fptr, int * status, char * obs, int * obs_flag)
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
int is_gbtstatus (fitsfile * fptr, int * status, char * obs, int * obs_flag)
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
int is_faststatus (fitsfile * fptr, int * status, char * obs, int * obs_flag)
//------------------------------------------------------------------------------
{

  char extname[FLEN_VALUE];
  fits_read_key(fptr, TSTRING, "EXTNAME", &extname, NULL, status);
  if (strcmp(extname, "FASTSTATUS") == 0)
  {
    if (!*obs_flag) 
    {
      strcpy(obs, "FAST");
      *obs_flag = 1;
    }
    return 1;
  }
  if (*status == KEY_NO_EXIST) *status=0;
  return 0; 
}

//------------------------------------------------------------------------------
int is_mrostatus (fitsfile * fptr, int * status, char * obs, int * obs_flag)
//------------------------------------------------------------------------------
{

  char extname[FLEN_VALUE];
  fits_read_key(fptr, TSTRING, "EXTNAME", &extname, NULL, status);
  if (strcmp(extname, "MROSTATUS") == 0)
  {
    if (!*obs_flag) 
    {
      strcpy(obs, "MRO");
      *obs_flag = 1;
    }
    return 1;
  }
  if (*status == KEY_NO_EXIST) *status=0;
  return 0; 
}

/*determine if the hdu we're looking at will include the hits or not*/
//------------------------------------------------------------------------------
int is_ethits (fitsfile * fptr, int * status)
//------------------------------------------------------------------------------
{
  char extname[16];
  fits_read_key(fptr, TSTRING, "EXTNAME", &extname, NULL, status);  
  if (*status == KEY_NO_EXIST) *status=0;
  if (strcmp(extname, "ETHITS") == 0) return 1;
  else return 0; 
}

//------------------------------------------------------------------------------
int is_ccpowers (fitsfile * fptr, int * status)
//------------------------------------------------------------------------------
{
  char extname[16];
  fits_read_key(fptr, TSTRING, "EXTNAME", &extname, NULL, status);  
  if (*status == KEY_NO_EXIST) *status=0;
  if (strcmp(extname, "CCPWRS") == 0) return 1;
  else return 0; 
}

/*convert from unix time to julian date*/
//------------------------------------------------------------------------------
double get_julian_from_unix (int unix_time)
//------------------------------------------------------------------------------
{
  return (double(unix_time) / 86400.0) + 2440587.5;
}

/*determines if we want to grab hits depending on the beampol/spectra provided
 * by the user. If no specifications are given (empty bors vector), then we'll
 * return 1 for everything*/
//------------------------------------------------------------------------------
int is_desired_bors (int bors, std::vector<int> desired_bors) 
//------------------------------------------------------------------------------
{
  if (desired_bors.empty()) return 1;
  if (find(bors, desired_bors)) return 1;
  else return 0;
}

/*same thing as above but for coarse channels*/
//------------------------------------------------------------------------------
int is_desired_coarchan (int cc, std::vector<int> desired_coarchan) 
//------------------------------------------------------------------------------
{
  if (desired_coarchan.empty()) return 1;
  if (find(cc, desired_coarchan)) return 1;
  else return 0;
}

/*might replace this with std::find since I have to include algorithm for sort
 * anyway*/
//------------------------------------------------------------------------------
int find (int val, std::vector<int> vec)
//------------------------------------------------------------------------------
{
  for (std::vector<int>::iterator it = vec.begin(); it != vec.end(); ++it)
  {
    if (val == * it) return 1;
  }
  return 0;
}

//------------------------------------------------------------------------------
void filter(s6dataspec_t * s6dataspec, std::vector<double> &rf_reference_vec) {
//------------------------------------------------------------------------------

	double rf_reference_mode, rf_reference_mode_test;
	unsigned long mode_count=1, mode_count_test=1;
	unsigned long i;

	BITMAP hit_filter;		// will be used by all filters

	if(!s6dataspec->filterby_rf_reference_mode) return;
    hit_filter.init((char*)calloc(1, BITMAP::get_size(s6dataspec->s6hits.size())));

	// start rf center mode filter
	std::sort(rf_reference_vec.begin(), rf_reference_vec.end());
	rf_reference_mode_test = rf_reference_vec.at(1);
	for(i=1; i < rf_reference_vec.size(); i++) {
		if (rf_reference_vec[i] == rf_reference_mode_test) {
			mode_count_test++;
		} else {
			if(mode_count_test > mode_count) { 	// we have a better mode
				mode_count 			= mode_count_test;			
				rf_reference_mode 		= rf_reference_vec[i-1];
				rf_reference_mode_test = rf_reference_vec[i];
			}
				mode_count_test = 1;
		}
	}
	if(mode_count_test > mode_count) {		// catch final rf_reference_vec[]
		rf_reference_mode = rf_reference_vec[i-1];
	}
	s6dataspec->rf_reference_mode = rf_reference_mode;
	// now go through s6dataspec.s6hits, filtering any element with rf_reference != mode
	for(unsigned long i=1; i < s6dataspec->s6hits.size(); i++) {
		if(s6dataspec->s6hits[i].rf_reference != rf_reference_mode) {
			hit_filter.set(i);
		}
	}
	// end rf center mode filter

	// Additional filters (with bit map setting) go here...
	
	
	// filtering complete, now copy any element not filtered to new vector
	std::vector<s6hits_t> filtered_hits;
	for(unsigned long i=1; i < s6dataspec->s6hits.size(); i++) {
		if(hit_filter.get(i) == 0) {		// this hit OK, ie not filtered
			filtered_hits.push_back(s6dataspec->s6hits[i]);
		}
	}
			
	// swap non-filtered for filtered
	s6dataspec->s6hits.swap(filtered_hits);
}

//------------------------------------------------------------------------------
void sort(s6dataspec_t * s6dataspec)
//------------------------------------------------------------------------------
{
  char const* sort_order[3] = {"", "", ""};
  if (s6dataspec->sortby_bors > 0) {
    sort_order[s6dataspec->sortby_bors] = "bors";
  }
  if (s6dataspec->sortby_time > 0) {
    sort_order[s6dataspec->sortby_time] = "time";
  } 
  if (s6dataspec->sortby_ifreq > 0) {
    sort_order[s6dataspec->sortby_ifreq] = "ifreq";
  }
  if (s6dataspec->sortby_rfreq > 0) {
    sort_order[s6dataspec->sortby_rfreq] = "rfreq";
  }
  for (int i=2; i >= 0; i--) 
  {
    if (strcmp(sort_order[i], "bors") == 0) 
      std::stable_sort(s6dataspec->s6hits.begin(), 
                       s6dataspec->s6hits.end(), 
                       cmpbors);    
    else if (strcmp(sort_order[i], "time") == 0) 
      std::stable_sort(s6dataspec->s6hits.begin(), 
                       s6dataspec->s6hits.end(), 
                       cmptime);    
    else if (strcmp(sort_order[i], "ifreq") == 0) 
      std::stable_sort(s6dataspec->s6hits.begin(), 
                       s6dataspec->s6hits.end(), 
                       cmpifreq);    
    else if (strcmp(sort_order[i], "rfreq") == 0) 
      std::stable_sort(s6dataspec->s6hits.begin(), 
                       s6dataspec->s6hits.end(), 
                       cmprfreq);    
    else {;}
  } 
}

//none of these comparison functions actually work
//------------------------------------------------------------------------------
bool cmpbors(const s6hits_t &lhs, const s6hits_t &rhs)
//------------------------------------------------------------------------------
{
  return lhs.bors < rhs.bors;
}

//------------------------------------------------------------------------------
void sort_bors(std::vector<s6hits_t> s6hits)
//------------------------------------------------------------------------------
{
  std::stable_sort(s6hits.begin(), s6hits.end(), cmpbors);
}

//------------------------------------------------------------------------------
bool cmptime(const s6hits_t &lhs, const s6hits_t &rhs)
//------------------------------------------------------------------------------
{
  return lhs.unix_time < rhs.unix_time;
}

//------------------------------------------------------------------------------
void sort_time(std::vector<s6hits_t> s6hits)
//------------------------------------------------------------------------------
{
  std::stable_sort(s6hits.begin(), s6hits.end(), cmptime);
}

//------------------------------------------------------------------------------
bool cmpifreq(const s6hits_t &lhs, const s6hits_t &rhs)
//------------------------------------------------------------------------------
{
  return lhs.ifreq < rhs.ifreq;
}

//------------------------------------------------------------------------------
void sort_ifreq(std::vector<s6hits_t> s6hits)
//------------------------------------------------------------------------------
{
  std::stable_sort(s6hits.begin(), s6hits.end(), cmpifreq);
}

//------------------------------------------------------------------------------
bool cmprfreq(const s6hits_t &lhs, const s6hits_t &rhs)
//------------------------------------------------------------------------------
{
  return lhs.rfreq < rhs.rfreq;
}

//------------------------------------------------------------------------------
void sort_rfreq(std::vector<s6hits_t> s6hits)
//------------------------------------------------------------------------------
{
  std::stable_sort(s6hits.begin(), s6hits.end(), cmprfreq);
}

//------------------------------------------------------------------------------
time_t get_time_over_file (char * filename) 
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
int get_hits_over_file (char * filename) 
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
void print_hits_structure (std::vector<s6hits_t> s6hits) 
//------------------------------------------------------------------------------
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
      printf("detected power: %20.0f\n", hit->detected_power);
      printf("mean power: %20.0f\n", hit->mean_power);
      printf("fine channel bin: %zu\n", (size_t) hit->fine_channel_bin);
      printf("coarse channel bin: %hu\n", hit->coarse_channel_bin);
      /*calculated frequencies*/ 
      printf("ifreq: %g\n", hit->ifreq);
      printf("rfreq: %g\n", hit->rfreq); 
      printf("\n");
    }
}

//------------------------------------------------------------------------------
void print_hits_table (std::vector<s6hits_t> s6hits)
//------------------------------------------------------------------------------
{
  printf("%15s", "Unix Time");
  printf("%20s", "Julian Date");
  printf("%10s", "RA");
  printf("%15s", "DEC");
  printf("%10s", "BorS");
  printf("%10s", "MSDPK");
  printf("%20s", "DETPOW");
  printf("%24s", "MEANPOW");
  printf("%22s", "SNR");
  printf("%13s", "FINECHAN");
  printf("%12s", "COARCHAN"); 
  printf("%11s", "ifreq");
  printf("%15s", "rfreq");
  printf("%17s\n", "rf_reference");
  for (std::vector<s6hits_t>::iterator hit = s6hits.begin() ; 
         hit != s6hits.end(); ++hit)
  {
    printf("%15d", (int) hit->unix_time);
    printf("%20f", hit->julian_date);
    printf("%10g", hit->ra);
    printf("%15g", hit->dec);
    printf("%10d", hit->bors);
    printf("%10d", hit->missedpk);
    printf("%22.0f", hit->detected_power);
    printf("%22.0f", hit->mean_power);
    printf("%22.0f", hit->detected_power/hit->mean_power);
    printf("%13d", (int)  hit->fine_channel_bin);
    printf("%10hu", hit->coarse_channel_bin);
    printf("%15f", hit->ifreq); 
    printf("%15f", hit->rfreq); 
    printf("%15f\n", hit->rf_reference); 
  }
}

//------------------------------------------------------------------------------
void print_ccpowers_table (std::vector<s6ccpowers_t> s6ccpowers)
//------------------------------------------------------------------------------
{
  printf("%15s", "Unix Time");
  printf("%20s", "Julian Date");
  printf("%10s", "RA");
  printf("%15s", "DEC");
  printf("%10s", "COARCHAN"); 
  printf("%20s", "POWER XPOL");
  printf("%20s\n", "POWER YPOL");
  //printf("%15s", "ifreq");
  //printf("%15s\n", "rfreq");
  for (std::vector<s6ccpowers_t>::iterator ccpowers = s6ccpowers.begin() ; 
         ccpowers != s6ccpowers.end(); ++ccpowers)
  {
    printf("%15d", (int) ccpowers->unix_time);
    printf("%20f", ccpowers->julian_date);
    printf("%10g", ccpowers->ra);
    printf("%15g", ccpowers->dec);
    printf("%10d", (int)ccpowers->coarse_channel_bin);
    printf("%20f", ccpowers->power_x);
    printf("%20f\n", ccpowers->power_y);
  }
}

//------------------------------------------------------------------------------
void print_hits_header_table (std::vector<s6hitsheader_t> s6hitsheaders)
//------------------------------------------------------------------------------
{
  printf("%15s", "Unix Time");
  printf("%20s", "Julian Date");
  printf("%10s", "RA");
  printf("%10s", "BorS");
  printf("%10s", "NHITS");
  printf("%15s", "DEC");
  printf("%10s\n", "MSDPK");
  for (std::vector<s6hitsheader_t>::iterator header = s6hitsheaders.begin() ; 
         header != s6hitsheaders.end(); ++header)
  {
    printf("%15d", (int) header->unix_time);
    printf("%20f", header->julian_date);
    printf("%10g", header->ra);
    printf("%10d", header->bors);
    printf("%10d", header->nhits);
    printf("%15g", header->dec);
    printf("%10d\n", header->missedpk);
  }
}

