#include <string.h>
#include <stdio.h>
#include "cfitsio/fitsio.h"

/*grabs and prints out keywords from a given fits files*/

int getPrimaryHeader(char* filename) {
	fitsfile *fptr;
	//char card[FLEN_CARD];
	char keyname[FLEN_KEYWORD];
	char value[FLEN_VALUE];
	char comment[FLEN_COMMENT];
	int status = 0;
	int single = 1, hdupos, nkeys, ii;
	
	if (!fits_open_file(&fptr, filename, READONLY, &status)) 
	{
		fits_get_hdu_num(fptr, &hdupos);
		
		fits_get_hdrspace(fptr, &nkeys, NULL, &status); /* gets number of keywords */

		printf("Header listing for HDU #%d:\n", hdupos);
		
		for (ii = 1; ii < nkeys; ii++) 
		{
			if (fits_read_keyn(fptr, ii, keyname, value, comment, &status)) break;	
			printf("%s: %s\n", keyname, value);
		}
		
	}
	status = 0;
	
	fits_close_file(fptr, &status);
	
	if (status) fits_report_error(stderr, status);
	
	return(status);
}

int printAllInstancesofKeyword(char* filename,char* keyword) {
    fitsfile *fptr;         /* FITS file pointer, defined in fitsio.h */
    char card[FLEN_CARD];   /* Standard string lengths defined in fitsio.h */
    int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
    int single = 0, hdupos, nkeys, ii;
    char value;
    
    if (!fits_open_file(&fptr, filename, READONLY, &status))
    {
    	fits_get_hdu_num(fptr, &hdupos);  /* Get the current HDU position */

    	/* List only a single header if a specific extension was given */ 
    	if (hdupos != 1 || strchr(filename, '[')) single = 1;

    	for (; !status; hdupos++)  /* Main loop through each extension */
    	{
		printf("Header listing for HDU #%d:\n", hdupos);
	
		fits_read_keyword(fptr, keyword, &value, NULL, &status);
		
		printf("%s", keyword);
		printf(" = %s\n", &value);
	
		if (status == KEY_NO_EXIST)  status = 0;
      
		fits_get_hdrspace(fptr, &nkeys, NULL, &status); /* <-- gets # of keywords.  Necessary for some as of yet undetermined reason. */

		fits_movrel_hdu(fptr, 1, NULL, &status);  /* try to move to next HDU */
    	}

    	if (status == END_OF_FILE)  status = 0; /* Reset after normal error */

    	fits_close_file(fptr, &status);
    
    }
    
    if (status) fits_report_error(stderr, status); /* print any error message */
    return(status);

}



int main(int argc, char *argv[])
{
	/*printf("testing just printing all instances of a keyword\n");
        printAllInstancesofKeyword(argv[1], argv[2]);*/
	printf("testing printing primaryHeader\n");
	getPrimaryHeader(argv[1]);
 	return 0;
}

