#include <string.h>
#include <stdio.h>
#include "cfitsio/fitsio.h"

/*cycles through and prints each keyword in a given header unit*/
int getHeaderbyHDU(char* filename, int hdunum) 
{
	fitsfile *fptr;
	//char card[FLEN_CARD];
	char keyname[FLEN_KEYWORD];
	char value[FLEN_VALUE];
	char comment[FLEN_COMMENT];
	int status = 0;
	int single = 1, nkeys, ii;
	
	if (!fits_open_file(&fptr, filename, READONLY, &status)) 
	{	
		fits_movabs_hdu(fptr, hdunum, IMAGE_HDU, &status);

		fits_get_hdrspace(fptr, &nkeys, NULL, &status); // gets number of keywords 

		printf("Header listing for HDU #%d:\n", hdunum);
		
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

/*finds and returns a particular keyword given an HDU*/
int getKeywordbyHDU(char* filename, int hdunum, char* keyword) 
{
	fitsfile *fptr;
	//char card[FLEN_CARD];
	int status = 0;
	int single = 1, nkeys, ii;
	char value[FLEN_VALUE];
	
	if (!fits_open_file(&fptr, filename, READONLY, &status)) 
	{	
		fits_movabs_hdu(fptr, hdunum, IMAGE_HDU, &status);

		fits_get_hdrspace(fptr, &nkeys, NULL, &status); // gets number of keywords 

		printf("Header listing for HDU #%d:\n", hdunum);
		
		fits_read_keyword(fptr, keyword, value, NULL, &status);	
		printf("%s: %s\n", keyword, value);
		
		
	}
	status = 0;
	
	fits_close_file(fptr, &status);
	
	if (status) fits_report_error(stderr, status);
	
	return(status);
}



/*simply returns the number of hits in each hdu for a given file and prints them*/
int get_hits(char* filename) {
    fitsfile *fptr;         /* FITS file pointer, defined in fitsio.h */
    char card[FLEN_CARD];   /* Standard string lengths defined in fitsio.h */
    int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
    int single = 0, hdupos, nkeys, ii;
    int nhits;

    if (!fits_open_file(&fptr, filename, READONLY, &status))
    {
    	fits_get_hdu_num(fptr, &hdupos);  /* Get the current HDU position */

    	/* List only a single header if a specific extension was given */ 
    	if (hdupos != 1 || strchr(filename, '[')) single = 1;

    	for (; !status; hdupos++)  /* Main loop through each extension */
    	{	
		fits_read_key(fptr, TINT, "NHITS", &nhits, NULL, &status);
		printf("%d\n", nhits);
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

int fits_to_ascii(char* filename) {
    fitsfile *fptr;
    char card[FLEN_CARD]; 
    int status = 0;  // CFITSIO status value MUST be initialized to zero! 
    int single = 0, hdupos, nkeys, ii;
    char value;
    const char* prim_header_keywords[] = {"OBSFREQ"};
    const char* header_keywords[] = {"TIME", "RA", "DEC", "BEAMPOL"};
    //const char* table_keywords[] = {"DETPOW", "MEANPOW", "COARCHID"};
    int num_prim_header_keywords = 1;	
    int num_header_keywords = 4; 
		//int num_table_keywords = 3;
    /*to be included: Compute Node*/
	
    if (!fits_open_file(&fptr, filename, READONLY, &status))
    {
    	fits_get_hdu_num(fptr, &hdupos);  // Get the current HDU position 

    	// List only a single header if a specific extension was given  
    	if (hdupos != 1 || strchr(filename, '[')) single = 1;

    	for (; !status; hdupos++)  // Main loop through each extension 
    	{
		/*for (int i = 0; i < num_prim_header_keywords; i++) 
		{
			fits_read_keyword(fptr, prim_header_keywords[i], &value, NULL, &status);
	
			printf("%s", prim_header_keywords[i]);
			printf(" = %s\n", &value);
	*/
//				fits_get_hdrspace(fptr, &nkeys, NULL, &status); // <-- gets # of keywords.  Necessary for some reason.*/
				for (int j = 0; j < 1; j++) 							
				{ 
					fits_read_keyword(fptr, header_keywords[j], &value, NULL, &status);
				
					printf("%s", header_keywords[j]);
					//printf("%s", "TIME");

					printf(" = %s\n", &value);
	
					if (status == KEY_NO_EXIST)  status = 0;
 				}

				fits_get_hdrspace(fptr, &nkeys, NULL, &status); // <-- gets # of keywords.  Necessary for some reason.*/
	
				fits_movrel_hdu(fptr, 1, NULL, &status);  // try to move to next HDU 

		//}
			}

    	if (status == END_OF_FILE)  status = 0; // Reset after normal error 

    	fits_close_file(fptr, &status);
    
    }
    
    if (status) fits_report_error(stderr, status); // print any error message
    return(status);

}


/*for each instance of the keyword that appears in the fits file, print the value on a new line.  Note that all keyword values are of type char* */
int printAllInstancesofKeyword(char* filename,char* keyword) 
{
    fitsfile *fptr;         /* FITS file pointer, defined in fitsio.h */
    //char card[FLEN_CARD];   /* Standard string lengths defined in fitsio.h */
    int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
    int single = 0, hdupos, nkeys, ii;
    char value[FLEN_VALUE];
    
    if (!fits_open_file(&fptr, filename, READONLY, &status))
    {
    	fits_get_hdu_num(fptr, &hdupos);  /* Get the current HDU position */
	
    	for (; !status; hdupos++)  /* Main loop through each extension */
    	{
				//printf("Header listing for HDU #%d:\n", hdupos);
				fits_get_hdrspace(fptr, &nkeys, NULL, &status); /* <-- gets # of keywords.  Necessary for some as of yet undetermined reason. */
				fits_read_keyword(fptr, keyword, value, NULL, &status);
				printf("%s = %s\n", keyword, &value);	

				if (status == KEY_NO_EXIST)  status = 0;      

				fits_movrel_hdu(fptr, 1, NULL, &status);  /* try to move to next HDU */
    	}

    	if (status == END_OF_FILE)  status = 0; /* Reset after normal error */

    	fits_close_file(fptr, &status);
    
    }
    
    if (status) fits_report_error(stderr, status); /* print any error message */
    return(status);
}

/*linked list definitions*/
typedef struct node 
{
	int *val;
	struct node *next;
} node_t;

void append(node_t* head, int val) 
{
	node_t* curr = head;

	while (curr->next != NULL) 
	{
		curr = curr->next;
	}
	
	curr->next = malloc(sizeof(node_t));
	curr->next->val = val;
	curr->next->next = NULL;
}

void print_list(node_t* head) 
{
	node_t* curr = head;
	
	while (curr != NULL) 
	{
		printf("%d\n", curr->val);
		curr = curr->next;
	}
}

/*saves all values of a single given keyword in a linked list.  Note this can only be used for keywords that return a type int*/
int saveAllInstancesofKeyword(char* filename,char* keyword, node_t* head) {
    fitsfile *fptr;         /* FITS file pointer, defined in fitsio.h */
    char card[FLEN_CARD];   /* Standard string lengths defined in fitsio.h */
    int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
    int single = 0, hdupos, nkeys, ii;
    int value;	

    if (!fits_open_file(&fptr, filename, READONLY, &status))
    {
    	fits_get_hdu_num(fptr, &hdupos);  /* Get the current HDU position */

    	/* List only a single header if a specific extension was given */ 
    	if (hdupos != 1 || strchr(filename, '[')) single = 1;

    	for (; !status; hdupos++)  /* Main loop through each extension */
    	{	
		fits_read_key(fptr, TINT, keyword, &value, NULL, &status);
		
		append(head, value);

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
	/*node_t* head = NULL;
	head = malloc(sizeof(node_t));
	saveAllInstancesofKeyword(argv[1], argv[2], head);
	print_list(head);*/
	//get_hits(argv[1]);
	//printf("testing just printing all instances of a keyword\n");
 	printAllInstancesofKeyword(argv[1], argv[2]);
	/*printf("testing printing primaryHeader\n");
	getPrimaryHeader(argv[1]);*/
	//getHeaderbyHDU(argv[1], atoi(argv[2]));
	//getKeywordbyHDU(argv[1], atoi(argv[2]), argv[3]);
	//fits_to_ascii(argv[1]);	
	return 0;
}

