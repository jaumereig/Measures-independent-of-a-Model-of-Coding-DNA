#include <stdio.h>
#include <math.h>
// gcc Fourier.c -o fffffff -lm
// ./fffffff 120 10 < ../../Measures_Jaume/fasta/MT066176.1.tbl > ../../Measures_Jaume/output/FOU_cpp_120 10_out
#define MAXWINDOW 10000
#define PERIOD 3
#define SEQLENGTH 1000000
#define IDLENGTH 300

main(int argc, char *argv[]){
  char Sequence[SEQLENGTH];
  long   LengthSequence;

  char id[IDLENGTH];

  int WindowLength,
      HalfWindowLength;
  int StepSize;

  double FourierExon (char *, int);

  long i;

  WindowLength=atoi(argv[1]);
  StepSize=atoi(argv[2]);

  HalfWindowLength=WindowLength/2;

  while (scanf("%s %s", id, Sequence) != EOF) {
    LengthSequence = strlen(Sequence);
    // printf("%d\n", LengthSequence);
    for (i=0; i<LengthSequence;i+=StepSize)
      printf ("%d %d %5.3f\n", i+1, i+WindowLength, FourierExon(Sequence+i,WindowLength));
      // printf ("%d %d %5.3f\n", i+1, i+WindowLength);
/*    printf ("\n"); */
  }
}

/*----------------------------------------------------------*
 * FOURIER
 * Fourier transform of the autocorrelation function
 *   for period 3. (adapted from Jim Fickett's code)
 *-----------------------------------------------------------*/
double FourierExon(char *window,int wlen)
{



   int diff;
   int end;
   int autoc[MAXWINDOW];
   int i;
   double fou;

  
  /*----- Autocorrelation function  and Fourier Transform*/
   fou=0.0;

   for ( diff=0; diff<wlen; ++diff ) {  
     autoc[diff] = 0;
     for ( end=diff; end<wlen; ++end )
       if ( window[end]==window[end-diff] )
	      autoc[diff]++;
     fou += autoc[diff]*cos(2*3.14159265*diff/PERIOD);
     //fou += autoc[diff]*(2*3.14159265*diff/PERIOD);
   }
   return fou/wlen;
}
