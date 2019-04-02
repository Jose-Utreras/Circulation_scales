#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include <string.h>

#ifndef max
    #define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
    #define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

char* concat(const char *s1, const char *s2)
  {
      char *result = malloc(strlen(s1)+strlen(s2)+1);
      strcpy(result, s1);
      strcat(result, s2);
      return result;
  }
int main(int argc, char* argv[]) {

int N_bins, j, lower, upper, N_datos, k, Number;
double c_min, c_max, c_sig, h, upper_c, lower_c, circ, pij;

char* name;
FILE *circ_vector;
FILE *histogram;
//sscanf(argv[1],"%s", &name);
name=argv[1];
char* circulation = concat("Files/", name);
circulation       = concat(circulation,"_circ_vector.txt");
char* hist_file   = concat("Files/", name);
hist_file         = concat(hist_file,"_hist_data.txt");
circ_vector  = fopen(circulation, "r");
histogram    = fopen(hist_file,"w");
sscanf(argv[2], "%lf", &c_min);
sscanf(argv[3], "%lf", &c_max);
sscanf(argv[4], "%d", &N_bins);
sscanf(argv[5], "%lf", &h);
sscanf(argv[6], "%lf", &c_sig);
sscanf(argv[7], "%d", &N_datos);
int *bins = malloc(sizeof(int) * N_bins);
float *errors = malloc(sizeof(float) * N_bins);
for(k=0 ; k < N_bins ; k++){
    bins[k]=0;
    errors[k]=0.0;
}
int l=0;
do{
    fscanf(circ_vector, "%lf", &circ);
    lower   = (int)((circ-c_min)/h);
    upper   = lower+1;
    bins[lower]++;
    Number++;
    for(k=0 ; k < N_bins ; k++){
        lower_c = c_min+(c_max-c_min)*k/(1.0*N_bins);
        upper_c = c_min+(c_max-c_min)*(k+1)/(1.0*N_bins);
        pij     = 0.5*erf((upper_c-circ)/(sqrt(2.0)*c_sig));
        pij     = pij-0.5*erf((lower_c-circ)/(sqrt(2.0)*c_sig));
        errors[k]+=pij*(1.0-pij);}

    N_datos--;
  }while(N_datos>0);
  for(k=0 ; k < N_bins ; k++){
      fprintf(histogram, "%d\t %f\n",bins[k],sqrt(errors[k]));
  }
fclose(circ_vector);
if(histogram==NULL)
	fclose(histogram);
return 0;
}
