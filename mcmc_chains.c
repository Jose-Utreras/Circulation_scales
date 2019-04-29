#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include <string.h>
#include <mpi.h>

#define steps 10
#define PI 3.14159265
#define tiny 1e-16
#define fmin 0.1
#define fmax 100.0
#define send_data_tag 2001
#define finish_data_tag 2002
#define pcyrtokms 0.9778

#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif


char* concat(const char *s1, const char *s2)
  {
      char *result = malloc(strlen(s1)+strlen(s2)+1);//+1 for the null-terminator
      //in real code you would check for errors in malloc here
      strcpy(result, s1);
      strcat(result, s2);
      return result;
  }

double uniform(double minimo, double maximo){
  return minimo+(maximo-minimo)*rand()/(1.0*RAND_MAX);
}

double walker(char *name, char *snapshot,double n1_min , double n1_max, double n2_min , double n2_max, double pc_min, double pc_max,
  double dv_min, double dv_max,
  double pmin, double pmax, int Nbins, double dx,float *sigma, float *delta, float *error, int nchain, int n_procs, int id){
    printf("This is the %03d\t sampler \n ", id);


    FILE *output;

    char outname[sizeof "_chain100.txt"];
    sprintf(outname, "_chain%03d.txt", id);
    char* outpar = concat("Chains/", name);
    outpar       = concat(outpar,outname);
    output       = fopen(outpar,"w");
		printf("%s\n",outpar);
    FILE *arx;
    char* parameters;
    char filename[sizeof "_sampler100.txt"];

    float pc,pc_b;
  	double n1  , n2 , fc  , dv, prob  ;
    double n1_b, n2_b, fc_b, dv_b, prob_b, accept_rate;
	  double alpha, rnd, numerator, denominator;
  	int i,counter,counter_b, ntrial, nacc, nsamp;
    int i1,i2,i3;
    int l1,l2,l3;
    int u1,u2,u3;
    int number_file;
		nsamp = steps/n_procs;
    double aux;
    double* tab_1 = malloc(Nbins*sizeof(double*));
    double* tab_2 = malloc(Nbins*sizeof(double*));
    double* tab_3 = malloc(Nbins*sizeof(double*));
    double* tab_4 = malloc(Nbins*sizeof(double*));
    double* tab_5 = malloc(Nbins*sizeof(double*));
    double* tab_6 = malloc(Nbins*sizeof(double*));
    double* tab_7 = malloc(Nbins*sizeof(double*));
    double* tab_8 = malloc(Nbins*sizeof(double*));
    double* temp  = malloc(Nbins*sizeof(double*));


    double mu1,mu2,mu3,mu4;
    double sd1,sd2,sd3,sd4;
    double n11,n21,pc1;
    double n12,n22,pc2;
    double n13,n23,pc3;
    double n14,n24,pc4;
    double n15,n25,pc5;
    double n16,n26,pc6;
    double n17,n27,pc7;
    double n18,n28,pc8;

    float v1,v2,v3,v4,v5,v6,v7,v8,vtot;

    char line[700];
    int iline;

	/***********************
	 * First random numbers
	************************/

    n1=uniform(n1_min,n1_max);
    n2=uniform(n2_min,n2_max);
    pc=uniform(pc_min,pc_max);
    dv=uniform(dv_min,dv_max);
    counter=0;
    ntrial=0;
    nacc=0;

    mu1=0.0*(n1_min+n1_max);
    mu2=0.0*(n2_min+n2_max);
    mu3=0.0*(pc_min+pc_max);
    mu4=0.0*(dv_min+dv_max);
    sd1=0.00*(n1_max-n1_min);
    sd2=0.00*(n2_max-n2_min);
    sd3=0.00*(pc_max-pc_min);
    sd4=0.00*(dv_max-dv_min);
    counter_b=-100;
    while(counter<nchain){

        ntrial++;
        if(counter<200){
            if(counter%4==0)n1=uniform(n1_min,n1_max);
            if(counter%4==1)n2=uniform(n2_min,n2_max);
            if(counter%4==2)pc=uniform(pc_min,pc_max);
            if(counter%4==3)dv=uniform(dv_min,dv_max);}

        else{
            do{n1=n1_b+0.6*sqrt(-2.0*sd1*log(uniform(0,1)))*cos(6.2831853*uniform(0,1));}while((n1<n1_min)||(n1>n1_max));
            do{n2=n2_b+0.6*sqrt(-2.0*sd2*log(uniform(0,1)))*cos(6.2831853*uniform(0,1));}while((n2<n2_min)||(n2>n2_max));
            do{pc=pc_b+0.6*sqrt(-2.0*sd3*log(uniform(0,1)))*cos(6.28318530718*uniform(0,1));}while((pc<pc_min)||(pc>pc_max));
            do{dv=dv_b+0.6*sqrt(-2.0*sd4*log(uniform(0,1)))*cos(6.2831853*uniform(0,1));}while((dv<dv_min)||(dv>dv_max));
            }
				//printf("%f\t %f\t %f\t %f\n",n1,n2,pc,dv);
				//printf("lims %f\t %f\t %f\t %f\n",pc_min,pc_max,dv_min,dv_max);
        aux= 1.0*(steps-1)*1.0*(n1-n1_min)/(n1_max-n1_min);

        l1 = (int) aux;
        u1 = l1+1;

        number_file = l1/nsamp;

				for(i=0;i<nsamp;i++){if(l1%nsamp==i) i1 = i*steps*steps;}

        aux = (steps-1)*1.0*(n2-n2_min)/(n2_max-n2_min);

        l2 = (int) aux;
        u2 = l2+1;

        i2 = l2*steps;
        aux = (steps-1)*1.0*(log(1.0*pc/pc_min)/log(1.0*pc_max/pc_min)) +1;

        l3 = (int) aux;
        i3 = l3;

				//printf("%d\t %d\t %d\t %d\t %d\t %d\n",number_file,l1,nsamp,i1,i2,i3);
        sprintf(filename, "_sampler%03d.txt", number_file);
        parameters = concat("Samplers/", snapshot);
        parameters       = concat(parameters,filename);
        arx = fopen(parameters,"r");

        iline = 0;
        while (fgets(line, sizeof(line), arx)) {
            iline++;
            //printf("ints %d\t %d\n",iline,i1+i2+i3);

            //// PRIMER NODO ////

            if(iline == i1+i2+i3){
//READSCALES
						sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
										&n11,&n21,&pc1
 										, &tab_1[0], &tab_1[1], &tab_1[2], &tab_1[3], &tab_1[4], &tab_1[5], &tab_1[6], &tab_1[7], &tab_1[8], &tab_1[9]
										, &tab_1[10], &tab_1[11], &tab_1[12], &tab_1[13], &tab_1[14], &tab_1[15], &tab_1[16], &tab_1[17], &tab_1[18], &tab_1[19]
										, &tab_1[20], &tab_1[21], &tab_1[22], &tab_1[23], &tab_1[24], &tab_1[25], &tab_1[26], &tab_1[27], &tab_1[28], &tab_1[29]
										, &tab_1[30], &tab_1[31], &tab_1[32], &tab_1[33], &tab_1[34], &tab_1[35], &tab_1[36], &tab_1[37], &tab_1[38], &tab_1[39]
										, &tab_1[40], &tab_1[41], &tab_1[42], &tab_1[43], &tab_1[44], &tab_1[45], &tab_1[46], &tab_1[47], &tab_1[48], &tab_1[49]);}
//SCALES
            //// SEGUNDO NODO ////

            if(iline == i1+i2+i3+1){
//READSCALES
						sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
										&n12,&n22,&pc2
 										, &tab_2[0], &tab_2[1], &tab_2[2], &tab_2[3], &tab_2[4], &tab_2[5], &tab_2[6], &tab_2[7], &tab_2[8], &tab_2[9]
										, &tab_2[10], &tab_2[11], &tab_2[12], &tab_2[13], &tab_2[14], &tab_2[15], &tab_2[16], &tab_2[17], &tab_2[18], &tab_2[19]
										, &tab_2[20], &tab_2[21], &tab_2[22], &tab_2[23], &tab_2[24], &tab_2[25], &tab_2[26], &tab_2[27], &tab_2[28], &tab_2[29]
										, &tab_2[30], &tab_2[31], &tab_2[32], &tab_2[33], &tab_2[34], &tab_2[35], &tab_2[36], &tab_2[37], &tab_2[38], &tab_2[39]
										, &tab_2[40], &tab_2[41], &tab_2[42], &tab_2[43], &tab_2[44], &tab_2[45], &tab_2[46], &tab_2[47], &tab_2[48], &tab_2[49]);}
//SCALES
            //// TERCER NODO ////

            if(iline == i1+i2+i3+steps){
//READSCALES
						sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
										&n13,&n23,&pc3
 										, &tab_3[0], &tab_3[1], &tab_3[2], &tab_3[3], &tab_3[4], &tab_3[5], &tab_3[6], &tab_3[7], &tab_3[8], &tab_3[9]
										, &tab_3[10], &tab_3[11], &tab_3[12], &tab_3[13], &tab_3[14], &tab_3[15], &tab_3[16], &tab_3[17], &tab_3[18], &tab_3[19]
										, &tab_3[20], &tab_3[21], &tab_3[22], &tab_3[23], &tab_3[24], &tab_3[25], &tab_3[26], &tab_3[27], &tab_3[28], &tab_3[29]
										, &tab_3[30], &tab_3[31], &tab_3[32], &tab_3[33], &tab_3[34], &tab_3[35], &tab_3[36], &tab_3[37], &tab_3[38], &tab_3[39]
										, &tab_3[40], &tab_3[41], &tab_3[42], &tab_3[43], &tab_3[44], &tab_3[45], &tab_3[46], &tab_3[47], &tab_3[48], &tab_3[49]);}
//SCALES
            //// CUARTO NODO ////

            if(iline == i1+i2+i3+steps+1){
//READSCALES
						sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
										&n14,&n24,&pc4
 										, &tab_4[0], &tab_4[1], &tab_4[2], &tab_4[3], &tab_4[4], &tab_4[5], &tab_4[6], &tab_4[7], &tab_4[8], &tab_4[9]
										, &tab_4[10], &tab_4[11], &tab_4[12], &tab_4[13], &tab_4[14], &tab_4[15], &tab_4[16], &tab_4[17], &tab_4[18], &tab_4[19]
										, &tab_4[20], &tab_4[21], &tab_4[22], &tab_4[23], &tab_4[24], &tab_4[25], &tab_4[26], &tab_4[27], &tab_4[28], &tab_4[29]
										, &tab_4[30], &tab_4[31], &tab_4[32], &tab_4[33], &tab_4[34], &tab_4[35], &tab_4[36], &tab_4[37], &tab_4[38], &tab_4[39]
										, &tab_4[40], &tab_4[41], &tab_4[42], &tab_4[43], &tab_4[44], &tab_4[45], &tab_4[46], &tab_4[47], &tab_4[48], &tab_4[49]);}
//SCALES

            if(iline == i1+i2+i3+steps+1)break;

            }

        fclose(arx);

        number_file = (l1+1)/nsamp;

				for(i=0;i<nsamp;i++){if((l1+1)%nsamp==i) i1 = i*steps*steps;}

        sprintf(filename, "_sampler%03d.txt", number_file);
        parameters = concat("Samplers/", snapshot);
        parameters = concat(parameters,filename);
        arx = fopen(parameters,"r");
        iline = 0;
				//printf("%d\t %d\t %d\t %d\t %d\t %d\n",number_file,l1,nsamp,i1,i2,i3);
        while (fgets(line, sizeof(line), arx)) {
            iline++;
            //printf("ints %d\t %d\n",iline,i1+i2+i3);

            //// PRIMER NODO ////

            if(iline == i1+i2+i3){
//READSCALES
						sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
										&n15,&n25,&pc5
 										, &tab_5[0], &tab_5[1], &tab_5[2], &tab_5[3], &tab_5[4], &tab_5[5], &tab_5[6], &tab_5[7], &tab_5[8], &tab_5[9]
										, &tab_5[10], &tab_5[11], &tab_5[12], &tab_5[13], &tab_5[14], &tab_5[15], &tab_5[16], &tab_5[17], &tab_5[18], &tab_5[19]
										, &tab_5[20], &tab_5[21], &tab_5[22], &tab_5[23], &tab_5[24], &tab_5[25], &tab_5[26], &tab_5[27], &tab_5[28], &tab_5[29]
										, &tab_5[30], &tab_5[31], &tab_5[32], &tab_5[33], &tab_5[34], &tab_5[35], &tab_5[36], &tab_5[37], &tab_5[38], &tab_5[39]
										, &tab_5[40], &tab_5[41], &tab_5[42], &tab_5[43], &tab_5[44], &tab_5[45], &tab_5[46], &tab_5[47], &tab_5[48], &tab_5[49]);}
//SCALES
            //// SEGUNDO NODO ////

            if(iline == i1+i2+i3+1){
//READSCALES
						sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
										&n16,&n26,&pc6
 										, &tab_6[0], &tab_6[1], &tab_6[2], &tab_6[3], &tab_6[4], &tab_6[5], &tab_6[6], &tab_6[7], &tab_6[8], &tab_6[9]
										, &tab_6[10], &tab_6[11], &tab_6[12], &tab_6[13], &tab_6[14], &tab_6[15], &tab_6[16], &tab_6[17], &tab_6[18], &tab_6[19]
										, &tab_6[20], &tab_6[21], &tab_6[22], &tab_6[23], &tab_6[24], &tab_6[25], &tab_6[26], &tab_6[27], &tab_6[28], &tab_6[29]
										, &tab_6[30], &tab_6[31], &tab_6[32], &tab_6[33], &tab_6[34], &tab_6[35], &tab_6[36], &tab_6[37], &tab_6[38], &tab_6[39]
										, &tab_6[40], &tab_6[41], &tab_6[42], &tab_6[43], &tab_6[44], &tab_6[45], &tab_6[46], &tab_6[47], &tab_6[48], &tab_6[49]);}
//SCALES
            //// TERCER NODO ////

            if(iline == i1+i2+i3+steps){
//READSCALES
						sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
										&n17,&n27,&pc7
 										, &tab_7[0], &tab_7[1], &tab_7[2], &tab_7[3], &tab_7[4], &tab_7[5], &tab_7[6], &tab_7[7], &tab_7[8], &tab_7[9]
										, &tab_7[10], &tab_7[11], &tab_7[12], &tab_7[13], &tab_7[14], &tab_7[15], &tab_7[16], &tab_7[17], &tab_7[18], &tab_7[19]
										, &tab_7[20], &tab_7[21], &tab_7[22], &tab_7[23], &tab_7[24], &tab_7[25], &tab_7[26], &tab_7[27], &tab_7[28], &tab_7[29]
										, &tab_7[30], &tab_7[31], &tab_7[32], &tab_7[33], &tab_7[34], &tab_7[35], &tab_7[36], &tab_7[37], &tab_7[38], &tab_7[39]
										, &tab_7[40], &tab_7[41], &tab_7[42], &tab_7[43], &tab_7[44], &tab_7[45], &tab_7[46], &tab_7[47], &tab_7[48], &tab_7[49]);}
//SCALES
            //// CUARTO NODO ////

            if(iline == i1+i2+i3+steps+1){
//READSCALES
						sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
										&n18,&n28,&pc8
 										, &tab_8[0], &tab_8[1], &tab_8[2], &tab_8[3], &tab_8[4], &tab_8[5], &tab_8[6], &tab_8[7], &tab_8[8], &tab_8[9]
										, &tab_8[10], &tab_8[11], &tab_8[12], &tab_8[13], &tab_8[14], &tab_8[15], &tab_8[16], &tab_8[17], &tab_8[18], &tab_8[19]
										, &tab_8[20], &tab_8[21], &tab_8[22], &tab_8[23], &tab_8[24], &tab_8[25], &tab_8[26], &tab_8[27], &tab_8[28], &tab_8[29]
										, &tab_8[30], &tab_8[31], &tab_8[32], &tab_8[33], &tab_8[34], &tab_8[35], &tab_8[36], &tab_8[37], &tab_8[38], &tab_8[39]
										, &tab_8[40], &tab_8[41], &tab_8[42], &tab_8[43], &tab_8[44], &tab_8[45], &tab_8[46], &tab_8[47], &tab_8[48], &tab_8[49]);}
//SCALES
            /////////////////////

            if(iline == i1+i2+i3+steps+1)break;
            }

        fclose(arx);
        v1=1.0/fabs((n1-n11)*(n2-n21)*(pc-pc1));
        v2=1.0/fabs((n1-n12)*(n2-n22)*(pc-pc2));
        v3=1.0/fabs((n1-n13)*(n2-n23)*(pc-pc3));
        v4=1.0/fabs((n1-n14)*(n2-n24)*(pc-pc4));
        v5=1.0/fabs((n1-n15)*(n2-n25)*(pc-pc5));
        v6=1.0/fabs((n1-n16)*(n2-n26)*(pc-pc6));
        v7=1.0/fabs((n1-n17)*(n2-n27)*(pc-pc7));
        v8=1.0/fabs((n1-n18)*(n2-n28)*(pc-pc8));
        vtot=v1+v2+v3+v4+v5+v6+v7+v8;
        v1/=vtot;
        v2/=vtot;
        v3/=vtot;
        v4/=vtot;
        v5/=vtot;
        v6/=vtot;
        v7/=vtot;
        v8/=vtot;

				/*
				printf("0th %f\t %f\t %f\n",n1,n2,pc);
				printf("1th %f\t %f\t %f\n",n11,n21,pc1);
				printf("2th %f\t %f\t %f\n",n12,n22,pc2);
				printf("3th %f\t %f\t %f\n",n13,n23,pc3);
				printf("4th %f\t %f\t %f\n",n14,n24,pc4);
				printf("5th %f\t %f\t %f\n",n15,n25,pc5);
				printf("6th %f\t %f\t %f\n",n16,n26,pc6);
				printf("7th %f\t %f\t %f\n",n17,n27,pc7);
				printf("8th %f\t %f\t %f\n",n18,n28,pc8);
				*/
        for(i=0;i<Nbins;i++)temp[i]=v1*tab_1[i]+v2*tab_2[i]+v3*tab_3[i]+v4*tab_4[i]+v5*tab_5[i]+v6*tab_6[i]+v7*tab_7[i]+v8*tab_8[i];

        numerator=0.0;
        denominator=0.0;

        for(i=0;i<Nbins;i++){
            numerator+=sigma[i]*sigma[i]/(error[i]*error[i]);
            denominator+=temp[i]*sigma[i]/(error[i]*error[i]);}
        fc=numerator/denominator;
        fc+=dv;

        prob=0.0;

        for(i=0;i<Nbins;i++){
            prob=prob-0.5*(pow((sigma[i]-fc*temp[i])/error[i],2)+2*log(error[i]));}

        if(counter==0)prob_b=prob-1.0e5;
        alpha=prob-prob_b;

        if (alpha>0.0){
            prob_b=prob;
            n1_b=n1;
            n2_b=n2;
            pc_b=pc;
            fc_b=fc;
            dv_b=dv;
            if(counter_b>=0){
            sd1 = (1.0*counter_b/(1.0*counter_b+1))*(sd1 +(n1_b-mu1)*(n1_b-mu1)/(1.0*counter_b+1));
            sd2 = (1.0*counter_b/(1.0*counter_b+1))*(sd2 +(n2_b-mu2)*(n2_b-mu2)/(1.0*counter_b+1));
            sd3 = (1.0*counter_b/(1.0*counter_b+1))*(sd3 +(pc_b-mu3)*(pc_b-mu3)/(1.0*counter_b+1));
            sd4 = (1.0*counter_b/(1.0*counter_b+1))*(sd4 +(dv_b-mu4)*(dv_b-mu4)/(1.0*counter_b+1));

            mu1 = (1.0*counter_b*mu1+n1_b)/(1.0*counter_b+1);
            mu2 = (1.0*counter_b*mu2+n2_b)/(1.0*counter_b+1);
            mu3 = (1.0*counter_b*mu3+pc_b)/(1.0*counter_b+1);
            mu4 = (1.0*counter_b*mu4+dv_b)/(1.0*counter_b+1);}
            counter_b++;
            counter++;
            nacc++;

            if(counter>1000)fprintf(output, "%f\t %f\t %10.8f\t %f\t %f\n", n1, n2, pc, pcyrtokms*fc, accept_rate);}
        else{
            rnd=log(uniform(0.0,1.0));
            if(rnd<alpha){
                prob_b=prob;
                n1_b=n1;
                n2_b=n2;
                pc_b=pc;
                fc_b=fc;
                dv_b=dv;

                if(counter_b>=0){
                sd1 = (1.0*counter_b/(1.0*counter_b+1))*(sd1 +(n1_b-mu1)*(n1_b-mu1)/(1.0*counter_b+1));
                sd2 = (1.0*counter_b/(1.0*counter_b+1))*(sd2 +(n2_b-mu2)*(n2_b-mu2)/(1.0*counter_b+1));
                sd3 = (1.0*counter_b/(1.0*counter_b+1))*(sd3 +(pc_b-mu3)*(pc_b-mu3)/(1.0*counter_b+1));
                sd4 = (1.0*counter_b/(1.0*counter_b+1))*(sd4 +(dv_b-mu4)*(dv_b-mu4)/(1.0*counter_b+1));

                mu1 = (1.0*counter_b*mu1+n1_b)/(1.0*counter_b+1);
                mu2 = (1.0*counter_b*mu2+n2_b)/(1.0*counter_b+1);
                mu3 = (1.0*counter_b*mu3+pc_b)/(1.0*counter_b+1);
                mu4 = (1.0*counter_b*mu4+dv_b)/(1.0*counter_b+1);}

                counter++;
                counter_b++;
                nacc++;
            if(counter>1000)fprintf(output, "%f\t %f\t %10.8f\t %f\t %f\n", n1, n2, pc, pcyrtokms*fc, accept_rate);

            }}

        accept_rate=1.0*nacc/(1.0*ntrial);
        //if(counter>200){
        //if(id==0)printf("acceptance rate %f\n",accept_rate);
        //if(id==0)printf("%f\t %f\t %.15f\t %f\t %f\t %d\n",sd1,sd2,sd3,sd4, pc_b,counter_b);}



    }

		fclose(output);

	//return 0.0;

}


main(int argc, char* argv[]) {
	MPI_Status status;
	int my_id, an_id, root_process, ierr, i, num_procs,rec_id, fin;
	root_process=0;
	fin=0;
	int N, Nres, Nbins, nchain, rmin,rmax;
	double L, out;
	double pmin,pmax,dx,n1_min,n1_max,n2_min,n2_max,pc_min,pc_max,dv_min,dv_max;
	char name[64];
	char line[2000];
	nchain=10000;
	sscanf(argv[1],"%d",&rmin);
	sscanf(argv[2],"%d",&rmax);

	//clock_t start, end;
	//double cpu_time_used;

	FILE *arx;
	arx = fopen("Input.txt","r");
	while (fgets(line, sizeof(line), arx)){
		sscanf(line,"%*s %63s %lf %d %*lf %d %lf %lf %lf %lf %lf %lf %*d",
		&name,&L,&N,&Nbins,&n1_min,&n1_max,&n2_min,&n2_max,&dv_min,&dv_max);
		break;}
	fclose(arx);


	char filename[sizeof "_00000_00000"];
	sprintf(filename, "_%05d_%05d", rmin, rmax);
	char* parameters = concat("Files/", name);
	parameters       = concat(parameters,filename);
	parameters       = concat(parameters,".txt");
	char* snapshot = concat(name,filename);
	arx = fopen(parameters,"r");

	float* delta = malloc(Nbins*sizeof(float*));
	float* sigma = malloc(Nbins*sizeof(float*));
	float* error = malloc(Nbins*sizeof(float*));

	while (fgets(line, sizeof(line), arx)){
	//READARRAYS
	sscanf(line," %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", &delta[0], &delta[1], &delta[2], &delta[3], &delta[4], &delta[5], &delta[6], &delta[7], &delta[8], &delta[9], &delta[10], &delta[11], &delta[12], &delta[13], &delta[14], &delta[15], &delta[16], &delta[17], &delta[18], &delta[19], &delta[20], &delta[21], &delta[22], &delta[23], &delta[24], &delta[25], &delta[26], &delta[27], &delta[28], &delta[29], &delta[30], &delta[31], &delta[32], &delta[33], &delta[34], &delta[35], &delta[36], &delta[37], &delta[38], &delta[39], &delta[40], &delta[41], &delta[42], &delta[43], &delta[44], &delta[45], &delta[46], &delta[47], &delta[48], &delta[49], &sigma[0], &sigma[1], &sigma[2], &sigma[3], &sigma[4], &sigma[5], &sigma[6], &sigma[7], &sigma[8], &sigma[9], &sigma[10], &sigma[11], &sigma[12], &sigma[13], &sigma[14], &sigma[15], &sigma[16], &sigma[17], &sigma[18], &sigma[19], &sigma[20], &sigma[21], &sigma[22], &sigma[23], &sigma[24], &sigma[25], &sigma[26], &sigma[27], &sigma[28], &sigma[29], &sigma[30], &sigma[31], &sigma[32], &sigma[33], &sigma[34], &sigma[35], &sigma[36], &sigma[37], &sigma[38], &sigma[39], &sigma[40], &sigma[41], &sigma[42], &sigma[43], &sigma[44], &sigma[45], &sigma[46], &sigma[47], &sigma[48], &sigma[49], &error[0], &error[1], &error[2], &error[3], &error[4], &error[5], &error[6], &error[7], &error[8], &error[9], &error[10], &error[11], &error[12], &error[13], &error[14], &error[15], &error[16], &error[17], &error[18], &error[19], &error[20], &error[21], &error[22], &error[23], &error[24], &error[25], &error[26], &error[27], &error[28], &error[29], &error[30], &error[31], &error[32], &error[33], &error[34], &error[35], &error[36], &error[37], &error[38], &error[39], &error[40], &error[41], &error[42], &error[43], &error[44], &error[45], &error[46], &error[47], &error[48], &error[49]);
	//ARRAYS
	break;}
	fclose(arx);


pmin=4.0/L;
dx=L/(1.0*N);
pmax=1.0/(4.0*dx);

pc_min=2*pmin;
pc_max=0.5*pmax;

ierr = MPI_Init(&argc, &argv);

ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);    // Get id of the current processor
ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);// Get number of processes
srand(time(NULL)+ my_id);

if(my_id==root_process){

	// Send respective id
	for(an_id=1 ; an_id<num_procs ; an_id++){
    		ierr = MPI_Send(&an_id, 1, MPI_INT, an_id, send_data_tag, MPI_COMM_WORLD);
	}

	out=walker(snapshot,name,n1_min , n1_max, n2_min, n2_max, pc_min, pc_max, dv_min, dv_max,
	  	pmin, pmax, Nbins, dx, sigma, delta, error, nchain,num_procs, 0);
	for(an_id=1 ; an_id<num_procs ; an_id++)
	ierr = MPI_Recv(&fin, 1, MPI_INT, MPI_ANY_SOURCE, finish_data_tag, MPI_COMM_WORLD, &status);
}

else{
    	ierr = MPI_Recv(&rec_id, 1, MPI_INT, 0, send_data_tag, MPI_COMM_WORLD, &status);
	out=walker(snapshot,name,n1_min , n1_max, n2_min, n2_max, pc_min, pc_max, dv_min, dv_max,
                pmin, pmax, Nbins, dx, sigma, delta, error, nchain,num_procs, rec_id);

	ierr = MPI_Send(&fin, 1, MPI_INT, 0, finish_data_tag, MPI_COMM_WORLD);
	}
ierr=MPI_Finalize();
return 0;
}
