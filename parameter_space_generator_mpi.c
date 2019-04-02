#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include <string.h>
#include <mpi.h>

#define PI 3.14159265
#define tiny 1e-16
#define fmin 0.1
#define fmax 100.0
#define send_data_tag 2001
#define finish_data_tag 2002
#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

double** two_column_file(char *filename, int Nrows){
  int i;
  int j;
  double** mat=malloc(Nrows*sizeof(double*));
  for(i=0;i<Nrows;++i)
  mat[i]=malloc(2*sizeof(double));
  FILE *file;
  file=fopen(filename, "r");

  for(i = 0; i < Nrows; i++){
    for(j = 0; j < 2; j++) {
       if (!fscanf(file, "%lf", &mat[i][j]))
          break;}}
  fclose(file);
  return mat;}

char* concat(const char *s1, const char *s2)
  {
      char *result = malloc(strlen(s1)+strlen(s2)+1);//+1 for the null-terminator
      //in real code you would check for errors in malloc here
      strcpy(result, s1);
      strcat(result, s2);
      return result;
  }

double sin_delta(double p, double q,double delta){
  double result;
  result=sin(PI*p*delta)*sin(PI*q*delta)/(delta*delta*PI*PI);
  result=result/((p+tiny)*(q+tiny));
  result=result*result;
  return result;
}

double gama(double p, double n1, double n2, double pc, double pmin, double pmax){
  if (p<pmin) return 0;
  if (p>pmax) return 0;
  if (p<pc) return 1.0/pow(p,n1);
  else return pow(pc,n2-n1)/pow(p,n2);
}

double func1(double r, double n1, double n2, double pmin, double pc ,double f1){
  double result;
  if (r<f1){
    result=pow(pow(pmin,2.0-2.0*n1)+(2.0-2.0*n1)*r,1.0/(2.0-2.0*n1));}
  else{
    result=pow(pc,2.0-2.0*n2)+(2.0-2.0*n2)*(r-f1)*pow(pc,2.0*(n1-n2));
    result=pow(result,1.0/(2.0-2.0*n2));}
  return result;
}

double func2(double r, double n1, double n2, double pmin, double pc ,double f1){
  double result;
  if (r<f1){
    result=pmin*exp(r);}
  else{
  result=pow(pc,2.0-2.0*n2)+(2.0-2.0*n2)*(r-f1)*pow(pc,2.0*(n1-n2));
  result=pow(result,1.0/(2.0-2.0*n2));}
  return result;
}

double func3(double r, double n1, double n2, double pmin, double pc ,double f1){
  double result;
  if (r<f1){
    result=pow(pow(pmin,2.0-2.0*n1)+(2.0-2.0*n1)*r,1.0/(2.0-2.0*n1));}
  else{
    result=pc*exp((r-f1)*pow(pc,2*(n1-n2)));}
  return result;}

double func4(double r, double n1, double n2, double pmin, double pc ,double f1){
  double result;
  if (r<f1){
    result=pmin*exp(r);}
  else{
    result=pc*exp((r-f1)*pow(pc,2*(n1-n2)));}
  return result;}

float *sigmas(double n1, double n2, double pc, double pmin, double pmax, int Nbins, float *delta, double dx){
  n1=n1-1;
  n2=n2-1;
  //srand(time(0));
  int Npoints, i, j;
  float random,angle;
  double f1,f2,g1,g2,A,B,p,q,sum;
  double maximo=1.0*RAND_MAX;
  float *result = malloc(sizeof(float) * Nbins);

  if (n1!=1.0 && n2!=1.0){
    f1=pow(pc,2.0-2.0*n1)/(2.0-2.0*n1) - pow(pmin,2.0-2.0*n1)/(2.0-2.0*n1);
    f2=pow(pmax,2.0-2.0*n2)-pow(pc,2.0-2.0*n2);
    f2=f2*pow(pc,2.0*(n2-n1))/(2.0-2.0*n2);

    //These are the condition of the index for the velocity field
    if (n1!=0.0 && n2!=0.0){
      g1=pow(pc,-2.0*n1)/(-2.0*n1) - pow(pmin,-2.0*n1)/(-2.0*n1);
      g2=pow(pmax,-2.0*n2)-pow(pc,-2.0*n2);
      g2=g2*pow(pc,2.0*(n2-n1))/(-2.0*n2);}
    else if(n1!=0.0){
      g1=pow(pc,-2.0*n1)/(-2.0*n1) - pow(pmin,-2.0*n1)/(-2.0*n1);
      g2=log(pmax/pc)*pow(pc,2.0*(n2-n1));}
    else if(n2!=0){
      g1=log(pc/pmin);
      g2=pow(pmax,-2.0*n2)-pow(pc,-2.0*n2);
      g2=g2*pow(pc,2.0*(n2-n1))/(-2.0*n2);}
    else {
      g1=log(pc/pmin);
      g2=log(pmax/pc)*pow(pc,2.0*(n2-n1));}

    A=f1+f2;
    B=g1+g2;


    for(j = 0; j < Nbins; j++){
      sum=0.0;
      Npoints= 10+(int) 15000*(pmax-pmin)*delta[j];
      for(i = 0; i < Npoints; i++){
        random = rand()/maximo;
        angle  = 0.5*PI*rand()/maximo;
        random = func1(random*A,n1,n2,pmin,pc,f1);
        p=random*cos(angle);
        q=random*sin(angle);
        sum=sum+sin_delta(p,q,delta[j]);}
      result[j]=2*PI*sqrt(A*sum/(1.0*B*Npoints))*dx*dx;}
  }
  else if (n1==1.0 && n2!=1.0){
    f1=log(pc/pmin);
    f2=pow(pmax,2.0-2.0*n2)-pow(pc,2.0-2.0*n2);
    f2=f2*pow(pc,2.0*(n2-n1))/(2.0-2.0*n2);
    g1=pow(pc,-2.0*n1)/(-2.0*n1) - pow(pmin,-2.0*n1)/(-2.0*n1);

    if (n2!=0.0){
      g2=pow(pmax,-2.0*n2)-pow(pc,-2.0*n2);
      g2=g2*pow(pc,2.0*(n2-n1))/(-2.0*n2);}
    else {g2=log(pmax/pc)*pow(pc,2.0*(n2-n1));}

    A=f1+f2;
    B=g1+g2;

    for(j = 0; j < Nbins; j++){
      sum=0.0;
      Npoints= 10+(int) 15000*(pmax-pmin)*delta[j];
      for(i = 0; i < Npoints; i++){
        random = rand()/maximo;
        angle  = 0.5*PI*rand()/maximo;
        random = func1(random*A,n1,n2,pmin,pc,f1);
        p=random*cos(angle);
        q=random*sin(angle);
        sum=sum+sin_delta(p,q,delta[j]);}
      result[j]=2*PI*sqrt(A*sum/(1.0*B*Npoints))*dx*dx;}

  }
  else if (n1!=1.0 && n2==1.0){
    f1=pow(pc,2.0-2.0*n1)/(2.0-2.0*n1) - pow(pmin,2.0-2.0*n1)/(2.0-2.0*n1);
    f2=log(pmax/pc)*pow(pc,2.0*(n2-n1));
    g2=pow(pmax,-2.0*n2)-pow(pc,-2.0*n2);
    g2=g2*pow(pc,2.0*(n2-n1))/(-2.0*n2);

    if (n1!=0.0){g1=pow(pc,-2.0*n1)/(-2.0*n1) - pow(pmin,-2.0*n1)/(-2.0*n1);}
    else {g1=log(pc/pmin);}

    A=f1+f2;
    B=g1+g2;

    for(j = 0; j < Nbins; j++){
      sum=0.0;
      Npoints= 10+(int) 15000*(pmax-pmin)*delta[j];
      for(i = 0; i < Npoints; i++){
        random = rand()/maximo;
        angle  = 0.5*PI*rand()/maximo;
        random = func1(random*A,n1,n2,pmin,pc,f1);
        p=random*cos(angle);
        q=random*sin(angle);
        sum=sum+sin_delta(p,q,delta[j]);}
      result[j]=2*PI*sqrt(A*sum/(1.0*B*Npoints))*dx*dx;}
  }
  else if (n1==1.0 && n2==1.0){
    f1=log(pc/pmin);
    f2=log(pmax/pc)*pow(pc,2.0*(n2-n1));
    g1=pow(pc,-2.0*n1)/(-2.0*n1) - pow(pmin,-2.0*n1)/(-2.0*n1);
    g2=pow(pmax,-2.0*n2)-pow(pc,-2.0*n2);
    g2=g2*pow(pc,2.0*(n2-n1))/(-2.0*n2);

    A=f1+f2;
    B=g1+g2;

    for( j = 0; j < Nbins; j++){
      sum=0.0;
      Npoints= 10+(int) 15000*(pmax-pmin)*delta[j];
      for( i = 0; i < Npoints; i++){
        random = rand()/maximo;
        angle  = 0.5*PI*rand()/maximo;
        random = func1(random*A,n1,n2,pmin,pc,f1);
        p=random*cos(angle);
        q=random*sin(angle);
        sum=sum+sin_delta(p,q,delta[j]);}
      result[j]=2*PI*sqrt(A*sum/(1.0*B*Npoints))*dx*dx;}}
  return result;
  }

double one_sigma(double n1, double n2, double pc, double pmin, double pmax, double dx){
    n1=n1-1;
    n2=n2-1;
    //srand(time(0));
    double f1,f2,g1,g2,A,B,p,q,sum;
    double result;

    if (n1!=1.0 && n2!=1.0){
      f1=pow(pc,2.0-2.0*n1)/(2.0-2.0*n1) - pow(pmin,2.0-2.0*n1)/(2.0-2.0*n1);
      f2=pow(pmax,2.0-2.0*n2)-pow(pc,2.0-2.0*n2);
      f2=f2*pow(pc,2.0*(n2-n1))/(2.0-2.0*n2);

      //These are the condition of the index for the velocity field
      if (n1!=0.0 && n2!=0.0){
        g1=pow(pc,-2.0*n1)/(-2.0*n1) - pow(pmin,-2.0*n1)/(-2.0*n1);
        g2=pow(pmax,-2.0*n2)-pow(pc,-2.0*n2);
        g2=g2*pow(pc,2.0*(n2-n1))/(-2.0*n2);}
      else if(n1!=0.0){
        g1=pow(pc,-2.0*n1)/(-2.0*n1) - pow(pmin,-2.0*n1)/(-2.0*n1);
        g2=log(pmax/pc)*pow(pc,2.0*(n2-n1));}
      else if(n2!=0){
        g1=log(pc/pmin);
        g2=pow(pmax,-2.0*n2)-pow(pc,-2.0*n2);
        g2=g2*pow(pc,2.0*(n2-n1))/(-2.0*n2);}
      else {
        g1=log(pc/pmin);
        g2=log(pmax/pc)*pow(pc,2.0*(n2-n1));}

      A=f1+f2;
      B=g1+g2;

      result=2*PI*sqrt(A/(1.0*B));
    }
    else if (n1==1.0 && n2!=1.0){
      f1=log(pc/pmin);
      f2=pow(pmax,2.0-2.0*n2)-pow(pc,2.0-2.0*n2);
      f2=f2*pow(pc,2.0*(n2-n1))/(2.0-2.0*n2);
      g1=pow(pc,-2.0*n1)/(-2.0*n1) - pow(pmin,-2.0*n1)/(-2.0*n1);

      if (n2!=0.0){
        g2=pow(pmax,-2.0*n2)-pow(pc,-2.0*n2);
        g2=g2*pow(pc,2.0*(n2-n1))/(-2.0*n2);}
      else {g2=log(pmax/pc)*pow(pc,2.0*(n2-n1));}

      A=f1+f2;
      B=g1+g2;

      result=2*PI*sqrt(A/(1.0*B));

    }
    else if (n1!=1.0 && n2==1.0){
      f1=pow(pc,2.0-2.0*n1)/(2.0-2.0*n1) - pow(pmin,2.0-2.0*n1)/(2.0-2.0*n1);
      f2=log(pmax/pc)*pow(pc,2.0*(n2-n1));
      g2=pow(pmax,-2.0*n2)-pow(pc,-2.0*n2);
      g2=g2*pow(pc,2.0*(n2-n1))/(-2.0*n2);

      if (n1!=0.0){g1=pow(pc,-2.0*n1)/(-2.0*n1) - pow(pmin,-2.0*n1)/(-2.0*n1);}
      else {g1=log(pc/pmin);}

      A=f1+f2;
      B=g1+g2;

      result=2*PI*sqrt(A/(1.0*B));
    }
    else if (n1==1.0 && n2==1.0){
      f1=log(pc/pmin);
      f2=log(pmax/pc)*pow(pc,2.0*(n2-n1));
      g1=pow(pc,-2.0*n1)/(-2.0*n1) - pow(pmin,-2.0*n1)/(-2.0*n1);
      g2=pow(pmax,-2.0*n2)-pow(pc,-2.0*n2);
      g2=g2*pow(pc,2.0*(n2-n1))/(-2.0*n2);

      A=f1+f2;
      B=g1+g2;

      result=2*PI*sqrt(A/(1.0*B));}
    return result*dx*dx;
    }

double uniform(double minimo, double maximo){
  return minimo+(maximo-minimo)*rand()/(1.0*RAND_MAX);
}

double walker(char *name,double n1_min , double n1_max, double n2_min , double n2_max, double pc_min, double pc_max,
  double pmin, double pmax, int Nbins, double dx, float *delta, int steps, int id, int ncores){
        printf("This is the %03d\t sampler \n ", id);

    FILE *arx;
    char filename[sizeof "_sampler100.txt"];
    sprintf(filename, "_sampler%03d.txt", id);
    char* parameters = concat("Files_model_2/", name);
    parameters       = concat(parameters,filename);
    arx = fopen(parameters,"w");
	
  	double n1  , n2,  pc  , fc  , r_factor;
  	float *temp;
  	int i,counter;
    int i1,i2,i3;
    int low, upp; 

    int cond;
    int div;
    div  = (int)steps/ncores; 
    cond = steps - id - ncores*div;

    if (cond > 0){
        low = id*(div+1);
        upp = (id+1)*(div+1)-1;}
    else{
        low = steps - (ncores-id)*div;
        upp = steps -1 - (ncores -id -1)*div;}

	
    r_factor=pow(pc_max/pc_min,1.0/(1.0*steps-1));

	/***********************
	 * First random numbers 	
	************************/
    printf("steps %d\n",steps);	
    for(i1=low;i1<=upp;i1++){
        printf("indice , id %d\t %d\n",i1,id);
    for(i2=0;i2<steps;i2++){
    for(i3=0;i3<steps;i3++){    
        n1   = n1_min + (n1_max-n1_min)*i1/(1.0*steps-1);
        n2   = n2_min + (n2_max-n2_min)*i2/(1.0*steps-1);
  		pc   = pc_min*pow(r_factor,i3);
        temp = sigmas(n1, n2, pc, pmin, pmax, Nbins, delta,dx);
    	fprintf(arx, "%f\t %f\t %f", n1, n2, pc);
        for(i=0;i<Nbins;i++){ 
            fprintf(arx, "\t %f", temp[i]);
            }
        fprintf(arx,"\n");
    }}}	
	printf("end sampler %03d\n",id);
	
	fclose(arx);
	return 0.0;

}


main(int argc, char* argv[]) {
MPI_Status status;
int my_id, an_id, root_process, ierr, i, num_procs,rec_id, fin;
root_process=0;
fin=0;
int N, Nres, Nbins, steps;
double L, out;
double pmin,pmax,dx,n1_min,n1_max,n2_min,n2_max,pc_min,pc_max;
char *name;
//clock_t start, end;
//double cpu_time_used;


name=argv[1];
sscanf(argv[2],"%lf",&L);
sscanf(argv[3],"%d",&N);
sscanf(argv[4],"%lf",&n1_min);
sscanf(argv[5],"%lf",&n1_max);
sscanf(argv[6],"%lf",&n2_min);
sscanf(argv[7],"%lf",&n2_max);
sscanf(argv[8],"%lf",&pc_min);
sscanf(argv[9],"%lf",&pc_max);

Nbins=argc-10;


steps=60;
float delta[Nbins];
for(i = 0; i < Nbins; i++){
 sscanf(argv[10+i],"%f",&delta[i]);}

pmin=4.0/L;
pmax=(1.0*N)/(4.0*L);
dx=L/(1.0*N);

ierr = MPI_Init(&argc, &argv);

ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);    // Get id of the current processor
ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);// Get number of processes
srand(time(NULL)+ my_id);
//srand(0);
if(my_id==0)printf("Nbins %d\n",Nbins);
if(my_id==root_process){

	// Send respective id
	for(an_id=1 ; an_id<num_procs ; an_id++){
    		ierr = MPI_Send(&an_id, 1, MPI_INT, an_id, send_data_tag, MPI_COMM_WORLD);	
	}

	out=walker(name,n1_min , n1_max, n2_min, n2_max, pc_min, pc_max,
	  	pmin, pmax, Nbins, dx, delta, steps, 0, num_procs);
	for(an_id=1 ; an_id<num_procs ; an_id++)
	ierr = MPI_Recv(&fin, 1, MPI_INT, MPI_ANY_SOURCE, finish_data_tag, MPI_COMM_WORLD, &status);
}

else{
    	ierr = MPI_Recv(&rec_id, 1, MPI_INT, 0, send_data_tag, MPI_COMM_WORLD, &status);
	out=walker(name,n1_min , n1_max, n2_min, n2_max, pc_min, pc_max,
                pmin, pmax, Nbins, dx, delta, steps, rec_id, num_procs);
   
	
	ierr = MPI_Send(&fin, 1, MPI_INT, 0, finish_data_tag, MPI_COMM_WORLD);
	}
ierr=MPI_Finalize();
return 0;
}
