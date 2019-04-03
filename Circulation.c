#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include <string.h>
#include <mpi.h>

#define PI 3.14159265
#define tiny 1e-16
#define send_data_tag 2001
#define finish_data_tag 2002

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


main(int argc, char* argv[]) {
    MPI_Status status;
    int my_id, an_id, root_process, ierr, num_procs,rec_id, fin, Ngrids;

		ierr = MPI_Init(&argc, &argv);

		ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);    // Get id of the current processor
		ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

		int i, j, k, knew, counter, pix_counter;
		double area, vort, omeg, radius;
		int scale_min, scale_max, scale_i, scale_aux, step, is, js;
		int N_tent, N_scales;
		double jump;
		FILE *arx;
		char line[700];
    root_process=0;
    char name[64];
    double L, Rmax;

		arx = fopen("Input.txt","r");
		while (fgets(line, sizeof(line), arx)){
			sscanf(line,"%*s %63s %lf %d %lf %d %*lf %*lf %*lf %*lf %*lf %*lf",&name,&L,&Ngrids,&Rmax,&N_tent);
			break;}
		fclose(arx);
		printf("%s\t %f\t %d\t %f\t %d\n",name,L,Ngrids,Rmax,N_tent);

		N_scales=N_tent;
		int scales[N_scales];

		float **vort_map = (float **)malloc(Ngrids * sizeof(float*));
		for(i = 0; i < Ngrids; i++) vort_map[i] = (float *)malloc(Ngrids * sizeof(float));

		float **omeg_map = (float **)malloc(Ngrids * sizeof(float*));
		for(i = 0; i < Ngrids; i++) omeg_map[i] = (float *)malloc(Ngrids * sizeof(float));


		printf("Saving maps\n");
		///// Saving vorticity map in matrix ////

    char* map_file	= concat("Maps/", name);
    map_file        = concat(map_file,"_vort.txt");
    arx 						= fopen(map_file,"r");

		for(i = 0; i < Ngrids; i++){
    		for (j = 0 ; j < Ngrids; j++){
      			fscanf(arx,"%f",&vort_map[i][j]);
    		}}
  	fclose(arx);
		//////////////////////////////////////////

		///// Saving average map in matrix ////
		map_file	= concat("Maps/", name);
		map_file  = concat(map_file,"_av_vort.txt");
		arx 			= fopen(map_file,"r");

		for(i = 0; i < Ngrids; i++){
				for (j = 0 ; j < Ngrids; j++){
						fscanf(arx,"%f",&omeg_map[i][j]);
				}}
		fclose(arx);
		//////////////////////////////////////////

		scale_min=1;
		scale_max=(int) Ngrids*(0.5*L - Rmax)/L;

		//////////// define ideal spacing ////
		printf("%d\t %d\t %d\t %d\n",scale_min,scale_max,N_tent,Ngrids);
		do{
		counter=1;
		jump=pow(10,log10(1.0*scale_max)/N_tent);
		scale_aux=scale_min;
		for(i = 0; i < N_tent; i++){
					scale_i=scale_min*pow(jump,i);
					if (scale_aux!=scale_i){
						counter++;
						scale_aux=scale_i;}}
				N_tent++;
				if (N_tent>Ngrids)break;}while(counter<N_scales);

		N_tent--;
		counter=0;
		jump=pow(10,log10(1.0*scale_max)/N_tent);
		for(i = 0; i < N_tent; i++){
					scale_i=scale_min*pow(jump,i);
					if (scale_aux!=scale_i){
						scales[counter]=scale_i;
						scale_aux=scale_i;
						counter++;}}


		//////////// Writing scales ///////////////

		map_file	= concat("Files/", name);
		map_file  = concat(map_file,"_scales.txt");
		arx       = fopen(map_file,"w");

		for(i=0; i<N_scales-1 ; i++)fprintf(arx,"%d\t",scales[i]);
		fprintf(arx,"%d\n",scales[N_scales-1]);
		fclose(arx);

		////// Name //////

    map_file	= concat("Files/", name);

		char outname[sizeof "_100.txt"];
		sprintf(outname, "_%03d.txt", my_id);
		map_file        = concat(map_file,"_vort");
		map_file        = concat(map_file,outname);
		arx       	    = fopen(map_file,"w");


		pix_counter=0;
		for(i = 0; i < Ngrids; i++){
				if(i*10%Ngrids==0)printf("%d\t %d\n", i, Ngrids);
    		for (j = 0 ; j < Ngrids; j++){
						radius = sqrt((1.0*i-0.5*Ngrids)*(1.0*i-0.5*Ngrids) + (1.0*j-0.5*Ngrids)*(1.0*j-0.5*Ngrids))*L/Ngrids;
						if (radius < Rmax){  // where circulation is going to be computed
						pix_counter++;
						/// which core ///
						if ((pix_counter-my_id)%num_procs==0){
						//printf("%d\t %d\t %d\t %d\t %f\n",i,j,my_id,pix_counter,radius);

						vort	=	vort_map[i][j];
						omeg	=	omeg_map[i][j];
						area	=	1.0;
						fprintf(arx, "%d\t %.4e\t %.3e\t %.3e", pix_counter, radius ,vort, omeg);
						knew	=	1;
						for(k = 1; k < N_scales; k++){

								scale_i = scales[k];
								step		= (int) (0.5*(scale_i-1.0)+0.5);


								do{knew++;
									 step		= (int) (0.5*(knew-1.0)+0.5);

								for(is = -step; is <= step; is++){
									if(abs(is)==step){
										for(js = -step; js <= step; js++){
											if(abs(js)==step){
												if(knew%2!=0){vort+=0.75*vort_map[i+is][j+js];
																			omeg+=0.75*omeg_map[i+is][j+js];
																			area+=0.75;}
												else{ vort+=0.25*vort_map[i+is][j+js];
															omeg+=0.25*omeg_map[i+is][j+js];
															area+=0.25;}
											}
											else{ vort+=0.5*vort_map[i+is][j+js];
														omeg+=0.5*omeg_map[i+is][j+js];
														area+=0.5;}
										}}

									else{	vort+=0.5*vort_map[i+is][j+step];
												omeg+=0.5*omeg_map[i+is][j+step];
												area+=0.5;

											 	vort+=0.5*vort_map[i+is][j-step];
												omeg+=0.5*omeg_map[i+is][j-step];
												area+=0.5;}

								}/*for is*/
								//if(k==20)printf("area %f\t %f\n", sqrt(area), area);
								}while(knew<scale_i);

								fprintf(arx, "\t %.3e \t %.3e" , vort/area, omeg/area);
						 }/*for k*/ fprintf(arx, "\n" );}/*if Rmax*/
					 }/* which core */
				}//if(pix_counter>1000)break;
				}/*for i,j*/

		//fclose(arx);
		ierr=MPI_Finalize();
		return 0;}
