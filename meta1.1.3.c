/*
METAPOPULATION MODEL
VIRAL REPLICATION HAPPENS IN MULTIPLE HOSTS WHO ARE CAPABLE OF TRANSMITTING THE VIRUS TO EACH OTHER IN A RANCOM NETWORK.
in 1.1.3, it uses simple one mutation accumulation equation.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>

// function pre-decleration
void mutate(double**** pop2, double**** pop1, int* curpop2, int* curpop1, double u, int kmax, int host_num, double* N2, double *N1, double factor[], int mutcap);
void reast(double**** pop2, int* curpop, int kmax, int host_num, double r, double* N2);
void repr(double**** pop2, double**** pop1, int* curpop2, int* curpop1, int kmax,
	int host_num, double s, double* N2, double* N1, double* N, double c, double* K, long* seed, double g);
float ran1(long *seed);
float gammln(float xx);
float poidev(float xm,long *idum);
float gasdev(long* idum);
double Nsum(int size,double N[]);
void record(double* N1, double* N2, double**** pop1, double**** pop2,
	int kmax, int host_num, int timestep, int krecord, int curpop1, int curpop2, int rep, int gen, FILE **fPointer);
void migrate(double**** pop2, double**** pop1, int* curpop2, int* curpop1, int kmax, int host_num, double* N2, double* N1, double* N, long* seed, double tr, double mig, double* K, int K0, double* kf, double* kthr, int gen, double kf0);
void evol(double**** pop2, double**** pop1, int* curpop2,
	int* curpop1, int kmax, int host_num, double* N2, double* N1, double* N, long* seed, double evolrate);
double fact(int num);
double C(int n, int k);
double poipmf(double l, int k);

// Macros
#define PRINTF 0 // print pop info at the end of generation
#define POPDIST 0 // print population distribution when XXPRINTF macros are 1.
#define MUTPRINTF 0 // print pop info before and after mutation
#define REAPRINTF 0 // print pop info before and after reassortment
#define REPPRINTF 0	// print pop info before and after reproduction
#define MIGPRINTF 0 // print pop info before and after migration
#define EVOPRINTF 0 // print pop info before and after evolution
#define SHORTCSV 1 //1= only prints out total mean fitness and population. 0= prints out each host's mean fitness and pop
#define EACHSTEPREC 0 // record after every step in a generation (mutate,reast,repr,migrate,evol). The population size info (N) is inaccurate after repr step. Only accurate after migration step.
#define DETERM 0// make the repr process deterministic.


int main(int argc, char *argv[])
{
	// clock start
	clock_t simbegin = clock();
	// prarmeter calling from os command
	char *description = argv[1];
	char *code_s = argv[2];
	char *destination = argv[3]; // directory to save data at
	char *timestep_s = argv[4]; // whether to record every generation
	char *krecord_s = argv[5]; // how to record k. (0=mean k, 1=min k, 2=histogram)
	char *untilext_s = argv[6]; // whether to run the simulation only until extinction (N = 0)
	char *rep_s = argv[7]; // number of repetition
	char *s_s = argv[8]; // selection coefficient
	char *N0_s = argv[9]; // base initial frequency
	char *K0_s = argv[10]; // carrying capacity in a host
	char *u_s = argv[11]; // mutation rate per segment
	char *gen_num_s = argv[12]; // number of generations
	char *c_s = argv[13]; // cost of multisegmentation
	char *r_s = argv[14]; // reassortment probability
	char *seed_s = argv[15]; // seed for random process
	char *host_num_s = argv[16]; // host number
	char *kmax_s = argv[17]; // maximum mutation amount in a segment
	char *pop2init_s = argv[18]; // information on the proportion of initial population to the N0 for 2 segments.
	char *pop2i_l_s = argv[19]; // subsidiary info for extracting pop2init info
	char *pop1init_s = argv[20]; // information on the proportion of initial population to the N0 for 1 segments.
	char *pop1i_l_s = argv[21]; // subsidiary info for extracting pop2init info
	char *tr_s = argv[22]; // transmission rate
	char *mig_s = argv[23]; // migration rate to the migration pool
	char *mutcap_s = argv[24]; // migration rate to the migration pool
	char *kf0_s = argv[25]; // right interval of the uniform distribution of the K fall rate
	char *kfsd_s = argv[26]; // NOT USED ANYMORE
	char *evolrate_s = argv[27]; // evolving rate of 2 segs from 1seg population
	char *g_s = argv[28]; // evolving rate of 2 segs from 1seg population
	char *end1;

	int code = (int) strtol(code_s,&end1,10);
	int timestep = (int) strtol(timestep_s,&end1,10);
	int krecord = (int) strtol(krecord_s, &end1, 10);
	int untilext = (int) strtol(untilext_s,&end1,10);
	int rep = (int) strtol(rep_s,&end1,10);
	double s = (double) strtof(s_s, NULL);
	int N0 = (int) strtol(N0_s,&end1,10);
	int K0 = (int) strtol(K0_s,&end1,10);
	double u = (double) strtof(u_s,NULL);
	int gen_num = (int) strtol(gen_num_s,&end1,10);
	double c = (double) strtof(c_s,NULL);
	double r = (double) strtof(r_s,NULL);
	long seed = strtol(seed_s,&end1,10);
	int host_num = (int) strtol(host_num_s,&end1,10);
	int kmax = (int) strtol(kmax_s,&end1,10);
	int pop2i_l = (int) strtol(pop2i_l_s,&end1,10);
	int pop1i_l = (int) strtol(pop1i_l_s,&end1,10);
	double tr = (double) strtof(tr_s, NULL);
	double mig = (double) strtof(mig_s, NULL);
	int mutcap = (int) strtol(mutcap_s,&end1,10);
	double kf0 = (double) strtof(kf0_s, NULL);
	double kfsd = (double) strtof(kfsd_s, NULL);
	double evolrate = (double) strtof(evolrate_s, NULL);
	double g = (double) strtof(g_s, NULL);

	int k,i,j,m;
	//double pop1init[host_num];
	double pop2init[host_num];
	double pop1init[host_num];
	int len = 0;

	char* dubstr = (char*) malloc(sizeof(char)*10);
	sprintf(dubstr,"");
	for (i=0; i < pop2i_l; i++)
	{

		if (pop2init_s[i] != '~')
		{
			sprintf(dubstr,"%s%c",dubstr,pop2init_s[i]);
		}
		else
		{
			pop2init[len] = (double) strtof(dubstr, NULL);
			len++;
			sprintf(dubstr,"");
		}
	}
	len = 0;
	for (i=0; i < pop1i_l; i++)
	{
		if (pop1init_s[i] != '~')
		{
			sprintf(dubstr,"%s%c",dubstr,pop1init_s[i]);
		}
		else
		{
			pop1init[len] = (double) strtof(dubstr, NULL);
			len++;
			sprintf(dubstr,"");
		}
	}
	free(dubstr);
	printf("destination=%s, timestep=%d, krecord=%d, hostnum=%d, untilext=%d, kamx=%d, rep=%d, s=%.2f, N0=%d, K0=%d, u=%.5f, gen_num=%d, c=%.7f, r=%.2f, tr=%.5f, mig=%.5f, mutcap=%d, seed=%ld, kf=%.2f, kfsd=%.2f, evolrate=%.7f,g=%.1f\n",destination,timestep,krecord,host_num,untilext,kmax,rep,s,N0,K0,u,gen_num,c,r,tr,mig,mutcap,seed,kf0,kfsd,evolrate,g);

	//check if the destination folder exists and if not, make one.
	char* dest2 = (char*) malloc(sizeof(char)*50);
	sprintf(dest2,"./data/%s",destination);
	struct stat st = {0};
	if (stat(dest2,&st) == -1)
	{
		mkdir(dest2, 0700);
	}
	char* paraminfo = (char*) malloc(sizeof(char)*1000);
	sprintf(paraminfo,"%d,%d,%d,%.3f,%d,%d,%.5f,%d,%.7f,%.2f,%d,%d,%.5f,%.5f,%.2f,%.2f,%.7f,%.1f",timestep,krecord,rep,s,N0,K0,u,gen_num,c,r,kmax,host_num,mig,tr,kf0,kfsd,evolrate,g);
	char* filename = (char*) malloc(sizeof(char)*1000);
	//sprintf(filename,"%s/m1.1.3s_%d,%d,%d,%.3f,%d,%d,%.5f,%d,%.2f,%.2f,%d,%d,%.5f,%.5f,%.2f,%.2f,%.2f(0).csv",dest2,timestep,krecord,rep,s,N0,K0,u,gen_num,c,r,kmax,host_num,mig,tr,kf0,kfsd,evolrate);
	sprintf(filename,"%s/m1.1.3s_%s(0).csv",dest2,paraminfo);
	int filenum = 0;

	while ( access(filename, F_OK) != -1) 
	{
		filenum += 1;
		sprintf(filename,"%s/m1.1.3s_%s(%d).csv",dest2,paraminfo,filenum);
	}
	FILE* fPointer;
	fPointer = fopen(filename,"w");
	free(dest2);
	free(filename);
	fprintf(fPointer,"#%s\n",description); // putting in comments on the csv
	fprintf(fPointer,"#%s.untilext=%d.SEED=%ld.code=%d\n",paraminfo,untilext,seed,code);
	char* str;


	if (SHORTCSV == 0)
	{
		if (krecord == 1 || krecord == 0)
		{
			str = (char*) malloc(sizeof(char)*(20*(host_num+1)*2) + 100); // (words + comma)*(host_num+1)*(seg1 and 2) + safety buffer
		}
		else if (krecord == 2)
		{
			str = (char*) malloc(sizeof(char)*(20*(2*kmax+1)*2) + 100); // (words + comma)*(host_num)*(seg1 and 2) + safety buffer
		}
		sprintf(str,"");
		if(krecord != 2)
		{
			for(i=0; i<=host_num; i++)
			{
				sprintf(str,"%s,pop1.%d,pop2.%d,k1.%d,k2.%d",str,i,i,i,i);	
			}
			sprintf(str,"%s,pop1.m,pop2.m,k1.m,k2.m",str);
		}
		else
		{
			sprintf(str,"%s,pop1.0,pop2.0",str);
			for (i=0; i<=kmax*2; i++)
			{
				sprintf(str,"%s,p1k%d,p2k%d",str,i,i);
			}
		}
		fprintf(fPointer,"rep,gen%s\n",str);
	}
	else
	{
		if (krecord == 1 || krecord == 0)
		{
			str = (char*) malloc(sizeof(char)*(20*(host_num+1)*2) + 100); // (words + comma)*(host_num+1)*(seg1 and 2) + safety buffer
		}
		else if (krecord == 2)
		{
			str = (char*) malloc(sizeof(char)*(20*(2*kmax+1)*2) + 100); // (words + comma)*(host_num)*(seg1 and 2) + safety buffer
		}
		if (krecord != 2)
		{
			fprintf(fPointer,"rep,gen,pop1.0,pop2.0,k1.0,k2.0\n");
		}
		else
		{
			sprintf(str,"%s,pop1.0,pop2.0",str);
			for (i=0; i<=kmax*2; i++)
			{
				sprintf(str,"%s,p1k%d,p2k%d",str,i,i);
			}
			fprintf(fPointer,"rep,gen%s\n",str);
		}
	}

	sprintf(str,"done");
	free(paraminfo);
	free(str);
	
	int gen,repe; //current generation and repetition.
	int curpop1, curpop2; //current population index ur working witih.
	double**** pop2; //population distribution array of 2 segments
	double**** pop1; //population distribution array of 1 segments
	double* N = (double*) malloc(sizeof(double)*(host_num+1)); // total population size. N[i] is a pop size of host i. N[0] is the total population size summing all hosts.
	double* N2 = (double*) malloc(sizeof(double)*(host_num+1)); // population size of 2segs. N2[i] is a segmented pop size of host i. N2[0] is the total population size of 2segs.
	double* N1 = (double*) malloc(sizeof(double)*(host_num+1)); // population size of 1segs. N1[i] is a non-segmented pop size of host i. N1[0] is the total population size of 1segs.
	double* K = (double*) malloc(sizeof(double)*(host_num));
	double* kf = (double*) malloc(sizeof(double)*(host_num));
	double* kthr = (double*) malloc(sizeof(double)*(host_num)); // array keeping track if the threshold for immunity process has been broken or not.
	for (i=0; i<host_num; i++)
	{
		kf[i] = ran1(&seed)*kf0;
		kthr[i] = 1.0;
	}
	double factor[2*kmax]; // probability of getting n mutations organized in an array

	// initialize population array
	// pop2: 2segmented population. pop[i][m][j][k] = number of individual with j mutation in 1st segment and k mutation in 2nd segment in host m.
	// i alternate between 0 and 1 after every step. In each step of a generation, the changed information is stored in a opposite index.
	// For example, before mutation step, pop[0][m][j][k] holds the information, but after the mutation step, pop[1][m][j][k] holds the changed information.
	// After reassortment step, pop[0][m][j][k] holds the latest information again.
	// host with index 0 is a common pool used in migration step.
	// pop1: same as pop2. Has 2 segments for the easement of transition into segmented virus during evolve step.
	// host 0 (pop1[i][0][j] or pop2[i][0][j][k]) is a migration pool.
	pop2 = (double****) malloc(sizeof(double***)*2);
	for (m=0; m<2; m++)
	{
		pop2[m] = (double***) malloc(sizeof(double**)*(host_num+1));
		for (i=0; i<=host_num; i++)
		{
			pop2[m][i] = (double**) malloc(sizeof(double*)*(kmax+1));
			for (j=0; j<=kmax; j++)
			{
				pop2[m][i][j] = (double*) malloc(sizeof(double)*(kmax+1));
			}
		}
	}


	pop1 = (double****) malloc(sizeof(double***)*2);
	for (m=0; m<2; m++)
	{
		pop1[m] = (double***) malloc(sizeof(double**)*(host_num+1));
		for (i=0; i<=host_num; i++)
		{
			pop1[m][i] = (double**) malloc(sizeof(double*)*(kmax+1));
			for (j=0; j<=kmax; j++)
			{
				pop1[m][i][j] = (double*) malloc(sizeof(double)*(kmax+1));
			}
		}
	}

	for (i=0; i<=2*kmax; i++){
		factor[i] = poipmf(2*u,i); // probability of choosing i amount of mutation in a generation step.
	}

	double mutate_time = 0;
	double reast_time = 0;
	double repr_time = 0;
	double migr_time = 0;
	double evol_time = 0;
	double reco_time = 0;
	clock_t begin, end;
	int poistart = 0; // initial population poisson distribution if 1. all k0 class if 0.
	int unistart = 0; // initial population uniform distribution if 1. all k0 class if 0.
	for (repe=0; repe < rep; repe++)
	{
		/* //print rep counts
		if (rep > 100)
		{
			if (repe % 100 == 0)
			{
				printf("\rREP = %d",repe);
				fflush(stdout);
			}
		}
		else
		{
			printf("\rREP = %d",repe);
			fflush(stdout);
		}
		*/
		// initializing N and pop
		if (poistart) // initialize pop with poisson distribution with rate of expected equilibria, 2u/s.
		{
			for(i=0; i<=host_num; i++)
			{
				N2[i] = 0;
				N1[i] = 0;
				N[i] = 0;
				for (j=0; j<=kmax; j++)
				{
					for (k=0; k<=kmax; k++)
					{
						if (i!=0)
						{
							//pop2[0][i][j][k] = 0;
							//pop1[0][i][j][k] = 0;
							if (j+k <= kmax)
							{
								pop1[0][i][j][k] = (double) N0*pop1init[i-1]*poipmf((2*u)/s,(j+k))/(j+k+1);
							}
							else
							{
								pop1[0][i][j][k] = (double) N0*pop1init[i-1]*poipmf((2*u)/s,(j+k))/(kmax+1-((j+k)-kmax));
							}
							//printf("pop1[0][%d][%d][%d]=%.3f\n",i,j,k,(double) N0*pop1init[i-1]*poipmf((2*u)/s,(j+k)));
							N1[i] += pop1[0][i][j][k];
							if (evolrate == 0)
							{
								if (j+k <= kmax)
								{
									pop2[0][i][j][k] = (double) N0*pop2init[i-1]*poipmf((2*u)/s,(j+k))/(j+k+1);
								}
								else
								{
									pop2[0][i][j][k] = (double) N0*pop2init[i-1]*poipmf((2*u)/s,(j+k))/(kmax+1-((j+k)-kmax));
								}
								//printf("pop2[0][%d][%d][%d]=%.3f\n",i,j,k,(double) N0*pop1init[i-1]*poipmf((2*u)/s,(j+k)));
								N2[i] += pop2[0][i][j][k];
							}
							else
							{
								pop2[0][i][j][k] = 0;
							}
						}
						else
						{
							pop2[0][i][j][k] = 0;
							pop1[0][i][j][k] = 0;
						}
					}
				}
				if (i!=0)
				{
					if (evolrate == 0)
					{
						N2[i] = 0;
					}
					N[i] = N1[i] + N2[i];
					K[i-1] = (double) K0;
				}
			}
		}
		else if (unistart) // initialize with uniform dist among classes
		{
			N2[i] = 0;
			N1[i] = 0;
			N[i] = 0;
			for(i=0; i<=host_num; i++)
			{
				for (j=0; j<=kmax; j++)
				{
					for (k=0; k<=kmax; k++)
					{
						if (i!=0)
						{
							pop2[0][i][j][k] = ((double)N0*pop2init[i-1])/((double)((kmax+1)*(kmax+1)));
							pop1[0][i][j][k] = ((double)N0*pop1init[i-1])/((double)((kmax+1)*(kmax+1)));
						}
					}
				}
				if (i!=0)
				{
					N1[i] = (double) N0*pop1init[i-1];
					N2[i] = (double) N0*pop2init[i-1];
					N[i] = N1[i] + N2[i];
					K[i-1] = (double) K0;
				}
			}	
		}
		else // default start. Only (0,0) class in the beginning.
		{
			for(i=0; i<=host_num; i++)
			{
				for (j=0; j<=kmax; j++)
				{
					for (k=0; k<=kmax; k++)
					{
						if (i!=0)
						{
							pop2[0][i][j][k] = 0;
							pop1[0][i][j][k] = 0;
						}
					}
				}
				if (i!=0)
				{
					pop1[0][i][0][0] = (double) N0*pop1init[i-1];
					N1[i] = (double) N0*pop1init[i-1];
					if (evolrate > 0)
					{
						N2[i] = 0;
					}
					else
					{
						pop2[0][i][0][0] = (double) N0*pop2init[i-1]; //N0 of virus with 0 mutations at initial condition.
						N2[i] = (double) N0*pop2init[i-1];
					}
					N[i] = N1[i] + N2[i];
					K[i-1] = (double) K0;
				}
			}
		}
		for (i=1; i<=host_num; i++)
		{
			//printf("N[%d]=%f, N1[i]=%f, N2[i]=%f\n",i,N[i],N1[i],N2[i]);
		}

		N1[0] = Nsum(host_num,N1);
		N2[0] = Nsum(host_num,N2);
		N[0] = N1[0] + N2[0];

		//N1[0] = Nsum(host_num,N1);
		curpop1 = 0;
		curpop2 = 0;
		for (gen=0; gen < gen_num; gen++)
		{
			begin = clock();
			mutate(pop2, pop1, &curpop2, &curpop1, u, kmax, host_num, N2, N1, factor, mutcap);
			if(EACHSTEPREC)
			{
				record(N1,N2,pop1,pop2,kmax,host_num,timestep,krecord,curpop1,curpop2,repe,gen,&fPointer);
			}
			end = clock();
			mutate_time += (double) (end - begin)/ CLOCKS_PER_SEC;
			begin = clock();
			reast(pop2, &curpop2, kmax, host_num, r, N2);
			if(EACHSTEPREC)
			{
				record(N1,N2,pop1,pop2,kmax,host_num,timestep,krecord,curpop1,curpop2,repe,gen,&fPointer);
			}
			end = clock();
			reast_time += (double) (end - begin)/ CLOCKS_PER_SEC;
			begin = clock();
			repr(pop2, pop1, &curpop2, &curpop1, kmax, host_num,s,N2,N1,N,c,K,&seed,g);
			if(EACHSTEPREC)
			{
				record(N1,N2,pop1,pop2,kmax,host_num,timestep,krecord,curpop1,curpop2,repe,gen,&fPointer);
			}
			end = clock();
			repr_time += (double) (end - begin)/ CLOCKS_PER_SEC;				
			begin = clock();
			migrate(pop2, pop1, &curpop2, &curpop1, kmax, host_num, N2, N1, N, &seed, tr, mig, K, K0, kf, kthr, gen, kf0);
			end = clock();
			migr_time += (double) (end - begin)/ CLOCKS_PER_SEC;
			if (evolrate > 0 && gen >= 300)
			{
				begin = clock();
				evol(pop2,pop1, &curpop2, &curpop1, kmax, host_num, N2, N1, N, &seed, evolrate);
				end = clock();
				evol_time += (double) (end - begin)/ CLOCKS_PER_SEC;
			}
			N[0] = Nsum(host_num,N);
			//if (gen > 300)
			//{
				//printf("gen=%d",gen);
				//printf("------------\n");
				//getchar();
			//}
			//printf("N[0]=%.3f\n",N[0]);
			//printf("N=%.3f N[1]=%.3f, N[2]=%.3f N[3]=%.3f N[4]=%.3f N[5]=%.3f\n",N[0],N[1],N[2],N[3],N[4],N[5]);


			if(PRINTF)
			{
				float count1, count2;
				count1 = 0;
				count2 = 0;
				printf("after repr outside the function\n");
				for (i=1; i<=host_num; i++)
				{
					for (j=0; j<=kmax; j++)
					{
						for (k=0; k<=kmax; k++)
						{
							//if (pop2[curpop2][i][j][k]>0)
							{
								printf("pop2[%d][%d][%d][%d]=%.3f\n",curpop2,i,j,k,pop2[curpop2][i][j][k]);
							}
							count2 += pop2[curpop2][i][j][k];
						}
					}
					for (j=0; j<=kmax; j++)
					{
						for (k=0; k<=kmax; k++)
						{
							if (pop1[curpop1][i][j][k]>0)
							{
								//printf("pop2[%d][%d][%d][%d]=%.3f\n",curpop2,i,j,k,pop1[curpop1][i][j][k]);
							}
							count1 += pop1[curpop1][i][j][k];
						}
					}
				}
				printf("count2=%.3f\n",count2);
				printf("count1=%.3f\n",count1);
				printf("N1[0]=%.3f\n",N1[0]);
				printf("N2[0]=%.3f\n",N2[0]);
			}
			
			//printf("gen=%d\n",gen);
			//printf("-----------------\n");
			// record it to the file
			if(timestep == 1)
			{
				begin = clock();
				record(N1,N2,pop1,pop2,kmax,host_num,timestep,krecord,curpop1,curpop2,repe,gen,&fPointer);
				end = clock();
				reco_time += (double) (end - begin)/ CLOCKS_PER_SEC;
			}
			// erase migration info.
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					pop1[curpop1][0][j][k] = 0;
					pop2[curpop2][0][j][k] = 0;
				}
			}


			if(untilext == 1)
			{
				if(evolrate == 0)
				{
						if (N2[0] == 0 || N1[0] == 0)
						{
							break;
						}
				}
				else
				{
						if (N1[0] == 0) // ** special edition for measuring waiting time. //N1[0] == 0 (float)N2[0]/(float)(N1[0]+N2[0]) >= 0.9
						{
							break;
						}
				}
			}
		}
		if(timestep == 0)
		{
			//begin = clock();
			record(N1,N2,pop1,pop2,kmax,host_num,timestep,krecord,curpop1,curpop2,repe,gen,&fPointer);
			//end = clock();
			//reco_time += (double) (end - begin)/ CLOCKS_PER_SEC;				
		}
		// erase pop data to start over simulation
		for (m=0; m<2; m++)
		{
			for (i=0; i<=host_num; i++)
			{
				for (j=0; j<=kmax; j++)
				{
					for (k=0; k<=kmax; k++)
					{
						pop2[m][i][j][k] = 0;
						pop1[m][i][j][k] = 0;
					}
				}
			}
		}
	}
	// freeing pop2
	for (m=0; m<2; m++)
	{
		for (i=0; i<=host_num; i++)
		{
			for (j=0; j<=kmax; j++)
			{
				free(pop2[m][i][j]);
				free(pop1[m][i][j]);
			}
			free(pop2[m][i]);
			free(pop1[m][i]);
			}
		free(pop2[m]);	
		free(pop1[m]);	
	}
	free(pop2);
	free(pop1);
	free(N2);
	free(N1);
	free(N);
	free(kf);
	free(K);
	// make mutation, recombination, reproduction into a single process.
	fclose(fPointer);
	clock_t simend = clock();
	double time_spent = (simend - simbegin)/ CLOCKS_PER_SEC;
	printf("\n");
	printf("mutate=%.3f, reast=%.3f, repr=%.3f, migr=%.3f\n",mutate_time, reast_time, repr_time, migr_time);
	if (evolrate > 0)
	{
		printf("evolve=%.3f\n",evol_time);
	}
	printf("record=%.3f\n",reco_time);
	printf("time spend was %.2f minutes\n", time_spent/60.0);
	return 0;
}


void mutate(double**** pop2, double**** pop1, int* curpop2, int* curpop1, double u, int kmax, int host_num, double* N2, double *N1, double factor[], int mutcap)
{
	int s2m,s2m2, s1m, s1m2, i,j,k,l,l2,l3;
	if (*curpop2 == 0)
	{
		s2m2 = *curpop2;
		s2m = 1;
		*curpop2 = s2m;
	}
	else 
	{
		s2m2 = *curpop2;
		s2m = 0;
		*curpop2 = s2m;
	}

	if (*curpop1 == 0)
	{
		s1m2 = *curpop1;
		s1m = 1;
		*curpop1 = s1m;
	}
	else 
	{
		s1m2 = *curpop1;
		s1m = 0;
		*curpop1 = s1m;
	}

	double count1, count2;
	double mean1, mean2;
	if(MUTPRINTF)
	{
		count2 = 0;
		count1 = 0;
		mean2 = 0;
		mean1 = 0;
		printf("before mutation inside the function\n");
		for (i=1; i<=host_num; i++)
		{
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					//if (pop2[s2m2][i][j][k] > 0)
					{
						if(POPDIST)
						{
							printf("pop2[%d][%d][%d][%d]=%f\n",s2m2,i,j,k,pop2[s2m2][i][j][k]);
						}
					}
					count2 += pop2[s2m2][i][j][k];
					mean2 += pop2[s2m2][i][j][k]*(j+k);
				}
			}
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					if (pop1[s1m2][i][j][k] > 0)
					{
						if (POPDIST)
						{
							printf("pop1[%d][%d][%d][%d]=%f\n",s1m2,i,j,k,pop1[s1m2][i][j][k]);
						}
					}
					count1 += pop1[s1m2][i][j][k];
					mean1 += pop1[s1m2][i][j][k]*(j+k);
				}
			}
		}
		printf("count1=%.3f\n",count1);
		printf("count2=%.3f\n",count2);
		printf("mean1=%.3f\n",mean1/count1);
		printf("mean2=%.3f\n",mean2/count2);
		printf("\n");
	}	

	double factor2 = 0;
	int left; // max of number of mutations that n(l,k) can give rise to.
	int cap;
	double lsum, binomprob;
	for (i=1; i<=host_num; i++)
	{
		// 2 seg mutation process
		if (N2[i] > 0) 
		{
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					// going to sum all the mutation_rate*n(l,k)
					left = 2*kmax - (j + k);
					pop2[s2m][i][j][k] += pop2[s2m2][i][j][k];
					if (left < mutcap) // if possible number of accumulation is smaller than mutcap, the max number of accumulation is the cap.
					{
						cap = left; 
					}
					else
					{
						cap = mutcap;
					}
					for(l=1; l<=cap; l++) // add different number of mutations with this loop
					{
						factor2 = factor[l]*pop2[s2m2][i][j][k];
						pop2[s2m][i][j][k] -= factor2;
						lsum = 0;
						binomprob = pow(0.5,l);
						for(l2=0; l2<=l; l2++) 
						{
							l3 = l - l2;
							if(l2 + j <= kmax && l3 + k <= kmax)
							{
								lsum += C(l,l2)*binomprob;
							}
						}
						for(l2=0; l2<=l; l2++) // different ways to give same number of mutations in this loop
						{ // l2= #of mut in one segment. l3= " in another segment
							l3 = l - l2;
							if(l2 + j <= kmax && l3 + k <= kmax)
							{
								pop2[s2m][i][j+l2][k+l3] += factor2*(C(l,l2)*binomprob/lsum);
							}
						}
					}
				}
			}
		}
		// 1seg mutation process
		if (N1[i] > 0) 
		{
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					// going to sum all the mutation_rate*n(l,k)
					left = 2*kmax - (j + k);
					pop1[s1m][i][j][k] += pop1[s1m2][i][j][k];
					if (left < mutcap) // if possible number of accumulation is smaller than mutcap, the max number of accumulation is the cap.
					{
						cap = left; 
					}
					else
					{
						cap = mutcap;
					}
					for(l=1; l<=cap; l++)
					{
						factor2 = factor[l]*pop1[s1m2][i][j][k];
						pop1[s1m][i][j][k] -= factor2;
						lsum = 0;
						binomprob = pow(0.5,l);
						for(l2=0; l2<=l; l2++) 
						{
							l3 = l - l2;
							if(l2 + j <= kmax && l3 + k <= kmax)
							{
								lsum += C(l,l2)*binomprob;
							}
						}
						for(l2=0; l2<=l; l2++) // different ways to give same number of mutations in this loop
						{ // l2= #of mut in one segment. l3= " in another segment
							l3 = l - l2;
							if(l2 + j <= kmax && l3 + k <= kmax)
							{
								pop1[s1m][i][j+l2][k+l3] += factor2*(C(l,l2)*binomprob/lsum);
							}
						}
					}
				}
			}
		}
	}

	if(MUTPRINTF)
	{
		count2 = 0;
		count1 = 0;
		mean2 = 0;
		mean1 = 0;
		printf("after mutation inside the function\n");
		for (i=1; i<=host_num; i++)
		{
			for (j=0; j<=kmax; j++)
			{	
				for (k=0; k<=kmax; k++)
				{
					if (pop2[s2m][i][j][k] > 0)
					{
						if(POPDIST)
						{
							printf("pop2[%d][%d][%d][%d]=%.2f\n",s2m,i,j,k,pop2[s2m][i][j][k]);
						}
					}
					//pop2[m2][i][j][k] = 0;
					count2 += pop2[s2m][i][j][k];
					mean2 += pop2[s2m][i][j][k]*(j+k);
				}
			}
			for (j=0; j<=kmax; j++)
			{	
				for (k=0; k<=kmax; k++)
				{
					if (pop1[s1m][i][j][k] > 0)
					{
						if(POPDIST)
						{
							printf("pop1[%d][%d][%d][%d]=%.2f\n",s1m,i,j,k,pop1[s1m][i][j][k]);
						}
					}
					//pop2[m2][i][j][k] = 0;
					count1 += pop1[s1m][i][j][k];
					mean1 += pop1[s1m][i][j][k]*(j+k);
				}
			}
		}
		printf("count1=%.3f\n",count1);
		printf("count2=%.3f\n",count2);
		printf("mean1=%.3f\n",mean1/count1);
		printf("mean2=%.3f\n",mean2/count2);
		printf("\n");
	}
}


void reast(double**** pop2, int* curpop, int kmax, int host_num, double r, double* N2)
{
	int m,m2,i,j,k;
	int jj;
	if (*curpop == 0)
	{
		m2 = *curpop;
		m = 1;
		*curpop = m;
	}
	else 
	{
		m2 = *curpop;
		m = 0;
		*curpop = m;
	}
	double count2;
	double mean2 = 0;
	double meanfit2 = 0; // mean fitness
	if(REAPRINTF)
	{
		count2 = 0;
		printf("before reast inside the function\n");
		for (i=1; i<=host_num; i++)
		{
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					if(POPDIST)
					{
						printf("pop2[%d][%d][%d][%d]=%.2f\n",m2,i,j,k,pop2[m2][i][j][k]);
					}
					count2 += pop2[m2][i][j][k];
					mean2 += pop2[m2][i][j][k]*(j+k);
					meanfit2 += pop2[m2][i][j][k]*pow(0.95,(j+k));
				}
			}
		}
		mean2 = mean2 / count2;
		meanfit2 = meanfit2 / count2;
		printf("count2=%.3f\n",count2);
		printf("mean2=%.3f\n",mean2);
		printf("meanfit2=%.3f\n",meanfit2);

	}	
	double jp[kmax+1];
	double kp[kmax+1];
	for (i=1; i<=host_num; i++)
	{
		if (N2[i] > 0)
		{
			// calculate all the proportions first
			for(j=0; j<=kmax; j++)
			{
				kp[j] = 0;
				jp[j] = 0;
				for(jj=0; jj<=kmax; jj++)
					{
						kp[j] += pop2[m2][i][jj][j];
						jp[j] += pop2[m2][i][j][jj];
					}
					kp[j] = kp[j]/N2[i];
					jp[j] = jp[j]/N2[i];
			}
			// those that didn't recombine remain, and the recombined ones get added.
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					pop2[m][i][j][k] = pop2[m2][i][j][k]*(1 - r);
					pop2[m][i][j][k] += N2[i]*kp[k]*jp[j]*r;
				}
			}
		}
	}


	if (REAPRINTF)
	{
		count2 = 0;
		mean2 = 0;
		meanfit2 = 0;
		printf("after reast inside the function\n");
		for (i=1; i<=host_num; i++)
		{
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					pop2[m2][i][j][k] = 0;
					if(POPDIST)
					{
						printf("pop2[%d][%d][%d][%d]=%.2f\n",m,i,j,k,pop2[m][i][j][k]);
					}
					count2 += pop2[m][i][j][k];
					mean2 += pop2[m][i][j][k]*(j+k);
					meanfit2 += pop2[m][i][j][k]*pow(0.95,(j+k));
				}
			}
		}
		mean2 = mean2/count2;
		meanfit2 = meanfit2/count2;
		printf("count2=%.3f\n",count2);
		printf("mean2=%.3f\n",mean2);
		printf("meanfit2=%.3f\n\n\n\n\n",meanfit2);
	}
}

void repr(double**** pop2, double**** pop1, int* curpop2, int* curpop1, int kmax, int host_num, double s, double* N2, double* N1, double* N, double c, double* K, long* seed, double g)
{
	int s2m,s2m2,s1m,s1m2,i,j,k;
	if (*curpop2 == 0)
	{
		s2m2 = *curpop2;
		s2m = 1;
		*curpop2 = s2m;
	}
	else 
	{
		s2m2 = *curpop2;
		s2m = 0;
		*curpop2 = s2m;
	}
	if (*curpop1 == 0)
	{
		s1m2 = *curpop1;
		s1m = 1;
		*curpop1 = s1m;
	}
	else 
	{
		s1m2 = *curpop1;
		s1m = 0;
		*curpop1 = s1m;
	}

	//double newN2, newN1;
	double count1, count2;
	double mean1, mean2;
	if (REPPRINTF)
	{
		count1 = 0;
		count2 = 0;
		mean1 = 0;
		mean2 = 0;
		printf("before repr inside the function\n");
		for (i=1; i<=host_num; i++)
		{
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					if (pop2[s2m2][i][j][k] > 0)
					{
						if(POPDIST)
						{
							printf("pop2[%d][%d][%d][%d]=%.2f\n",s2m2,i,j,k,pop2[s2m2][i][j][k]);
						}
						count2 += pop2[s2m2][i][j][k];
						mean2 += pop2[s2m2][i][j][k]*(j+k);
					}
				}
			}
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					if (pop1[s1m2][i][j][k] > 0)
					{
						if(POPDIST)
						{
							printf("pop1[%d][%d][%d][%d]=%.2f\n",s1m2,i,j,k,pop1[s1m2][i][j][k]);
						}
						count1 += pop1[s1m2][i][j][k];
						mean1 += pop1[s1m2][i][j][k]*(j+k);
					}
				}
			}
		}
		printf("count1=%.3f\n",count1);
		printf("count2=%.3f\n",count2);
		printf("mean1=%.3f\n",mean1/count1);
		printf("mean2=%.3f\n",mean2/count2);
	}
	int determ = DETERM; // choose weather to run the model deterministic or stochastic.
	// stochastic mode
	double poirate;
	//double g = 5;
	if (determ == 0)
	{
		for (i=1; i<=host_num; i++)
		{
			if(N2[i] > 0)
			{
				//newN2 = 0; //reset N2 to 0
				for (j=0; j<=kmax; j++)
				{
					for (k=0; k<=kmax; k++)
					{
						// if total mutation number = 2*kmax the reproduction rate is 0.
						if (k + j == 2*kmax)
						{
							pop2[s2m][i][j][k] = poidev(0,seed);
						}
						else
						{
							poirate = pop2[s2m2][i][j][k] * pow((1 - s), (k + j)) * (1 - c) * ((double)(1 + g) / (1.0 + g*(N[i]/K[i-1])));
							pop2[s2m][i][j][k] = poidev(poirate,seed);
						}
					}
				}
			}
			if(N1[i] > 0)
			{
				for (j=0; j<=kmax; j++)
				{
					for (k=0; k<=kmax; k++)
					{
						// if total mutation number = 2*kmax the reproduction rate is 0.
						if (k + j == 2*kmax)
						{
							pop1[s1m][i][j][k] = poidev(0,seed);
						}
						else
						{
							poirate = pop1[s1m2][i][j][k] * pow((1 - s), (k + j)) * ((double)(1 + g) / (1.0 + g*(N[i]/K[i-1])));
							pop1[s1m][i][j][k] = poidev(poirate,seed);
						}
					}
				}
			}
		}
	}
	else // FOR DETERMINISTIC SYSTEM
	{
		double numerator = 0;
		double numerator2 = 0;
		for (j=0; j<=kmax; j++)
		{
			for (k=0; k<=kmax; k++)
			{
				// if total mutation number = 2*kmax the reproduction rate is 0.
				numerator += pop1[s1m2][1][j][k]*pow((1-s),(j+k));
				numerator2 += pop2[s2m2][1][j][k]*pow((1-s),(j+k));
			}
		}
		for (j=0; j<=kmax; j++)
		{
			for (k=0; k<=kmax; k++)
			{
				pop1[s1m][1][j][k] = pop1[s1m2][1][j][k]*pow((1-s),(j+k))/numerator * N1[1];
				pop2[s2m][1][j][k] = pop2[s2m2][1][j][k]*pow((1-s),(j+k))/numerator2 * N2[1];
			}
		}
	}
	
	if(REPPRINTF)
	{
		count1 = 0;
		count2 = 0;
		mean1 = 0;
		mean2 = 0;
		printf("after repr inside the function\n");
	}
	
	// erase m2 info in order to make mutate function work on the next generation.
	for (i=0; i<=host_num; i++)
	{
		for (j=0; j<=kmax; j++)
		{
			for (k=0; k<=kmax; k++)
			{
				//pop2[s2m2][i][j][k] = 0;
				//if (pop2[m][i][j][k] >0 )
				{
					if (REPPRINTF)
					{
						if (pop2[s2m][i][j][k] > 0)
						{
							if(POPDIST)
							{
								printf("pop2[%d][%d][%d][%d]=%.2f\n",s2m,i,j,k,pop2[s2m][i][j][k]);
							}
							count2 += pop2[s2m][i][j][k];
							mean2 += pop2[s2m][i][j][k]*(j+k);
						}
					}
				}
			}
		}
		for (j=0; j<=kmax; j++)
		{
			for (k=0; k<=kmax; k++)
			{
				if (REPPRINTF)
				{
					if (pop1[s1m][i][j][k] > 0)
					{
						if(POPDIST)
						{
							printf("pop1[%d][%d][%d][%d]=%.2f\n",s1m,i,j,k,pop1[s1m][i][j][k]);
						}
						count1 += pop1[s1m][i][j][k];
						mean1 += pop1[s1m][i][j][k]*(j+k);
					}
				}
			}
		}

	}
	if (REPPRINTF)
	{
		printf("count1=%.3f\n",count1);
		printf("count2=%.3f\n",count2);
		printf("mean1=%.3f\n",mean1/count1);
		printf("mean2=%.3f\n",mean2/count2);
	}
}

void migrate(double**** pop2, double**** pop1, int* curpop2, int* curpop1, int kmax, int host_num, double* N2, double* N1, double* N, long* seed, double tr, double mig, double* K, int K0, double* kf, double* kthr, int gen, double kf0)
{
	int s2m,s2m2,s1m,s1m2,i,j,k;
	double count1, count2;
	if (host_num > 1)
	{
		if (*curpop2 == 0)
		{
			s2m2 = *curpop2;
			s2m = 1;
			*curpop2 = s2m;
		}
		else 
		{
			s2m2 = *curpop2;
			s2m = 0;
			*curpop2 = s2m;
		}
		if (*curpop1 == 0)
		{
			s1m2 = *curpop1;
			s1m = 1;
			*curpop1 = s1m;
		}
		else 
		{
			s1m2 = *curpop1;
			s1m = 0;
			*curpop1 = s1m;
		}
	}
	else
	{
		s2m = *curpop2;
		s1m = *curpop1;
		if(*curpop1 == 0)
		{
			s1m2 = 1;
		}
		else
		{
			s1m2 = 0;
		}
		if(*curpop2 == 0)
		{
			s2m2 = 1;
		}
		else
		{
			s2m2 = 0;
		}
	}	
	double newN2, newN1;

	if (MIGPRINTF)
	{
		if (host_num > 1)
		{
			printf("before migration\n");
			count1 = 0;
			count2 = 0;
			for (i=0; i<=host_num; i++)
			{
				for (j=0; j<=kmax; j++)
				{
					for (k=0; k<=kmax; k++)
					{
						if (pop2[s2m2][i][j][k] > 0)
						{
							//if(POPDIST)
							{
								printf("pop2[%d][%d][%d][%d]=%.2f\n",s2m2,i,j,k,pop2[s2m2][i][j][k]);
							}
							count2 += pop2[s2m2][i][j][k];
						}
					}
				}
				for (j=0; j<=kmax; j++)
				{
					for (k=0; k<=kmax; k++)
					{
						if (pop1[s1m2][i][j][k] > 0)
						{
							//if(POPDIST)
							{
								printf("pop1[%d][%d][%d][%d]=%.2f\n",s1m2,i,j,k,pop1[s1m2][i][j][k]);
							}
							count1 += pop1[s1m2][i][j][k];
						}
					}
				}
			}
			printf("count2=%.3f\n",count2);
			printf("count1=%.3f\n",count1);
		}
		else
		{
			printf("only 1 host no migration\n");
		}
	}
	if (host_num > 1)
	{
		// going into the migration pool (pop[m][0])
		for(i=1; i<= host_num; i++)
		{
			if (N2[i] > 0)
			{
				for (j=0; j<=kmax; j++)
				{
					for (k=0; k<=kmax; k++)
					{
						pop2[s2m][i][j][k] = pop2[s2m2][i][j][k];
						pop2[s2m][i][j][k] -= pop2[s2m2][i][j][k] * mig;
						if (i==1)
						{
							pop2[s2m][0][j][k] = 0;
						}
						pop2[s2m][0][j][k] += pop2[s2m2][i][j][k] * mig;
						//printf("pop2[s2m2][%d][%d][%d]=%.2f,mig=%.2f\n",i,j,k,pop2[s2m2][i][j][k],mig);
					}
				}
			}
			if (N1[i] > 0)
			{
				for (j=0; j<=kmax; j++)
				{
					for (k=0; k<=kmax; k++)
					{
						pop1[s1m][i][j][k] = pop1[s1m2][i][j][k];
						pop1[s1m][i][j][k] -= pop1[s1m2][i][j][k] * mig;
						pop1[s1m][0][j][k] += pop1[s1m2][i][j][k] * mig;
					}
				}
			}
		}
	}
	// transmission to hosts from migration pool
	for(i=1; i<= host_num; i++)
	{
		newN2 = 0; //reset N2 to 0
		for (j=0; j<=kmax; j++)
		{
			for (k=0; k<=kmax; k++)
			{
				if (host_num > 1)
				{
					pop2[s2m][i][j][k] += poidev(pop2[s2m][0][j][k]/host_num*tr,seed);
				}
				newN2 += pop2[s2m][i][j][k];
			}
		}
		N2[i] = newN2;
		newN1 = 0;
		for (j=0; j<=kmax; j++)
		{
			for (k=0; k<=kmax; k++)
			{
				if (host_num > 1)
				{
					pop1[s1m][i][j][k] += poidev(pop1[s1m][0][j][k]/host_num*tr,seed);
				}
				newN1 += pop1[s1m][i][j][k];

			}
		}
		N1[i] = newN1;
		N[i] = N2[i] + N1[i];
		//immunization process
		if (N[i] >= 0.5*K0 && kf[i-1] > 0)
		{
			kthr[i-1] = 0;
		}

		if (kthr[i-1] == 0) // decrease only if kthr[i] value is greater than 0.
		{
			if (K[i-1] > 0.00000001) // carrying capacity has to be greater than 0 to decrease
			{
				K[i-1] -= kf[i-1];
				//printf("kf=%.3f\n",kf);
				//printf("K[%d]inside=%.2f\n",i-1,K[i-1]);
			}
			if (K[i-1] <= 0) // host carrying capacity goes to 0.
			{
				K[i-1] = 0.00000001;
			}
			if (N[i] == 0) // new host setup
			{
				kf[i-1] = kf0*ran1(seed);
				K[i-1] = (double) K0;
				kthr[i-1] = 1;
			}
		}
	}
	//double Kmean=0, Kvar=0;
	//int Ncount = 0;
	//if (gen == 200 || gen == 300 || gen == 400)
	/*{
		for (i=0; i< host_num; i++)
		{
			/
			if (N[i] > 0)
			{
				Ncount++;
			}
			if (gen == 200 || gen == 300 || gen == 400)
			{
				printf("K[%d]=%.0f\n",i,K[i]);
			}
			//Kmean += K[i];
		}
		if (Ncount == host_num && gen < 70)
		{
			printf("all hosts infected at gen=%d\n",gen);
		}

		Kmean = Kmean/host_num;
		for (i=0; i< host_num; i++)
		{
			Kvar += pow((K[i] - Kmean),2);
		}
		Kvar = Kvar/host_num;
		printf("Kmean = %.2f, Kvar = %.2f\n",Kmean,Kvar);
	}*/


	N2[0] = Nsum(host_num,N2);
	N1[0] = Nsum(host_num,N1);
	
	// erase migration info (** it's erased after record since 12/17/2018)

	if(MIGPRINTF)
	{
		printf("migration pool\n");
	}
	for (j=0; j<=kmax; j++)
	{
		for (k=0; k<=kmax; k++)
		{
			if(MIGPRINTF)
			{
				//if(POPDIST)
				{
					if(pop1[s1m][0][j][k]>0)
					{
						printf("pop1[%d][0][%d][%d]=%.2f\n",s1m,j,k,pop1[s1m][0][j][k]);
					}
				}
			}
			//pop1[s1m][0][j][k] = 0;
		}
	}
	for (j=0; j<=kmax; j++)
	{
		for (k=0; k<=kmax; k++)
		{
			if(MIGPRINTF)
			{
				//if (POPDIST)
				{
					if(pop2[s2m][0][j][k]>0)
					{
						printf("pop2[%d][0][%d][%d]=%.2f\n",s2m,j,k,pop2[s2m][0][j][k]);
					}
				}
			}
			//pop2[s2m][0][j][k] = 0;
		}
	}

	if(MIGPRINTF)
	{
		count1 = 0;
		count2 = 0;
		printf("after migr inside the function\n");
	}
	// erase m2 info in order to make mutate function work on the next generation.
	for (i=1; i<=host_num; i++)
	{
		for (j=0; j<=kmax; j++)
		{
			for (k=0; k<=kmax; k++)
			{
				pop2[s2m2][i][j][k] = 0;
				pop1[s1m2][i][j][k] = 0;
				//if (pop2[m][i][j][k] >0 )
				{
					if (MIGPRINTF)
					{
						if(pop2[s2m][i][j][k] > 0)
						{
							if(POPDIST)
							{
								printf("pop2[%d][%d][%d][%d]=%.2f\n",s2m,i,j,k,pop2[s2m][i][j][k]);
							}
							count2 += pop2[s2m][i][j][k];
						}
						if(pop1[s1m][i][j][k] > 0)
						{
							if (POPDIST)
							{
								printf("pop1[%d][%d][%d][%d]=%.2f\n",s1m,i,j,k,pop1[s1m][i][j][k]);
							}
							count1 += pop1[s1m][i][j][k];
						}
					}
				}
			}
		}
	}
	if(MIGPRINTF)
	{
		printf("count1=%.2f\n",count1);
		printf("count2=%.2f\n",count2);
	}
}

void evol(double**** pop2, double**** pop1, int* curpop2, int* curpop1, int kmax, int host_num, double* N2, double* N1, double* N, long* seed, double evolrate)
{
	// every 1seg individual have chance of evolrate of turning into 2 seg.
	// of times in a poisson distribution with lambda = evolrate
	int s2,s1,i,j,k;
	double count1, count2;
	s2 = *curpop2;
	s1 = *curpop1;	
	if (EVOPRINTF)
	{
		count1 = 0;
		count2 = 0;
		printf("before evol inside the function\n");
		for (i=1; i<=host_num; i++)
		{
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					if (pop2[s2][i][j][k] > 0)
					{
						if(POPDIST)
						{
							printf("pop2[%d][%d][%d][%d]=%.2f\n",s2,i,j,k,pop2[s2][i][j][k]);
						}
						count2 += pop2[s2][i][j][k];
					}
				}
			}
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					if (pop1[s1][i][j][k] > 0)
					{
						if(POPDIST)
						{
							printf("pop1[%d][%d][%d][%d]=%.2f\n",s1,i,j,k,pop1[s1][i][j][k]);
						}
						count1 += pop1[s1][i][j][k];
					}
				}
			}
		}
		printf("count2=%.3f\n",count2);
		printf("count1=%.3f\n",count1);
	}
	double p; // amount added to pop2 and subtracted from pop1
	for (i=1; i<=host_num; i++)
	{
		if (N1[i] > 0)
		{
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					//p = poidev(evolrate*pop1[s1][i][j][k],seed); // stochastic version
					p = evolrate*pop1[s1][i][j][k]; // amount of seg1 virus turned into 2 seg in class (j,k)
					pop2[s2][i][j][k] += p;
					N2[i] += p;
					N2[0] += p;
					pop1[s1][i][j][k] -= p;
					N1[i] -= p;
					N1[0] -= p;
				}
			}
		}
	}
	if (EVOPRINTF)
	{
		count1 = 0;
		count2 = 0;
		printf("after evol inside the function\n");
		for (i=1; i<=host_num; i++)
		{
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					if (pop2[s2][i][j][k] > 0)
					{
						if(POPDIST)
						{
							printf("pop2[%d][%d][%d][%d]=%.2f\n",s2,i,j,k,pop2[s2][i][j][k]);
						}
						count2 += pop2[s2][i][j][k];
					}
				}
			}
			for (j=0; j<=kmax; j++)
			{
				for (k=0; k<=kmax; k++)
				{
					if (pop1[s1][i][j][k] > 0)
					{
						if(POPDIST)
						{
							printf("pop1[%d][%d][%d][%d]=%.2f\n",s1,i,j,k,pop1[s1][i][j][k]);
						}							
						count1 += pop1[s1][i][j][k];
					}
				}
			}
		}
		printf("count2=%.3f\n",count2);
		printf("count1=%.3f\n",count1);
	}
}


double Nsum(int size,double N[])
{
  int i;
  double sum = 0;
  for(i=1; i<=size; i++){
    sum += N[i];
  }
  return sum;
}

void record(double* N1, double* N2, double**** pop1, double**** pop2, int kmax, int host_num, int timestep, int krecord, int curpop1, int curpop2, int rep, int gen, FILE **fPointer)
{
	int i,j,k;
	double krecord1, krecord2;
	double krecord1t, krecord2t;
	double p1hist[2*kmax+1], p2hist[2*kmax+1];
	
	for (j=0; j<=2*kmax; j++) // pop1
	{
		p2hist[j] = 0;
		p1hist[j] = 0;
	}
	char* str;

	if (krecord == 0 || krecord == 1)
	{
		str = (char*) malloc(sizeof(char)*(10+10+2*2*20*(host_num+1)));
	}
	else if (krecord == 2)
	{
		str = (char*) malloc(sizeof(char)*((2*kmax+1)*2*20 + 50)); //(number of mutations)*(1&2seg)*(digits for pop size + commas) + buffer
 	}
	sprintf(str,"");
	if(krecord == 0) // record mean
	{
		krecord1t = 0; // record of entire population for 1 seg
		krecord2t = 0; // `` for 2 seg
		for (i=1; i<=host_num; i++)
		{
			krecord1 = 0; // record for 1 host for 1seg
			krecord2 = 0; // for 2 seg
			if (N2[i] > 0)
			{
				for (j=0; j<=kmax; j++)
				{
					for (k=0; k<=kmax; k++)
					{
						krecord2 += pop2[curpop2][i][j][k]/N2[i] * (j + k);
					}
				}
				krecord2t += krecord2*N2[i]/N2[0];
			}
			else
			{
				krecord2 = -1.0;
			}
			if (N1[i] > 0)
			{
				for (j=0; j<=kmax; j++)
				{
					for (k=0; k<=kmax; k++)
					{
						krecord1 += pop1[curpop1][i][j][k]/N1[i] * (j + k);
					}
				}
				krecord1t += krecord1*N1[i]/N1[0];
			}
			else
			{
				krecord1 = -1.0;
			}
			sprintf(str,"%s,%.2f,%.2f,%.2f,%.2f",str,N1[i],N2[i],krecord1,krecord2);
		}
		i=0;
		krecord1 = 0;
		krecord2 = 0;
		float migpop1 = 0;
		float migpop2 = 0;
		for (j=0; j<=kmax; j++)
		{
			for (k=0; k<=kmax; k++)
			{
				krecord2 += pop2[curpop2][i][j][k] * (j + k);
				migpop2 += pop2[curpop2][i][j][k];
			}
		}
		for (j=0; j<=kmax; j++)
		{
			for (k=0; k<=kmax; k++)
			{
				krecord1 += pop1[curpop1][i][j][k] * (j + k);
				migpop1 += pop1[curpop1][i][j][k];
			}
		}
		krecord1 = krecord1/migpop1;
		krecord2 = krecord2/migpop2;
		sprintf(str,"%s,%.2f,%.2f,%.2f,%.2f",str,migpop1,migpop2,krecord1,krecord2);
	}
	else if(krecord == 1)// record minimum
	{
		krecord2 = kmax*2 + 1;
		krecord1 = kmax*2 + 1;
		for (i=1; i<=host_num; i++)
		{
			if (N2[i] > 0)
			{
				for (j=0; j<=kmax; j++)
				{
					for (k=0; k<=kmax; k++)
					{
						if((j + k) < krecord2)
						{
							if(pop2[curpop2][i][j][k] > 0)
							{
								krecord2 = j + k;
							}
						}
					}
				}
				if (krecord2t > krecord2)
				{
					krecord2t = krecord2;
				}
			}
			else
			{
				krecord2 = -1.0;
			}
			if (N1[i] > 0)
			{
				for (j=0; j<=kmax; j++)
				{
					for (k=0; k<=kmax; k++)
					{
						if((j + k) < krecord1)
						{
							if(pop1[curpop1][i][j][k] > 0)
							{
								krecord1 = j + k;
							}
						}
					}
				}
				if (krecord1t > krecord1)
				{
					krecord1t = krecord1;
				}
			}
			else
			{
				krecord1 = -1.0;
			}
			sprintf(str,"%s,%.2f,%.2f,%.2f,%.2f",str,N1[i],N2[i],krecord1,krecord2);
		}
	}
	else // histogram (krecord == 2)
	{
		for (i=1; i<=host_num; i++)
		{
			for (j=0; j<=kmax; j++) // pop2
			{
				for (k=0; k<=kmax; k++)
				{
					p2hist[(j+k)] += pop2[curpop2][i][j][k];
				}
			}
			for (j=0; j<=kmax; j++) // pop1
			{
				for (k=0; k<=kmax; k++)
				{
					p1hist[(j+k)] += pop1[curpop1][i][j][k];
				}
			}
		}
	}
	if(timestep == 1)
	{

		if (N1[0] == 0)
		{
			krecord1t = -1.0;
		}
		if (N2[0] == 0)
		{
			krecord2t = -1.0;
		}
		if (krecord != 2) // not histogram
		{
			if (SHORTCSV)
			{
				fprintf(*fPointer,"%d,%d,%.2f,%.2f,%.8f,%.8f\n",rep+1,gen+1,N1[0],N2[0],krecord1t, krecord2t);
			}
			else
			{
				fprintf(*fPointer,"%d,%d,%.2f,%.2f,%.8f,%.8f%s\n",rep+1,gen+1,N1[0],N2[0],krecord1t, krecord2t, str);
			}
		}
		else // histogram
		{
			for (j=0; j<=2*kmax; j++)
			{
				sprintf(str,"%s,%.8f,%.8f",str,p1hist[j],p2hist[j]);
			}
			fprintf(*fPointer,"%d,%d,%.2f,%.2f%s\n",rep+1,gen+1,N1[0],N2[0], str);
		}
	}
	else
	{
		if(N1[0] == 0)
		{
			krecord1t = -1.0;
		}
		if(N2[0] == 0)
		{
			krecord2t = -1.0;
		}
		if (krecord != 2)
		{
			if (SHORTCSV)
			{
				fprintf(*fPointer,"%d,%d,%.2f,%.2f,%.8f,%.8f\n",rep+1,gen+1,N1[0],N2[0],krecord1t, krecord2t);
			}
			else
			{
				fprintf(*fPointer,"%d,%d,%.2f,%.2f,%.8f,%.8f%s\n",rep+1,gen+1,N1[0],N2[0],krecord1t, krecord2t, str);
			}
		}
		else
		{
			for (j=0; j<=2*kmax; j++)
			{
				sprintf(str,"%s,%.8f,%.8f",str,p1hist[j],p2hist[j]);
			}
			fprintf(*fPointer,"%d,%.2f,%.2f%s\n",rep+1,N1[0],N2[0], str);
		}
	}
	free(str); 
}

// random number generating functions ran1=uniform [0,1], bnldev= binomial. gammln is needed for bnldev.
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)	
float ran1(long *idum)
{
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

float gammln(float xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

#define PI 3.141592654

float poidev(float xm,long *idum)
{
	float gammln(float xx);
	float ran1(long *idum);
	static float sq,alxm,g,oldm=(-1.0);
	float em,t,y;

	if (xm < 12.0) {
		if (xm != oldm) {
			oldm=xm;
			g=exp(-xm);
		}
		em = -1;
		t=1.0;
		do {
			++em;
			t *= ran1(idum);
		} while (t > g);
	} else {
		if (xm != oldm) {
			oldm=xm;
			sq=sqrt(2.0*xm);
			alxm=log(xm);
			g=xm*alxm-gammln(xm+1.0);
		}
		do {
			do {
				y=tan(PI*ran1(idum));
				em=sq*y+xm;
			} while (em < 0.0);
			em=floor(em);
			t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
		} while (ran1(idum) > t);
	}
	return em;
}

#undef PI

double fact(int num)
{
	// factorial
	double val = 1;
	int i;
	for (i=1; i<=num; i++)
	{
		val *= i;
	}
	return val;
}

double C(int n, int k)
{
	// combination
	double val = fact(n)/(fact(k)*fact(n-k));
	return val;
}

double poipmf(double l, int k)
{
	//pmf function of poisson
	double val = (pow(l,k)*exp(-1*l))/fact(k);
	return val;
}

float gasdev(idum)
long *idum;
{
	static int iset=0;
	static float gset;
	float fac,r,v1,v2;
	float ran1();

	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			r=v1*v1+v2*v2;
		} while (r >= 1.0 || r == 0.0);
		fac=sqrt(-2.0*log(r)/r);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}