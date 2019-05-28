/*
 Influenza Competition model 1.3
 competition model btw 1seg and 2seg in a WF model
 Some key differences from 1.2's:
 - WF model (no demographic factor or carrying capcacity for that matter.)
 Reproduction process:
 1. pick 2 random parents
 2. If either one is 1segment, only one of them gets replicated.
 3. If both of them are 2segement, it goes through recombination with 
    a probability of r, and one of the progeny is created.
 4. the created progeny from either proces  2. or 3. are assessed 
  whether they will survive by their fitness factor w.
 6. If a random number ([0,1]) exceeds w, the progeny dies and 
    has to go through the process 2.-5. again.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
/*
 Parameters
 back = whether the simulation allows back mutation (1 if there's one. 0 if there isn't)
 timestep = whether a simulation wants to record every generation or only the end generation per repetition
 krecord = how you want to record k in the output csv. 0= mean k. 1=all k values of indivs in an array. 2= minimum k.
 N0 = Population size
 K = Carrying capacity
 L = sequence length
 s = fitness decrease from deleterious mutation
 mu = mutation rate per site
 gen_num = generation amount 
 r = reassortment rate
 rep = repetition amount
 N1r = ratio of 1segment virus
*/
// virus basic structure
struct virus {
	int id;
	int k1;
	int k2;
	int k;
};

//https://stackoverflow.com/questions/2620146/how-do-i-return-multiple-values-from-a-function-in-c
float ran1(long *seed);
float gammln(float xx);
float bnldev(float pp, int n, long *idum);
void mutate(long *seed, int back, int N0, double mu, int L, struct virus popop[]);
struct virus *step(long *seed, int rep, int t, double cost, int N0, int L, int timestep, int krecord, double s, int K, double mu, double r,struct virus popop[],struct virus *next_gen_p,FILE **fPointer, int* N1, int* N2, int gen, int gen_num, double q, double a, double b, int type);
int intmin(int argc,int array[]); //min value of an integer array
int intsum(int size,int a[]);
float survp(int type, double s, double q, double a, double b, int k, int hascost, double cost);

int main(int argc, char *argv[]) {
	// set progress bar and initiate timer
	clock_t begin = clock();

	// call in parameters

	char *description = argv[1];
	char *destination = argv[3];
	char *back_s = argv[4];
	char *timestep_s = argv[5];
	char *krecord_s = argv[6];
	char *rep_s = argv[7];
	char *L_s = argv[8];
	char *s_s = argv[9];
	char *N0_s = argv[10];
	char *K_s = argv[11];
	char *mu_s = argv[12];
	char *gen_num_s = argv[13];
	char *cost_s = argv[14];
	char *r_s = argv[15];
	char *N1r_s = argv[16];
	char *seed_s = argv[17];
	char *untilext_s = argv[18];
	char *q_s = argv[19];
	char *a_s = argv[20];
	char *b_s = argv[21];
	char *type_s = argv[22]; // type of mutation interactions. 0 = independent, 1 = intermediate eqn1, 2 = intermediate eqn2
	char *end1;

	int back = (int) strtol(back_s,&end1,10);
	int timestep = (int) strtol(timestep_s,&end1,10);
	int krecord = (int) strtol(krecord_s, &end1, 10);
	int rep = (int) strtol(rep_s,&end1,10);
	int L = (int) strtol(L_s,&end1,10);
	double s = (double) strtof(s_s, NULL);
	int N0 = (int) strtol(N0_s,&end1,10);
	int K = (int) strtol(K_s,&end1,10);
	double mu = (double) strtof(mu_s,NULL);
	int gen_num = (int) strtol(gen_num_s,&end1,10);
	double cost = (double) strtof(cost_s,NULL);
	double r = (double) strtof(r_s,NULL);
	double N1r = (double) strtof(N1r_s,NULL);
	long seed = strtol(seed_s,&end1,10);
	int untilext = (int) strtol(untilext_s,&end1,10);
	double q = (double) strtof(q_s,NULL);
	double a = (double) strtof(a_s,NULL);
	double b = (double) strtof(b_s,NULL);
	int type = (int) strtol(type_s, &end1,10);

	printf("destination=%s, type=%d, back=%d, timestep=%d, untilext=%d, krecord=%d, rep=%d, L=%d, s=%.2f, N0=%d, K=%d, mu=%.5f, gen_num=%d, cost=%.2f, r=%.2f, N1r=%.2f, q=%.2f, a=%.2f, b=%.2f\n", destination, type, back, timestep, untilext, krecord, rep, L, s, N0, K, mu, gen_num, cost, r, N1r, q, a, b);

	//initiate csv file
	//// set up folder
	char *dest2 = (char*) malloc(50*sizeof(char)); 
	sprintf(dest2, "./data/%s", destination);
	struct stat st = {0};
	if (stat(dest2, &st) == -1) { // if the destination folder doesn't exist, make one with that name.
    	mkdir(dest2, 0700);
	}
	//// initiate file
	char* paraminfo = (char*) malloc(sizeof(char)*1000);
	sprintf(paraminfo,"%d,%.3f,%.3f,%.3f,%d,%d,%d,%.2f,%d,%d,%.5f,%d,%.7f,%.2f,%.2f,%d,%d",type,q,a,b,back,rep,L,s,N0,K,mu,gen_num,cost,r,N1r,timestep,krecord);
	char *filename = (char*) malloc(1000*sizeof(char));
	sprintf(filename,"%s/c1.3s_%s(0).csv",dest2,paraminfo);
	//sprintf(filename,"%s/c1.3s_%d,%d,%d,%d,%.2f,%d,%d,%.5f,%d,%.2f,%.2f,%.2f(0).csv",dest2,type,back,rep,L,s,N0,K,mu,gen_num,cost,r,N1r);
	int filenum  = 0;
	while( access( filename, F_OK ) != -1 ) { // check if file exists and change the file number if it exists
	    filenum += 1;
		sprintf(filename,"%s/c1.2s_%s(%d).csv",dest2,paraminfo,filenum);
	}
	FILE * fPointer;
	fPointer = fopen(filename,"w");
	free(dest2);
	free(filename);
	fprintf(fPointer,"#%s\n",description); // putting in comments on the csv
	fprintf(fPointer,"#%s\n",paraminfo);

	// simulation start
	struct virus pop[N0];
	int N1 = N0*N1r; // initial 1seg pop
	int N2 = N0*(1-N1r); // initial 2seg pop
	struct virus next_gen[N0];
	struct virus *pop2;

	if (timestep && krecord != 1) 
	{
		fprintf(fPointer,"rep,t,pop1,pop2,k1,k2\n");
	} 
	else if (timestep && krecord == 1)
	{
		char *header1 = (char*) malloc(650*sizeof(char));
		sprintf(header1, "rep,t,pop1,pop2,");
		int m;
		for( m=0; m<51; m++) 
		{
			sprintf(header1,"%s,k1.%d",header1,m);			
		}
		for( m=0; m<51; m++)
		{
			sprintf(header1,"%s,k2.%d",header1,m);			
		}
		sprintf(header1,"%s\n",header1);
		printf("%s",header1);
		free(header1);
	}
	else
	{
		fprintf(fPointer,"pop1,pop2,k1,k2\n");
	}
	int i,repe,gen;
	for (repe=0; repe<rep; repe++)
	{	
		N1 = N0*N1r; // initial 1seg pop
		N2 = N0*(1-N1r); // initial 2seg pop
		// initialize pop
		for (i=0;i<N1;i++)
		{
			pop[i].id = 1;
			pop[i].k1 = 0;
			pop[i].k2 = 0;
			pop[i].k = 0;
		}
		for (i=0;i<N2;i++)
		{
			pop[i+N1].id = 2;
			pop[i+N1].k1 = 0;
			pop[i+N1].k2 = 0;
			pop[i+N1].k = 0;
		}
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
		for (gen=0; gen<gen_num; gen++)
		{
			// run through generation
			//printf("GEN=%d/%d\n",gen,gen_num);
			if (untilext == 1 )
			{
				if (N1 == 0 || N2 == 0)
				{
					break;	
				}
				
			}
			mutate(&seed,back,N0,mu,L,pop);
			pop2 = step(&seed,(repe+1),(gen+1),cost,N0,L,timestep,krecord,s,K,mu,r,pop,next_gen,&fPointer,&N1,&N2,gen,gen_num,q,a,b,type);
			memcpy(pop,pop2,sizeof(struct virus)*N0); // cycle between pop and pop2 to continue looping.
		}
	}

	// close file and timer
	fclose(fPointer);
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("time spent was %.2f minutes\n", time_spent/60.0);

	return 0;
}


struct virus *step(long *seed,int rep, int t, double cost, int N0, int L, int timestep, int krecord, double s, int K, double mu, double r,struct virus popop[],struct virus *next_gen_p,FILE **fPointer,int* N1, int* N2, int gen, int gen_num, double q, double a, double b, int type) {
	// goes through reproduction process
	// the process is depicted at the top of the script.
	// input: pop struct array
	// output: next generation's pop struct array
	int l = 0; // next gen length
	int ks1[N0];
	int ks1l = 0;
	int ks2[N0];
	int ks2l = 0;
	int s1, s2;
	while (l < N0) 
	{
		s1 = floor(ran1(seed)*N0); // sample 1
		s2 = floor(ran1(seed)*N0); // sample 2
		//printf("s1:%d k1:%d k2:%d k:%d\n",s1,popop[s1].k1,popop[s1].k2,popop[s1].k);
		//printf("s2:%d k1:%d k2:%d k:%d\n",s2,popop[s2].k1,popop[s2].k2,popop[s2].k);
		if (popop[s1].id == 1 || popop[s2].id == 1) 
		{ 
			// either parents is segment 1
			if (ran1(seed) < 0.5) 
			{	
				// pick s1
				if (ran1(seed) < survp(type,s,q,a,b,popop[s1].k,0,cost))
				{
					next_gen_p[l].id = popop[s1].id;
					next_gen_p[l].k1 = popop[s1].k1;
					next_gen_p[l].k2 = popop[s1].k2;
					next_gen_p[l].k = popop[s1].k;
					if (next_gen_p[l].id == 1) 
					{
						ks1[ks1l] = next_gen_p[l].k;
						ks1l += 1;
					} 
					else 
					{
						ks2[ks2l] = next_gen_p[l].k;
						ks2l += 1;
					}
					l += 1;
				}
			} 
			else 
			{ 
				// pick s2
				if (ran1(seed) < survp(type,s,q,a,b,popop[s2].k,0,cost))
				{
					next_gen_p[l].id = popop[s2].id;
					next_gen_p[l].k1 = popop[s2].k1;
					next_gen_p[l].k2 = popop[s2].k2;
					next_gen_p[l].k = popop[s2].k;
					if (next_gen_p[l].id == 1) 
					{
						ks1[ks1l] = next_gen_p[l].k;
						ks1l += 1;
					} 
					else 
					{
						ks2[ks2l] = next_gen_p[l].k;
						ks2l += 1;
					}					
					l += 1;
				}
			}
		}
		else
		{ // both parents segment 2
			if (ran1(seed)<r) 
			{ // recombination happens
				if (ran1(seed) < 0.5) 
				{ // pick k1 from s1 and k2 from s2
					if (ran1(seed) < survp(type,s,q,a,b,(popop[s1].k1 + popop[s2].k2),1,cost))
					{
						next_gen_p[l].id = popop[s1].id;
						next_gen_p[l].k1 = popop[s1].k1;
						next_gen_p[l].k2 = popop[s2].k2;
						next_gen_p[l].k = popop[s1].k1 + popop[s2].k2;
						ks2[ks2l] = next_gen_p[l].k;
						ks2l += 1;
						l += 1;
					}
				}
				else
				{ // pick k1 from s2 and k2 from s1
					if (ran1(seed) < survp(type,s,q,a,b,(popop[s1].k2 + popop[s2].k1),1,cost))
					{
						next_gen_p[l].id = popop[s2].id;
						next_gen_p[l].k1 = popop[s2].k1;
						next_gen_p[l].k2 = popop[s1].k2;
						next_gen_p[l].k = popop[s1].k2 + popop[s2].k1;
						ks2[ks2l] = next_gen_p[l].k;
						ks2l += 1;
						l += 1;
					}
				}
			}
			else 
			{ // recombintion doesn't happen
				if (ran1(seed) < 0.5) 
				{ // pick s1
					if (ran1(seed) < survp(type,s,q,a,b,popop[s1].k,1,cost))
					{
						next_gen_p[l].id = popop[s1].id;
						next_gen_p[l].k1 = popop[s1].k1;
						next_gen_p[l].k2 = popop[s1].k2;
						next_gen_p[l].k = popop[s1].k;
						ks2[ks2l] = next_gen_p[l].k;
						ks2l += 1;
						l += 1;
					}
				}
				else 
				{ // pick s2
					if (ran1(seed) < survp(type,s,q,a,b,popop[s2].k,1,cost))
					{
						next_gen_p[l].id = popop[s2].id;
						next_gen_p[l].k1 = popop[s2].k1;
						next_gen_p[l].k2 = popop[s2].k2;
						next_gen_p[l].k = popop[s2].k;
						ks2[ks2l] = next_gen_p[l].k;
						ks2l += 1;				
						l += 1;
					}
				}
			}
		}	
	}
	*N1 = ks1l;
	*N2 = ks2l;
	if (krecord == 0)
	{
		float mk1;
		float mk2;
		if (ks1l == 0) 
		{
			mk1 = -1;
		} 
		else 
		{
			mk1 = (float)intsum(ks1l,ks1)/ks1l;				
		}
		if (ks2l == 0)
		{
			mk2 = -1;
		} 
		else 
		{
			mk2 = (float)intsum(ks2l,ks2)/ks2l;
		}
		if (timestep)
		{
			fprintf(*fPointer,"%d,%d,%d,%d,%.2f,%.2f\n",rep,t,ks1l,ks2l,mk1,mk2);			
		}
		else if (gen == (gen_num -1) || ks1l == 0 || ks2l == 0)
		{
			fprintf(*fPointer,"%d,%d,%.2f,%.2f\n",ks1l,ks2l,mk1,mk2);
		}
	} 
	else if (krecord == 1) 
	{
		int k1list[51];
		int k2list[51];
		int i;
		for (i=0; i<ks1l; i++) 
		{
			if (ks1[i] > 50) 
			{
				printf("WARNINGWARNING THERE'S INDIVIDUAL WITH K HIGHER THAN 50!!!!\n");
			}
			k1list[ks1[i]] += 1; 
		}
		for (i=0; i<ks2l; i++)
		{
			if (ks1[i] > 50) 
			{
				printf("WARNINGWARNING THERE'S INDIVIDUAL WITH K HIGHER THAN 50!!!!\n");
			}
			k2list[ks2[i]] += 1;
		}
		char *line = (char*) malloc(10100*sizeof(char));
		sprintf(line,"%d,%d,%d,%d\n",rep,t,ks1l,ks2l);
		for (i=0; i<51; i++)
		{
			sprintf(line,"%s,%d",line,k1list[i]);
		}
		for (i=0; i<51; i++)
		{
			sprintf(line,"%s,%d",line,k2list[i]);
		}
		fprintf(*fPointer,"%s\n",line);
		free(line);
	}
	else
	{
		int kmin1;
		int kmin2;
		if (ks1l == 0) 
		{
			kmin1 = -1;
		}
		else 
		{
			kmin1 = intmin(ks1l,ks1);
		} 
		if (ks2l == 0) 
		{
			kmin2 = -1;
		} 
		else 
		{
			kmin2 = intmin(ks2l,ks2);
		}
		if (timestep) 
		{
			fprintf(*fPointer,"%d,%d,%d,%d,%d,%d\n",rep,t,ks1l,ks2l,kmin1,kmin2);
		}
		else if (gen == (gen_num -1) || ks1l == 0 || ks2l == 0)
		{

			fprintf(*fPointer,"%d,%d,%d,%d\n",ks1l,ks2l,kmin1,kmin2);
		}		
	}	
	return next_gen_p;
}


void mutate(long *seed, int back, int N0, double mu, int L, struct virus popop[]) {
	// goes through population and make every inidividual go through 
	// mutation process.
	// input: population struct array.
	// output: population struct with updated k values according to
	//   mutation rate.
	int i,j,mut_num,breakpt;
	for (i=0;i<N0; i++)
	{
		mut_num = bnldev(mu,L,seed); // binomial pick of number of mutation based on mu.
		if (back == 1)
		{ // back mutation
			if (popop[i].id == 1)
			{ // individual is 1segment
				for (j=0; j<mut_num; j++) 
				{ // add or subtract k
					if (ran1(seed)*L < popop[i].k) 
					{ // k/L is the probability of back mutating.
						popop[i].k -= 1;
					}
					else
					{
						popop[i].k += 1;
					}
				}
			} 
			else
			{ // individual is 2segment
				for (j=0; j<mut_num; j++) 
				{ // add or subtract k1 or k2, and k.
					if (ran1(seed)*L < popop[i].k) 
					{
						if (popop[i].k1==0) 
						{
							popop[i].k2 -= 1;
							popop[i].k -= 1;
						} 
						else if (popop[i].k2 == 0) 
						{
							popop[i].k1 -= 1;
							popop[i].k -= 1;
						} 
						else 
						{
							if (ran1(seed) < (float)popop[i].k1/popop[i].k) 
							{
								popop[i].k1 -= 1;
								popop[i].k -= 1;
							} 
							else 
							{
								popop[i].k2 -= 1;
								popop[i].k -= 1;
							}
						}
					}
					else 
					{
						if (ran1(seed) < 0.5) 
						{
							popop[i].k1 += 1;
							popop[i].k += 1;
						}
						else
						{
							popop[i].k2 += 1;
							popop[i].k += 1;
						}
					}
				}
			}
		}
		else 
		{ // no back mutation
			int dist = 1; // 0=uniform 1=binom
			if (popop[i].id == 1) 
			{
				popop[i].k += mut_num;
			}
			else 
			{
				//binom
				if (dist == 1) //binom
				{
					breakpt = bnldev(0.5, mut_num, seed);
				}
				else // uniform
				{
					breakpt = floor(ran1(seed)*(mut_num+1));
				}				
				popop[i].k1 += breakpt;
				popop[i].k2 += mut_num - breakpt;
				popop[i].k += mut_num;
			}
		}
	}
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
float bnldev(float pp, int n,long *idum)
{
	float gammln(float xx);
	float ran1(long *idum);
	int j;
	static int nold=(-1);
	float am,em,g,angle,p,bnl,sq,t,y;
	static float pold=(-1.0),pc,plog,pclog,en,oldg;

	p=(pp <= 0.5 ? pp : 1.0-pp);
	am=n*p;
	if (n < 25) {
		bnl=0.0;
		for (j=1;j<=n;j++)
			if (ran1(idum) < p) bnl += 1.0;
	} else if (am < 1.0) {
		g=exp(-am);
		t=1.0;
		for (j=0;j<=n;j++) {
			t *= ran1(idum);
			if (t < g) break;
		}
		bnl=(j <= n ? j : n);
	} else {
		if (n != nold) {
			en=n;
			oldg=gammln(en+1.0);
			nold=n;
		} if (p != pold) {
			pc=1.0-p;
			plog=log(p);
			pclog=log(pc);
			pold=p;
		}
		sq=sqrt(2.0*am*pc);
		do {
			do {
				angle=PI*ran1(idum);
				y=tan(angle);
				em=sq*y+am;
			} while (em < 0.0 || em >= (en+1.0));
			em=floor(em);
			t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)
				-gammln(en-em+1.0)+em*plog+(en-em)*pclog);
		} while (ran1(idum) > t);
		bnl=em;
	}
	if (p != pp) bnl=n-bnl;
	return bnl;
}
#undef PI

int intmin(int argc, int array[]){
  int min_num = array[0];
  int i;
  for(i=1;i<argc;i++)
  {
    if(array[i]<min_num)
    {
      min_num = array[i];
    }
  }
  return min_num;
}

int intsum(int size,int a[]){
  int i;
  int sum = 0;
  for(i=0; i<size; i++)
  {
    sum += a[i];
  }
  return sum;
}

float survp(int type, double s, double q, double a, double b, int k, int hascost, double cost)
{
	// survival probability calculation depending on the type
	float val;
	if (type == 0)
	{
		val = pow(1.0-s,k);
	}
	else if (type == 1)
	{
		val = pow(1.0-s,pow(k,q));
	}
	else
	{
		val = 1.0 - a*k - b*pow(k,2);
	}
	if (hascost)
	{
		val = val*(1.0-cost);
	}
	return val;
}
