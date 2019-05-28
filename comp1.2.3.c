/*
Influenza Competition model 1.2
Same concept as comp model 1 but more efficient. 
Some key differences:
- back mutation, tracking info, and program cut strategy parameterized.
- No sequence information using array. Only keeps track of mutant allele amount
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
	int progeny;
	double w;
};

float ran1(long *seed);
float gammln(float xx);
float bnldev(float pp, int n, long *idum);
float poidev(float xm,long *idum);
void mutate(long *seed, int back, int N0, double mu, int L, struct virus popop[], int type, double q, double a, double b, double s);
struct virus *step(long *seed, int rep, int* N1, int* N2, int* N,int L, int timestep, int krecord,
	double s, int K, double mu, double r, struct virus pop[], struct virus *next_gen_p, FILE **fPointer, int gen_num, int untilext, int gen, double cost);
int *intsample(int *array, int l, int n, int replacement, long *seed, int *receive);
int *intdel(int *array, int l, int ran, int *receive);
int intmin(int argc,int array[]); //min value of an integer array
int intsum(int size,int a[]);
float survp(int type, double s, double q, double a, double b, int k);

int main(int argc, char *argv[]) {
	// set progress bar and initiate timer
	clock_t begin = clock();
	
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
	/*
	int back = 0;
	int timestep = 1;
	int krecord = 0;
	int untilext = 0;
	int rep = 1;
	int L = 300;
	double s = 0.05;
	int N0 = 10000;
	int K = 10000;
	double mu = 0.06;
	int gen_num = 2;
	double cost = 0.0;
	double r = 0.5;
	double N1r = 0.5;
	char *destination = "ctest";
	long seed = -2389;
	*/
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
	sprintf(filename,"%s/c1.2s_%s(0).csv",dest2,paraminfo);
	//sprintf(filename,"%s/c1.3s_%d,%d,%d,%d,%.2f,%d,%d,%.5f,%d,%.2f,%.2f,%.2f(0).csv",dest2,type,back,rep,L,s,N0,K,mu,gen_num,cost,r,N1r);
	int filenum  = 0;
	while( access( filename, F_OK ) != -1 ) { // check if file exists and change the file number if it exists
	    filenum += 1;
		sprintf(filename,"%s/c1.2s_%s(%d).csv",dest2,paraminfo,filenum);
		//sprintf(filename,"%s/c1.3s_%d,%d,%d,%d,%.2f,%d,%d,%.5f,%d,%.2f,%.2f,%.2f(%d).csv",dest2,type,back,rep,L,s,N0,K,mu,gen_num,cost,r,N1r,filenum);
	}
	FILE * fPointer;
	fPointer = fopen(filename,"w");
	free(dest2);
	free(filename);
	fprintf(fPointer,"#%s\n",description); // putting in comments on the csv
	fprintf(fPointer,"#%s\n",paraminfo);
	
	// simulation start
	
	if (timestep) {
		fprintf(fPointer,"rep,t,pop1,pop2,k1,k2\n");
	} else {
		fprintf(fPointer,"pop1,pop2,k1,k2\n");
	}
	
	// initiate some variables that will change over generations
	int N1; // seg1 population
	int N2; // seg2 population
	int N; // N1 + N2
	struct virus next_gen[K+K/2];
	struct virus *popb;
	struct virus *pop;
	int repe,i,gen;
	for (repe=0; repe<rep; repe++){	
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
		// initialize pop (generation 0)
		N1 = N0*N1r; // initial 1seg pop
		N2 = N0*(1-N1r); // initial 2seg pop
		N = N0;
		pop = (struct virus*) malloc(sizeof(struct virus)*N);
		
		for (i=0;i<N1;i++){
			pop[i].id = 1;
			pop[i].k1 = 0;
			pop[i].k2 = 0;
			pop[i].k = 0;
			pop[i].progeny = 0;
			pop[i].w = 1;
		}
		for (i=0;i<N2;i++){
			pop[i + N1].id = 2;
			pop[i + N1].k1 = 0;
			pop[i + N1].k2 = 0;
			pop[i + N1].k = 0;
			pop[i + N1].progeny = 0;
			pop[i + N1].w = 1;
		}
		
		// make first generation
		mutate(&seed,back,N,mu,L,pop,type,q,a,b,s); // seg1's mutate
		popb = step(&seed,(repe+1),&N1,&N2,&N,L,timestep,krecord,s,K,mu,r,pop,next_gen,&fPointer,gen_num,untilext,1,cost); // reproduction of seg1
		/*for (int i=0; i<N; i++) {
			printf("popb%d has k = %d\n",i,popb[i].k);
		}*/
		free(pop);

		for (gen=1; gen<gen_num; gen++){ // run through generation
			//printf("GEN=%d/%d\n",gen+1,gen_num);
			// cycle btw pop and popb to continue looping.
			pop = (struct virus*) malloc(sizeof(struct virus)*N);
			memcpy(pop,popb,sizeof(struct virus)*N); 
			/*for (int i=0; i<N; i++) {
				if (pop[i].k != popb[i].k){
					printf("not same(%d)\n",i);
				}
				//printf("pop%d has k = %d\n",i,popb[i].k);
			}*/
			mutate(&seed,back,N,mu,L,pop,type,q,a,b,s); // seg1's mutate
			if (N1 == 0 || N2 == 0)
			{
				if (untilext)
				{
					break;
				}
			}
			popb = step(&seed,(repe+1),&N1,&N2,&N,L,timestep,krecord,s,K,mu,r,pop,next_gen,&fPointer,gen_num,untilext,gen+1,cost); // reproduction of seg1
			free(pop);
		}
	}
	// close file and timer
	fclose(fPointer);
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("time spent was %.2f minutes\n", time_spent/60.0);
	return 0;
}

struct virus *step(long *seed, int rep, int* N1, int* N2, int* N,int L, int timestep, int krecord,
	double s, int K, double mu, double r, struct virus pop[], struct virus *next_gen_p, FILE **fPointer, int gen_num, int untilext, int gen, double cost) {
	// this function goes through reassortment and reproduction process for 2 segmented viruses. 
	// every reproduction from a pair of parents give 2 progenies.
	// Reproduction process:
	// 1. decide how many progeny each virus is having.
	// 2. Pick random 2 parents whose progeny number is not 0.
	// 3. Make 2 progenies from the two parents and deduct there progeny number by 1.
	// 4. If they are both 2 segs, with the probability of r(reassortment prob), the two progenies have reassorted genes
	//    (one progeny have 1st segment from parent 1 and 2nd segment from parent 2, and another have
	//    1st segment from parent 2 and 2nd from parent 1.)
	int remaining[*N]; // list of viruses who still have progenise to give.
	int remainingl = 0;
	int i;
	int l = 0; // length of next_gen_p


	for (i=0;i<*N;i++) // decide how many youngs each virus will have and fill the 'remaining' array.
	{
		if (pop[i].id == 1) 
		{
			pop[i].progeny = poidev(pop[i].w*(2/(1+(float)*N/K)),seed);
		}
		else
		{
			pop[i].progeny = poidev(pop[i].w*(2/(1+(float)*N/K)*(1-cost)),seed);
		}
		if (pop[i].progeny > 0)
		{
			remaining[remainingl] = i;
			remainingl++;
		}
	}

	//printf("remainings: ");
	//for (i=0; i<remainingl; i++) {
	//	printf("%d, ",remaining[i]);
	//}
	int *samp;
	int sampreceive[2];
	int temp;
	int *temparray;
	int remainreceive[remainingl];
	int sampindex; //samp index in remaining;
	float ks1 = 0;
	int ks1l = 0;
	float ks2 = 0;
	int ks2l = 0;	
	
	while (remainingl >= 2)
	{
		samp = intsample(remaining,remainingl, 2, 0, seed, sampreceive); // sample 2 parent viruses
		if (samp[0] < samp[1]) {
			temp = samp[0];
			samp[0] = samp[1];
			samp[1] = temp;
		}
		//printf("sample[0]=%d\n",samp[0]);	
		//printf("sample[1]=%d\n",samp[1]);
		if ((pop[samp[0]].id != pop[samp[1]].id) || (pop[samp[0]].id == 1 && pop[samp[1]].id == 1))
		{
			next_gen_p[l].id = pop[samp[0]].id; 
			next_gen_p[l].k1 = pop[samp[0]].k1; 
			next_gen_p[l].k2 = pop[samp[0]].k2;
			next_gen_p[l].k = pop[samp[0]].k;
			next_gen_p[l].progeny = 0;
			next_gen_p[l].w = pop[samp[0]].w; // supposed to calculate w newly, but doesn't matter cuz it will be changed after mutation
			l++;
			if (pop[samp[0]].id == 1)
			{
				if (krecord == 0)
				{
					ks1 += pop[samp[0]].k;
				}
				else if (krecord == 2)
				{
					if (pop[samp[0]].k < (int) ks1)
					{
						ks1 = pop[samp[0]].k;
					}
				}
				ks1l++;
			}
			else
			{
				if (krecord == 0)
				{
					ks2 += pop[samp[0]].k;
				}
				else if (krecord == 2)
				{
					if (pop[samp[0]].k < (int) ks2)
					{
						ks2 = pop[samp[0]].k;
					}
				}
				ks2l++;
			}
			next_gen_p[l].id = pop[samp[1]].id; 
			next_gen_p[l].k1 = pop[samp[1]].k1;
			next_gen_p[l].k2 = pop[samp[1]].k2;
			next_gen_p[l].k = pop[samp[1]].k;
			next_gen_p[l].progeny = 0;
			next_gen_p[l].w = pop[samp[1]].w;
			l++;
			if (pop[samp[0]].id == 1)
			{
				if (krecord == 0)
				{
					ks1 += pop[samp[0]].k;
				}
				else if (krecord == 2)
				{
					if (pop[samp[0]].k < (int) ks1)
					{
						ks1 = pop[samp[0]].k;
					}
				}
				ks1l++;
			}
			else
			{
				if (krecord == 0)
				{
					ks2 += pop[samp[0]].k;
				}
				else if (krecord == 2)
				{
					if (pop[samp[0]].k < (int) ks2)
					{
						ks2 = pop[samp[0]].k;
					}
				}
				ks2l++;
			}		
		}
		else
		{
			if (ran1(seed) > r) // no recombination
			{
				next_gen_p[l].id = pop[samp[0]].id; 
				next_gen_p[l].k1 = pop[samp[0]].k1; 
				next_gen_p[l].k2 = pop[samp[0]].k2;
				next_gen_p[l].k = pop[samp[0]].k;
				next_gen_p[l].progeny = 0;
				next_gen_p[l].w = pop[samp[0]].w; // supposed to calculate w newly, but doesn't matter cuz it will be changed after mutation
				l++;
				if (krecord == 0)
				{
					ks2 += pop[samp[0]].k;
				}
				else if (krecord == 2)
				{
					if (pop[samp[0]].k < (int) ks2)
					{
						ks2 = pop[samp[0]].k;
					}
				}
				ks2l++;
				next_gen_p[l].id = pop[samp[1]].id; 
				next_gen_p[l].k1 = pop[samp[1]].k1;
				next_gen_p[l].k2 = pop[samp[1]].k2;
				next_gen_p[l].k = pop[samp[1]].k;
				next_gen_p[l].progeny = 0;
				next_gen_p[l].w = pop[samp[1]].w;
				l++;
				if (krecord == 0)
				{
					ks2 += pop[samp[0]].k;
				}
				else if (krecord == 2)
				{
					if (pop[samp[0]].k < (int) ks2)
					{
						ks2 = pop[samp[0]].k;
					}
				}
				ks2l++;
			}
			else
			{
				// parent 1 give 1st segment
				next_gen_p[l].id = pop[samp[0]].id; 
				next_gen_p[l].k1 = pop[samp[0]].k1;
				next_gen_p[l].k2 = pop[samp[1]].k2;
				next_gen_p[l].k = next_gen_p[l].k1 + next_gen_p[l].k2;
				next_gen_p[l].progeny = 0;
				next_gen_p[l].w = pop[samp[1]].w; 
				if (krecord == 0)
				{
					ks2 += next_gen_p[l].k;
				}
				else if (krecord == 2)
				{
					if (next_gen_p[l].k < (int) ks2)
					{
						ks2 = next_gen_p[l].k;
					}
				}
				ks2l++;
				l++;
				// parent 2 give 1st segment
				next_gen_p[l].id = pop[samp[1]].id; 
				next_gen_p[l].k1 = pop[samp[1]].k1;
				next_gen_p[l].k2 = pop[samp[0]].k2;
				next_gen_p[l].k = next_gen_p[l].k1 + next_gen_p[l].k2;
				next_gen_p[l].progeny = 0;
				next_gen_p[l].w = pop[samp[1]].w;
				if (krecord == 0)
				{
					ks2 += next_gen_p[l].k;
				}
				else if (krecord == 2)
				{
					if (next_gen_p[l].k < (int) ks2)
					{
						ks2 = next_gen_p[l].k;
					}
				}
				ks2l++;
				l++;
			}
		}
		// deduct one progeny from the sampled individual's progeny number.
		pop[samp[1]].progeny -= 1;
		pop[samp[0]].progeny -= 1;
		// delete those with progeny num 0 from remaining list
		if (pop[samp[0]].progeny == 0)
		{
			for (i= (remainingl-1); i>-1; i--)
			{
				if (samp[0] == remaining[i])
				{
					sampindex = i;
				}

			}
			temparray = intdel(remaining,remainingl,sampindex,remainreceive);
			//printf("sample[0] %d is erased!\n",samp[0]);
			remainingl--;
			memcpy(remaining, temparray, sizeof(int)*remainingl);
		}
		if (pop[samp[1]].progeny == 0)
		{
			for (i=0; i<remainingl; i++)
			{
				if (samp[1] == remaining[i])
				{
					sampindex = i;
				}

			}
			temparray = intdel(remaining,remainingl,sampindex,remainreceive);
			//printf("sample[1] %d is erased!\n",samp[1]);
			remainingl--;
			memcpy(remaining, temparray, sizeof(int)*remainingl);
		}
		/*
		printf("remainings: ");
		for (i=0; i<remainingl; i++) {
			printf("%d, ",remaining[i]);
		}
		printf("\n");
		*/
	}

	if (l > K+K/2) 
	{
		printf("WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING\n");
		printf("WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING\n");
		printf("WARNING: NEXT_GEN ALLOCATION WAS NOT ENOUGH. ALLOCATE MORE MEMORIES TO RUN THIS SIMULATION\n");
		printf("WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING\n");
		printf("WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING\n");
	}

	if (remainingl > 0)
	{		
		for (i=0; i<pop[remaining[0]].progeny; i++) 
		{
			next_gen_p[l].id = pop[remaining[0]].id;
			next_gen_p[l].k1 = pop[remaining[0]].k1;
			next_gen_p[l].k2 = pop[remaining[0]].k2;
			next_gen_p[l].k = pop[remaining[0]].k;
			next_gen_p[l].progeny = 0;
			next_gen_p[l].w = pop[remaining[0]].w;
			if (pop[remaining[0]].id == 1)
			{
				if (krecord == 0)
				{
					ks1 += pop[remaining[0]].k;
				}
				else if (krecord == 2)
				{
					if (pop[samp[0]].k < (int) ks1)
					{
						ks1 = pop[remaining[0]].k;
					}
				}
			}
			else
			{
				if (krecord == 0)
				{
					ks2 += pop[remaining[0]].k;
				}
				else if (krecord == 2)
				{
					if (pop[samp[0]].k < (int) ks2)
					{
						ks2 = pop[remaining[0]].k;
					}
				}
			}
				
			l++;		
		}
	}

	/*for (i=0; i<*N2; i++) {
		printf("pop%d has k=%d\n",i,next_gen_p[i].k);
	}*/

	if (krecord == 0)
	{
		if (ks1l == 0)
		{
			ks1 = -1.0;
		}
		else
		{
			ks1 = ks1/ks1l;
		}
		if (ks2l == 0) 
		{
			ks2 = -1.0;
		} 
		else 
		{
			ks2 = ks2/ks2l;
		}
		if (timestep)
		{
			fprintf(*fPointer,"%d,%d,%d,%d,%.2f,%.2f\n",rep,gen,ks1l,ks2l,ks1,ks2);			
		}
		else if (gen == gen_num || ks1l == 0 || ks2l == 0)
		{
			fprintf(*fPointer,"%d,%d,%.2f,%.2f\n",ks1l,ks2l,ks1,ks2);
		}


	}
	else if (krecord == 1) 
	{
		// WILL FINISH WHEN NECESSARY.
		/*if (ks1l == 0){
			char *k1str = "'NA'"
		} else {
			char *k1str	= (char *) malloc(N0*4*sizeof(char));
		}
		if (ks2l == 0){
			char *k2str = "'NA'"
		} else {
			char *k2str = (char *) malloc(N0*4*sizeof(char));
			ks2strl = 0;
			for (int i=0; i<ks2l; i++) {
				char *int2str = (char *) malloc(4);
				sprintf(int2str,"%d",ks2[i]);
				k2str[i] 
			}
		}
		fprintf(fPointer,"",);*/
	}
	else 
	{
		if (ks1l == 0) {
			ks1 = -1.0;
		}
		if (ks2l == 0) {
			ks2 = -1.0;
		}
		if (timestep)
		{
			fprintf(*fPointer,"%d,%d,%d,%d,%.2f,%.2f\n",rep,gen,ks1l,ks2l,ks1,ks2);
		}
		else if (gen == gen_num || ks1l == 0 || ks2l == 0)
		{
			fprintf(*fPointer,"%d,%d,%.2f,%.2f\n",ks1l,ks2l,ks1,ks2);
		}

	}
	*N1 = ks1l;
	*N2 = ks2l;
	*N = l; // update N2	
	return next_gen_p;
}

int *intsample(int *array, int l, int n, int replacement,long *seed, int *receive) {
	// sampling from an integer array
	// array = list of integer to sapmle from
	// l = list length
	// n = sample size (don't put more than the length of list if no replacement)
	// replacement = true or false for replacement
	int i, ll, ran;
	ll = l;
	if (replacement)
	{
		for( i=0; i<n; i++)
		{
			ran = floor(ran1(seed)*ll);
			receive[i] = array[ran];
		}
	}
	else
	{
		int receive2[ll];
		for(i=0; i<n; i++)
		{
			ran = floor(ran1(seed)*ll);
			receive[i] = array[ran];
			array = intdel(array,ll,ran,receive2);
			ll -= 1;
		}
	}
	return receive;
}

int *intdel(int *array, int l, int ran, int *receive)
{
	// deletes an specified element in an array and returns the revised one.
	// array = an array where an element will be deleted
	// l = array length
	// ran = element index in an array

	int i;
	for (i=0; i<ran; i++)
	{
		receive[i] = array[i];
	}
	for (i=ran+1; i<l; i++)
	{
		receive[i-1] = array[i];
	}
	return receive;
}

void mutate(long *seed, int back, int N0, double mu, int L, struct virus popop[], int type, double q, double a, double b, double s) {

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
		popop[i].w = survp(type,s,q,a,b,popop[i].k);
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

int intmin(int argc, int array[]){
  int min_num = array[0];
  int i;
  for(i=1;i<argc;i++){
    if(array[i]<min_num){
      min_num = array[i];
    }
  }
  return min_num;
}

int intsum(int size,int a[]){
  int i;
  int sum = 0;
  for(i=0; i<size; i++){
    sum += a[i];
  }
  return sum;
}

float survp(int type, double s, double q, double a, double b, int k)
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
	return val;
}
