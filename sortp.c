#include <time.h>		// for CPU time
#include <sys/time.h>	// for gettimeofday
#include <stdio.h>
#include <mpi.h>

#define max_len 100000
#define LENGTH 40

int main(int argc, char *argv[])
{
	// my variables
	int me, nproc, lenme;
	double bme[max_len+1], cme[max_len+1], low, high, dif, lowerBound, upperBound;
	float cpuTime12, cpuTime03, maxCpuTime12, maxCpuTime03;
	double maxTime12, maxTime03;
	//
	
	int i=1, len, ind[max_len+1], j, cur, prev;
	double b[max_len+1], c[max_len+1], new, cnew;
	char name[LENGTH] = "sort.txt", line[LENGTH];
	FILE *fp;
	clock_t cpu0, cpu1, cpu2, cpu3;				// clock_t defined in <time.h> and <sys/types.h> as int
	struct timeval time0, time1, time2, time3;	// for wall clock in s and us
	double dtime12, dtime03;					// for wall clock in s (real number)
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	
	// timing 0
	cpu0 = clock();
	gettimeofday(&time0, NULL);
	
	if (me == 0)
	{
		fp = fopen(name, "r");
		
		while (1)
		{
			if (fgets(line, LENGTH, fp) == NULL)
				break;
			
			if (sscanf(line, "%lf %lf", &b[i], &c[i]) == -1)
				break;
			
			i++;
		}
		
		len = i - 1;
		fclose(fp);
		printf("Number of items to sort: %i\n", len);
		
		low = b[0];
		high = b[0];
		for (i=0; i<len; i++)
		{
			if (b[i] < low)
				low = b[i];
			
			if (b[i] > high)
				high = b[i];
		}
		
		dif = high - low;
	}
	
	MPI_Bcast(&len,		1,			MPI_INT,	0, MPI_COMM_WORLD);
	MPI_Bcast(&b,		max_len+1, 	MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&c, 		max_len+1, 	MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&high,	1, 			MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&low, 	1, 			MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&dif, 	1,			MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	// timing 1
	cpu1 = clock();
	gettimeofday(&time1, NULL);
	
	// split numbers
	lowerBound = (double)((me * (dif / nproc)) + low);
	upperBound = (double)((me * (dif / nproc)) + dif / nproc + low);
	
	if (me == 0)
		lowerBound = low;
	
	if (me == nproc-1)
		upperBound = high;
	
	lenme = 1;
	
	if (me != 0)
	{
		for (j=1; j<=len; j++)
		{
			if ((lowerBound < b[j]) & (b[j] <= upperBound))
			{
				bme[lenme] = b[j];
				cme[lenme] = c[j];
				lenme++;
			}
		}
		
		lenme = lenme - 1;
	}
	else
	{
		for (j=1; j<=len; j++)
		{
			if ((lowerBound <= b[j]) & (b[j] <= upperBound))
			{
				bme[lenme] = b[j];
				cme[lenme] = c[j];
				lenme++;
			}
		}
		
		lenme = lenme - 1;
	}
	
	printf("Thread %i has %i numbers to sort\n", me, lenme);
	
	// sort numbers in each thread
	ind[0] = 1;
	for (j=2; j<=lenme; j++)	// start sorting with the second item
	{
		new = bme[j];
		cnew = cme[j];
		cur = 0;
		
		for (i=1; i<j; i++)
		{
			prev = cur;
			cur = ind[cur];
			
			if (new == bme[cur])
				printf("Equal numbers %lf\n", new);
			
			if ((new < bme[cur]) | ((new == bme[cur]) & (cnew < cme[cur])))
			{
				ind[prev] = j;
				ind[j] = cur;
				goto loop;
			}
		}
		
		// new number is the largest so far
		ind[cur] = j;
		loop: ;
	}
	
	// timing 2
	cpu2 = clock();
	gettimeofday(&time2, NULL);
	
	cpuTime12 = (float) (cpu2 - cpu1)/CLOCKS_PER_SEC;
	dtime12 = ((time2.tv_sec - time1.tv_sec) + (time2.tv_usec - time1.tv_usec) / 1e6);
	
	// write sorted numbers
	for (j=0; j<=nproc-1; j++)
	{
		if (me == j)
		{
			if (me == 0)
				fp = fopen("sortedp.txt", "w");
			else
				fp = fopen("sortedp.txt", "a");
			
			cur = 0;
			for (i=1; i<=lenme; i++)
			{
				cur = ind[cur];
				fprintf(fp, "%lf %lf\n", bme[cur], cme[cur]);
			}
			
			fclose(fp);
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
	}
	
	// timing 3
	cpu3 = clock();
	gettimeofday(&time3, NULL);
	
	cpuTime03 = (float) (cpu3 - cpu0)/CLOCKS_PER_SEC;
	dtime03 = ((time3.tv_sec - time0.tv_sec) + (time3.tv_usec - time0.tv_usec) / 1e6);
	
	// max timing
	MPI_Reduce(&cpuTime12,	&maxCpuTime12,	1, MPI_FLOAT,  MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&cpuTime03,	&maxCpuTime03,	1, MPI_FLOAT,  MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&dtime12,	&maxTime12,		1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&dtime03,	&maxTime03,		1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	
	if (me == 0)
	{
		printf("\nElapsed wall time sorting: %f \n",	maxTime12);
		printf("Elapsed CPU time sorting: %f \n",		maxCpuTime12);
		printf("Elapsed wall time complete: %f \n",		maxTime03);
		printf("Elapsed CPU time complete: %f \n\n",		maxCpuTime03);
	}
	
	MPI_Finalize();
}