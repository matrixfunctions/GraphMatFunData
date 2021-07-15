#include<stdio.h>
#include<time.h>
#include<math.h>
#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

LOC_INC

#define NOF_SAMPLES 10


int compare (const void * a, const void * b)
{
  return (*(double*)a > *(double*)b) ? 1 : (*(double*)a < *(double*)b) ? -1:0 ;
}

double median (double * values, int nof_values)
{
  qsort (values, nof_values, sizeof(double), compare);
  return values[nof_values/2];
}


int main (int argc, char *argv[]){
  double *A;
  double *E;
  double timing[NOF_SAMPLES];
  int n,i,j,k;
  //int n0;
  //clock_t tv1, tv2;
  //double time;

  long seconds,nanoseconds;
  double elapsed;
  struct timespec begin, end;

  //double timingv[NOF_SAMPLES];

  double sum=0;
  clock_gettime(CLOCK_REALTIME, &begin);

  n = MATSIZE;
  A = malloc(n * n * sizeof(*A));
  E = malloc(n * n * sizeof(*A));

  // START REPEATED CODE
  elapsed=0;
  sum=0;
  printf("* GRAPH: NAME\n    TIMINGS:  ");
  for (k=0;k<NOF_SAMPLES;k++){
    // Recreate the matrix every time since it is overwritten.
    for(i=0; i<n; i++){
      for(j=0; j<n; j++){
	A[i+n*j]=0.0;
	if (i!=j){ // Normalize at the same time
	  A[i+n*j]=1.0 / (123.57039475388511*sqrt(fabs(1.0*(i-j))));
	}
      }
    }
    A[n-1]=0.02236627204212922;

    // Run the computation with timing
    clock_gettime(CLOCK_REALTIME, &begin);
    dFUNCTION(A, n, E);
    clock_gettime(CLOCK_REALTIME, &end);

    seconds = end.tv_sec - begin.tv_sec;
    nanoseconds = end.tv_nsec - begin.tv_nsec;
    elapsed = seconds + nanoseconds*1e-9;

    timing[k]=elapsed;
    printf("%6.3f ",elapsed);
    sum += elapsed;
    sleep(1);
  }
  printf("\n");
  printf("    MEAN: %.3f, MEDIAN: %.3f\n",sum / (1.0*NOF_SAMPLES), median(timing,NOF_SAMPLES));

  sleep(2);  // Let CPU rest / let OS do mem handling
  // END REPEATED CODE

  return 0;

}
