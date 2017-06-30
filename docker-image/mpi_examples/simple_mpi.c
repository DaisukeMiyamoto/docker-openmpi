#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "mpi.h"

#define DATASIZE 10000

int main(int argc, char **argv)
{
  int myid;
  int p;
  int i,j;
  double ans[DATASIZE];
  double sum = 0.0;
  double s_time, e_time;

  MPI_Status status;
  char hostname[1024];
  hostname[1023] = '\0';
  gethostname(hostname, 1023);


  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  printf("Hello from %s (%d/%d)\n", hostname, myid, p);

  MPI_Barrier(MPI_COMM_WORLD);
  s_time = MPI_Wtime();

  for (i=0; i<DATASIZE; i++)
    {
      ans[i] = 0.0;
    }


  for (i=0; i<DATASIZE; i++)
    {
      for (j=0; j<DATASIZE; j++)
        {
          ans[i]  += exp(i* 0.00000005 + j*0.000002);
        }
    }

  for (i=0; i<DATASIZE; i++)
    {
      sum += ans[i];
    }

  MPI_Barrier(MPI_COMM_WORLD);
  e_time = MPI_Wtime();

  if (myid == 0)
    {
      printf ("time=%f\n", e_time - s_time);
    }

  MPI_Finalize();

}
