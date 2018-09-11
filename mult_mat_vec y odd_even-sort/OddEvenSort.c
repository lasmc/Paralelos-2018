#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#define tam 3200

void PrintVector(int a[], int n)
{
	printf("El vector contiene:\n");
	for (int i=0;i<n;i++)
		printf("%7d ",a[i]);
}

int cmp(const void* ap, const void* bp) {
  int a = * ((const int*) ap);
  int b = * ((const int*) bp);

  if (a < b) 
  {
    return -1;
  } 
  else if (a > b) 
  {
    return 1;
  } 
  else 
  {
    return 0;
  }
}

int max_index(int* data) 
{
  int i, max = data[0], maxi = 0;
  for (i = 1; i < tam; i++) 
  {
    if (data[i] > max) 
    {
      max = data[i];
      maxi = i;
    }
  }
  return maxi;
}

int min_index(int* data) 
{
  int i, min = data[0], mini = 0;
  for (i = 1; i < tam; i++) 
  {
    if (data[i] < min) 
    {
      min = data[i];
      mini = i;
    }
  }
  return mini;
}

void parallel_sort(int* data, int rank, int size) {
  int i;
  int local_vec[tam];
  for (i = 0; i < size; i++) {
	qsort(data, tam, sizeof(int), &cmp);
    int part;
    if (i % 2 == 0) {
      if (rank % 2 == 0) {
        part = rank + 1;
      } else {
        part = rank - 1;
      }
    } else {
      if (rank % 2 == 0) {
        part = rank - 1;
      } else {
        part = rank + 1;
      }
    }

    if (part < 0 || part >= size) {
      continue;
    }
    if (rank % 2 == 0) {
      MPI_Send(data, tam, MPI_INT, part, 0, MPI_COMM_WORLD);
      MPI_Recv(local_vec, tam, MPI_INT, part, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    } else {
      MPI_Recv(local_vec, tam, MPI_INT, part, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Send(data, tam, MPI_INT, part, 0, MPI_COMM_WORLD);
    }
    if (rank < part) {
      while (1) {
        int mini = min_index(local_vec);
        int maxi = max_index(data);
        if (local_vec[mini] < data[maxi]) {
          int temp = local_vec[mini];
          local_vec[mini] = data[maxi];
          data[maxi] = temp;
        } else {
          break;
        }
      }
    } else {
      while (1) {
        int maxi = max_index(local_vec);
        int mini = min_index(data);
        if (local_vec[maxi] > data[mini]) {
          int temp = local_vec[maxi];
          local_vec[maxi] = data[mini];
          data[mini] = temp;
        } else {
          break;
        }
      }
    }
  }
}

int main(int argc, char** argv) 
{
	int * vec;
	int comn_sz, rank;
  double startwtime, endwtime;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comn_sz);
  vec=malloc(tam*sizeof(int));

  startwtime = MPI_Wtime();
  for (int i=0;i<tam;i++)
   	*(vec + i)=1+(int)(10*random()/(RAND_MAX+1.0));
  //PrintVector(vec,tam);
  parallel_sort(vec, rank, comn_sz);
  //PrintVector(vec,tam);
  endwtime = MPI_Wtime();
  MPI_Finalize();
  printf("wall clock time = %f\n",(endwtime - startwtime)*1000);
  return 0;
}