#include <iostream>
#include <time.h>

using namespace std;
__global__
void matMultKernel(float *d_M, float *d_N, float *d_P, int Width){
  int Row = blockIdx.y*blockDim.y + threadIdx.y;
  int Col = blockIdx.x*blockDim.x + threadIdx.x;
  int k = 0;
  if(Row < Width && Col < Width){
      float Pvalue = 0;
      for(k = 0; k < Width; ++k){
          Pvalue += d_M[Row*Width + k] * d_N[k*Width+Col];
      }
      d_P[Row*Width+Col] = Pvalue;
  }
}

void matMult(float* A, float* B, float* C, int n){
  int size = n*n*sizeof(float);
  float *d_A, *d_B, *d_C;

  cudaMalloc((void **) &d_A, size);
  cudaMemcpy(d_A,A,size,cudaMemcpyHostToDevice);
  cudaMalloc((void **) &d_B, size);
  cudaMemcpy(d_B,B,size,cudaMemcpyHostToDevice);
  cudaMalloc((void **) &d_C, size);

  dim3 dimGrid(ceil(n/1024.0),ceil(n/1024.0),1);
  dim3 dimBlock(1024,1024,1);  /*Dimensiónde 2 X 2, eje Z desactivado con 1*/
  matMultKernel<<<dimGrid, dimBlock>>>(d_A,d_B,d_C,n);
  
  cudaMemcpy(C,d_C,size,cudaMemcpyDeviceToHost);

  cudaFree(d_A); cudaFree(d_B); cudaFree(d_C);
}


int main(){
  int n,i,j;
  float *h_A,*h_B,*h_C;
  n = 10000;
  /******************Time Variables*********************************/
  clock_t time;
  cout<<"El Tam de la matriz Cuadrada es: "<<n<<" X "<<n<<endl;
  h_A = new float[n*n];
  h_B = new float[n*n];
  h_C = new float[n*n];

  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++)
      h_A[i*n+j] = 1;
  }
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++)
      h_B[i*n+j] = 1;
  }
  time = clock();
  matMult(h_A,h_B,h_C,n);
  time = clock() - time;
  cout<<"El Tiempo  es: "<<(((float)time)/CLOCKS_PER_SEC)<<endl;


  /*for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
        cout<<h_C[i*n+j]<<" ; ";
    }
    cout<<endl;
  }
    cout<<endl;*/
  return 0;
}