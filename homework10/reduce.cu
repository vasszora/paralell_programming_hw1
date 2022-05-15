#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <chrono>
#include <iostream>

//This is a little wrapper that checks for error codes returned by CUDA API calls
#define cudaCheck(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

__global__ void reduce_atomic(double *c, double *result, int n) {
    int tid = blockIdx.x*blockDim.x+threadIdx.x;
 
    // Make sure we do not go out of bounds
    if (tid < n)
        atomicAdd(result, c[tid]);
}

__global__ void reduce_shared1(double *c, double *result, int n) {
    int tid = blockIdx.x*blockDim.x+threadIdx.x;
    __shared__ double local[1024];
 
    // Make sure we do not go out of bounds
    if (tid < n)
      local[threadIdx.x] = c[tid];
    else
      local[threadIdx.x] = 0.0;

    __syncthreads();
    if (threadIdx.x == 0) {
      double sum = 0.0;
      for (int i = 0; i < 1024; i++) sum += local[i];
      atomicAdd(result, sum);
    }
}

__global__ void reduce_shared2(double *c, double *result, int n) {
    int tid = blockIdx.x*blockDim.x+threadIdx.x;
    __shared__ double local[1024];
 
    // Make sure we do not go out of bounds
    if (tid < n)
      local[threadIdx.x] = c[tid];
    else
      local[threadIdx.x] = 0.0;

    for (int d = blockDim.x >> 1; d >= 1; d >>= 1) {
      __syncthreads();
      if (threadIdx.x < d) local[threadIdx.x] += local[threadIdx.x+d];
    }

    if (threadIdx.x == 0) {
      atomicAdd(result, local[0]);
    }
}


// CUDA kernel. Each thread takes care of one element of c
__global__ void first(double *c, int n)
{
    // Get our global thread ID
    int tid = blockIdx.x*blockDim.x+threadIdx.x;
 
    // Make sure we do not go out of bounds
    if (tid < n)
        c[tid] = tid;
}
 
int main( int argc, char* argv[] )
{
    // Size of vectors
    int n = 1<<30;
    //Host vector
    double *h_c;
    double h_result;
    //Device output vector
    double *d_c;
    double *d_result;
    // Size, in bytes, of each vector
    size_t bytes = n*sizeof(double);
    // Allocate memory on host
    h_c = (double*)malloc(bytes);
    // Allocate memory on GPU
    // Note how we use the cudaCheck wrapper to check for error codes returned
    cudaCheck(cudaMalloc(&d_c, bytes));
    cudaCheck(cudaMalloc(&d_result, sizeof(double)));
    cudaCheck(cudaMemset(d_result, 0, sizeof(double)));
    // Copy host vectors to device
    int blockSize, gridSize;
    // Number of threads in each thread block
    blockSize = 1024;
    // Number of thread blocks in grid
    gridSize = (int)ceil((float)n/blockSize);
    // Execute the kernel
    first<<<gridSize, blockSize>>>(d_c, n);
    // Synchronize
    cudaCheck(cudaDeviceSynchronize());
    auto t1 = std::chrono::high_resolution_clock::now();
    reduce_atomic<<<gridSize,blockSize>>>(d_c, d_result, n);
    cudaCheck(cudaDeviceSynchronize());
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";
    cudaCheck(cudaMemset(d_result, 0, sizeof(double)));
    t1 = std::chrono::high_resolution_clock::now();
    reduce_shared1<<<gridSize,blockSize>>>(d_c, d_result, n);
    cudaCheck(cudaDeviceSynchronize());
    t2 = std::chrono::high_resolution_clock::now();
    std::cout << "took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";

    cudaCheck(cudaMemset(d_result, 0, sizeof(double)));
    t1 = std::chrono::high_resolution_clock::now();
    reduce_shared2<<<gridSize,blockSize>>>(d_c, d_result, n);
    cudaCheck(cudaDeviceSynchronize());
    t2 = std::chrono::high_resolution_clock::now();
    std::cout << "took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
              << " milliseconds\n";




    // Copy array back to host
    cudaCheck(cudaMemcpy( &h_result, d_result, sizeof(double), cudaMemcpyDeviceToHost ));
    // Print resulting array sequentially on the GPU
    printf("%g\n", h_result);
 
    // Release device memory
    cudaFree(d_c);
    cudaFree(d_result);
 
    // Release host memory
    free(h_c);
 
    return 0;
}

