#include <cstdio>
#include <cstdlib>
#include <cuda.h>

__global__ void mandel_kernel(int *output, float lower_x, float lower_y, float step_x, float step_y, int res_x, int res_y, int max_iterations)
{
    int thisX = blockIdx.x * blockDim.x + threadIdx.x;
    int thisY = blockIdx.y * blockDim.y + threadIdx.y;

    if (thisX < res_x && thisY < res_y) {
        float c_re = lower_x + thisX * step_x;
        float c_im = lower_y + thisY * step_y;

        float z_re = c_re;
        float z_im = c_im;

        int i;
        float z_re2 = z_re * z_re;
        float z_im2 = z_im * z_im;

        #pragma unroll 4
        for (i = 0; i < max_iterations; ++i) {
            if (z_re2 + z_im2 > 4.f) {
                break;
            }
            float new_re = z_re2 - z_im2;
            float new_im = 2.f * z_re * z_im;

            z_re = c_re + new_re;
            z_im = c_im + new_im;

            z_re2 = z_re * z_re;
            z_im2 = z_im * z_im;
        }
        int index = thisY * res_x + thisX;
        output[index] = i;
    }
}

// Host front-end function that allocates the memory and launches the GPU kernel
void host_fe(float upper_x,
             float upper_y,
             float lower_x,
             float lower_y,
             int *img,
             int res_x,
             int res_y,
             int max_iterations)
{
    float step_x = (upper_x - lower_x) / (float)res_x;
    float step_y = (upper_y - lower_y) / (float)res_y;

    int size_bytes = res_x * res_y * sizeof(int);
    int *h_buffer;
    cudaHostAlloc((void**)&h_buffer, size_bytes, cudaHostAllocDefault);

    int *d_img;
    cudaMalloc((void **)&d_img, size_bytes);

    dim3 threadsPerBlock(32, 8);
    dim3 numBlocks((res_x + threadsPerBlock.x - 1) / threadsPerBlock.x, (res_y + threadsPerBlock.y - 1) / threadsPerBlock.y);

    mandel_kernel<<<numBlocks, threadsPerBlock>>>(d_img, lower_x, lower_y, step_x, step_y, res_x, res_y, max_iterations);
    cudaDeviceSynchronize();
    
    cudaMemcpy(h_buffer, d_img, size_bytes, cudaMemcpyDeviceToHost);

    memcpy(img, h_buffer, size_bytes);

    cudaFree(d_img);
    cudaFreeHost(h_buffer);
}
