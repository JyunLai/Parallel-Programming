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
        for (i = 0; i < max_iterations; ++i) {
            if (z_re * z_re + z_im * z_im > 4.f) {
                break;
            }
            float new_re = z_re * z_re - z_im * z_im;
            float new_im = 2.f * z_re * z_im;
            z_re = c_re + new_re;
            z_im = c_im + new_im;
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

    int size = res_x * res_y * sizeof(int);
    int *h_buffer = new int[res_x * res_y];

    int *d_img;
    cudaMalloc((void **)&d_img, size);

    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks((res_x + threadsPerBlock.x - 1) / threadsPerBlock.x, (res_y + threadsPerBlock.y - 1) / threadsPerBlock.y);

    mandel_kernel<<<numBlocks, threadsPerBlock>>>(d_img, lower_x, lower_y, step_x, step_y, res_x, res_y, max_iterations);
    cudaDeviceSynchronize();
    
    cudaMemcpy(h_buffer, d_img, size, cudaMemcpyDeviceToHost);

    for (int i = 0; i < res_x * res_y; i++) {
        img[i] = h_buffer[i];
    }
    cudaFree(d_img);
    delete[] h_buffer;
}
