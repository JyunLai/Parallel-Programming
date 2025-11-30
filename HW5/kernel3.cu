#include <cstdio>
#include <cstdlib>
#include <cuda.h>

#define GROUP_SIZE 4

__global__ void mandel_kernel(int *output, size_t pitch, float lower_x, float lower_y, float step_x, float step_y, int res_x, int res_y, int max_iterations)
{
    int startX = (blockIdx.x * blockDim.x + threadIdx.x) * GROUP_SIZE;
    int thisY = blockIdx.y * blockDim.y + threadIdx.y;

    if (thisY < res_y) {
        int *row = (int*)((char*)output + thisY * pitch);
        
        for (int i = 0; i < GROUP_SIZE; i++) {
            int currentX = startX + i;

            if (currentX < res_x) {
                float c_re = lower_x + currentX * step_x;
                float c_im = lower_y + thisY * step_y;

                float z_re = c_re;
                float z_im = c_im;
                int count;
                for (count = 0; count < max_iterations; ++count) {
                    if (z_re * z_re + z_im * z_im > 4.f) {
                        break;
                    }
                    float new_re = z_re * z_re - z_im * z_im;
                    float new_im = 2.f * z_re * z_im;
                    z_re = c_re + new_re;
                    z_im = c_im + new_im;
                }
                row[currentX] = count;
            }
        }
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
    size_t pitch;
    cudaMallocPitch((void **)&d_img, &pitch, res_x * sizeof(int), res_y);

    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks((res_x / GROUP_SIZE + threadsPerBlock.x - 1) / threadsPerBlock.x, (res_y + threadsPerBlock.y - 1) / threadsPerBlock.y);

    mandel_kernel<<<numBlocks, threadsPerBlock>>>(d_img, pitch, lower_x, lower_y, step_x, step_y, res_x, res_y, max_iterations);
    cudaDeviceSynchronize();
    
    cudaMemcpy2D(h_buffer, res_x * sizeof(int), d_img, pitch, res_x * sizeof(int), res_y, cudaMemcpyDeviceToHost);

    memcpy(img, h_buffer, size_bytes);

    cudaFree(d_img);
    cudaFreeHost(h_buffer);
}
