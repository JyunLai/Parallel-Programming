#include "host_fe.h"
#include "helper.h"
#include <stdio.h>
#include <stdlib.h>

void host_fe(int filter_width,
             float *filter,
             int image_height,
             int image_width,
             float *input_image,
             float *output_image,
             cl_device_id *device,
             cl_context *context,
             cl_program *program)
{
    static cl_command_queue queue = NULL;
    static cl_kernel kernel_3x3 = NULL;
    static cl_kernel kernel_5x5 = NULL;
    static cl_kernel kernel_7x7 = NULL;

    static cl_mem d_input = NULL;
    static cl_mem d_output = NULL;
    static cl_mem d_filter = NULL;

    static int current_img_size = 0;
    static int current_filter_size = 0;
    static float *prev_input_ptr = NULL;

    if (queue == NULL) {
        queue = clCreateCommandQueue(*context, *device, 0, NULL);
        kernel_3x3 = clCreateKernel(*program, "convolution_3x3", NULL);
        kernel_5x5 = clCreateKernel(*program, "convolution_5x5", NULL);
        kernel_7x7 = clCreateKernel(*program, "convolution_7x7", NULL);
    }

    cl_kernel kernel;
    if (filter_width == 3) {
        kernel = kernel_3x3;
    }
    else if (filter_width == 5) {
        kernel = kernel_5x5;
    }
    else {
        kernel = kernel_7x7;
    }

    int req_img_size = image_height * image_width * sizeof(float);
    int req_filter_size = filter_width * filter_width *sizeof(float);
    int input_reallocated = 0;
    if (req_img_size > current_img_size) {
        if (d_input) {
            clReleaseMemObject(d_input);
        }
        if (d_output) {
            clReleaseMemObject(d_output);
        }
        d_input = clCreateBuffer(*context, CL_MEM_READ_ONLY, req_img_size, NULL, NULL);
        d_output = clCreateBuffer(*context, CL_MEM_WRITE_ONLY, req_img_size, NULL, NULL);
        current_img_size = req_img_size;
        input_reallocated = 1;
        prev_input_ptr = NULL;
    }
    if (req_filter_size > current_filter_size) {
        if (d_filter) {
            clReleaseMemObject(d_filter);
        }
        d_filter = clCreateBuffer(*context, CL_MEM_READ_ONLY, req_filter_size, NULL, NULL);
        current_filter_size = req_filter_size;
    }

    if (input_reallocated || input_image != prev_input_ptr) {
        clEnqueueWriteBuffer(queue, d_input, CL_FALSE, 0, req_img_size, input_image, 0, NULL, NULL);
        prev_input_ptr = input_image;
    }
    clEnqueueWriteBuffer(queue, d_filter, CL_FALSE, 0, req_filter_size, filter, 0, NULL, NULL);

    clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&d_input);
    clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&d_output);
    clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&d_filter);
    clSetKernelArg(kernel, 3, sizeof(int), (void *)&filter_width);
    clSetKernelArg(kernel, 4, sizeof(int), (void *)&image_height);
    clSetKernelArg(kernel, 5, sizeof(int), (void *)&image_width);

    size_t global_work_size[2];
    global_work_size[0] = (image_width + 31) / 32 * 32;
    global_work_size[1] = (image_height + 15) / 16 * 16;
    size_t local_work_size[2] = {32, 16};

    clEnqueueNDRangeKernel(queue, kernel, 2, NULL, global_work_size, local_work_size, 0, NULL, NULL);
    clEnqueueReadBuffer(queue, d_output, CL_TRUE, 0, req_img_size, output_image, 0, NULL, NULL);
}
