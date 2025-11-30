__kernel void convolution_3x3(__global float *input,
                          __global float *output,
                          __constant float *filter,
                          int filterWidth,
                          int imageHeight,
                          int imageWidth) 
{
    int gx = get_global_id(0);
    int gy = get_global_id(1);
    int lx = get_local_id(0);
    int ly = get_local_id(1);

    __local float tile[18][34];
    
    if (gy < imageHeight && gx < imageWidth) {
        tile[ly+1][lx+1] = input[gy * imageWidth + gx];
    }
    else {
        tile[ly+1][lx+1] = 0.0f;
    }

    int l_idx = ly * 32 + lx;
    if (l_idx < 100) {
        int t_r, t_c;
        int tmp;

        if (l_idx < 34) {
            t_r = 0;
            t_c = l_idx;
        }
        else if (l_idx < 68) {
            tmp = l_idx - 34;
            t_r = 17;
            t_c = tmp;
        }
        else if (l_idx < 84) {
            tmp = l_idx - 68;
            t_r = tmp + 1;
            t_c = 0;
        }
        else {
            tmp = l_idx - 84;
            t_r = tmp + 1;
            t_c = 33;
        }

        int gh_r = (gy - ly - 1) + t_r;
        int gh_c = (gx - lx - 1) + t_c;
        if (gh_r >= 0 && gh_r < imageHeight && gh_c >= 0 && gh_c < imageWidth) {
            tile[t_r][t_c] = input[gh_r * imageWidth + gh_c];
        }
        else {
            tile[t_r][t_c] = 0.0f;
        }
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    if (gx >= imageWidth || gy >= imageHeight) {
        return;
    }
    float sum = 0.0f;
    int filterIdx = 0;
    #pragma unroll
    for (int k = 0; k < 3; k++) {
        #pragma unroll
        for (int l = 0; l < 3; l++) {
            sum += tile[ly+k][lx+l] * filter[filterIdx++];
        }
    }
    output[gy * imageWidth + gx] = sum;
}

__kernel void convolution_5x5(__global float *input,
                          __global float *output,
                          __constant float *filter,
                          int filterWidth,
                          int imageHeight,
                          int imageWidth) 
{
    int gx = get_global_id(0);
    int gy = get_global_id(1);
    int lx = get_local_id(0);
    int ly = get_local_id(1);

    __local float tile[20][36];
    
    if (gy < imageHeight && gx < imageWidth) {
        tile[ly+2][lx+2] = input[gy * imageWidth + gx];
    }
    else {
        tile[ly+2][lx+2] = 0.0f;
    }
    
    int l_idx = ly * 32 + lx;
    if (l_idx < 208) {
        int t_r, t_c;
        int tmp;

        if (l_idx < 36) {
            t_r = 0;
            t_c = l_idx;
        }
        else if (l_idx < 72) {
            t_r = 1;
            t_c = l_idx - 36;
        }
        else if (l_idx < 108) {
            t_r = 18;
            t_c = l_idx - 72;
        }
        else if (l_idx < 144) {
            t_r = 19;
            t_c = l_idx - 108;
        }
        else {
            tmp = l_idx - 144;
            if (tmp < 32) {
                t_r = 2 + (tmp >> 1);
                t_c = tmp & 1;
            }
            else {
                tmp -= 32;
                t_r = 2 + (tmp >> 1);
                t_c = 34 + (tmp & 1);
            }
        }

        int gh_r = (gy - ly - 2) + t_r;
        int gh_c = (gx - lx - 2) + t_c;
        if (gh_r >= 0 && gh_r < imageHeight && gh_c >= 0 && gh_c < imageWidth) {
            tile[t_r][t_c] = input[gh_r * imageWidth + gh_c];
        }
        else {
            tile[t_r][t_c] = 0.0f;
        }
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    if (gx >= imageWidth || gy >= imageHeight) {
        return;
    }
    float sum = 0.0f;
    int filterIdx = 0;
    #pragma unroll
    for (int k = 0; k < 5; k++) {
        #pragma unroll
        for (int l = 0; l < 5; l++) {
            sum += tile[ly+k][lx+l] * filter[filterIdx++];
        }
    }
    output[gy * imageWidth + gx] = sum;
}

__kernel void convolution_7x7(__global float *input,
                          __global float *output,
                          __constant float *filter,
                          int filterWidth,
                          int imageHeight,
                          int imageWidth) 
{
    int gx = get_global_id(0);
    int gy = get_global_id(1);
    int lx = get_local_id(0);
    int ly = get_local_id(1);

    __local float tile[22][38];
    
    if (gy < imageHeight && gx < imageWidth) {
        tile[ly+3][lx+3] = input[gy * imageWidth + gx];
    }
    else {
        tile[ly+3][lx+3] = 0.0f;
    }
    
    int l_idx = ly * 32 + lx;
    if (l_idx < 324) {
        int t_r, t_c;
        int tmp;

        if (l_idx < 38) {
            t_r = 0;
            t_c = l_idx;
        }
        else if (l_idx < 76) {
            t_r = 1;
            t_c = l_idx - 38;
        }
        else if (l_idx < 114) {
            t_r = 2;
            t_c = l_idx - 76;
        }
        else if (l_idx < 152) {
            t_r = 19;
            t_c = l_idx - 114;
        }
        else if (l_idx < 190) {
            t_r = 20;
            t_c = l_idx - 152;
        }
        else if (l_idx < 228) {
            t_r = 21;
            t_c = l_idx - 190;
        }
        else {
            tmp = l_idx - 228;
            if (tmp < 48) {
                t_r = 3 + tmp / 3;
                t_c = tmp % 3;
            }
            else {
                tmp -= 48;
                t_r = 3 + tmp / 3;
                t_c = 35 + tmp % 3;
            }
        }

        int gh_r = (gy - ly - 3) + t_r;
        int gh_c = (gx - lx - 3) + t_c;
        if (gh_r >= 0 && gh_r < imageHeight && gh_c >= 0 && gh_c < imageWidth) {
            tile[t_r][t_c] = input[gh_r * imageWidth + gh_c];
        }
        else {
            tile[t_r][t_c] = 0.0f;
        }
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    if (gx >= imageWidth || gy >= imageHeight) {
        return;
    }
    float sum = 0.0f;
    int filterIdx = 0;
    #pragma unroll
    for (int k = 0; k < 7; k++) {
        #pragma unroll
        for (int l = 0; l < 7; l++) {
            sum += tile[ly+k][lx+l] * filter[filterIdx++];
        }
    }
    output[gy * imageWidth + gx] = sum;
}
