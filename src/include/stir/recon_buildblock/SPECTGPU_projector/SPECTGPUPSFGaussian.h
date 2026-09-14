#pragma once

#include <cuda_runtime.h>

__global__
void GaussianConvolutionKernel_pull(
        const float* in_im,
        float* out_im,
        int3 image_dim,
        float3 spacing,
        float sigma0,
        float slope);



__global__
void GaussianConvolutionKernel_push(
        const float* in_im,
        float* out_im,
        int3 image_dim,
        float3 spacing,
        float sigma0,
        float slope);
