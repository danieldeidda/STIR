#pragma once

#include <cuda_runtime.h>

__global__
void forwardKernel(
        const float* in_image,
        const float* in_umap,
        float* sino,
        int3 image_dim,
        float3 spacing,
        bool do_atten);

__global__
void backwardKernel(
        const float* in_sino,
        float* image,
        const float* in_umap,
        int3 image_di,
        float3 spacing,
        bool do_atten);
