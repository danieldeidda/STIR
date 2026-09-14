#include "stir/recon_buildblock/SPECTGPU_projector/SPECTGPUBackwardProjectorCUDA.h"

#include <cuda_runtime.h>

#include "stir/recon_buildblock/SPECTGPU_projector/SPECTGPURotateAndGaussianInterpolate.h"
#include "stir/recon_buildblock/SPECTGPU_projector/SPECTGPUPSFGaussian.h"
#include "stir/recon_buildblock/SPECTGPU_projector/SPECTGPUProjection.h"

START_NAMESPACE_STIR

void run_backward_projection_cuda(
        const RelatedViewgrams<float>& stir_sino,
        DiscretisedDensity<3,float>& stir_image,
        const DiscretisedDensity<3,float>& stir_umap,
        bool do_atten,
        float coll_sigma0_cm,
        float coll_slope,
        int num_views,
        unsigned int block_x,
        unsigned int block_y,
        unsigned int block_z,
        unsigned int grid_x,
        unsigned int grid_y,
        unsigned int grid_z,
        float spacing_x,
        float spacing_y,
        float spacing_z,
        float origin_x,
        float origin_y,
        float origin_z,
        int dim_x,
        int dim_y,
        int dim_z,
        int min_z,
        int min_y,
        int min_x)
{
//    std::cout << "ENTER BP CUDA" << std::endl;
    dim3 cuda_block_dim(
        block_x,
        block_y,
        block_z);

    dim3 cuda_grid_dim(
        grid_x,
        grid_y,
        grid_z);

    float* dev_image;
    cudaMalloc(
        &dev_image,
        stir_image.size_all() * sizeof(float));


//    cudaMemset(
//        dev_image,
//        0,
//        stir_image.size_all() * sizeof(float));
//    this is different than Forward as STIR calls actual_backproject() for every view
    array_to_device(dev_image, stir_image);

    float* rotated_im;
    cudaMalloc(
        &rotated_im,
        stir_image.size_all() * sizeof(float));

    float* rotated_umap;
    float* dev_umap;

    if (do_atten)
    {
        cudaMalloc(
                    &dev_umap,
                    stir_image.size_all() * sizeof(float));
        array_to_device(dev_umap, stir_umap);

        cudaMalloc(
            &rotated_umap,
            stir_image.size_all() * sizeof(float));
//        array_to_device(rotated_umap, stir_umap);
    }




    float3 spacing = make_float3(spacing_x,
                                 spacing_y,
                                 spacing_z);

    float3 origin = make_float3(origin_x,
                                origin_y,
                                origin_z
                                );

    int3 image_dim = make_int3(dim_x, dim_y, dim_z);
    int3 min_indeces = make_int3(min_x, min_y, min_z);

    float* blurred_im;
    if(coll_sigma0_cm>=0 && coll_slope>=0)
        cudaMalloc(&blurred_im, stir_image.size_all() * sizeof(float));


    auto vg_iter = stir_sino.begin();
    //    const Viewgram<float>& vg0 = *vg_iter;
    //    const auto sino_size = vg0.size_all();

//    for (auto vg_iter = stir_sino.begin();
//         vg_iter != stir_sino.end();
//         ++vg_iter)
//    {
        const Viewgram<float>& vg = *vg_iter;
        const auto sino_size = vg.size_all();


        float angle_rad = -vg.get_view_num() * 2.f * M_PI / num_views;

        if (do_atten)
        {
            rotateKernel_pull<<<cuda_grid_dim, cuda_block_dim>>>(
                                                                   dev_umap,
                                                                   rotated_umap,
                                                                   image_dim,
                                                                   spacing,
                                                                   origin,
                                                                   min_indeces,
                                                                   angle_rad);

            cudaDeviceSynchronize();

            auto err0 = cudaGetLastError();
            if (err0 != cudaSuccess)
                error(cudaGetErrorString(err0));
        }

        float* dev_sino;
        cudaMalloc(
            &dev_sino,
            sino_size * sizeof(float));

        array_to_device(dev_sino, vg);

        cudaMemset(rotated_im,
                   0,
                   stir_image.size_all() * sizeof(float));

        //Actual BP
        backwardKernel<<<cuda_grid_dim,cuda_block_dim>>>(
                                                           dev_sino,
                                                           rotated_im,
                                                           rotated_umap,
                                                           image_dim,
                                                           spacing,
                                                           do_atten);

        cudaDeviceSynchronize();

        auto err = cudaGetLastError();
        if (err != cudaSuccess)
            error(cudaGetErrorString(err));

        //        PSF
        if (coll_sigma0_cm>=0 && coll_slope>=0)
        {

            cudaMemset(
                blurred_im,
                0,
                stir_image.size_all() * sizeof(float));

            GaussianConvolutionKernel_push<<<cuda_grid_dim, cuda_block_dim>>>(
                                                                                rotated_im,
                                                                                blurred_im,
                                                                                image_dim,
                                                                                spacing,
                                                                                coll_sigma0_cm,
                                                                                coll_slope);

            cudaDeviceSynchronize();

            auto errpsf_f0 = cudaGetLastError();
            if (errpsf_f0 != cudaSuccess)
                error(cudaGetErrorString(errpsf_f0));


            //Rotation+Interpolation
            rotateKernel_push<<<cuda_grid_dim, cuda_block_dim>>>(
                                                                   blurred_im,
                                                                   dev_image,
                                                                   image_dim,
                                                                   spacing,
                                                                   origin,
                                                                   min_indeces,
                                                                   angle_rad);

            cudaDeviceSynchronize();

            auto errpsf_r = cudaGetLastError();
            if (errpsf_r != cudaSuccess)
                error(cudaGetErrorString(errpsf_r));
        }
        else
        {
            //Rotation+Interpolation
            rotateKernel_push<<<cuda_grid_dim, cuda_block_dim>>>(
                                                                   rotated_im,
                                                                   dev_image,
                                                                   image_dim,
                                                                   spacing,
                                                                   origin,
                                                                   min_indeces,
                                                                   angle_rad);

            cudaDeviceSynchronize();

            auto err1 = cudaGetLastError();
            if (err1 != cudaSuccess)
                error(cudaGetErrorString(err1));
        }

    array_to_host(stir_image, dev_image);
    cudaFree(dev_sino);

    cudaFree(rotated_im);
    if(coll_sigma0_cm>=0 && coll_slope>=0)
        cudaFree(blurred_im);
    cudaFree(dev_image);

    if (do_atten)
    {
        cudaFree(rotated_umap);
        cudaFree(dev_umap);
    }
}

END_NAMESPACE_STIR
