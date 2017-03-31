#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include "image.h"

#define TILE_WIDTH 16

int main (int argc, char **argv)
{
	int imageWidth, imageHeight, matSize;
	int *hostInputImage, *hostOutputImage;
	int *deviceInputImage, *deviceOutputImage;
	
	clock_t gpu_start, gpu_end;

	/*if (argc != 5)
	{
		printf("Usage: ./lab2 <input-image> <output-image-name> <blur/gaussian/emboss/sharp> <1d/2d kernel>\n");
		exit(1);
	}*/

	// Read in image and convert to readable format
	read_image_template<int>(argv[1], &hostInputImage, &imageWidth, &imageHeight);

	// Set image size information
	int img_size = imageWidth * imageHeight * sizeof(int);

	// Allocate memory for image on GPU
	//cudaMalloc((void **)&deviceInputImage, img_size);
	//cudaMalloc((void **)&deviceOutputImage, img_size);
		
	// Copy image to device
	//cudaMemcpy( deviceInputImage, hostInputImage, img_size, cudaMemcpyHostToDevice );
	//cudaMemcpy( deviceMatrix, hostMatrix, sizeof(double)*matSize*matSize, cudaMemcpyHostToDevice );
	
	
	/*****Pre-Processing*****/
	
	
	
	/*****Transforming*******/
	
	
	/*****Quantization*******/
	
	
	/********Encoding********/
	

	return 0;
}