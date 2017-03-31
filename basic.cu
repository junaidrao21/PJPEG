#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include "image.h"

#define TILE_WIDTH 16

//Code written for 1D arrays ahead of time
//for CUDA implementation. 
//
//Assumptions: the 'arr' has been malloc'd
//			   before the function call.
//
//Functionality: linearizes the traversal of
//			     the 8x8 block before the Huffman
//				 encoding scheme is used.
void Traverse(char *block, char *arr)
{
	int count = 0;
	int r = 0;
	int c = 0;

	while(count < 64)
	{
		//FIRST HALF OF TRAVERSE
		if(c < 7)
		{
			//MOVE RIGHT
			arr[count++] = block[r*8 + (c++)];
			
			//ALGORITHM ALWAYS ENDS HERE
			if(count == 64)
				break;
		}	
		else
		{
			//MOVE DOWN
			arr[count++] = block[(r++)*8 + c];
		}

		//MOVE DOWN AND LEFT	
		while((r>0) && (r<7) && (c>0) && (c<7))
		{
			arr[count++] = block[(r++)*8 + (c--)];
		}



		//SECOND HALF OF TRAVERSE
		if(r < 7)
		{
			//MOVE DOWN
			arr[count++] = block[(r--)*8 + c];
		}
		else
		{
			//MOVE RIGHT
			arr[count++] = block[r*8 + (c++)];
		}

		//MOVE UP AND RIGHT	
		while((r>0) && (r<7) && (c>0) && (c<7))
		{
			arr[count++] = block[(r--)*8 + (c++)];
		}
	}
}



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