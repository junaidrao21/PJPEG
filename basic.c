#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <cuda.h>
#include <cuda_runtime.h>
//#include "image.h"

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
		while((r<7) && (c>0))
		{
		    arr[count++] = block[(r++)*8 + (c--)];
		}



		//SECOND HALF OF TRAVERSE
		if(r < 7)
		{
			//MOVE DOWN
			arr[count++] = block[(r++)*8 + c];
		}
		else
		{
			//MOVE RIGHT
			arr[count++] = block[r*8 + (c++)];
		}

		//MOVE UP AND RIGHT	
		while((r>0) && (c<7))
		{
			arr[count++] = block[(r--)*8 + (c++)];
		}
	}
}



int main (int argc, char **argv)
{
//	int imageWidth, imageHeight, matSize;
//	int *hostInputImage, *hostOutputImage;
//	int *deviceInputImage, *deviceOutputImage;

	char test[64] = {  1,  2,  6,  7, 15, 16, 28, 29,
		          3,  5,  8, 14, 17, 27, 30, 43,
	                  4,  9, 13, 18, 26, 31, 42, 44,
			 10, 12, 19, 25, 32, 41, 45, 54,
			 11, 20, 24, 33, 40, 46, 53, 55,
			 21, 23, 34, 39, 47, 52, 56, 61,
	 		 22, 35, 38, 48, 51, 57, 60, 62,
			 36, 37, 49, 50, 58, 59, 63, 64};


	char *A = (char *)calloc(64, sizeof(char));
	char *B = (char *)calloc(64, sizeof(char));

	memcpy(A, test, 64*sizeof(char));

	Traverse(A, B);

	int i;
	for(i=0; i<64; i++)
		printf("B[%d] = (%d).\n", i, (int)B[i]);

    /*
	
	clock_t gpu_start, gpu_end;

	*if (argc != 5)
	{
		printf("Usage: ./lab2 <input-image> <output-image-name> <blur/gaussian/emboss/sharp> <1d/2d kernel>\n");
		exit(1);
	}*/

	// Read in image and convert to readable format
	//read_image_template<int>(argv[1], &hostInputImage, &imageWidth, &imageHeight);

	// Set image size information
	//int img_size = imageWidth * imageHeight * sizeof(int);

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
