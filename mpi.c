#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
// #include <cuda.h>
// #include <cuda_runtime.h>
// extern "C" {
#include "lab4.h"
// }
// #define COEFFS(Cu,Cv,u,v) {\
// 			if (u==0) Cu = 1/sqrt(2); else Cu = 1.0;\
// 			if (v==0) Cv = 1/sqrt(2); else Cv = 1.0;}

//Function prototypes
void Traverse(char *block, char *arr, int col);
void Inverse(char *block, char *arr, int col);

//MPI DCT Function
void parallel_dct(int *ky, int *kcb, int *kcr, double *kdy,
			                        double *kdcb, double *kdcr, int row, int col)
{
	// int idx = (threadIdx.x + blockIdx.x * blockDim.x)*8;
	// int idy = (threadIdx.y + blockIdx.y * blockDim.y)*8;
	int rank, numprocs;
  // MPI_Status status;
  // MPI_Comm_rank( MPI_COMM_WORLD, &rank);
  // MPI_Comm_size( MPI_COMM_WORLD, &numprocs);
	int i,j,x,m,n;
	for (m =0; m < row; m+=8)
		for(n=0; n < col; n+=8)
		{
			double temp_dcty[64];
			double temp_dctcb[64];
			double temp_dctcr[64];
			double DCT[8][8] =  {0.3536,  0.3536,  0.3536,  0.3536,  0.3536,  0.3536,  0.3536,  0.3536,
													 0.4904,  0.4157,  0.2778,  0.0975, -0.0975, -0.2778, -0.4157, -0.4904,
													 0.4619,  0.1913, -0.1913, -0.4619, -0.4619, -0.1913,  0.1913,  0.4619,
													 0.4157, -0.0975, -0.4904, -0.2778,  0.2778,  0.4904,  0.0975, -0.4157,
													 0.3536, -0.3536, -0.3536,  0.3536,  0.3536, -0.3536, -0.3536,  0.3536,
													 0.2778, -0.4904,  0.0975,  0.4157, -0.4157, -0.0975,  0.4904, -0.2778,
													 0.1913, -0.4619,  0.4619, -0.1913, -0.1913,  0.4619, -0.4619,  0.1913,
													 0.0975, -0.2778,  0.4157, -0.4904,  0.4904, -0.4157,  0.2788, -0.0975};
				for(i = 0; i < 8; i++)
					for(j = 0; j < 8; j++)
					{
						int index = (m+x) * col + (n+j);
						if(index < row*col){
							temp_dcty[i*8+j] = 0;
							temp_dctcb[i*8+j] = 0;
							temp_dctcr[i*8+j] = 0;
							for(x = 0; x < 8; x++)
							{
								temp_dcty[i*8+j] += DCT[i][x] * (double)ky[(m+x) * col + (n+j)];
								temp_dctcb[i*8+j] += DCT[i][x] * (double)kcb[(m+x) * col + (n+j)];
								temp_dctcr[i*8+j] += DCT[i][x] * (double)kcr[(m+x) * col + (n+j)];
							}
						}
					}
				for(i = 0; i < 8; i++)
					for(j = 0; j < 8; j++)
					{
						int index = (m+x) * col + (n+j);
							if(index < row*col){
							kdy[(m+i) * col + (n+j)] = 0;
							kdcb[(m+i) * col + (n+j)] = 0;
							kdcr[(m+i) * col + (n+j)] = 0;
							for(x = 0; x < 8; x++)
							{
								kdy[(m+i) * col + (n+j)] += temp_dcty[i*8+x] * DCT[j][x];
								kdcb[(m+i) * col + (n+j)] += temp_dctcb[i*8+x] * DCT[j][x];
								kdcr[(m+i) * col + (n+j)] += temp_dctcr[i*8+x] * DCT[j][x];
							}
						}
					}
		}
}

//MPI Inverse DCT Function
void parallel_idct(int *ky, int *kcb, int *kcr, double *kdy,
			                        double *kdcb, double *kdcr, int row, int col)
{
	// int idx = (threadIdx.x + blockIdx.x * blockDim.x)*8;
	// int idy = (threadIdx.y + blockIdx.y * blockDim.y)*8;
	int m,n;
	for (m =0; m < row; m+=8)
		for(n=0; n < col; n+=8)
		{
			double temp_dcty[64];
			double temp_dctcb[64];
			double temp_dctcr[64];
			int i,j,x;
			double IDCT[8][8] = {0.3535,  0.4905,  0.4620,  0.4157,  0.3535,  0.2778,  0.1914,  0.0975,
													 0.3536,  0.4156,  0.1914, -0.0975, -0.3536, -0.4903, -0.4621, -0.2777,
													 0.3534,  0.2780, -0.1914, -0.4905, -0.3534,  0.0973,  0.4622,  0.4156,
													 0.3537,  0.0973, -0.4619, -0.2778,  0.3533,  0.4160, -0.1916, -0.4903,
													 0.3533, -0.0973, -0.4621,  0.2778,  0.3537, -0.4160, -0.1911,  0.4903,
													 0.3537, -0.2780, -0.1913,  0.4905, -0.3537, -0.0973,  0.4618, -0.4156,
													 0.3534, -0.4156,  0.1913,  0.0975, -0.3534,  0.4903, -0.4619,  0.2777,
													 0.3535, -0.4905,  0.4620, -0.4157,  0.3535, -0.2778,  0.1913, -0.0975};
			for(i = 0; i < 8; i++)
				for(j = 0; j < 8; j++)
				{
					int index = (m+x) * col + (n+j);
					if(index < row*col){
						temp_dcty[i*8+j] = 0;
						temp_dctcb[i*8+j] = 0;
						temp_dctcr[i*8+j] = 0;
						for(x = 0; x < 8; x++)
						{
							temp_dcty[i*8+j] += IDCT[i][x] * kdy[(m+x)*col + (n+j)];
							temp_dctcb[i*8+j] += IDCT[i][x] * kdcb[(m+x)*col + (n+j)];
							temp_dctcr[i*8+j] += IDCT[i][x] * kdcr[(m+x)*col + (n+j)];
						}
					}
				}
			for(i = 0; i < 8; i++)
				for(j = 0; j < 8; j++)
				{
					int index = (m+x) * col + (n+j);
					if(index < row*col){
						ky[(m+i) * col + (n+j)] = 0;
						kcb[(m+i) * col + (n+j)] = 0;
						kcr[(m+i) * col + (n+j)] = 0;
						for(x = 0; x < 8; x++)
						{
							ky[(m+i) * col + (n+j)] += (int)(temp_dcty[i*8+x] * IDCT[j][x]);
							kcb[(m+i) * col + (n+j)] += (int)(temp_dctcb[i*8+x] * IDCT[j][x]);
							kcr[(m+i) * col + (n+j)] += (int)(temp_dctcr[i*8+x] * IDCT[j][x]);
						}
					}
				}
		}
}

int main (int argc, char **argv)
{
	//Timing variables
	struct timeval startrgb, endrgb, startdct, enddct, startquant, endquant, starthuff, endhuff, startcmptot, endcmptot;
	struct timeval startirgb, endirgb, startidct, endidct, startiquant, endiquant, startihuff, endihuff, startdectot, enddectot;
	struct timeval starttot, endtot;
	int row, col, char_val, orig_row, orig_col, rank, numprocs;
	int c, i, j, k, m, n;
	char img_type[16];
	int counter = 0;

	unsigned char *img_c;// = (unsigned char *)calloc(row*col*3, sizeof(unsigned char));
	int *img_y;//  = (int *)calloc(img_size*2, sizeof(int));
	int *img_cb; // = (int *)calloc(img_size*2, sizeof(int));
	int *img_cr;
	double *img_dy;//  = (double *)calloc(img_size*2, sizeof(double));
	double *img_dcb;// = (double *)calloc(img_size*2, sizeof(double));
	double *img_dcr;
	char *img_qy;//  = (char *)calloc(img_size*2, sizeof(char));
	char *img_qcb;// = (char *)calloc(img_size*2, sizeof(char));
	char *img_qcr;// = (char *)calloc(img_size*2, sizeof(char));

	//Temp arrays for Traverse()
	char *trav_arr_qy;//  = (char *)calloc(64, sizeof(char));
	char *trav_arr_qcb;// = (char *)calloc(64, sizeof(char));
	char *trav_arr_qcr;// = (char *)calloc(64, sizeof(char));
	char *huff;// = (char *)calloc(row*col*3, sizeof(char));
	FILE *f_in, *f_out, *g_in, *g_out;//  = fopen(argv[1], "r");
	int img_size, per_proc;

	//Quantization Matrix
	unsigned char Q[8][8] = {	16,  11,  10,  16,  24,  40,  51,  61,
														12,  12,  14,  19,  26,  58,  60,  55,
														14,  13,  16,  24,  40,  57,  69,  56,
														14,  17,  22,  29,  51,  87,  80,  62,
														18,  22,  37,  56,  68, 109, 103,  77,
														24,  35,  55,  64,  81, 104, 113,  92,
														49,  64,  78,  87, 103, 121, 120, 101,
														72,  92,  95,  98, 112, 100, 103,  99};

	MPI_Init( &argc,&argv);
	MPI_Comm_rank( MPI_COMM_WORLD, &rank);
	MPI_Comm_size( MPI_COMM_WORLD, &numprocs);
	if(rank == 0){
		gettimeofday(&starttot, NULL);
		gettimeofday(&startcmptot, NULL);

		//Declaration of variables
		f_in  = fopen(argv[1], "r");
		f_out = fopen("output.ppm", "w");
		//double temp_y, temp_cb, temp_cr;
		//double Ci, Cj;

		//Parse header
		fscanf(f_in, "%s\n", img_type);
		fscanf(f_in, "%d %d\n", &orig_col, &orig_row);
		fscanf(f_in, "%d\n", &char_val);

		//Pad row and col if necessary
		col = orig_col;
		row = orig_row;
		if(orig_col % 8 != 0)
		col = orig_col + (8-(orig_col % 8));
		if(orig_row % 8 != 0)
		row = orig_row + (8-(orig_row % 8));

		img_size = col*row;
		//Full RGB matrix
		img_c = (unsigned char *)calloc(row*col*3, sizeof(unsigned char));

		//Separate YCbCr matrices
		img_y  = (int *)calloc(img_size*2, sizeof(int));
		img_cb = (int *)calloc(img_size*2, sizeof(int));
		img_cr = (int *)calloc(img_size*2, sizeof(int));
		//Kernel YCbCr matrices
		// int *ky, *kcb, *kcr;
		// cudaMalloc((void **)&ky, img_size*2*sizeof(int));
		// cudaMalloc((void **)&kcb, img_size*2*sizeof(int));
		// cudaMalloc((void **)&kcr, img_size*2*sizeof(int));

		//Discrete cosine transform matrices
		img_dy  = (double *)calloc(img_size*2, sizeof(double));
		img_dcb = (double *)calloc(img_size*2, sizeof(double));
		img_dcr = (double *)calloc(img_size*2, sizeof(double));
		//Kernel DCT matrices
		// double *kdy, *kdcb, *kdcr;
		// cudaMalloc((void **)&kdy, img_size*2*sizeof(double));
		// cudaMalloc((void **)&kdcb, img_size*2*sizeof(double));
		// cudaMalloc((void **)&kdcr, img_size*2*sizeof(double));

		//Quantization matrices
		img_qy  = (char *)calloc(img_size*2, sizeof(char));
		img_qcb = (char *)calloc(img_size*2, sizeof(char));
		img_qcr = (char *)calloc(img_size*2, sizeof(char));

		//Temp arrays for Traverse()
		trav_arr_qy  = (char *)calloc(64, sizeof(char));
		trav_arr_qcb = (char *)calloc(64, sizeof(char));
		trav_arr_qcr = (char *)calloc(64, sizeof(char));

		//Temp arrays for DCT matrix mult
		// double *temp_dcty  = (double *)calloc(64, sizeof(double));
		// double *temp_dctcb = (double *)calloc(64, sizeof(double));
		// double *temp_dctcr = (double *)calloc(64, sizeof(double));

		//Rearranged matrix for huffman (1D)
		huff = (char *)calloc(row*col*3, sizeof(char));

		//Read in pixel data
		for(i=0; i<orig_row; i++)
			for(j=0; j<orig_col*3; j++)
				fscanf(f_in, "%c", &img_c[i*col*3 + j]);

		// for(m=0;m<64;m++){
		// 	if(m%8==0)
		// 		printf("\n");
		// 	printf("%6d",img_c[m]);
		// }

			//RBG -> YCbCr
		gettimeofday(&startrgb, NULL);
		for(i=0; i<row; i++)
			for(j=0,k=0; j<col*3; j+=3,k++)
			{
				img_y[i * col + k]  = (0.299)*img_c[i*col*3 + j] + (0.587)*img_c[i*col*3 + (j+1)] + (0.114)*img_c[i*col*3 + (j+2)];
			  img_cb[i * col + k] = 128 - (0.168736)*img_c[i*col*3 + j] - (0.331264)*img_c[i*col*3 + (j+1)] + (0.5)*img_c[i*col*3 + (j+2)];
			  img_cr[i * col + k] = 128 + (0.5)*img_c[i*col*3 + j] - (0.418688)*img_c[i*col*3 + (j+1)] - (0.081312)*img_c[i*col*3 + (j+2)];
			}

			//Center
		for(i=0; i<row; i++)
			for(j=0; j<col; j++)
			{
				img_y[i * col + j]  -= 128;
				img_cb[i * col + j] -= 128;
				img_cr[i * col + j] -= 128;
			}
		gettimeofday(&endrgb, NULL);

		per_proc = row/numprocs;
		// if(per_proc%8){
		// 	per_proc += 8-(per_proc%8);
		// }
	}

		// printf("Before DCT\n");
		// for(m=0;m<8;m++){
		// 	for(n=0;n<8;n++)
		// 		printf("%4d", img_cr[m*col+n]);
		// 	printf("\n");
		// }

		//Copying YCbCr to Device
		// cudaMemcpy( ky, img_y, img_size*2*sizeof(int), cudaMemcpyHostToDevice );
		// cudaMemcpy( kcb, img_cb, img_size*2*sizeof(int), cudaMemcpyHostToDevice );
		// cudaMemcpy( kcr, img_cr, img_size*2*sizeof(int), cudaMemcpyHostToDevice );
		// dim3 dimGrid(ceil((double)col / 8),
		// 	     ceil((double)row / 8),1);
		// dim3 dimBlock(8, 8, 1);
		//Discrete Cosine Transform
	MPI_Bcast(&per_proc, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&col, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int *sub_img_y  = (int *)calloc(per_proc*col, sizeof(int));
	int *sub_img_cb = (int *)calloc(per_proc*col, sizeof(int));
	int *sub_img_cr = (int *)calloc(per_proc*col, sizeof(int));
	double *sub_img_dy  = (double *)calloc(per_proc*col, sizeof(double));
	double *sub_img_dcb = (double *)calloc(per_proc*col, sizeof(double));
	double *sub_img_dcr = (double *)calloc(per_proc*col, sizeof(double));

	gettimeofday(&startdct, NULL);
	//Scatter the matrices all process
	int snd_count = per_proc*col;
	MPI_Scatter((void *)img_y, snd_count, MPI_INT, (void *)sub_img_y,
              snd_count, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter((void *)img_cb, snd_count, MPI_INT, (void *)sub_img_cb,
              snd_count, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter((void *)img_cr, snd_count, MPI_INT, (void *)sub_img_cr,
	            snd_count, MPI_INT, 0, MPI_COMM_WORLD);
	parallel_dct(sub_img_y, sub_img_cb, sub_img_cr, sub_img_dy, sub_img_dcb, sub_img_dcr, per_proc, col);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gather(sub_img_dy, per_proc*col, MPI_DOUBLE, img_dy, per_proc*col, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(sub_img_dcb, per_proc*col, MPI_DOUBLE, img_dcb, per_proc*col, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(sub_img_dcr, per_proc*col, MPI_DOUBLE, img_dcr, per_proc*col, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		//MPI_Barrier(MPI_COMM_WORLD);
		//MPI_Finalize();
	gettimeofday(&enddct, NULL);
		//Copying results back from device
		// cudaMemcpy( img_dy, kdy, img_size*2*sizeof(double), cudaMemcpyDeviceToHost );
		// cudaMemcpy( img_dcb, kdcb, img_size*2*sizeof(double), cudaMemcpyDeviceToHost );
		// cudaMemcpy( img_dcr, kdcr, img_size*2*sizeof(double), cudaMemcpyDeviceToHost );

		// cudaFree(ky);
		// cudaFree(kcb);
		// cudaFree(kcr);
		// cudaFree(kdy);
		// cudaFree(kdcb);
		// cudaFree(kdcr);

		// printf("After DCT\n");
		// for(m=0;m<8;m++){
		// 	for(n=0;n<8;n++)
		// 		printf("%6.2f", img_dcr[m*col+n]);
		// 	printf("\n");
		// }

	if(rank == 0){
		//Quantization
		gettimeofday(&startquant, NULL);
		for(m=0; m<row; m+=8)
			for(n=0; n<col; n+=8)
				for(i=0; i<8; i++)
					for(j=0; j<8; j++)
					{
						img_qy[(m+i) * col + (n+j)]  = (char)rint((img_dy[(m+i)*col + (n+j)]/Q[i][j]));
						img_qcb[(m+i) * col + (n+j)] = (char)rint((img_dcb[(m+i)*col + (n+j)]/Q[i][j]));
						img_qcr[(m+i) * col + (n+j)] = (char)rint((img_dcr[(m+i)*col + (n+j)]/Q[i][j]));
					}
		gettimeofday(&endquant, NULL);

		// printf("After Quant\n");
		// for(m=0;m<8;m++){
		// 	for(n=0;n<8;n++)
		// 		printf("%6d", img_qcr[m*row+n]);
		// 	printf("\n");
		// }

		//Linearization of each 8x8 block before compression
		for(m=0; m<row; m+=8)
			for(n=0; n<col; n+=8)
			{
				//Linearization
				Traverse(img_qy+(m*col+n), trav_arr_qy, col);
				Traverse(img_qcb+(m*col+n), trav_arr_qcb, col);
				Traverse(img_qcr+(m*col+n), trav_arr_qcr, col);

				//Combination into single Huffman array
				for(c=0; c<64; c++, counter++)
				{
					huff[counter] = trav_arr_qy[c];
					huff[col*row+counter] = trav_arr_qcb[c];
					huff[col*row*2+counter] = trav_arr_qcr[c];
				}
			}

		// printf("After Traverse");
		// for(m=row*col*2;m<row*col*2+64;m++){
		// 	if(m%8 == 0)
		// 		printf("\n");
		// 	printf("%6d", huff[m]);
		// }
		// printf("\n");

		//Write out combined matrix to "output.ppm"
		fprintf(f_out, "%s\n", img_type);
		fprintf(f_out, "%d %d\n", orig_col, orig_row);
		fprintf(f_out, "%d\n", char_val);
		for(m=0; m<row*col*3; m++)
			fprintf(f_out, "%c", huff[m]);

		fclose(f_in);
		fclose(f_out);

		//Huffman Compression
		gettimeofday(&starthuff, NULL);
		//system("./huff -c output.ppm output.ppm.huf");
		char arg1[] = "c";
		char arg2[] = "output.ppm";
		char arg3[] = "output.ppm.huf";
		lab4(arg1, arg2, arg3);
		gettimeofday(&endhuff, NULL);

		gettimeofday(&endcmptot, NULL);

		//Timing outputs in milliseconds
		// printf("file: %s\n",argv[1]);
		// printf("\n");
		// printf("%6.3f,", (double)(endrgb.tv_usec - startrgb.tv_usec) / 1000000 + (endrgb.tv_sec - startrgb.tv_sec));
		// printf("%6.3f,", (double)(enddct.tv_usec - startdct.tv_usec) / 1000000 + (enddct.tv_sec - startdct.tv_sec));
		// printf("%6.3f,", (double)(endquant.tv_usec - startquant.tv_usec) / 1000000 + (endquant.tv_sec - startquant.tv_sec));
		// printf("%6.3f,", (double)(endhuff.tv_usec - starthuff.tv_usec) / 1000000 + (endhuff.tv_sec - starthuff.tv_sec));
		// printf("%6.3f\n", (double)(endcmptot.tv_usec - startcmptot.tv_usec) / 1000000 + (endcmptot.tv_sec - startcmptot.tv_sec));



		gettimeofday(&startdectot, NULL);
		//printf("Decompressing...\n");

		//Huffman Decompression
		gettimeofday(&startihuff, NULL);
		char arg4[] = "d";
		char arg5[] = "output.ppm.uhuf";
		lab4(arg4, arg3, arg5);
		gettimeofday(&endihuff, NULL);


		//End of JPEG encoding



		//////////////////////////////////////////////////////////////////////////////
		//  BEGIN DECODING OF UNCOMPRESSED JPEG IMAGE
		//////////////////////////////////////////////////////////////////////////////

		//Open uncompressed huffman file
		g_in  = fopen("output.ppm.uhuf", "r");
		g_out = fopen("result.ppm", "w");

		//Parse header
		fscanf(g_in, "%s\n", img_type);
		fscanf(g_in, "%d %d\n", &orig_col, &orig_row);
		fscanf(g_in, "%d\n", &char_val);

		col = orig_col;
		row = orig_row;
		if(orig_col % 8 != 0)
		 	col = orig_col + (8-(orig_col % 8));
		if(orig_row % 8 != 0)
		 	row = orig_row + (8-(orig_row % 8));


		//Prepare huffman array
		memset(huff, 0, row*col*3*sizeof(char));

		//printf("Reading in file\n");
		//Read in file contents
		for(m=0; m<row*col*3; m++)
				fscanf(g_in, "%c", &huff[m]);

		// printf("Before ITraverse");
		// for(m=row*col*2;m<row*col*2+64;m++){
		// 	if(m%8 == 0)
		// 		printf("\n");
		// 	printf("%6d", huff[m]);
		// }
		// printf("\n");

		//printf("Reverse Linear\n");
		//Reverse the linearization
		counter = 0;
		for(i=0; i<3; i++)
			for(m=0; m<row; m+=8)
				for(n=0; n<col; n+=8, counter++)
					if(i==0)
						Inverse(&img_qy[m * col + n], (huff + counter*64), col);
					else if(i==1)
						Inverse(&img_qcb[m * col + n], (huff + counter*64), col);
					else
						Inverse(&img_qcr[m * col + n], (huff + counter*64), col);

		// printf("Before IQuant\n");
		// for(m=0;m<8;m++){
		// 	for(n=0;n<8;n++)
		// 		printf("%6d", img_qcr[m*row+n]);
		// 	printf("\n");
		// }

		//printf("IQuant\n");
		gettimeofday(&startiquant, NULL);
		//Inverse Quantization
		for(m=0; m<row; m+=8)
			for(n=0; n<col; n+=8)
				for (i=0; i<8; i++)
					for (j=0; j<8; j++)
					{
								img_dy[(m+i)*col + (n+j)] = (img_qy[(m+i)*col + (n+j)]*Q[i][j]);
								img_dcb[(m+i)*col + (n+j)] = (img_qcb[(m+i)*col + (n+j)]*Q[i][j]);
								img_dcr[(m+i)*col + (n+j)] = (img_qcr[(m+i)*col + (n+j)]*Q[i][j]);
					}
		gettimeofday(&endiquant, NULL);
	}

		// printf("Before IDCT\n");
		// for(m=0;m<8;m++){
		// 	for(n=0;n<8;n++)
		// 		printf("%6.2f", img_dcr[m*col+n]);
		// 	printf("\n");
		// }

		//Iniatializing Device matrices
		// cudaMalloc((void **)&kdy, img_size*2*sizeof(double));
		// cudaMalloc((void **)&kdcb, img_size*2*sizeof(double));
		// cudaMalloc((void **)&kdcr, img_size*2*sizeof(double));
		// cudaMalloc((void **)&ky, img_size*2*sizeof(int));
		// cudaMalloc((void **)&kcb, img_size*2*sizeof(int));
		// cudaMalloc((void **)&kcr, img_size*2*sizeof(int));

		//Copying DCT matrices to device memory
		// cudaMemcpy( kdy, img_dy, img_size*2*sizeof(double), cudaMemcpyHostToDevice );
		// cudaMemcpy( kdcb, img_dcb, img_size*2*sizeof(double), cudaMemcpyHostToDevice );
		// cudaMemcpy( kdcr, img_dcr, img_size*2*sizeof(double), cudaMemcpyHostToDevice );
		MPI_Barrier(MPI_COMM_WORLD);
		gettimeofday(&startidct, NULL);
		// IDCT
		MPI_Scatter((void *)img_dy, snd_count, MPI_DOUBLE, (void *)sub_img_dy,
	              snd_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatter((void *)img_dcb, snd_count, MPI_DOUBLE, (void *)sub_img_dcb,
	              snd_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatter((void *)img_dcr, snd_count, MPI_DOUBLE, (void *)sub_img_dcr,
		            snd_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		parallel_idct(sub_img_y, sub_img_cb, sub_img_cr, sub_img_dy, sub_img_dcb, sub_img_dcr, per_proc, col);

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Gather(sub_img_y, per_proc*col, MPI_INT, img_y, per_proc*col, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Gather(sub_img_cb, per_proc*col, MPI_INT, img_cb, per_proc*col, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Gather(sub_img_cr, per_proc*col, MPI_INT, img_cr, per_proc*col, MPI_INT, 0, MPI_COMM_WORLD);


		// cudaDeviceSynchronize();

		gettimeofday(&endidct, NULL);

		// cudaMemcpy( img_y, ky, img_size*2*sizeof(int), cudaMemcpyDeviceToHost );
		// cudaMemcpy( img_cb, kcb, img_size*2*sizeof(int), cudaMemcpyDeviceToHost );
		// cudaMemcpy( img_cr, kcr, img_size*2*sizeof(int), cudaMemcpyDeviceToHost );

		// printf("After IDCT\n");
		// for(m=0;m<8;m++){
		// 	for(n=0;n<8;n++)
		// 		printf("%6d", img_cr[m*col+n]);
		// 	printf("\n");
		// }

	if(rank == 0){
		gettimeofday(&startirgb, NULL);
		//Un-Center
		for(i=0; i<row; i++)
			for(j=0; j<col; j++)
			{
				img_y[i * col + j]  += 128;
				img_cb[i * col + j] += 128;
				img_cr[i * col + j] += 128;
			}
			// printf("YCbCr[0][0]\n");
			// printf("%8d %8d %8d\n",img_y[0][0], img_cb[0][0], img_cr[0][0]);


		int tempr, tempg, tempb;

		// printf("YCbCr->RGB\n");
		//YCbCr back to RGB
		for(m=0; m<row; m++)
			for(n=0, j=0; n<col; n++, j+=3)
			{
				tempr = (img_y[m * col + n] + 1.40200 * (img_cr[m * col + n] - 128));
				tempg = (img_y[m * col + n] - 0.34414 * (img_cb[m * col + n] - 128) - 0.71414 * (img_cr[m * col + n] - 128));
				tempb = (img_y[m * col + n] + 1.77200 * (img_cb[m * col + n] - 128));

				if(tempr > 255)
					tempr = 255;
				if(tempr < 0)
					tempr = 0;
				if(tempg > 255)
					tempg = 255;
				if(tempg < 0)
					tempg = 0;
				if(tempb > 255)
					tempb = 255;
				if(tempb < 0)
					tempb = 0;

				img_c[m*col*3 + j] = (unsigned char)tempr;
				img_c[m*col*3 + (j+1)] = (unsigned char)tempg;
				img_c[m*col*3 + (j+2)] = (unsigned char)tempb;
			}
			gettimeofday(&endirgb, NULL);

			// printf("RGB[0][0]\n");
			// printf("%8d %8d %8d\n",img_c[0][0], img_c[0][1], img_c[0][2]);

		// printf("Write out\n");
		//Write out the final, reconstructed RGB image
		fprintf(g_out, "%s\n", img_type);
		fprintf(g_out, "%d %d\n", orig_col, orig_row);
		fprintf(g_out, "%d\n", char_val);
		for(i=0; i<orig_row; i++)
			for(j=0; j<orig_col*3; j++)
				fprintf(g_out, "%c", img_c[i*col*3 + j]);

				// for(m=0;m<64;m++){
				// 	if(m%8==0)
				// 		printf("\n");
				// 	printf("%6d",img_c[m]);
				// }

		//Clean up
		fclose(g_in);
		fclose(g_out);

		gettimeofday(&enddectot, NULL);
		// printf("\n");
		// printf("%6.3f,", (double)(endirgb.tv_usec - startirgb.tv_usec) / 1000000 + (endirgb.tv_sec - startirgb.tv_sec));
		// printf("%6.3f,", (double)(endidct.tv_usec - startidct.tv_usec) / 1000000 + (endidct.tv_sec - startidct.tv_sec));
		// printf("%6.3f,", (double)(endiquant.tv_usec - startiquant.tv_usec) / 1000000 + (endiquant.tv_sec - startiquant.tv_sec));
		// printf("%6.3f,", (double)(endihuff.tv_usec - startihuff.tv_usec) / 1000000 + (endihuff.tv_sec - startihuff.tv_sec));
		// printf("%6.3f\n", (double)(enddectot.tv_usec - startdectot.tv_usec) / 1000000 + (enddectot.tv_sec - startdectot.tv_sec));

		// printf("freeing memory\n");
		//Free allocated memory
		// cudaFree(ky);
		// cudaFree(kcb);
		// cudaFree(kcr);
		// cudaFree(kdy);
		// cudaFree(kdcb);
		// cudaFree(kdcr);
		free(img_c);
		free(img_y);
		free(img_cb);
		free(img_cr);
		free(img_dy);
		free(img_dcb);
		free(img_dcr);
		// free(img_qy);
		// free(img_qcb);
		// free(img_qcr);
		free(trav_arr_qy);
		free(trav_arr_qcb);
		free(trav_arr_qcr);
		free(huff);

		// system("rm output.ppm.huf output.ppm.uhuf");

		gettimeofday(&endtot, NULL);
			printf("%s,%d,%d,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f\n",
							argv[1],orig_col,orig_row,
							(double)(endrgb.tv_usec - startrgb.tv_usec) / 1000000 + (endrgb.tv_sec - startrgb.tv_sec),
							(double)(enddct.tv_usec - startdct.tv_usec) / 1000000 + (enddct.tv_sec - startdct.tv_sec),
							(double)(endquant.tv_usec - startquant.tv_usec) / 1000000 + (endquant.tv_sec - startquant.tv_sec),
							(double)(endhuff.tv_usec - starthuff.tv_usec) / 1000000 + (endhuff.tv_sec - starthuff.tv_sec),
							(double)(endcmptot.tv_usec - startcmptot.tv_usec) / 1000000 + (endcmptot.tv_sec - startcmptot.tv_sec),
							(double)(endirgb.tv_usec - startirgb.tv_usec) / 1000000 + (endirgb.tv_sec - startirgb.tv_sec),
							(double)(endidct.tv_usec - startidct.tv_usec) / 1000000 + (endidct.tv_sec - startidct.tv_sec),
							(double)(endiquant.tv_usec - startiquant.tv_usec) / 1000000 + (endiquant.tv_sec - startiquant.tv_sec),
							(double)(endihuff.tv_usec - startihuff.tv_usec) / 1000000 + (endihuff.tv_sec - startihuff.tv_sec),
							(double)(enddectot.tv_usec - startdectot.tv_usec) / 1000000 + (enddectot.tv_sec - startdectot.tv_sec),
							(double)(endtot.tv_usec - starttot.tv_usec) / 1000000 + (endtot.tv_sec - starttot.tv_sec));
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	// printf("%s,%d,%d,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f\n",
	// 				argv[1],orig_col,orig_row,
	// 				(double)(endrgb.tv_usec - startrgb.tv_usec) / 1000000 + (endrgb.tv_sec - startrgb.tv_sec),
	// 				(double)(enddct.tv_usec - startdct.tv_usec) / 1000000 + (enddct.tv_sec - startdct.tv_sec),
	// 				(double)(endquant.tv_usec - startquant.tv_usec) / 1000000 + (endquant.tv_sec - startquant.tv_sec),
	// 				(double)(endhuff.tv_usec - starthuff.tv_usec) / 1000000 + (endhuff.tv_sec - starthuff.tv_sec),
	// 				(double)(endcmptot.tv_usec - startcmptot.tv_usec) / 1000000 + (endcmptot.tv_sec - startcmptot.tv_sec),
	// 				(double)(endirgb.tv_usec - startirgb.tv_usec) / 1000000 + (endirgb.tv_sec - startirgb.tv_sec),
	// 				(double)(endidct.tv_usec - startidct.tv_usec) / 1000000 + (endidct.tv_sec - startidct.tv_sec),
	// 				(double)(endiquant.tv_usec - startiquant.tv_usec) / 1000000 + (endiquant.tv_sec - startiquant.tv_sec),
	// 				(double)(endihuff.tv_usec - startihuff.tv_usec) / 1000000 + (endihuff.tv_sec - startihuff.tv_sec),
	// 				(double)(enddectot.tv_usec - startdectot.tv_usec) / 1000000 + (enddectot.tv_sec - startdectot.tv_sec),
	// 				(double)(endtot.tv_usec - starttot.tv_usec) / 1000000 + (endtot.tv_sec - starttot.tv_sec));


	//Scott is n00b
	return 0;
}

////////////////////////////////////////////////////////////////////////////////
//  FUNCTION DEFINITIONS of Traverse() and Inverse()
////////////////////////////////////////////////////////////////////////////////
void Traverse(char *block, char *arr, int col)
{
	int count = 0;
	int r = 0;
	int c = 0;

	while(count < 64)
	{
		if(c < 7)
		{
			arr[count++] = block[(r*col) + (c++)];
			if(count == 64)
				break;
		}
		else
			arr[count++] = block[(r++)*col + c];

		while((r<7) && (c>0))
		    arr[count++] = block[(r++)*col + (c--)];

		if(r < 7)
			arr[count++] = block[(r++)*col + c];
		else
			arr[count++] = block[(r*col) + (c++)];

		while((r>0) && (c<7))
			arr[count++] = block[(r--)*col + (c++)];
	}
}

void Inverse(char *block, char *arr, int col)
{
	int count = 0;
	int r = 0;
	int c = 0;

	while(count < 64)
	{
		if(c < 7)
		{
			block[(r*col) + (c++)] = arr[count++];

			if(count == 64)
				break;
		}
		else
			block[(r++)*col + c] = arr[count++];

		while((r<7) && (c>0))
		    block[(r++)*col + (c--)] = arr[count++];

		if(r < 7)
			block[(r++)*col + c] = arr[count++];
		else
			block[(r*col) + (c++)] = arr[count++];

		while((r>0) && (c<7))
			block[(r--)*col + (c++)] = arr[count++];
	}
}
