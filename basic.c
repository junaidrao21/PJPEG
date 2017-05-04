#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "lab4.h"

#define COEFFS(Cu,Cv,u,v) {\
			if (u==0) Cu = 1/sqrt(2); else Cu = 1.0;\
			if (v==0) Cv = 1/sqrt(2); else Cv = 1.0;}

//Function prototypes
void Traverse(char *block, char *arr, int col);
void Inverse(char *block, char *arr, int col);


int main (int argc, char **argv)
{
	//Timing variables
	struct timeval startrgb, endrgb, startdct, enddct, startquant, endquant, starthuff, endhuff, startcmptot, endcmptot;
	struct timeval startirgb, endirgb, startidct, endidct, startiquant, endiquant, startihuff, endihuff, startdectot, enddectot;
	struct timeval starttot, endtot;
	gettimeofday(&starttot, NULL);
	gettimeofday(&startcmptot, NULL);

	//Declaration of variables
	FILE *f_in  = fopen(argv[1], "r");
	FILE *f_out = fopen("output.ppm", "w");
	char img_type[16];
	int row, col, char_val, orig_row, orig_col;
	int c, i, j, k, m, n, x, y;
	int counter = 0;
	double temp_y, temp_cb, temp_cr;
	double Ci, Cj;

	//Quantization Matrix
	unsigned char Q[8][8] = {	16,  11,  10,  16,  24,  40,  51,  61,
	 													12,  12,  14,  19,  26,  58,  60,  55,
														14,  13,  16,  24,  40,  57,  69,  56,
														14,  17,  22,  29,  51,  87,  80,  62,
														18,  22,  37,  56,  68, 109, 103,  77,
														24,  35,  55,  64,  81, 104, 113,  92,
														49,  64,  78,  87, 103, 121, 120, 101,
														72,  92,  95,  98, 112, 100, 103,  99};

	//DCT cosine matrix
	double DCT[8][8] =  {0.3536,  0.3536,  0.3536,  0.3536,  0.3536,  0.3536,  0.3536,  0.3536,
											 0.4904,  0.4157,  0.2778,  0.0975, -0.0975, -0.2778, -0.4157, -0.4904,
											 0.4619,  0.1913, -0.1913, -0.4619, -0.4619, -0.1913,  0.1913,  0.4619,
											 0.4157, -0.0975, -0.4904, -0.2778,  0.2778,  0.4904,  0.0975, -0.4157,
											 0.3536, -0.3536, -0.3536,  0.3536,  0.3536, -0.3536, -0.3536,  0.3536,
											 0.2778, -0.4904,  0.0975,  0.4157, -0.4157, -0.0975,  0.4904, -0.2778,
											 0.1913, -0.4619,  0.4619, -0.1913, -0.1913,  0.4619, -0.4619,  0.1913,
											 0.0975, -0.2778,  0.4157, -0.4904,  0.4904, -0.4157,  0.2788, -0.0975};

	//IDCT cosine matrix
	double IDCT[8][8] = {0.3535,  0.4905,  0.4620,  0.4157,  0.3535,  0.2778,  0.1914,  0.0975,
											 0.3536,  0.4156,  0.1914, -0.0975, -0.3536, -0.4903, -0.4621, -0.2777,
											 0.3534,  0.2780, -0.1914, -0.4905, -0.3534,  0.0973,  0.4622,  0.4156,
											 0.3537,  0.0973, -0.4619, -0.2778,  0.3533,  0.4160, -0.1916, -0.4903,
											 0.3533, -0.0973, -0.4621,  0.2778,  0.3537, -0.4160, -0.1911,  0.4903,
											 0.3537, -0.2780, -0.1913,  0.4905, -0.3537, -0.0973,  0.4618, -0.4156,
											 0.3534, -0.4156,  0.1913,  0.0975, -0.3534,  0.4903, -0.4619,  0.2777,
											 0.3535, -0.4905,  0.4620, -0.4157,  0.3535, -0.2778,  0.1913, -0.0975};

	//Parse header
	fscanf(f_in, "%s\n", img_type);
	fscanf(f_in, "%d %d\n", &orig_col, &orig_row);
	fscanf(f_in, "%d\n", &char_val);


	// If greyscale, say it doesn't work and exit
	if(strcmp(img_type,"P5") == 0){
		printf("This program does not currently work on greyscale (P5) images.\nChoose a color ppm image (P6).\n");
		exit(1);
	}
	else if(strcmp(img_type,"P6") != 0){
		printf("This image is not a color ppm (P6) image. Try again with a P6 image.\n");
		exit(1);
	}

	//Pad row and col if necessary
	col = orig_col;
	row = orig_row;
	if(orig_col % 8 != 0)
	col = orig_col + (8-(orig_col % 8));
	if(orig_row % 8 != 0)
	row = orig_row + (8-(orig_row % 8));
	// printf("%d %d\n",col, row);

	//Full RGB matrix
	unsigned char *img_c = (unsigned char *)calloc(row*col*3, sizeof(unsigned char));

	//Separate YCbCr matrices
	int *img_y  = (int *)calloc(row*col*2, sizeof(int));
	int *img_cb = (int *)calloc(row*col*2, sizeof(int));
	int *img_cr = (int *)calloc(row*col*2, sizeof(int));

	//Discrete cosine transform matrices
	double *img_dy  = (double *)calloc(row*col*2, sizeof(double));
	double *img_dcb = (double *)calloc(row*col*2, sizeof(double));
	double *img_dcr = (double *)calloc(row*col*2, sizeof(double));

	//Quantization matrices (1D)
	char *img_qy  = (char *)calloc(row*col*2, sizeof(char));
	char *img_qcb = (char *)calloc(row*col*2, sizeof(char));
	char *img_qcr = (char *)calloc(row*col*2, sizeof(char));

	//Temp arrays for Traverse()
	char *trav_arr_qy  = (char *)calloc(64, sizeof(char));
	char *trav_arr_qcb = (char *)calloc(64, sizeof(char));
	char *trav_arr_qcr = (char *)calloc(64, sizeof(char));

	//Temp arrays for DCT matrix mult
	double *temp_dcty  = (double *)calloc(64, sizeof(double));
	double *temp_dctcb = (double *)calloc(64, sizeof(double));
	double *temp_dctcr = (double *)calloc(64, sizeof(double));

	//Rearranged matrix for huffman (1D)
	char *huff = (char *)calloc(row*col*3, sizeof(char));

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

	// printf("Before DCT\n");
	// for(m=0;m<8;m++){
	// 	for(n=0;n<8;n++)
	// 		printf("%4d", img_cr[m*col+n]);
	// 	printf("\n");
	// }

	//Discrete Cosine Transform
	gettimeofday(&startdct, NULL);

	for(m=0; m<row; m+=8)
		for(n=0; n<col; n+=8)
		{
			for(i = 0; i < 8; i++)
				for(j = 0; j < 8; j++)
				{
					temp_dcty[i*8+j] = 0;
					temp_dctcb[i*8+j] = 0;
					temp_dctcr[i*8+j] = 0;
					for(x = 0; x < 8; x++)
					{
						temp_dcty[i*8+j] += DCT[i][x] * (double)img_y[(m+x) * col + (n+j)];
						temp_dctcb[i*8+j] += DCT[i][x] * (double)img_cb[(m+x) * col + (n+j)];
						temp_dctcr[i*8+j] += DCT[i][x] * (double)img_cr[(m+x) * col + (n+j)];
					}
				}
			for(i = 0; i < 8; i++)
				for(j = 0; j < 8; j++)
				{
					img_dy[(m+i) * col + (n+j)] = 0;
					img_dcb[(m+i) * col + (n+j)] = 0;
					img_dcr[(m+i) * col + (n+j)] = 0;
					for(x = 0; x < 8; x++)
					{
						img_dy[(m+i) * col + (n+j)] += temp_dcty[i*8+x] * DCT[j][x];
						img_dcb[(m+i) * col + (n+j)] += temp_dctcb[i*8+x] * DCT[j][x];
						img_dcr[(m+i) * col + (n+j)] += temp_dctcr[i*8+x] * DCT[j][x];
					}
				}
		}

	gettimeofday(&enddct, NULL);

	// printf("After DCT\n");
	// for(m=0;m<8;m++){
	// 	for(n=0;n<8;n++)
	// 		printf("%6.2f", img_dcr[m*col+n]);
	// 	printf("\n");
	// }

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
	int count = 1;
	// printf("img_qcr: start = (%p)\n", &img_qcr[row*col-1]);
	for(m=0; m<row; m+=8)
		for(n=0; n<col; n+=8)
		{
			//Linearization
			Traverse(img_qy+(m*col+n), trav_arr_qy, col);
			Traverse(img_qcb+(m*col+n), trav_arr_qcb, col);
			Traverse(img_qcr+(m*col+n), trav_arr_qcr, col);
			// printf("%d   ",count);
			// printf("%d\n", img_qcr[m*col+n]);
			count++;
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
	lab4("c", "output.ppm", "output.ppm.huf");
	gettimeofday(&endhuff, NULL);

	gettimeofday(&endcmptot, NULL);

	//Timing outputs in milliseconds
	// printf("\n");
	// printf("RGB->YCbCr = %6.3f\n", (double)(endrgb.tv_usec - startrgb.tv_usec) / 1000000 + (endrgb.tv_sec - startrgb.tv_sec));
	// printf("DCT =        %6.3f\n", (double)(enddct.tv_usec - startdct.tv_usec) / 1000000 + (enddct.tv_sec - startdct.tv_sec));
	// printf("Quant =      %6.3f\n", (double)(endquant.tv_usec - startquant.tv_usec) / 1000000 + (endquant.tv_sec - startquant.tv_sec));
	// printf("Huffman =    %6.3f\n", (double)(endhuff.tv_usec - starthuff.tv_usec) / 1000000 + (endhuff.tv_sec - starthuff.tv_sec));
	// printf("Compr Tot =  %6.3f\n", (double)(endcmptot.tv_usec - startcmptot.tv_usec) / 1000000 + (endcmptot.tv_sec - startcmptot.tv_sec));



	gettimeofday(&startdectot, NULL);
	// printf("Decompressing...\n");

	//Huffman Decompression
	gettimeofday(&startihuff, NULL);
	lab4("d", "output.ppm.huf", "output.ppm.uhuf");
	gettimeofday(&endihuff, NULL);


	//End of JPEG encoding



	//////////////////////////////////////////////////////////////////////////////
	//  BEGIN DECODING OF UNCOMPRESSED JPEG IMAGE
	//////////////////////////////////////////////////////////////////////////////

	//Open uncompressed huffman file
	FILE *g_in  = fopen("output.ppm.uhuf", "r");
	FILE *g_out = fopen("result.ppm", "w");

	//Parse header
	fscanf(g_in, "%s\n", img_type);
	fscanf(g_in, "%d %d\n", &orig_col, &orig_row);
	fscanf(g_in, "%d\n", &char_val);

	//Pad row and col if necessary
	col = orig_col;
	row = orig_row;
	if(orig_col % 8 != 0)
	col = orig_col + (8-(orig_col % 8));
	if(orig_row % 8 != 0)
	row = orig_row + (8-(orig_row % 8));
	// printf("%d %d\n",col, row);


	//Prepare huffman array
	memset(huff, 0, row*col*3*sizeof(char));

	// printf("Reading in file\n");
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

	// printf("Reverse Linear\n");
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

	// printf("Before IDCT\n");
	// for(m=0;m<8;m++){
	// 	for(n=0;n<8;n++)
	// 		printf("%6.2f", img_dcr[m*col+n]);
	// 	printf("\n");
	// }

gettimeofday(&startidct, NULL);
	// IDCT
	for(m=0; m<row; m+=8)
		for(n=0; n<col; n+=8)
		{
			for(i = 0; i < 8; i++)
				for(j = 0; j < 8; j++)
				{
					temp_dcty[i*8+j] = 0;
					temp_dctcb[i*8+j] = 0;
					temp_dctcr[i*8+j] = 0;
					for(x = 0; x < 8; x++)
					{
						temp_dcty[i*8+j] += IDCT[i][x] * img_dy[(m+x)*col + (n+j)];
						temp_dctcb[i*8+j] += IDCT[i][x] * img_dcb[(m+x)*col + (n+j)];
						temp_dctcr[i*8+j] += IDCT[i][x] * img_dcr[(m+x)*col + (n+j)];
					}
				}
			for(i = 0; i < 8; i++)
				for(j = 0; j < 8; j++)
				{
					img_y[(m+i) * col + (n+j)] = 0;
					img_cb[(m+i) * col + (n+j)] = 0;
					img_cr[(m+i) * col + (n+j)] = 0;
					for(x = 0; x < 8; x++)
					{
						img_y[(m+i) * col + (n+j)] += (int)(temp_dcty[i*8+x] * IDCT[j][x]);
						img_cb[(m+i) * col + (n+j)] += (int)(temp_dctcb[i*8+x] * IDCT[j][x]);
						img_cr[(m+i) * col + (n+j)] += (int)(temp_dctcr[i*8+x] * IDCT[j][x]);
					}
				}
		}
gettimeofday(&endidct, NULL);

	// printf("After IDCT\n");
	// for(m=0;m<8;m++){
	// 	for(n=0;n<8;n++)
	// 		printf("%6d", img_cr[m*col+n]);
	// 	printf("\n");
	// }

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

	//printf("YCbCr->RGB\n");
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

	//printf("Write out\n");
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
	// printf("YCbCr->RGB = %6.3f\n", (double)(endirgb.tv_usec - startirgb.tv_usec) / 1000000 + (endirgb.tv_sec - startirgb.tv_sec));
	// printf("IDCT =       %6.3f\n", (double)(endidct.tv_usec - startidct.tv_usec) / 1000000 + (endidct.tv_sec - startidct.tv_sec));
	// printf("IQuant =     %6.3f\n", (double)(endiquant.tv_usec - startiquant.tv_usec) / 1000000 + (endiquant.tv_sec - startiquant.tv_sec));
	// printf("IHuffman =   %6.3f\n", (double)(endihuff.tv_usec - startihuff.tv_usec) / 1000000 + (endihuff.tv_sec - startihuff.tv_sec));
	// printf("Decom Tot =  %6.3f\n", (double)(enddectot.tv_usec - startdectot.tv_usec) / 1000000 + (enddectot.tv_sec - startdectot.tv_sec));

	// printf("freeing memory\n");
	//Free allocated memory
	free(img_c);
	free(img_y);
	free(img_cb);
	free(img_cr);
	free(img_dy);
	free(img_dcb);
	free(img_qy);
	free(img_qcb);
	free(img_qcr);
	free(trav_arr_qy);
	free(trav_arr_qcb);
	free(trav_arr_qcr);
	free(huff);

	system("rm output.ppm output.ppm.huf output.ppm.uhuf");

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
