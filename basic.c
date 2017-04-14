#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

//Function prototypes
void Traverse(char *block, char *arr, int row);
void Inverse(char *block, char *arr, int row);

int main (int argc, char **argv)
{
	//Timing variables
	struct timeval startrgb, endrgb, startdct, enddct, startquant, endquant, starthuff, endhuff, starttot, endtot;
	gettimeofday(&starttot, NULL);

	//Declaration of variables
	FILE *f_in  = fopen(argv[1], "r");
	FILE *f_out = fopen("output.ppm", "w");
	char img_type[16];
	int row, col, char_val, orig_row, orig_col;
	int c, i, j, k, m, n, x, y;
	int counter = 0;
	double temp_y, temp_cb, temp_cr;

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
	// double DCT[8][8] =  {       sqrt(2)/4,         sqrt(2)/4,         sqrt(2)/4,         sqrt(2)/4,         sqrt(2)/4,         sqrt(2)/4,         sqrt(2)/4,          sqrt(2)/4,
	// 	    	         		   cos(M_PI/16)/2,  cos(3*M_PI/16)/2,  cos(5*M_PI/16)/2,  cos(7*M_PI/16)/2,  cos(9*M_PI/16)/2, cos(11*M_PI/16)/2, cos(13*M_PI/16)/2,  cos(15*M_PI/16)/2,
	//                      cos(2*M_PI/16)/2,  cos(6*M_PI/16)/2, cos(10*M_PI/16)/2, cos(14*M_PI/16)/2, cos(18*M_PI/16)/2, cos(22*M_PI/16)/2, cos(26*M_PI/16)/2,  cos(30*M_PI/16)/2,
	//                      cos(3*M_PI/16)/2,  cos(9*M_PI/16)/2, cos(15*M_PI/16)/2, cos(21*M_PI/16)/2, cos(27*M_PI/16)/2, cos(33*M_PI/16)/2, cos(39*M_PI/16)/2,  cos(45*M_PI/16)/2,
	//                      cos(4*M_PI/16)/2, cos(12*M_PI/16)/2, cos(20*M_PI/16)/2, cos(28*M_PI/16)/2, cos(36*M_PI/16)/2, cos(44*M_PI/16)/2, cos(52*M_PI/16)/2,  cos(60*M_PI/16)/2,
	//                      cos(5*M_PI/16)/2, cos(15*M_PI/16)/2, cos(25*M_PI/16)/2, cos(35*M_PI/16)/2, cos(45*M_PI/16)/2, cos(55*M_PI/16)/2, cos(65*M_PI/16)/2,  cos(75*M_PI/16)/2,
	//                      cos(6*M_PI/16)/2, cos(18*M_PI/16)/2, cos(30*M_PI/16)/2, cos(42*M_PI/16)/2, cos(54*M_PI/16)/2, cos(66*M_PI/16)/2, cos(78*M_PI/16)/2,  cos(90*M_PI/16)/2,
	//                      cos(7*M_PI/16)/2, cos(21*M_PI/16)/2, cos(35*M_PI/16)/2, cos(49*M_PI/16)/2, cos(63*M_PI/16)/2, cos(77*M_PI/16)/2, cos(91*M_PI/16)/2, cos(105*M_PI/16)/2};

	double DCT[8][8] =  {0.3536,  0.3536,  0.3536,  0.3536,  0.3536,  0.3536,  0.3536,  0.3536,
											 0.4904,  0.4157,  0.2778,  0.0975, -0.0975, -0.2778, -0.4157, -0.4904,
											 0.4619,  0.1913, -0.1913, -0.4619, -0.4619, -0.1913,  0.1913,  0.4619,
											 0.4157, -0.0975, -0.4904, -0.2778,  0.2778,  0.4904,  0.0975, -0.4157,
											 0.3536, -0.3536, -0.3536,  0.3536,  0.3536, -0.3536, -0.3536,  0.3536,
											 0.2778, -0.4904,  0.0975,  0.4157, -0.4157, -0.0975,  0.4904, -0.2778,
											 0.1913, -0.4619,  0.4619, -0.1913, -0.1913,  0.4619, -0.4619,  0.1913,
											 0.0975, -0.2778,  0.4157, -0.4904,  0.4904, -0.4157,  0.2788, -0.0975};

	//Parse header
	fscanf(f_in, "%s\n", img_type);
	fscanf(f_in, "%d %d\n", &orig_col, &orig_row);
	fscanf(f_in, "%d\n", &char_val);

	//Pad row and col if necessary
	row = orig_row + (orig_row % 8);
	col = orig_col + (orig_col % 8);

	//Full RGB matrix
	unsigned char **img_c = (unsigned char **)calloc(row, sizeof(char *));

	//Separate YCbCr matrices
	int **img_y  = (int **)calloc(row, sizeof(int *));
	int **img_cb = (int **)calloc(row, sizeof(int *));
	int **img_cr = (int **)calloc(row, sizeof(int *));

	//Discrete cosine transform matrices
	double **img_dy  = (double **)calloc(row, sizeof(double *));
	double **img_dcb = (double **)calloc(row, sizeof(double *));
	double **img_dcr = (double **)calloc(row, sizeof(double *));

	//Quantization matrices (1D)
	char *img_qy  = (char *)calloc(row*col, sizeof(char));
	char *img_qcb = (char *)calloc(row*col, sizeof(char));
	char *img_qcr = (char *)calloc(row*col, sizeof(char));

	//Temp arrays for Traverse()
	char *trav_arr_qy  = (char *)calloc(64, sizeof(char));
	char *trav_arr_qcb = (char *)calloc(64, sizeof(char));
	char *trav_arr_qcr = (char *)calloc(64, sizeof(char));

	//Rearranged matrix for huffman (1D)
	char *huff = (char *)calloc(row*col*3, sizeof(char));

	//Allocate 2D arrays
	for(i=0; i<col; i++)
	{
		img_c[i]   = (unsigned char *)calloc(col*3, sizeof(char));
	  img_y[i]   = (int *)calloc(col, sizeof(int));
		img_cb[i]  = (int *)calloc(col, sizeof(int));
		img_cr[i]  = (int *)calloc(col, sizeof(int));
		img_dy[i]  = (double *)calloc(col, sizeof(double));
		img_dcb[i] = (double *)calloc(col, sizeof(double));
		img_dcr[i] = (double *)calloc(col, sizeof(double));
	}

	//Read in pixel data
	for(i=0; i<row; i++)
		for(j=0; j<col*3; j++)
			fscanf(f_in, "%c", &img_c[i][j]);

	//RBG -> YCbCr
	gettimeofday(&startrgb, NULL);
	for(i=0; i<row; i++)
		for(j=0,k=0; j<col*3; j+=3,k++)
		{
			img_y[i][k]  = (0.299)*img_c[i][j] + (0.587)*img_c[i][j+1] + (0.114)*img_c[i][j+2];
		  img_cb[i][k] = 128 - (0.168736)*img_c[i][j] - (0.331264)*img_c[i][j+1] + (0.5)*img_c[i][j+2];
		  img_cr[i][k] = 128 + (0.5)*img_c[i][j] - (0.418688)*img_c[i][j+1] - (0.081312)*img_c[i][j+2];
		}

	//Center
	for(i=0; i<row; i++)
		for(j=0; j<col; j++)
		{
			img_y[i][j]  -= 127;
			img_cb[i][j] -= 127;
			img_cr[i][j] -= 127;
		}
	gettimeofday(&endrgb, NULL);

	//Discrete Cosine Transform
	gettimeofday(&startdct, NULL);
	for(m=0; m<row; m+=8)
		for(n=0; n<col; n+=8)
			for(i = 0; i < 8; i++)
				for(j = 0; j < 8; j++)
				{
					temp_y  = 0.0;
					temp_cb = 0.0;
					temp_cr = 0.0;
				  for (x=0; x<8; x++)
				  	for (y=0; y<8; y++)
						{
							temp_y  += DCT[x][i] * DCT[y][j] * img_y[m+x][n+y];
							temp_cb += DCT[x][i] * DCT[y][j] * img_cb[m+x][n+y];
							temp_cr += DCT[x][i] * DCT[y][j] * img_cr[m+x][y+n];
		 				}
					if(i==0)
					{
		 				temp_y  /= sqrt(2);
						temp_cb /= sqrt(2);
						temp_cr /= sqrt(2);
					}
					if(j==0)
					{
						temp_y  /= sqrt(2);
						temp_cb /= sqrt(2);
						temp_cr /= sqrt(2);
					}
					temp_y  /= 4;
					temp_cb /= 4;
					temp_cr /= 4;

				  img_dy[m+i][n+j]  = temp_y;
					img_dcb[m+i][n+j] = temp_cb;
					img_dcr[m+i][n+j] = temp_cr;
				}
	gettimeofday(&enddct, NULL);

	//Quantization
	gettimeofday(&startquant, NULL);
	for(m=0; m<row; m+=8)
		for(n=0; n<col; n+=8)
			for(i=0; i<8; i++)
				for(j=0; j<8; j++)
				{
					img_qy[(m+i)*row + (n+j)]  = (char)rint((img_dy[m+i][n+j]/Q[i][j]));
					img_qcb[(m+i)*row + (n+j)] = (char)rint((img_dcb[m+i][n+j]/Q[i][j]));
					img_qcr[(m+i)*row + (n+j)] = (char)rint((img_dcr[m+i][n+j]/Q[i][j]));
				}
	gettimeofday(&endquant, NULL);

	//Linearization of each 8x8 block before compression
	for(m=0; m<row; m+=8)
		for(n=0; n<col; n+=8)
		{
			//Linearization
			Traverse(img_qy+(m*row+n), trav_arr_qy, row);
			Traverse(img_qcb+(m*row+n), trav_arr_qcb, row);
			Traverse(img_qcr+(m*row+n), trav_arr_qcr, row);

			//Combination into single Huffman array
			for(c=0; c<64; c++, counter++)
			{
				huff[counter] = trav_arr_qy[c];
				huff[col*row+counter] = trav_arr_qcb[c];
				huff[col*row*2+counter] = trav_arr_qcr[c];
			}
		}


	//Write out combined matrix to "output.ppm"
	fprintf(f_out, "%s\n", img_type);
	fprintf(f_out, "%d %d\n", col, row);
	fprintf(f_out, "%d\n", char_val);
	for(m=0; m<row; m++)
		for(n=0; n<col; n++)
			fprintf(f_out, "%c", huff[m*row+n]);

	//Huffman Compression
	gettimeofday(&starthuff, NULL);
	char arg_comp[100] = {"./huff -c output.ppm output.ppm.huf"};
	system(arg_comp);
	gettimeofday(&endhuff, NULL);

	//Huffman Decompression
	char arg_decomp[100] = {"./huff -d output.ppm.huf output.ppm.uhuf"};
	system(arg_decomp);
	char dif[100] = {"diff output.ppm output.ppm.uhuf"};

	//Check to ensure Huffman was successful
	system(dif);
	gettimeofday(&endtot, NULL);

	//Timing calculations
	double delta_us_rgb   = (double)(endrgb.tv_usec - startrgb.tv_usec) / 1000000 + (endrgb.tv_sec - startrgb.tv_sec);
	double delta_us_dct   = (double)(enddct.tv_usec - startdct.tv_usec) / 1000000 + (enddct.tv_sec - startdct.tv_sec);
	double delta_us_quant = (double)(endquant.tv_usec - startquant.tv_usec) / 1000000 + (endquant.tv_sec - startquant.tv_sec);
	double delta_us_huff  = (double)(endhuff.tv_usec - starthuff.tv_usec) / 1000000 + (endhuff.tv_sec - starthuff.tv_sec);
	double delta_us_tot   = (double)(endtot.tv_usec - starttot.tv_usec) / 1000000 + (endtot.tv_sec - starttot.tv_sec);

	//Timing outputs in milliseconds
	printf("RGB->YCbCr = %6.3f\n", delta_us_rgb);
	printf("DCT =        %6.3f\n", delta_us_dct);
	printf("Quant =      %6.3f\n", delta_us_quant);
	printf("Huffman =    %6.3f\n", delta_us_huff);
	printf("Total      = %6.3f\n", delta_us_tot);

	//End of JPEG encoding
	fclose(f_in);
	fclose(f_out);


	//////////////////////////////////////////////////////////////////////////////
	//  BEGIN DECODING OF UNCOMPRESSED JPEG IMAGE
	//////////////////////////////////////////////////////////////////////////////

	//Open uncompressed huffman file
	f_in  = fopen("output.ppm.uhuf", "r");
	f_out = fopen("result.ppm", "w");

	//Parse header
	fscanf(f_in, "%s\n", img_type);
	fscanf(f_in, "%d %d\n", &col, &row);			//Why orig_col and orig_row here?
	fscanf(f_in, "%d\n", &char_val);

	//Prepare huffman array
	memset(huff, 0, row*col*3*sizeof(char));		//Will it clear old memory?

	//Read in file contents
	for(m=0; m<row; m++)
		for(n=0; n<col*3; n++)
			fscanf(f_in, "%c", &huff[m * row + n]);

	//Reverse the linearization
	counter = 0;
	for(i=0; i<3; i++)
		for(m=0; m<row; m+=8)
			for(n=0; n<col; n+=8, counter++)
				if(i==0)
					Inverse(&img_qy[m * row + n], (huff + counter*64), row);
				else if(i==1)
					Inverse(&img_qcb[m * row + n], (huff + counter*64), row);
				else
					Inverse(&img_qcr[m * row + n], (huff + counter*64), row);


	//Inverse Quantization
	for(m=0; m<row; m+=8)
		for(n=0; n<col; n+=8)
			for (i=0; i<8; i++)
				for (j=0; j<8; j++)
				{
							img_dy[m+i][n+j] = rint((img_qy[(m+i)*row + (n+j)]*Q[i][j]));
							img_dcb[m+i][n+j] = rint((img_qcb[(m+i)*row + (n+j)]*Q[i][j]));
							img_dcr[m+i][n+j] = rint((img_qcr[(m+i)*row + (n+j)]*Q[i][j]));
				}

	//Inverse Discrete Cosine Transform
	for(m=0; m<row; m+=8)
		for(n=0; n<col; n+=8)
			for(i=0; i<8; i++)
				for(j=0; j<8; j++)
				{
					temp_y = 0.0;
					temp_cb = 0.0;
					temp_cr = 0.0;
				  for (x=0; x<8; x++)
				  	for (y=0; y<8; y++)
						{
							temp_y  += DCT[x][i] * DCT[y][j] * img_dy[m+x][n+y];
							temp_cb += DCT[x][i] * DCT[y][j] * img_dcb[m+x][n+y];
							temp_cr += DCT[x][i] * DCT[y][j] * img_dcr[m+x][y+n];
		 				}
					if(i != 0)
					{
		 				temp_y  *= 2;
						temp_cb *= 2;
						temp_cr *= 2;
					}
					if(j != 0)
					{
						temp_y  *= 2;
						temp_cb *= 2;
						temp_cr *= 2;
					}

				  img_y[m+i][n+j]  = temp_y;
					img_cb[m+i][n+j] = temp_cb;
					img_cr[m+i][n+j] = temp_cr;
				}

	//Un-Center
	for(i=0; i<row; i++)
		for(j=0; j<col; j++)
		{
			img_y[i][j]  += 127;
			img_cb[i][j] += 127;
			img_cr[i][j] += 127;
		}

	//YCbCr back to RGB
	for(m=0; m<row; m++)
		for(n=0, j=0; n<col; n++, j+=n*3)
		{
			img_c[m][j]   = (char)(img_y[m][n] + 1.40200 * (img_cr[m][n] - 0x80));
			img_c[m][j+1] = (char)(img_y[m][n] - 0.34414 * (img_cb[m][n] - 0x80) - 0.71414 * (img_cr[m][n] - 0x80));
			img_c[m][j+2] = (char)(img_y[m][n] + 1.77200 * (img_cb[m][n] - 0x80));
		}

	//Write out the final, reconstructed RGB image
	fprintf(f_out, "%s\n", img_type);
	fprintf(f_out, "%d %d\n", col, row);
	fprintf(f_out, "%d\n", char_val);
	for(m=0; m<row; m++)
		for(n=0; n<col*3; n++)
			fprintf(f_out, "%c", img_c[m][n]);

	//Clean up
	fclose(f_in);
	fclose(f_out);

	//Free allocated memory
	for(i=0; i<col; i++)
	{
		free(img_c[i]);
		free(img_y[i]);
		free(img_cb[i]);
		free(img_cr[i]);
		free(img_dy[i]);
		free(img_dcb[i]);
		free(img_dcr[i]);
	}

	free(img_c);
	free(img_y);
	free(img_cb);
	free(img_cr);
	free(img_dy);
	free(img_dcb);
	free(img_cr);
	free(img_qy);
	free(img_qcb);
	free(img_qcr);
	free(trav_arr_qy);
	free(trav_arr_qcb);
	free(trav_arr_qcr);
	free(huff);

	//Scott is n00b
	return 0;
}

//////////////////////////////////////////////////////////////////////////////
//  FUNCTION DEFINITIONS of Traverse() and Inverse()
//////////////////////////////////////////////////////////////////////////////
void Traverse(char *block, char *arr, int row)
{
	int count = 0;
	int r = 0;
	int c = 0;

	while(count < 64)
	{
		if(c < 7)
		{
			arr[count++] = block[(r*row) + (c++)];
			if(count == 64)
				break;
		}
		else
			arr[count++] = block[(r++)*row + c];

		while((r<7) && (c>0))
		    arr[count++] = block[(r++)*row + (c--)];

		if(r < 7)
			arr[count++] = block[(r++)*row + c];
		else
			arr[count++] = block[(r*row) + (c++)];

		while((r>0) && (c<7))
			arr[count++] = block[(r--)*row + (c++)];
	}
}

void Inverse(char *block, char *arr, int row)
{
	int count = 0;
	int r = 0;
	int c = 0;

	while(count < 64)
	{
		if(c < 7)
		{
			block[(r*row) + (c++)] = arr[count++];

			if(count == 64)
				break;
		}
		else
			block[(r++)*row + c] = arr[count++];

		while((r<7) && (c>0))
		    block[(r++)*row + (c--)] = arr[count++];

		if(r < 7)
			block[(r++)*row + c] = arr[count++];
		else
			block[(r*row) + (c++)] = arr[count++];

		while((r>0) && (c<7))
			block[(r--)*row + (c++)] = arr[count++];
	}
}
