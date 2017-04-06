#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>


void Traverse(unsigned char *block, unsigned char *arr)
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
			arr[count++] = block[(r++)*8 + c];
		//MOVE DOWN AND LEFT
		while((r<7) && (c>0))
		    arr[count++] = block[(r++)*8 + (c--)];

		//SECOND HALF OF TRAVERSE
		if(r < 7) //MOVE DOWN
			arr[count++] = block[(r++)*8 + c];
		else //MOVE RIGHT
			arr[count++] = block[r*8 + (c++)];
		//MOVE UP AND RIGHT
		while((r>0) && (c<7))
			arr[count++] = block[(r--)*8 + (c++)];
	}
}



int main (int argc, char **argv)
{
//	int imageWidth, imageHeight, matSize;
//	int *hostInputImage, *hostOutputImage;
//	int *deviceInputImage, *deviceOutputImage;

	FILE *f_in  = fopen(argv[1], "r");
	FILE *f_out = fopen("output.ppm", "w");
	char img_type[16];
	int row, col, char_val, orig_row, orig_col;
	int i, j, k, x, y;

	//Quantization Matrix
	unsigned char Q[64] = {	16,  11,  10,  16,  24,  40,  51,  61,
									12,  12,  14,  19,  26,  58,  60,  55,
									14,  13,  16,  24,  40,  57,  69,  56,
									14,  17,  22,  29,  51,  87,  80,  62,
									18,  22,  37,  56,  68, 109, 103,  77,
									24,  35,  55,  64,  81, 104, 113,  92,
									49,  64,  78,  87, 103, 121, 120, 101,
									72,  92,  95,  98, 112, 100, 103,  99};

	double DCT[64] =  {         sqrt(2)/4,         sqrt(2)/4,         sqrt(2)/4,         sqrt(2)/4,         sqrt(2)/4,         sqrt(2)/4,         sqrt(2)/4,          sqrt(2)/4,
		    	         		 cos(M_PI/16)/2,  cos(3*M_PI/16)/2,  cos(5*M_PI/16)/2,  cos(7*M_PI/16)/2,  cos(9*M_PI/16)/2, cos(11*M_PI/16)/2, cos(13*M_PI/16)/2,  cos(15*M_PI/16)/2,
	                   cos(2*M_PI/16)/2,  cos(6*M_PI/16)/2, cos(10*M_PI/16)/2, cos(14*M_PI/16)/2, cos(18*M_PI/16)/2, cos(22*M_PI/16)/2, cos(26*M_PI/16)/2,  cos(30*M_PI/16)/2,
	                   cos(3*M_PI/16)/2,  cos(9*M_PI/16)/2, cos(15*M_PI/16)/2, cos(21*M_PI/16)/2, cos(27*M_PI/16)/2, cos(33*M_PI/16)/2, cos(39*M_PI/16)/2,  cos(45*M_PI/16)/2,
	                   cos(4*M_PI/16)/2, cos(12*M_PI/16)/2, cos(20*M_PI/16)/2, cos(28*M_PI/16)/2, cos(36*M_PI/16)/2, cos(44*M_PI/16)/2, cos(52*M_PI/16)/2,  cos(60*M_PI/16)/2,
	                   cos(5*M_PI/16)/2, cos(15*M_PI/16)/2, cos(25*M_PI/16)/2, cos(35*M_PI/16)/2, cos(45*M_PI/16)/2, cos(55*M_PI/16)/2, cos(65*M_PI/16)/2,  cos(75*M_PI/16)/2,
	                   cos(6*M_PI/16)/2, cos(18*M_PI/16)/2, cos(30*M_PI/16)/2, cos(42*M_PI/16)/2, cos(54*M_PI/16)/2, cos(66*M_PI/16)/2, cos(78*M_PI/16)/2,  cos(90*M_PI/16)/2,
	                   cos(7*M_PI/16)/2, cos(21*M_PI/16)/2, cos(35*M_PI/16)/2, cos(49*M_PI/16)/2, cos(63*M_PI/16)/2, cos(77*M_PI/16)/2, cos(91*M_PI/16)/2, cos(105*M_PI/16)/2};

	unsigned char *A = (char *)calloc(64, sizeof(char));
	unsigned char *B = (char *)calloc(64, sizeof(char));

	//memcpy(A, test, 64*sizeof(char));

	//Traverse(A, B);

	//for(i=0; i<64; i++)
	//	printf("B[%d] = (%d).\n", i, (int)B[i]);

	  fscanf(f_in, "%s\n", img_type);
    fscanf(f_in, "%d %d\n", &orig_row, &orig_col);
    fscanf(f_in, "%d\n", &char_val);

		row = orig_row+orig_row%8;
		col = orig_col+orig_col%8;
		//full rgb matrix
		unsigned char **img_c   = (unsigned char **)calloc(row, sizeof(char *));
		// separate ycbcr matrices
    int **img_y   = (int **)calloc(row, sizeof(int *));
	  int **img_cb   = (int **)calloc(row, sizeof(int *));
	  int **img_cr   = (int **)calloc(row, sizeof(int *));
		// discrete cosine transform matrices
    double **img_dy = (double **)calloc(row, sizeof(double *));
		double **img_dcb = (double **)calloc(row, sizeof(double *));
		double **img_dcr = (double **)calloc(row, sizeof(double *));
		// quantization matrices
		unsigned char **img_qy = (unsigned char **)calloc(row, sizeof(char *));
		unsigned char **img_qcb = (unsigned char **)calloc(row, sizeof(char *));
		unsigned char **img_qcr = (unsigned char **)calloc(row, sizeof(char *));
    for(i=0; i<row; i++)
    {
    	img_c[i] = (unsigned char *)calloc(col*3, sizeof(char));
    	img_y[i] = (int *)calloc(col, sizeof(int));
		  img_cb[i] = (int *)calloc(col, sizeof(int));
		  img_cr[i] = (int *)calloc(col, sizeof(int));
    	img_dy[i] = (double *)calloc(col, sizeof(double));
			img_dcb[i] = (double *)calloc(col, sizeof(double));
			img_dcr[i] = (double *)calloc(col, sizeof(double));
			img_qy[i] = (unsigned char *)calloc(col, sizeof(char));
			img_qcb[i] = (unsigned char *)calloc(col, sizeof(char));
			img_qcr[i] = (unsigned char *)calloc(col, sizeof(char));
    }

    //READ IMAGE
    for(i=0; i<row; i++)
	    for(j=0; j<col*3; j++)
	         fscanf(f_in, "%c", &img_c[i][j]);

	//RBG -> YCbCr
	for(i=0; i<row; i++)
	    for(j=0,k=0; j<col*3; j+=3,k++)
	    {
	    	img_y[i][k]   = (0.299)*img_c[i][j] + (0.587)*img_c[i][j+1] + (0.114)*img_c[i][j+2];
	    	img_cb[i][k] = 128 - (0.168736)*img_c[i][j] - (0.331264)*img_c[i][j+1] + (0.5)*img_c[i][j+2];
	    	img_cr[i][k] = 128 + (0.5)*img_c[i][j] - (0.418688)*img_c[i][j+1] - (0.081312)*img_c[i][j+2];
	    }

	//CENTER
	for(i=0; i<row; i++)
	    for(j=0; j<col; j++){
	    	img_y[i][j] -= 127;
			  img_cb[i][j] -= 127;
			  img_cr[i][j] -= 127;
		}

	int m,n;
	for(m=0; m<row; m+=8){
		for(n=0; n<col; n+=8){
			  // Discrete Cosine Transform
				double temp_y, temp_cb, temp_cr;
				for (i = 0; i < 8; i++) {
			    for (j = 0; j < 8; j++) {
			        temp_y = 0.0;
							temp_cb = 0.0;
							temp_cr = 0.0;
			        for (x = 0; x < 8; x++) {
			            for (y = 0; y < 8; y++) {
										temp_y += cos((2*x+1)*i*M_PI/16) * cos((2*y+1)*j*M_PI/16) * img_y[m+x][n+y];
										temp_cb += cos((2*x+1)*i*M_PI/16) * cos((2*y+1)*j*M_PI/16) * img_cb[m+x][n+y];
										temp_cr += cos((2*x+1)*i*M_PI/16) * cos((2*y+1)*j*M_PI/16) * img_cr[m+x][yn+];
	 								}
							}
							if(i==0){
	 							temp_y /= sqrt(2);
								temp_cb /= sqrt(2);
								temp_cr /= sqrt(2);
							}
							if(j==0){
								temp_y /= sqrt(2);
								temp_cb /= sqrt(2);
								temp_cr /= sqrt(2);
							}
							temp_y /= 4;
							temp_cb /= 4;
							temp_cr /= 4;

			        img_dy[m+i][n+j] = temp_y;
							img_dcb[m+i][n+j] = temp_cb;
							img_dcr[m+i][n+j] = temp_cr;
							//Quantization
							img_qy[m+i][n+j] = (unsigned char)(img_dy[m+i][n+j]/Q[i][j]);
							img_qcb[m+i][n+j] = (unsigned char)(img_dcb[m+i][n+j]/Q[i][j]);
							img_qcr[m+i][n+j] = (unsigned char)(img_dcr[m+i][n+j]/Q[i][j]);
						}
			    }
			  }
			}
	fprintf(f_out, "%s\n", img_type);
    fprintf(f_out, "%d %d\n", row, col);
    fprintf(f_out, "%d\n", char_val);
	for(i=0; i<row; i++)
	    for(j=0; j<col; j++){
	         fprintf(f_out, "%c", img_r[i][j]);
	         fprintf(f_out, "%c", img_g[i][j]);
	         fprintf(f_out, "%c", img_b[i][j]);
		}

	//DISCRETE COSINE TRANSFORM




	return 0;
}
