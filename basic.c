#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>


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

	FILE *f_in  = fopen(argv[1], "r");
	FILE *f_out = fopen("output.ppm", "w");
	char img_type[16];
	int row, col, char_val;
	int i, j, k;  

	//Quantization Matrix
	char Q[64] = {
		16,  11,  10,  16,  24,  40,  51,  61, 
		12,  12,  14,  19,  26,  58,  60,  55,
		14,  13,  16,  24,  40,  57,  69,  56,
		14,  17,  22,  29,  51,  87,  80,  62,
		18,  22,  37,  56,  68, 109, 103,  77,
		24,  35,  55,  64,  81, 104, 113,  92,
		49,  64,  78,  87, 103, 121, 120, 101,
		72,  92,  95,  98, 112, 100, 103,  99};	

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

	//for(i=0; i<64; i++)
	//	printf("B[%d] = (%d).\n", i, (int)B[i]);
	
	fscanf(f_in, "%s\n", img_type);
    fscanf(f_in, "%d %d\n", &row, &col);
    fscanf(f_in, "%d\n", &char_val); 

    int **img_r   = (int **)calloc(row, sizeof(int *));
	int **img_g   = (int **)calloc(row, sizeof(int *));
	int **img_b   = (int **)calloc(row, sizeof(int *));
    char  **img_c   = (char **)calloc(row, sizeof(char *));
    double **img_d = (double **)calloc(row, sizeof(double *));
    for(i=0; i<row; i++)
    {
    	img_r[i] = (int *)calloc(col, sizeof(int));
		img_g[i] = (int *)calloc(col, sizeof(int));
		img_b[i] = (int *)calloc(col, sizeof(int));
    	img_c[i] = (char *)calloc(col*3, sizeof(char));
    	img_d[i] = (double *)calloc(col*3, sizeof(double));
    }

    //READ IMAGE
    for(i=0; i<row; i++)
	    for(j=0; j<col*3; j++)
	         fscanf(f_in, "%c", &img_c[i][j]);

	//RBG -> YCbCr
	for(i=0; i<row; i++)
	    for(j=0,k=0; j<col*3; j+=3,k++)
	    {
	    	img_r[i][k]   = (0.299)*img_c[i][j] + (0.587)*img_c[i][j+1] + (0.114)*img_c[i][j+2];
	    	img_g[i][k] = 128 - (0.168736)*img_c[i][j] - (0.331264)*img_c[i][j+1] + (0.5)*img_c[i][j+2];
	    	img_b[i][k] = 128 + (0.5)*img_c[i][j] - (0.418688)*img_c[i][j+1] - (0.081312)*img_c[i][j+2];
	    }	

	//CENTER
	for(i=0; i<row; i++)
	    for(j=0; j<col; j++){
	    	img_r[i][j] -= 127;
			img_g[i][j] -= 127;
			img_b[i][j] -= 127;
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
