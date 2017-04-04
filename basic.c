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
	char Q[64] = {	16,  11,  10,  16,  24,  40,  51,  61, 
					12,  12,  14,  19,  26,  58,  60,  55,
					14,  13,  16,  24,  40,  57,  69,  56,
					14,  17,  22,  29,  51,  87,  80,  62,
					18,  22,  37,  56,  68, 109, 103,  77,
					24,  35,  55,  64,  81, 104, 113,  92,
					49,  64,  78,  87, 103, 121, 120, 101,
					72,  92,  95,  98, 112, 100, 103,  99};	

	char DCT[64] =  {         sqrt(2)/4,         sqrt(2)/4,         sqrt(2)/4,         sqrt(2)/4,         sqrt(2)/4,         sqrt(2)/4,         sqrt(2)/4,          sqrt(2)/4,
		    	         cos(M_PI/16)/2,  cos(3*M_PI/16)/2,  cos(5*M_PI/16)/2,  cos(7*M_PI/16)/2,  cos(9*M_PI/16)/2, cos(11*M_PI/16)/2, cos(13*M_PI/16)/2,  cos(15*M_PI/16)/2,
	                   cos(2*M_PI/16)/2,  cos(6*M_PI/16)/2, cos(10*M_PI/16)/2, cos(14*M_PI/16)/2, cos(18*M_PI/16)/2, cos(22*M_PI/16)/2, cos(26*M_PI/16)/2,  cos(30*M_PI/16)/2,
	                   cos(3*M_PI/16)/2,  cos(9*M_PI/16)/2, cos(15*M_PI/16)/2, cos(21*M_PI/16)/2, cos(27*M_PI/16)/2, cos(33*M_PI/16)/2, cos(39*M_PI/16)/2,  cos(45*M_PI/16)/2,
	                   cos(4*M_PI/16)/2, cos(12*M_PI/16)/2, cos(20*M_PI/16)/2, cos(28*M_PI/16)/2, cos(36*M_PI/16)/2, cos(44*M_PI/16)/2, cos(52*M_PI/16)/2,  cos(60*M_PI/16)/2,
	                   cos(5*M_PI/16)/2, cos(15*M_PI/16)/2, cos(25*M_PI/16)/2, cos(35*M_PI/16)/2, cos(45*M_PI/16)/2, cos(55*M_PI/16)/2, cos(65*M_PI/16)/2,  cos(75*M_PI/16)/2,
	                   cos(6*M_PI/16)/2, cos(18*M_PI/16)/2, cos(30*M_PI/16)/2, cos(42*M_PI/16)/2, cos(54*M_PI/16)/2, cos(66*M_PI/16)/2, cos(78*M_PI/16)/2,  cos(90*M_PI/16)/2,
	                   cos(7*M_PI/16)/2, cos(21*M_PI/16)/2, cos(35*M_PI/16)/2, cos(49*M_PI/16)/2, cos(63*M_PI/16)/2, cos(77*M_PI/16)/2, cos(91*M_PI/16)/2, cos(105*M_PI/16)/2};

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
