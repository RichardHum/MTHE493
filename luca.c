
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

long double NeighNrg(int** img, int rows, int cols, int x, int y);
long double origNrg(int** img, int** imgCopy, int rows, int cols);
long double origNNrg(int** img, int rows, int cols);
long double SiteNrg(int** img, int** imgCopy, int rows, int cols, int x, int y);
int flip(int** img, int**imgCopy, int rows, int cols, double b, int x, int y);

int main(int argc, char* argv[])
{
    // --------------------------- I/O stuff ---------------------------
    if (argc != 5)
    {
        printf("Invalid number of arguments\n%s [INFILE] [OUTFILE] [COLS] [ROWS]\n",argv[0]);
        return 1;
    }
    int i=0, j=0, k=0, l=0;
    char* fName = argv[1];
    char* oName = argv[2];
    int cols = atoi(argv[3]);
    int rows = atoi(argv[4]);

    size_t count;
    char buffer[5000];
    char *line;
    char *record;

    int** img;
    img = malloc(rows * sizeof *img);
    for (i=0;i<rows;i++)
        img[i] = malloc(cols * sizeof *img[i]);
    int** origImg;
    origImg = malloc(rows * sizeof *origImg);
    for (i=0;i<rows;i++)
        origImg[i] = malloc(cols * sizeof *origImg[i]);
    FILE *file;
    file = fopen(fName,"r");
    i = 0;
    while ((line = fgets(buffer,sizeof(buffer),file)) != NULL){
        record = strtok(line,",");
        j = 0;
        while(record != NULL){
            img[i][j++] = atoi(record);
            record = strtok(NULL,",");
        }
        i++;
    }
    fclose(file);

	// Initialize prng
    srand(time(NULL));

    // copying image
    int** imgCopy;
    imgCopy = malloc(rows * sizeof(*imgCopy));
    for (i=0;i<rows;i++)
        imgCopy[i] = malloc(cols * sizeof(*imgCopy[i]));
    for (i=0;i<rows;i++){
        for (j=0;j<cols;j++)
            imgCopy[i][j] = img[i][j];
    }

    //--------------------------main loop-------------------------
    double b;
    double h;
    long r;
    int p;

    printf("\n");
    // k iterations
	for(r=1;r<100;r=r+1){
        // increase beta as function of r, oscilliation = 100554
        b = 1/(rows * cols * 100554) * log(r*10000) + 0.0000001*r;
        // for printing energy
        //if(p%5 == 1){
        //    h=0;
        //    h = origNrg(img, imgCopy, rows, cols) + origNNrg(img, rows, cols);
        //    printf("H = %f\n", h);
        //}
        //p++;
        // loop to traverse sites
		for (i=0;i<rows;i++){
	        for (j=0;j<cols;j++)
	        	img[i][j]=flip(img, imgCopy, rows, cols, b, i, j);
		}
	}
    // ------------------------More I/O stuff-----------------------
	
    file = fopen(oName,"w");
    for (i=0;i<rows;i++){
        for (j=0;j<cols;j++)
            fprintf(file,"%d,",img[i][j]);
        fprintf(file,"\n");
    }
    fclose(file);
    for (i=0;i<rows;i++)
        free(origImg[i]);
    free(origImg);

    fclose(file);
    for (i=0;i<rows;i++)
        free(img[i]);
    free(img);
    return 0;
}

int flip(int** img, int**imgCopy, int rows, int cols, double b, int x, int y){

	long double z=0;
	long double p[256];
	double randomNumber=0;

    long double oSNrg=0;
    long double oNNrg=0;

    long double sNrg=0;
    long double nNrg=0;

    long double tNrg=0;
    int k;
    int OMG=0;

    // calculate energy for full image
    oSNrg = origNrg(img, imgCopy, rows, cols);
    oNNrg = origNNrg(img, rows, cols);
    OMG = img[x][y];
    // loop through different values for a particular site
	for (k=0;k<256;k++){
        // calculate energy associated with specific site before changing its value
        sNrg = SiteNrg(img, imgCopy, rows, cols, x, y);
        nNrg = NeighNrg(img, rows, cols, x, y);
        // change the site's value
		img[x][y] = k;
        // re-calculate energy using optimization that only current site's neighbouring energy changes
        tNrg = oSNrg - sNrg + SiteNrg(img, imgCopy, rows, cols, x, y) + oNNrg - nNrg + NeighNrg(img, rows, cols, x, y);
        // calculate numerator of gibb's distribution;
        p[k] = exp(-1.0 * b * tNrg);
        // start summing all these numerators for factor of z in gibb's distribution
	    z += p[k];
	}
    // apply the factor of 1/z to all numerator's in gibb's distribution
	for(k=0;k<256;k++)
		p[k] = p[k] / z; 

    // creating a CDF
	for (k=1;k<256;k++)
        p[k] = p[k-1] + p[k];

    // generate a random number uniform(0,1)
    randomNumber = rand();
    randomNumber = randomNumber/RAND_MAX;
    // recurse through CDF array until you reach the right k value
    k = 0;
    while (1){
        if (randomNumber < p[k])
            break;
        if (k == 255){
            //k = OMG;
            break;
        } 
        k++;
    }
    return k;
}
// calculates energy factor from neighbouring sites for full image
long double origNNrg(int** img, int rows, int cols){
    long double n=0;
    int i,j;
    for (i=0;i<rows;i++){
            for (j=0;j<cols;j++)
                n += NeighNrg(img, rows, cols, i, j);
    }
    return n;
}
// calculates energy factor between observed and current image for full image
long double origNrg(int** img, int** imgCopy, int rows, int cols){
    long double n=0;
    int i,j;
    for (i=0;i<rows;i++){
            for (j=0;j<cols;j++)
                n += (img[i][j] - imgCopy[i][j]) * (img[i][j] - imgCopy[i][j]);
    }
    return n;
}
// calculates energy factor between and observed image for just specific site
long double SiteNrg(int** img, int** imgCopy, int rows, int cols, int x, int y){
    double n=0;
    n = (img[x][y] - imgCopy[x][y]) * (img[x][y] - imgCopy[x][y]);
    return n;
}
// calculates energy factor from neighbouring sites for just specific site
long double NeighNrg(int** img, int rows, int cols, int x, int y){
    long double n=0;
    // all cases to take care of edges
    if(x == 0 && y == 0)
        n += (img[x+1][y]-img[x][y])*(img[x+1][y]-img[x][y]) + (img[x][y+1]-img[x][y])*(img[x][y+1]-img[x][y]);

    else if(x == 0 && y == cols-1)
        n += (img[x+1][y]-img[x][y])*(img[x+1][y]-img[x][y]) + (img[x][y-1]-img[x][y])*(img[x][y-1]-img[x][y]);

    else if(x == rows-1 && y == cols-1)
        n += (img[x-1][y]-img[x][y])*(img[x-1][y]-img[x][y]) + (img[x][y-1]-img[x][y])*(img[x][y-1]-img[x][y]);

    else if(x == rows-1 && y == 0)
        n += (img[x-1][y]-img[x][y])*(img[x-1][y]-img[x][y]) + (img[x][y+1]-img[x][y])*(img[x][y+1]-img[x][y]); 

    else if(x == 0 && y!=0 && y!=cols-1)
        n += (img[x+1][y]-img[x][y])*(img[x+1][y]-img[x][y]) + (img[x][y+1]-img[x][y])*(img[x][y+1]-img[x][y]) + (img[x][y-1]-img[x][y])*(img[x][y-1]-img[x][y]); 

    else if(x == rows-1 && y!=0 && y!=cols-1)
        n += (img[x][y+1]-img[x][y])*(img[x][y+1]-img[x][y]) + (img[x][y-1]-img[x][y])*(img[x][y-1]-img[x][y]) + (img[x-1][y]-img[x][y])*(img[x-1][y]-img[x][y]); 

    else if(y == 0 && x!=0 && x!=rows-1)
        n += (img[x+1][y]-img[x][y])*(img[x+1][y]-img[x][y]) + (img[x][y+1]-img[x][y])*(img[x][y+1]-img[x][y]) + (img[x-1][y]-img[x][y])*(img[x-1][y]-img[x][y]); 

    else if(y == cols-1 && x!=0 && x!=rows-1)
        n += (img[x-1][y]-img[x][y])*(img[x-1][y]-img[x][y]) + (img[x+1][y]-img[x][y])*(img[x+1][y]-img[x][y]) + (img[x][y-1]-img[x][y])*(img[x][y-1]-img[x][y]); 

    else
        n += (img[x-1][y]-img[x][y])*(img[x-1][y]-img[x][y]) + (img[x+1][y]-img[x][y])*(img[x+1][y]-img[x][y]) + (img[x][y-1]-img[x][y])*(img[x][y-1]-img[x][y]) + (img[x][y+1]-img[x][y])*(img[x][y+1]-img[x][y]);
    
    // ising model generalized to more than binary
    /*if(img[x-1][y] == img[x][y])
        n += 1;
    if(img[x+1][y] == img[x][y])
        n += 1;
    if(img[x][y+1] == img[x][y])
        n += 1;
    if(img[x][y-1] == img[x][y])
        n += 1;*/

    return n;
}
