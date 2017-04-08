#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

double absf(double x);
double min(double x, double y);
double energy(int rows, int cols, int** img, int** origImg);
double nbhdEnergy(int rows, int cols, int** img, int** origImg, int x, int y);
int sampleProbability(double beta, int rows, int cols, int** img, int** origImg, double* initEnergy, int x, int y);
int metropolis(double beta, int rows, int cols, int** img, int** origImg, double* initEnergy, int x, int y);
int testFcn(int rows, int cols, int** imgIn, int** imgOut);

double lambda = 1000;
double sigma = 1000;
int** imgCopy1;


int main(int argc, char* argv[])
{
    if (argc != 4)
    {
        printf("Invalid number of arguments\n%s [BASENAME] [COLS] [ROWS]\n",argv[0]);
        return 1;
    }
    int i=0, j=0, k=0, l=0;
    char fName[255];
    char oName[255];
    char suffix[255];
    char* basename = argv[1];
    int cols = atoi(argv[2]);
    int rows = atoi(argv[3]);

    //size_t count;
    char buffer[5000];
    char *line;
    char *record;

    int** img;
    img = malloc(rows * sizeof *img);
    for (i=0;i<rows;i++)
    {
        img[i] = malloc(cols * sizeof *img[i]);
    }

    int** origImg;
    origImg = malloc(rows * sizeof *origImg);
    for (i=0;i<rows;i++)
    {
        origImg[i] = malloc(cols * sizeof *origImg[i]);
    }

    strcpy(fName,basename);
    strcat(fName,"-noisy.csv");
    FILE *file;
    file = fopen(fName,"r");

    i = 0;
    while ((line = fgets(buffer,sizeof(buffer),file)) != NULL)
    {
        record = strtok(line,",");
        j = 0;
        while(record != NULL)
        {
            img[i][j++] = atoi(record);
            record = strtok(NULL,",");
        }
        i++;
    }
    fclose(file);

//////////////////////////////////////////////////////////////////////////////////

    // Initialize prng
    srand(time(NULL));
    // Calculate oscilation
    double osc = 0;
    double locOsc = 0;
    double current = 0;
    double delta;

    int** imgCopy;
    imgCopy = malloc(rows * sizeof(*imgCopy));
    for (i=0;i<rows;i++)
    {
        imgCopy[i] = malloc(cols * sizeof(*imgCopy[i]));
    }

    imgCopy1 = malloc(rows * sizeof(*imgCopy1));
    for (i=0;i<rows;i++)
    {
        imgCopy1[i] = malloc(cols * sizeof(*imgCopy1[i]));
    }


    for (i=0;i<rows;i++)
    {
        for (j=0;j<cols;j++)
            imgCopy[i][j] = img[i][j];
    }
	
	for (i=0;i<rows;i++)
    {
        for (j=0;j<cols;j++)
            origImg[i][j] = img[i][j];
    }

    double initEnergy = energy(rows,cols,img,origImg);
    printf("initEnergy: %f\n", initEnergy);
    double initNbhdEnergy = 0;
    double curNbhdEnergy = 0;
    l = 0;

    //printf("Oscillation: %f\n",delta);

    double beta;
    double curEnergy;

    int mod;

    float** lambdas;
    lambdas = malloc(noFilts *sizeof(*lambdas));
    for(i=0;i<noFilts;i++)
        lambdas[i] = malloc(8*sizeof(*lambdas[i]));

    float** lambdaTemp[][] = 
        {{0.0440000000000000,0.438800000000000,-0.372400000000000,0.398400000000000,-0.422000000000000,0.000400000000000000}
        {0,-0.00600000000000000,0.00520000000000000,-0.0108000000000000,0.00520000000000000,0.00240000000000000}
        {0,0.0212000000000000,0.00360000000000000,0.0168000000000000,0.00200000000000000,0.0176000000000000}
        {0,-0.00480000000000000,0.0120000000000000,-0.00320000000000000,0.0172000000000000,0.00720000000000000}
        {0,-0.000800000000000000,0.00800000000000000,-0.000400000000000000,-0.00520000000000000,0.00800000000000000}
        {0,0,0.00360000000000000,-0.00280000000000000,0.0100000000000000,0.0108000000000000}
        {0,0,0.469600000000000,0.00200000000000000,0.00640000000000000,0.0100000000000000}
        {-0.0439999999999999,-0.448400000000000,-0.129600000000000,-0.400000000000000,0.386400000000000,-0.0564000000000000}};

    for (i=0;i<noFilts;i++)
    {
        for (j=0;j<8;j++)
        {
            lambdas[i][j] = lambdaTemp[i][j];
        }
    }

    //Main loop
    for (k=1;k<=10000;k+=1)
    {
        //beta = 3 / log(k+1);
		beta = (rows * cols * delta)/log(k+1);
        //beta = 1;
        curEnergy = energy(rows,cols,img,origImg);
        for (i=0;i<rows;i++)
        {
            for (j=0;j<cols;j++)
            {
                img[i][j]=metropolis(beta,rows,cols,img,origImg,&curEnergy,j,i);
                imgCopy1[i][j] = img[i][j];
            }
        }

        mod = 100;
        if (k%mod == 0)
        {
            strcpy(oName,basename);
            sprintf(suffix,"-smooth-%d.csv",k/mod);
            strcat(oName,suffix);
            file = fopen(oName,"w");

            for (i=0;i<rows;i++)
            {
                for (j=0;j<cols;j++)
                    fprintf(file,"%d,",img[i][j]);
                fprintf(file,"\n");
            }
            fclose(file);
        }
    }

    for (i=0;i<rows;i++)
    {
        free(imgCopy[i]);
    }
    free(imgCopy);

    for (i=0;i<rows;i++)
    {
        free(imgCopy1[i]);
    }
    free(imgCopy1);


/////////////////////////////////////////////////////////////////////////////////

    strcpy(oName,basename);
    strcat(oName,"-smooth.csv");
    file = fopen(oName,"w");

    for (i=0;i<rows;i++)
    {
        for (j=0;j<cols;j++)
            fprintf(file,"%d,",img[i][j]);
        fprintf(file,"\n");
    }

    fclose(file);

    for (i=0;i<rows;i++)
    {
        free(origImg[i]);
    }
    free(origImg);

    fclose(file);
    for (i=0;i<rows;i++)
    {
        free(img[i]);
    }
    free(img);
    return 0;
}

int testFcn(int rows, int cols, int** imgIn, int** imgOut)
{
    int count = 0;
    int i,j;

    for (i=1;i<rows-1;i++)
    {
        for (j=1;j<cols-1;j++)
        {
            count = count + imgIn[i][j];
            imgOut[i][j] = (imgIn[i][j] + imgIn[(i+1)][j] + imgIn[(i-1)][j] + imgIn[i][(j+1)] + imgIn[i][(j-1)])/5;
        }
    }

    for (i=1;i<rows-1;i++)
    {
        imgOut[i][0] = (imgIn[i][0] + imgIn[(i+1)][0] + imgIn[(i-1)][0] + imgIn[i][(1)])/4;
        imgOut[i][cols-1] = (imgIn[i][cols-1] + imgIn[(i+1)][cols-1] + imgIn[(i-1)][cols-1] + imgIn[i][(cols-2)])/4;
    }

    for (j=1;j<cols-1;j++)
    {
        imgOut[0][j] = (imgIn[0][j] + imgIn[1][j] + imgIn[0][(j+1)] + imgIn[0][(j-1)])/4;
        imgOut[rows-1][j] = (imgIn[rows-1][j] + imgIn[rows-2][j] + imgIn[rows-1][(j+1)] + imgIn[rows-1][(j-1)])/4;
    }

    return count;
}

double energy(int rows, int cols, int** img, int** origImg)
{
    int i,j,k;


    double ans = 0;

    double pt1 = 0;
    double pt2 = 0;

    double diff = 0;

    for (i=0;i<rows;i++)
    {
        for (j=0;j<cols;j++)
        {
            diff = origImg[i][j] - img[i][j];
            pt1 += abs(diff) * abs(diff) / (255.0*255.0);
        }
    }
    
    pt1 = pt1 * 1/(2.0 * sigma * sigma);
    
    for (i=0;i<rows;i++)
    {
        for (j=0;j<cols;j++)
        {
            if (j+1 < cols)
            {
                diff = img[i][j] - img[i][j+1];
                pt2 += diff * diff/(255.0*255.0);
            }
            if (j-1 >=  0)
            {
                diff = img[i][j] - img[i][j-1];
                pt2 += diff * diff/(255.0*255.0);
            }
            if (i+1 < rows)
            {
                diff = img[i][j] - img[i+1][j];
                pt2 += diff * diff/(255.0*255.0);
            }
            if (i-1 >= 0)
            {
                diff = img[i][j] - img[i-1][j];
                pt2 += diff * diff/(255.0*255.0);
            }

        }
    }

    pt2 = pt2 * lambda;

    ans = pt1 + pt2;

    //printf("pt1: %f pt2: %f\n",pt1,pt2);
    return ans;
}

double nbhdEnergy(int rows, int cols, int** img, int** origImg, int x, int y)
{
    int i,j;
    int xCoord, yCoord;

    double pt1 = 0;
    double pt2 = 0;
    double ans = 0;

    double diff = 0;


    for (i=-1;i<=1;i++)
    {
        for (j=-1;j<=1;j++)
        {
            yCoord = y+i;
            xCoord = x+j;

            if (xCoord < cols && xCoord >= 0 && yCoord < rows && yCoord >= 0)
            {
                pt1 += abs(origImg[yCoord][xCoord] - img[yCoord][xCoord]) * abs(origImg[yCoord][xCoord] - img[yCoord][xCoord])/(255.0*255.0);
				
            }
        }
    }
    pt1 = pt1 * 1/(2.0 * sigma * sigma);

    for (i=-1;i<=1;i++)
    {
        for (j=-1;j<=1;j++)
        {
            yCoord = y+i;
            xCoord = x+j;

            if (xCoord < cols && xCoord >= 0 && yCoord < rows && yCoord >= 0)
            {
				if (xCoord+1 < cols)
                {
                    diff = img[yCoord][xCoord] - img[yCoord][xCoord+1];
                    pt2 += diff * diff/(255.0*255.0);
                }
                if (xCoord-1 >= 0)
                {
                    diff = img[yCoord][xCoord] - img[yCoord][xCoord-1];
                    pt2 += diff * diff/(255.0*255.0);
                }
                if (yCoord+1 < rows)
                {
                    diff = img[yCoord][xCoord] - img[yCoord+1][xCoord];
                    pt2 += diff * diff/(255.0*255.0);
                }
                if (yCoord-1 >= 0)
                {
                    diff = img[yCoord][xCoord] - img[yCoord-1][xCoord];
                    pt2 += diff * diff/(255.0*255.0);
                }
            }
        }
    }

    pt2 = pt2 * lambda;
    ans = pt1 + pt2;
    return ans;
}

int metropolis(double beta, int rows, int cols, int** img, int** origImg, double* initEnergy, int x, int y)
{
    double probs[256];
    double energies[256];
    double sum = 0;
    double randomNumber;
    double initNbhdEnergy,curNbhdEnergy,current;
    double max;
    double prob;
    int randInt;
    int backup;

    int i,j,k;


    /*for (i=0;i<rows;i++)
    {
        for (j=0;j<cols;j++)
            imgCopy1[i][j] = img[i][j];
    }
    */


    initNbhdEnergy = nbhdEnergy(rows,cols,img,origImg,x,y);

    randomNumber = rand();
    randomNumber = (randomNumber/RAND_MAX) * 255;
    randInt = randomNumber;

    imgCopy1[y][x] = randInt;
    curNbhdEnergy = nbhdEnergy(rows,cols,imgCopy1,origImg,x,y);
    current = *initEnergy - initNbhdEnergy + curNbhdEnergy;

    if (current < *initEnergy)
    {
        *initEnergy = current;
        return randInt;
    }
    else
    {
        randomNumber = rand();
        randomNumber = randomNumber/RAND_MAX;
        prob = exp(-(1/beta) * (current - *initEnergy));
        printf("%f %f %f\n ",beta,prob,(current-*initEnergy));
        //printf("%f %f %f %f %f %f\n",beta,prob,curNbhdEnergy,initNbhdEnergy,*initEnergy,current);
        if (randomNumber < prob)
        {
            *initEnergy = current;
            return randInt;
        }
        else
        {
            return img[y][x];
        }
    }
}


double absf(double x)
{
    return x >= 0 ? x : -1 * x;
}

double min(double x, double y)
{
    return x > y ? x : y;
}

