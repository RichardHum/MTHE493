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

//long double origNrg(int** img, int** imgCopy, int rows, int cols);
//long double SiteNrg(int** img, int** imgCopy, int rows, int cols, int x, int y);
//long double origNNrg(int** img, int rows, int cols);
//long double NeighNrg(int** img, int rows, int cols, int x, int y);

double lambda = 10;
double sigma = 0.25;
int** imgCopy1;


int main(int argc, char* argv[])
{
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
    //double initEnergy = origNNrg(img,rows,cols);
    printf("initEnergy: %f\n", initEnergy);
    double initNbhdEnergy = 0;
    double curNbhdEnergy = 0;
    l = 0;

    /*for (i=0; i<rows; i++)
    {
        for (j=0; j<cols; j++)
        {
            initNbhdEnergy = nbhdEnergy(rows,cols,img,origImg,i,j);
            for (k=0;k<256;k++)
            {
                imgCopy[i][j] = k;
                curNbhdEnergy = nbhdEnergy(rows,cols,imgCopy,origImg,i,j);
                current = initEnergy - initNbhdEnergy + curNbhdEnergy;
                locOsc = current > locOsc ? current : locOsc;
            }
            imgCopy[i][j] = img[i][j];
            osc = locOsc > osc ? locOsc : osc;
        }
    }*/

    //delta = osc;
    //printf("Oscillation: %f\n",delta);

    double beta;
    double curEnergy;


    //Main loop
    for (k=1;k<10000;k+=1)
    {
        beta = 3 / log(k+1);
		//beta = (rows * cols * delta)/log(k+1);
        //beta = 1;
        curEnergy = energy(rows,cols,img,origImg);
    	//curEnergy = origNNrg(img,rows,cols);
        for (i=0;i<rows;i++)
        {
            for (j=0;j<cols;j++)
            {
                //img[i][j]=sampleProbability(beta,rows,cols,img,origImg,&curEnergy,j,i);
                img[i][j]=metropolis(beta,rows,cols,img,origImg,&curEnergy,j,i);
            }
        }
		
    }


    //beta = 100;
    //printf("%d ",img[3][3]);
    //img[3][3]=sampleProbability(beta,rows,cols,img,origImg,&curEnergy,3,3);
    //printf("%d\n",img[3][3]);
    

    int count = 0;

    //for (i=0;i<rows;i++)
    //{
    //    for (j=0;j<cols;j++)
    //    {
    //        count = count + origImg[i][j];
    //    }
    //}
    //printf("%d\n",count);

    //int** tempImg;
    //testFcn(rows,cols,img,origImg);

    //img = origImg;

    //for (i=0;i<10;i++)
    //{
    //    testFcn(rows,cols,origImg,img);
    //    tempImg = img;
    //    img = origImg;
    //    origImg = tempImg;
    //}
    ////printf("%d\n",count);

    //double endEnergy = energy(rows,cols,img,imgCopy);
    //printf("endEnergy: %f\n", endEnergy);

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

    printf("pt1: %f pt2: %f\n",pt1,pt2);
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

int sampleProbability(double beta, int rows, int cols, int** img, int** origImg, double* initEnergy, int x, int y)
{
    double probs[256];
    double energies[256];
    double sum = 0;
    double randomNumber;
    double initNbhdEnergy,curNbhdEnergy,current;
    double max;

    int i,j,k;


    for (i=0;i<rows;i++)
    {
        for (j=0;j<cols;j++)
            imgCopy1[i][j] = img[i][j];
    }


    initNbhdEnergy = nbhdEnergy(rows,cols,img,origImg,x,y);
    for (k=0;k<256;k++)
    {
        imgCopy1[y][x] = k;
        curNbhdEnergy = nbhdEnergy(rows,cols,imgCopy1,origImg,x,y);
        current = *initEnergy - initNbhdEnergy + curNbhdEnergy;
        energies[k] = current;
    }

    for (k=0;k<256;k++)
    {
        probs[k] = exp((-1 * (1/beta) * (energies[k] - max)));// - (-1 * (1/beta) * *initEnergy));
        printf("%f %f %f %f %f %f\n",beta,energies[k]-max,probs[k],curNbhdEnergy,initNbhdEnergy,*initEnergy);
        sum += probs[k];
    }

    for (k=0;k<256;k++)
    {
        probs[k] = probs[k] / sum;
    }

    for (k=1;k<256;k++)
    {
 //       printf("%f\n",probs[k]);
        probs[k] = probs[k-1] + probs[k];
    }
    
    randomNumber = rand();
    randomNumber = randomNumber/RAND_MAX;

    k = 0;
    while (1)
    {
        if (randomNumber < probs[k])
            break;
        k++;
    }

    if (k>255)
        printf("k: %d sum: %f probs[255]: %f\t",k,sum,probs[255]);

    *initEnergy = energies[k];
    return k;
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

    int i,j,k;


    for (i=0;i<rows;i++)
    {
        for (j=0;j<cols;j++)
            imgCopy1[i][j] = img[i][j];
    }


    initNbhdEnergy = nbhdEnergy(rows,cols,img,origImg,x,y);
    //initNbhdEnergy = NeighNrg(imgCopy1,rows,cols,x,y);

    randomNumber = rand();
    randomNumber = (randomNumber/RAND_MAX) * 255;
    randInt = randomNumber;

    imgCopy1[y][x] = randInt;
    curNbhdEnergy = nbhdEnergy(rows,cols,imgCopy1,origImg,x,y);
    //curNbhdEnergy = NeighNrg(imgCopy1,rows,cols,x,y);
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



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// calculates energy factor between observed and current image for full image
//long double origNrg(int** img, int** imgCopy, int rows, int cols){
//    long double n=0;
//    int i,j;
//    for (i=0;i<rows;i++){
//            for (j=0;j<cols;j++)
//                n += SiteNrg(img, imgCopy, rows, cols, i, j)/255;
//    }
//    return n;
//}
//// calculates energy factor between and observed image for just specific site
//long double SiteNrg(int** img, int** imgCopy, int rows, int cols, int x, int y){
//    long double n=0;
//    n = pow((img[x][y] - imgCopy[x][y]),2);
//    return n;
//}
//
//// calculates energy factor from neighbouring sites for full image
//long double origNNrg(int** img, int rows, int cols){
//    long double n=0;
//    int i,j;
//    for (i=0;i<rows;i++){
//            for (j=0;j<cols;j++)
//                n += NeighNrg(img, rows, cols, i, j);
//    }
//    return n;
//}
//// calculates energy factor from neighbouring sites for just specific site
//long double NeighNrg(int** img, int rows, int cols, int x, int y){
//    long double n=0;
//    // all cases to take care of edges
//    if(x == 0 && y == 0)
//        n = (img[x+1][y]-img[x][y])*(img[x+1][y]-img[x][y]) + (img[x][y+1]-img[x][y])*(img[x][y+1]-img[x][y]);
//
//    else if(x == 0 && y == cols-1)
//        n = (img[x+1][y]-img[x][y])*(img[x+1][y]-img[x][y]) + (img[x][y-1]-img[x][y])*(img[x][y-1]-img[x][y]);
//
//    else if(x == rows-1 && y == cols-1)
//        n = (img[x-1][y]-img[x][y])*(img[x-1][y]-img[x][y]) + (img[x][y-1]-img[x][y])*(img[x][y-1]-img[x][y]);
//
//    else if(x == rows-1 && y == 0)
//        n = (img[x-1][y]-img[x][y])*(img[x-1][y]-img[x][y]) + (img[x][y+1]-img[x][y])*(img[x][y+1]-img[x][y]); 
//
//    else if(x == 0 && y!=0 && y!=cols-1)
//        n = (img[x+1][y]-img[x][y])*(img[x+1][y]-img[x][y]) + (img[x][y+1]-img[x][y])*(img[x][y+1]-img[x][y]) + (img[x][y-1]-img[x][y])*(img[x][y-1]-img[x][y]); 
//
//    else if(x == rows-1 && y!=0 && y!=cols-1)
//        n = (img[x][y+1]-img[x][y])*(img[x][y+1]-img[x][y]) + (img[x][y-1]-img[x][y])*(img[x][y-1]-img[x][y]) + (img[x-1][y]-img[x][y])*(img[x-1][y]-img[x][y]); 
//
//    else if(y == 0 && x!=0 && x!=rows-1)
//        n = (img[x+1][y]-img[x][y])*(img[x+1][y]-img[x][y]) + (img[x][y+1]-img[x][y])*(img[x][y+1]-img[x][y]) + (img[x-1][y]-img[x][y])*(img[x-1][y]-img[x][y]); 
//
//    else if(y == cols-1 && x!=0 && x!=rows-1)
//        n = (img[x-1][y]-img[x][y])*(img[x-1][y]-img[x][y]) + (img[x+1][y]-img[x][y])*(img[x+1][y]-img[x][y]) + (img[x][y-1]-img[x][y])*(img[x][y-1]-img[x][y]); 
//
//    else
//        n = (img[x-1][y]-img[x][y])*(img[x-1][y]-img[x][y]) + (img[x+1][y]-img[x][y])*(img[x+1][y]-img[x][y]) + (img[x][y-1]-img[x][y])*(img[x][y-1]-img[x][y]) + (img[x][y+1]-img[x][y])*(img[x][y+1]-img[x][y]);
//
//    return n;
//}
