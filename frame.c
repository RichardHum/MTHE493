#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

void localConvolve(int imgRows, int imgCols, int filtRows, int filtCols, int** imgIn, int** imgOut, float** filter,int x, int y);
void convolve(int imgRows, int imgCols, int filtRows, int filtCols, int** imgIn, int** imgOut, float** filter);
void computeHistogram(int imgRows, int imgCols, int** img, int* histogram);
double energy(int noFilts, float lambdas[][256], int* histogram);

int** imgCopy1;

int main(int argc, char* argv[])
{
    if (argc != 4)
    {
        printf("Invalid number of arguments\n%s [BASENAME] [COLS] [ROWS]\n",argv[0]);
        return 1;
    }
    int i=0, j=0, k=0, l=0, m=0;
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
    strcat(fName,".csv");
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

    int** filteredImg;
    filteredImg = malloc(rows * sizeof(*filteredImg));
    for (i=0;i<rows;i++)
    {
        filteredImg[i] = malloc(cols * sizeof(*filteredImg[i]));
    }

    int** synthImg;
    synthImg = malloc(rows * sizeof(*synthImg));
    for (i=0;i<rows;i++)
    {
        synthImg[i] = malloc(cols * sizeof(*synthImg[i]));
    }

    imgCopy1 = malloc(rows * sizeof(*imgCopy1));
    for (i=0;i<rows;i++)
    {
        imgCopy1[i] = malloc(cols * sizeof(*imgCopy1[i]));
    }

    for (i=0;i<rows;i++)
    {
        for (j=0;j<cols;j++)
            filteredImg[i][j] = img[i][j];
    }
	
	for (i=0;i<rows;i++)
    {
        for (j=0;j<cols;j++)
            origImg[i][j] = img[i][j];
    }

	for (i=0;i<rows;i++)
    {
        for (j=0;j<cols;j++)
            imgCopy1[i][j] = img[i][j];
    }
////////////////////////////////////////////////////////////////////////////////

    int filtRows = 5;
    int filtCols = 5;
    int noFilts = 1;
    long epsilon = 255*rows*cols;
    //long epsilon = 700;
    int delta = 1;
    long diff;
    double sum;
    double prob;
    float** filter;
    float curEnergy;
    float beta;

    double probs[256];
    double energies[256];
    double randomNumber;
    double initNbhdEnergy,curNbhdEnergy,current;
    int randInt;
    int backup;

    filter = malloc(filtRows * sizeof(*filter));
    for (i=0;i<filtRows;i++)
    {
        filter[i] = malloc(filtCols * sizeof(*filter[i]));
    }

    //Gaussian
    //filter[0][0] = (1/273.0);
    //filter[0][1] = (4/273.0);
    //filter[0][2] = (7/273.0);
    //filter[0][3] = (4/273.0);
    //filter[0][4] = (1/273.0);
    //filter[1][0] = (4/273.0);
    //filter[1][1] = (16/273.0);
    //filter[1][2] = (26/273.0);
    //filter[1][3] = (16/273.0);
    //filter[1][4] = (4/273.0);
    //filter[2][0] = (7/273.0);
    //filter[2][1] = (26/273.0);
    //filter[2][2] = (41/273.0);
    //filter[2][3] = (26/273.0);
    //filter[2][4] = (7/273.0);
    //filter[3][0] = (4/273.0);
    //filter[3][1] = (16/273.0);
    //filter[3][2] = (26/273.0);
    //filter[3][3] = (16/273.0);
    //filter[3][4] = (4/273.0);
    //filter[4][0] = (1/273.0);
    //filter[4][1] = (4/273.0);
    //filter[4][2] = (7/273.0);
    //filter[4][3] = (4/273.0);
    //filter[4][4] = (1/273.0);

    //Gabor
	filter[0][0] = 0.158635;
    filter[0][1] = 0.757461;
    filter[0][2] = 0.193165;
    filter[0][3] = -0.668173;
    filter[0][4] = -0.431051;
    filter[1][0] = -0.304412;
    filter[1][1] = 0.706700;
    filter[1][2] = 0.745220;
    filter[1][3] = -0.362229;
    filter[1][4] = -0.746554;
    filter[2][0] = -0.705139;
    filter[2][1] = 0.231120;
    filter[2][2] = 1.000000;
    filter[2][3] = 0.231120;
    filter[2][4] = -0.705139;
    filter[3][0] = -0.746554;
    filter[3][1] = -0.362229;
    filter[3][2] = 0.745220;
    filter[3][3] = 0.706700;
    filter[3][4] = -0.304412;
    filter[4][0] = -0.431051;
    filter[4][1] = -0.668173;
    filter[4][2] = 0.193165;
    filter[4][3] = 0.757461;
    filter[4][4] = 0.158635;

    sum = 0;
	for (i=0;i<filtRows;i++)
    {
        for (j=0;j<filtCols;j++)
        {
            sum += (filter[i][j] * filter[i][j]);
        }
    }
	for (i=0;i<filtRows;i++)
    {
        for (j=0;j<filtCols;j++)
        {
            filter[i][j] /= sum;
        }
    }
    printf("sum: %f\n",sum);
    sum = 0;
	for (i=0;i<filtRows;i++)
    {
        for (j=0;j<filtCols;j++)
        {
            sum += (filter[i][j] * filter[i][j]);
        }
    }
    printf("sum: %f\n",sum);

    int obsHistogram[256];
    int synthHistogram[256];

    convolve(rows,cols,filtRows,filtCols,img,filteredImg,filter);    
    
    float max;
    float min;
    max = filteredImg[0][0];
    min = filteredImg[0][0];
    for (i=0;i<rows;i++)
    {
        for (j=0;j<cols;j++)
        {
            if (filteredImg[i][j] > max)
                max = filteredImg[i][j];
            if (filteredImg[i][j] < min)
                min = filteredImg[i][j];
        }
    }

    if (max-min != 0)
    {
        for (i=0;i<rows;i++)
        {
            for (j=0;j<cols;j++)
            {
                filteredImg[i][j] += -1 * min;
                filteredImg[i][j] *= 255;
                filteredImg[i][j] = (int)(filteredImg[i][j]/(max - min));
            }
        }
    }

    printf("min: %f max: %f\n",min,max);


    computeHistogram(rows,cols,filteredImg,obsHistogram);

    float lambdas[noFilts][256];

    for (i=0;i<noFilts;i++)
    {
        for (j=0;j<256;j++)
            lambdas[i][j] = 0;
    }

    for (i=0;i<rows;i++)
    {
        for (j=0;j<cols;j++)
        {
            synthImg[i][j] = (rand()/RAND_MAX) * 255;
        }
    }

    int count = 0;
    int loopNo = 0;
    do
    {
        loopNo++;
        convolve(rows,cols,filtRows,filtCols,synthImg,filteredImg,filter);
        max = filteredImg[0][0];
        min = filteredImg[0][0];
        for (i=0;i<rows;i++)
        {
            for (j=0;j<cols;j++)
            {
                if (filteredImg[i][j] > max)
                    max = filteredImg[i][j];
                if (filteredImg[i][j] < min)
                    min = filteredImg[i][j];
            }
        }

        if (max-min != 0)
        {
            for (i=0;i<rows;i++)
            {
                for (j=0;j<cols;j++)
                {
                    filteredImg[i][j] += -1 * min;
                    filteredImg[i][j] *= 255;
                    filteredImg[l][m] = (int)(filteredImg[l][m]/(max - min));
                }
            }
        }
        printf("min: %f max: %f\n",min,max);
        computeHistogram(rows,cols,filteredImg,synthHistogram);

        for (i=0;i<256;i++)
        {
            lambdas[0][i] += delta * (synthHistogram[i] - obsHistogram[i]);
        }

        curEnergy = energy(noFilts,lambdas,synthHistogram);

        prob = exp(-1 * curEnergy);

        
        count = 0;
        for (k=1;k<=100;k+=1)
        {
            count++;
            beta = 3 / log(k+1);
            //beta = (rows * cols * delta)/log(k+1);
            //beta = 1;
            curEnergy = energy(noFilts,lambdas,synthHistogram);
            for (i=0;i<rows;i++)
            {
                if ((i+1)%rows == 0)
                    printf("count: %d\n",count);
                for (j=0;j<cols;j++)
                {
                        //metropolis start
                        //initNbhdEnergy = nbhdEnergy(rows,cols,img,origImg,x,y);

                        randomNumber = rand();
                        randomNumber = (randomNumber/RAND_MAX) * 255;
                        randInt = randomNumber;

                        imgCopy1[i][j] = randInt;
                        //curNbhdEnergy = nbhdEnergy(rows,cols,imgCopy1,origImg,x,y);
                        //current = *initEnergy - initNbhdEnergy + curNbhdEnergy;
                        localConvolve(rows,cols,filtRows,filtCols,imgCopy1,filteredImg,filter,j,i);

                        max = filteredImg[0][0];
                        min = filteredImg[0][0];
                        for (l=0;l<rows;l++)
                        {
                            for (m=0;m<cols;m++)
                            {
                                if (filteredImg[l][m] > max)
                                    max = filteredImg[i][j];
                                if (filteredImg[l][m] < min)
                                    min = filteredImg[i][j];
                            }
                        }

                        if (max-min != 0)
                        {
                            for (l=0;l<rows;l++)
                            {
                                for (m=0;m<cols;m++)
                                {
                                    filteredImg[l][m] += -1 * min;
                                    filteredImg[l][m] *= 255;
                                    filteredImg[l][m] = (int)(filteredImg[l][m]/(max - min));
                                }
                            }
                        }

                        computeHistogram(rows,cols,filteredImg,synthHistogram);
                        current = energy(noFilts,lambdas,synthHistogram);


                        if (current < curEnergy)
                        {
                            curEnergy = current;
                            synthImg[i][j] = randInt;
                        }
                        else
                        {
                            randomNumber = rand();
                            randomNumber = randomNumber/RAND_MAX;
                            prob = exp(-(1/beta) * (current - curEnergy));
                            //printf("%f %f %f %f %f %f\n",beta,prob,curNbhdEnergy,initNbhdEnergy,*initEnergy,current);
                            if (randomNumber < prob)
                            {
                                curEnergy = current;
                                synthImg[i][j] = randInt;
                            }
                            else
                            {
                                synthImg[i][j] = synthImg[i][j];
                            }
                        }
                        //metropolis end
                        imgCopy1[i][j] = synthImg[i][j];
                }
            }

            //mod = 100;
            //if (k%mod == 0)
            //{
            //    strcpy(oName,basename);
            //    sprintf(suffix,"-smooth-%d.csv",k/mod);
            //    strcat(oName,suffix);
            //    file = fopen(oName,"w");

            //    for (i=0;i<rows;i++)
            //    {
            //        for (j=0;j<cols;j++)
            //            fprintf(file,"%d,",img[i][j]);
            //        fprintf(file,"\n");
            //    }
            //    fclose(file);
            //}
        }
        printf("Loop %d done\n",loopNo);
        diff = 0;

        computeHistogram(rows,cols,filteredImg,synthHistogram);
        
        for (i=0;i<256;i++)
        {
            diff += abs(obsHistogram[i] - synthHistogram[i]);
        }
        diff = diff/2;
        printf("diff: %d\n",diff);
    } while (diff > epsilon);





    for (i=0;i<rows;i++)
    {
        for (j=0;j<cols;j++)
        {
            img[i][j] = filteredImg[i][j];
            //img[i][j] = synthImg[i][j];
        }
    }

////////////////////////////////////////////////////////////////////////////////
    for (i=0;i<filtRows;i++)
    {
        free(filter[i]);
    }
    free(filter);

    for (i=0;i<rows;i++)
    {
        free(filteredImg[i]);
    }
    free(filteredImg);

    for (i=0;i<rows;i++)
    {
        free(synthImg[i]);
    }
    free(synthImg);

    for (i=0;i<rows;i++)
    {
        free(imgCopy1[i]);
    }
    free(imgCopy1);
/////////////////////////////////////////////////////////////////////////////////

    strcpy(oName,basename);
    strcat(oName,"-filtered.csv");
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
    return 1;
}

void localConvolve(int imgRows, int imgCols, int filtRows, int filtCols, int** imgIn, int** imgOut, float** filter,int x, int y)
{
    int i,j,k,l;
    int pixTemp;
    int iCoord, jCoord;

    for (i=-filtRows;i<filtRows;i++)
    {
        for (j=-filtCols;j<filtCols;j++)
        {
            pixTemp = 0;

            iCoord = i+y;
            jCoord = j+x;
            
            if (iCoord < imgRows && iCoord >= 0 && jCoord < imgCols && jCoord >= 0)
            {
                for (k=0;k<filtRows;k++)
                {
                    for (l=0;l<filtCols;l++)
                    {
                        if (iCoord+k-(filtRows/2) < imgRows && iCoord+k-(filtRows/2) >= 0 && jCoord+l-(filtCols/2) < imgCols && jCoord+l-(filtCols/2) >= 0)
                            pixTemp += imgIn[iCoord+k-(filtRows/2)][jCoord+l-(filtCols/2)] * filter[k][l];
                    }
                }
                imgOut[iCoord][jCoord] = pixTemp;
            }
        }
    }
}

void convolve(int imgRows, int imgCols, int filtRows, int filtCols, int** imgIn, int** imgOut, float** filter)
{
    int i,j,k,l;
    int pixTemp;

    for (i=0;i<imgRows;i++)
    {
        for (j=0;j<imgCols;j++)
        {
            pixTemp = 0;
            for (k=0;k<filtRows;k++)
            {
                for (l=0;l<filtCols;l++)
                {
                    if (i+k-(filtRows/2) < imgRows && i+k-(filtRows/2) >= 0 && j+l-(filtCols/2) < imgCols && j+l-(filtCols/2) >= 0)
                        pixTemp += imgIn[i+k-(filtRows/2)][j+l-(filtCols/2)] * filter[k][l];
                }
            }
//            if (pixTemp > 255 || pixTemp < 0)
//                printf("i: %d j: %d: pixTemp: %d\n",i,j,pixTemp);
            imgOut[i][j] = pixTemp;
        }
    }
}

void computeHistogram(int imgRows, int imgCols, int** img, int* histogram)
{
    int i,j,k;

    for (i=0;i<256;i++)
    {
        histogram[i] = 0;
    }

    for (i=0;i<imgRows;i++)
    {
        for (j=0;j<imgCols;j++)
        {
            if (img[i][j] > 255 || img[i][j] < 0)
                printf("i: %d j: %d: img[i][j]: %d\n",i,j,img[i][j]);
            histogram[img[i][j]]++;
        }
    }
}

double energy(int noFilts, float lambdas[][256], int* histogram)
{
    int i,j;
    float ans=0;

    for (i=0;i<noFilts;i++)
        for (j=0;j<256;j++)
    {
        {
            ans += lambdas[i][j] * histogram[j];
        }
    }
    return ans;
}

