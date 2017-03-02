#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

void convolve(int imgRows, int imgCols, int filtRows, int filtCols, int** imgIn, int** imgOut, float** filter);
void computeHistogram(int imgRows, int imgCols, int** img, int* histogram);

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
////////////////////////////////////////////////////////////////////////////////

    int filtRows = 5;
    int filtCols = 5;
    int noFilts = 1;
    int epsilon = 10;
    int diff;
    float** filter;
    filter = malloc(filtRows * sizeof(*filter));
    for (i=0;i<filtRows;i++)
    {
        filter[i] = malloc(filtCols * sizeof(*filter[i]));
    }

    filter[0][0] = (1/273.0);
    filter[0][1] = (4/273.0);
    filter[0][2] = (7/273.0);
    filter[0][3] = (4/273.0);
    filter[0][4] = (1/273.0);
    filter[1][0] = (4/273.0);
    filter[1][1] = (16/273.0);
    filter[1][2] = (26/273.0);
    filter[1][3] = (16/273.0);
    filter[1][4] = (4/273.0);
    filter[2][0] = (7/273.0);
    filter[2][1] = (26/273.0);
    filter[2][2] = (41/273.0);
    filter[2][3] = (26/273.0);
    filter[2][4] = (7/273.0);
    filter[3][0] = (4/273.0);
    filter[3][1] = (16/273.0);
    filter[3][2] = (26/273.0);
    filter[3][3] = (16/273.0);
    filter[3][4] = (4/273.0);
    filter[4][0] = (1/273.0);
    filter[4][1] = (4/273.0);
    filter[4][2] = (7/273.0);
    filter[4][3] = (4/273.0);
    filter[4][4] = (1/273.0);

    int obsHistogram[256];
    int synthHistogram[256];

    convolve(rows,cols,filtRows,filtCols,img,filteredImg,filter);    

    computeHistogram(rows,cols,filteredImg,obsHistogram);

    int lambdas[noLambdas][256];

    for (i=0;i<noLambdas;i++)
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

    do
    {
        convolve(rows,cols,filtRows,filtCols,synthImg,filteredImg,filter);
        computeHistogram(rows,cols,filteredImg,synthHistogram);

        for (i=0;i<256;i++)
        {
            lambdas[1][i] += (synthHistogram - obsHistogram);
        }

        for (i=0;i<256;i++)
        {
            
        }

    } while (diff > epsilon);





    for (i=0;i<rows;i++)
    {
        for (j=0;j<cols;j++)
        {
            img[i][j] = filteredImg[i][j];
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
    return 0;
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
                    if (i+k < imgRows && j+l < imgCols)
                        pixTemp += imgIn[i+k][j+l] * filter[k][l];
                }
            }
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
            histogram[img[i][j]]++;
        }
    }
}
