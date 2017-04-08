#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

void convolve(int imgRows, int imgCols, int filtRows, int filtCols, int** imgIn, int** imgOut, double** filter);
void localHistogramOut(int* histogram, int** img, int imgRows, int imgCols, int a, int b);
void localHistogramIn(int* histogram, int** img, int imgRows, int imgCols, int a, int b);
void localConvolve(int imgRows, int imgCols, int filtRows, int filtCols, int** imgIn, int** imgOut, double** filter,int x, int y);
void computeHistogram(int imgRows, int imgCols, int** img, int* histogram);
double energy(int noFilts, double* lambdas, double* histogram);
void gaussian(double** filter, double sigma, int filtRows, int filtCols);
void gabor(double** filter, double sigma, int theta, int filtRows, int filtCols);
double absf(double x);
void norm(int imgRows, int imgCols, int* histogram, double* histOut);

int main(int argc, char* argv[])
{
    if (argc != 4)
    {
        printf("Invalid number of arguments\n%s [BASENAME] [COLS] [ROWS]\n",argv[0]);
        return 1;
    }
    int i=0, j=0, k=0, l=0, m=0;
    char fName[50];
    char oName[50];
    char suffix[50];
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
    for (i=0;i<rows;i++){
        for (j=0;j<cols;j++){
            img[i][j] = (int)(img[i][j]/32);
        }
    }

//////////////////////////////////////////////////////////////////////////////////

    // Initialize prng
    srand(time(NULL));

    int** synthImg;
    synthImg = malloc(rows * sizeof(*synthImg));
    for (i=0;i<rows;i++)
    {
        synthImg[i] = malloc(cols * sizeof(*synthImg[i]));
    }

    int** imgCopy1;
    imgCopy1 = malloc(rows * sizeof(*imgCopy1));
    for (i=0;i<rows;i++)
    {
        imgCopy1[i] = malloc(cols * sizeof(*imgCopy1[i]));
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

    double min = 10000000;
    double diff;
    double curEnergy;
    int noFilts = 7;
    double randomNumber;
    double current;
    int randInt;
    int count = 0;
    int loopNo = 0;

    int** obsHistogram;
    obsHistogram = malloc(noFilts * sizeof(*obsHistogram));
    for (i=0; i<noFilts; i++)
        obsHistogram[i] = malloc(8 * sizeof(*obsHistogram[i]));

    int** synthHistogram;
    synthHistogram = malloc(noFilts * sizeof(*synthHistogram));
    for (i=0; i<noFilts; i++)
        synthHistogram[i] = malloc(8 * sizeof(*synthHistogram[i]));

    double** normObsHistogram;
    normObsHistogram = malloc(noFilts * sizeof(*normObsHistogram));
    for (i=0; i<noFilts; i++)
        normObsHistogram[i] = malloc(8 * sizeof(*normObsHistogram[i]));

    double** normSynthHistogram;
    normSynthHistogram = malloc(noFilts * sizeof(*normSynthHistogram));
    for (i=0; i<noFilts; i++)
        normSynthHistogram[i] = malloc(8 * sizeof(*normSynthHistogram[i]));

    // filter stuff
    int filtSize = 5;
    int filtRows[noFilts];
    int filtCols[noFilts];

    filtRows[0] = 1;
    filtRows[1] = 5;
    filtRows[2] = 5;
    filtRows[3] = 5;
    filtRows[4] = 5;
    filtRows[5] = 5;
    filtRows[6] = 5;

    for (i=0;i<noFilts;i++){
        filtCols[i] = filtRows[i];
    }

    double*** filter;
    filter = malloc(noFilts * sizeof(double**));
    for (i=0;i<noFilts;i++){
        filter[i] = malloc(filtRows[i] * sizeof(double*));
        for (j=0;j<filtRows[i];j++){
            filter[i][j] = malloc(filtCols[i] * sizeof(double));
        }
    }

    int*** filteredImg;
    filteredImg = malloc(noFilts * sizeof(int**));
    for (i=0;i<noFilts;i++){
        filteredImg[i] = malloc(rows * sizeof(int*));
        for (j=0;j<rows;j++){
            filteredImg[i][j] = malloc(cols * sizeof(int));
        }
    }

    filter[0][0][0] = 1;
    gabor(filter[1],2,45,filtRows[1],filtCols[1]);
    gabor(filter[2],2,135,filtRows[2],filtCols[2]);
    gaussian(filter[3],2,filtRows[3],filtCols[3]);
    gabor(filter[4],2,0,filtRows[4],filtCols[4]);
    gabor(filter[5],2,90,filtRows[5],filtCols[5]);
    gaussian(filter[6],1,filtRows[6],filtCols[6]);
    

    // lambda stuff
    double** lambdas;
    lambdas = malloc(noFilts *sizeof(*lambdas));
    for(i=0;i<noFilts;i++)
        lambdas[i] = malloc(8*sizeof(*lambdas[i]));

    // main code
    computeHistogram(rows,cols,img,obsHistogram[0]);
    norm(rows,cols,obsHistogram[0], normObsHistogram[0]);
    for (i=1;i<noFilts;i++){
        convolve(rows,cols,filtRows[i],filtCols[i],img,filteredImg[i],filter[i]);
        computeHistogram(rows,cols,filteredImg[i],obsHistogram[i]);
        norm(rows,cols,obsHistogram[i], normObsHistogram[i]);
    }

    for(i=0;i<noFilts;i++){
        for (j=0;j<8;j++)
            lambdas[i][j] = 0;
    }
    for (i=0;i<rows;i++){
        for (j=0;j<cols;j++){
            synthImg[i][j] = rand()%8;
        }
    }
    for (i=0;i<rows;i++){
        for (j=0;j<cols;j++)
            imgCopy1[i][j] = synthImg[i][j];
    }

    double temp = 32;
    double dt = 0.05;
    int x,y;
    double probs[8];
    double nrg[8];
    double sum;
    double nrgN[8];
    double minNrg;
    double maxNrg;

    computeHistogram(rows,cols,synthImg,synthHistogram[0]);
    norm(rows,cols,synthHistogram[0], normSynthHistogram[0]);
    for (m=1;m<noFilts;m++){
        convolve(rows,cols,filtRows[m],filtCols[m],synthImg,filteredImg[m],filter[m]);
        computeHistogram(rows,cols,filteredImg[m],synthHistogram[m]);
        norm(rows,cols,synthHistogram[m], normSynthHistogram[m]);
    }

    do{
        loopNo++;
        for (m=0;m<noFilts;m++){
            for (i=0;i<8;i++){
                lambdas[m][i] = lambdas[m][i] + temp * (normSynthHistogram[m][i] - normObsHistogram[m][i]);
            }
        }

        if(temp>1)
            temp = temp * (1-dt);
        else
            temp = 1;

        for (k=1;k<=4*(rows-4)*(cols-4);k++){

            sum = 0;
            minNrg = 999999999;
            maxNrg = -999999999;
            x = (rand()%(rows-4))+2;
            y = (rand()%(cols-4))+2;


            for(i=0;i<8;i++){

                nrg[i] = 0;

                localHistogramOut(synthHistogram[0], imgCopy1, rows, cols, x, y);
                imgCopy1[x][y] = i;
                localHistogramIn(synthHistogram[0], imgCopy1, rows, cols, x, y);
                norm(rows,cols,synthHistogram[0], normSynthHistogram[0]);

                for (m=1;m<noFilts;m++){
                    localHistogramOut(synthHistogram[m], filteredImg[m], rows, cols, x, y);
                    localConvolve(rows,cols,filtRows[m],filtCols[m],imgCopy1,filteredImg[m],filter[m],x,y);
                    localHistogramIn(synthHistogram[m], filteredImg[m], rows, cols, x, y);
                    norm(rows,cols,synthHistogram[m], normSynthHistogram[m]);
                }
                for (m=0;m<noFilts;m++){
                    nrg[i] += energy(noFilts,lambdas[m],normSynthHistogram[m]);
                }

                if(nrg[i] < minNrg)
                    minNrg = nrg[i];
                if(nrg[i] > maxNrg)
                    maxNrg = nrg[i];
            }

            for(i=0;i<8;i++){
                nrgN[i] = (nrg[i] - minNrg)/(maxNrg - minNrg);
                probs[i] = exp(-nrgN[i]*50);
                sum += probs[i];
            }

            for(i=0;i<8;i++)
                probs[i] = probs[i]/sum;

            randomNumber = rand();
            randomNumber = randomNumber/RAND_MAX;

            //if(x==100 && y == 100){
            //    printf("p(0) is %f  ", probs[0]);
            //    printf("p(7) is %f  ", probs[7]);
            //    printf("rand is %f\n", randomNumber);
            //}

            for (i=1;i<8;i++)
                probs[i] = probs[i-1] + probs[i];

            i = 0;
            while (1){
                if (randomNumber < probs[i])
                    break;
                if (i == 7)
                    break;
                i++;
            }

            //if(x==100 && y == 100){
            //    printf("Changed pixel from %d  ",synthImg[x][y]);
            //    printf("to %d\n",i);
            //}

            synthImg[x][y] = i;
            imgCopy1[x][y] = 7;

            localHistogramOut(synthHistogram[0], imgCopy1, rows, cols, x, y);
            imgCopy1[x][y] = i;    
            localHistogramIn(synthHistogram[0], imgCopy1, rows, cols, x, y);
            norm(rows,cols,synthHistogram[0], normSynthHistogram[0]);

            for (m=1;m<noFilts;m++){
                localHistogramOut(synthHistogram[m], filteredImg[m], rows, cols, x, y);
                localConvolve(rows,cols,filtRows[m],filtCols[m],imgCopy1,filteredImg[m],filter[m],x,y);
                localHistogramIn(synthHistogram[m], filteredImg[m], rows, cols, x, y);
                norm(rows,cols,synthHistogram[m], normSynthHistogram[m]);
            }
        }
        
        diff = 0;
        for (m=0;m<noFilts;m++){
            for (i=0;i<8;i++){
                diff += absf(normObsHistogram[m][i] - normSynthHistogram[m][i]);
            }
        }

        if(loopNo < 50 && loopNo > 0){
            if(diff<min){
                min = diff;
                printf("\nBANG NEW MIN\n");
            }
        }

        printf("\n diff = %f and min = %f and temp = %f\n",diff, min, temp);
        printf("------Loop %d -----\n",loopNo);


        if(loopNo%10 == 0){
            for (i=0;i<rows;i++){
                for (j=0;j<cols;j++){
                    img[i][j] = (int)(synthImg[i][j]*32);
                }
            }
            strcpy(oName,basename);
            sprintf(suffix,"-output-%d.csv",loopNo);
            strcat(oName,suffix);
            file = fopen(oName,"w");

            for (i=0;i<rows;i++){
                for (j=0;j<cols;j++)
                    fprintf(file,"%d,",img[i][j]);
                fprintf(file,"\n");
            }
            fclose(file);
        }

    } while ((diff > min || loopNo<=50) && loopNo<50);

    for (i=0;i<rows;i++){
        for (j=0;j<cols;j++){
            img[i][j] = synthImg[i][j];
        }
    }

/////////////////////////////////////////////////////////////////////////////////

    for (i=0;i<rows;i++)
    {
        for (j=0;j<cols;j++)
        {
            img[i][j] = (int)(img[i][j]*32);
        }
    }
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

    strcpy(oName,basename);
    strcat(oName,"-lambdas.csv");
    file = fopen(oName,"w");

    for (i=0;i<noFilts;i++)
    {
        for (j=0;j<8;j++)
            fprintf(file,"%f,",lambdas[i][j]);
        fprintf(file,"\n");
    }

    fclose(file);

}

void norm(int imgRows, int imgCols, int* histogram, double* histOut){

    int i;
    for(i=0;i<8;i++){
        histOut[i] = histogram[i] * (1.0/((imgRows-4)*(imgCols-4)));
    }
}

void localHistogramOut(int* histogram, int** img, int imgRows, int imgCols, int a, int b){

    int x,y;
    int m = 2;
    int n = 2;

    for(x=a-m; x<=a+m; x++){
        for(y=b-n; y<=b+n; y++){
            if(img[x][y] > 7)
                histogram[7] = histogram[7] - 1;
            else if(img[x][y] < 0)
                histogram[0] = histogram[0] - 1;
            else
                histogram[img[x][y]] = histogram[img[x][y]] - 1;
        }
    } 
}

void localHistogramIn(int* histogram, int** img, int imgRows, int imgCols, int a, int b){

    int x,y;
    int m = 2;
    int n = 2;

    for(x=a-m; x<=a+m; x++){
        for(y=b-n; y<=b+n; y++){
            if(img[x][y] > 7)
                histogram[7] = histogram[7] + 1;
            else if(img[x][y] < 0)
                histogram[0] = histogram[0] + 1;
            else
                histogram[img[x][y]] = histogram[img[x][y]] + 1;
        }
    } 
}

void localConvolve(int imgRows, int imgCols, int filtRows, int filtCols, int** imgIn, int** imgOut, double** filter,int a, int b){

    int x,y,i,j;
    double sum;
    int m = (filtRows-1)/2;
    int n = (filtCols-1)/2;

    for(x=a-m; x<=a+m; x++){
        for(y=b-n; y<=b+n; y++){
            if(x<m || y<n || x>=imgRows-m || y>=imgCols-n)
                continue;
            sum = 0;
            for(i=-m; i<=m; i++){
                for(j=-n; j<=n; j++){
                    sum = sum + filter[j+n][i+m]*imgIn[x-j][y-i]; 
                }
            }
        imgOut[x][y] = sum;
        }
    } 
}

double energy(int noFilts, double *lambdas, double* histogram){
    int i;
    double ans=0;
        for (i=0;i<8;i++){
            ans += lambdas[i]*histogram[i];
        }
    return ans;
}

void computeHistogram(int imgRows, int imgCols, int** img, int* histogram){
    int i,j,k;

    for (i=0;i<8;i++)
        histogram[i] = 0;

    for (i=2;i<imgRows-2;i++){
        for (j=2;j<imgCols-2;j++){
            if (img[i][j] < 0)
                histogram[0]++;
            else if (img[i][j] > 7)
                histogram[7]++;
            else
                histogram[img[i][j]]++;
        }
    }
}

void convolve(int imgRows, int imgCols, int filtRows, int filtCols, int** imgIn, int** imgOut, double** filter){

    int x,y,i,j;
    double sum;
    int m = (filtRows-1)/2;
    int n = (filtCols-1)/2;

    for(x=m; x<imgRows-m; x++){
        for(y=n; y<imgCols-m; y++){
            sum = 0;
            for(i=-m; i<=m; i++){
                for(j=-n; j<=n; j++){
                    sum = sum + filter[j+n][i+m]*imgIn[x-j][y-i];  
                }
            }
        imgOut[x][y] = sum;
        }
    } 
}

void gaussian(double** filter, double sigma, int filtRows, int filtCols){

    int i;

    double r;
    double t = 2.0 * sigma * sigma;
    int a = filtRows/2;
    int b = filtCols/2;
    // sum is for normalization
    double sum = 0.0;
    // generate kernel
    for (int x = -a; x <= a; x++){
        for(int y = -b; y <= b; y++){
            r = sqrt(x*x + y*y);
            filter[x + a][y + b] = (exp(-(r*r)/t))/(M_PI * t);
            sum += filter[x + a][y + b];
        }
    }
    // normalize the Kernel
    for(int i = 0; i < filtRows; ++i){
        for(int j = 0; j < filtCols; ++j){
            filter[i][j] /= sum;
        }
    }
}

void gabor(double** filter, double sigma, int theta, int filtRows, int filtCols){

    int x,y;
    double xp,yp,t,s;
    double m = (filtRows)/2;
    double n = (filtCols)/2;

    // radian conversion
    t = theta * M_PI/180;
    s = sigma;

    for(x=0; x<filtRows; x++){
        for(y=0; y<filtCols; y++){
            xp = (x-m)*cos(t) + (y-n)*sin(t);
            yp = -(x-m)*sin(t) + (y-n)*cos(t);
            filter[x][y] = exp(-(1/(4*s*s))*(4*(xp*xp+yp*yp)))*cos(2*M_PI*xp*(1/(sqrt(2)*s)));
        }
    }
}
double absf(double x)
{
    return x >= 0 ? x : -1 * x;
}
