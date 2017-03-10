#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

void localConvolve(int imgRows, int imgCols, int filtRows, int filtCols, int** imgIn, int** imgOut, double** filter,int x, int y);
void convolve(int imgRows, int imgCols, int filtRows, int filtCols, int** imgIn, int** imgOut, double** filter);
void computeHistogram(int imgRows, int imgCols, int** img, int* histogram);
double energy(int noFilts, float* lambdas, int* histogram);
void gaussian(double** filter, double sigma, int filtRows, int filtCols);
void gabor(double** filter, double sigma, int theta, int filtRows, int filtCols);

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

    long epsilon = 700;
    double delta = 0.5;
    long diff;
    double sum;
    double prob;
    float curEnergy;
    float beta;
    int noFilts = 1;

    double probs[256];
    double energies[256];
    double randomNumber;
    double initNbhdEnergy,curNbhdEnergy,current;
    double max;
    int randInt;
    int backup;

    int obsHistogram[256];
    int synthHistogram[256];
    

    // filter stuff
    int filtSize = 5;
    int filtRows = filtSize;
    int filtCols = filtSize;
    double** filter3;
    filter3 = malloc(filtRows * sizeof(*filter3));
    for (i=0; i<filtCols; i++)
        filter3[i] = malloc(filtRows * sizeof(*filter3[i]));
    gabor(filter3, 3, 90, filtRows, filtCols);


    // main code
    convolve(rows,cols,filtRows,filtCols,img,filteredImg,filter3);  
    computeHistogram(rows,cols,filteredImg,obsHistogram);

    float lambdas[256];

    for (i=0;i<256;i++)
            lambdas[i] = 0;

    for (i=0;i<rows;i++){
        for (j=0;j<cols;j++)
            synthImg[i][j] = rand()%256;
    }
    int count = 0;
    int loopNo = 0;
    do
    {
        loopNo++;
        convolve(rows,cols,filtRows,filtCols,synthImg,filteredImg,filter3);
        computeHistogram(rows,cols,filteredImg,synthHistogram);

        for (i=0;i<256;i++)
            lambdas[i] += delta * (synthHistogram[i] - obsHistogram[i]);

        curEnergy = energy(noFilts,lambdas,synthHistogram);
        prob = exp(-1 * curEnergy);

        
        count = 0;
        for (k=1;k<=100;k+=1){
            count++;
            beta = 3 / log(k+1);
            curEnergy = energy(noFilts,lambdas,synthHistogram);
            for (i=0;i<rows;i++){
                if ((i+1)%rows == 0)
                    printf("count: %d\n",count);
                for (j=0;j<cols;j++){
                        //metropolis start
                        randomNumber = rand()% 256;
                        randInt = randomNumber;
                        imgCopy1[i][j] = randInt;
                        localConvolve(rows,cols,filtRows,filtCols,imgCopy1,filteredImg,filter3,j,i);
                        //convolve(rows,cols,filtRows,filtCols,imgCopy1,filteredImg,filter3);
                        computeHistogram(rows,cols,filteredImg,synthHistogram);
                        current = energy(noFilts,lambdas,synthHistogram);
                        if (current < curEnergy){
                            curEnergy = current;
                            synthImg[i][j] = randInt;
                        }
                        else{
                            randomNumber = rand()%256;
                            prob = exp(-(1/beta) * (current - curEnergy));
                            if (randomNumber < prob){
                                curEnergy = current;
                                synthImg[i][j] = randInt;
                            }
                            else{
                                synthImg[i][j] = synthImg[i][j];
                            }
                        }
                        //metropolis end
                        imgCopy1[i][j] = synthImg[i][j];
                }
            }
        }
        printf("Loop %d done\n",loopNo);
        diff = 0;

        computeHistogram(rows,cols,filteredImg,synthHistogram);
        
        for (i=0;i<256;i++){
            diff += abs(obsHistogram[i] - synthHistogram[i]);
        }
        diff = diff/2;
        printf("diff: %d\n",diff);
    } while (loopNo < 12);
    //} while (diff > epsilon);

    for (i=0;i<rows;i++){
        for (j=0;j<cols;j++){
            img[i][j] = synthImg[i][j];
        }
    }
    
////////////////////////////////////////////////////////////////////////////////

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
    return 0;
}

double energy(int noFilts, float *lambdas, int* histogram){
    int i;
    float ans=0;

        for (i=0;i<256;i++){
            ans += lambdas[i] * histogram[i];
        }
    return ans;
}

void computeHistogram(int imgRows, int imgCols, int** img, int* histogram){
    int i,j,k;

    for (i=0;i<256;i++)
        histogram[i] = 0;

    for (i=0;i<imgRows;i++){
        for (j=0;j<imgCols;j++){
            histogram[img[i][j]]++;
            if(img[i][j] > 255 || img[i][j]<0)
                printf("%d\n", img[i][j] );
        }
    }
}

void localConvolve(int imgRows, int imgCols, int filtRows, int filtCols, int** imgIn, int** imgOut, double** filter,int a, int b){

    int x,y,i,j;
    double sum;
    int m = (filtRows-1)/2;
    int n = (filtCols-1)/2;

    for(x=a-m; x<a+m; x++){
        for(y=b-n; y<b+n; y++){
            sum = 0;
            if(x<m || x>=imgRows-m || y<n || y>=imgCols-n)
                continue;
            for(i=-m; i<=m; i++){
                for(j=-n; j<=n; j++){
                    sum = sum + filter[j+n][i+m]*imgIn[x-j][y-i];  
                }
            }
        imgOut[x][y] = sum;
        }
    } 

    double min = imgOut[0][0];
    double max = imgOut[0][0];
    for(x=m; x<imgRows-m; x++){
        for(y=n; y<imgCols-m; y++){
            if(imgOut[x][y]>max)
                max = imgOut[x][y];
            if(imgOut[x][y]<min)
                min = imgOut[x][y];
        }
    }  
    for(x=m; x<imgRows-m; x++){
        for(y=n; y<imgCols-m; y++){
            imgOut[x][y] = (imgOut[x][y]-min)*255/(max-min);
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

    double min = imgOut[0][0];
    double max = imgOut[0][0];
    for(x=m; x<imgRows-m; x++){
        for(y=n; y<imgCols-m; y++){
            if(imgOut[x][y]>max)
                max = imgOut[x][y];
            if(imgOut[x][y]<min)
                min = imgOut[x][y];
        }
    }  
    for(x=m; x<imgRows-m; x++){
        for(y=n; y<imgCols-m; y++){
            imgOut[x][y] = (imgOut[x][y]-min)*255/(max-min);
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
            printf("%f\n", filter[i][j]);
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

    for(x=0 ; x<filtRows; x++){
        for(y=0; y<filtCols; y++){

            xp = (x-m)*cos(t) + (y-n)*sin(t);
            yp = -(x-m)*sin(t) + (y-n)*cos(t);
            filter[x][y] = exp(-(1/(4*s*s))*(4*(xp*xp+yp*yp)))*cos(2*M_PI*xp*(1/(sqrt(2)*s)));
        }
    }
}
