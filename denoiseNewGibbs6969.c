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
double localEnergyMatch(int x, int y, int** img, int** origImg);
double energy(int noFilts, double* lambdas, double* histogram);
double energyMatch(int rows, int cols, int** img, int** origImg);
void gaussian(double** filter, double sigma, int filtRows, int filtCols);
void gabor(double** filter, double sigma, int theta, int filtRows, int filtCols);
double absf(double x);
void norm(int imgRows, int imgCols, int* histogram, double* histOut);

float sigma = 1;
float lambda = 300;

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
    for (i=0;i<rows;i++)
    {
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

    //Theorem-big
    float tempLambdas[7][8] = {{-5.172199,1.870683,1.511740,2.176885,1.897370,2.939844,-0.305859,-4.918465},
                                {-0.537380,-4.751398,-1.556484,0.154660,1.053504,1.102417,1.463649,3.071032},
                                {-1.113643,-5.065499,-2.549878,-0.068455,1.355884,1.476000,1.927394,4.038198},
                                {-1.706341,5.232651,3.227347,-0.127444,-4.631689,-6.412724,4.418200,0.000000},
                                {-5.768399,-0.798955,0.625114,1.169057,1.303987,0.906752,0.222656,2.339788},
                                {-1.445146,0.556661,0.653496,0.251219,0.986082,-0.126777,-1.916005,1.040471},
                                {-2.591420,-2.237303,2.340083,3.134268,6.502385,4.316854,-3.638152,-7.826716}};


    //Theorem
    //float tempLambdas[7][8] = {{-4.772002,1.259613,0.740862,0.884719,0.918161,1.790169,1.257460,-2.078982},
    //                            {-3.216437,-2.364406,-0.787150,0.073240,0.786722,1.176355,1.638008,2.693668},
    //                            {-2.693116,-1.964676,-0.396044,0.156014,1.402432,1.354778,0.935025,1.205587},
    //                            {0.719284,1.319837,1.229935,2.697305,-2.981215,-1.771714,-1.213431,0.000000},
    //                            {-3.676325,-2.070627,-0.702331,0.430259,1.052437,1.277627,1.489661,2.199299},
    //                            {-0.566914,0.752498,1.244862,0.697884,0.833054,-1.146222,-0.899837,-0.915325},
    //                            {0.485283,-2.032143,-6.059077,1.433645,6.504421,7.798379,2.157284,-10.287793}};

    //Snell
    //float tempLambdas[7][8] = {{-10.104918,2.756894,1.398456,2.319907,3.360725,2.482131,2.823906,-5.037102},
    //                            {-4.157557,-4.281057,-1.191598,1.490794,1.696462,1.508549,1.809487,3.124919},
    //                            {-4.525711,-4.186759,-0.677246,1.299550,1.835748,1.679962,1.346148,3.228308},
    //                            {2.328334,-0.153508,-7.217595,0.386008,1.407728,2.382952,0.866082,0.000000},
    //                            {-7.639743,1.658043,2.022225,1.898639,1.490888,0.134796,-0.154279,0.589432},
    //                            {-1.324878,1.032792,0.135508,0.597000,0.942994,0.202214,-0.739702,-0.845928},
    //                            {-4.093105,-2.928330,-4.230888,5.617408,9.579729,4.535682,0.280546,-8.761042}};

                                

    for (i=0;i<noFilts;i++)
    {
        for (j=0;j<8;j++)
        {
            lambdas[i][j] = tempLambdas[i][j];
        }
    }

    // main code
    computeHistogram(rows,cols,img,obsHistogram[0]);
    norm(rows,cols,obsHistogram[0], normObsHistogram[0]);
    for (i=1;i<noFilts;i++){
        convolve(rows,cols,filtRows[i],filtCols[i],img,filteredImg[i],filter[i]);
        computeHistogram(rows,cols,filteredImg[i],obsHistogram[i]);
        norm(rows,cols,obsHistogram[i], normObsHistogram[i]);
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

    float beta;
    double initEnergy = energyMatch(rows, cols, imgCopy1, img);
    double pt1 = 0;
    double tempEnergy;
    for (k=1;k<=200*(rows-4)*(cols-4);k++){
        if (!(k%((rows-4)*(cols-4))))
            printf("---Loop: %d---\n",k/((rows-4)*(cols-4)));
        beta = 3 / log(k+1);

        sum = 0;
        minNrg = 999999999;
        maxNrg = -999999999;
        x = (rand()%(rows-4))+2;
        y = (rand()%(cols-4))+2;


        pt1 = initEnergy - localEnergyMatch(x,y,imgCopy1,img);
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
            tempEnergy = pt1;
            tempEnergy += localEnergyMatch(x,y,imgCopy1,img);

            //printf("tempEnergy: %f\tnrg[i]: %f\n",tempEnergy,nrg[i]);
            nrg[i] =  tempEnergy + lambda * nrg[i];

            if(nrg[i] < minNrg)
                minNrg = nrg[i];
            if(nrg[i] > maxNrg)
                maxNrg = nrg[i];
        }

        for(i=0;i<8;i++){
            nrgN[i] = (nrg[i] - minNrg)/(maxNrg - minNrg);
            probs[i] = exp(-nrgN[i]*50*beta);
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

        initEnergy = pt1 + localEnergyMatch(x,y,imgCopy1,img);

        for (m=1;m<noFilts;m++){
            localHistogramOut(synthHistogram[m], filteredImg[m], rows, cols, x, y);
            localConvolve(rows,cols,filtRows[m],filtCols[m],imgCopy1,filteredImg[m],filter[m],x,y);
            localHistogramIn(synthHistogram[m], filteredImg[m], rows, cols, x, y);
            norm(rows,cols,synthHistogram[m], normSynthHistogram[m]);
        }
    }
    
    //if(loopNo%10 == 0){
    //    for (i=0;i<rows;i++){
    //        for (j=0;j<cols;j++){
    //            img[i][j] = (int)(synthImg[i][j]*32);
    //        }
    //    }
    //    strcpy(oName,basename);
    //    sprintf(suffix,"-output-%d.csv",loopNo);
    //    strcat(oName,suffix);
    //    file = fopen(oName,"w");

    //    for (i=0;i<rows;i++){
    //        for (j=0;j<cols;j++)
    //            fprintf(file,"%d,",img[i][j]);
    //        fprintf(file,"\n");
    //    }
    //    fclose(file);
    //}


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
    strcat(oName,"-restored.csv");
    file = fopen(oName,"w");

    for (i=0;i<rows;i++)
    {
        for (j=0;j<cols;j++)
            fprintf(file,"%d,",img[i][j]);
        fprintf(file,"\n");
    }

    fclose(file);
}

double localEnergyMatch(int x, int y, int** img, int** origImg)
{
    double diff, pt1;
    diff = origImg[x][y] - img[x][y];
    pt1 = abs(diff) * abs(diff)/(8.0*8.0);
    pt1 = pt1 * 1/(2.0 * sigma * sigma);
    return pt1;
}
double energyMatch(int rows, int cols, int** img, int** origImg)
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
            pt1 += abs(diff) * abs(diff) / (8.0*8.0);
        }
    }
    
    pt1 = pt1 * 1/(2.0 * sigma * sigma);
    
    ans = pt1;

    //printf("pt1: %f pt2: %f\n",pt1,pt2);
    return ans;
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
