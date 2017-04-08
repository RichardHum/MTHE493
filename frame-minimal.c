#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

void localConvolve(int imgRows, int imgCols, int filtRows, int filtCols, int** imgIn, int** imgOut, double** filter,int x, int y);
void convolve(int imgRows, int imgCols, int filtRows, int filtCols, int** imgIn, int** imgOut, double** filter);
void localHistogram(int* histogram, int** img, int** img2, int filtSize, int imgRows, int imgCols, int a, int b);
void computeHistogram(int imgRows, int imgCols, int** img, int* histogram);
void normalizeHistogram(int imgRows, int imgCols, int* histIn, float* histOut);
double energy(int noFilts, float* lambdas, float* histogram);
void gaussian(double** filter, double sigma, int filtRows, int filtCols);
void gabor(double** filter, double sigma, int theta, int filtRows, int filtCols);
double absf(double x);

void testFcn(int** filter,int filtRows,int filtCols)
{
    int i,j;
    for (i=0;i<filtRows;i++)
    {
        for (j=0;j<filtCols;j++)
        {
            printf("%0.3f ",filter[i][j]);
        }
        printf("\n");
    }
}
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

    for (i=0;i<rows;i++)
    {
        for (j=0;j<cols;j++)
        {
            img[i][j] = (int)(img[i][j]/8);
        }
    }

//////////////////////////////////////////////////////////////////////////////////

    // Initialize prng
    srand(time(NULL));

    int** filteredImg;
    filteredImg = malloc(rows * sizeof(*filteredImg));
    for (i=0;i<rows;i++)
    {
        filteredImg[i] = malloc(cols * sizeof(*filteredImg[i]));
    }

    int** oldfilteredImg;
    oldfilteredImg = malloc(rows * sizeof(*oldfilteredImg));
    for (i=0;i<rows;i++)
    {
        oldfilteredImg[i] = malloc(cols * sizeof(*oldfilteredImg[i]));
    }

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
            filteredImg[i][j] = img[i][j];
    }

    for (i=0;i<rows;i++)
    {
        for (j=0;j<cols;j++)
            oldfilteredImg[i][j] = img[i][j];
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

    float min = 999999999;
    float minEnergy;
    double delta = 0.01;
    float diff;
    double sum;
    double prob;
    float curEnergy;
    float beta;
    int noFilts = 6;

    double probs[32];
    double energies[32];
    double randomNumber;
    double current;
    double max;
    int randInt;
    int backup;
    int count = 0;
    int loopNo = 0;

    int** obsHistogram;
    obsHistogram = malloc(noFilts * sizeof(*obsHistogram));
    for (i=0; i<noFilts; i++)
        obsHistogram[i] = malloc(32 * sizeof(*obsHistogram[i]));

    float** normObsHistogram;
    normObsHistogram = malloc(noFilts * sizeof(*normObsHistogram));
    for (i=0; i<noFilts; i++)
        normObsHistogram[i] = malloc(32 * sizeof(*normObsHistogram[i]));

    int** synthHistogram;
    synthHistogram = malloc(noFilts * sizeof(*synthHistogram));
    for (i=0; i<noFilts; i++)
        synthHistogram[i] = malloc(32 * sizeof(*synthHistogram[i]));

    int** CsynthHistogram;
    CsynthHistogram = malloc(noFilts * sizeof(*CsynthHistogram));
    for (i=0; i<noFilts; i++)
        CsynthHistogram[i] = malloc(32 * sizeof(*CsynthHistogram[i]));

    float** normSynthHistogram;
    normSynthHistogram = malloc(noFilts * sizeof(*normSynthHistogram));
    for (i=0; i<noFilts; i++)
        normSynthHistogram[i] = malloc(32 * sizeof(*normSynthHistogram[i]));

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

    for (i=0;i<noFilts;i++)
    {
        filtCols[i] = filtRows[i];
    }

    double*** filter;
    filter = malloc(noFilts * sizeof(double**));
    for (i=0;i<noFilts;i++)
    {
        filter[i] = malloc(filtRows[i] * sizeof(double*));
        for (j=0;j<filtRows[i];j++)
        {
            filter[i][j] = malloc(filtCols[i] * sizeof(double));
        }
    }

    filter[0][0][0] = 1;
    gabor(filter[1],2,0,filtRows[1],filtCols[1]);
    gabor(filter[2],2,90,filtRows[2],filtCols[2]);
    gabor(filter[3],2,45,filtRows[3],filtCols[3]);
    gabor(filter[4],2,135,filtRows[4],filtCols[4]);
    gaussian(filter[5],2,filtRows[5],filtCols[5]);
    //gaussian(filter[3],2,filtRows[3],filtCols[3]);

    // lambda stuff
    float** lambdas;
    lambdas = malloc(noFilts *sizeof(*lambdas));
    for(i=0;i<noFilts;i++)
        lambdas[i] = malloc(32*sizeof(*lambdas[i]));

    // main code
    computeHistogram(rows,cols,filteredImg,obsHistogram[0]);
    normalizeHistogram(rows,cols,obsHistogram[0],normObsHistogram[0]);
    for (i=1;i<noFilts;i++)
    {
        convolve(rows,cols,filtRows[i],filtCols[i],img,filteredImg,filter[i]);
        computeHistogram(rows,cols,filteredImg,obsHistogram[i]);
        normalizeHistogram(rows,cols,obsHistogram[i],normObsHistogram[i]);
    }

    for(i=0;i<noFilts;i++){
        for (j=0;j<32;j++)
            lambdas[i][j] = 0;
    }
    for (i=0;i<rows;i++){
        for (j=0;j<cols;j++){
            synthImg[i][j] = rand()%32;
        }
    }
    for (i=0;i<rows;i++){
        for (j=0;j<cols;j++)
            imgCopy1[i][j] = synthImg[i][j];
    }

    for (m=0;m<noFilts;m++)
    {
        for (i=0;i<32;i++)
        {
            synthHistogram[m][i] = 0;
        }
    }
    
    computeHistogram(rows,cols,synthImg,synthHistogram[0]);
    normalizeHistogram(rows,cols,synthHistogram[0],normSynthHistogram[0]);
    for (i=1;i<noFilts;i++)
    {
        convolve(rows,cols,filtRows[i],filtCols[i],synthImg,filteredImg,filter[i]);
        computeHistogram(rows,cols,filteredImg,synthHistogram[i]);
        normalizeHistogram(rows,cols,synthHistogram[i],normSynthHistogram[i]);
    }

    for (m=0;m<noFilts;m++)
    {
        sum = 0;
        for (i=0;i<32;i++)
        {
            printf("m: %d i: %d\tObs: %d\t Synth: %d\n",m,i,obsHistogram[m][i], synthHistogram[m][i]);
            sum += synthHistogram[m][i];
            //printf("m: %d i: %d\tObs: %f\t Synth: %f\n",m,i,normObsHistogram[m][i], normSynthHistogram[m][i]);
        }
        printf("sum: %f\n\n",sum);
    }

    count = 0;
    float temp = 32;
    float dTemp = 0.05;

    convolve(rows,cols,filtRows[4],filtCols[4],img,filteredImg,filter[4]);
    //do{
    //    loopNo++;
    //    for (m=0;m<noFilts;m++)
    //    {
    //        for (i=0;i<32;i++){
    //            //lambdas[m][i] = lambdas[m][i] + delta * (synthHistogram[m][i] - obsHistogram[m][i]);
    //            //lambdas[m][i] = lambdas[m][i] + delta * (normSynthHistogram[m][i] - normObsHistogram[m][i]);
    //            lambdas[m][i] = lambdas[m][i] + temp * (normSynthHistogram[m][i] - normObsHistogram[m][i]);
    //        }
    //    }

    //    if (temp > 1)
    //        temp = temp * (1-dTemp);
    //    else
    //        temp = 1;
    //    
    //    curEnergy = 0;
    //    for (m=0;m<noFilts;m++)
    //    {
    //        curEnergy += energy(noFilts,lambdas[m],normSynthHistogram[m]);
    //    }

    //    prob = exp(-1 * curEnergy);

    //    for (k=1;k<=4;k+=1){
    //        count += 1;
    //        beta = 3 / log(count+1);

    //        curEnergy = 0;
    //        for (m=0;m<noFilts;m++)
    //        {
    //            curEnergy += energy(noFilts,lambdas[m],normSynthHistogram[m]);
    //        }
    //        
    //        for (i=0;i<rows;i++){
    //            for (j=0;j<cols;j++){

    //                    //GIBBS
    //                    for (m=0;m<noFilts;m++)
    //                    {
    //                        for(l=0;l<32;l++){
    //                            CsynthHistogram[m][l] = synthHistogram[m][l];
    //                        }
    //                    }

    //                    for (l=0;l<32;l++)
    //                    {
    //                        imgCopy1[i][j] = l;
    //                        //localHistogram(synthHistogram[0],synthImg,imgCopy1,filtSize,rows,cols,i,j);
    //                        computeHistogram(rows,cols,filteredImg,synthHistogram[0]);
    //                        normalizeHistogram(rows,cols,synthHistogram[0],normSynthHistogram[0]);
    //                        for (m=1;m<noFilts;m++)
    //                        {
    //                            localConvolve(rows,cols,filtRows[m],filtCols[m],synthImg,oldfilteredImg,filter[m],i,j);
    //                            localConvolve(rows,cols,filtRows[m],filtCols[m],imgCopy1,filteredImg,filter[m],i,j);
    //                            //localHistogram(synthHistogram[m], oldfilteredImg, filteredImg, filtRows[m], rows, cols, i, j);
    //                            computeHistogram(rows,cols,filteredImg,synthHistogram[m]);
    //                            normalizeHistogram(rows,cols,synthHistogram[m],normSynthHistogram[m]);
    //                        }

    //                        current = 0;
    //                        for (m=0;m<noFilts;m++)
    //                        {
    //                            current += energy(noFilts,lambdas[m],normSynthHistogram[m]);
    //                        }
    //                        energies[l] = current;
    //                    }


    //                    minEnergy = energies[0];
    //                    for (l=1;l<32;l++)
    //                    {
    //                        if (energies[l] < minEnergy)
    //                            minEnergy = energies[l];
    //                    }

    //                    sum = 0;
    //                    for (l=0;l<32;l++)
    //                    {
    //                        probs[l] = exp(-(1/beta) * (energies[l]-minEnergy));
    //                        sum = sum + probs[l];
    //                    }

    //                    for (l=0;l<32;l++)
    //                    {
    //                        probs[l] = probs[l]/sum;
    //                    }
    //                    //if (i==23 && j==23)
    //                    //{
    //                    //    printf("beta: %f\n",beta);
    //                    //    for (l=0;l<32;l++)
    //                    //        printf("energies[%d] = %f\tprobs[%d] = %f\n",l,energies[l]-minEnergy,l,probs[l]);
    //                    //    printf("\n");
    //                    //}

    //                    for (l=1;l<32;l++)
    //                    {
    //                        probs[l] = probs[l-1] + probs[l];
    //                    }

    //                    randomNumber = (float)rand()/(float)RAND_MAX;
    //                    //randomNumber = randomNumber/RAND_MAX;

    //                    l = 0;
    //                    while (1)
    //                    {
    //                        if (randomNumber < probs[l])
    //                            break;
    //                        if (l == 31)
    //                            break;
    //                        l++;
    //                    }

    //                    imgCopy1[i][j] = l;
    //                    //localHistogram(synthHistogram[0],synthImg,imgCopy1,filtSize,rows,cols,i,j);
    //                    computeHistogram(rows,cols,filteredImg,synthHistogram[0]);
    //                    normalizeHistogram(rows,cols,synthHistogram[0],normSynthHistogram[0]);
    //                    for (m=1;m<noFilts;m++)
    //                    {
    //                        localConvolve(rows,cols,filtRows[m],filtCols[m],synthImg,oldfilteredImg,filter[m],i,j);
    //                        localConvolve(rows,cols,filtRows[m],filtCols[m],imgCopy1,filteredImg,filter[m],i,j);
    //                        //localHistogram(synthHistogram[m], oldfilteredImg, filteredImg, filtRows[m], rows, cols, i, j);
    //                        computeHistogram(rows,cols,filteredImg,synthHistogram[m]);
    //                        normalizeHistogram(rows,cols,synthHistogram[m],normSynthHistogram[m]);
    //                    }
    //                    synthImg[i][j] = l;




    //                    /////////////////////////////////////////////////////////////////////////////////////////////////
    //                    ////METROPOLIS
    //                    //randInt = rand()%32;
    //                    //imgCopy1[i][j] = randInt;

    //                    //for (m=0;m<noFilts;m++)
    //                    //{
    //                    //    for(l=0;l<32;l++){
    //                    //        CsynthHistogram[m][l] = synthHistogram[m][l];
    //                    //    }
    //                    //}

    //                    //localHistogram(synthHistogram[0], synthImg, imgCopy1, filtSize, rows, cols, i, j);
    //                    //normalizeHistogram(rows,cols,synthHistogram[0],normSynthHistogram[0]);
    //                    //for (m=1;m<noFilts;m++)
    //                    //{
    //                    //    localConvolve(rows,cols,filtRows[m],filtCols[m],synthImg,oldfilteredImg,filter[m],i,j);
    //                    //    localConvolve(rows,cols,filtRows[m],filtCols[m],imgCopy1,filteredImg,filter[m],i,j);
    //                    //    localHistogram(synthHistogram[m], oldfilteredImg, filteredImg, filtRows[m], rows, cols, i, j);
    //                    //    normalizeHistogram(rows,cols,synthHistogram[m],normSynthHistogram[m]);
    //                    //}

    //                    //current = 0;
    //                    //for (m=0;m<noFilts;m++)
    //                    //{
    //                    //    current += energy(noFilts,lambdas[m],normSynthHistogram[m]);
    //                    //}
    //        
    //                    //if (current < curEnergy){
    //                    //    curEnergy = current;
    //                    //    synthImg[i][j] = randInt;
    //                    //}
    //                    //else{
    //                    //    randomNumber = rand();
    //                    //    randomNumber = randomNumber/RAND_MAX;
    //                    //    prob = exp(-(1/beta) * (current - curEnergy));
    //                    //    //prob = exp(curEnergy) / exp(current + curEnergy);
    //                    //    //if (i == 0 && j == 0)
    //                    //    //    printf("%f %f %f \t %f\n ",beta,prob,current-curEnergy,(curEnergy));

    //                    //    if (randomNumber < prob){
    //                    //        curEnergy = current;
    //                    //        synthImg[i][j] = randInt;
    //                    //    }
    //                    //    else{
    //                    //        for (m=0;m<noFilts;m++)
    //                    //        {
    //                    //            for(l=0;l<32;l++){
    //                    //                synthHistogram[m][l] = CsynthHistogram[m][l];
    //                    //            }
    //                    //            normalizeHistogram(rows,cols,synthHistogram[m],normSynthHistogram[m]);
    //                    //        }
    //                    //        imgCopy1[i][j]=synthImg[i][j];                           
    //                    //    }
    //                    //}     
    //            }
    //        }
    //    }
    //    diff = 0;
    //    
    //    if(loopNo%1==0){
    //        for (m=0;m<noFilts;m++)
    //        {
    //            for (i=0;i<32;i++){
    //                diff += absf(normObsHistogram[m][i] - normSynthHistogram[m][i]);
    //            }
    //        }
    //    }

    //    for (m=0;m<noFilts;m++)
    //    {
    //        sum = 0;
    //        for (i=0;i<32;i++)
    //        {
    //            printf("m: %d i: %d\tObs: %d\t Synth: %d\n",m,i,obsHistogram[m][i], synthHistogram[m][i]);
    //            sum += synthHistogram[m][i];
    //            //printf("m: %d i: %d\tObs: %f\t Synth: %f\n",m,i,normObsHistogram[m][i], normSynthHistogram[m][i]);
    //        }
    //        printf("sum: %f\n\n",sum);
    //    }

    //    if(loopNo < 1000 && loopNo > 50){
    //        if(diff<min)
    //            min = diff;
    //    }

    //    printf("\n diff = %f and min = %f and temp = %f and beta = %f and prob = %f\n",diff, min,temp,beta,prob);
    //    printf("------Loop %d -----\n",loopNo);

    //} while ((diff > min || loopNo<=1000) && loopNo<200);

    for (i=0;i<rows;i++){
        for (j=0;j<cols;j++){
            img[i][j] = filteredImg[i][j];
            //img[i][j] = synthImg[i][j];
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

    for (i=0;i<rows;i++)
    {
        for (j=0;j<cols;j++)
        {
            img[i][j] = (int)(img[i][j]*8);
        }
    }
    strcpy(oName,basename);
    strcat(oName,"-filtered-new.csv");
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

double energy(int noFilts, float *lambdas, float* histogram){
    int i;
    float ans=0;
        for (i=0;i<32;i++){
            ans += lambdas[i]*histogram[i];
        }
    return ans;
}

void localHistogram(int* histogram, int** img1, int** img2, int filtSize, int imgRows, int imgCols, int a, int b){

    int x,y;
    int filtRows = filtSize;
    int filtCols = filtSize;
    int tempHist[32];
    int m = (filtRows-1)/2;
    int n = (filtCols-1)/2;

    for (x=0;x<32;x++)
        tempHist[x] = 0;
    
    for(x=a-m; x<a+m; x++){
        for(y=b-n; y<b+n; y++){
            if(x<m || x>=imgRows-m || y<n || y>=imgCols-n)
                continue;
            else if (img1[x][y] < 0)
                tempHist[0]++;
            else if (img1[x][y] > 31)
                tempHist[31]++;
            else
                tempHist[img1[x][y]]++;
            //if(x<m || x>=imgRows-m || y<n || y>=imgCols-n)
            //    continue;
            //if(img1[x][y]<0)
            //    histogram[0]--;
            //else if(img1[x][y]>31)
            //    histogram[31]--;
            //else if(histogram[img1[x][y]] != 0)
            //    histogram[img1[x][y]]--;

            //if(img2[x][y]<0)
            //    histogram[0]++;
            //else if(img2[x][y]>31)
            //    histogram[31]++;
            //else
            //    histogram[img2[x][y]]++;
        }
    } 
    int o;
    for (x=0;x<32;x++)
    {
        o = histogram[x];
        histogram[x] = histogram[x] - tempHist[x];
        if (histogram[x] < 0 || histogram[x] > imgRows * imgCols)
            printf("histogram[x]: %d o: %d tempHist[x]: %d x: %d\n",histogram[x],o,tempHist[x],x);
    }
    for(x=a-m; x<a+m; x++){
        for(y=b-n; y<b+n; y++){
            //if(x<m || x>=imgRows-m || y<n || y>=imgCols-n)
            //    continue;
            //if(img1[x][y]<0)
            //    histogram[0]--;
            //else if(img1[x][y]>31)
            //    histogram[31]--;
            //else if(histogram[img1[x][y]] != 0)
            //    histogram[img1[x][y]]--;

            if(x<m || x>=imgRows-m || y<n || y>=imgCols-n)
                continue;
            else if (img2[x][y] < 0)
                histogram[0]++;
            else if (img2[x][y] > 31)
                histogram[31]++;
            else
                histogram[img2[x][y]]++;
        }
    } 
}

void computeHistogram(int imgRows, int imgCols, int** img, int* histogram){
    int i,j,k;

    for (i=0;i<32;i++)
        histogram[i] = 0;

    for (i=0;i<imgRows;i++){
        for (j=0;j<imgCols;j++){
            if (img[i][j] < 0)
                histogram[0]++;
            else if (img[i][j] > 31)
                histogram[31]++;
            else
                histogram[img[i][j]]++;
        }
    }
}

void normalizeHistogram(int imgRows, int imgCols, int* histIn, float* histOut)
{
    int i;
    for (i=0;i<32;i++)
    {
        histOut[i] = histIn[i] / (float)(imgRows * imgCols);
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

    /*double min = filter[0][0];
    double sum = 0;

    for(x=0; x<filtRows; x++){
        for(y=0; y<filtCols; y++){
            if(filter[x][y]< min)
                min = filter[x][y];
        }
    }
    for(x=0 ; x<filtRows; x++){
        for(y=0; y<filtCols; y++){
            filter[x][y] = filter[x][y]-min;
            sum += filter[x][y];
        }
    }
    for(x=0 ; x<filtRows; x++){
        for(y=0; y<filtCols; y++){
            filter[x][y] = filter[x][y]/sum;
        }
    }*/
}

double absf(double x)
{
    return x >= 0 ? x : -1 * x;
}
