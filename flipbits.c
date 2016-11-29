#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char* argv[])
{
    if (argc < 5)
        return 1;
    int i = 0, j = 0;
    char* fName = argv[1];
    char* oName = argv[2];
    int x = atoi(argv[3]);
    int y = atoi(argv[4]);

    size_t count;
    char buffer[5000];
    char *line;
    char *record;

    int img[y][x];

    FILE *file;
    file = fopen(fName,"r");

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

    for (i=0; i<y; i++)
    {
        for (j=0; j<x; j++)
            img[i][j] = 255-img[i][j];
    }

    file = fopen(oName,"w");

    for (i=0;i<y;i++)
    {
        for (j=0;j<x;j++)
            fprintf(file,"%d,",img[i][j]);
        fprintf(file,"\n");
    }

    fclose(file);
    return 0;
}
