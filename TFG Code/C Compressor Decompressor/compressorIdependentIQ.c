#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

const int MAXBIT = 128, MAXDICT = 2048, MAXFLOAT = 32;

void compressSignal(char *, char *, char *);
void splitSample(char *, char *, char *);
void getSubString(char *, char *, int, int);
int getDict(float *, char *, char *);
void quantitzationIndependent(float *, char *, float, float, int, char *);
void strCpyArray(char *, char *, int, int);
int sampletobin(char *, char *, char *, int, FILE *);

int main(int argc, char *argv[])
{
	char fileName[128],fileNameRes[128],fileNameDict[128];
	
	if(argc == 4){
		
		strcpy(fileName, argv[1]);
		printf("Singnal file: %s\n", fileName);
		strcpy(fileNameRes, argv[2]);
		printf("Results file: %s\n", fileNameRes);
	    strcpy(fileNameDict,argv[3]);
		printf("Dictionary fiel: %s\n", fileNameDict);
	    compressSignal(fileName,fileNameRes,fileNameDict);

	}else{
		printf("error");
	}
    
    return 0;
}

void getSubString(char *dest, char *string, int start, int end)
{
	char result[MAXBIT];
	int j = 0;
	for (int i = start; i < end; ++i)
	{
		result[j++] = string[i];
	}

	result[j] = '\0';
    strcpy(dest,result);
}

void compressSignal(char *fileSig, char *filenameResult, char *fileNameDict)
{
    
    int dictSize, i = 0, leftintbin = 0;
   	char sampleR[MAXFLOAT],sampleI[MAXFLOAT], line[64], bitsCode[MAXDICT][MAXBIT], samplestr[2*MAXBIT], binstr[32], leftstrbin = '\0';
   	float values[MAXDICT], real, imag;
    FILE *file, *fileResult, *fileDict;

	printf("\n");
    printf("--------------------------------------------------------------\n");
    printf("-------------------  READING DICTIONARY  ---------------------\n");
    printf("--------------------------------------------------------------\n");
    printf("\n");

    dictSize = getDict(&values[0], bitsCode, fileNameDict);
    for(int j = 0; j < dictSize; j++)
    {
    	printf("Value float: \"%f\", Value binary: \"%s\"\n",values[j], bitsCode[j]);
    }
    fileResult = fopen(filenameResult, "w");

    printf("\n");
    printf("---------------------------------------------------------------\n");
    printf("----------------------  WRITTING FLAG  ------------------------\n");
    printf("---------------------------------------------------------------\n");
    printf("\n");

    leftintbin = sampletobin(binstr, bitsCode[dictSize-1], &leftstrbin, leftintbin, fileResult);

    printf("\n");
    printf("---------------------------------------------------------------\n");
    printf("-------------------  PROCESSING SAMPLES  ----------------------\n");
    printf("---------------------------------------------------------------\n");
    printf("\n");

    file = fopen(fileSig, "r");
    while(fgets(line, sizeof(line), file))
    {
        splitSample(line,sampleI,sampleR);
        imag = strtof(sampleI,NULL);
        real = strtof(sampleR,NULL);
        printf("\n");
		printf("sample: %d\n", i+1);
        printf("Real sample: %f, Imag sample: %f\n", real, imag);

        /*quantizer + huffman*/

        quantitzationIndependent(&values[0], bitsCode, real, imag, dictSize, samplestr);

        /*string to binary*/
        leftintbin = sampletobin(binstr, samplestr, &leftstrbin, leftintbin, fileResult);
        //printf("Result value: %s, Value ramining bit: %d, Num. bits remaining: %d\n", binstr, leftstrbin, leftintbin);
        /*-----------------------------------------------------------------------------------------*/
        //jointSample(fileResult, imag, real);
		i++;
		/*if(i == 10)
		{
			break;
		}*/
    }
    leftintbin = sampletobin(binstr, bitsCode[dictSize-1], &leftstrbin, leftintbin, fileResult);
    printf("\n");
    printf("---------------------------------------------------------------\n");
    printf("--------------------  FILLING LAST BYTE  ----------------------\n");
    printf("---------------------------------------------------------------\n");
    printf("\n");
    if(leftintbin != 0)
    {
    	leftstrbin = leftstrbin << (8-leftintbin);
    	printf("Final value filling last byte: \"%c\"(%d)\n", leftstrbin, leftstrbin);
    	fprintf(fileResult, "%c", leftstrbin);
    }
    fclose(file);
    fclose(fileResult);
}

void splitSample(char *sampleStr, char *imagS, char *realS)
{

	int count = strlen(sampleStr);
	int i = 0;
	char samChar;
	if(sampleStr[0] == '-')
	{
		i = 1;
	}
	do
	{
		samChar = sampleStr[i++];
	}
	while(i <= count && samChar != '-' && samChar != '+');
    getSubString(realS,sampleStr,0,i-1);
	if(sampleStr[i-1] == '+')
		i++;
	getSubString(imagS,sampleStr,i-1,count-1);
}

int getDict(float *values, char *bitsCode, char *fileNameDict)
{
	FILE *fileDict;
	int j = 0;
	char sampleBit[MAXBIT],sampleVal[MAXFLOAT], line[128];
	int count;

	fileDict = fopen(fileNameDict, "r");
	while(fgets(line, sizeof(line), fileDict)){
        count = strlen(line);
		int i = 0;
		char dictChar;
		do
		{
			dictChar = line[i++];
		}
		while(i <= count && dictChar != ' ');
	    getSubString(sampleVal,line,0,i-1);
		getSubString(sampleBit,line,i,count);
		values[j] = strtof(sampleVal,NULL);
		strCpyArray(bitsCode, sampleBit, j*MAXBIT, 0);
		j++;
		/*printf("%s\n", sampleBit);
		printf("value: %f bit: %c\n", values[0], bitsCode[0]);
		printf("%s\n", line);*/
		
    }
    fclose(fileDict);

    return j;
}

void quantitzationIndependent(float *values, char *bitsCode, float sampleR, float sampleI, int dictSize, char *sampleStr)
{

	float differencesR, differencesI, minR =1000, minI = 1000;
	int idxI, idxR;
	char bitI[MAXBIT], bitR[MAXBIT];

	for(int i = 0; i < dictSize; i++)
	{
		differencesR = fabs(sampleR - values[i]);
		differencesI = fabs(sampleI - values[i]);

		//printf("%f\n", values[i]);
		//printf("Real difference: %f, Imag difference %f\n", differencesR, differencesI);

		if(differencesR < minR)
		{
			idxR = i;
			minR = differencesR;
		}
		if(differencesI < minI)
		{
			idxI = i;
			minI = differencesI;
		}
	}

	strCpyArray(bitR, bitsCode, 0, idxR*MAXBIT);
	strCpyArray(bitI, bitsCode, 0, idxI*MAXBIT);
	printf("===============================================================================\n");
	printf("Real: %d(%f), Imag: %d(%f)\n", idxR, values[idxR], idxI, values[idxI]);
	printf("Length bit real: %d, Length bit imag: %d\n", strlen(bitR), strlen(bitI));
	printf("Real bit: \"%s\", Imag bit: \"%s\"\n", bitR, bitI);
	printf("===============================================================================\n");

	strcat(bitR,bitI);
	strcpy(sampleStr,bitR);
}

void strCpyArray(char *dest, char *src, int startDest, int startSrc)
{

	int k = 0;
	while(src[startSrc+k] == '1' || src[startSrc+k] == '0')
	{
		dest[startDest+k] = src[startSrc+k];
		k++;
	}
	

	dest[startDest+k] = '\0';
}

int sampletobin(char *binstr, char *samplestr,char *leftstrbin, int leftintbin, FILE *fileResult)
{
	int lensample = strlen(samplestr), i = 0, j = 0, k = 0;
	char byte = *leftstrbin, tempstr[32];
	printf("===============================================================================\n");
	printf("Value current sample: \"%s\", Value ramining bit: \"%d\", Num. bits remaining: \"%d\"\n", samplestr, *leftstrbin, leftintbin);

	while(lensample-i > 0)
	{

		for(j = leftintbin; j < 8 && i < lensample; j++)
		{
			

			byte = byte << 1;
			if(samplestr[i] == '1')
			{
				byte++;
			}

			printf("Value byte: %d, Value actual bin: %c, Bit num: %d, Position current byte: %d\n", byte, samplestr[i], i, j+1);
			i++;

		}	
		
		leftintbin = j;
		if(j == 8)
		{
			tempstr[k++] = byte;
			leftintbin = 0;
			*leftstrbin = '\0';
			printf("Bit complete value: %d\n", tempstr[k-1]);
			fprintf(fileResult, "%c", byte);
			
		}else
		{
			*leftstrbin = byte;
		}

		byte = '\0';
	}

	strcpy(binstr, tempstr);
	printf("Result value: \"%s\", Value ramining bit: \"%d\", Num. bits remaining: \"%d\"\n", binstr, *leftstrbin, leftintbin);
	printf("===============================================================================\n");
	return leftintbin;
}