#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

const int MAXBIT = 128, MAXDICT = 2048, MAXFLOAT = 32;

void getSubString(char *, char *, int, int);
void jointSample(FILE *, float ,float );
int getDict(float *, char *, char *);
void strCpyArray(char *, char *, int, int);
void decompressSignal(char *, char *, char*);
void bytetostr(char *, int);
int bytetosignalvalue(char *, char *, char *, int);

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
	    decompressSignal(fileName,fileNameRes,fileNameDict);

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

void jointSample(FILE *fileResult, float imag, float real)
{
	if(imag >= 0)
	{
		fprintf(fileResult, "%f+%f\n", real,imag);
	}
	else
	{
		fprintf(fileResult, "%f%f\n", real,imag);
	}
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

void decompressSignal(char *fileSig, char *filenameResult, char *fileNameDict)
{
	int flagEnd = 0, ch, n = 0, idx = -1, dictSize, flagIQ = 0;
   	char bitsCode[MAXDICT][MAXBIT], newStr[MAXBIT], strrem[MAXBIT];
   	unsigned byte;
   	float values[MAXDICT], real, imag;
    FILE *file, *fileResult, *fileDict;

    strrem[0] = '\0';

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
    file = fopen(fileSig, "r");

    printf("\n");
    printf("---------------------------------------------------------------\n");
    printf("------------------------  GETTING FLAG  -----------------------\n");
    printf("---------------------------------------------------------------\n");
    printf("\n");

    while (ch = fgetc(file)) 
    {
    	printf("Char value: %x\n", ch);
	    bytetostr(newStr, ch);
	    printf("===============================================================================\n");
	    printf("Decimal value: \"%d\", Binary value: \"%s\"\n", ch, newStr);
	    idx = bytetosignalvalue(strrem, newStr, bitsCode, dictSize);
	    n++;

	    if (-1 != idx)
	    {
	    	break;
	    }
   	}

   	printf("\n");
	printf("---------------------------------------------------------------\n");
	printf("---------------------  PROCESSING BYTES  ----------------------\n");
	printf("---------------------------------------------------------------\n");
	printf("\n");

    
    fileResult = fopen(filenameResult, "w");
    while ((ch = fgetc(file))) 
    {
    	printf("Char value: %x\n", ch);
    	bytetostr(newStr, ch);
    	idx = bytetosignalvalue(strrem, newStr, bitsCode, dictSize);

    	while(idx != -1 && idx != dictSize-1)
    	{
    		if(flagIQ == 0 && idx != dictSize-1)
    		{
    			flagIQ = 1;
    			real = values[idx];
    		}
    		 else if(flagIQ == 1 && idx != -1)
    		{
    			flagIQ = 0;
    			imag = values[idx];
    			printf("Real part: \"%f\", Imaginary part: \"%f\"\n", real, imag);
    			jointSample(fileResult, imag, real);
    		}
     	 	idx = bytetosignalvalue(strrem, "\0", bitsCode, dictSize);
    	}
    	
    	if (idx == dictSize-1)
    	{
    		break;
    	}

    	n++;
   	}
}

void bytetostr(char *dest, int byte)
{
	char strtemp[8];
	int i = 0;

	while(byte != 0)
	{

		if(byte % 2 == 0)
		{
			strtemp[i++] = '0';
			
		}
		else
		{
			strtemp[i++] = '1';
			byte--;
		}
		byte /= 2;

	}

	while(strlen(strtemp) < 8)
	{
		strtemp[i++] = '0';
	}

	strtemp[i] = '\0';
	strrev(strtemp);
	strcpy(dest,strtemp);
}

int bytetosignalvalue(char *strrem, char *currentstr, char *bitsCode, int dictSize)
{
	int startBit = strlen(strrem), j = 1, i;
	char strtemp[MAXBIT], tempbit[MAXBIT];

	strcat(strrem, currentstr);
	printf("Actual bit code: \"%s\", Length: \"%d\"\n", strrem, strlen(strrem));


	while(j - strlen(strrem) != 0)
	{
		getSubString(strtemp,strrem,0,j);

		for(i = 0; i < dictSize; i++)
		{
			strCpyArray(tempbit, bitsCode, 0, i*MAXBIT);
			if(strcmp(strtemp,tempbit) == 0)
			{
				break;
			}
		}

		if(i < dictSize)
		{
			getSubString(strrem,strrem,j,strlen(strrem));	
			printf("Dictionary index: \"%d\", Bit code: \"%s\", Bits remaining: \"%s\"\n", i+1, tempbit, strrem);
			printf("===============================================================================\n");
			return i;
		}
		j++;
	}
	printf("===============================================================================\n");
	return -1;	
}