#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>	// includes sin
#include <string>
#include <fstream>
#include "functions.h"

//code from numericalrecipes in c textbook provided by professor
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

//  Standard sample rate in Hz
#define SAMPLE_RATE     	44100.0

//  Standard sample size in bits
#define BITS_PER_SAMPLE		16

// Standard sample size in bytes		
#define BYTES_PER_SAMPLE	(BITS_PER_SAMPLE/8)

// Rescaling factor to convert between 16-bit shorts and doubles between -1 and 1
#define MAX_SHORT_VALUE		32768

// Number of channels
#define MONOPHONIC			1

// Offset of the fmt chunk in the WAV header
#define FMT_OFFSET			12

using namespace std;


/*
Writes the header for a WAV file with the given attributes to 
 the provided filestream
*/
//code from testtone.cpp, given to us by the professor with permission to use 
void writeWavFileHeader(int channels, int numberSamples, double outputRate, FILE *outputFile) {
    // Note: channels is not currently used. You will need to add this functionality
	// yourself for the bonus part of the assignment
	
	/*  Calculate the total number of bytes for the data chunk  */
    int dataChunkSize = numberSamples * BYTES_PER_SAMPLE;
	
    /*  Calculate the total number of bytes for the form size  */
    int formSize = 36 + dataChunkSize;
	
    /*  Calculate the total number of bytes per frame  */
    short int frameSize = channels * BYTES_PER_SAMPLE;
	
    /*  Calculate the byte rate  */
    int bytesPerSecond = (int)ceil(outputRate * frameSize);

    /*  Write header to file  */
    /*  Form container identifier  */
    fputs("RIFF", outputFile);
      
    /*  Form size  */
    fwriteIntLSB(formSize, outputFile);
      
    /*  Form container type  */
    fputs("WAVE", outputFile);

    /*  Format chunk identifier (Note: space after 't' needed)  */
    fputs("fmt ", outputFile);
      
    /*  Format chunk size (fixed at 16 bytes)  */
    fwriteIntLSB(16, outputFile);

    /*  Compression code:  1 = PCM  */
    fwriteShortLSB(1, outputFile);

    /*  Number of channels  */
    fwriteShortLSB((short) MONOPHONIC, outputFile);

    /*  Output Sample Rate  */
    fwriteIntLSB((int)outputRate, outputFile);

    /*  Bytes per second  */
    fwriteIntLSB(bytesPerSecond, outputFile);

    /*  Block alignment (frame size)  */
    fwriteShortLSB(frameSize, outputFile);

    /*  Bits per sample  */
    fwriteShortLSB(BITS_PER_SAMPLE, outputFile);

    /*  Sound Data chunk identifier  */
    fputs("data", outputFile);

    /*  Chunk size  */
    fwriteIntLSB(dataChunkSize, outputFile);
}


/*
Creates a WAV file with the contents of the provided outputArray as the samples, and writes
it to the given filename
 */
//code from testtone.cpp, given to us by the professor with permission to use 
void writeWavFile(double *outputArray, int outputArraySize, int channels, char *filename) {
    // Note: channels is not currently used. You will need to add this functionality
	// yourself for the bonus part of the assignment

  //open a binary output file stream for writing
    FILE *outputFileStream = fopen(filename, "wb");
    if (outputFileStream == NULL) {
      printf("File %s cannot be opened for writing\n", filename);
        return;
    }

    //create an 16-bit integer array to hold rescaled samples
    short *intArray = new short[outputArraySize];

    //find the largest entry and uses that to rescale all other
    // doubles to be in the range (-1, 1) to prevent 16-bit integer overflow
    double largestDouble = 1;
    for (int i=0; i< outputArraySize; i++) {
		if (outputArray[i] > largestDouble) {
			largestDouble = outputArray[i];
		}
    }

    for (int i=0; i<outputArraySize; i++) {
		intArray[i] = (short) ((outputArray[i]/largestDouble)*MAX_SHORT_VALUE);
    }
	
    int numSamples = outputArraySize;

	// actual file writing
    writeWavFileHeader(channels, numSamples, SAMPLE_RATE, outputFileStream);
    fwrite(intArray, sizeof(short), outputArraySize, outputFileStream);
    
    //clear memory from heap
    delete[] intArray;
}


//writes an integer to the provided stream in little-endian form
//code from testtone.cpp, given to us by the professor with permission to use 
size_t fwriteIntLSB(int data, FILE *stream) {
    unsigned char array[4];

    array[3] = (unsigned char)((data >> 24) & 0xFF);
    array[2] = (unsigned char)((data >> 16) & 0xFF);
    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 4, stream);
}


//reads an integer from the provided stream in little-endian form
//code from testtone.cpp, given to us by the professor with permission to use 
int freadIntLSB(FILE *stream) {
    unsigned char array[4];

    fread(array, sizeof(unsigned char), 4, stream);
    
    int data;
    data = array[0] | (array[1] << 8) | (array[2] << 16) | (array[3] << 24);

    return data;
}


//writes a short integer to the provided stream in little-endian form
//code from testtone.cpp, given to us by the professor with permission to use 
size_t fwriteShortLSB(short int data, FILE *stream) {
    unsigned char array[2];

    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 2, stream);
}


//reads a short integer from the provided stream in little-endian form
//code from testtone.cpp, given to us by the professor with permission to use 
short int freadShortLSB(FILE *stream) {
    unsigned char array[2];

    fread(array, sizeof(unsigned char), 2, stream);
    
    int data;
    data = array[0] | (array[1] << 8);

    return data;
}

//code from the T04 lecture recording
void readWavFileHeader(int *channels, int *numSamples, FILE *inputFile){

    int sampleRate;
    int bytesPerSecond;
    int dataChunkSize;

    unsigned char buffer[64];
    fread(buffer, sizeof(unsigned char), FMT_OFFSET, inputFile);

    freadIntLSB(inputFile);
    int fmtSize = freadIntLSB(inputFile);
    freadShortLSB(inputFile);

    *channels = freadShortLSB(inputFile);
    sampleRate = freadIntLSB(inputFile);
    bytesPerSecond = freadIntLSB(inputFile);

    int frameSize = freadShortLSB(inputFile);
    int bitRate = freadShortLSB(inputFile);

    if(bitRate != BITS_PER_SAMPLE){
        printf("Error: bit rate provided is not 16. exiting.");
        exit(-1);
    }

    if(sampleRate != SAMPLE_RATE){
        printf("Error: sample rate is not 44.1 khz. exiting.");
        exit(-1);
    }

    fread(buffer, sizeof(unsigned char), fmtSize-12, inputFile);

    dataChunkSize = freadIntLSB(inputFile);
    printf("Data chunk size: %d\n", dataChunkSize);

    *numSamples = dataChunkSize / (BYTES_PER_SAMPLE * (*channels));

}

//code from the t04 lecture recording
double* readWavFile(int *arraySize, int *channels, char *filename){

    double *array;
    FILE *inputFileStream = fopen(filename, "rb");

    if(inputFileStream == NULL){
        printf("File %s could not be opened for reading\n", filename);
        exit(-1);
    }

    int numSamples; 

    readWavFileHeader(channels, &numSamples, inputFileStream);

    printf("Channels: %d\n", *channels);
    printf("Samples: %d\n", numSamples);

    if(numSamples <= 0){
        printf("the file %s does not have any samples. exiting.\n", filename);
        exit(0);
    }

    *arraySize = numSamples * (*channels);

    array = new double[*arraySize];
    short *intArray = new short[*arraySize];

    int count = fread(intArray, BYTES_PER_SAMPLE, numSamples, inputFileStream);


    int largest  = 0;

    for (int i = 0; i < *arraySize; i++)
    {
        if(abs(intArray[i]) > largest){
            largest = intArray[i];
        }
    }
    for (int i = 0; i < *arraySize; i++)
    {
        array[i] = ((double) intArray[i]) / largest;
    }

    delete [] intArray;

    return array;
    
}

void convolve(char *inputFileName, char *irFileName, char *outputFileName){

    int inputArraySize;
    int irArraySize;
    int outputArraySize;
    int channels;

    double *inputArray;
    double *irArray;
    //reading
    inputArray = readWavFile(&inputArraySize, &channels, inputFileName);
    irArray = readWavFile(&irArraySize, &channels, irFileName);

    outputArraySize = inputArraySize + irArraySize -1;

    double *output = new double[outputArraySize];
    //writing
    for(int i = 0; i < irArraySize; i++){
        for(int j = 0; j < inputArraySize; j++){
            output[i+j] += irArray[i] * inputArray[j];
        }
    }

    double largest = 0.0;
    //scaling
    for (int i = 0; i < outputArraySize; i++)
    {
        if(abs(output[i]) > largest){
            largest = output[i];
        }
    }
    for (int i = 0; i < outputArraySize; i++)
    {
        output[i] = ((double) output[i]) / largest;
    }


    writeWavFile(output, outputArraySize, channels, outputFileName);

}
//code from numericalrecipes in c textbook provided by professor
void four1(double data[], unsigned long nn, int isign){

    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    //printf("inside four1\n");

    n = nn << 1;
    j = 1;

    for(i=1; i < n; i+=2 ){
        if(j > i){
            SWAP(data[j],data[i]);
            SWAP(data[j+1],data[i+1]);
        }
        m = nn;
        while(m >= 2 && j > m){
            j -=m;
            m >>= 1;
        }
        j+=m;
    }
    mmax = 2;
    while(n > mmax){
        istep = mmax << 1;
        theta=isign*(6.28318530717959/mmax);
        wtemp=sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi=sin(theta);
        wr=1.0;
        wi=0.0;
        for (m=1;m<mmax;m+=2){
            for (i=m;i<=n;i+=istep){
                j=i+mmax;
                tempr=wr*data[j]-wi*data[j+1];
                tempi=wr*data[j+1]+wi*data[j];
                data[j]=data[i]-tempr;
                data[j+1]=data[i+1]-tempi;
                data[i] += tempr;
                data[i+1] += tempi;

            }
            wr=(wtemp=wr)*wpr-wi*wpi+wr;
            wi=wi*wpr+wtemp*wpi+wi;
        }
        mmax=istep;
    }
    //printf("last line of four1\n");
}
unsigned int log2(unsigned int x){
    if(x < 2) {return 0;}
    else if(x < 4) {return 1;}
    else if(x < 8) {return 2;}
    else if(x < 16) {return 3;}
    else if(x < 32) {return 4;}
    else if(x < 64) {return 5;}
    else if(x < 128) {return 6;}
    else if(x < 256) {return 7;}
    else if(x < 512) {return 8;}
    else if(x < 1024) {return 9;}
    else if(x < 2048) {return 10;}
    else if(x < 4096) {return 11;}
    else if(x < 8192) {return 12;}
    else if(x < 16384) {return 13;}
    else if(x < 32768) {return 14;}
    else if(x < 65536) {return 15;}
    else if(x < 131072) {return 16;}
    else if(x < 262144) {return 17;}
    else if(x < 524288) {return 18;}
    else if(x < 1048576) {return 19;}
    else if(x < 2097152) {return 20;}
    else if(x < 4194304) {return 21;}
    else if(x < 8388608) {return 22;}
    else if(x < 16777216) {return 23;}
    else if(x < 33554432) {return 24;}
    else if(x < 67108864) {return 25;}
    else if(x < 134217728) {return 26;}
    else if(x < 268435456) {return 27;}
    else if(x < 536870912) {return 28;}
    else if(x < 1073741824) {return 29;}
    else if(x < 2147483648) {return 30;}
    else {return 31;}

}



void convolveFFT(char *inputFileName, char *irFileName, char *outputFileName){
    int inputArraySize;
    int irArraySize;
    int outputArraySize;
    int channels;
    int complexInputArraySize;
    int complexIrArraySize;


    double *inputArray;
    double *irArray;
    //reading
    inputArray = readWavFile(&inputArraySize, &channels, inputFileName);
    irArray = readWavFile(&irArraySize, &channels, irFileName);
    //printf("finished reading\n");

    complexInputArraySize = 2*inputArraySize;
    complexIrArraySize = 2*irArraySize;

    double *complexInputArray = new double[complexInputArraySize];
    double *complexIrArray = new double[complexIrArraySize];
    //printf("before for loop\n");
    //creating complex arrays
    for(int i = 0; i < complexInputArraySize; i++){
       // printf("in for loop\n");
        complexInputArray[i] = inputArray[i/2];
        complexInputArray[i+1] = 0.0;
    }
    //printf("after first for loop\n");
    for (int i = 0; i < complexIrArraySize; i++)
    {
        //printf("in seconf for loop\n");
       complexIrArray[i] = irArray[i/2];
       complexIrArray[i+1] = 0.0;
    }
    //printf("finished making complex arrays.\n");

    int next;
    //padding
    if(complexInputArraySize > complexIrArraySize){
         next = pow(2, log2(complexInputArraySize));
    }else{
        next = pow(2, log2(complexIrArraySize));
    }

    //printf("next= %d \n", next);
    double *complexFinalInputArray = new double[next];
    double *complexFinalIrArray = new double[next];

    for (int i = 0; i < complexInputArraySize; i++)
    {
        complexFinalInputArray[i] = complexInputArray[i];
    }

    for (int i = 0; i < complexIrArraySize; i++)
    {
        complexFinalIrArray[i] = complexIrArray[i];
    }
    
    //printf("got to line 419\n");
    //performing fft
    four1(complexFinalInputArray, next/2, 1);
    four1(complexFinalIrArray, next/2, 1);

    //printf("got past four1 functions\n");
    //writing
    double *outputArrayComplex = new double[next];
    for(int i = 0; i < next; i++){

        outputArrayComplex[i] = complexFinalInputArray[i] * complexFinalIrArray[i];
    }

    //printf("got to inverse four1 call\n");
    //inverse fft 
    four1(outputArrayComplex, next/2, -1);
    //printf("got past inverse four1 call\n");
    double *outputArrayComplex2 = new double[next/2];
    //unpadding
    for(int i = 0; i < next/2; i++){

        outputArrayComplex2[i] = outputArrayComplex[2*i];
    }

    outputArraySize = inputArraySize + irArraySize - 1;

    double *outputArray = new double[outputArraySize];
    double largest = 0.0;

    for(int i = 0; i < outputArraySize; i++){

        outputArray[i] =  outputArrayComplex2[i];
        if(abs(outputArray[i]) > largest){
            largest = outputArray[i];
        }
    }
    
    //scaling
    for (int i = 0; i < outputArraySize; i++)
    {
        outputArray[i] = ((double) outputArray[i]) / largest;
    }

    writeWavFile(outputArray, outputArraySize, channels, outputFileName);

}

int main(int argc, char** argv){

    char *inputFileName;
    char *irFileName;
    char *outputFileName;

    if(argc != 4){
        printf("Please enter as input the input file name, ir file name, and output file name; in that order.");
        exit(-1);
    }
    inputFileName = argv[1];
    irFileName = argv[2];
    outputFileName = argv[3];

    printf("Convolving.\n");
    //uncomment the line below and comment out line 500 to run the convolve program without the fft
    //convolve(inputFileName, irFileName, outputFileName);
    convolveFFT(inputFileName, irFileName, outputFileName);
    printf("finished.");
    exit(0);
}