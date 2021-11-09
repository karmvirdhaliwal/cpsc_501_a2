#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>	// includes sin
#include <string>
#include <fstream>
#include "functions.h"


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
size_t fwriteIntLSB(int data, FILE *stream) {
    unsigned char array[4];

    array[3] = (unsigned char)((data >> 24) & 0xFF);
    array[2] = (unsigned char)((data >> 16) & 0xFF);
    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 4, stream);
}


//reads an integer from the provided stream in little-endian form
int freadIntLSB(FILE *stream) {
    unsigned char array[4];

    fread(array, sizeof(unsigned char), 4, stream);
    
    int data;
    data = array[0] | (array[1] << 8) | (array[2] << 16) | (array[3] << 24);

    return data;
}


//writes a short integer to the provided stream in little-endian form
size_t fwriteShortLSB(short int data, FILE *stream) {
    unsigned char array[2];

    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 2, stream);
}


//reads a short integer from the provided stream in little-endian form
short int freadShortLSB(FILE *stream) {
    unsigned char array[2];

    fread(array, sizeof(unsigned char), 2, stream);
    
    int data;
    data = array[0] | (array[1] << 8);

    return data;
}

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

    inputArray = readWavFile(&inputArraySize, &channels, inputFileName);
    irArray = readWavFile(&irArraySize, &channels, irFileName);

    outputArraySize = inputArraySize + irArraySize -1;

    double *output = new double[outputArraySize];

    for(int i = 0; i < irArraySize; i++){
        for(int j = 0; j < inputArraySize; j++){
            output[i+j] += irArray[i] * inputArray[j];
        }
    }

    double largest = 0.0;

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

    printf("Convolving.");
    convolve(inputFileName, irFileName, outputFileName);
    printf("finished.");
    exit(0);
}