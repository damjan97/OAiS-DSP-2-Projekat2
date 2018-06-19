#include "ImageFilter.h"


void convolve2D (uchar Y_buff[], int xSize, int ySize, double filterCoeff, int N)
{
	//TO DO
}

void extendBorders(const uchar input[], int xSize, int ySize, uchar output[], int delta)
{
	memset(output, 0, (xSize + 2 * delta) * (ySize + 2 * delta) * 3);
	
	for (int j = delta; j < ySize + delta; j++)
	{
		for (int i = delta; i < xSize + delta; i++)
		{
			output[j * (xSize + 2 * delta) * 3 + i * 3]     = input[(j - delta)*xSize * 3 + (i - delta) * 3];
			output[j * (xSize + 2 * delta) * 3 + i * 3 + 1] = input[(j - delta)*xSize * 3 + (i - delta) * 3 + 1];
			output[j * (xSize + 2 * delta) * 3 + i * 3 + 2] = input[(j - delta)*xSize * 3 + (i - delta) * 3 + 2];
		}
	}
}
	
void performNFFilter (uchar input[], int xSize, int ySize)
{
	//TO DO
}

void performVFFilter (uchar input[], int xSize, int ySize)
{
	//TO DO
}

void performSuccessiveVFFilter (uchar input[], int xSize, int ySize, int stages)
{
	//TO DO
}

void performSobelEdgeDetection(uchar input[], int xSize, int ySize, uchar threshold)
{
	//TO DO
}

void performNFplusSobelEdgeDetection(uchar input[], int xSize, int ySize, int stages, uchar threshold)
{
	//TO DO
}
