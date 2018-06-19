#include "ImageInterpolation.h"
#include "ColorSpaces.h"
#include "ImageFilter.h"
#include <math.h>
#include <ctime>

void clip(int* x, int* y, int xMax, int yMax)
{
	if (*x > xMax)
		*x = xMax;
	if (*y > yMax)
		*y = yMax;
}

void sampleAndHold(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	std::clock_t t_start = std::clock();
	double x_f = (double)newXSize / xSize;
	double y_f = (double)newYSize / ySize;
	int i, j;

	for (int q = 0; q < newYSize; q++)
	{
		j = floor(q / y_f + 0.5);
		for (int p = 0; p < newXSize; p++)
		{
			i = floor(p / x_f + 0.5);
			clip(&i, &j, xSize - 1, ySize - 1);
			output[q * newXSize * 3 + p * 3]     = input[j * xSize * 3 + i * 3];
			output[q * newXSize * 3 + p * 3 + 1] = input[j * xSize * 3 + i * 3 + 1];
			output[q * newXSize * 3 + p * 3 + 2] = input[j * xSize * 3 + i * 3 + 2];
		}
	}
	std::clock_t t_end = std::clock();
	qDebug("Execution time for sample and hold: %.5f ms", (1000.0 * (t_end - t_start)) / CLOCKS_PER_SEC);
}

void bilinearInterpolate(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	std::clock_t t_start = std::clock();
	double x_f = (double)newXSize / xSize;
	double y_f = (double)newYSize / ySize;
	double a, b;
	int i, j;

	for (int q = 0; q < newYSize; q++)
	{
		j = q / y_f;
		for (int p = 0; p < newXSize; p++)
		{
			a = p / x_f - floor(p / x_f);
			b = q / y_f - floor(q / y_f);
			i = p / x_f;
			clip(&i, &j, xSize - 2, ySize - 2);
			output[q * newXSize * 3 + p * 3]     = (1 - a) * (1 - b) * input[j * xSize * 3 + i * 3]     + (1 - a) * b * input[(j + 1) * xSize * 3 + i * 3]     + a * (1 - b) * input[j * xSize * 3 + (i + 1) * 3]     + a * b * input[(j + 1) * xSize * 3 + (i + 1) * 3];
			output[q * newXSize * 3 + p * 3 + 1] = (1 - a) * (1 - b) * input[j * xSize * 3 + i * 3 + 1] + (1 - a) * b * input[(j + 1) * xSize * 3 + i * 3 + 1] + a * (1 - b) * input[j * xSize * 3 + (i + 1) * 3 + 1] + a * b * input[(j + 1) * xSize * 3 + (i + 1) * 3 + 1];
			output[q * newXSize * 3 + p * 3 + 2] = (1 - a) * (1 - b) * input[j * xSize * 3 + i * 3 + 2] + (1 - a) * b * input[(j + 1) * xSize * 3 + i * 3 + 2] + a * (1 - b) * input[j * xSize * 3 + (i + 1) * 3 + 2] + a * b * input[(j + 1) * xSize * 3 + (i + 1) * 3 + 2];
		}
	}
	std::clock_t t_end = std::clock();
	qDebug("Execution time for bilinear interpolation: %.5f ms", (1000.0 * (t_end - t_start)) / CLOCKS_PER_SEC);
}

double wFunction(double d)
{
	double ret = 0.0;
	double fabs_tmp = fabs(d);
	if (fabs_tmp < 1)
		ret = 1.5 * fabs_tmp * d * d - 2.5 * d * d + 1;
	else if (fabs_tmp >= 1 && fabs_tmp < 2)
		ret = -0.5 * fabs_tmp * d * d + 2.5 * d * d - 4 * fabs_tmp + 2;
	return ret;
}

uchar cubicInterpolate(uchar x[4], double d)
{
	double y = 0;
	for (int i = 0; i < 4; i++)
		y += x[i] * wFunction((i < 2) ? (1 - i + d) : (i - 1 - d));
	if (y > 255)
		return 255;
	else if (y < 0)
		return 0;
	else
		return (uchar)y;
}

void bicubicInterpolate(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	std::clock_t t_start = std::clock();
	double x_f = (double)newXSize / xSize;
	double y_f = (double)newYSize / ySize;
	double a, b;
	int i, j, ind;
	uchar x[12];
	uchar y[12];
	uchar* temp = new uchar[(xSize + 4) * (ySize + 4) * 3];
	extendBorders(input, xSize, ySize, temp, 2);

	for (int q = 0; q < newYSize; q++)
	{
		j = q / y_f;
		for (int p = 0; p < newXSize; p++)
		{
			a = p / x_f - floor(p / x_f);
			b = q / y_f - floor(q / y_f);
			i = p / x_f;
			ind = 0;
			for (int k = j - 1; k < j + 3; k++)
			{
				for (int l = -1; l < 3; l++)
				{
					x[l + 1] = temp[(k + 2) * (xSize + 4) * 3 + (i + l + 2) * 3];
					x[l + 5] = temp[(k + 2) * (xSize + 4) * 3 + (i + l + 2) * 3 + 1];
					x[l + 9] = temp[(k + 2) * (xSize + 4) * 3 + (i + l + 2) * 3 + 2];
				}
				y[ind]       = cubicInterpolate(x, a);
				y[ind + 4]   = cubicInterpolate(&(x[4]), a);
				y[ind++ + 8] = cubicInterpolate(&(x[8]), a);
			}
			output[q * newXSize * 3 + p * 3]     = cubicInterpolate(y, b);
			output[q * newXSize * 3 + p * 3 + 1] = cubicInterpolate(&(y[4]), b);
			output[q * newXSize * 3 + p * 3 + 2] = cubicInterpolate(&(y[8]), b);
		}
	}
	std::clock_t t_end = std::clock();
	qDebug("Execution time for bicubic interpolation: %.5f ms", (1000.0 * (t_end - t_start)) / CLOCKS_PER_SEC);
	delete[] temp;
}

void imageTransform(const uchar input[], int xSize, int ySize, uchar output[], double k1, double k2)
{
	std::clock_t t_start = std::clock();
	int m = xSize / 2;
	int n = ySize / 2;
	double i_norm, j_norm;
	int i_prim, j_prim;
	double r;

	for (int j = 0; j < ySize; j++)
	{
		j_norm = (double)(j - n) / ySize;
		for (int i = 0; i < xSize; i++)
		{
			i_norm = (double)(i - m) / xSize;
			r = i_norm * i_norm + j_norm * j_norm;
			i_prim = floor(m + (i - m) * (1 + k1 * r + k2 * r * r) + 0.5);
			j_prim = floor(n + (j - n) * (1 + k1 * r + k2 * r * r) + 0.5);

			if (i_prim >= 0 && i_prim < xSize && j_prim >= 0 && j_prim < ySize)
			{
				output[j * xSize * 3 + i * 3]     = input[j_prim * xSize * 3 + i_prim * 3];
				output[j * xSize * 3 + i * 3 + 1] = input[j_prim * xSize * 3 + i_prim * 3 + 1];
				output[j * xSize * 3 + i * 3 + 2] = input[j_prim * xSize * 3 + i_prim * 3 + 2];
			}
			else
			{
				output[j * xSize * 3 + i * 3]     = 0;
				output[j * xSize * 3 + i * 3 + 1] = 0;
				output[j * xSize * 3 + i * 3 + 2] = 0;
			}
		}
	}
	std::clock_t t_end = std::clock();
	qDebug("Execution time for fisheye: %.5f ms", (1000.0 * (t_end - t_start)) / CLOCKS_PER_SEC);
}

void imageTransformBilinear(const uchar input[], int xSize, int ySize, uchar output[], double k1, double k2)
{
	std::clock_t t_start = std::clock();
	int m = xSize / 2;
	int n = ySize / 2;
	double i_norm, j_norm, r, a, b;
	int i_prim, j_prim;

	for (int j = 0; j < ySize; j++)
	{
		j_norm = (double)(j - n) / ySize;
		for (int i = 0; i < xSize; i++)
		{
			i_norm = (double)(i - m) / xSize;
			r = i_norm * i_norm + j_norm * j_norm;
			a = (i - m) * (1 + k1 * r + k2 * r * r) - floor((i - m) * (1 + k1 * r + k2 * r * r));
			b = (j - n) * (1 + k1 * r + k2 * r * r) - floor((j - n) * (1 + k1 * r + k2 * r * r));
			i_prim = m + (i - m) * (1 + k1 * r + k2 * r * r);
			j_prim = n + (j - n) * (1 + k1 * r + k2 * r * r);

			if (i_prim >= 0 && i_prim < (xSize - 1) && j_prim >= 0 && j_prim < (ySize - 1))
			{
				output[j * xSize * 3 + i * 3]     = (1 - a) * (1 - b) * input[j_prim * xSize * 3 + i_prim * 3]     + (1 - a) * b * input[(j_prim + 1) * xSize * 3 + i_prim * 3]     + a * (1 - b) * input[j_prim * xSize * 3 + (i_prim + 1) * 3]     + a * b * input[(j_prim + 1) * xSize * 3 + (i_prim + 1) * 3];
				output[j * xSize * 3 + i * 3 + 1] = (1 - a) * (1 - b) * input[j_prim * xSize * 3 + i_prim * 3 + 1] + (1 - a) * b * input[(j_prim + 1) * xSize * 3 + i_prim * 3 + 1] + a * (1 - b) * input[j_prim * xSize * 3 + (i_prim + 1) * 3 + 1] + a * b * input[(j_prim + 1) * xSize * 3 + (i_prim + 1) * 3 + 1];
				output[j * xSize * 3 + i * 3 + 2] = (1 - a) * (1 - b) * input[j_prim * xSize * 3 + i_prim * 3 + 2] + (1 - a) * b * input[(j_prim + 1) * xSize * 3 + i_prim * 3 + 2] + a * (1 - b) * input[j_prim * xSize * 3 + (i_prim + 1) * 3 + 2] + a * b * input[(j_prim + 1) * xSize * 3 + (i_prim + 1) * 3 + 2];
			}
			else
			{
				output[j * xSize * 3 + i * 3] = 0;
				output[j * xSize * 3 + i * 3 + 1] = 0;
				output[j * xSize * 3 + i * 3 + 2] = 0;
			}
		}
	}
	std::clock_t t_end = std::clock();
	qDebug("Execution time for fisheye bilinear: %.5f ms", (1000.0 * (t_end - t_start)) / CLOCKS_PER_SEC);
}

void imageInverseTransform(const uchar input[], int xSize, int ySize, uchar output[], double k1, double k2)
{
	std::clock_t t_start = std::clock();
	int m = xSize / 2;
	int n = ySize / 2;
	double i_norm, j_norm;
	int i_prim, j_prim;
	double r;

	memset(output, 0, xSize * ySize * 3);

	for (int j = 0; j < ySize; j++)
	{
		j_norm = (double)(j - n) / ySize;
		for (int i = 0; i < xSize; i++)
		{
			i_norm = (double)(i - m) / xSize;
			r = i_norm * i_norm + j_norm * j_norm;
			i_prim = floor(m + (i - m) * (1 + k1 * r + k2 * r * r) + 0.5);
			j_prim = floor(n + (j - n) * (1 + k1 * r + k2 * r * r) + 0.5);

			if (i_prim >= 0 && i_prim < xSize && j_prim >= 0 && j_prim < ySize)
			{
				output[j_prim * xSize * 3 + i_prim * 3]     = input[j * xSize * 3 + i * 3];
				output[j_prim * xSize * 3 + i_prim * 3 + 1] = input[j * xSize * 3 + i * 3 + 1];
				output[j_prim * xSize * 3 + i_prim * 3 + 2] = input[j * xSize * 3 + i * 3 + 2];
			}
		}
	}
	std::clock_t t_end = std::clock();
	qDebug("Execution time for inverse fisheye: %.5f ms", (1000.0 * (t_end - t_start)) / CLOCKS_PER_SEC);
}

void imageInverseTransformBilinear(const uchar input[], int xSize, int ySize, uchar output[], double k1, double k2)
{
	std::clock_t t_start = std::clock();
	int m = xSize / 2;
	int n = ySize / 2;
	double i_norm, j_norm, r, a, b;
	int i_prim, j_prim;

	memset(output, 0, xSize * ySize * 3);
	
	for (int j = 0; j < ySize; j++)
	{
		j_norm = (double)(j - n) / ySize;
		for (int i = 0; i < xSize; i++)
		{
			i_norm = (double)(i - m) / xSize;
			r = i_norm * i_norm + j_norm * j_norm;
			a = (i - m) * (1 + k1 * r + k2 * r * r) - floor((i - m) * (1 + k1 * r + k2 * r * r));
			b = (j - n) * (1 + k1 * r + k2 * r * r) - floor((j - n) * (1 + k1 * r + k2 * r * r));
			i_prim = m + (i - m) * (1 + k1 * r + k2 * r * r);
			j_prim = n + (j - n) * (1 + k1 * r + k2 * r * r);

			if (i_prim >= 0 && i_prim < xSize && j_prim >= 0 && j_prim < ySize)
			{
				output[j_prim * xSize * 3 + i_prim * 3]     = (1 - a) * (1 - b) * input[j * xSize * 3 + i * 3]     + (1 - a) * b * input[(j + 1) * xSize * 3 + i * 3]     + a * (1 - b) * input[j * xSize * 3 + (i + 1) * 3]     + a * b * input[(j + 1) * xSize * 3 + (i + 1) * 3];
				output[j_prim * xSize * 3 + i_prim * 3 + 1] = (1 - a) * (1 - b) * input[j * xSize * 3 + i * 3 + 1] + (1 - a) * b * input[(j + 1) * xSize * 3 + i * 3 + 1] + a * (1 - b) * input[j * xSize * 3 + (i + 1) * 3 + 1] + a * b * input[(j + 1) * xSize * 3 + (i + 1) * 3 + 1];
				output[j_prim * xSize * 3 + i_prim * 3 + 2] = (1 - a) * (1 - b) * input[j * xSize * 3 + i * 3 + 2] + (1 - a) * b * input[(j + 1) * xSize * 3 + i * 3 + 2] + a * (1 - b) * input[j * xSize * 3 + (i + 1) * 3 + 2] + a * b * input[(j + 1) * xSize * 3 + (i + 1) * 3 + 2];
			}
		}
	}
	std::clock_t t_end = std::clock();
	qDebug("Execution time for inverse fisheye bilinear: %.5f ms", (1000.0 * (t_end - t_start)) / CLOCKS_PER_SEC);
}