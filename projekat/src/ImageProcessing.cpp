
#include "ImageProcessing.h"
#include "ImageInterpolation.h"

#include <QDebug>

void imageProcessingFun(const QString& progName, QImage* const outImgs, const QImage* const inImgs, const QVector<double>& params) 
{
	int X_SIZE = inImgs->width();
	int Y_SIZE = inImgs->height();

	/* NOTE: Calculate output image resolution and construct output image object */

	if(progName == "Sample and hold") 
	{	
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		/* Vertical scale factor is params[0] */
		/* Horizontal scale factor is params[1] */
		int newXSize = X_SIZE * params[1];
		int newYSize = Y_SIZE * params[0];

		/* TO DO: Calculate output image resolution and construct output image object */
		int x_moduo;
		x_moduo = newXSize % 4;
		if (x_moduo != 0)
		{
			newXSize = newXSize + 4 - x_moduo;
		}

		/* Create empty output image */
		*outImgs = *(new QImage(newXSize, newYSize, inImgs->format()));

		/* TO DO: Perform Sample and hold interpolation  */
		sampleAndHold(inImgs->bits(), X_SIZE, Y_SIZE, outImgs->bits(), newXSize, newYSize);

	}
	else if (progName == "Bilinear") 
	{
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		/* Vertical scale factor is params[0] */
		/* Horizontal scale factor is params[1] */
		int newXSize = X_SIZE * params[1];
		int newYSize = Y_SIZE * params[0];

		/* TO DO: Calculate output image resolution and construct output image object */
		int x_moduo;
		x_moduo = newXSize % 4;
		if (x_moduo != 0)
		{
			newXSize = newXSize + 4 - x_moduo;
		}

		/* Create empty output image */
		*outImgs = *(new QImage(newXSize, newYSize, inImgs->format()));

		/* TO DO: Perform Bilinear interpolation  */
		bilinearInterpolate(inImgs->bits(), X_SIZE, Y_SIZE, outImgs->bits(), newXSize, newYSize);
	}
	else if (progName == "Bicubic")
	{
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		/* Vertical scale factor is params[0] */
		/* Horizontal scale factor is params[1] */
		int newXSize = X_SIZE * params[1];
		int newYSize = Y_SIZE * params[0];

		/* TO DO: Calculate output image resolution and construct output image object */
		int x_moduo;
		x_moduo = newXSize % 4;
		if (x_moduo != 0)
		{
			newXSize = newXSize + 4 - x_moduo;
		}

		/* Create empty output image */
		*outImgs = *(new QImage(newXSize, newYSize, inImgs->format()));

		/* TO DO: Perform Bilinear interpolation  */
		bicubicInterpolate(inImgs->bits(), X_SIZE, Y_SIZE, outImgs->bits(), newXSize, newYSize);
	}
	else if(progName == "Transform") 
	{	
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		/* k1 and k2 parameters are given as params[0] and params[1]*/
		double k1 = params[0];
		double k2 = params[1];

		/* TO DO: Construct output image object */
		*outImgs = *(new QImage(X_SIZE, Y_SIZE, inImgs->format()));

		/* TO DO: Perform image transformation */
		imageTransform(inImgs->bits(), X_SIZE, Y_SIZE, outImgs->bits(), k1, k2);
	}
	else if (progName == "Transform Bilinear") 
	{
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		/* k1 and k2 parameters are given as params[0] and params[1]*/
		double k1 = params[0];
		double k2 = params[1];

		/* TO DO: Construct output image object */
		*outImgs = *(new QImage(X_SIZE, Y_SIZE, inImgs->format()));

		/* TO DO: Perform image transformation with bilinear interpolation */
		imageTransformBilinear(inImgs->bits(), X_SIZE, Y_SIZE, outImgs->bits(), k1, k2);
	}
	else if (progName == "Inverse Transform")
	{
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		/* k1 and k2 parameters are given as params[0] and params[1]*/
		double k1 = params[0];
		double k2 = params[1];

		/* TO DO: Construct output image object */
		*outImgs = *(new QImage(X_SIZE, Y_SIZE, inImgs->format()));

		/* TO DO: Perform image transformation */
		imageInverseTransform(inImgs->bits(), X_SIZE, Y_SIZE, outImgs->bits(), k1, k2);
	}
	else if (progName == "Inverse Transform Bilinear")
	{
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		/* k1 and k2 parameters are given as params[0] and params[1]*/
		double k1 = params[0];
		double k2 = params[1];

		/* TO DO: Construct output image object */
		*outImgs = *(new QImage(X_SIZE, Y_SIZE, inImgs->format()));

		/* TO DO: Perform image transformation with bilinear interpolation */
		imageInverseTransformBilinear(inImgs->bits(), X_SIZE, Y_SIZE, outImgs->bits(), k1, k2);
	}
}

