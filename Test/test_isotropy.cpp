#include <stdlib.h>
#include <ASTex/Perioblend/perioblend.h>
#include <ASTex/easy_io.h>
#include <ASTex/rpn_utils.h>
#include "Algo/ProceduralNoiseFiltering/color_map.h"

int main(int argc, char **argv)
{
	if(argc < 2)
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " <in_texture> [in_texture2]..." << std::endl;
		return EXIT_FAILURE;
	}
	using ImageType = ImageRGBd;
	ImageType im_in, im_out;
	IO::loadu8_in_01(im_in, argv[1]);
	TexturePool<ImageType> texturePool;
	texturePool.addTexture(im_in, true);
	Eigen::Matrix3d R(Eigen::Matrix3d::Identity());
	Eigen::Rotation2D<double> rotateMatrix(M_PI/2.0);
	Eigen::Matrix2d rotationMatrix2d(rotateMatrix.matrix());
	R.row(0) << rotationMatrix2d(0, 0), rotationMatrix2d(0, 1), 0;
	R.row(1) << rotationMatrix2d(1, 0), rotationMatrix2d(1, 1), 0;
	R.row(2) << 0, 0, 1;
	texturePool.addTransformationMatrix(R);
	texturePool.generate();
	ImageGrayd rIn, gIn, bIn;
	ImageGrayd r, g, b;
	ImageGrayd autocovariance1, autocovariance2;
	extract3Channels(texturePool[0].texture, rIn, gIn, bIn);
	extract3Channels(texturePool[1].texture, r, g, b);
	Fourier::truePeriodicStationaryAutocovariance(rIn, autocovariance1, true, false);
	Fourier::truePeriodicStationaryAutocovariance(r, autocovariance2, true, false);
	auto absLambda = [&] (ImageGrayd::PixelType &pix, int x, int y)
	{
		//pix = std::abs(pix);
		pix = pix;
	};
	autocovariance1.for_all_pixels(absLambda);
	autocovariance2.for_all_pixels(absLambda);
	auto substractAbs = [&] (ImageGrayd &image1, ImageGrayd &image2) -> ImageGrayd
	{
		ImageGrayd output;
		output.initItk(image1.width(), image1.height());
		output.for_all_pixels([&] (ImageGrayd::PixelType &pix, int x, int y)
		{
			pix = std::abs( image1.pixelAbsolute(x, y) - image2.pixelAbsolute(x, y) );
		});
		return output;
	};
	ImageGrayd diff = substractAbs(autocovariance1, autocovariance2);

	auto diffToColor = [&] (ImageGrayd::PixelType &pix) -> ImageRGBd::PixelType
	{
		ImageRGBd::PixelType colorPix;
		for(int i=0; i<3; ++i)
			colorPix[i] = 0.0;
		if(pix<0.1)
		{
			colorPix[2] = (pix*10.0);
		}
		else
		{
			if(pix<0.4)
			{
				colorPix[2] = 1.0-colorPix[2];
				colorPix[0] = ((pix - 0.1)/3.0)*10.0;
			}
			else
			{
				colorPix[0] = 1.0;
			}
		}
		return colorPix;
	};

	ImageRGBd colorDiffMap;
	colorDiffMap.initItk(diff.width(), diff.height());
	colorDiffMap.for_all_pixels([&] (ImageRGBd::PixelType &pix, int x, int y)
	{
		pix = diffToColor(diff.pixelAbsolute(x, y));
	});

	IO::save01_in_u8(autocovariance1, "/home/nlutz/autocovariance1.png");
	IO::save01_in_u8(autocovariance2, "/home/nlutz/autocovariance2.png");
	IO::save01_in_u8(colorDiffMap, "/home/nlutz/colorDiffMap.png");
	IO::save01_in_u8(r, "/home/nlutz/r.png");
	IO::save01_in_u8(rIn, "/home/nlutz/rIn.png");

//	Fourier::truePeriodicStationaryAutocovariance(gIn, autocovariance1, true);
//	Fourier::truePeriodicStationaryAutocovariance(g, autocovariance1, true);
	autocovariance1.for_all_pixels(absLambda);
	autocovariance2.for_all_pixels(absLambda);


//	Fourier::truePeriodicStationaryAutocovariance(bIn, autocovariance1, true);
//	Fourier::truePeriodicStationaryAutocovariance(b, autocovariance1, true);
	autocovariance1.for_all_pixels(absLambda);
	autocovariance2.for_all_pixels(absLambda);


	return 0;
}
