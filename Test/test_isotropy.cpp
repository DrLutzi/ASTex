#include <stdlib.h>
#include <ASTex/Perioblend/perioblend.h>
#include <ASTex/easy_io.h>
#include <ASTex/rpn_utils.h>

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
	texturePool.addTexture(im_in);
	texturePool.addTransformationMatrix(Eigen::Rotation2D<double>(M_PI).matrix());
	texturePool.generate();
	ImageGrayd rIn, gIn, bIn;
	ImageGrayd r, g, b;
	ImageGrayd autocovariance1, autocovariance2;
	extract3Channels(im_in, rIn, gIn, bIn);
	extract3Channels(texturePool[1].texture, r, g, b);
	Fourier::truePeriodicStationaryAutocovariance(rIn, autocovariance1, true);
	Fourier::truePeriodicStationaryAutocovariance(r, autocovariance1, true);
	IO::save01_in_u8(autocovariance1, "/home/nlutz/autocovariance1.png");
	IO::save01_in_u8(autocovariance2, "/home/nlutz/autocovariance2.png");
	IO::save01_in_u8(r, "/home/nlutz/r.png");

	Fourier::truePeriodicStationaryAutocovariance(gIn, autocovariance1, true);
	Fourier::truePeriodicStationaryAutocovariance(g, autocovariance1, true);


	Fourier::truePeriodicStationaryAutocovariance(bIn, autocovariance1, true);
	Fourier::truePeriodicStationaryAutocovariance(b, autocovariance1, true);


	return 0;
}
