#include <stdlib.h>
#include <vector>
#include <algorithm>
#include "ASTex/fourier.h"
#include "ASTex/utils.h"
#include "ASTex/easy_io.h"

using namespace ASTex;

typedef std::pair<itk::Index<2>, ImageGrayd::PixelType> FullTexel;

class CompareTexelValue
{
public:

	CompareTexelValue(){}

	typedef FullTexel value_type;
	bool operator()(FullTexel &tex1, FullTexel &tex2)
	{
		return tex1.second < tex2.second;
	};
};

int main(int argc, char **argv)
{
	if(argc < 2)
	{
		std::cerr << "Usage:" << argv[0] << "<texture file.png>" << std::endl;
		return EXIT_FAILURE;
	}

	ImageGrayd img, acorr;
	ASTex::IO::loadu8_in_01(img, argv[1]);

	std::string name_file = IO::remove_path(argv[1]);
	std::string name_noext = IO::remove_ext(name_file);

	std::cout << Fourier::stationaryMean(img) << std::endl;
	std::cout << Fourier::stationaryVariance(img) << std::endl;

	acorr.initItk(img.width(), img.height());
	bool exists;
	exists = Fourier::truePeriodicStationaryAutocovariance(img, acorr, true);

	if(exists)
	{
		acorr.for_all_pixels([&] (ImageGrayd::PixelType &pix)
		{
			pix = pix < 0 ? 0 : pix;
		});
		IO::save01_in_u8(acorr, std::string("/home/nlutz/acorr_") + name_noext + ".png");

//		std::vector<FullTexel> texelVector;
//		texelVector.reserve(acorr.width()*acorr.height());
//		acorr.for_all_pixels([&] (ImageGrayd::PixelType &pix, int x, int y)
//		{
//			FullTexel tex;
//			tex.first = itk::Index<2>();
//			tex.first[0] = x;
//			tex.first[1] = y;
//			tex.second = pix;
//		});
//		CompareTexelValue compare;
//		std::sort(texelVector.begin(), texelVector.end(), compare);
		//std::rand()%texelVector.size()/6;
	}

	return 0;
}
