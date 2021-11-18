#include <stdlib.h>
#include <vector>
#include <algorithm>
#include "ASTex/fourier.h"
#include "ASTex/utils.h"
#include "ASTex/easy_io.h"
#include "Eigen/Eigen"
#include "ASTex/CSN/csn_texture.h"
#include "ASTex/Stamping/stamper.h"

using ImageType = ImageRGBd;
using PcaImageType = CSN::CSN_Texture<ImageType>::PcaImageType;

ImageGrayd acorr;

PcaImageType spotNoiseProcedure(const PcaImageType &image)
{
	float spotNoiseOutputScale = 1.0;
	PcaImageType centeredImage;
	centeredImage.initItk(image.width(), image.height());
	centeredImage.copy_pixels(image);

	HistogramRGBBase<typename PcaImageType::DataType> histo(centeredImage);
	PcaImageType::PixelType mean = histo.meanPixelType();
	PcaImageType meanImage;
	meanImage.initItk(image.width(), image.height(), true);
	meanImage.for_all_pixels([&] (typename PcaImageType::PixelType &pix)
	{
		pix = mean;
	});

	centeredImage.for_all_pixels([&] (PcaImageType::PixelType &pix, int x, int y)
	{
		pix[0] -= meanImage.pixelAbsolute(x, y)[0];
		pix[1] -= meanImage.pixelAbsolute(x, y)[1];
		pix[2] -= meanImage.pixelAbsolute(x, y)[2];
		pix = pix * (1.0/std::sqrt(centeredImage.width()*centeredImage.height()));
	});
	double alpha = 0.0;
	//normalization
	if(alpha>0.0)
	{
		std::cout << "applying smooth transition function" << std::endl;
		PcaImageType::PixelType norm = PcaImageType::zero(), normPhi = PcaImageType::zero();
		centeredImage.for_all_pixels([&] (PcaImageType::PixelType &pix)
		{
			for(unsigned i=0; i<3; ++i)
			{
				norm[i] += pix[i] * pix[i];
			}
		});
		for(unsigned i=0; i<3; ++i)
		{
			norm[i] = sqrt(norm[i]);
		}

		centeredImage.for_all_pixels([&] (PcaImageType::PixelType &pix, int x, int y)
		{
			double distanceToBorder;
			double distanceX = double(std::min(x, centeredImage.width()-1-x));
			double distanceY = double(std::min(y, centeredImage.height()-1-y));
			distanceToBorder = std::min(distanceX/centeredImage.width(), distanceY/centeredImage.height());
			for(unsigned i=0; i<3; ++i)
			{
				pix[i] = pix[i]/norm[i];
				if(distanceToBorder<alpha)
				{
					pix[i] = pix[i] * sqrt(distanceToBorder/alpha);
				}
				normPhi[i] += pix[i]*pix[i];
			}
		});
		for(unsigned i=0; i<3; ++i)
		{
			normPhi[i] = sqrt(normPhi[i]);
		}
		centeredImage.for_all_pixels([&] (PcaImageType::PixelType &pix, int x, int y)
		{
			for(unsigned i=0; i<3; ++i)
			{
				pix[i] /= normPhi[i];
				pix[i] *= norm[i];
			}
		});
	}

	Stamping::StampDiscrete<PcaImageType> stamp(centeredImage);
	stamp.setInterpolationRule(Stamping::StampDiscrete<PcaImageType>::BILINEAR_PERIODIC);
	Stamping::SamplerImportance sampler(acorr);
	sampler.setNbPoints(100);
	Stamping::StamperTexton<PcaImageType> stamper(&sampler, &stamp);
	stamper.setPeriodicity(true);
	stamper.setSpot(false);
	centeredImage = stamper.generate(	image.width()*spotNoiseOutputScale,
										image.height()*spotNoiseOutputScale);
	centeredImage.for_all_pixels([&] (PcaImageType::PixelType &pix, int x, int y)
	{
		itk::Index<2> index;
		index[0] = x%meanImage.width();
		index[1] = y%meanImage.height();
		for(int i=0; i<3; ++i)
		{
			pix[i] += meanImage.pixelAbsolute(index)[i];
		}
	});
	return centeredImage;
};

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

	ImageType img;
	ASTex::IO::loadu8_in_01(img, argv[1]);

	std::string name_file = IO::remove_path(argv[1]);
	std::string name_noext = IO::remove_ext(name_file);

	CSN::CSN_Texture<ImageType> csn;
	csn.setTexture(img);

	ImageGrayd img_r, img_g, img_b;
	extract3Channels(img, img_r, img_g, img_b);

	acorr.initItk(img.width(), img.height());
	bool exists;
	exists = Fourier::truePeriodicStationaryAutocovariance(img_r, acorr, true);

	if(exists)
	{
		acorr.for_all_pixels([&] (ImageGrayd::PixelType &pix)
		{
			pix = pix < 0 ? 0 : pix;
		});
		IO::save01_in_u8(acorr, std::string("/home/nlutz/") + name_noext + "_acorr.png");

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

//	//The following is how I produced the figure showing the distance map
//	ImageRGBd cycleEvaluationMap = csn.debug_cycleEvaluationMap(129, 129, Eigen::Vector2d(cyclePair.vectors[0][0], cyclePair.vectors[1][1]), 0.1);
//	IO::save01_in_u8(cycleEvaluationMap, std::string("/home/nlutz/cycleEvaluationMap129_") + textureName + ".png");

	csn.setCycles(Eigen::Vector2d(), Eigen::Vector2d());
	csn.setUseCycles(false);
	csn.setGamma(1.0);
	csn.setUsePca(true);
	csn.setUseGaussianTransfer(true);
	csn.setUseYCbCr(false);
	csn.setUseCyclicTransfer(false);
	csn.setUVScale(1.0);
	csn.setProceduralBlendingSubstitute(spotNoiseProcedure);
	ImageType output = csn.synthesize(img.width(), img.height());
	IO::save01_in_u8(output, std::string("/home/nlutz/" + name_noext + "_output.png"));
//	output.for_all_pixels([&] (ImageType::PixelType &pix)
//	{
//		for(int i=0; i<3; ++i)
//		{
//			pix[i] = pix[i] > 1.0 ? 1.0 : (pix[i] < 0.0 ? 0.0 : pix[i]);
//		}
//	});

	std::cout << Fourier::stationaryMean(img_r) << std::endl;
	std::cout << Fourier::stationaryVariance(img_r) << std::endl;

	ImageGrayu8 pdf;
	pdf.initItk(acorr.width(), acorr.height(), true);
	Stamping::SamplerImportance sampler(acorr);
	sampler.setNbPoints(10000);
	std::vector<Eigen::Vector2f> points = sampler.generate();
	for(std::vector<Eigen::Vector2f>::const_iterator cit = points.begin(); cit!=points.end(); ++cit)
	{
		ImageGrayu8::PixelType &p = pdf.pixelAbsolute(int((*cit)[0]*(pdf.width()-1)), int((*cit)[1]*(pdf.height()-1)));
		p = p == 255? 255 : p+15;
	}
	pdf.save(std::string("/home/nlutz/") + name_noext + "_pdf.png");

	return 0;
}
