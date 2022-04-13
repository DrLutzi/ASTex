#ifndef PERIOBLEND_PERIOBLEND_H_
#define PERIOBLEND_PERIOBLEND_H_

#include <Eigen/Eigen>
#include "Algo/ProceduralNoiseFiltering/gaussian_transfer.h"
#include "ASTex/utils.h"
#include "ASTex/pca.h"
#include <cmath>
#include "ASTex/colorspace_filters.h"
#include "ASTex/easy_io.h"
#include "ASTex/histogram.h"
#include "ASTex/texture_pool.h"
#include "ASTex/CSN/csn_texture.h"
#include "ASTex/Perioblend/affineTransformations.h"

namespace ASTex
{

/**
 * @brief The TandBFunctionType class -- see Perioblend::addFunction for all the infos
 */
class TandBFunctionType
{
public:
	TandBFunctionType() {}

	using BlendingFunctionType = std::function<double (double u, double v)>;
	using TilingFunctionType = std::function<Eigen::Vector2i (double u, double v)>;

	TilingFunctionType tilingFunction;
	BlendingFunctionType blendingFunction;

	ImageRGBd visualizeTiling(unsigned int width, unsigned int height, double uMax, double vMax) const;
	ImageGrayd visualizeBlending(unsigned int width, unsigned int height, double uMax, double vMax) const;
};

ImageRGBd TandBFunctionType::visualizeTiling(unsigned int width, unsigned int height, double uMax, double vMax) const
{
	ImageRGBd visualization;
	visualization.initItk(width, height, true);
	Eigen::Vector2i min, max;
	min = tilingFunction(0, 0);
	max = tilingFunction(0, 0);
	visualization.for_all_pixels([&] (ImageRGBd::PixelType &pix, int x, int y)
	{
		(void) pix;
		double u, v;
		u = double(x)/(visualization.width()-1);
		u = u*uMax;
		v = double(y)/(visualization.height()-1);
		v = v*vMax;
		Eigen::Vector2i value = tilingFunction(u, v);
		min[0] = value[0] < min[0] ? value[0] : min[0];
		min[1] = value[1] < min[1] ? value[1] : min[1];
		max[0] = value[0] > max[0] ? value[0] : max[0];
		max[1] = value[1] > max[1] ? value[1] : max[1];
	});
	Eigen::Vector2d minDouble = min.cast<double>(), maxDouble = max.cast<double>();
	std::cout << minDouble[0] << ", " << minDouble[1] << std::endl;
	std::cout << maxDouble[0] << ", " << maxDouble[1] << std::endl;
	visualization.for_all_pixels([&] (ImageRGBd::PixelType &pix, int x, int y)
	{
		double u, v;
		u = double(x)/(visualization.width()-1);
		u = u*uMax;
		v = double(y)/(visualization.height()-1);
		v = v*vMax;
		Eigen::Vector2i value = tilingFunction(u, v);
		Eigen::Vector2d valueDouble = value.cast<double>();
		valueDouble[0] = (valueDouble[0] - minDouble[0])/(maxDouble[0] - minDouble[0]);
		valueDouble[1] = (valueDouble[1] - minDouble[1])/(maxDouble[1] - minDouble[1]);
		pix[0] = valueDouble[0];
		pix[1] = valueDouble[1];
	});
	return visualization;
}

ImageGrayd TandBFunctionType::visualizeBlending(unsigned int width, unsigned int height, double uMax, double vMax) const
{
	ImageGrayd visualization;
	visualization.initItk(width, height, true);
	double min=0;
	double max=0;
	visualization.for_all_pixels([&] (ImageGrayd::PixelType &pix, int x, int y)
	{
		(void) pix;
		double u, v;
		u = double(x)/(visualization.width()-1);
		u = u*uMax;
		v = double(y)/(visualization.height()-1);
		v = v*vMax;
		double value = blendingFunction(u, v);
		min = value < min ? value : min;
		max = value > max ? value : max;
	});
	visualization.for_all_pixels([&] (ImageGrayd::PixelType &pix, int x, int y)
	{
		double u, v;
		u = double(x)/(visualization.width()-1);
		u = u*uMax;
		v = double(y)/(visualization.height()-1);
		v = v*vMax;
		double value = blendingFunction(u, v);
		pix = (value - min)/(max - min);
	});

	return visualization;
}


template<typename I>
/**
 * @brief The TilingAndBlending class represents a generic tiling and blending process for texture synthesis.
 */
class TilingAndBlending
{
public:

	typedef enum{STATIONARY, CYCLOSTATIONARY, PERIODIC} ProcessType;

	using ImageType					= I;
	using ImageNoiseType			= ImageGrayd;
	using PixelType					= typename ImageType::PixelType;
	using DataType					= typename ImageType::DataType;
	using PixelPosType				= itk::Index<2>;
	using PcaImageType				= ImageRGB<DataType>;
	using PcaPixelType				= typename PcaImageType::PixelType;
	using PcaType					= PCA<DataType>;
	using GaussianTransferType		= Gaussian_transfer<PcaImageType>;
	using LutType					= PcaImageType;
	using ProceduralBlendingType	= std::function<PcaImageType(const PcaImageType &)>;
	using PtrImageType				= ImageGrayu64;

	TilingAndBlending();

	/**
	 * @brief synthesize_stationary uses the first texture added to synthesize a texture using a by-example stationary ergodic synthesis.
	 * You need to have added at least two TandBFunctionType using addFunction for it to produce good results.
	 * @return the synthesized texture.
	 */
	ImageType synthesize_stationary();

	/**
	 * @brief synthesize_periodic uses all textures added to synthesize a texture using a by-examples periodic non ergodic synthesis.
	 * You need to have added at least two TandBFunctionType using addFunction for it to produce good results.
	 * @return the synthesized texture.
	 */
	ImageType synthesize_periodic();

	/**
	 * @brief addTexture adds a texture to the pool. Only useful to be called more than once for periodic, non ergodic synthesis.
	 * @param image examplar texture.
	 */
	void addTexture(const ImageType &image);

	/**
	 * @brief setProcessType unused yet
	 */
	void setProcessType(ProcessType);

	/**
	 * @brief addFunction add a user-defind blending + tiling function.
	 * The blending function must reach 0 where the tiling function changes its returned index, otherwise the results will have seams.
	 * You can add as many as you want, although this class was re-purposed to usually have only two or three.
	 * Different functions are meant to overlap.
	 * @param tabft user-defind blending + tiling function
	 */
	void addFunction(const TandBFunctionType &tabft);

	/**
	 * @brief setSingularityBlendingFunction sets a blending function to hide singularities.
	 * @param sbf user-defined blending function.
	 */
	void setSingularityBlendingFunction(const TandBFunctionType::BlendingFunctionType &sbf);
	void setWidth(unsigned int width);
	void setHeight(unsigned int height);
	void setUVScale(double uScale, double vScale);

	/**
	 * @brief setUseHistogramTransfer if the histogram transfer of [HN18] is used.
	 * @param useTransfer boolean
	 */
	void setUseHistogramTransfer(bool useTransfer);

	/**
	 * @brief setRandomAffineTransform the random affine transform generator used to generate random rotations in each tile.
	 * @param randomAffineTransform -- see the related class.
	 */
	void setRandomAffineTransform(const RandomAffineTransform &randomAffineTransform);

	ImageGrayd visualizeBlendingSum(unsigned int width, unsigned int height, double uMax, double vMax);
	ImageRGBd visualizeSynthesis(unsigned int width, unsigned int height);

private:

	Eigen::Vector2d hash(const Eigen::Vector2d &p) const;
	/**
	 * @brief cantorPairingFn attributes one int from two integers, used for seeding.
	 * @param v vector of two integers
	 * @return one integer
	 */
	int cantorPairingFunction(const Eigen::Vector2i &v) const;

	TexturePool<ImageType>	m_texturePool;
	std::list<TandBFunctionType> m_functions;
	TandBFunctionType::BlendingFunctionType	m_singularityBlendingFunction;
	ProcessType				m_processType;
	unsigned int			m_width;
	unsigned int			m_height;
	double					m_uScale;
	double					m_vScale;
	bool					m_useTransfer;
	RandomAffineTransform	m_randomAffineTransform;
};

template<typename I>
TilingAndBlending<I>::TilingAndBlending() :
	m_texturePool(),
	m_functions(),
	m_singularityBlendingFunction(nullptr),
	m_processType(STATIONARY),
	m_width(1024),
	m_height(1024),
	m_uScale(1.0),
	m_vScale(1.0),
	m_useTransfer(true),
	m_randomAffineTransform()
{}

template<typename I>
typename TilingAndBlending<I>::ImageType TilingAndBlending<I>::synthesize_stationary()
{
	m_texturePool.generate();
	const ImageType &texture = m_texturePool[0].texture;
	PcaImageType pcaTexture = CSN::CSN_Texture<ImageType>::toPcaImageType(texture);

	//Forward PCA
	PcaImageType outputPCA;
	outputPCA.initItk(m_width, m_height);
	PcaType pca(pcaTexture);
	MaskBool mb_alwaysTrue(pcaTexture.width(), pcaTexture.height());
	mb_alwaysTrue |= [] (int, int) {return true;};
	pca.computePCA(mb_alwaysTrue);
	pca.project(pcaTexture);

	//Forward Histogram transfer
	PcaImageType pcaGaussianTexture;
	LutType lut;
	lut.initItk(128, 1);
	GaussianTransferType gtt;
	pcaGaussianTexture.initItk(pcaTexture.width(), pcaTexture.height());
	PcaPixelType mean = PcaImageType::zero();
	if(m_useTransfer)
	{
		gtt.ComputeTinput(pcaTexture, pcaGaussianTexture);
		gtt.ComputeinvT(pcaTexture, lut);
		for(int i = 0; i<3; ++i)
		{
			mean[i] = 0.5;
		}
	}
	else
	{
		pcaGaussianTexture = pcaTexture;
		pcaTexture.for_all_pixels([&] (PcaPixelType &pix)
		{
			for(int i=0; i<3; ++i)
			{
				mean[i] += pix[i];
			}
		});
		for(int i=0; i<3; ++i)
		{
			mean[i]/=pcaTexture.width()*pcaTexture.height();
		}
	}

	//Computation
	outputPCA.for_all_pixels([&] (PcaPixelType &pix, int x, int y)
	{
		double u = double(x)/(pcaTexture.width()-1);
		double v = double(y)/(pcaTexture.height()-1);
		double sum = 0;
		double sumSquare = 0;
		for(auto const &fonctions : m_functions)
		{
			double blendingValue = fonctions.blendingFunction(u * m_uScale, v * m_vScale);
			sum += blendingValue;
		}
		if(m_singularityBlendingFunction != nullptr)
		{
			double blendingValue = m_singularityBlendingFunction(u * m_uScale, v * m_vScale);
			sum += blendingValue;
		}
		for(auto const &fonctions : m_functions)
		{
			Eigen::Vector2i tile = fonctions.tilingFunction(u * m_uScale, v * m_vScale);
			Eigen::Vector2d uvTexture = hash(tile.cast<double>());
			uvTexture[0] += u;
			uvTexture[1] += v;
			AffineTransform affT = m_randomAffineTransform.generate(true, cantorPairingFunction(tile));
			Eigen::Vector2d uvTextureFract = CSN::CSN_Texture<I>::fract(affT.transform()*uvTexture);
			PixelPosType xyTexture;
			xyTexture[0] = int(uvTextureFract[0] * pcaTexture.width()) % pcaTexture.width();
			xyTexture[1] = int(uvTextureFract[1] * pcaTexture.height()) % pcaTexture.height();
			double blendingValue = fonctions.blendingFunction(u * m_uScale, v * m_vScale) / sum;
			pix += pcaGaussianTexture.pixelAbsolute(xyTexture) * blendingValue;
			sumSquare += blendingValue*blendingValue;
		}
		if(m_singularityBlendingFunction != nullptr)
		{
			double blendingValue = m_singularityBlendingFunction(u * m_uScale, v * m_vScale) / sum;
			sumSquare += blendingValue*blendingValue;
		}
		sumSquare = std::sqrt(sumSquare);
		for(unsigned int i=0; i<3; ++i)
		{
			pix[i] -= mean[i];
			pix[i] /= sumSquare;
			pix[i] += mean[i];
		}
	});

	//Inverse Histogram transfer
	if(m_useTransfer)
	{
		outputPCA = GaussianTransferType::invT(outputPCA, lut);
	}

	//Inverse PCA
	pca.back_project(outputPCA, pcaTexture);

	ImageType output = CSN::CSN_Texture<ImageType>::fromPcaImageType(pcaTexture);
	return output;
}

template<typename I>
typename TilingAndBlending<I>::ImageType TilingAndBlending<I>::synthesize_periodic()
{
	m_texturePool.generate();
	TexturePool<PcaImageType> texturePool;
	for(unsigned int i=0; i<m_texturePool.size(); ++i)
	{
		assert(	m_texturePool[i].texture.width() == m_texturePool[0].texture.width() &&
				m_texturePool[i].texture.height() == m_texturePool[0].texture.height());
		PcaImageType pcaTexture = CSN::CSN_Texture<ImageType>::toPcaImageType(m_texturePool[i].texture);
		texturePool.addTexture(pcaTexture);
	}
	texturePool.generate();

	PcaImageType outputPCA;
	outputPCA.initItk(m_width, m_height);
	unsigned width = m_texturePool[0].texture.width();
	unsigned height = m_texturePool[0].texture.height();

	PcaImageType massivePCATexture;
	massivePCATexture.initItk(width * m_texturePool.size(), height);
	massivePCATexture.for_all_pixels([&] (PcaPixelType &pix, int x, int y)
	{
		int index = x/width;
		int xMod = x%width;
		pix = texturePool[index].texture.pixelAbsolute(xMod, y);
	});
	PcaType pca(massivePCATexture);
	MaskBool mb_alwaysTrue(massivePCATexture.width(), massivePCATexture.height());
	mb_alwaysTrue |= [] (int, int) {return true;};
	pca.computePCA(mb_alwaysTrue);
	for(unsigned int i=0; i<m_texturePool.size(); ++i)
	{
		PcaImageType projection;
		pca.project(texturePool[i].texture, projection);
		texturePool[i].texture = projection;
		IO::save01_in_u8(projection, std::string("/home/nlutz/im") + std::to_string(i) + ".png");
	}

	outputPCA.for_all_pixels([&] (PcaPixelType &pix, int x, int y)
	{
		double u = double(x)/(width-1);
		double v = double(y)/(height-1);
		Eigen::Vector2d uvTexture(u, v);
		double sum = 0;
		double sumSquare = 0;
		for(auto const &fonctions : m_functions)
		{
			double blendingValue = fonctions.blendingFunction(u * m_uScale, v * m_vScale);
			sum += blendingValue;
		}
		for(auto const &fonctions : m_functions)
		{
			Eigen::Vector2i transition = fonctions.tilingFunction(u * m_uScale, v * m_vScale);
			int indexTexture1D = int(transition[0]*41 + transition[1]*97)%texturePool.size();
			Eigen::Vector2d uvTextureFract = CSN::CSN_Texture<I>::fract(uvTexture);
			PixelPosType xyTexture;
			xyTexture[0] = int(uvTextureFract[0] * width) % width;
			xyTexture[1] = int(uvTextureFract[1] * height) % height;
			double blendingValue = fonctions.blendingFunction(u * m_uScale, v * m_vScale) / sum;
			const PcaImageType &texture = texturePool[indexTexture1D].texture;
			pix += texture.pixelAbsolute(xyTexture) * blendingValue;
			sumSquare += blendingValue*blendingValue;
		}
		sumSquare = std::sqrt(sumSquare);
		for(unsigned int i=0; i<3; ++i)
		{
//			pix[i] -= 0.5;
			pix[i] /= sumSquare;
//			pix[i] += 0.5;
		}
	});

	//Inverse PCA
	PcaImageType pcaTexture;
	pcaTexture.initItk(outputPCA.width(), outputPCA.height());
	pca.back_project(outputPCA, pcaTexture);

	ImageType output = CSN::CSN_Texture<ImageType>::fromPcaImageType(pcaTexture);
	return output;
}

template<typename I>
void TilingAndBlending<I>::addTexture(const ImageType &texture)
{
	m_texturePool.addTexture(texture);
}

template<typename I>
void TilingAndBlending<I>::setProcessType(ProcessType processType)
{
	m_processType = processType;
}

template<typename I>
void TilingAndBlending<I>::addFunction(const TandBFunctionType &tabft)
{
	m_functions.push_back(tabft);
}

template<typename I>
void TilingAndBlending<I>::setSingularityBlendingFunction(const TandBFunctionType::BlendingFunctionType &sbf)
{
	m_singularityBlendingFunction = sbf;
}

template<typename I>
void TilingAndBlending<I>::setWidth(unsigned int width)
{
	m_width	= width;
}

template<typename I>
void TilingAndBlending<I>::setHeight(unsigned int height)
{
	m_height = height;
}

template<typename I>
void TilingAndBlending<I>::setUVScale(double uScale, double vScale)
{
	m_uScale = uScale;
	m_vScale = vScale;
}

template<typename I>
void TilingAndBlending<I>::setUseHistogramTransfer(bool useTransfer)
{
	m_useTransfer = useTransfer;
}

template<typename I>
void TilingAndBlending<I>::setRandomAffineTransform(const RandomAffineTransform &randomAffineTransform)
{
	m_randomAffineTransform = randomAffineTransform;
}

template<typename I>
Eigen::Vector2d TilingAndBlending<I>::hash(const Eigen::Vector2d &p) const
{
	Eigen::Matrix2d hashMat;
	hashMat << 127.1, 269.5, 311.7, 183.3;
	Eigen::Vector2d q = hashMat * p;
	q[0] = sin(q[0]);
	q[1] = sin(q[1]);
	return CSN::CSN_Texture<I>::fract ( q * 43758.5453 );
}

template<typename I>
int TilingAndBlending<I>::cantorPairingFunction(const Eigen::Vector2i &v) const
{
	return (v[0] + v[1]) * (v[0] + v[1] + 1) / 2 + v[0];
}

template<typename I>
ImageGrayd TilingAndBlending<I>::visualizeBlendingSum(unsigned int width, unsigned int height, double uMax, double vMax)
{
	ImageGrayd visualization;
	visualization.initItk(width, height, true);
	double min=0;
	double max=0;
	visualization.for_all_pixels([&] (ImageGrayd::PixelType &pix, int x, int y)
	{
		(void) pix;
		double u, v;
		u = double(x)/(visualization.width()-1);
		u = u*uMax;
		v = double(y)/(visualization.height()-1);
		v = v*vMax;
		double value = 0;
		for(auto const &fonctions : m_functions)
		{
			value += fonctions.blendingFunction(u, v);
		}
		min = value < min ? value : min;
		max = value > max ? value : max;
	});
	visualization.for_all_pixels([&] (ImageGrayd::PixelType &pix, int x, int y)
	{
		double u, v;
		u = double(x)/(visualization.width()-1);
		u = u*uMax;
		v = double(y)/(visualization.height()-1);
		v = v*vMax;
		double value = 0;
		for(auto const &fonctions : m_functions)
		{
			value += fonctions.blendingFunction(u, v);
		}
		pix = (value - min)/(max - min);
	});

	return visualization;
}

template<typename I>
ImageRGBd TilingAndBlending<I>::visualizeSynthesis(unsigned int width, unsigned int height)
{
	ImageRGBd visualization;
	visualization.initItk(width, height, true);
	visualization.for_all_pixels([&] (ImageRGBd::PixelType &pix, int x, int y)
	{
		double u = double(x)/(visualization.width()-1);
		double v = double(y)/(visualization.height()-1);
		double sum = 0;
		double sumSquare = 0;
		for(auto const &fonctions : m_functions)
		{
			double blendingValue = fonctions.blendingFunction(u * m_uScale, v * m_vScale);
			sum += blendingValue;
		}
		for(auto const &fonctions : m_functions)
		{
			Eigen::Vector2i transition = fonctions.tilingFunction(u * m_uScale, v * m_vScale);
			Eigen::Vector2d uvTexture = hash(transition.cast<double>());
			uvTexture[0] += u * m_uScale;
			uvTexture[1] += v * m_vScale;
			Eigen::Vector2d uvTextureFract = CSN::CSN_Texture<I>::fract(uvTexture);
			double blendingValue = fonctions.blendingFunction(u * m_uScale, v * m_vScale) / sum;
			pix[0] += uvTextureFract[0] * blendingValue;
			pix[1] += uvTextureFract[1] * blendingValue;
		}
	});
	return visualization;
}

}

#endif
