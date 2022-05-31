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
#include "ASTex/Stamping/sampler.h"

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
		u = double(x)/(visualization.width());
		u = u*uMax;
		v = double(y)/(visualization.height());
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
		u = double(x)/(visualization.width());
		u = u*uMax;
		v = double(y)/(visualization.height());
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
		u = double(x)/(visualization.width());
		u = u*uMax;
		v = double(y)/(visualization.height());
		v = v*vMax;
		double value = blendingFunction(u, v);
		min = value < min ? value : min;
		max = value > max ? value : max;
	});
	visualization.for_all_pixels([&] (ImageGrayd::PixelType &pix, int x, int y)
	{
		double u, v;
		u = double(x)/(visualization.width());
		u = u*uMax;
		v = double(y)/(visualization.height());
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

	using Vec2i = Eigen::Vector2i;
	using Vec2 = Eigen::Vector2d;

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
	ImageRGBd visualizePrimalAndDual();

	/**
	 * @brief cantorPairingFn attributes one int from two integers, used for seeding.
	 * @param v vector of two integers
	 * @return one integer
	 */
	static int cantorPairingFunction(const Eigen::Vector2i &v);


	//Cyclostationary DLC

	void setCSPolyphaseComponentSamplesNumber(Vec2i samplesNumber);
	void setCSPolyphaseComponentSamplesNumber(int samplesNumber);
	void setCSCycles(Vec2 t0, Vec2 t1);

	/**
	 * @brief synthesize_stationary uses the first texture added to synthesize a texture using a by-example cyclostationary cycloergodic synthesis.
	 * You need to have added at least two TandBFunctionType using addFunction for it to produce good results.
	 * Cyclostationary synthesis requires setting additional parameters (cycles and polyphase component samples).
	 * @return the synthesized texture.
	 */
	ImageType synthesize_cyclostationary();

	//Importance sampling DLC
	void setImportanceSampler(Stamping::SamplerImportance *sampler);


private:

	Vec2 fract(Vec2 v) const;

	Vec2 hash(const Vec2 &p) const;
	Vec2 cyclicHash(Vec2 p) const;
	Vec2 importanceHash(const Vec2 &p) const;

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

	//Cyclostationary DLC
	Vec2i m_cs_polyphaseComponentSamplesNumber;
	Vec2 m_cs_t0;
	Vec2 m_cs_t1;

	PtrImageType			m_transferPtrImage;
	PtrImageType			m_pcaPtrImage;
	unsigned				m_largestCycleProduct;

	//Importance sampling DLC
	Stamping::SamplerImportance *m_sampler;
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
	m_randomAffineTransform(),
	m_sampler(nullptr)
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
		//The mean is literally 0 but at least this demonstrates what you would do WITHOUT a PCA
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
			Vec2 uvTexture;
			if(m_sampler != nullptr)
				uvTexture = hash(tile.cast<double>());
			else
				uvTexture = importanceHash(tile.cast<double>());
			uvTexture[0] += u;
			uvTexture[1] += v;
			AffineTransform affT = m_randomAffineTransform.generate(true, cantorPairingFunction(tile));
			Vec2 uvTextureFract = CSN::CSN_Texture<I>::fract(affT.transform()*uvTexture);
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
			pix += mean * blendingValue;
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
ImageRGBd TilingAndBlending<I>::visualizePrimalAndDual()
{
	m_texturePool.generate();
	const ImageType &texture = m_texturePool[0].texture;

	ImageRGBd output;
	output.initItk(m_width, m_height);
	output.for_all_pixels([&] (ImageRGBd::PixelType &pix, int x, int y)
	{
		double u = double(x)/(texture.width()-1);
		double v = double(y)/(texture.height()-1);
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
		unsigned int i = 0;
		for(auto const &fonctions : m_functions)
		{
			Eigen::Vector2i tile = fonctions.tilingFunction(u * m_uScale, v * m_vScale);
			Vec2 uvTexture = hash(tile.cast<double>());
			uvTexture[0] += u;
			uvTexture[1] += v;
			AffineTransform affT = m_randomAffineTransform.generate(true, cantorPairingFunction(tile));
			Vec2 uvTextureFract = CSN::CSN_Texture<I>::fract(affT.transform()*uvTexture);
			PixelPosType xyTexture;
			xyTexture[0] = int(uvTextureFract[0] * texture.width()) % texture.width();
			xyTexture[1] = int(uvTextureFract[1] * texture.height()) % texture.height();
			double blendingValue = fonctions.blendingFunction(u * m_uScale, v * m_vScale) / sum;
			pix[i++] = blendingValue*blendingValue;
			sumSquare += blendingValue*blendingValue;
		}
		if(m_singularityBlendingFunction != nullptr)
		{
			double blendingValue = m_singularityBlendingFunction(u * m_uScale, v * m_vScale) / sum;
			pix[i++] = blendingValue*blendingValue;
			sumSquare += blendingValue*blendingValue;
		}
		sumSquare = std::sqrt(sumSquare);
		pix[0] += 0.2;
		pix[1] += 0.2;
		for(unsigned int i=0; i<3; ++i)
		{
			pix[i] /= sumSquare;
			pix[i] = std::max(std::min(1.0, pix[i]), 0.0);
		}
	});
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
		Vec2 uvTexture(u, v);
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
			Vec2 uvTextureFract = CSN::CSN_Texture<I>::fract(uvTexture);
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
			pix[i] /= sumSquare;
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
typename TilingAndBlending<I>::Vec2 TilingAndBlending<I>::hash(const Vec2 &p) const
{
	Eigen::Matrix2d hashMat;
	hashMat << 127.1, 269.5, 311.7, 183.3;
	Vec2 q = hashMat * p;
	q[0] = sin(q[0]);
	q[1] = sin(q[1]);
	return CSN::CSN_Texture<I>::fract ( q * 43758.5453 );
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
		u = double(x)/(visualization.width());
		u = u*uMax;
		v = double(y)/(visualization.height());
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
		u = double(x)/(visualization.width());
		u = u*uMax;
		v = double(y)/(visualization.height());
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
	m_texturePool.generate();
	ImageType &exemplar = m_texturePool[0].texture;
	double naturalScaleU = double(width)/exemplar.width(), naturalScaleV = double(height)/exemplar.height();
	ImageRGBd visualization;
	visualization.initItk(width, height, true);
	visualization.for_all_pixels([&] (ImageRGBd::PixelType &pix, int x, int y)
	{
		double u = double(x)/(visualization.width()) * naturalScaleU;
		double v = double(y)/(visualization.height()) * naturalScaleV;
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
			Vec2 uvTexture = hash(transition.cast<double>());
			uvTexture[0] += u * m_uScale;
			uvTexture[1] += v * m_vScale;
			Vec2 uvTextureFract = CSN::CSN_Texture<I>::fract(uvTexture);
			double blendingValue = fonctions.blendingFunction(u * m_uScale, v * m_vScale) / sum;
			double colorR = uvTextureFract[0] > 0.5 ? 1.0 - ((uvTextureFract[0] - 0.5) * 2.0) : uvTextureFract[0] * 2.0;
			double colorG = uvTextureFract[1] > 0.5 ? 1.0 - ((uvTextureFract[1] - 0.5) * 2.0) : uvTextureFract[1] * 2.0;
			pix[0] += colorR * blendingValue;
			pix[1] += colorG * blendingValue;
		}
	});
	return visualization;
}

template<typename I>
int TilingAndBlending<I>::cantorPairingFunction(const Eigen::Vector2i &v)
{
	return (v[0] + v[1]) * (v[0] + v[1] + 1) / 2 + v[0];
}

template<typename I>
void TilingAndBlending<I>::setCSPolyphaseComponentSamplesNumber(Vec2i samplesNumber)
{
	m_cs_polyphaseComponentSamplesNumber = samplesNumber;
}

template<typename I>
void TilingAndBlending<I>::setCSPolyphaseComponentSamplesNumber(int samplesNumber)
{
	m_cs_polyphaseComponentSamplesNumber[0] = samplesNumber;
	m_cs_polyphaseComponentSamplesNumber[1] = samplesNumber;
}

template<typename I>
void TilingAndBlending<I>::setCSCycles(Vec2 t0, Vec2 t1)
{
	m_cs_t0 = t0;
	m_cs_t1 = t1;
}

template<typename I>
typename TilingAndBlending<I>::ImageType TilingAndBlending<I>::synthesize_cyclostationary()
{
	std::cout << "Init" << std::endl;
	m_texturePool.generate();
	const ImageType &texture = m_texturePool[0].texture;
	PcaImageType pcaTexture = CSN::CSN_Texture<ImageType>::toPcaImageType(texture);
	int exemplarWidth = pcaTexture.width(), exemplarHeight = pcaTexture.height();

	//Input treatment -- can be written in another function if you want to separate pre-computations and synthesis

	//////////////////////////////////////////////////////////////////////////////
	///Forward PCA -- 1 per polyphase component. I picked it from CSN and don't remember why there are so many 0.5 in the code
	PcaImageType pcaOutput;
	pcaOutput.initItk(m_width, m_height);
	int pcaWidth, pcaHeight, GWidth, GHeight;
	//Representing the size of the image for which a single PCA is computed, i.e. one polyphase component
	pcaWidth = int(std::floor(1.0/m_cs_t0[0]));
	pcaHeight = int(std::floor(1.0/m_cs_t1[1]));
	GWidth = pcaWidth*m_cs_polyphaseComponentSamplesNumber[0];
	GHeight = pcaHeight*m_cs_polyphaseComponentSamplesNumber[1];
	PcaImageType pcaSubTexture, pcaGaussianSubTexture, pcaInputCopy;
	pcaInputCopy.initItk(pcaTexture.width(), pcaTexture.height());
	pcaInputCopy.copy_pixels(pcaTexture);
	pcaSubTexture.initItk(GWidth, GHeight);
	pcaGaussianSubTexture.initItk(GWidth, GHeight);
	m_pcaPtrImage.initItk(	std::round((m_cs_t0[0])*(pcaTexture.width()) + 0.5),
							std::round((m_cs_t1[1])*(pcaTexture.height()) + 0.5) );
	std::cout << "PcaPtrImage allocated" << std::endl;
	m_pcaPtrImage.for_all_pixels([&] (PtrImageType::PixelType &ptr, int inX, int inY)
	{
		//std::cout << "trying new sub PCA: inX=" << inX << ", inY=" << inY << std::endl;
		auto getColor = [&] (int subX, int subY, int wx, int wy) -> PcaPixelType
		{
			double x, y;
			x = std::max(std::round(inX + (subX*m_cs_t0[0] + subY*m_cs_t1[0])*(pcaInputCopy.width())) + (wx - m_cs_polyphaseComponentSamplesNumber[0]/2.0 + 0.5)*1.0, 0.0);
			y = std::max(std::round(inY + (subX*m_cs_t0[1] + subY*m_cs_t1[1])*(pcaInputCopy.height())) + (wy - m_cs_polyphaseComponentSamplesNumber[1]/2.0 + 0.5)*1.0, 0.0);
			return bilinear_interpolation(pcaInputCopy, x, y, true);
		};
		for(int subX=0; subX<pcaWidth; ++subX)
			for(int subY=0; subY<pcaHeight; ++subY)
				for(unsigned wx=0; wx<m_cs_polyphaseComponentSamplesNumber[0]; ++wx)
					for(unsigned wy=0; wy<m_cs_polyphaseComponentSamplesNumber[1]; ++wy)
					{
						pcaSubTexture.pixelAbsolute(subX*m_cs_polyphaseComponentSamplesNumber[0]+wx, subY*m_cs_polyphaseComponentSamplesNumber[1]+wy) = getColor(subX, subY, int(wx), int(wy));
					}
		PcaType *subPca = new PcaType(pcaSubTexture);
		ptr = reinterpret_cast<PtrImageType::PixelType>(subPca);
		MaskBool mb_alwaysTrue(pcaSubTexture.width(), pcaSubTexture.height());
		mb_alwaysTrue |= [] (int, int) {return true;};
		subPca->computePCA(mb_alwaysTrue);
		subPca->project(pcaSubTexture);

		for(int subX=0; subX<pcaWidth; ++subX)
			for(int subY=0; subY<pcaHeight; ++subY)
			{
				PixelPosType inCoordinates;
				inCoordinates[0] = int(std::round(inX + (subX*m_cs_t0[0] + subY*m_cs_t1[0])*(pcaInputCopy.width())))%pcaInputCopy.width();
				inCoordinates[1] = int(std::round(inY + (subX*m_cs_t0[1] + subY*m_cs_t1[1])*(pcaInputCopy.height())))%pcaInputCopy.height();
				pcaTexture.pixelAbsolute(inCoordinates) = pcaSubTexture.pixelAbsolute(	subX*m_cs_polyphaseComponentSamplesNumber[0]+int(m_cs_polyphaseComponentSamplesNumber[0]/2.0-0.5),
																						subY*m_cs_polyphaseComponentSamplesNumber[1]+int(m_cs_polyphaseComponentSamplesNumber[1]/2.0-0.5));
			}
	});
	std::cout << "Pca ptr image initialized" << std::endl;
	//////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////
	///Histogram transfer -- 1 per polyphase component again. And again with the 0.5s
	PcaImageType pcaGaussianTexture;
	GaussianTransferType gtt;
	//Mean computation (in reality all transfers set the mean to 0.5 for visualization purposes so we take that in account, sure, why not)
	PcaPixelType mean = PcaImageType::zero();
	if(m_useTransfer)
	{
		for(int i = 0; i<3; ++i)
		{
			mean[i] = 0.5;
		}
	}
	pcaGaussianTexture.initItk(pcaTexture.width(), pcaTexture.height());
	if(m_useTransfer)
	{
		//Same as PCA; if using cycles, there is one Gaussianization (and one histogram) for each polyphase component.
		int transferWidth, transferHeight, GWidth, GHeight;
		transferWidth = int(std::floor(1.0/m_cs_t0[0]));
		transferHeight = int(std::floor(1.0/m_cs_t1[1]));
		GWidth = transferWidth*m_cs_polyphaseComponentSamplesNumber[0];
		GHeight = transferHeight*m_cs_polyphaseComponentSamplesNumber[1];
		PcaImageType pcaSubTexture, pcaGaussianSubTexture, pcaInputCopy;
		pcaInputCopy.initItk(pcaTexture.width(), pcaTexture.height());
		pcaInputCopy.copy_pixels(pcaTexture);
		pcaSubTexture.initItk(GWidth, GHeight);
		pcaGaussianSubTexture.initItk(GWidth, GHeight);
		m_transferPtrImage.initItk(	std::round((m_cs_t0[0])*(pcaTexture.width()) +0.5),
									std::round((m_cs_t1[1])*(pcaTexture.height()) +0.5) );
		m_transferPtrImage.for_all_pixels([&] (PtrImageType::PixelType &ptr, int inX, int inY)
		{
			LutType *subLut = new LutType;
			subLut->initItk(256, 1);
			ptr = reinterpret_cast<PtrImageType::PixelType>(subLut);
			auto getColor = [&] (int subX, int subY, int wx, int wy) -> PcaPixelType
			{
				double x, y;
				x = std::max(std::round(inX + (subX*m_cs_t0[0] + subY*m_cs_t1[0])*(pcaInputCopy.width())) + (wx - m_cs_polyphaseComponentSamplesNumber[0]/2.0 + 0.5)*1.0, 0.0);
				y = std::max(std::round(inY + (subX*m_cs_t0[1] + subY*m_cs_t1[1])*(pcaInputCopy.height())) + (wy - m_cs_polyphaseComponentSamplesNumber[1]/2.0 + 0.5)*1.0, 0.0);
				return bilinear_interpolation(pcaInputCopy, x, y, true);
			};
			for(int subX=0; subX<transferWidth; ++subX)
				for(int subY=0; subY<transferHeight; ++subY)
					for(unsigned wx=0; wx<m_cs_polyphaseComponentSamplesNumber[0]; ++wx)
						for(unsigned wy=0; wy<m_cs_polyphaseComponentSamplesNumber[1]; ++wy)
						{
							pcaSubTexture.pixelAbsolute(subX*m_cs_polyphaseComponentSamplesNumber[0]+wx, subY*m_cs_polyphaseComponentSamplesNumber[1]+wy) = getColor(subX, subY, wx, wy);
						}
			gtt.ComputeTinput(pcaSubTexture, pcaGaussianSubTexture);
			gtt.ComputeinvT(pcaSubTexture, *subLut);
			for(int subX=0; subX<transferWidth; ++subX)
				for(int subY=0; subY<transferHeight; ++subY)
				{
					PixelPosType inCoordinates;
					inCoordinates[0] = int(std::round(inX + (subX*m_cs_t0[0] + subY*m_cs_t1[0])*(pcaInputCopy.width())))%pcaInputCopy.width();
					inCoordinates[1] = int(std::round(inY + (subX*m_cs_t0[1] + subY*m_cs_t1[1])*(pcaInputCopy.height())))%pcaInputCopy.height();
					pcaGaussianTexture.pixelAbsolute(inCoordinates) = pcaGaussianSubTexture.pixelAbsolute(subX*m_cs_polyphaseComponentSamplesNumber[0]+int(m_cs_polyphaseComponentSamplesNumber[0]/2.0-0.5),
																										  subY*m_cs_polyphaseComponentSamplesNumber[1]+int(m_cs_polyphaseComponentSamplesNumber[1]/2.0-0.5));
				}
		});
	}
	else
	{
		pcaGaussianTexture.copy_pixels(pcaTexture);
	}
	std::cout << "Transfer ptr image initialized" << std::endl;

	//Actual synthesis, using pre-computed data
	pcaOutput.for_all_pixels([&] (PcaPixelType &pix, int x, int y)
	{
		double u = double(x)/(pcaTexture.width());
		double v = double(y)/(pcaTexture.height());
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
			Vec2 uvTexture = cyclicHash(tile.cast<double>()); //The only thing that really changes in the synthesis function with respect to the stationary version
			uvTexture[0] += u;
			uvTexture[1] += v;
			AffineTransform affT = m_randomAffineTransform.generate(true, cantorPairingFunction(tile));
			Vec2 uvTextureFract = CSN::CSN_Texture<I>::fract(affT.transform()*uvTexture);
			PixelPosType xyTexture;
			xyTexture[0] = int(std::round(uvTextureFract[0] * (pcaTexture.width()))) % pcaTexture.width();
			xyTexture[1] = int(std::round(uvTextureFract[1] * (pcaTexture.height()))) % pcaTexture.height();
			double blendingValue = fonctions.blendingFunction(u * m_uScale, v * m_vScale) / sum;
			pix += pcaGaussianTexture.pixelAbsolute(xyTexture) * blendingValue;
			sumSquare += blendingValue*blendingValue;
		}
		if(m_singularityBlendingFunction != nullptr)
		{
			double blendingValue = m_singularityBlendingFunction(u * m_uScale, v * m_vScale) / sum;
			pix += mean * blendingValue;
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

	/////////////////////
	//Output post-process
	/////////////////////

	///Inverse transfer
	if(m_useTransfer)
	{
		unsigned transferWidth, transferHeight, GWidth, GHeight;
		transferWidth = unsigned(std::ceil(double(m_width)/exemplarWidth) * std::floor(1.0/m_cs_t0[0]));
		transferHeight = unsigned(std::ceil(double(m_height)/exemplarHeight) * std::floor(1.0/m_cs_t1[1]));
		GWidth = transferWidth*m_cs_polyphaseComponentSamplesNumber[0];
		GHeight = transferHeight*m_cs_polyphaseComponentSamplesNumber[1];
		PcaImageType pcaSubTexture, pcaGaussianSubTexture, pcaOutputCopy;
		pcaOutputCopy.initItk(pcaOutput.width(), pcaOutput.height());
		pcaOutputCopy.copy_pixels(pcaOutput);
		pcaSubTexture.initItk(GWidth, GHeight);
		m_transferPtrImage.for_all_pixels([&] (PtrImageType::PixelType &ptr, int inX, int inY)
		{
			LutType *subLut = reinterpret_cast<LutType *>(ptr);
			pcaSubTexture.for_all_pixels([&] (typename PcaImageType::PixelType &pix, int subX, int subY)
			{
				PixelPosType inCoordinates;
//				if(unsigned(std::ceil(double(m_width)/pcaTexture.width()) * std::floor(1.0/m_cs_t0[0])))
				inCoordinates[0] = int(std::round(inX + (subX*m_cs_t0[0] + subY*m_cs_t1[0])*(exemplarWidth)))%pcaOutputCopy.width();
				inCoordinates[1] = int(std::round(inY + (subX*m_cs_t0[1] + subY*m_cs_t1[1])*(exemplarHeight)))%pcaOutputCopy.height();
				pix = pcaOutputCopy.pixelAbsolute(inCoordinates);
			});
			pcaGaussianSubTexture = GaussianTransferType::invT(pcaSubTexture, *subLut);
			pcaGaussianSubTexture.for_all_pixels([&] (const typename PcaImageType::PixelType &pix, int subX, int subY)
			{
				PixelPosType inCoordinates;
				inCoordinates[0] = int(std::round(inX + (subX*m_cs_t0[0] + subY*m_cs_t1[0])*(exemplarWidth)))%pcaOutputCopy.width();
				inCoordinates[1] = int(std::round(inY + (subX*m_cs_t0[1] + subY*m_cs_t1[1])*(exemplarHeight)))%pcaOutputCopy.height();
				pcaOutput.pixelAbsolute(inCoordinates) = pix;
			});
		});
	}
//	pcaOutput.for_all_pixels([&] (typename PcaImageType::PixelType &pix)
//	{
//		pix[0] += 0.5;
//		pix[1] += 0.5;
//		pix[2] += 0.5;
//	});
	IO::save01_in_u8(pcaOutput, "/home/nlutz/pca.png");
	///Inverse PCA
	pcaTexture.initItk(pcaOutput.width(), pcaOutput.height());
	pcaWidth = unsigned(std::ceil(double(m_width)/exemplarWidth * std::floor(1.0/m_cs_t0[0])));
	pcaHeight = unsigned(std::ceil(double(m_height)/exemplarHeight * std::floor(1.0/m_cs_t1[1])));
	GWidth = pcaWidth*m_cs_polyphaseComponentSamplesNumber[0];
	GHeight = pcaHeight*m_cs_polyphaseComponentSamplesNumber[1];
	PcaImageType pcaUnfoldedSubTexture, pcaOutputCopy;
	pcaOutputCopy.initItk(pcaOutput.width(), pcaOutput.height());
	pcaOutputCopy.copy_pixels(pcaOutput);
	pcaSubTexture.initItk(GWidth, GHeight);
	pcaUnfoldedSubTexture.initItk(pcaOutput.width(), pcaOutput.height());
	m_pcaPtrImage.for_all_pixels([&] (PtrImageType::PixelType &ptr, int inX, int inY)
	{
		PcaType *subPca = reinterpret_cast<PcaType *>(ptr);
		pcaSubTexture.for_all_pixels([&] (typename PcaImageType::PixelType &pix, int subX, int subY)
		{
			PixelPosType inCoordinates;
			inCoordinates[0] = int(std::round(inX + (subX*m_cs_t0[0] + subY*m_cs_t1[0])*(exemplarWidth)))%pcaOutputCopy.width();
			inCoordinates[1] = int(std::round(inY + (subX*m_cs_t0[1] + subY*m_cs_t1[1])*(exemplarHeight)))%pcaOutputCopy.height();
			pix = pcaOutputCopy.pixelAbsolute(inCoordinates);
		});

		subPca->back_project(pcaSubTexture, pcaUnfoldedSubTexture);

		pcaUnfoldedSubTexture.for_all_pixels([&] (const typename PcaImageType::PixelType &pix, int subX, int subY)
		{
			PixelPosType inCoordinates;
			inCoordinates[0] = int(std::round(inX + (subX*m_cs_t0[0] + subY*m_cs_t1[0])*(exemplarWidth)))%pcaOutputCopy.width();
			inCoordinates[1] = int(std::round(inY + (subX*m_cs_t0[1] + subY*m_cs_t1[1])*(exemplarHeight)))%pcaOutputCopy.height();
			pcaTexture.pixelAbsolute(inCoordinates) = pix;
		});
		std::cout << pcaUnfoldedSubTexture.width() << std::endl;
	});
	if(m_useTransfer)
	{
		m_transferPtrImage.for_all_pixels([&] (PtrImageType::PixelType &pix)
		{
			GaussianTransferType *gtt = reinterpret_cast<GaussianTransferType *>(pix);
			delete gtt;
		});
	}
	m_pcaPtrImage.for_all_pixels([&] (PtrImageType::PixelType &pix)
	{
		PcaType *pca = reinterpret_cast<PcaType *>(pix);
		delete pca;
	});
	ImageType output =  CSN::CSN_Texture<ImageType>::fromPcaImageType(pcaTexture);
	return output;
}

template<typename I>
typename TilingAndBlending<I>::Vec2 TilingAndBlending<I>::fract(Vec2 v) const
{
	Eigen::Vector2d w;
	w[0] = v[0]-std::floor(v[0]);
	w[1] = v[1]-std::floor(v[1]);
	return w;
}

template<typename I>
typename TilingAndBlending<I>::Vec2 TilingAndBlending<I>::cyclicHash(Vec2 p) const
{
	int randMax = m_largestCycleProduct;
	Eigen::Matrix2d hashMat;
	hashMat << 127.1, 269.5, 311.7, 183.3;
	Eigen::Vector2d q = hashMat * p;
	q[0] = sin(q[0]);
	q[1] = sin(q[1]);
	Eigen::Vector2d h = fract ( q * 43758.5453 );
	int cycle1 = int(h[0]*randMax);
	int cycle2 = int(h[1]*randMax);
	return double(cycle1)*m_cs_t0 + double(cycle2)*m_cs_t1;
}

template<typename I>
typename TilingAndBlending<I>::Vec2 TilingAndBlending<I>::importanceHash(const Vec2 &p) const
{
	assert(m_sampler && "importanceHash: no sampler is set.");
	const ImageType &texture = m_texturePool[0].texture;
	srand(int(p[1]*(1<<14)*texture.width()) + int(p[0]*(1<<14)));
	Eigen::Vector2f pf = m_sampler->next();
	return Eigen::Vector2d(double(pf[0]), double(pf[1]));
}

}

#endif
