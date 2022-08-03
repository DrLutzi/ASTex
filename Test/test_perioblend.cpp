#include <stdlib.h>
#include <ASTex/Perioblend/perioblend.h>
#include <ASTex/easy_io.h>
#include <ASTex/fourier.h>
#include <ASTex/utils.h>
#include <ASTex/rpn_utils.h>
#include "itkImage.h"
#include "itkScaleTransform.h"
#include "itkResampleImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#define SQRT3DIV2 0.86602540378
#define SQRT3 1.73205080757
#define SQRT2DIV2 0.70710678118

typedef struct
{
	Eigen::Vector2d vectors[2];
} CyclePair;

using CycleMapType = std::map<std::string, CyclePair>;

/**
 * @brief loadCycles load a cycles file.
 * Format of the file (don't copy the = symbols):
 * [<texture name>					//texture name with no extension or directory
 * <denominator> <denominator>		//two denomionators of the first cycle (cycle is [1/denominator]^2)
 * <denominator> <denominator>]*	//two denomionators of the second cycle
 * =======
 * Exemple:
 * =======
 * bricks5
 * 14	20
 * 7	0
 *
 * herringbone2
 * 7	0
 * 0	26
 * =======
 */
CycleMapType loadCycles(std::string filename)
{
	CycleMapType knownCycles;
	std::string name;
	CyclePair cycles;
	std::ifstream ifs(filename);
	if(!ifs)
	{
		return knownCycles;
	}
	while(!ifs.eof())
	{
		ifs >> name;
		if(!ifs.eof()) //I hate C++ sometimes
		{
			double readNumber1, readNumber2;
			ifs >> readNumber1 >> readNumber2;
			if((readNumber1>0 && readNumber1<1) || (readNumber2>0 && readNumber2<1))
				cycles.vectors[0] = Eigen::Vector2d(readNumber1, readNumber2);
			else
				cycles.vectors[0] = Eigen::Vector2d(readNumber1 == 0 ? 0 : 1.0/readNumber1,
													readNumber2 == 0 ? 0 : 1.0/readNumber2);
			ifs >> readNumber1 >> readNumber2;
			if((readNumber1>0 && readNumber1<1) || (readNumber2>0 && readNumber2<1))
				cycles.vectors[1] = Eigen::Vector2d(readNumber1, readNumber2);
			else
				cycles.vectors[1] = Eigen::Vector2d(readNumber1 == 0 ? 0 : 1.0/readNumber1,
													readNumber2 == 0 ? 0 : 1.0/readNumber2);
			knownCycles.insert({name, cycles});
		}
	}
	return knownCycles;
}

typedef enum {NO_ISOTROPY=0, FULL_ISOTROPY=1, PI_ISOTROPY=2, HALF_PI_ISOTROPY=3, UNSPECIFIED_ISOTROPY=4} Isotropy_enum_t;
typedef enum {SQUARE_SQUARE=0, ALTSQUARE_TRIANGULAR=1, ROTSQUARE_ROTSQUARE=2, HEXAGONAL_TRIANGULAR=3, FOUR_SQUARES=4} Tiling_enum_t;

typedef struct
{
	bool useHistogramTransfer;			//if a histogram transfer should be used
	unsigned int width;					//output width
	unsigned int height;				//output height
	Isotropy_enum_t isotropy;			//type of affine invariance (used to be isotropy only)
	double uScale;						//scale in u (only for synthesis structures)
	double vScale;						//scale in v (only for synthesis structures)
	double uTranslation;				//translation in u (only for synthesis structures)
	double vTranslation;				//translation in v (only for synthesis structures)
	bool useSingularityHider;			//if a singularity hider should be used, but it's only implemented for SqSq tiles
	Tiling_enum_t tiling;				//the type of primal/dual
	bool isCyclostationary;				//if the process must be assumed to be cyclostationary (in which case periods will be seeked)
	bool useImportanceSampling;			//(UNUSED) if importance sampling should be used
	bool useAutocovarianceSampling;		//if importance sampling with the AC fonction should be used
	bool compareAutocovariance;			//if the autocovariance of the output must be compared to that of the exemplar
	bool cheatAutocovariance;			//if the autocovariance function should be inverted in the importance sampler to force some of the worst possible outputs
	bool tryToUseMask;					//if a mask should be loaded if possible (mask image should be texture_mask.png in the same directory)
	bool useCustomImportanceMap;
	double importanceMapCutThreshold;
	bool computeOnlyOutputAutocovariance;
	int contentMaskTest;
} ArgumentsType;

ArgumentsType loadArguments(std::string argFilename)
{
	std::ifstream ifs(argFilename);
	assert(ifs);

	ArgumentsType arguments;

	std::string line;
	while( std::getline(ifs, line) )
	{
		std::istringstream is_line(line);
		std::string key;
		if( std::getline(is_line, key, '=') )
		{
			std::string value;
			if( std::getline(is_line, value) )
			{
				if(key == "useHistogramTransfer")
				{
					arguments.useHistogramTransfer = std::stoi(value);
				}
				else if(key == "width")
				{
					arguments.width = std::stoi(value);
				}
				else if(key == "height")
				{
					arguments.height = std::stoi(value);
				}
				else if(key == "isotropy")
				{
					int vInt = std::stoi(value);
					if(vInt <= 4 && vInt >= 0)
					{
						arguments.isotropy = Isotropy_enum_t(vInt);
					}
					else
					{
						arguments.isotropy = UNSPECIFIED_ISOTROPY;
					}
				}
				else if(key == "uScale")
				{
					arguments.uScale = std::stod(value);
				}
				else if(key == "vScale")
				{
					arguments.vScale = std::stod(value);
				}
				else if(key == "useSingularityHider")
				{
					arguments.useSingularityHider = std::stoi(value);
				}
				else if(key == "tiling")
				{
					arguments.tiling = Tiling_enum_t(std::stoi(value));
				}
				else if(key == "uTranslation")
				{
					arguments.uTranslation = std::stod(value);
				}
				else if(key == "vTranslation")
				{
					arguments.vTranslation = std::stod(value);
				}
				else if(key == "isCyclostationary")
				{
					arguments.isCyclostationary = bool(std::stoi(value));
				}
				else if(key == "useImportanceSampling")
				{
					arguments.useImportanceSampling = bool(std::stoi(value));
				}
				else if(key == "useAutocovarianceSampling")
				{
					arguments.useAutocovarianceSampling = bool(std::stoi(value));
				}
				else if(key == "compareAutocovariance")
				{
					arguments.compareAutocovariance = bool(std::stoi(value));
				}
				else if(key == "cheatAutocovariance")
				{
					arguments.cheatAutocovariance = bool(std::stoi(value));
				}
				else if(key == "tryToUseMask")
				{
					arguments.tryToUseMask = bool(std::stoi(value));
				}
				else if(key == "useCustomImportanceMap")
				{
					arguments.useCustomImportanceMap = bool(std::stoi(value));
				}
				else if(key == "importanceMapCutThreshold")
				{
					arguments.importanceMapCutThreshold = std::stod(value);
				}
				else if(key == "computeOnlyOutputAutocovariance")
				{
					arguments.computeOnlyOutputAutocovariance = bool(std::stoi(value));
				}
				else if(key == "contentMaskTest")
				{
					arguments.contentMaskTest = std::stoi(value);
				}
			}
		}
	}
	return arguments;
}

int mod3( int n )
{
	return (n<0) ? 2-((2-n)%3) : n%3;
}

Eigen::Vector2d floor(Eigen::Vector2d q)
{
	return Eigen::Vector2d(std::floor(q.x()), std::floor(q.y()));
}

Eigen::Vector3d floor(Eigen::Vector3d q)
{
	return Eigen::Vector3d(std::floor(q.x()), std::floor(q.y()), std::floor(q.z()));
}

Eigen::Vector2d fract(Eigen::Vector2d q)
{
	return Eigen::Vector2d(q.x() - std::floor(q.x()), q.y() - std::floor(q.y()));
}

Eigen::Vector3d fract(Eigen::Vector3d q)
{
	return Eigen::Vector3d(q.x() - std::floor(q.x()), q.y() - std::floor(q.y()), q.z() - std::floor(q.z()));
}

Eigen::Vector2i hexagonGetID(Eigen::Vector2d p )
{
	Eigen::Vector2d  q = Eigen::Vector2d( p.x(), p.x()*0.5+p.y()*SQRT3DIV2);

	Eigen::Vector2i i = floor(q).cast<int>();
	Eigen::Vector2d  f =       fract(q);

	int v = mod3(i.x()+i.y());
	Eigen::Vector2i id = i + Eigen::Vector2i(v, v); //TODO : v, 0 ou v, v ?
	if( v==2 ) id -= (f.x()>f.y())?Eigen::Vector2i(1,2):Eigen::Vector2i(2,1);

	return Eigen::Vector2i( id.x(), (2*id.y()-id.x())/3 );
}

// return the center of an hexagon
Eigen::Vector2d hexagonCenFromID(Eigen::Vector2i id)
{
	return Eigen::Vector2d(float(id.x()),float(id.y())*SQRT3);
}
//////////////////////////////////////////////////////////////////////////

float hexagonDistanceFromEdge(Eigen::Vector2i id, Eigen::Vector2d p)
{
	Eigen::Vector2d c = hexagonCenFromID(id);
	Eigen::Matrix2d rotate60; rotate60 << 0.5, -SQRT3DIV2, SQRT3DIV2, 0.5;
	Eigen::Matrix2d rotate120; rotate120 << -0.5, -SQRT3DIV2, SQRT3DIV2, -0.5;
	Eigen::Vector2d cen = p - c;
	Eigen::Vector2d cen60 = rotate60*(p - c);
	Eigen::Vector2d cen120 = rotate120*(p - c);
	float dist1 = cen.x()*cen.x();
	float dist2 = cen60.x()*cen60.x();
	float dist3 = cen120.x()*cen120.x();
	return 1.0 - std::max(std::max(dist1, dist2), dist3);
}

Eigen::Vector3d triangleCoordinates(Eigen::Vector2d UV)
{
	Eigen::Vector3d t;
	Eigen::Matrix2d rotate60; rotate60 << 0.5, -SQRT3DIV2, SQRT3DIV2, 0.5;
	Eigen::Matrix2d rotate30; rotate30 << SQRT3DIV2, -0.5, 0.5, SQRT3DIV2;
	t.x() = UV.y();
	t.y() = (rotate30 * UV).x();
	t.z() = (rotate60 * UV).y();
	return t;
}

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " <in_texture> [in_texture2]..." << std::endl;
		return EXIT_FAILURE;
	}
	ArgumentsType arguments=loadArguments(argv[1]);
	unsigned int visualizationWidth = 256;
	unsigned int visualizationHeight = 128;
	double visualizationU = 4.0;
	double visualizationV = 2.0;

	using ImageType = ASTex::ImageRGBd;
	auto divideBySum = [&] (const ASTex::ImageGrayd &blendingFunction, const ASTex::ImageGrayd &sumFunction) -> ImageGrayd
	{
		ASTex::ImageGrayd output_dbs;
		output_dbs.initItk(blendingFunction.width(), blendingFunction.height());
		output_dbs.for_all_pixels([&] (ASTex::ImageGrayd::PixelType &pix, int x, int y)
		{
			pix = blendingFunction.pixelAbsolute(x, y) / sumFunction.pixelAbsolute(x, y);
		});
		return output_dbs;
	};


	auto noTransition = [&] (double u, double v) -> Eigen::Vector2i
	{
		Eigen::Vector2i vec;
		return vec;
	};
	auto noTransition2 = [&] (double u, double v) -> Eigen::Vector2i
	{
		Eigen::Vector2i vec;
		vec[0] = 1;
		vec[1] = 1;
		return vec;
	};

//	BlendingFunctionType bft1;
//	bft1.blendingFunction = [&] (double u, double v)
//	{
//		return (std::sin(u*M_PI)+1) * (std::sin(v*M_PI)+1);
//	};
//	bft1.transitionFunction = [&] (double u, double v)
//	{
//		Eigen::Vector2i vec;
//		vec[0] = int((u+0.5)/2)*63;
//		vec[1] = int((v+0.5)/2)*63;
//		return vec;
//	};
	TandBFunctionType bft1;
	bft1.blendingFunction = [&] (double u, double v)
	{
		u += arguments.uTranslation;
		v += arguments.vTranslation;
		float sinU = std::sin(u*M_PI);
		float sinV = std::sin(v*M_PI);
		return sqrt(sinU*sinU * sinV*sinV);
	};
	bft1.tilingFunction = [&] (double u, double v)
	{
		u += arguments.uTranslation;
		v += arguments.vTranslation;
		Eigen::Vector2i vec;
		vec[0] = int(u)*63;
		vec[1] = int(v)*63;
		return vec;
	};
	ImageGrayd bft1Blending = bft1.visualizeBlending(visualizationWidth, visualizationHeight, visualizationU, visualizationV);
	ImageRGBd bft1Transition = bft1.visualizeTiling(visualizationWidth, visualizationHeight, visualizationU, visualizationV);
	//IO::save01_in_u8(bft1Blending, std::string("/home/nlutz/bft1Blending.png"));
	//IO::save01_in_u8(bft1Transition, std::string("/home/nlutz/bft1Transition.png"));

	TandBFunctionType bft2;
	bft2.blendingFunction = [&] (double u, double v)
	{
		u += arguments.uTranslation;
		v += arguments.vTranslation;
		u += 0.5;
		v += 0.5;
		float sinU = std::sin(u*M_PI);
		float sinV = std::sin(v*M_PI);
		return sqrt(sinU*sinU * sinV*sinV);
	};
	bft2.tilingFunction = [&] (double u, double v)
	{
		u += arguments.uTranslation;
		v += arguments.vTranslation;
		Eigen::Vector2i vec;
		vec[0] = int(floor(u + 0.5)*127);
		vec[1] = int(floor(v + 0.5)*127);
		return vec;
	};

	ImageGrayd bft2Blending = bft2.visualizeBlending(visualizationWidth, visualizationHeight, visualizationU, visualizationV);
	ImageRGBd bft2Transition = bft2.visualizeTiling(visualizationWidth, visualizationHeight, visualizationU, visualizationV);
	//IO::save01_in_u8(bft2Blending, std::string("/home/nlutz/bft2Blending.png"));
	//IO::save01_in_u8(bft2Transition, std::string("/home/nlutz/bft2Transition.png"));

	TandBFunctionType bft3;
	bft3.blendingFunction = [&] (double u, double v)
	{
		u += arguments.uTranslation;
		v += arguments.vTranslation;
		u+=0.5;
		float sinU = std::sin(u*M_PI);
		float sinV = std::sin(v*M_PI);
		return sqrt(sinU*sinU * sinV*sinV);
	};
	bft3.tilingFunction = [&] (double u, double v)
	{
		u += arguments.uTranslation;
		v += arguments.vTranslation;
		u+=0.5;
		Eigen::Vector2i vec;
		vec[0] = int(u)*63;
		vec[1] = int(v)*63;
		return vec;
	};
	ImageGrayd bft3Blending = bft3.visualizeBlending(visualizationWidth, visualizationHeight, visualizationU, visualizationV);
	ImageRGBd bft3Transition = bft3.visualizeTiling(visualizationWidth, visualizationHeight, visualizationU, visualizationV);
	//IO::save01_in_u8(bft3Blending, std::string("/home/nlutz/bft3Blending.png"));
	//IO::save01_in_u8(bft3Transition, std::string("/home/nlutz/bft3Transition.png"));

	TandBFunctionType bft4;
	bft4.blendingFunction = [&] (double u, double v)
	{
		u += arguments.uTranslation;
		v += arguments.vTranslation;
		v+=0.5;
		float sinU = std::sin(u*M_PI);
		float sinV = std::sin(v*M_PI);
		return sqrt(sinU*sinU * sinV*sinV);
	};
	bft4.tilingFunction = [&] (double u, double v)
	{
		u += arguments.uTranslation;
		v += arguments.vTranslation;
		v+=0.5;
		Eigen::Vector2i vec;
		vec[0] = int(u)*63;
		vec[1] = int(v)*63;
		return vec;
	};
	ImageGrayd bft4Blending = bft4.visualizeBlending(visualizationWidth, visualizationHeight, visualizationU, visualizationV);
	ImageRGBd bft4Transition = bft4.visualizeTiling(visualizationWidth, visualizationHeight, visualizationU, visualizationV);
	//IO::save01_in_u8(bft4Blending, std::string("/home/nlutz/bft4Blending.png"));
	//IO::save01_in_u8(bft4Transition, std::string("/home/nlutz/bft4Transition.png"));

	TandBFunctionType bftTri;
	bftTri.blendingFunction = [&] (double u, double v)
	{
		u += arguments.uTranslation;
		v += arguments.vTranslation;
		Eigen::Vector2d UV;
		UV[0] = u;
		UV[1] = v;
		Eigen::Vector3d t = fract(triangleCoordinates(UV));
		Eigen::Vector3d tmix = (t - Eigen::Vector3d(0.5, 0.5, 0.5))*2.0;
		tmix.x() = 1.0 - std::pow(tmix.x(), 2.0);
		tmix.y() = 1.0 - std::pow(tmix.y(), 2.0);
		tmix.z() = 1.0 - std::pow(tmix.z(), 2.0);
		return std::min(std::min(tmix.x(), tmix.y()), tmix.z());
	};
	bftTri.tilingFunction = [&] (double u, double v)
	{
		u += arguments.uTranslation;
		v += arguments.vTranslation;
		Eigen::Vector2d UV;
		UV[0] = u;
		UV[1] = v;
		Eigen::Vector3d t = floor(triangleCoordinates(UV));
		Eigen::Vector2d txy;
		txy.x() = t.x();
		txy.y() = t.y();
		return Eigen::Vector2i(TilingAndBlending<ImageType>::cantorPairingFunction(txy.cast<int>()), t.z());
	};

	ImageGrayd bftTriBlending = bftTri.visualizeBlending(visualizationWidth, visualizationHeight, visualizationU, visualizationV);
	ImageRGBd bftTriTransition = bftTri.visualizeTiling(visualizationWidth, visualizationHeight, visualizationU, visualizationV);
	//IO::save01_in_u8(bftTriBlending, std::string("/home/nlutz/bftTriBlending.png"));
	//IO::save01_in_u8(bftTriTransition, std::string("/home/nlutz/bftTriTransition.png"));

	TandBFunctionType bftAltTri;
	bftAltTri.blendingFunction = [&] (double u, double v)
	{
		u += arguments.uTranslation;
		v += arguments.vTranslation;
		u += 0.5;
		v += 0.5;
		Eigen::Vector2d UV;
		UV[0] = u / SQRT3DIV2;
		UV[1] = v; //+0.5
		Eigen::Vector3d t = fract(triangleCoordinates(UV));
		Eigen::Vector3d tmix = (t - Eigen::Vector3d(0.5, 0.5, 0.5))*2.0;
		tmix.x() = 1.0 - std::pow(tmix.x(), 2.0);
		tmix.y() = 1.0 - std::pow(tmix.y(), 2.0);
		tmix.z() = 1.0 - std::pow(tmix.z(), 2.0);
		return std::min(std::min(tmix.x(), tmix.y()), tmix.z());
	};
	bftAltTri.tilingFunction = [&] (double u, double v)
	{
		u += arguments.uTranslation;
		v += arguments.vTranslation;
		u += 0.5;
		v += 0.5;
		Eigen::Vector2d UV;
		UV[0] = u / SQRT3DIV2;
		UV[1] = v; //+0.5
		Eigen::Vector3d t = floor(triangleCoordinates(UV));
		Eigen::Vector2d txy;
		txy.x() = t.x();
		txy.y() = t.y();
		return Eigen::Vector2i(TilingAndBlending<ImageType>::cantorPairingFunction(txy.cast<int>()), t.z());
	};

	ImageGrayd bftAltTriBlending = bftAltTri.visualizeBlending(visualizationWidth, visualizationHeight, visualizationU, visualizationV);
	ImageRGBd bftAltTriTransition = bftAltTri.visualizeTiling(visualizationWidth, visualizationHeight, visualizationU, visualizationV);
	//IO::save01_in_u8(bftAltTriBlending, std::string("/home/nlutz/bftAltTriBlending.png"));
	//IO::save01_in_u8(bftAltTriTransition, std::string("/home/nlutz/bftAltTriTransition.png"));

	TandBFunctionType bftHex;
	bftHex.blendingFunction = [&] (double u, double v)
	{
		u += arguments.uTranslation;
		v += arguments.vTranslation;
		Eigen::Vector2d UV;
		UV[0] = u;
		UV[1] = v;
		UV *= SQRT3;
		Eigen::Vector2i H = hexagonGetID(UV);
		return hexagonDistanceFromEdge(H, UV);
	};
	bftHex.tilingFunction = [&] (double u, double v)
	{
		u += arguments.uTranslation;
		v += arguments.vTranslation;
		Eigen::Vector2d UV;
		UV[0] = u;
		UV[1] = v;
		UV *= SQRT3;
		Eigen::Vector2i H = hexagonGetID(UV);
		return H;
	};

	ImageGrayd bftHexBlending = bftHex.visualizeBlending(visualizationWidth, visualizationHeight, visualizationU, visualizationV);
	ImageRGBd bftHexTransition = bftHex.visualizeTiling(visualizationWidth, visualizationHeight, visualizationU, visualizationV);
	//IO::save01_in_u8(bftHexBlending, std::string("/home/nlutz/bftHexBlending.png"));
	//IO::save01_in_u8(bftHexTransition, std::string("/home/nlutz/bftHexTransition.png"));

	TandBFunctionType bftAlteratingSquares;
	bftAlteratingSquares.blendingFunction = [&] (double u, double v) -> float
	{
		u += arguments.uTranslation;
		v += arguments.vTranslation;
		int vInt = std::floor(v);
		if(vInt%2 == 1)
		{
			u += 0.5;
		}
		float sinU = std::sin(u*M_PI);
		float sinV = std::sin(v*M_PI);
		return sqrt(sinU*sinU * sinV*sinV);
	};
	bftAlteratingSquares.tilingFunction = [&] (double u, double v)
	{
		u += arguments.uTranslation;
		v += arguments.vTranslation;
		int vInt = std::floor(v);
		if(vInt%2 == 1)
		{
			u += 0.5;
		}
		Eigen::Vector2i vec;
		vec[0] = int(floor(u)*127);
		vec[1] = int(floor(v)*127);
		return vec;
	};

	ImageGrayd bftAltSquaresBlending = bftAlteratingSquares.visualizeBlending(visualizationWidth, visualizationHeight, visualizationU, visualizationV);
	ImageRGBd bftAltSquaresTransition = bftAlteratingSquares.visualizeTiling(visualizationWidth, visualizationHeight, visualizationU, visualizationV);
	//IO::save01_in_u8(bftAltSquaresBlending, std::string("/home/nlutz/bftAltSquaresBlending.png"));
	//IO::save01_in_u8(bftAltSquaresTransition, std::string("/home/nlutz/bftAltSquaresTransition.png"));

	TandBFunctionType bftRotSquare;
	bftRotSquare.blendingFunction = [&] (double u, double v) -> float
	{
		u += arguments.uTranslation;
		v += arguments.vTranslation;
		Eigen::Matrix2d rotate45; rotate45 << SQRT2DIV2, -SQRT2DIV2, SQRT2DIV2, SQRT2DIV2;
		//Eigen::Matrix2d scaleSqrt2; scaleSqrt2 << 1.0, 0, 0, 1.0;
		Eigen::Matrix2d scaleSqrt2; scaleSqrt2 << sqrt(2.0), 0, 0, sqrt(2);
		Eigen::Vector2d UV;
		UV[0] = u;
		UV[1] = v;
//		UV[0] += 0.25;
//		UV[1] += 0.25;
		UV = rotate45 * scaleSqrt2 * UV;
		UV[0] += 0.0;
		UV[1] += 0.0;
		float sinU = std::sin(UV[0]*M_PI);
		float sinV = std::sin(UV[1]*M_PI);
		return sqrt(sinU*sinU * sinV*sinV);
	};
	bftRotSquare.tilingFunction = [&] (double u, double v)
	{
		u += arguments.uTranslation;
		v += arguments.vTranslation;
		Eigen::Matrix2d rotate45; rotate45 << SQRT2DIV2, -SQRT2DIV2, SQRT2DIV2, SQRT2DIV2;
		//Eigen::Matrix2d scaleSqrt2; scaleSqrt2 << 1.0, 0, 0, 1.0;
		Eigen::Matrix2d scaleSqrt2; scaleSqrt2 << sqrt(2.0), 0, 0, sqrt(2);
		Eigen::Vector2d UV;
		UV[0] = u;
		UV[1] = v;
//		UV[0] += 0.25;
//		UV[1] += 0.25;
		UV = rotate45 * scaleSqrt2 * UV;
		UV[0] += 0.0;
		UV[1] += 0.0;
		Eigen::Vector2i vec;
		vec[0] = int(UV[0]+10)*63;
		vec[1] = int(UV[1]+10)*63;
		return vec;
	};

	ImageGrayd bftRotSquareBlending = bftRotSquare.visualizeBlending(visualizationWidth, visualizationHeight, visualizationU, visualizationV);
	ImageRGBd bftRotSquareTransition = bftRotSquare.visualizeTiling(visualizationWidth, visualizationHeight, visualizationU, visualizationV);
	//IO::save01_in_u8(bftRotSquareBlending, std::string("/home/nlutz/bftRotSquareBlending.png"));
	//IO::save01_in_u8(bftRotSquareTransition, std::string("/home/nlutz/bftRotSquareTransition.png"));

	TandBFunctionType bftAltRotSquare;
	bftAltRotSquare.blendingFunction = [&] (double u, double v) -> float
	{
		u += arguments.uTranslation;
		v += arguments.vTranslation;
		Eigen::Matrix2d rotate45; rotate45 << SQRT2DIV2, -SQRT2DIV2, SQRT2DIV2, SQRT2DIV2;
		//Eigen::Matrix2d scaleSqrt2; scaleSqrt2 << 1.0, 0, 0, 1.0;
		Eigen::Matrix2d scaleSqrt2; scaleSqrt2 << sqrt(2.0), 0, 0, sqrt(2);
		Eigen::Vector2d UV;
		UV[0] = u;
		UV[1] = v;
//		UV[0] += 0.25;
//		UV[1] += 0.25;
		UV = rotate45 * scaleSqrt2 * UV;
		UV[0] += 0.5;
		UV[1] += 0.5;
		float sinU = std::sin(UV[0]*M_PI);
		float sinV = std::sin(UV[1]*M_PI);
		return sqrt(sinU*sinU * sinV*sinV);
	};
	bftAltRotSquare.tilingFunction = [&] (double u, double v)
	{
		u += arguments.uTranslation;
		v += arguments.vTranslation;
		Eigen::Matrix2d rotate45; rotate45 << SQRT2DIV2, -SQRT2DIV2, SQRT2DIV2, SQRT2DIV2;
		//Eigen::Matrix2d scaleSqrt2; scaleSqrt2 << 1.0, 0, 0, 1.0;
		Eigen::Matrix2d scaleSqrt2; scaleSqrt2 << sqrt(2.0), 0, 0, sqrt(2);
		Eigen::Vector2d UV;
		UV[0] = u;
		UV[1] = v;
//		UV[0] += 0.25;
//		UV[1] += 0.25;
		UV = rotate45 * scaleSqrt2 * UV;
		UV[0] += 0.5;
		UV[1] += 0.5;
		Eigen::Vector2i vec;
		vec[0] = int(UV[0]+100)*256;
		vec[1] = int(UV[1]+100)*256;
		return vec;
	};

	ImageGrayd bftAltRotSquareBlending = bftAltRotSquare.visualizeBlending(visualizationWidth, visualizationHeight, visualizationU, visualizationV);
	ImageRGBd bftAltRotSquareTransition = bftAltRotSquare.visualizeTiling(visualizationWidth, visualizationHeight, visualizationU, visualizationV);
	//IO::save01_in_u8(bftAltRotSquareBlending, std::string("/home/nlutz/bftAltRotSquareBlending.png"));
	//IO::save01_in_u8(bftAltRotSquareTransition, std::string("/home/nlutz/bftAltRotSquareTransition.png"));

	TandBFunctionType::BlendingFunctionType blendingSingularityFunction = [&] (double u, double v) -> double
	{
		u += arguments.uTranslation;
		v += arguments.vTranslation;
		u -= 0.5f;
		v -= 0.5f;
		u*=1.41421356237;
		v*=1.41421356237;
		const float a = M_PI/4.0;
		float cosA = cos(a);
		float sinA = sin(a);
		Eigen::Matrix2d rotateA;
		rotateA << cosA, sinA,
				-sinA, cosA;
		Eigen::Vector2d V = rotateA*Eigen::Vector2d(u, v);
		float cosu = sin(V[0]*M_PI);
		float cosv = sin(V[1]*M_PI);
		return 0.05f*cosu*cosv*cosu*cosv;
	};

	RandomAffineTransform rAffT;
	if(arguments.isotropy == FULL_ISOTROPY)
	{
		rAffT.addAngle(RandomAffineTransform::AngleRangeType(std::make_pair(-M_PI, M_PI)));
	}
	else if(arguments.isotropy != UNSPECIFIED_ISOTROPY)
	{
		rAffT.addAngle(0);
		if(arguments.isotropy == PI_ISOTROPY || arguments.isotropy == HALF_PI_ISOTROPY)
		{
			rAffT.addAngle(M_PI);
		}
		if(arguments.isotropy == HALF_PI_ISOTROPY)
		{
			rAffT.addAngle(M_PI/2.0);
			rAffT.addAngle(3.0*M_PI / 2.0);
		}
	}
	else
	{
		//custom program
		rAffT.addAngle(std::make_pair(-M_PI/4.0, M_PI/4.0));
	}

//	BlendingFunctionType bftNeutral;
//	bftNeutral.blendingFunction = [&] (double u, double v)
//	{
//		return 1.0/* * cos(v*2.0)*/;
//	};
//	bftNeutral.transitionFunction = [&] (double u, double v)
//	{
//		Eigen::Vector2i vec;
//		return vec;
//	};

	std::string filename_cycles = std::string(argv[2]);
	std::string name_exemplar = IO::remove_path(argv[3]);
	std::string name_noext = IO::remove_ext(name_exemplar);
	ImageType im_out;
	TilingAndBlending<ImageType> perioBlend;
	for(int i=3; i<argc; ++i)
	{
		ImageType im_in, im_out;
		IO::loadu8_in_01(im_in, argv[i]);
		perioBlend.addTexture(im_in);
	}
	TilingAndBlending<ImageType>::MaskType *im_mask = nullptr;
	if(arguments.tryToUseMask)
	{
		im_mask = new TilingAndBlending<ImageType>::MaskType;
		if(argc > 3)
		{
			IO::loadu8_in_01(*im_mask, argv[4]);
		}
		perioBlend.setContentWeightingMask(im_mask);
		perioBlend.setWeightingMaskTest(arguments.contentMaskTest);
	}
	std::string suffix = name_noext + "_"
			+ std::to_string(arguments.width) + "x"
			+ std::to_string(arguments.height) +  "_uScale="
			+ std::to_string(arguments.uScale) + "_vScale="
			+ std::to_string(arguments.vScale) + (arguments.useHistogramTransfer ? "_HT_" : "_noHT_")
			+ (arguments.isotropy == NO_ISOTROPY ? "noIso" :
			   (arguments.isotropy == FULL_ISOTROPY ? "fullIso" :
			   (arguments.isotropy == PI_ISOTROPY ? "piIso" :
			   (arguments.isotropy == HALF_PI_ISOTROPY ? "halfPiIso" : "unspecIso"))))
			+ (arguments.tiling == SQUARE_SQUARE ? "_SqSq" :
				(arguments.tiling == ALTSQUARE_TRIANGULAR ? "_ASqTr" :
				(arguments.tiling == ROTSQUARE_ROTSQUARE ? "_rotSqSq" :
				(arguments.tiling == HEXAGONAL_TRIANGULAR ? "_hexTri" :
				(arguments.tiling == FOUR_SQUARES ? "_4Sq" : "_unspecTiling")))))
			+ (arguments.useSingularityHider ? "_SH" : "_noSH")
			+ "_uTr=" + std::to_string(arguments.uTranslation)
			+ "_vTr=" + std::to_string(arguments.vTranslation)
			+ (arguments.useAutocovarianceSampling ? (arguments.useCustomImportanceMap ? "_customSampling" : "_ACSampling") : "")
			+ (arguments.cheatAutocovariance && arguments.useAutocovarianceSampling ? "_cheatedAC" : "")
			+ (im_mask != nullptr ? "_withMask" : "")
			+ (im_mask != nullptr ? "_m" + std::to_string(arguments.contentMaskTest) : "")
			+".png";

	CycleMapType loadedCycles = loadCycles(filename_cycles);
	std::string textureName = IO::remove_ext(IO::remove_path(name_exemplar));
	CyclePair cyclePair;
	CycleMapType::const_iterator cit = loadedCycles.find(textureName);
	if(cit != loadedCycles.end())
		cyclePair = (*cit).second;

	switch (arguments.tiling)
	{
	case SQUARE_SQUARE:
		perioBlend.addFunction(bft1);
		perioBlend.addFunction(bft2);
		break;
	case ALTSQUARE_TRIANGULAR:
		perioBlend.addFunction(bftAlteratingSquares);
		perioBlend.addFunction(bftAltTri);
		break;
	case ROTSQUARE_ROTSQUARE:
		perioBlend.addFunction(bftRotSquare);
		perioBlend.addFunction(bftAltRotSquare);
		break;
	case HEXAGONAL_TRIANGULAR:
		perioBlend.addFunction(bftTri);
		perioBlend.addFunction(bftHex);
		break;
	case FOUR_SQUARES:
		perioBlend.addFunction(bft1);
		perioBlend.addFunction(bft2);
		perioBlend.addFunction(bft3);
		perioBlend.addFunction(bft4);
		break;
	}


	if(arguments.useSingularityHider)
	{
		perioBlend.setSingularityBlendingFunction(blendingSingularityFunction);
	}
	unsigned int outputWidth = arguments.width;
	unsigned int outputHeight = arguments.height;
	perioBlend.setWidth(arguments.width);
	perioBlend.setHeight(arguments.height);
	perioBlend.setUVScale(arguments.uScale, arguments.vScale);
	perioBlend.setUseHistogramTransfer(arguments.useHistogramTransfer);
	perioBlend.setRandomAffineTransform(rAffT);

	ImageGrayd bftSumBlending = perioBlend.visualizeBlendingSum(visualizationWidth, visualizationHeight, visualizationU, visualizationV);
	//IO::save01_in_u8(bftSumBlending, std::string("/home/nlutz/bftSumBlending.png"));

	perioBlend.setCSCycles(cyclePair.vectors[0], cyclePair.vectors[1]);
	perioBlend.setCSPolyphaseComponentSamplesNumber(Eigen::Vector2i(1, 1));

	unsigned int indexCustomMap = 4;
	if(arguments.useCustomImportanceMap && arguments.tryToUseMask)
	{
		indexCustomMap = 5;
	}
	ImageGrayd customMap;
	if(argc > indexCustomMap)
	{
		IO::loadu8_in_01(customMap, argv[indexCustomMap]);
	}
	else
	{
		arguments.useCustomImportanceMap = false;
		arguments.useAutocovarianceSampling = false;
	}

	ImageGrayd L, a, b, autocorrelation;
	if(arguments.isCyclostationary)
	{
		im_out = perioBlend.synthesize_cyclostationary();
		im_out.for_all_pixels([&] (ImageType::PixelType &pix)
		{
			for(unsigned i=0; i<3; ++i)
				pix[i] = std::max(std::min(1.0, pix[i]), 0.0);
		});
		IO::save01_in_u8(im_out, std::string("/home/nlutz/") + suffix);
	}
	else
	{
		if(arguments.useAutocovarianceSampling || arguments.compareAutocovariance)
		{
			bool rescaled = false;
			double widthScale = 0;
			double heightScale = 0;
			if(!arguments.useCustomImportanceMap)
			{
				ImageType im_in;
				IO::loadu8_in_01(im_in, argv[3]);
				if(im_in.width()>512 || im_in.height()>512)
				{ //I just want to scale an exemplar, why does it take like a million lines with itk
					rescaled = true;
					ImageType::ItkImg::SizeType inputSize  = im_in.itk()->GetLargestPossibleRegion().GetSize();
					ImageType::ItkImg::SizeType outputSize;
					outputSize[0] = std::min(512, im_in.width());
					outputSize[1] = std::min(512, im_in.height());
					using TransformType = itk::IdentityTransform <double, 2>;
					using ResampleImageFilterType = itk::ResampleImageFilter<ImageType::ItkImg, ImageType::ItkImg>;

					widthScale = (double(inputSize[0]) / double(outputSize[0]));
					heightScale = (double(inputSize[1]) / double(outputSize[1]));
					ImageType::ItkImg::SpacingType outputSpacing;
					outputSpacing[0] = im_in.itk()->GetSpacing()[0] * widthScale;
					outputSpacing[1] = im_in.itk()->GetSpacing()[1] * heightScale;

					ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();
					resample->SetInput(im_in.itk());
					resample->SetSize(outputSize);
					resample->SetOutputSpacing(outputSpacing);
					resample->SetTransform(TransformType::New());
					resample->UpdateOutputInformation();
					resample->Update();
					im_in.itk() = resample->GetOutput();
				}
				ImageType res;
				PCA<double> pca(im_in);
				MaskBool mb(im_in.width(), im_in.height());
				mb |= [] (int, int) {return true;};
				pca.computePCA(mb);
				pca.project(res);
				extract3Channels(res, L, a, b);

				//extractLab(im_in, L, a, b, true);
				IO::save01_in_u8(L, "/home/nlutz/L.png");
				if(!arguments.useCustomImportanceMap && !arguments.computeOnlyOutputAutocovariance)
				{
					Fourier::truePeriodicStationaryAutocovariance(L, autocorrelation, true, false);
					autocorrelation.for_all_pixels([&] (ImageGrayd::PixelType &pix)
					{
						if(arguments.cheatAutocovariance)
						{
							pix = pix > 0 ? 0 : -pix;
						}
						else
						{
							pix = pix < 0 ? 0 : pix;
						}
					});
					IO::save01_in_u8(autocorrelation, "/home/nlutz/autocorrelation_" + suffix);
				}
			}
			if(arguments.useAutocovarianceSampling)
			{
				Stamping::SamplerImportance *si = nullptr;
				if(arguments.useCustomImportanceMap)
				{
					if(arguments.importanceMapCutThreshold > 0)
					{
						customMap.for_all_pixels([&] (ImageGrayd::PixelType &pix, int x, int y)
						{
							pix = pix < arguments.importanceMapCutThreshold ? 0 : arguments.importanceMapCutThreshold;
						});
						IO::save01_in_u8(customMap, "/home/nlutz/adjustedImportanceMap_" + suffix);
					}
					si = new Stamping::SamplerImportance(customMap);
				}
				else
				{
					if(arguments.importanceMapCutThreshold > 0)
					{
						autocorrelation.for_all_pixels([&] (ImageGrayd::PixelType &pix, int x, int y)
						{
							pix = pix < arguments.importanceMapCutThreshold ? 0 : arguments.importanceMapCutThreshold;
						});
						IO::save01_in_u8(autocorrelation, "/home/nlutz/adjustedImportanceMap_" + suffix);
					}
					si = new Stamping::SamplerImportance(autocorrelation);
				}
				perioBlend.setImportanceSampler(si);
				im_out = perioBlend.synthesize_stationary();
				delete si;
				perioBlend.setImportanceSampler(nullptr);
			}
			else
			{
				im_out = perioBlend.synthesize_stationary();
			}
			im_out.for_all_pixels([&] (ImageType::PixelType &pix)
			{
				for(unsigned i=0; i<3; ++i)
					pix[i] = std::max(std::min(1.0, pix[i]), 0.0);
			});
			IO::save01_in_u8(im_out, std::string("/home/nlutz/") + suffix);
			if(arguments.compareAutocovariance)
			{
				if(rescaled)
				{ //I just want to scale an exemplar, why does it take like a million lines with itk
					using TransformType = itk::IdentityTransform <double, 2>;
					using ResampleImageFilterType = itk::ResampleImageFilter<ImageType::ItkImg, ImageType::ItkImg>;
					ImageType::ItkImg::SizeType outputSize;
					outputSize[0] = im_out.width() * widthScale;
					outputSize[1] = im_out.height() * heightScale;

					ImageType::ItkImg::SpacingType outputSpacing;
					outputSpacing[0] = im_out.itk()->GetSpacing()[0] * widthScale;
					outputSpacing[1] = im_out.itk()->GetSpacing()[1] * heightScale;

					ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();
					resample->SetInput(im_out.itk());
					resample->SetSize(outputSize);
					resample->SetOutputSpacing(outputSpacing);
					resample->SetTransform(TransformType::New());
					resample->UpdateOutputInformation();
					resample->Update();
					im_out.itk() = resample->GetOutput();
				}
				ImageGrayd autocorrelation_out;
				extractLab(im_out, L, a, b, true);

				Fourier::trueBoundedStationaryAutocovariance(L, autocorrelation_out, true, false, 512, 512);
				autocorrelation_out.for_all_pixels([&] (ImageGrayd::PixelType &pix)
				{
					pix = pix < 0 ? 0 : pix;
				});
				IO::save01_in_u8(autocorrelation_out, "/home/nlutz/autocorrelation_out_" + suffix);
				assert(autocorrelation.width() == autocorrelation_out.width());
				if(!arguments.computeOnlyOutputAutocovariance)
				{
					ImageGrayd norm2diff;
					norm2diff.initItk(autocorrelation.width(), autocorrelation.height());
					double totalDiff;
					norm2diff.for_all_pixels([&] (ImageGrayd::PixelType &pix, int x, int y)
					{
						ImageGrayd::PixelType &pixACIn = autocorrelation.pixelAbsolute(x, y);
						ImageGrayd::PixelType &pixACOut = autocorrelation_out.pixelAbsolute(x, y);
						pix = std::sqrt((pixACIn - pixACOut) * (pixACIn - pixACOut));
						totalDiff += pix;
					});
					std::cout << "autocovariance diff=" << totalDiff/(norm2diff.width()*norm2diff.height()) << std::endl;
					IO::save01_in_u8(norm2diff, "/home/nlutz/autocorrelationNorm2Diff" + suffix);
				}
			}
		}
		else
		{
			im_out = perioBlend.synthesize_stationary();
			im_out.for_all_pixels([&] (ImageType::PixelType &pix)
			{
				for(unsigned i=0; i<3; ++i)
					pix[i] = std::max(std::min(1.0, pix[i]), 0.0);
			});
			IO::save01_in_u8(im_out, std::string("/home/nlutz/") + suffix);
		}
	}


	//ImageRGBd im_visualizeSynthesis = perioBlend.visualizeSynthesis(outputWidth, outputHeight);

//	ImageGrayd bft1BlendingNormalized = divideBySum(bft1Blending, bftSumBlending);
//	ImageGrayd bft2BlendingNormalized = divideBySum(bft2Blending, bftSumBlending);

//	IO::save01_in_u8(bft1BlendingNormalized, std::string("/home/nlutz/bft1BlendingNormalized.png"));
//	IO::save01_in_u8(bft2BlendingNormalized, std::string("/home/nlutz/bft2BlendingNormalized.png"));

//	im_visualizeSynthesis.for_all_pixels([&] (ImageType::PixelType &pix)
//	{
//		for(unsigned i=0; i<3; ++i)
//			pix[i] = std::max(std::min(1.0, pix[i]), 0.0);
//	});

	//IO::save01_in_u8(im_visualizeSynthesis, std::string("/home/nlutz/visualization_") + suffix);
	IO::save01_in_u8(perioBlend.visualizePrimalAndDual(), std::string("/home/nlutz/blending_") + suffix);
	if(im_mask != nullptr)
		delete im_mask;
	return 0;
}
