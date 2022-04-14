#include <stdlib.h>
#include <ASTex/Perioblend/perioblend.h>
#include <ASTex/easy_io.h>

typedef enum {NO_ISOTROPY=0, FULL_ISOTROPY=1, PI_ISOTROPY=2, HALF_PI_ISOTROPY=3, UNSPECIFIED_ISOTROPY=4} Isotropy_enum_t;

typedef struct
{
	bool useHistogramTransfer;
	unsigned int width;
	unsigned int height;
	Isotropy_enum_t isotropy;
	double uScale;
	double vScale;
	bool useSingularityHider;
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
			}
		}
	}
	return arguments;
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
		return std::min(std::min(u - std::floor(u), std::ceil(u) - u), std::min(v - std::floor(v), std::ceil(v) - v));
	};
	bft1.blendingFunction = [&] (double u, double v)
	{
		float sinU = std::sin(u*M_PI);
		float sinV = std::sin(v*M_PI);
		return sqrt(sinU*sinU * sinV*sinV);
	};
	bft1.tilingFunction = [&] (double u, double v)
	{
		Eigen::Vector2i vec;
		vec[0] = int(u)*63;
		vec[1] = int(v)*63;
		return vec;
	};
	ImageGrayd bft1Blending = bft1.visualizeBlending(visualizationWidth, visualizationHeight, visualizationU, visualizationV);
	ImageRGBd bft1Transition = bft1.visualizeTiling(visualizationWidth, visualizationHeight, visualizationU, visualizationV);
	IO::save01_in_u8(bft1Blending, std::string("/home/nlutz/bft1Blending.png"));
	IO::save01_in_u8(bft1Transition, std::string("/home/nlutz/bft1Transition.png"));

	TandBFunctionType bft2;
	bft2.blendingFunction = [&] (double u, double v)
	{
		u += 0.5;
		v += 0.5;
		return std::min(std::min(u - std::floor(u), std::ceil(u) - u), std::min(v - std::floor(v), std::ceil(v) - v));
	};
	bft2.blendingFunction = [&] (double u, double v)
	{
		u += 0.5f;
		v += 0.5f;
		float sinU = std::sin(u*M_PI);
		float sinV = std::sin(v*M_PI);
		return sqrt(sinU*sinU * sinV*sinV);
	};
	bft2.tilingFunction = [&] (double u, double v)
	{
		Eigen::Vector2i vec;
		vec[0] = int(u + 0.5)*127;
		vec[1] = int(v + 0.5)*127;
		return vec;
	};
	ImageGrayd bft2Blending = bft2.visualizeBlending(visualizationWidth, visualizationHeight, visualizationU, visualizationV);
	ImageRGBd bft2Transition = bft2.visualizeTiling(visualizationWidth, visualizationHeight, visualizationU, visualizationV);
	IO::save01_in_u8(bft2Blending, std::string("/home/nlutz/bft2Blending.png"));
	IO::save01_in_u8(bft2Transition, std::string("/home/nlutz/bft2Transition.png"));

	TandBFunctionType::BlendingFunctionType blendingSingularityFunction = [&] (double u, double v) -> double
	{
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

	std::string name_file = IO::remove_path(argv[2]);
	std::string name_noext = IO::remove_ext(name_file);
	ImageType im_out;
	TilingAndBlending<ImageType> perioBlend;
	for(int i=2; i<argc; ++i)
	{
		ImageType im_in, im_out;
		IO::loadu8_in_01(im_in, argv[i]);
		perioBlend.addTexture(im_in);
	}

	perioBlend.addFunction(bft1);
	perioBlend.addFunction(bft2);
	perioBlend.setSingularityBlendingFunction(blendingSingularityFunction);
	unsigned int outputWidth = 2048;
	unsigned int outputHeight = 1024;
	perioBlend.setWidth(arguments.width);
	perioBlend.setHeight(arguments.height);
	perioBlend.setUVScale(arguments.uScale, arguments.vScale);
	perioBlend.setUseHistogramTransfer(arguments.useHistogramTransfer);
	perioBlend.setRandomAffineTransform(rAffT);

	ImageGrayd bftSumBlending = perioBlend.visualizeBlendingSum(visualizationWidth, visualizationHeight, visualizationU, visualizationV);
	IO::save01_in_u8(bftSumBlending, std::string("/home/nlutz/bftSumBlending.png"));
	im_out = perioBlend.synthesize_stationary();
	ImageRGBd im_visualizeSynthesis = perioBlend.visualizeSynthesis(outputWidth, outputHeight);

	ImageGrayd bft1BlendingNormalized = divideBySum(bft1Blending, bftSumBlending);
	ImageGrayd bft2BlendingNormalized = divideBySum(bft2Blending, bftSumBlending);

	IO::save01_in_u8(bft1BlendingNormalized, std::string("/home/nlutz/bft1BlendingNormalized.png"));
	IO::save01_in_u8(bft2BlendingNormalized, std::string("/home/nlutz/bft2BlendingNormalized.png"));

	im_visualizeSynthesis.for_all_pixels([&] (ImageType::PixelType &pix)
	{
		for(unsigned i=0; i<3; ++i)
			pix[i] = std::max(std::min(1.0, pix[i]), 0.0);
	});

	im_out.for_all_pixels([&] (ImageType::PixelType &pix)
	{
		for(unsigned i=0; i<3; ++i)
			pix[i] = std::max(std::min(1.0, pix[i]), 0.0);
	});

	IO::save01_in_u8(im_out, std::string("/home/nlutz/") + name_noext + "_"
					 + std::to_string(arguments.width) + "x"
					 + std::to_string(arguments.height) +  "_uScale="
					 + std::to_string(arguments.uScale) + "_vScale="
					 + std::to_string(arguments.vScale) + (arguments.useHistogramTransfer ? "_HT_" : "_noHT_")
					 + (arguments.isotropy == NO_ISOTROPY ? "noIso" :
						(arguments.isotropy == FULL_ISOTROPY ? "fullIso" :
						(arguments.isotropy == PI_ISOTROPY ? "piIso" :
						(arguments.isotropy == HALF_PI_ISOTROPY ? "halfPiIso" : "unspecifiedIso")))) + ".png");
	IO::save01_in_u8(im_visualizeSynthesis, std::string("/home/nlutz/visualization_out.png"));
	return 0;
}
