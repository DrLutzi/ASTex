#include <stdlib.h>
#include "ASTex/ContentExchange/atlas.h"
#include "ASTex/pca.h"
#include "ASTex/rpn_utils.h"

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0] << " <in_texture> <out_directory> [save rendering pack]" << std::endl;
		return EXIT_FAILURE;
	}

	ImageRGBu8 im_in;
	im_in.load(std::string(argv[1]));

	std::string out_dir = argv[2];
	std::string name_file = IO::remove_path(argv[1]);
	std::string name_noext = IO::remove_ext(name_file);
	create_directory(out_dir);

    ContentExchange::PatchProcessor<ImageRGBu8> pProcessor(im_in);
	pProcessor.setFilteringMode(NO_FILTER);
	pProcessor.setNbContentsPerPatch(5);
	pProcessor.setSeed(7);
	pProcessor.patches_initRandom(32);
	pProcessor.patches_dilate();
	pProcessor.contents_initDefault();
	pProcessor.contents_initRandom();
	//pProcessor.fullProcess_oldMethod();

	Stamping::SamplerImportance *si = nullptr;
	if(argc>3)
	{
		ImageGrayd importanceMap;
		IO::loadu8_in_01(importanceMap, argv[3]);
		si = new Stamping::SamplerImportance(importanceMap);
		pProcessor.setImportanceSampler(si);
	}

	pProcessor.debug_savePatchMap(out_dir + "/" + name_noext + "_patchMap" + ".png");

	if(argc>3 && std::atoi(argv[3])!=0)
	{
		create_directory(out_dir);
		std::string renderingPackDir = out_dir + "/" + name_noext + "_renderingPack";
		create_directory(renderingPackDir);
		pProcessor.saveRenderingPack(renderingPackDir);
	}
	pProcessor.setOutputSize(2*im_in.width(), 1*im_in.height());
	ImageRGBu8 output = pProcessor.generate().texture();
	output.save(out_dir + "/" + name_noext + "_contentExchange" + ".png");

//	ImageRGBd outputDouble;
//	outputDouble.initItk(output.width(), output.height());
//	outputDouble.for_all_pixels([&] (ImageRGBd::PixelType &pix, int x, int y)
//	{
//		for(int i=0; i<3; ++i)
//		{
//			pix[i] = output.pixelAbsolute(x, y)[i]/255.0;
//		}
//	});

//	//delete this
//	ImageRGBd res;
//	PCA<double> pca(outputDouble);
//	MaskBool mb(outputDouble.width(), outputDouble.height());
//	mb |= [] (int, int) {return true;};
//	pca.computePCA(mb);
//	pca.project(res);
//	ImageGrayd L, a, b, autocorrelation_out;
//	extract3Channels(res, L, a, b);

//	Fourier::trueBoundedStationaryAutocovariance(L, autocorrelation_out, true, false, 512, 512);
//	autocorrelation_out.for_all_pixels([&] (ImageGrayd::PixelType &pix)
//	{
//			pix = pix > 1.0 ? 1.0 : (pix < 0.0 ? 0.0 : pix);
//	});
//	IO::save01_in_u8(autocorrelation_out, "/home/nlutz/autocorrelation_out_contex.png");

	if(si)
		delete si;
	return 0;
}
