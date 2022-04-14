#include <stdlib.h>
#include <ASTex/Perioblend/perioblend.h>
#include <ASTex/easy_io.h>

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
	FouriertexturePool[0]
	return 0;
}
