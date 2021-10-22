#include <stdlib.h>
#include <vector>
#include <algorithm>
#include "ASTex/fourier.h"
#include "ASTex/utils.h"
#include "ASTex/easy_io.h"
#include "Eigen/Eigen"

template <typename Predicate> int FindInterval(int size,
		const Predicate &pred) {
	int first = 0, len = size;
	while (len > 0) {
		int half = len >> 1, middle = first + half;
		   if (pred(middle)) {
			   first = middle + 1;
			   len -= half + 1;
		   } else
			   len = half;

	}
	return ASTex::clamp_scalar<float, float>(first - 1, 0, size - 2);
}

struct Distribution1D {
	   Distribution1D(const float *f, int n)
		   : func(f, f + n), cdf(n + 1) {
			  cdf[0] = 0;
			  for (int i = 1; i < n + 1; ++i)
				  cdf[i] = cdf[i - 1] + func[i - 1] / n;

			  funcInt = cdf[n];
			  if (funcInt == 0) {
				  for (int i = 1; i < n + 1; ++i)
					  cdf[i] = float(i) / float(n);
			  } else {
				  for (int i = 1; i < n + 1; ++i)
					  cdf[i] /= funcInt;
			  }

	   }
	   int Count() const { return func.size(); }
	   float SampleContinuous(float u, float *pdf, int *off = nullptr) const {
			  int offset = FindInterval(cdf.size(),
				  [&](int index) { return cdf[index] <= u; });

		   if (off) *off = offset;
			  float du = u - cdf[offset];
			  if ((cdf[offset + 1] - cdf[offset]) > 0)
				  du /= (cdf[offset + 1] - cdf[offset]);
			  if (pdf) *pdf = func[offset] / funcInt;
			  return (offset + du) / Count();

	   }
	   int SampleDiscrete(float u, float *pdf = nullptr,
			   float *uRemapped = nullptr) const {
			  int offset = FindInterval(cdf.size(),
				  [&](int index) { return cdf[index] <= u; });

		   if (pdf) *pdf = func[offset] / (funcInt * Count());
		   if (uRemapped)
			   *uRemapped = (u - cdf[offset]) / (cdf[offset + 1] - cdf[offset]);
		   return offset;
	   }
	   float DiscretePDF(int index) const {
		   return func[index] / (funcInt * Count());
	   }

	   std::vector<float> func, cdf;
	   float funcInt;

};

class Distribution2D {
public:
	   Distribution2D(const float *data, int nu, int nv);
	   Eigen::Vector2f SampleContinuous(const Eigen::Vector2f &u, float *pdf) const {
		   float pdfs[2];
		   int v;
		   float d1 = pMarginal->SampleContinuous(u[1], &pdfs[1], &v);
		   float d0 = pConditionalV[v]->SampleContinuous(u[0], &pdfs[0]);
		   *pdf = pdfs[0] * pdfs[1];
		   return Eigen::Vector2f(d0, d1);
	   }
	   float Pdf(const Eigen::Vector2f &p) const {
		   int iu = ASTex::clamp_scalar(int(p[0] * pConditionalV[0]->Count()),
						  0, pConditionalV[0]->Count() - 1);
		   int iv = ASTex::clamp_scalar(int(p[1] * pMarginal->Count()),
						  0, pMarginal->Count() - 1);
		   return pConditionalV[iv]->func[iu] / pMarginal->funcInt;
	   }

private:
	   std::vector<std::unique_ptr<Distribution1D>> pConditionalV;
	   std::unique_ptr<Distribution1D> pMarginal;

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
