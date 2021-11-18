#ifndef __SAMPLER__H__
#define __SAMPLER__H__

#include <vector>
#include <iostream>
#include <fstream>
#include <memory.h>

#include <cmath>
#include <Eigen/Core>
#include <random>
#include <ctime>
#include "ASTex/image_gray.h"
#include "ASTex/utils.h"

namespace ASTex
{

namespace Stamping
{

//TODO: add iterator for samplers

/**
 * @brief The SamplerBase class is an interface for every Samplers.
 * the main idea is for them to provide an array of coordinates between (0, 0) and (1, 1),
 * and setters to modify the parameters of the generative process.
 */
class SamplerBase
{
public:

	/**
	 * @brief SamplerBase constructor for SamplerBase.
	 */
	SamplerBase() {}
	virtual ~SamplerBase() {}

	/**
	 * @brief generate a function which yields an array of floating point coordinates.
	 * @return the coordinates array.
	 */
	virtual std::vector<Eigen::Vector2f> generate()=0;
};

/**
 * @brief The SamplerOrigin class is an override of SamplerBase,
 * supposed to yield an array of points all located at the origin.
 */
class SamplerOrigin : public SamplerBase
{
public:
	SamplerOrigin(unsigned nbPoints=1):
		SamplerBase(),
		m_nbPoints(nbPoints)
	{}

	void setNbPoints(unsigned nbPoints) {m_nbPoints = nbPoints;}

	std::vector<Eigen::Vector2f> generate();

private:
	unsigned int m_nbPoints;
};

/**
 * @brief The SamplerRegular class is an override of SamplerBase,
 * supposed to yield an array of points sampled on a regular grid.
 * Behavior: The first point is always at (0,0),
 * but the last point will be before (M, N). (could change?)
 */
class SamplerRegular : public SamplerBase
{
public:
	/**
	 * @brief SamplerRegular constructor for SamplerRegular.
	 * @param nbPoints is the default number of points the generate() function yields.
	 */
	SamplerRegular(unsigned nbPoints=0) :
		SamplerBase(),
		m_nbPoints(nbPoints)
	{}

	/**
	 * @param nbPoints is the number of points per dimension the generate() function will produce.
	 * For example, giving 3 will produce 9 points.
	 */
	void setNbPoints(unsigned nbPoints) {m_nbPoints = nbPoints;}

	/**
	 * @brief generate yields an array of floating point coordinates disposed on a regular grid.
	 * @return the coordinates array.
	 */
	std::vector<Eigen::Vector2f> generate();

private:
	unsigned int m_nbPoints;
};

/**
 * @brief The SamplerUniform class is an override of SamplerBase,
 * supposed to yield an array of points in respect to a uniform random process.
 */
class SamplerUniform : public SamplerBase
{
public:
	/**
	 * @brief SamplerUniform constructor for SamplerUniform.
	 * @param nbPoints is the default number of points the generate() function yields.
	 */
	SamplerUniform(unsigned nbPoints=0) :
		SamplerBase(),
		m_nbPoints(nbPoints)
	{}

	/**
	 * @param nbPoints is the total number of points the generate() function will produce.
	 */
	void setNbPoints(unsigned nbPoints) {m_nbPoints = nbPoints;}

	/**
	 * @brief generate yields an array of floating point coordinates
	 * randomly distributed according to 2 uniform laws between 0 and 1.
	 * @return the coordinates array.
	 */
	std::vector<Eigen::Vector2f> generate();

private:
	unsigned int m_nbPoints;
};

/**
 * @brief The SamplerPoisson class is an override of SamplerBase,
 * supposed to yield an array of points in respect to a Poisson random process.
 */
class SamplerPoissonGrid: public SamplerBase
{
public:
	/**
	 * @brief SamplerPoisson constructor for SamplerPoisson.
	 * @param nbPoints the default number of points the generate() function yields.
	 */
	SamplerPoissonGrid(unsigned nbPoints=0) :
		SamplerBase(),
		m_nbPoints(nbPoints),
		m_generator(),
		m_newPointsCount(30),
		m_generateInCircle(false),
		m_minDistance(-1.0f)
	{}

	/**
	 * @param nbPoints is the total number of points the generate() function will produce.
	 */
	void setNbPoints(unsigned nbPoints)     {m_nbPoints = nbPoints;}
	/**
	 * @param count defines the number of points the process tries to pick in a single loop, 30 by default.
	 * Refer to bridson-siggraph07-poissondisk.pdf for details (the value 'k')
	 */
	void setNewPointsCount(int count)       {m_newPointsCount = std::max(0, count);}
	/**
	 * @param b true for generating the poisson process within a circle of diameter 1.
	 * false for generating the process within a unity square.
	 */
	void setGenerateInCircle(bool b)        {m_generateInCircle = b;}
	/**
	 * @param d the minimal distance between points.
	 * If under 0, will be determined according to the number of points.
	 */
	void setMinDistance(float d)            {m_minDistance = d;}

	/**
	 * @brief generate yields an array of floating point coordinates
	 * randomly distributed according to a poisson law between (0,0) and (1,1).
	 * @return the coordinates array.
	 */
	std::vector<Eigen::Vector2f> generate();

private:
	/**
	 * @brief The DefaultPRNG class a pseudo random number generator class
	 * based on the std PRNG classes.
	 * Note for developpement: we should include a PRNG in every other sampler class.
	 */
	class DefaultPRNG
	{
	public:
		DefaultPRNG()
		: m_Gen( std::random_device()() )
		, m_Dis( 0.0f, 1.0f )
		{
			// prepare PRNG
			m_Gen.seed( time( nullptr ) );
		}

		explicit DefaultPRNG( uint32_t seed )
		: m_Gen( seed )
		, m_Dis( 0.0f, 1.0f )
		{
		}

		float RandomFloat()
		{
			return static_cast<float>( m_Dis( m_Gen ) );
		}

		int RandomInt( int Max )
		{
			std::uniform_int_distribution<> DisInt( 0, Max );
			return DisInt( m_Gen );
		}

	private:
		std::mt19937 m_Gen;
		std::uniform_real_distribution<float> m_Dis;
	};

	struct sGridPoint
	{
		sGridPoint( int X, int Y )
			: x( X )
			, y( Y )
		{}
		int x;
		int y;
	};


	struct sGrid
	{
		sGrid( int W, int H, float CellSize )
			: m_W( W )
			, m_H( H )
			, m_CellSize( CellSize )
		{
			m_Grid.resize( m_H );

			for ( auto i = m_Grid.begin(); i != m_Grid.end(); i++ ) { i->resize( m_W ); }
		}
		sGridPoint ImageToGrid( const Eigen::Vector2f& P, float CellSize )
		{
			return sGridPoint( ( int )( P[0] / CellSize ), ( int )( P[1] / CellSize ) );
		}
		void Insert( const Eigen::Vector2f& P )
		{
			sGridPoint G = ImageToGrid( P, m_CellSize );
			m_Grid[ G.x ][ G.y ] = P;
		}
		bool IsInNeighbourhood( Eigen::Vector2f Point, float MinDist, float CellSize )
		{
			sGridPoint G = ImageToGrid( Point, CellSize );

			// number of adjucent cells to look for neighbour points
			const int D = 5;

			// scan the neighbourhood of the point in the grid
			for ( int i = G.x - D; i < G.x + D; i++ )
			{
				for ( int j = G.y - D; j < G.y + D; j++ )
				{
					if ( i >= 0 && i < m_W && j >= 0 && j < m_H )
					{
						Eigen::Vector2f P = m_Grid[ i ][ j ];

						Eigen::Vector2f np = P-Point;
						float X = np[0];
						float Y = np[1];
						float d = std::sqrt( X*X + Y*Y );

						if ( d < MinDist ) { return true; }
					}
				}
			}


			return false;
		}

	private:
		int m_W;
		int m_H;
		float m_CellSize;

		std::vector< std::vector<Eigen::Vector2f> > m_Grid;
	};

	bool IsInRectangle(Eigen::Vector2f p) const
	{
		return p[0] >= 0 && p[1] >= 0 && p[0] <= 1 && p[1] <= 1;
	}
	//
	bool IsInCircle(Eigen::Vector2f p) const
	{
		float fx = p[0] - 0.5f;
		float fy = p[1] - 0.5f;
		return ( fx*fx + fy*fy ) <= 0.25f;
	}

	Eigen::Vector2f PopRandom( std::vector<Eigen::Vector2f>& Points, DefaultPRNG& Generator )
	{
		const int Idx = Generator.RandomInt( Points.size()-1 );
		const Eigen::Vector2f P = Points[ Idx ];
		Points.erase( Points.begin() + Idx );
		return P;
	}

	Eigen::Vector2f GenerateRandomPointAround( const Eigen::Vector2f& P, float MinDist, DefaultPRNG& Generator )
	{
		// start with non-uniform distribution
		float R1 = Generator.RandomFloat();
		float R2 = Generator.RandomFloat();

		// radius should be between MinDist and 2 * MinDist
		float Radius = MinDist * ( R1 + 1.0f );

		// random angle
		float Angle = 2 * 3.141592653589f * R2;

		// the new point is generated around the point (x, y)
		float X = P[0] + Radius * cos( Angle );
		float Y = P[1] + Radius * sin( Angle );

		return Eigen::Vector2f( X, Y );
	}

	unsigned int m_nbPoints;
	DefaultPRNG m_generator;
	int m_newPointsCount;
	int m_generateInCircle;
	float m_minDistance;
};

/**
 * @brief The SamplerCycles class is an override of SamplerBase,
 * supposed to yield an array of points in respect to a uniform random process.
 */
class SamplerCycles : public SamplerBase
{
public:
	/**
	 * @brief SamplerCycles constructor for SamplerUniform.
	 * @param nbPoints is the default number of points the generate() function yields.
	 */
	SamplerCycles(unsigned nbPoints=0) :
		SamplerBase(),
		m_nbPoints(nbPoints),
		m_cycles()
	{}

	void setCycles(const Eigen::Vector2f &firstCycle, const Eigen::Vector2f &secondCycle)
	{
		m_cycles[0] = firstCycle;
		m_cycles[1] = secondCycle;
	}

	/**
	 * @param nbPoints is the total number of points the generate() function will produce.
	 */
	void setNbPoints(unsigned nbPoints) {m_nbPoints = nbPoints;}

	/**
	 * @brief generate yields an array of floating point coordinates
	 * randomly distributed according to 2 uniform laws between 0 and 1.
	 * @return the coordinates array.
	 */
	std::vector<Eigen::Vector2f> generate();

private:
	unsigned int		m_nbPoints;
	Eigen::Vector2f		m_cycles[2];
};

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

/**
 * @brief The SamplerCycles class is an override of SamplerBase,
 * supposed to yield an array of points in respect to a uniform random process.
 */
class SamplerImportance : public SamplerBase
{
public:
	/**
	 * @brief SamplerCycles constructor for SamplerUniform.
	 * @param nbPoints is the default number of points the generate() function yields.
	 */
	SamplerImportance(ImageGrayd importanceFunction, unsigned nbPoints=0);

	~SamplerImportance()
	{
		delete m_functionFloatPointer;
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

		Distribution2D() {}

		void init(const float *func, int nu, int nv)
		{
			for (int v = 0; v < nv; ++v)
			{
				pConditionalV.emplace_back(new Distribution1D(&func[v * nu], nu));
			}
			std::vector<float> marginalFunc;
			for (int v = 0; v < nv; ++v)
				marginalFunc.push_back(pConditionalV[v]->funcInt);
			pMarginal.reset(new Distribution1D(&marginalFunc[0], nv));
		}


		Eigen::Vector2f SampleContinuous(const Eigen::Vector2f &u, float *pdf) const
		{
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

	/**
	 * @param nbPoints is the total number of points the generate() function will produce.
	 */
	void setNbPoints(unsigned nbPoints) {m_nbPoints = nbPoints;}

	/**
	 * @brief generate yields an array of floating point coordinates
	 * randomly distributed according to 2 uniform laws between 0 and 1.
	 * @return the coordinates array.
	 */
	std::vector<Eigen::Vector2f> generate();

private:
	ImageGrayd			m_importanceFunction;
	unsigned int		m_nbPoints;
	Distribution2D		m_distribution2D;
	float *				m_functionFloatPointer;
};

} //namespace Stamping

} //namespace ASTex


#endif //__SAMPLER__H__
