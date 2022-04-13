#include "affineTransformations.h"
#include <random>

using namespace ASTex;

AffineTransform::AffineTransform() :
	m_transform()
{}

AffineTransform::AffineTransform(const Eigen::Matrix2d &transform) :
	m_transform(transform)
{}

const Eigen::Matrix2d &AffineTransform::transform() const
{
	return m_transform;
}



RandomAffineTransform::RandomAffineTransform() :
	m_rotationIsSample(true),
	m_angleRanges(),
	m_totalWeight(0.0)
{}

void RandomAffineTransform::addAngle(double angle)
{
	m_angleRanges.push_back(std::make_pair(std::make_pair(angle, angle), 1.0));
}

void RandomAffineTransform::addAngle(AngleRangeType angleRange)
{
	m_angleRanges.push_back(std::make_pair(angleRange, 1.0));
	m_rotationIsSample = false;
}

AffineTransform RandomAffineTransform::generate(bool ignoreWeight, int seed) const
{
	if(m_angleRanges.size() == 0)
	{
		return AffineTransform(Eigen::Matrix2d::Identity());
	}
	(void) ignoreWeight;
	//Sadly no C++17 for me yet
	static std::random_device rd;
	static std::mt19937 gen;
	if(seed != -1)
	{
		gen.seed(seed);
	}
	else
	{
		gen.seed(rd());
	}
	std::default_random_engine engine;
	AngleRangeRandomContainerType::const_iterator cit = m_angleRanges.begin();
	std::uniform_int_distribution<> dis(0, std::distance(cit, m_angleRanges.end()) - 1);
	std::advance(cit, dis(gen));
	const AngleRangeRandomType &arrt = *cit;
	const AngleRangeType &art = arrt.first;
	double angle;
	if(m_rotationIsSample)
	{
		angle = art.first;
	}
	else
	{
		std::uniform_real_distribution<double> unifd(art.first, art.second);
		angle = unifd(gen);
	}
	Eigen::Rotation2D<double> rot2D(angle);
	return AffineTransform(rot2D.matrix());
}
