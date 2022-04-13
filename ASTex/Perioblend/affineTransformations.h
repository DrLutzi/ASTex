#ifndef _AFFINE_TRANSFORMATIONS_H_
#define _AFFINE_TRANSFORMATIONS_H_

#include <Eigen/Eigen>

namespace ASTex
{

/**
 * @brief The AffineTransform class is a warper for a 2x2 matrix.
 */
class AffineTransform
{
public:
	AffineTransform();
	AffineTransform(const Eigen::Matrix2d &transform);

	const Eigen::Matrix2d &transform() const;

private:
	Eigen::Matrix2d m_transform;
};

/**
 * @brief The RandomAffineTransform class manages random affine transforms through a random 2x2 matrix generator.
 */
class RandomAffineTransform
{
public:
	using AngleRangeType = std::pair<double, double>;
	using AngleRangeRandomType = std::pair<AngleRangeType, double>;
	RandomAffineTransform();

	/**
	 * @brief addAngle adds a rotation angle to the random transform pool.
	 * Don't forget to add 0 if you need it.
	 */
	void addAngle(double);

	/**
	 * @brief addAngle adds a rotation angle range to the random transform pool and enables ranges.
	 */
	void addAngle(AngleRangeType);

	/**
	 * @brief generate generates a random affine transform from the data provided
	 * @param ignoreWeight unimplemented feature
	 * @return an AffineTransformation containing the random affine transform
	 */
	AffineTransform generate(bool ignoreWeight = true, int seed = -1) const;

private:

	using AngleRangeRandomContainerType = std::vector<AngleRangeRandomType>;
	//Rotations
	bool m_rotationIsSample;
	AngleRangeRandomContainerType m_angleRanges;

	//weight consideration is unimplemented yet (the idea is that each range has a random chance of being picked at generation).
	double m_totalWeight;

	//I don't need other transformations do it yourself if you need them <3
};

};

#endif
