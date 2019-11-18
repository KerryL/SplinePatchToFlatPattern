// File:  main.cpp
// Date:  11/17/2019
// Auth:  K. Loux
// Desc:  Entry point for SplinePatchToFlatPattern application.

// optimization headers
#include "optimization/nelderMead.h"

// Eigen headers
#include <Eigen/Eigen>

// Standard C++ headers
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <cassert>

bool ParseToken(const std::string& token, double& value)
{
	std::istringstream ss(token);
	return !(ss >> value).fail();
}

bool ParseLine(const std::string& line, std::vector<Eigen::Vector3d>& curve1, std::vector<Eigen::Vector3d>& curve2)
{
	std::istringstream ss(line);
	Eigen::Vector3d p1, p2;
	std::string token;
	if (!std::getline(ss, token, ','))
		return false;
	
	if (!ParseToken(token, p1(0)))
		return false;
		
	if (!std::getline(ss, token, ','))
		return false;
	
	if (!ParseToken(token, p1(1)))
		return false;
	p1(1) = fabs(p1(1));
		
	if (!std::getline(ss, token, ','))
		return false;
	
	if (!ParseToken(token, p1(2)))
		return false;
		
	if (!std::getline(ss, token, ','))
		return false;
	
	if (!ParseToken(token, p2(0)))
		return false;
		
	if (!std::getline(ss, token, ','))
		return false;
	
	if (!ParseToken(token, p2(1)))
		return false;
	p2(1) = fabs(p2(1));
		
	if (!std::getline(ss, token, ','))
		return false;
	
	if (!ParseToken(token, p2(2)))
		return false;

	curve1.push_back(p1);
	curve2.push_back(p2);
		
	return true;
}

bool ReadInputFile(const std::string& fileName, std::vector<Eigen::Vector3d>& curve1, std::vector<Eigen::Vector3d>& curve2)
{
	std::ifstream file(fileName);
	if (!file.is_open() || !file.good())
	{
		std::cerr << "Failed to open '" << fileName << "' for input\n";
		return false;
	}
	
	std::string line;
	unsigned int lineCount(0);
	while (std::getline(file, line))
	{
		++lineCount;
		if (!ParseLine(line, curve1, curve2))
		{
			std::cerr << "Failed to parse line " << line << '\n';
			return false;
		}
	}
	
	return true;
}

class Spline
{
public:
	void AddPoint(const Eigen::Vector3d& p, const Eigen::Vector3d& c)
	{
		intersectionPoints.push_back(p);
		controlVectors.push_back(c);
	}

	void SetControlVector(const unsigned int& i, const Eigen::Vector3d& v) { controlVectors[i] = v; }

	unsigned int GetSegmentCount() const { return intersectionPoints.size() - 1; }
	Eigen::Vector3d GetIntersectionPoint(const unsigned int& i) const { return intersectionPoints[i]; }
	Eigen::Vector3d GetControlVector(const unsigned int& i) const { return controlVectors[i]; }
	
private:
	std::vector<Eigen::Vector3d> intersectionPoints;
	std::vector<Eigen::Vector3d> controlVectors;
};

std::vector<Eigen::Vector3d> ComputeSpline(const Spline& s, const unsigned int& segmentResolution)
{
	const auto segments(s.GetSegmentCount());
	std::vector<Eigen::Vector3d> points(segments * segmentResolution);
	for (unsigned int i = 0; i < segments; ++i)
	{
		double t(0.0);
		const double tStep(1.0 / segmentResolution);
		for (unsigned int j = 0; j < segmentResolution; ++j)
		{
			const Eigen::Vector3d p0(s.GetIntersectionPoint(i));
			const Eigen::Vector3d p1([&i, &s]() -> Eigen::Vector3d
			{
				if (i == 0)
					return s.GetIntersectionPoint(i) + s.GetControlVector(i);
				return s.GetIntersectionPoint(i) - s.GetControlVector(i);
			}());
			const Eigen::Vector3d p2(s.GetIntersectionPoint(i + 1) + s.GetControlVector(i + 1));
			const Eigen::Vector3d p3(s.GetIntersectionPoint(i + 1));

			points[i * segmentResolution + j] = pow(1.0 - t, 3) * p0 + 3.0 * pow(1.0 - t, 2) * t * p1 + 3.0 * (1. - t) * t * t * p2 + pow(t, 3) * p3;
			t += tStep;
		}
	}

	return points;
}

double ComputeError(const Spline& s, const std::vector<Eigen::Vector3d>& goalPoints)
{
	const unsigned int resolution(1000);
	const auto sPoints(ComputeSpline(s, resolution));
	double e(0.0);
	for (const auto& p : sPoints)
	{
		double minDistance(std::numeric_limits<double>::max());
		for (const auto& v : goalPoints)
		{
			const auto distance((p - v).norm());
			if (distance < minDistance)
				minDistance = distance;
		}
		e += minDistance;
	}

	return e;
}

std::vector<Eigen::Vector3d> BuildControlVectors(const Eigen::VectorXd& x)
{
	std::vector<Eigen::Vector3d> controlVectors;
	controlVectors.push_back(Eigen::Vector3d(0.0, fabs(x(0)), 0.0));

	int i;
	for (i = 1; i < x.size() - 1; i += 3)
		controlVectors.push_back(Eigen::Vector3d(x(i), x(i + 1), x(i + 2)));

	controlVectors.push_back(Eigen::Vector3d(0.0, fabs(x(i)), 0.0));

	return controlVectors;
}

struct SplineFitArgs : public Optimizer::AdditionalArgs
{
	SplineFitArgs(const std::vector<Eigen::Vector3d>& goalPoints,
		const std::vector<Eigen::Vector3d>& intersectionPoints) : goalPoints(goalPoints), intersectionPoints(intersectionPoints) {}

	const std::vector<Eigen::Vector3d>& goalPoints;
	const std::vector<Eigen::Vector3d>& intersectionPoints;
};

Eigen::VectorXd DoIteration(const Eigen::VectorXd& guess, const Optimizer::AdditionalArgs* args)
{
	const auto& arguments(*dynamic_cast<const SplineFitArgs*>(args));
	const auto controlVectors(BuildControlVectors(guess));
	Spline s;
	for (unsigned int i = 0; i < controlVectors.size(); ++i)
		s.AddPoint(arguments.intersectionPoints[i], controlVectors[i]);

	return Eigen::VectorXd(guess.size()).setOnes() * ComputeError(s, arguments.goalPoints);
}

void FitSplineToPoints(const std::vector<Eigen::Vector3d>& points, Spline& spline)
{
	constexpr unsigned int splineSegmentCount(3);// Assume that we'll get a good fit if we choose three segments.
	assert(points.size() > splineSegmentCount);
	Eigen::VectorXd initialGuess((splineSegmentCount - 1) * 3 + 2, 1);
	initialGuess.setOnes();
	
	std::vector<Eigen::Vector3d> intersectionPoints;
	spline.AddPoint(points.front(), Eigen::Vector3d(0.0, 1.0, 0.0));
	intersectionPoints.push_back(points.front());

	const unsigned int i1(points.size() / splineSegmentCount);
	for (unsigned int a = 1; a < splineSegmentCount; ++a)
	{
		spline.AddPoint(points[i1 * a], points[i1 * a - 1] - points[i1 * a + 1]);
		intersectionPoints.push_back(points[i1 * a]);
		initialGuess((a - 1) * 3 + 1) = spline.GetControlVector(a)(0);
		initialGuess((a - 1) * 3 + 2) = spline.GetControlVector(a)(1);
		initialGuess((a - 1) * 3 + 3) = spline.GetControlVector(a)(2);
	}

	spline.AddPoint(points.back(), Eigen::Vector3d(0.0, 1.0, 0.0));
	intersectionPoints.push_back(points.back());
	
	SplineFitArgs arguments(points, intersectionPoints);
	const unsigned int iterationLimit(10000);
	NelderMead<(splineSegmentCount - 1) * 3 + 2> optimizer(DoIteration, iterationLimit, &arguments);
	optimizer.SetInitialGuess(initialGuess * 4.0);
	const auto x(optimizer.Optimize());
	const auto newControlVectors(BuildControlVectors(x));

	for (unsigned int i = 0; i < newControlVectors.size(); ++i)
		spline.SetControlVector(i, newControlVectors[i]);
}

double ComputeLength(const std::vector<Eigen::Vector3d>& p)
{
	double length(0.0);
	for (unsigned int i = 1; i < p.size(); ++i)
		length += (p[i] - p[i - 1]).norm();
	return length;
}

bool FindIntersectionOfTwoCircles(const Eigen::Vector2d& c1, const double& r1,
	const Eigen::Vector2d& c2, const double& r2, Eigen::Vector2d& isect1, Eigen::Vector2d& isect2)
{
	const double distance((c1 - c2).norm());
	if (distance > r1 + r2 || distance < fabs(r1 -r2) || (distance == 0.0 && r1 == r2))// If there are no solutions, or infinite solutions, we cannot proceed
		return false;

	const double a((r1 * r1 - r2 * r2 + distance * distance) / (2.0 * distance));
	const double h(sqrt(r1 * r1 - a * a));
	const Eigen::Vector2d p(c1 + a * (c2 - c1) / distance);

	isect1(0) = p(0) + h * (c2(1) - c1(1)) / distance;
	isect1(1) = p(0) - h * (c2(0) - c1(0)) / distance;

	isect2(0) = p(0) - h * (c2(1) - c1(1)) / distance;
	isect2(1) = p(0) + h * (c2(0) - c1(0)) / distance;

	return true;
}

Eigen::Vector2d ChooseBestIntersection(const Eigen::Vector2d& isect1, const Eigen::Vector2d& isect2, const std::vector<Eigen::Vector2d>& c)
{
	if (c.size() < 2)
		return isect1;

	if ((c.back() - isect1).norm() > (c.back() - isect2).norm())
		return isect1;
	return isect2;
}

bool GenerateFlatPattern(const Spline& s1, const Spline& s2, const double& stepTarget, std::vector<Eigen::Vector2d>& flatPatternPoints)
{
	const unsigned int resolution(1000);
	const auto c1(ComputeSpline(s1, resolution));
	const auto c2(ComputeSpline(s2, resolution));

	const double s1Length(ComputeLength(c1));
	const double s2Length(ComputeLength(c2));

	const double step1(s1Length > s2Length ? stepTarget : stepTarget * s1Length / s2Length);
	const double step2(s2Length > s1Length ? stepTarget : stepTarget * s2Length / s1Length);

	std::vector<Eigen::Vector2d> curve1, curve2;
	double d((c1.front() - c2.front()).norm());
	curve1.push_back(Eigen::Vector2d(0.0, 0.0));
	curve2.push_back(Eigen::Vector2d(d, 0.0));

	unsigned int i1(1), i2(1);
	unsigned int i1Last(0), i2Last(0);
	while (i1Last < c1.size() && i2Last < c2.size())
	{
		for (; i1 < curve1.size() - 1; ++i1)
		{
			if ((c1[i1] - c1[i1Last]).norm() > step1)
				break;
		}

		for (; i2 < curve2.size() - 1; ++i2)
		{
			if ((c2[i2] - c2[i2Last]).norm() > step2)
				break;
		}

		const double d1From1((c1[i1] - c1[i1Last]).norm());
		const double d1From2((c1[i1] - c2[i2Last]).norm());
		const double d2From1((c2[i2] - c1[i1Last]).norm());
		const double d2From2((c2[i2] - c2[i2Last]).norm());
		i1Last = i1;
		i2Last = i2;

		Eigen::Vector2d isect11, isect12, isect21, isect22;
		if (!FindIntersectionOfTwoCircles(curve1.back(), d1From1, curve2.back(), d1From2, isect11, isect12))
			return false;
		if (!FindIntersectionOfTwoCircles(curve1.back(), d2From1, curve2.back(), d2From2, isect21, isect22))
			return false;

		curve1.push_back(ChooseBestIntersection(isect11, isect12, curve1));
		curve2.push_back(ChooseBestIntersection(isect21, isect22, curve2));
	}

	flatPatternPoints = curve1;
	std::reverse(curve2.begin(), curve2.end());
	flatPatternPoints.insert(flatPatternPoints.end(), curve2.begin(), curve2.end());

	return true;
}

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		std::cout << "Usage:  " << argv[0] << " <input file>\n"
			<< "  Input file must be comma-delimited and must contain four columns.\n"
			<< "  The first three columns are (x,y,z) for a series of points\n"
			<< "  describing one spline, and columns 4-5 are (x,y,z) for a series of\n"
			<< "  points describing the second spline.  Points should only be included\n"
			<< "  for half of each curve (i.e. positive y-ordinates only).  It is\n"
			<< "  assumed that x-z plane symmetry is desired, and curves are\n"
			<< "  constrained to have slopes parallel to the y-axis where the curves\n"
			<< "  meet the x-z plane\n" << std::endl;
	}

	std::vector<Eigen::Vector3d> curve1, curve2;
	if (!ReadInputFile(argv[1], curve1, curve2))
		return 1;

	Spline spline1, spline2;
	FitSplineToPoints(curve1, spline1);
	FitSplineToPoints(curve2, spline2);

	/*const unsigned int res(30);
	const auto c1(ComputeSpline(spline1, res));
	const auto c2(ComputeSpline(spline2, res));
	std::ofstream o("out.csv");
	for (unsigned int i = 0; i < c1.size(); ++i)
		o << c1[i](0) << ',' << c1[i](1) << ',' << c1[i](2) << ',' << c2[i](0) << ',' << c2[i](1) << ',' << c2[i](2) << '\n';//*/

	const double distanceResolution(0.01);
	std::vector<Eigen::Vector2d> flatPatternPoints;
	if (!GenerateFlatPattern(spline1, spline2, distanceResolution, flatPatternPoints))
		return 1;

	std::ofstream o("flatPattern.csv");
	o.precision(10);
	for (const auto& p : flatPatternPoints)
		o << std::fixed << p(0) << ',' << p(1) << '\n';
	
	return 0;
}
