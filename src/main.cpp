// File:  main.cpp
// Date:  11/17/2019
// Auth:  K. Loux
// Desc:  Entry point for SplinePatchToFlatPattern application.

// Local headers


// Standard C++ headers
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <cassert>

struct Point
{
	Point() = default;
	Point(const double& xIn, const double& yIn, const double& zIn) : x(xIn), y(yIn), z(zIn) {}
	
	double x;
	double y;
	double z;
	
	Point operator-(const Point& p) const
	{
		Point r;
		r.x = x - p.x;
		r.y = y - p.y;
		r.z = z - p.z;
		return r;
	}
};

bool ParseToken(const std::string& token, double& value)
{
	std::istringstream ss(token);
	return !(ss >> value).fail();
}

bool ParseLine(const std::string& line, std::vector<Point>& curve1, std::vector<Point>& curve2)
{
	std::istringstream ss(line);
	Point p1, p2;
	std::string token;
	if (!std::getline(ss, token, ','))
		return false;
	
	if (!ParseToken(token, p1.x))
		return false;
		
	if (!std::getline(ss, token, ','))
		return false;
	
	if (!ParseToken(token, p1.y))
		return false;
		
	if (!std::getline(ss, token, ','))
		return false;
	
	if (!ParseToken(token, p1.z))
		return false;
		
	if (!std::getline(ss, token, ','))
		return false;
	
	if (!ParseToken(token, p2.x))
		return false;
		
	if (!std::getline(ss, token, ','))
		return false;
	
	if (!ParseToken(token, p2.y))
		return false;
		
	if (!std::getline(ss, token, ','))
		return false;
	
	if (!ParseToken(token, p2.z))
		return false;
		
	return true;
}

bool ReadInputFile(const std::string& fileName, std::vector<Point>& curve1, std::vector<Point>& curve2)
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
	void AddPoint(const Point& p, const Point& tangent);
	double GetLength() const;
	
private:
	std::vector<Point> points;
	std::vector<Point> tangents;
};

bool FitSplineToPoints(const std::vector<Point>& points, Spline& spline)
{
	const unsigned int splineControlPointCount(4);// Assume that we'll get a good fit if we choose four points.
	assert(points.size() > splineControlPointCount);
	
	spline.AddPoint(points.front(), Point(0.0, 1.0, 0.0));
	
	unsigned int i1(points.size() / (splineControlPointCount - 1));
	for (unsigned int a = 0; a < splineControlPointCount - 2; ++a)
		spline.AddPoint(points[i1 * a], points[i1 * a + 1] - points[i1 * a - 1]);
	
	spline.AddPoint(points.back(), Point(0.0, 1.0, 0.0));
	
	// TODO:  Best fit optimization?  variables are y-ordinate of end-point slopes and x, y, and z of slopes for all other points.
	
	return false;
}

bool GenerateFlatPattern(const Spline& s1, const Spline& s2, const unsigned int& divisions, std::vector<Point>& flatPatternPoints)
{
	return false;
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

	std::vector<Point> curve1, curve2;
	if (!ReadInputFile(argv[2], curve1, curve2))
		return 1;

	Spline spline1, spline2;
	if (!FitSplineToPoints(curve1, spline1) ||
		!FitSplineToPoints(curve2, spline2))
		return 1;

	const unsigned int divisions(100);
	std::vector<Point> flatPatternPoints;
	if (!GenerateFlatPattern(spline1, spline2, divisions, flatPatternPoints))
		return 1;
	
	return 0;
}
