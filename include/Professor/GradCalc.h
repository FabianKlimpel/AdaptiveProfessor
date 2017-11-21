#ifndef __GRADCALC__H
#define __GRADCALC__H

#include <vector>
#include "Professor/LinAlg.h"
#include "Professor/ParamPoints.h"

using namespace std;

namespace Professor{class ParamPoints;}

class GradCalc
{
public:
	GradCalc(){}

	void initStructure(const Professor::ParamPoints& pts);
	const vector<double> getGradVector(const size_t i, const std::vector<double>& bfp, Professor::ParamPoints& pts, const int order);
	const vector<vector<double>> getAllGradVectors(const std::vector<double>& bfp, Professor::ParamPoints& pts, const int order);

private:
	const double extendStructure(const size_t i, const size_t j, const size_t numFitParams, Professor::ParamPoints& pts, const int order);
	const double addMonomial(const size_t i, const size_t j, const size_t k, Professor::ParamPoints& pts, const int order);
	void setStructure(const size_t i, const size_t j, const vector<double> vec){_structure[i][j] = vec;}

	vector<vector<vector<double>>> _structure;
};
#endif
