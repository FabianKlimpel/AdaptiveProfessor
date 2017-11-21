#ifndef __HYPERCUBEIPOL__H
#define __HYPERCUBEIPOL__H

#include <vector>
#include "Professor/LinAlg.h"
#include "Professor/ParamPoints.h"
#include "Professor/QRHandler.h"

using namespace std;

class HyperCubeIpol
{
public:
	//Constructor
	HyperCubeIpol();

	// Setter for the indices referring to a hypercube surrounding a point with index @i
	const vector< size_t>& gethypercube(const size_t i, const vector<vector<double>>& pts);

	// Getter of @_hypercubes
	const vector<vector<size_t>>& hypercubes() {return _hypercubes;}
	const vector<vector<double>>& getAllFitParams(const vector<vector<double>>& pts, const vector<double>& ptvals);
	const vector<double>& getFitParams(const size_t i, const vector<vector<double>>& pts, const vector<double>& ptvals);
private:
	void buildHyperCube(const size_t i, const vector<vector<double>>& pts);
	void calcFitParams(const size_t i, const vector<vector<double>>& pts, const vector<double>& ptvals);

	// list of indices for hypercubes of different points
	vector<vector<size_t>> _hypercubes;
	vector<vector<double>> _fitparams;
};

#endif

