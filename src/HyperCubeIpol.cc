#include "Professor/HyperCubeIpol.h"

vector<double> rescale(const size_t i, const vector<vector<double>>& pts){
	
	vector<double> result;
	double min = numeric_limits<double>::infinity(), max = -numeric_limits<double>::infinity();
	
	for(size_t j = 0; j < pts[0].size(); j++)
	{
		for(size_t k = 0; k < pts.size(); k++)
		{
			if(pts[k][j] < min)
				min = pts[k][j];
			if(pts[k][j] > max)
				max = pts[k][j];
		}
		result.push_back((pts[i][j] - min) / (max - min));
		min = numeric_limits<double>::infinity();
		max = -numeric_limits<double>::infinity();
	}
	return result;	
}

HyperCubeIpol::HyperCubeIpol(){}

void HyperCubeIpol::buildHyperCube(const size_t i, const vector<vector<double>>& pts){
  
	//setting up @center and the resultvector
	vector<double> center = pts[i];
	vector<size_t> result;
	
	//if a point isn't set, it's value is #number of anchor points + 1, so that these points can be found and won't be used for further calculations
	//this is important for points at the border of the sampleregion
	result.assign(pow(2., pts[0].size()), pts.size());
	
	//The same procedure as before is done for the distances, so that a check is possible if the smallest possible hypercube cann be constructed. Therefore it's initial values are set to infinity.
	vector<double> distances;
	distances.assign(pow(2., pts[0].size()), numeric_limits<double>::infinity());	
	size_t bin = 0;
	
	for(size_t j = 0; j < pts.size(); j++)
		//if a datapoint matches the @center, it will be skipped
		if(i == j)
			continue;
		else
		{
			//component wise check, if the components are bigger than the center point
			for(size_t k = 0; k < pts[0].size(); k++)
				if(center[k] >= pts[j][k])
					//All 'if' checks were binary questions, therefore a binary representation as summary can be found.
					//@bin is a value that can directly connect to the overall situation of the test point due to its value
					bin += pow(2., k);

			//check up, if the new trial point is closer than the stored one in the category specified by @bin
			if(LinAlg::getDistanceOfVectors(center, pts[j]) < distances[bin])
			{
				//if closer than set it as new point
				distances[bin] = LinAlg::getDistanceOfVectors(center, pts[j]);
				result[bin] = j;				
			}	
			bin = 0;
		}
	//store the result
	_hypercubes[i] = result;
	_hypercubes[i].push_back(i);
}
  
/**
 * This function calculates a list of indices referring to the hypercube surrounding the point with index @i
 * @i: number of the point around which a hypercube should be constructed
 * @center: temporary storage of the coordinates of point @i
 * @result: storage of the hypercube indices
 * @distances: list of the shortest distances found in every direction
 * @bin: this parameter is used as an indicator for the direction of the proposal point in a binary coding
 */
const vector<size_t>& HyperCubeIpol::gethypercube(const size_t i, const vector<vector<double>>& pts){
	
	//If the size of the hypercubes does not fit, it will be resized.
	//This construction is meant as initializer of @_hypercubes
	if(_hypercubes.size() != pts.size())
		_hypercubes.resize(pts.size());
		
	//if the hypercube is already calculated, return it
	if(!_hypercubes[i].empty())
		return _hypercubes[i];
	
	buildHyperCube(i, pts);	
	return _hypercubes[i];
}

const vector<vector<double>>& HyperCubeIpol::getAllFitParams(const vector<vector<double>>& pts, const vector<double>& ptvals){
	if(_fitparams.empty())
		_fitparams.resize(_hypercubes.size());
		
	for(size_t i = 0; i < _fitparams.size(); i++)
		if(_fitparams[i].empty())
			calcFitParams(i, pts, ptvals);
	return _fitparams;
}

const vector<double>& HyperCubeIpol::getFitParams(const size_t i, const vector<vector<double>>& pts, const vector<double>& ptvals){

	if(_hypercubes.empty())
		for(size_t j = 0; j < pts.size(); j++)
			gethypercube(j, pts);
	if(_fitparams.empty())
		_fitparams.resize(_hypercubes.size());

	if(_fitparams[i].empty())
		calcFitParams(i, pts, ptvals);

	return _fitparams[i];
}

void HyperCubeIpol::calcFitParams(const size_t i, const vector<vector<double>>& pts, const vector<double>& ptvals){
	
	vector<vector<double>> tmp_points, tmp_points_scaled;
	vector<double> tmp_ptvals;
	QRHandler qrh;
	for(size_t j = 0; j < gethypercube(i, pts).size(); j++)
		if(gethypercube(i, pts)[j] < pts.size())
		{
			tmp_points.push_back(pts[gethypercube(i, pts)[j]]);
			tmp_points_scaled.push_back(rescale(gethypercube(i, pts)[j], pts));
			tmp_ptvals.push_back(ptvals[gethypercube(i, pts)[j]]);
		}
	Professor::ParamPoints tmp_pts(tmp_points);	
	tmp_pts.setpointsScaled(tmp_points_scaled);
	
	qrh.init(tmp_pts, tmp_ptvals);
	for(size_t j = 0; j < pts[0].size(); j++)
		qrh.iterate(tmp_pts);
	
	vector<double> a;
	a.resize(pts[0].size() + 1);
	LinAlg::collinearity(a, qrh.getR(), qrh.getB(), 10e-10, 10e-10);
		
	_fitparams[i] = LinAlg::getBestFitParameters(qrh.getR(), a, qrh.getBprime());
		
	for(size_t j = 0; j < _fitparams[i].size(); j++)
		_fitparams[i][j] *= LinAlg::getAbs(qrh.getB()) / LinAlg::getAbs(LinAlg::getCol(qrh.getM(), j));
}











