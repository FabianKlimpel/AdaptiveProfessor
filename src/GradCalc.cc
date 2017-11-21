#include "Professor/GradCalc.h"

void GradCalc::initStructure(const Professor::ParamPoints& pts){
	if(_structure.empty())
	{
			_structure.resize(pts.numPoints());
			for(size_t j = 0; j < pts.numPoints(); j++)
				_structure[j].resize(pts.dim());
	}	
}

const double GradCalc::addMonomial(const size_t i, const size_t j, const size_t k, Professor::ParamPoints& pts, const int order){

	double tmp = 1;
	//extracting the respective parameter values
	for(size_t l = 0; l < pts.dim(); l++)
		//if the parameter is the one that is derived in the component of the gradient, then its derivative is used
		if(j == l)
			//if the value of the parameter is 0 & the power of the parameter to derive is != 1, the whole monomial will be 0 but the special case in the derivative of 0^0 = 1
			//if the power is 0, the whole monomial will be 0 after derived 
			if((pts.pointScaled(i)[l] == 0 && pts.getPower(order)[k][l] != 1) || pts.getPower(order)[k][l] == 0)
			{	
				tmp = 0;
				continue;
			}							
			else			
				//multiply the derived factor
				tmp *= pow(pts.pointScaled(i)[l], pts.getPower(order)[k][l] - 1) * pts.getPower(order)[k][l];	
		else
			//if at least one of the not derived parameters is 0 and its power is !=0, the whole monomial become 0
			if(pts.pointScaled(i)[l] == 0 && pts.getPower(order)[k][l] != 0)
			{
				tmp = 0;
				continue;
			}
			else
				tmp *= pow(pts.pointScaled(i)[l], pts.getPower(order)[k][l]);
	return tmp;
}

const double GradCalc::extendStructure(const size_t i, const size_t j, const size_t numFitParams, Professor::ParamPoints& pts, const int order){

	double tmp = 0;
	vector<double> tmpstruc = _structure[i][j];
	tmpstruc.resize(numFitParams);

	//start the calculation after the last calculated component and walk up to the new maximum needed
	for(size_t k = _structure[i][j].size(); k < tmpstruc.size(); k++)
		//put the new monomial part to the @tmpstruc at the specific point in the list
		tmpstruc[k] = addMonomial(i, j, k, pts, order);

	setStructure(i, j, tmpstruc);	

	return tmp;
}

const vector<double> GradCalc::getGradVector(const size_t i, const vector<double>& bfp, Professor::ParamPoints& pts, const int order){
	
	initStructure(pts);

	//set up the gradient and resize to the number of dimensions, initialized as 0
	vector<double> functiongradient;
	functiongradient.assign(_structure[0].size(), 0);

	//calculating the components of the gradient
	for(size_t j = 0; j < functiongradient.size(); j++)
		//ifenough monomials were already calculated, they can be taken directly
		if(_structure[i][j].size() > bfp.size()) 
			//walking over every monomial and calculate the contribution to the gradient vector
			for(size_t k = 1; k < bfp.size(); k++)
				functiongradient[j] += bfp[k] * _structure[i][j][k];
		else
		{
			extendStructure(i, j, bfp.size(), pts, order);
			//walking over every monomial, calculate what is already available
			for(size_t k = 1; k < _structure[i][j].size(); k++)
				functiongradient[j] += bfp[k] * _structure[i][j][k];
		}
	
	//normalizing the gradient
	double tmp = LinAlg::getAbs(functiongradient);
	
	//in order to avoid NaN's, the gradient is only normalized if it isn't a zerovector, else it's normalized vector is the vector itself
	if(tmp != 0)
		for(size_t i = 0; i < functiongradient.size(); i++)
			functiongradient[i] /= tmp;
		
	return functiongradient;
}

const vector<vector<double>> GradCalc::getAllGradVectors(const std::vector<double>& bfp, Professor::ParamPoints& pts, const int order){
	vector<vector<double>> result;
	for(size_t i = 0; i < _structure.size(); i++)
		result.push_back(getGradVector(i, bfp, pts, order));
	return result;	
}
