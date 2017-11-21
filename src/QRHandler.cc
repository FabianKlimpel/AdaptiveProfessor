#include "Professor/QRHandler.h"
 
void QRHandler::init(Professor::ParamPoints& pts, const vector<double>& ptvals){
	_iterationcounter = 0;
	_max = 0;
	_b = ptvals;
	_power = pts.getPower(0);
	//setting up the matrices and @_a
	increaseM(pts);
	makeMprime();
	initQR();
	iterateb();		
}

void QRHandler::iterateM(Professor::ParamPoints& pts){

	//if every power mentioned in @_power is already in use in @_m, new components of a higher power needs to be calculated
	if(_m[0].size() == _power.size())
	{
		//@_max is the order of the polynom. It will be increased and the new powers will be calculated
		_max++;
		_power = pts.getPower(_max);
	}

	//setting up the new member variables
	increaseM(pts);	
	_iterationcounter++;
	makeMprime();
	expandQR();
}

void QRHandler::iterateb(){
	_bprime = LinAlg::normalizeVec(_b);
	_bprime = LinAlg::multMatVec(LinAlg::transpose(_q), _bprime);
}

void QRHandler::iterate(Professor::ParamPoints& pts, const bool walkthrough){
	iterateM(pts);
	cout << "M calculated in iteration " << _iterationcounter;
	if(!walkthrough){
		iterateb(); cout << "\tb calculated in iteration " << _iterationcounter << endl;}
	else
		cout << endl;
}

const bool QRHandler::initQR(){
	if(_q.empty() || _r.empty())
	{
		_q.resize(_mprime.size());
		_r.resize(_mprime[0].size());
		
		for(size_t k = 0; k < _q.size(); k++)
			_q[k].resize(_mprime[0].size());
		
		for(size_t k = 0; k < _r.size(); k++)
			_r[k].resize(_mprime[0].size());
		
		//calculating the first elements by using http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5342430&tag=1 on top of eq. 5
		for(size_t k = 0; k < _q.size(); k++)
			_q[k][0] = LinAlg::getCol(_mprime,0)[k] / LinAlg::getAbs(LinAlg::getCol(_mprime,0));

		_r[0][0] = LinAlg::getAbs(LinAlg::getCol(_mprime,0));
		return true;
	}
	return false;
}

void QRHandler::resizeQR(){
	//in a later iteration, the size should be adapted according to the current size of @_mprime 
	_q.resize(_mprime.size());
	_r.resize(_mprime[0].size());
	for(size_t k = 0; k < _q.size(); k++)
		_q[k].resize(_mprime[0].size());
	for(size_t k = 0; k < _r.size(); k++)
		_r[k].resize(_mprime[0].size());
}

/**
 * This function builds the QR-decomposition or increase the matrices, if they already exist
 * @i: represents the iterationstep
 * @_q, @_r: the respective matrices of the QR-decomposition
 * @_mprime: coloumnwise rescaled @_m, c.f. @FitHandler::makemprime() 
 * @tmp: temporary storage for the increase of @_r
 * @sum: temporary storage for the increase of @_q
 */
void QRHandler::expandQR(){

	if(initQR()) return;

	resizeQR();
	
	//the following calculations add the new elements to @_q & @_r
	double tmp;
	for(size_t l = 0; l < _iterationcounter; l++)
	{
		tmp = 0;
		for(size_t k = 0 ; k < _q.size(); k++)
			tmp += _q[k][l] * _mprime[k][_iterationcounter];
		_r[l][_iterationcounter] = tmp;
	}

	vector<double> sum;
	sum.assign(_q.size(), 0);
	for(size_t l = 0; l < sum.size(); l++)
		for(size_t k = 0; k < _iterationcounter; k++)
			sum[l] += _r[k][_iterationcounter] * _q[l][k];
			
	for(size_t k = 0; k < _q.size(); k++)
		_q[k][_iterationcounter] = _mprime[k][_iterationcounter] - sum[k];

	_r[_iterationcounter][_iterationcounter] = LinAlg::getAbs(LinAlg::getCol(_q, _iterationcounter));

	tmp = LinAlg::getAbs(LinAlg::getCol(_q, _iterationcounter));
	for(size_t k = 0; k < _q.size(); k++)
	{
		_q[k][_iterationcounter] /= tmp;
		
		//If the absolut value of @q[k][i + 1] is 0, the result would become nan. Setting it instead to 0 "stabilizes" the calculations.
		if(std::isnan(_q[k][_iterationcounter])) 
			_q[k][_iterationcounter] = 0;
	}
}

/**
 * This function normalizes every coloumn of @_m by using the absolut value of the respective coloumn.
 * @coloumn: represents the coloumn of interest
 * @abs: absolut value of the coloumn
 */
void QRHandler::makeMprime(){

	//setting @_mprime to the same size as @_m
	_mprime.resize(_m.size());
	for(size_t i = 0; i < _mprime.size(); i++)
		_mprime[i].resize(_m[i].size());

	//setting the absolut value of the coloumn, because while assigning it row wise, it changes after every assigned value
	const double abs = LinAlg::getAbs(LinAlg::getCol(_m, _iterationcounter));
	
	//setting the normalization
	for(size_t i = 0; i < _m.size(); i++)
		_mprime[i][_iterationcounter] = _m[i][_iterationcounter] / abs;

}

/**
 * This function adds a new coloumn to @_m
 * @pts: container for the anchor points
 * @size: old number of coloumns in @_m
 * @tmp: temporary storage of the product of the different powers of the values
 */
void QRHandler::increaseM(Professor::ParamPoints& pts){

	//set the size of @_m, if not done yet
	if(_m.empty())
		_m.resize(pts.numPoints());
	
	//adding a new coloumn to @_m
	size_t size = _m[0].size();
	for(size_t i = 0; i < _m.size(); i++)
		_m[i].resize(size + 1);

	double tmp = 1;
	//walking over every row of @_m
	for(size_t j = 0; j < _m.size(); j++)
	{	
		//in every row the product of the values and the respective powers is calculated
		for(size_t i = 0; i < pts.dim(); i++)
			tmp *= pow(pts.pointScaled(j)[i], _power[_m[0].size() - 1][i]);
			
		//the new terms will be added to the last coloumn in every row and @tmp is resetted
		_m[j][size] = tmp;
		tmp = 1;
	}
}

void QRHandler::load(const int order, Professor::ParamPoints& pts, const size_t numFitParams){
	_power = pts.getPower(order);
	_max = order;
	
	//set the size of @_m, if not done yet
	if(_m.empty())
		_m.resize(pts.numPoints());
		
	//set up @_m
	while(_m[0].size() < numFitParams)
		increaseM(pts);
}

void QRHandler::reset(){
	_m.clear();
	_r.clear();
	_q.clear();
	_mprime.clear();
}
