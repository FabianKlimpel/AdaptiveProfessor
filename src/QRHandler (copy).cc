#include "Professor/QRHandler.h"
 
/**
 * This function initializes the object. Every necessary is getting calculated and the first iteration is performed.
 * @pts: Storage of the anchor points and the powers
 * @ptvals: Bin values of the anchor points
 */
void QRHandler::init(Professor::ParamPoints& pts, const vector<double>& ptvals){
	_iterationcounter = 0;
	_max = 0;
	_b = ptvals;
	_power = pts.getPower(0);
	//setting up the matrices and iterate
	increaseM(pts);
	makeMprime();
	initQR();
	iterateb();		
}

/**
 * This function performs an iteration of the matrix M and forwards the new entries to Q and R
 * @pts: Storage of the anchor points and the powers
 */
void QRHandler::iterateM(Professor::ParamPoints& pts){
	
	//if every power in @_power is already in use in @_m, new components of a higher power need to be calculated
	if(_m[0].size() == _power.size())
	{
		//@_max is the order of the polynomial function. It will be increased and the new powers will be calculated
		_max++;
		_power = pts.getPower(_max);
	}
	//iterate
	increaseM(pts);	
	_iterationcounter++;
	makeMprime();
	expandQR();
}

/**
 * This function updates @_bprime
 */
void QRHandler::iterateb(){
	_bprime = LinAlg::normalizeVec(_b);
	_bprime = LinAlg::multMatVec(LinAlg::transpose(_q), _bprime);
}

/**
 * This function function performs an iteration of all matrices and vectors
 * @pts: Storage of the anchor points and powers
 * @walkthrough: Flag that decides, if @_bprime needs to be calculated at every iteration
 */
void QRHandler::iterate(Professor::ParamPoints& pts, const bool walkthrough){
	iterateM(pts);
	//If the flag is set, @_bprime will not get calculated.
	if(!walkthrough)
		iterateb();
}

/**
 * This function initializes the matrices Q and R and sets the first values.
 */
const bool QRHandler::initQR(){
	//set the sizes of @_q and @_r
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

/**
 * This function updates the sizes of Q and R.
 */
void QRHandler::resizeQR(){
	//after the first iteration, the sizes are adapted according to the current size of @_mprime 
	_q.resize(_mprime.size());
	_r.resize(_mprime[0].size());
	for(size_t k = 0; k < _q.size(); k++)
		_q[k].resize(_mprime[0].size());
	for(size_t k = 0; k < _r.size(); k++)
		_r[k].resize(_mprime[0].size());
}

/**
 * This function builds the QR-decomposition or increase the matrices, if they already exist
 * @tmp: Temporary storage for the increase of @_r
 * @sum: Temporary storage for the increase of @_q
 */
void QRHandler::expandQR(){

	//if @_q and @_r are not initialize yet, it will be done here and the function quits
	if(initQR()) return;
	//update the sizes of @_q and @_r
	resizeQR();
	
	//the following calculations add the new elements to @_q and @_r
	double tmp;
	for(size_t l = 0; l < _iterationcounter; l++)
	{
		tmp = 0;
		for(size_t k = 0 ; k < _q.size(); k++)
			tmp += _q[k][l] * _mprime[k][_iterationcounter];
		_r[l][_iterationcounter] = tmp;
	}
	std::cout << "1" << std::endl;
	vector<double> sum;
	sum.assign(_q.size(), 0);
	for(size_t l = 0; l < sum.size(); l++)
		for(size_t k = 0; k < _iterationcounter; k++)
			sum[l] += _r[k][_iterationcounter] * _q[l][k];
	std::cout << "2" << std::endl;		
	for(size_t k = 0; k < _q.size(); k++)
		_q[k][_iterationcounter] = _mprime[k][_iterationcounter] - sum[k];
	std::cout << "3\t" << _iterationcounter << std::endl;
	for(size_t m = 0; m < _q.size(); m++){
		for(size_t m1 = 0; m1 < _q[m].size(); m1++)
			std::cout << _q[m][m1] << "\t";
		std::cout << std::endl;}
	
	_r[_iterationcounter][_iterationcounter] = LinAlg::getAbs(LinAlg::getCol(_q, _iterationcounter));
	std::cout << "4" << std::endl;
	tmp = LinAlg::getAbs(LinAlg::getCol(_q, _iterationcounter));
	for(size_t k = 0; k < _q.size(); k++)
	{
		_q[k][_iterationcounter] /= tmp;
		
		//If @tmp is 0, the result would become nan. Setting it instead to 0 "stabilizes" the calculations.
		if(std::isnan(_q[k][_iterationcounter])) 
			_q[k][_iterationcounter] = 0;
	}
		std::cout << "5" << std::endl;
}

/**
 * This function normalizes every column of @_m by using the absolut value of the respective column.
 * @abs: Absolut value of the column
 */
void QRHandler::makeMprime(){

	//setting @_mprime to the same size as @_m
	_mprime.resize(_m.size());
	for(size_t i = 0; i < _mprime.size(); i++)
		_mprime[i].resize(_m[i].size());

	//getting the absolut value of the column
	const double abs = LinAlg::getAbs(LinAlg::getCol(_m, _iterationcounter));

	//setting the normalization
	for(size_t i = 0; i < _m.size(); i++)
		_mprime[i][_iterationcounter] = _m[i][_iterationcounter] / abs;
}

/**
 * This function adds a new column to @_m
 * @pts: Storage of the anchor points and powers
 * @size: Old number of columns in @_m
 * @tmp: Temporary storage of the product of the parts of a monomial
 */
void QRHandler::increaseM(Professor::ParamPoints& pts){

	//set the size of @_m, if not done yet
	if(_m.empty())
		_m.resize(pts.numPoints());
	
	//adding a new column to @_m
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
			
		//the new terms will be added to the last column in every row and @tmp is resetted
		_m[j][size] = tmp;
		tmp = 1;
	}
}

/**
 * This function sets the member variables to a certain iteration
 * @order: Order of the polynomial function
 * @pts: Storage of the anchor points and powers
 * @numFitParams: Number of fit parameters
 */
void QRHandler::load(const int order, Professor::ParamPoints& pts, const size_t numFitParams){
	_power = pts.getPower(order);
	_max = order;
	_iterationcounter = 0;
	
	//set the size of @_m, if not done yet
	if(_m.empty())
		_m.resize(pts.numPoints());
	
	//set up @_m
	while(_m[0].size() < numFitParams){std::cout << "bsize: " << _bprime.size() << std::endl;
		iterate(pts, true); std::cout << "iteration\t" << _iterationcounter << "\t_msize: " << _m[0].size() << "\tnumfitparams: " << numFitParams << std::endl;}
	iterate(pts);
	std::cout << "bsize: " << _bprime.size() << std::endl;
}

/**
 * This function empties all vectors in this object
 */
void QRHandler::reset(){
	_m.clear();
	_r.clear();
	_q.clear();
	_mprime.clear();
}
