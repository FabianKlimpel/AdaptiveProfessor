#ifndef __QRHANDLER__H
#define __QRHANDLER__H

#include <iostream>
#include <vector>
#include "Professor/LinAlg.h"
#include <limits>
#include "Professor/ParamPoints.h"
#include <cmath>

using namespace std;
namespace Professor{class ParamPoints;};

class QRHandler
{

public:
	QRHandler(){}

	void init(Professor::ParamPoints& pts, const vector<double>& ptvals);
	void iterate(Professor::ParamPoints& pts, bool walkthrough = false);
	void load(int order, Professor::ParamPoints& pts, const size_t numFitParams);
	void reset();

	const vector<vector<double>>& getM() const {return _m;}
	const vector<vector<double>>& getR() const {return _r;}
	const vector<double>& getBprime() const {return _bprime;}
	const vector<double>& getB() const {return _b;}
	const int getMaxPower() const {return _max;}
	const size_t getIterationCounter() const {return _iterationcounter;}
	const vector<vector<int>>& getPower() const {return _power;}

private:
	/**
	 * @_m: storage of the matrix M
	 * @_q, @r: storage of the matrices Q & R of the QR-decomposition
	 * @_d: identity matrix
	 * @_mprime: coloumnwise normalized M
	 * @_b: reference values
	 * @_a: vector for the fitparameters; is used as indicator of RR constrained parameters
	 * @_bprime: normalized @_b
	 * @_bfp, @_bfperr: best fit parameters and corresponding errors
	 * @_sigma: error of the reference values
	 * @_max: maximum of powers
	 * @_iterationcounter: number of iterations in the fitting process
	 * @_power: list that contains the powers of the variables at the terms of the fitting function
	 * @_dsmooth_current: Storage of the dsmooth value in a certain step in order to avoid recalculation. The default value indicates an unset value.
	 */
	vector<vector<double>> _m, _r, _q, _mprime;
	vector<double> _b, _bprime;
	vector<vector<int>> _power;
	int _max; 
	size_t _iterationcounter;

	const bool initQR();
	void iterateM(Professor::ParamPoints& pts);
	void iterateb();
	void resizeQR();

	//adding components to @_q & @_r
	void expandQR();

	//calculate @_mprime from @_m
	void makeMprime();

	//adding elements to @_m
	void increaseM(Professor::ParamPoints& pts);	
};

#endif

