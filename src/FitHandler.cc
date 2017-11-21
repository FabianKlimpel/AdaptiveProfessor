#include "Professor/FitHandler.h"

/**
 * Constructor
 */
FitHandler::FitHandler(){}

/**
 * Main constructor.
 * @fitparams: vector containing the fitparameters
 * @pterrs: vector containing the uncertainties of the data
 * @pts: Object containing everything about the anchor points
 * @order: order of the polynomial function
 * @num_ipol: bin number
 */
FitHandler::FitHandler(const vector<double>& fitparams, const vector<double>& pterrs, Professor::ParamPoints& pts, const int order, const int num_ipol){

	//local storage of parameters
	_bfp = fitparams;
	_sigma = pterrs;
	_qrh.load(order, pts, getNumFitParams());

	//calculate the fit errors
	setFitErrors(num_ipol);
}

/**
 * This function rescales the best fit parameters due to normalization of @_m and @_b
 * @vec: vector, containing the best fit parameters
 */
void FitHandler::rescaleBestFitParameters(){

	//backwards calculating of the norms that influenzed @_m and @_b in order to get the regular best fit parameters
	for(size_t i = 0; i < _bfp.size(); i++)
		_bfp[i] *= LinAlg::getAbs(_qrh.getB()) / LinAlg::getAbs(LinAlg::getCol(_qrh.getM(), i));
}

/**
 * This function calculates the Chi^2 value of the fit.
 * @result: resulting Chi^2
 * @functionvalue: functionvalue of the fit at a certain point
 */
const double FitHandler::getChi2() const{

	double result = 0, functionvalue = 0;
	cout << "sizes: " << _qrh.getM().size() << "\t" << _qrh.getB().size() << endl;
	//walk over every polynomial
	for(size_t i = 0; i < _qrh.getM().size(); i++)
	{
		//if the bin entry is 0, it can be a simple misscalculation in Rivet etc.
		if(_qrh.getB()[i] == 0)
			continue;
		cout << "sizes2: " << _qrh.getM()[i].size() << "\t" << _bfp.size() << endl;
		for(size_t j = 0; j < _qrh.getM()[i].size(); j++)
			functionvalue += _qrh.getM()[i][j] * _bfp[j];
			
		//get the difference between the functionvalue and the datapoint		
		functionvalue -= _qrh.getB()[i];
		cout << "functionvalue: " << functionvalue << endl;
		//if the uncertainty is 0, the calculation would rise a nan
		//this is prevented by setting it to the arbitrary value of 1e-10
		if(_sigma[i] == 0)
			result += (functionvalue * functionvalue) / (1e-10 * 1e-10);
		else
			result += (functionvalue * functionvalue) / (_sigma[i] * _sigma[i]);

		cout << "result: " << result << endl;
		functionvalue = 0;
	}
	return result;
}

/**
 * This function returns Dsmooth
 * @pts: container for the anchor points
 * @gh: container of the gradient vectors of the anchor points
 * @result: result of the calculation that will be returned
 */
const double FitHandler::getDsmooth(Professor::ParamPoints& pts){

	//if the smoothness was not calculated in the iteration ...
	if(_dsmooth_current == 2)
	{
		double result = 0;
		//walking over all anchor points and adding the dotproduct of the gradient vectors
		for(size_t i = 0; i < pts.numPoints(); i++)
			result += LinAlg::dotProduct(pts.getGradient(i), _gc.getGradVector(i, _bfp, pts, getMaxPower()));
		_dsmooth_current = result / (double) pts.numPoints();
	}
	//returning the mean
	return _dsmooth_current;
}

/**
 * This function calculates the uncertainty of the fit parameters.
 * @num_ipol: number of the bin
 * @mat: storage of the precision matrix, later of the covariance matrix
 * @tmp: temporary storage of precision matrix elements
 */
void FitHandler::setFitErrors(const size_t num_ipol){

	//setting up variables
	MatrixXd mat(_bfp.size(), _bfp.size());
	double tmp;
	//fill the matrix; walking over the upper triangle
	for(size_t row = 0; row < _bfp.size(); row++)
		for(size_t col = row; col < _bfp.size(); col++)
		{	
			//pull an element
			tmp = getInvCovMatElement(row, col);
			if(row == col)
				//if the element is on the diagonal then set it
				mat(row, col) = tmp;
			else
			{
				//using the symmetry, setting off-diagonal elements twice
				mat(row, col) = tmp;
				mat(col, row) = tmp;
			}
		}

	//calculating the inverse of the matrix
	mat = mat.fullPivHouseholderQr().inverse();

	//If an element is +/-infinity, it will be set to the maximum of double instead. That way one keeps regular numerical eleements
	for(size_t i = 0; i < _bfp.size(); i++)
		if(mat(i, i) == std::numeric_limits<double>::infinity())
			_bfperr.push_back(std::numeric_limits<double>::max());
		else
			if(-mat(i, i) == std::numeric_limits<double>::infinity())
				_bfperr.push_back(-std::numeric_limits<double>::max());
			else
				_bfperr.push_back(sqrt(mat(i, i)));
		
	//writing the covariance matrix to file
	OutputHandler oh;
	oh.writeCovMat(mat, num_ipol);
}

/**
 * This function calculates an element of the precision matrix
 * @i, @j: row/coloumn of the matrix
 * @resulg: numerical value of the matrix element
 */
const double FitHandler::getInvCovMatElement(const size_t i, const size_t j) const{
	double result = 0;

	//walk over the fit function
	for(size_t k = 0; k < _qrh.getM().size(); k++)
	{ 
		//If the uncertainty is zero, the result would be inf. Those elements will be skipped. It does not change the overall result, since it will be left out for every element.
		if(_sigma[k] == 0)
			continue;		
		//bin wise calculating d^2\chi^2/da_ida_j with the fit parameters \vec{a}
		result += _qrh.getM()[k][i] * _qrh.getM()[k][j] / (_sigma[k] * _sigma[k]);	
	}
	return result;
}

/**
 * This function calculates the next iterationstep in a bin
 * @pts: container of the anchor points
 * @threshold: threshold for the RR constraint of the fit
 * @kappa: shift in case of applied RR constraint
 */
void FitHandler::nextStep(Professor::ParamPoints& pts, const double threshold, const double kappa){
	
	_gc.initStructure(pts);
	//reset the smoothness, so that will be calculated again
	_dsmooth_current = 2;
	
	_qrh.iterate(pts);
	_a.resize(_qrh.getM()[0].size());
	LinAlg::collinearity(_a, _qrh.getR(), _qrh.getB(), threshold, kappa);

	cout << "_a: ";
	for(size_t i = 0; i < _a.size(); i++)
		cout << _a[i] << "\t";
	cout << endl;
	
	//getting the fit and rescale it
	_bfp = LinAlg::getBestFitParameters(_qrh.getR(), _a, _qrh.getBprime());	
	
	cout << "_bfp: ";
	for(size_t i = 0; i < _bfp.size(); i++)
		cout << _bfp[i] << "\t";
	cout << endl;
	
	rescaleBestFitParameters();
	
	cout << "rescaled _bfp: ";
	for(size_t i = 0; i < _bfp.size(); i++)
		cout << _bfp[i] << "\t";
	cout << endl;
}

/**
 * This function sets up the object for a new bin
 * @pts: container of the anchor points
 * @ptvals: bin values for each anchor point
 * @pterrs: error of the bin values for each anchor point
 * @num_ipol: number of the bin
 * @threshold: threshold for the RR constraint of the fit
 * @kappa: shift in case of applied RR constraint
 */
void FitHandler::startBin(Professor::ParamPoints& pts, const vector<double>& ptvals, const vector<double>& pterrs, const double threshold, const double kappa){

	_sigma = pterrs;
	_qrh.init(pts, ptvals);
	_a.resize(_qrh.getM()[0].size());
	
	LinAlg::collinearity(_a, _qrh.getR(), _qrh.getB(), threshold, kappa);
	
	//getting the fit and rescaling it
	_bfp = LinAlg::getBestFitParameters(_qrh.getR(), _a, _qrh.getBprime());
	rescaleBestFitParameters();	
}

/**
 * This function is a getter for the size of @_bfp
 * @zeros: counter of 0's as fit parameters
 */
const size_t FitHandler::getNumFitParams() const{
	size_t zeros = 0;

	//count the number of 0's at the end of the fit parameters list
	for(size_t i = _bfp.size(); i > 0; i--)
		if(_bfp[i - 1] == 0)
			zeros++;
		else
			break;
		
	//return the number of fit parameters that are != 0	
	return _bfp.size() - zeros;
}

/**
 * This function sets the object up for a certain bin and calculates a given number of iterations
 * @pts: container of the anchor points
 * @ptvals: bin values for each anchor point
 * @pterrs: error of the bin values for each anchor point
 * @num_ipol: number of the bin
 * @threshold: threshold for the RR constraint of the fit
 * @kappa: shift in case of applied RR constraint
 * @bestiteration: number of iterations that will be performed
 */
void FitHandler::setToIteration(Professor::ParamPoints& pts, const vector<double>& ptvals, const vector<double>& pterrs, const int num_ipol, const double threshold, const double kappa, const int bestiteration){

	_bfp.clear();
	_bfperr.clear();
	_qrh.reset();
	_d.clear();
	_a.clear();
	//setting up the bin
	startBin(pts, ptvals, pterrs, threshold, kappa);
	if(bestiteration > 0)
	{
		//Fast walking the number of iterations along.
		for(size_t i = 0; i < bestiteration - 1; i++) 
			_qrh.iterate(pts, true);
		
		//Perform the whole calculation only for the last step. That way the Object is indifferent to a regular stepping.
		nextStep(pts, threshold, kappa);
	}	
}

/**
 * This function calculates the dot product of all normalized gradients at every anchor point with every other
 * @pts: container for the anchor points
 * @gh: container of the gradient vectors of the anchor points
 * @result: vector containing all dot products
 * @tmp: temporary storage of the gradient vectors
 */
const vector<double> FitHandler::getAllGradDotProducts(Professor::ParamPoints& pts){
	
	vector<double> result;
	vector<vector<double>> tmp;
	tmp.resize(pts.numPoints());
	
	//walk over every possible combination of anchor points
	for(size_t i = 0; i < pts.numPoints(); i++)
	{
		//if a gradient was not calculated yet, do it
		if(tmp[i].empty())
			tmp[i] = _gc.getGradVector(i, _bfp, pts, getMaxPower());
		for(size_t j = 0; j < pts.numPoints(); j++)
		{
			//if a gradient was not calculated yet, do it
			if(tmp[j].empty())
				tmp[j] = _gc.getGradVector(j, _bfp, pts, getMaxPower());
			//only calculate the dot product, if the points aren't the same
			if(i != j)
				result.push_back(LinAlg::dotProduct(tmp[i], tmp[j]));
		}
	}
	return result;
}

/**
 * This function adds 0's to the list of fit parameters and their errors in order to fill up a certain order. Afterwards, the parameters will be reordered for Professor 2.2.1 usage.
 * That is needed in order to become usable for Professor 2.2.1.
 * @match_index: mapping list that connects the current order to the order needed by Professor 2.2.1
 * @bfp_sorted, @bfperr_sorted: helper to store the new ordering 
 */
void FitHandler::sortFitParams(const vector<int>& match_index){
	
	//push back 0's until there are enough 0's that it matches a certain polynomial order
    for(size_t i = _bfp.size(); i < _qrh.getPower().size(); i++)
		_bfp.push_back(0);
	for(size_t i = _bfperr.size(); i < _qrh.getPower().size(); i++)
		_bfperr.push_back(0);

	//set up the helper
	vector<double> bfp_sorted, bfperr_sorted;
	bfp_sorted.resize(_bfp.size());
	bfperr_sorted.resize(_bfperr.size());
	
	//filling the helper
	for(size_t i = 0; i < match_index.size(); i++)
	{
		bfp_sorted[i] = _bfp[match_index[i]];
		bfperr_sorted[i] = _bfperr[match_index[i]];
	}
	
	//set the ordered list as the entries of the member variables
	_bfp = bfp_sorted;
	_bfperr = bfperr_sorted;
}