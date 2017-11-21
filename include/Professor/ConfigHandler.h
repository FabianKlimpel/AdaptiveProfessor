#ifndef __CONFIGHANDLER__H
#define __CONFIGHANDLER__H

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

/**
* This class handles a config file provided by the user.
* It handles the reading of the file and stores all necessary information in the object.
*/
class ConfigHandler
{
public:
	//Constructor. If no filename is given, default parameters are set
	ConfigHandler(const string configfile);

	const int getSummaryFlag(){return _summaryflag;}
	const int getOutDotFlag(){return _outdotflag;}
	const int getCovmatFlag(){return _covmat;}
	const double getThresholdFit(){return _thresholdfit;}
	const double getThresholdData(){return _thresholddata;}
	const double getThresholdErr(){return _thresholderr;}
	const double getChi2Mean(){return _chi2mean;}
	const double getKappa(){return _kappa;}
	const double getExponent(){return _exponent;}

private:
	//setter for the member variables
	void readThresholdFit(const string line);
	void readThresholdData(const string line);
	void readThresholdErr(const string line);
	void readChi2Mean(const string line);
	void readKappa(const string line);
	void readExponent(const string line);
	void readSummaryFlag(const string line);
	void readOutDotFlag(const string line);
	//~ void readRunCombs(const string line);
	//~ void readLeaveOut(const string line);
	//~ void readRngSeed(const string line);
	void readCovMat(const string line);

	/**
	* @_summaryflag, @_outdotflag: flags for writing additional summaries
	* @_covmat: flag for writing covariance matrices
	* @_thresholdfit: threshold for the RR constrain in the fitting
	* @_thresholddata: threshold for the RR constrain in the hypercube-fitting of the data
	* @_thresholderr: threshold for the RR constrain in getter of the fitting errors
	* @_chi2mean: number of modified chi2's to store in order to state a best value regarding the modified chi2
	* @_kappa: shifting parameter if the RR constrain is applied
	* @_exponent: exponent for the distance weighting of the hypercube
	*/
	int _summaryflag = 1, _outdotflag = 0, _covmat = 1;
	double _thresholdfit = 1e-10, _thresholddata = 1e-10, _thresholderr = 1e-10, _chi2mean = 100, _kappa = 1e-10, _exponent = 1;
};

#endif
