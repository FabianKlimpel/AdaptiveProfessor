#include "Professor/OutputHandler.h"

/**
 * Constructor
 * Here, the sizes of the vectors that will later contain fitparameters etc. will be set. Additionally, the distances and centers will be calculated
 * @ch: carries the flag for the dotproduct-summary
 * @rh: handler of the reference data
 * @minparval: minimum of the parametervalues
 * @maxparval: maximum of the parametervalues
 */
OutputHandler::OutputHandler(const Professor::ParamPoints& pts, const int& outdotflag, const int& summaryflag, const size_t& num_ipol){	
	//if the dotproduct-summary should be written, the distances will be calculated
	if(outdotflag)
		//Calculate the distances between the anchor points
		setDistances(pts);
	if(summaryflag && (num_ipol == 0))
		setupSummary();
}	

OutputHandler::OutputHandler(){}

/**
 * This function calculates the distances between the anchors points
 * @ch: carries the number of MC runs
 * @rh: handles the reference data
 * @skips: no distances between an anchor point and itself is calculated, so the index shift needs to be registered by this variable
 */
void OutputHandler::setDistances(const Professor::ParamPoints& pts){

	//Resizing the @distances vectors. No distances between an anchor point and itself will be calculated, therefore the size is modified accordingly.
	distances.resize(pts.numPoints() * pts.numPoints() - pts.numPoints());
	
	for(size_t i = 0; i < pts.numPoints(); i++)
	{		
		size_t skips = 0;
		//compare every point with every other
		for(size_t j = 0; j < pts.numPoints(); j++)
			//if both points are the same, the iteration is skipped
			if(i == j)
				skips--;				
			else
				//calculate the distance between the vectors and store them
				distances[i * pts.numPoints() - i - 1 + j + skips] = LinAlg::getDistanceOfVectors(pts.pointScaled(i), pts.pointScaled(j));
	}
}

/**
 * This function creates the summary file and writes its header
 * @outsummary: delivers the output to file functionality
 */
void OutputHandler::setupSummary() const{	
	//set up the file and write the header
	ofstream outsummary;
	outsummary.open("summary");
	outsummary << "Chi^2" << "\t" << "Chi2^2,red" << "\t" << "Iterations" << "\t" << "Dsmooth" << endl;
	outsummary.close();
}

/**
 * This function write the result of the fit for a bin to the console
 * @num_ipol: # of the bin
 * @ch: container for the observable name and the analysis name
 * @rh: handler of the reference data
 * @fh: handler of the fit
 * @store: counter to identify the right analysis/observable
 */
void OutputHandler::writeBinResult(const size_t num_ipol, Professor::ParamPoints& pts, FitHandler& fh) const{
	cout << endl << "Result for bin " << num_ipol << ":" << endl;

	//the constrained monomial numbers will be printed to the terminal
	cout << "RR constraint:\t";
	for(size_t i = 0; i < fh.getNumFitParams(); i++)
		if(!std::isnan(fh.get_a()[i]))
			cout << i << "\t";
	cout << endl;
	
	//write further summary variables will be printed to the terminal
	cout << "Dsmooth:\t\t" << fh.getDsmooth(pts) << endl;
	cout << "chi2:\t\t\t" << fh.getChi2() << endl;
	cout << "iterationcounter:\t" << fh.getIterationCounter() << endl;
	cout << "max. power:\t\t" << fh.getMaxPower() << endl;
	cout << "-------------------------" << endl;	
}

/**
 * This function writes the dotproduct-summary
 * @num_ipol: # of the bin
 * @ch: contains the number of runs
 * @rh: handles the reference data
 * @fh: handles the fit
 * @outdot: delivers the output to file functionality
 * @rhdot, @fhdot: stores the dot products of the reference data and the fit at every anchor point
 */
void OutputHandler::writeDotProduct(const size_t num_ipol, Professor::ParamPoints& pts, FitHandler& fh) const{
	//set up the output
	ofstream outdot;
	outdot.open(("dotproduct" + to_string(num_ipol)).c_str());
	
	//calculate all dot products
	vector<double> rhdot = pts.getAllGradDotProducts(), fhdot = fh.getAllGradDotProducts(pts);
	
	//write some summary parameters
	outdot << fh.getDsmooth(pts) << "\t" << fh.getChi2() << "\t" << fh.getIterationCounter() << endl;
		
	//write the dot products of the reference data, the fit and the distances between the anchor points used for the dot product
	for(size_t k = 0; k < rhdot.size(); k++)
		outdot << rhdot[k] << "\t" << fhdot[k] << "\t" << distances[k] << "\t";
	outdot << endl;
						
	outdot.close();
}
	
/**
 * This function writes a summary of a fit to the summary file
 * @fh: handles the fit
 * @rh: handles the reference data
 * @outsummary: delivers the output to file functionality
 */
void OutputHandler::writeSummary(FitHandler& fh, Professor::ParamPoints& pts) const{
	//open the file and continue writing at its end
	ofstream outsummary;
	outsummary.open("summary", ofstream::out | ofstream::app);
	//write out the resulting Chi^2, the reduced Chi^2, the number of iterations needed and the smoothness
	outsummary << fh.getChi2() << "\t" << fh.getIterationCounter() << "\t" << fh.getDsmooth(pts) << endl;
	outsummary.close();
}

void OutputHandler::writeCovMat(const MatrixXd& mat, const size_t num_ipol) const{
	ofstream outcovmat;
	outcovmat.open(("covmat_" + to_string(num_ipol)).c_str());
	
	for(size_t row = 0; row < (size_t) mat.rows(); row++)
	{
		for(size_t col = 0; col < (size_t) mat.cols(); col++)
			if(mat(row, col) == std::numeric_limits<double>::infinity())
				outcovmat << std::numeric_limits<double>::max() << " ";
			else
				if(-mat(row, col) == std::numeric_limits<double>::infinity())
					outcovmat << -std::numeric_limits<double>::max() << " ";
				else
					outcovmat << mat(row, col) << " ";
		outcovmat << "\n";
	}
	outcovmat.close();
}
