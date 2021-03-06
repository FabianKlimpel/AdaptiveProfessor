#include "Professor/Ipol.h"
#include "Eigen/SVD"
#include <sstream>
#include <cassert>
#include <cmath>
#include <set>

namespace Professor {

  using namespace std;
  using namespace Eigen;


  namespace { //< hide this symbol, since not in API

    // Scaling function to map x from [a,b] into [1,2].
    // NB. Target range does not touch 0, e.g. [0,1] to avoid raising very small numbers to large powers.
    double map_prange(double x, double a, double b) {
      return (x-a)/(b-a);
    }
  }


  // NB. Not a member function
  int numCoeffs(int dim, int order) {
    int ntok = 1;
    int r = min(order, dim);
    for (int i = 0; i < r; ++i) {
      ntok = ntok*(dim+order-i)/(i+1);
    }
    return ntok;
  }


  // NB. Not a member function
  std::vector<double> calcCoeffs(const ParamPoints& pts, const vector<double>& vals, int order,
                                 double threshold, const vector<vector<int> >& structure) {

    // Early exit if this is a trivial 0th order polynomial
    vector<double> rtn;
    if (order == 0) {
      rtn.push_back(vals[0]);
      return rtn;
    }

    // Check the inputs
    if (pts.numPoints() != vals.size())
      throw IpolError("pts.numPoints() != vals.size() in calcCoeffs");
    const int ncoeff = numCoeffs(pts.dim(), order);
    if (ncoeff > pts.numPoints()) {
      stringstream ss;
      ss << "Ipol: not enough (" << ncoeff << " vs. " << pts.numPoints() << ") anchor points "
         << "for interpolating with " << pts.dim() << " params at order " << order;
      for (unsigned int i_order=1;i_order<order;i_order++) {
        if (numCoeffs(pts.dim(), i_order)<=pts.numPoints())
          ss << "\n Order " << i_order  << " requires " << numCoeffs(pts.dim(), i_order) << " anchors";
      }
      throw IpolError(ss.str());
    }

    // Create Eigen objects for the SVD solving
    MatrixXd DP = MatrixXd(pts.numPoints(), ncoeff);
    VectorXd MC = VectorXd(pts.numPoints());

    // The parameter scaling business
    std::vector<std::vector<double> > origpoints = pts.points();
    std::vector<std::vector<double> > scaledpoints;
    std::vector<double> minPV = pts.ptmins();
    std::vector<double> maxPV = pts.ptmaxs();

    for (int p = 0; p < origpoints.size(); ++p) {
      std::vector<double> temp;
      for (int i = 0; i < pts.dim(); ++i) {
        temp.push_back(map_prange(origpoints[p][i], minPV[i], maxPV[i]));
      }
      scaledpoints.push_back(temp);
    }


    // Populate the matrix to be inverted
    vector<double> tempLV;
    for (int a = 0; a < pts.numPoints(); ++a) {
      tempLV = mkLongVector(scaledpoints[a], order, structure);
      for (size_t i = 0; i < tempLV.size(); ++i) {
        DP(a, i) = tempLV[i];
      }
      // The vector of values (corresponding to anchors)
      MC[a] = vals[a];
    }
    JacobiSVD<MatrixXd> svd = DP.jacobiSvd(ComputeThinU|ComputeThinV);
    #if EIGEN_WORLD_VERSION >= 3 && EIGEN_MAJOR_VERSION >= 2 && EIGEN_MINOR_VERSION >= 1
    svd.setThreshold(threshold); // Needed TODO find transform for dependence on stuff
    #endif

    // Check for singular values, i.e. fully correlated parameters
    /// @todo Maybe figure out how to use Eigen's setThreshold better?
    VectorXd svals = svd.singularValues();
    for (unsigned int i = 0; i < svd.nonzeroSingularValues();++i) {
      if (fabs(svals[i]) < threshold) {
        std::cout << "Singular value encountered, aborting" << std::endl;
        abort();
      }
    }

    // Solve for coefficients
    VectorXd co = svd.solve(MC);

    // Populate the coefficient std::vector and return
    for (size_t i = 0; i < ncoeff; ++i) rtn.push_back(co[i]);
    return rtn;
  }


  double calcValue(const vector<double>& params,
                   const vector<double>& coeffs, int order,
                   const vector< vector<int> >& structure) {
    const vector<double> lv = mkLongVector(params, order, structure);
    return calcValue(lv, coeffs);
  }


  double calcValue(const vector<double>& paramslongvector,
                   const vector<double>& coeffs) {
    // Dot product of params long-vector with coeffs -> value
    assert(paramslongvector.size() == coeffs.size());
    double v = 0.0;
    for (size_t i = 0; i < paramslongvector.size(); ++i) {
      //cout << i << ": " << coeffs[i] << " * " << paramslongvector[i] << " = " << coeffs[i]*paramslongvector[i] << endl;
      v += coeffs[i] * paramslongvector[i];
    }
    return v;
  }


  vector<vector<int> > mkStructure(int dim, int order) {
    if (order < 0)
      throw IpolError("Polynomial order " + to_string(order) + " not implemented");

    const vector<int> zero(dim, 0);
    set< vector<int> > rtn;
    rtn.insert(zero);



    if (order>0) {
    // The set of parameter base vectors
      vector<vector<int> > BS;
      for (unsigned int d = 0; d < dim; ++d) {
        vector<int> p(dim,0);
        p[d]=1;
        rtn.insert(p);
        BS.push_back(p);
      }


      auto temp = BS;
      vector<vector<int> > temp2;
      vector<int> e(dim,0);

      // Recursively add base vectors
      for (unsigned int o = 1; o < order; ++o) {
        temp2.clear();
       
        for ( auto const & t : temp) {
          for (auto const & bs : BS) {
            // Create a new element
            for (unsigned int d = 0; d < dim; d++) {
              e[d] = t[d] + bs[d];
            }
            temp2.push_back(e);
          }
        }
        temp=temp2; // For the next order, we want to add base vectors 
                    // to each element of the current order
        for (auto const &v : temp2) {
          rtn.insert(v); // The set takes care of not having duplicates.
        }
      }
    }

    // Convert back to vector<vector<int> > TODO rewrite everything else such that it takes set<vector<int> > ?
    vector<vector<int> > rtn2;
    for (auto const &v : rtn) rtn2.push_back(v);
    return rtn2;
  }


  vector<double> mkLongVector(const vector<double>& params, int order, const vector< vector<int> >& structure) {
    if (order < 0)
      throw IpolError("Polynomial order " + to_string(order) + " not implemented");

    vector<double> rtn;
    for (const vector<int>& v : structure) {
      double prod = 1.0;
      for (size_t i = 0; i < v.size(); ++i) {
        if (v[i] == 0) continue;
        /// @todo Can be speeded with (precomputable?) integer powers / exp-by-doubling?
        prod *= std::pow(params[i], v[i]);
      }
      rtn.push_back(prod);
    }
    return rtn;
  }


  /// @todo Why the min/maxPV args?
  /// @todo Expose to API
  vector<double> mkLongVectorDerivative(const vector<double>& params, int order,
                                        const vector<double>& minPV, const vector<double>& maxPV,
                                        const vector<vector<int> >& structure) {
    if (order < 0)
      throw IpolError("Polynomial order " + to_string(order) + " not implemented");

    vector<double> rtn;
    bool firstItem = true;
    for (const vector<int>& s : structure) {

      if (firstItem) {
        rtn.push_back(0.0); // Derivative of constant term
        firstItem = false;
        continue;
      }
      double part = 0.0;
      // Differentiate x^a*y^b*z^c*...
      for (unsigned int c = 0; c < s.size(); c++) { // d/dx, d/dy, d/dz, ...

        double temp2 = 1.0;
        for (unsigned int i = 0; i <s.size(); i++) { // x, y, z
          if (c==i) {  // d/dx x*y*z
            temp2 *= s[i];
            if (s[c] == 0) continue;
            else temp2 *= std::pow(params[i], s[i]-1)/(maxPV[i]- minPV[i]); // Jacobian factor: 'd map_prange / dx' = 1./(b-a)
          } else {
            temp2 *= std::pow(params[i], s[i] );
          }
        }
        part += temp2;
      }
      rtn.push_back(part);
    }

    return rtn;
  }


  /// @todo Why the min/maxPV args?
  /// @todo Expose to API
  vector<double> mkLongVectorGradient(const vector<double>& params, int coord, int order,
                                      const vector<double>& minPV, const vector<double>& maxPV,
                                      const vector<vector<int> >& structure) {
    if (order < 0)
      throw IpolError("Polynomial order " + to_string(order) + " not implemented");

    vector<double> rtn;
    bool firstItem = true;
    for (const vector<int>& s : structure) {
      if (firstItem) {
        rtn.push_back(0.0); // Derivative of constant term
        firstItem = false;
        continue;
      }

      if (s[coord] == 0) {
        rtn.push_back(0);
        continue;
      }
      double temp = 1.0;
      for (unsigned int i = 0; i <s.size(); i++) { // x, y, z
        if (i == coord) {  // d/dx x*y*z
          temp *= s[i];  // d/dx  x^a = a*x^(a-1)
          temp *= std::pow(params[i], s[i]-1)/(maxPV[i]- minPV[i]); // Jacobian factor: 'd map_prange / dx' = 1./(b-a)
        } else {
          temp *= std::pow(params[i], s[i] );
        }
      }
      rtn.push_back(temp);
    }

    return rtn;
  }


  ///////////////////////////////////////////////////////

/**
 * Constructor for adaptive interpolation
 * @pts: Anchor points
 * @ptvlas: Response of each point in @pts
 * @pterrs: Uncertainty of each response in @ptvals
 * @num_ipol: Index of the bin, used for output
 * @doipol: Boolean flag. If true then an interpolation will be calculated, otherwise the uncertainty of the fit parameters is performed
 * @order: If a fit function is already calculated, the order of the polynomial is passed to improve the uncertainty calculation speed
 * @fitparams: If a fit function is already calculated, this parameter represents the fit parameters for an uncertainty calculation
 * @configfile: Path to a file that contains additional configuration information
 */
Ipol::Ipol(ParamPoints& pts, const std::vector<double>& ptvals, const vector<double>& pterrs, const int num_ipol, const string histname, const bool doipol, const int order, const std::vector<double>& fitparams, const string configfile){
	
	//rescale parameter ranges
	pts.rescale();
	
	//either perform an interpolation or an uncertainty calculation
	if(doipol)
		calcipol(pts, ptvals, pterrs, num_ipol, configfile, histname);
	else
		calcerr(pts, pterrs, order, fitparams, num_ipol, histname);
	
	//clean memory
	pts.clearAll();
}
    
/**
 * This function calculates the interpolation.
 * @pts: Anchor points
 * @ptvlas: Response of each point in @pts
 * @pterrs: Uncertainty of each response in @ptvals
 * @num_ipol: Index of the bin, used for output
 * @configfile: Path to a file that contains additional configuration information
 * @ch: Object that loads and stores the infromation in the configuration file
 * @oh: Object that handles the output
 * @fh: Object that handles the interpolation
 * @qp: Structure that stores the last couple quality parameters and tracks the currently best iteration
 */ 
void Ipol::calcipol(ParamPoints& pts, const std::vector<double>& ptvals, const std::vector<double>& pterrs, const int num_ipol, const std::string configfile, const std::string histname) {

	cout << "start with bin " << num_ipol << " in " << histname << endl;
		
	//Crate and initialize objects
	ConfigHandler* ch = new ConfigHandler(configfile);
	OutputHandler oh(pts, ch->getOutDotFlag(), ch->getSummaryFlag(), num_ipol);
	FitHandler fh;	
	fh.startBin(pts, ptvals, pterrs, ch->getThresholdFit(), ch->getKappa());
	qualityparams qp;
	qp.bestiteration = fh.getIterationCounter();
	qp.bestF = fh.getChi2();
	qp.Fvalues.push_back(fh.getChi2());
	
	//calculate the new normal vectors of the data
	pts.setGradients(ptvals);
	cout << "start iterating for bin " << num_ipol << endl;
	//continuously iterate while a fit becomes better with more iterations
	while(iterate(fh, pts, ch, qp)){}

	cout << "fit complete for bin " << num_ipol << " in " << histname << endl;
	
	//recalculate the best iteration
	fh.setToIteration(pts, ptvals, pterrs, num_ipol, ch->getThresholdFit(), ch->getKappa(), qp.bestiteration);			
	
	//write out the result
	oh.writeBinResult(num_ipol, pts, fh);
	if(ch->getOutDotFlag())
		oh.writeDotProduct(num_ipol, pts, fh);
	if(ch->getSummaryFlag())
		oh.writeSummary(fh, pts);
		
	delete(ch);
	qp.Fvalues.clear();
	//set parameters used by Professor
	setparams(pts, fh);
}

/**
 * This function performs an iteration of the adaptive interpolation
 * @fh: Object that handles the interpolation
 * @pts: Anchor points
 * @ch: Object that loads and stores the infromation in the configuration file
 * @qp: Structure that stores the last couple quality parameters and tracks the currently best iteration
 * @mean: Storage of the mean of the F-values stored in @qp.
 */
const bool Ipol::iterate(FitHandler& fh, ParamPoints& pts, ConfigHandler* const ch, qualityparams& qp) const{
	//iterate
	fh.nextStep(pts, ch->getThresholdFit(), ch->getKappa());
	
	//store quality parameter
	qp.Fvalues.push_back(fh.getChi2() * (1 - fh.getDsmooth(pts)) / (1 + fh.getDsmooth(pts)));
	
	//if enough iterations were made, the convergence check will be performed
	if(qp.Fvalues.size() >= ch->getChi2Mean())
	{
		double mean = 0.;
		//calculate the mean of @qp.Fvalues
		for(size_t i = 0; i < qp.Fvalues.size(); i++)
			mean += qp.Fvalues[i];
		mean /= qp.Fvalues.size();

		//if the current value is equal or worse the sliding mean, the loop ends
		if(qp.Fvalues[qp.Fvalues.size() - 1] >= mean)
			return false;

		//delete the first entry in @qp.Fvalues -> list length remains constant
		qp.Fvalues.erase(qp.Fvalues.begin());
	}

	//store 
	if(qp.Fvalues.back() < qp.bestF)
	{
		qp.bestF = qp.Fvalues.back();
		qp.bestiteration = fh.getIterationCounter();
	}

	return true;
}

/**
 * This function sets the member parameters such that Professor can work without additional interfaces
 * @pts: Anchor points
 * @fh: Object that handles the interpolation
 * @match_index: List of indices; Used for re-ordering the fit parameters to fit Professor 2.2.1's order
 */
void Ipol::setparams(ParamPoints& pts, FitHandler& fh){
	
	_dim = pts.dim();	
	_minPV = pts.ptmins();
	_maxPV = pts.ptmaxs();
	_name = "";
	_order = fh.getMaxPower();
	
	//sort fit parameters
	std::vector<int> match_index = sort_strucs(pts);
	fh.sortFitParams(match_index);	
	match_index.clear();

	_coeffs = fh.getFitParams();
	_structure = mkStructure(_dim, _order);
}

/**
 * This function performs a re-ordering of the order of the fit parameters in order to fit Professor 2.2.1's order.
 * @pts: Anchor points
 * @profstruc: The structure of the fit parameters of Professor 2.2.1
 * @result: List of mapping numbers for later usage
 */
const std::vector<int> Ipol::sort_strucs(ParamPoints& pts) const{

	//calculate the Professor 2.2.1 structure
	std::vector<std::vector<int>> profstruc = mkStructure(_dim, _order);
	std::vector<int> result;

	//walk over both list, if there is a match, write down the number in @result
	for(size_t pi = 0; pi < profstruc.size(); pi++)
		for(size_t i = 0; i < pts.getPower(_order).size(); i++)
			if(pts.getPower(_order)[i] == profstruc[pi])
				result.push_back(i);
	return result;	
}

void Ipol::calcerr(ParamPoints& pts, const std::vector<double>& pterrs, const int order, const std::vector<double>& fitparams, const int num_ipol, const string histname){

	std::cout << "calculating error for bin " << num_ipol << " in " << histname << endl;
		
	FitHandler fh(fitparams, pterrs, pts, order, num_ipol, histname);
	setparams(pts, fh);
	
	cout << "error calculation complete for bin " << num_ipol << " in " << histname << endl;
}

  ///////////////////////////////////////////////////////

  string Ipol::toString(const string& name) const {
    stringstream ss;
    if (!name.empty()) ss << name << ": ";
    else if (!_name.empty()) ss << _name << ": ";
    ss << this->dim() << " ";
    ss << this->order() << " ";
    for (const double& a : coeffs())
      ss << a << " ";
    return ss.str();
  }
  
  /// TODO: How do we want to read in the MinMaxValues here?
  void Ipol::fromString(const string& s) {
    // Extract a name if given at the start of the string
    _name = (s.find(":") != std::string::npos) ? s.substr(0, s.find(":")) : "";
    // Load the rest of the string into a stringstream and load into numerical variables
    istringstream numss( (s.find(":") != std::string::npos) ? s.substr(s.find(":")+1) : s );
    numss >> _dim;
    numss >> _order;
    double tmp; int ncoeffs = 0;
    while (numss >> tmp) {
      _coeffs.push_back(tmp);
      ncoeffs += 1;
    }
    assert(ncoeffs == numCoeffs(dim(),order()));
    _structure = mkStructure(dim(), order());
  }


  vector<double> Ipol::sparams(const vector<double>& params) const {
    if (params.size() != dim()) {
      stringstream ss;
      ss << "Incorrect number of parameters given ("
         << dim() << " params required, " << params.size() << " supplied)";
      throw IpolError(ss.str());
    }

    // Param scaling into [0,1] ranges defined by sampling limits (if set)
    vector<double> sparams = params;
    if (!_minPV.empty() && !_maxPV.empty()) {
      for (size_t i = 0; i < dim(); ++i) {
        sparams[i] = map_prange(params[i], _minPV[i], _maxPV[i]);
      }
    }
    return sparams;
  }


  double Ipol::value(const vector<double>& params) const {
    return calcValue(sparams(params), coeffs(), order(), _structure);
  }


  double Ipol::derivative(const vector<double>& params) const {
    /// @todo Extract into a standalone calc function
    if (params.size() != dim()) {
      stringstream ss;
      ss << "Incorrect number of parameters passed to Ipol::derivative ("
         << dim() << " params required, " << params.size() << " supplied)";
      throw IpolError(ss.str());
    }

    // Param scaling into [0,1] ranges defined by sampling limits (if set)
    vector<double> sparams = params;
    if (!_minPV.empty() && !_maxPV.empty()) {
      for (size_t i = 0; i < dim(); ++i) {
        sparams[i] = map_prange(params[i], _minPV[i], _maxPV[i]);
      }
    }

    // Dot product for value
    const vector<double> lv = mkLongVectorDerivative(sparams, order(), _minPV, _maxPV, _structure);
    assert(lv.size() == coeffs().size());
    double v = 0.0;
    for (size_t i = 1; i < lv.size(); ++i) {
      v += lv[i] * coeff(i);
    }
    return v;
  }


  vector<double> Ipol::gradient(const vector<double>& params) const {
    /// @todo Extract into a standalone calc function
    if (params.size() != dim()) {
      stringstream ss;
      ss << "Incorrect number of parameters passed to Ipol::gradient ("
         << dim() << " params required, " << params.size() << " supplied)";
      throw IpolError(ss.str());
    }

    vector<double> grad;

    // Param scaling into [0,1] ranges defined by sampling limits (if set)
    vector<double> sparams = params;
    if (!_minPV.empty() && !_maxPV.empty()) {
      for (size_t i = 0; i < dim(); ++i) {
        sparams[i] = map_prange(params[i], _minPV[i], _maxPV[i]);
      }
    }

    for (int c=0; c< params.size(); c++) {
      // Dot product for value
      const vector<double> lv = mkLongVectorGradient(sparams, c, order(), _minPV, _maxPV, _structure);
      assert(lv.size() == coeffs().size());
      double v = 0.0;
      for (size_t i = 1; i < lv.size(); ++i) {
        v += lv[i] * coeff(i);
      }
      grad.push_back(v);

    }

    return grad;
  }
}
