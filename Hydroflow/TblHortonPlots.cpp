/*IlwisObjects is a framework for analysis, processing and visualization of remote sensing and gis data
Copyright (C) 2018  52n North

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

#include <functional>
#include <future>
#include "kernel.h"
#include "raster.h"
#include "symboltable.h"
#include "ilwisoperation.h"
#include "blockiterator.h"
#include "pixeliterator.h"

#include "featurecoverage.h"
#include "basetable.h"
#include "table.h"
#include "flattable.h"
#include "feature.h"
#include "featureiterator.h"
#include "conventionalcoordinatesystem.h"
#include "geos/geom/CoordinateArraySequence.h"
#include "geos/geom/CoordinateSequence.h"
#include "geos/geom/CoordinateSequenceFactory.h"
#include "geos/geom/GeometryFactory.h"
#include "geos/geom/LineString.h"
#include "geos/geom/Point.h"
#include "geos/geom/Polygon.h"
#include "coordinatedomain.h"
#include "ellipsoid.h"


#include "TblHortonPlots.h"


using namespace Ilwis;
using namespace Hydroflow;
using namespace std;

REGISTER_OPERATION(TableHortonPlots)


////////////////////////////////////////////////////////////

const double rDefaultEarthRadius = 6371007.0;

static double rSphericalDistance(double rRadius, const LatLon& ll_1, const LatLon& ll_2)
{
	if (ll_1.Lat() == rUNDEF || ll_1.Lon() == rUNDEF || ll_2.Lat() == rUNDEF || ll_2.Lon() == rUNDEF)
		return rUNDEF;
	double phi1 = ll_1.Lat() * M_PI / 180.0; //conversion to radians
	double lam1 = ll_1.Lon() * M_PI / 180.0;
	double phi2 = ll_2.Lat() * M_PI / 180.0; ;
	double lam2 = ll_2.Lon() * M_PI / 180.0; ;
	double sinhalfc = fabs(sin((phi2 - phi1) / 2) * sin((phi2 - phi1) / 2) +
		cos(phi1) * cos(phi2) * sin((lam2 - lam1) / 2) * sin((lam2 - lam1) / 2));
	sinhalfc = sqrt(sinhalfc);
	double c; // the shortest spherical arc
	if (sinhalfc < sqrt(2.0) / 2)
		c = 2.0 * asin(sinhalfc);
	else
	{
		phi2 = -phi2;
		lam2 = M_PI + lam2;
		sinhalfc = fabs(sin((phi2 - phi1) / 2) * sin((phi2 - phi1) / 2) +
			cos(phi1) * cos(phi2) * sin((lam2 - lam1) / 2) * sin((lam2 - lam1) / 2));
		sinhalfc = sqrt(sinhalfc);
		c = M_PI - 2.0 * asin(sinhalfc);
	}
	return c * rRadius;
}


// ####################

CVector LeastSquaresFit::LeastSquares(const RealMatrix& A, const CVector& b)
{
	assert(A.iRows() >= A.iCols());
	assert(A.iRows() == b.iRows());
	int i, j, k;
	int M, N;
	M = A.iCols();
	N = A.iRows();

	RealMatrix AtA(M);
	for (i = 0; i < M; ++i)
		for (j = 0; j < M; ++j) {
			double rSum = 0;
			for (k = 0; k < N; ++k)
				rSum += A(k, i) * A(k, j);
			AtA(i, j) = rSum;
		}

	CVector Atb(M);
	for (i = 0; i < M; ++i) {
		double rSum = 0;
		for (j = 0; j < N; ++j)
			rSum += A(j, i) * b(j);
		Atb(i) = rSum;
	}

	CVector vRes(M);
	AtA.Solve(Atb, vRes);
	//  vRes.SetError(AtA.err());
	return vRes;
}

LeastSquaresFit::LeastSquaresFit(const CVector& cvXData, const CVector& cvYData, int iTerms)
: cvX(cvXData), cvY(cvYData), 
  cvSolution(iTerms), mtBasis(cvXData.iRows(), iTerms),
  cvW(cvXData), cvZ(cvYData),
  iNrTerms(iTerms), 
  iNrPoints(cvXData.iRows())
{
}

LeastSquaresFit::~LeastSquaresFit()
{
}

static void ErrorTooLessPoints()
{
  throw ErrorObject(TR("Not enough points"));
}

static void ErrorTooLessTerms()
{
  throw ErrorObject(TR("Not enough terms"));
}

static void ErrorMoreTermsThanPoints()
{
  throw ErrorObject(TR("Points terms mismatch"));
}

static void ErrorNoPossibleSolution()
{
  throw ErrorObject(TR("No possible solution"));
}

void LeastSquaresFit::Compute()
{
  InitializeAndFormBasisVectors();
  CreateAndSolveEquations();
  
  vector<vector<double>> mb(iNrPoints, vector<double>(iNrTerms, 0.0));

  for (int index = 0; index < iNrPoints; ++index)
  {
    for (int term = 0; term < iNrTerms; ++term) 
	{
		mb[index][term] = mtBasis(index,term);
	}
  }

  FinalizeSolution();
   double a1 = cvSolution(0);
  double a2 = cvSolution(1);
  int a =0;
}

void LeastSquaresFit::InitializeAndFormBasisVectors()
{
  if (iNrPoints < 2)
    ErrorTooLessPoints();
  if (iNrTerms <= 1)
    ErrorTooLessTerms();
  if (iNrTerms >= iNrPoints)
    ErrorMoreTermsThanPoints();
  for (int i = 0; i < cvSolution.iRows(); ++i)
    cvSolution(i) = 0;
  Transform();
  CreateBasisFunctions();      
}

void LeastSquaresFit::CreateAndSolveEquations()
{
  cvSolution = LeastSquares(mtBasis, cvZ);
  double a1 = cvSolution(0);
  double a2 = cvSolution(1);
  int a =0;

}

void LeastSquaresFit::FinalizeSolution()
{
  double rSum = 0;
  for (int index = 0; index < iNrPoints; ++index) {
    double rYFit = 0;
    for (int term = 0; term < iNrTerms; ++term) 
      rYFit += cvSolution(term) * mtBasis(index, term);
    InverseTransform(rYFit);
    rSum += (rYFit - cvY(index)) * (rYFit - cvY(index));
  }
  if (iNrPoints == iNrTerms)
    _rStdDev = 0;
  else 
    _rStdDev = sqrt(rSum / (iNrPoints - iNrTerms));  
    
  TransformSolution();
}

void LeastSquaresFit::Transform()
{
  // default implementation is empty
  // because constructor already copies cvX to cvW and cvY to cvZ
}

void LeastSquaresFit::TransformSolution()
{
  // default implementation is empty
}

void LeastSquaresFit::InverseTransform(double rYFit)
{
  // default implementation is empty
}


LeastSquaresFitPolynomial::LeastSquaresFitPolynomial(const CVector& cvXData, const CVector& cvYData, int order)
: LeastSquaresFit(cvXData, cvYData, order)
{
  Compute();
}

void LeastSquaresFitPolynomial::WritePoly(int iTerm, ostream& os) const
{
  double r = cvSolution(iTerm);
  if (abs(r) < EPS10)  r = 0.0;
  os << r;
  iTerm += 1;
  if (iTerm == iNrTerms)
    return;
  os << " + X * (";
  WritePoly(iTerm, os);
  os << ")";  
}

QString LeastSquaresFitPolynomial::sFormula() const
{
//  String sText = String("Equation of Fitting Polynomial Function:\r\n");
	QString sText = QString("Y =");
  int iTerm = 0;
  double r = cvSolution(0);
  if (abs(r) < EPS10)  r = 0.0;
  QString sPolyn;
  sPolyn.sprintf("  %.6f", r);
  for ( iTerm = 1; iTerm < iNrTerms; iTerm++) {
    r = cvSolution(iTerm);
    if (abs(r) < EPS10)  r = 0.0;
	QString temp;
    if (iTerm == 1)
	  temp.sprintf("  + %.6f * X", r);
	else
	  temp.sprintf("  + %.6f * X^%d", r, iTerm);
	sPolyn += temp;
  }
  /*char sBuf[400];
  ostrstream str(sBuf, 200);
  WritePoly(0,str);
  String s(str.str());
  */
  sText += sPolyn;
  return sText;
}

double LeastSquaresFitPolynomial::rFuncPoly(int iTerm, double rX) const
{
  double a = cvSolution(iTerm);
  if (iTerm + 1 == iNrTerms)
    return cvSolution(iTerm);
  else
    return cvSolution(iTerm) + rX * rFuncPoly(iTerm+1, rX);  
}

double LeastSquaresFitPolynomial::rCalc(double rX) const
{
  return rFuncPoly(0, rX);
}

void LeastSquaresFitPolynomial::Transform()
{
  int i;
  double rXMin = cvX(0);
  double rXMax = cvX(0);
  for (i = 1; i < iNrPoints; ++i) {
    if (rXMin > cvX(i))
      rXMin = cvX(i);
    else if (rXMax < cvX(i))
      rXMax = cvX(i);
  }    
  rMultiplier = 2.0 / (rXMax - rXMin);
  rConstant = - rMultiplier * (rXMax + rXMin) / 2.0;
  std::vector<double> cw;
  for (i = 0; i < iNrPoints; ++i)
  {
    cvW(i) = rMultiplier * cvX(i) + rConstant;
	cw.push_back(rMultiplier * cvX(i) + rConstant);
  }
  int a=0;
//  ZDAta := YData;
}

void LeastSquaresFitPolynomial::CreateBasisFunctions()
{
  int i, j;
  for (i = 0; i < iNrPoints; ++i) {
    mtBasis(i, 0) = 1;
    mtBasis(i, 1) = cvW(i);
    for (j = 2; j < iNrTerms; ++j)
      mtBasis(i, j) = 2 * cvW(i) * mtBasis(i, j - 1) - mtBasis(i, j - 2);
  }    
}

void LeastSquaresFitPolynomial::TransformSolution()
{
  int index, term;
  double rSum;
  CVector cvOldConversion(iNrTerms), cvConversion(iNrTerms);
  CVector cvOldSolution(cvSolution);

  double a1  = cvSolution(0);
  double a2  = cvSolution(1);

  double a3  = cvOldSolution(0);
  double a4  = cvOldSolution(1);

  for (int i = 0; i < iNrTerms; ++i)
    cvOldConversion(i) = 0;
  for (index = 0; index < iNrTerms; ++index) {
    rSum = 0;
    if (index > 0)
      cvConversion(index - 1) = 0;
    for (term = index; term < iNrTerms; ++term) {
      if (term == 0)
        cvConversion(term) = 1.0;
      else if (term == 1)
        if (index == 0)
          cvConversion(term) = rConstant;
        else
          cvConversion(term) = rMultiplier;
      else
        cvConversion(term) = 2 * rMultiplier * cvOldConversion(term - 1)
                           + 2 * rConstant * cvConversion(term - 1)
                           - cvConversion(term - 2);
      rSum += cvConversion(term) * cvOldSolution(term);
    }  
    cvSolution(index) = rSum;
    cvOldConversion = cvConversion;
  }
}


// ####################




/*---------------------------------------------------------
   Chebyshev polynomials  T_k(t)
   t ∈ [-1, 1]，not original x
----------------------------------------------------------*/
double ChebyshevFitter::Chebyshev(int k, double t) {
	if (k == 0) return 1.0;
	if (k == 1) return t;

	double Tkm2 = 1.0;
	double Tkm1 = t;
	double Tk;

	for (int i = 2; i <= k; i++) {
		Tk = 2.0 * t * Tkm1 - Tkm2;
		Tkm2 = Tkm1;
		Tkm1 = Tk;
	}
	return Tk;
}

/*---------------------------------------------------------
   mapping [xmin, xmax] to [-1,1]
   t = (2*x - (xmin + xmax)) / (xmax - xmin)
----------------------------------------------------------*/
double ChebyshevFitter::Map_to_chebyshev_domain(double x, double xmin, double xmax) {
	return (2.0 * x - (xmin + xmax)) / (xmax - xmin);
}

/*---------------------------------------------------------
	Chebyshev fit：LeastSquaresFit
----------------------------------------------------------*/
vector<double> ChebyshevFitter::Chebyshev_fit( const vector<double>& xx,const vector<double>& yy, int n_terms)
{
	vector<double> x, y;

	long iNrValidPnts = 0;
	for (int i = 0; i < xx.size(); i++) 
	{
		if (xx[i] == rUNDEF || yy[i] == rUNDEF)
			continue;
		x.push_back(xx[i]);
		y.push_back(yy[i]);
		iNrValidPnts++;     // count all valid X,Y points
	}
	
	m_x = x;
	m_y = y;
	int n = x.size();
	m_iNrPoints = n;
	m_iNrTerms = n_terms;

	double xmin = *min_element(x.begin(), x.end());
	double xmax = *max_element(x.begin(), x.end());

	m_rMultiplier = 2.0 / (xmax - xmin);
	m_rConstant = -m_rMultiplier * (xmax + xmin) / 2.0;

	// Construct A
	vector<vector<double>> A(n, vector<double>(n_terms, 0.0));
	for (int i = 0; i < n; i++) {
		double t = Map_to_chebyshev_domain(x[i], xmin, xmax);
		for (int j = 0; j < n_terms; j++) {
			A[i][j] = Chebyshev(j, t);
		}
	}

	mtBasis = A;
	// ATA *a = ATy
	vector<vector<double>> ATA(n_terms, vector<double>(n_terms, 0.0));
	vector<double> ATy(n_terms, 0.0);

	for (int i = 0; i < n_terms; i++) {
		for (int j = 0; j < n_terms; j++) {
			for (int k = 0; k < n; k++) {
				ATA[i][j] += A[k][i] * A[k][j];
			}
		}
		for (int k = 0; k < n; k++) {
			ATy[i] += A[k][i] * y[k];
		}
	}

	// Gaussian elimination solution
	vector<double> a(n_terms, 0.0);

	for (int i = 0; i < n_terms; i++) {
		double pivot = ATA[i][i];
		for (int j = i; j < n_terms; j++)
			ATA[i][j] /= pivot;
		ATy[i] /= pivot;

		for (int k = i + 1; k < n_terms; k++) {
			double factor = ATA[k][i];
			for (int j = i; j < n_terms; j++)
				ATA[k][j] -= factor * ATA[i][j];
			ATy[k] -= factor * ATy[i];
		}
	}

	// coefficients
	for (int i = n_terms - 1; i >= 0; i--) {
		a[i] = ATy[i];
		for (int j = i + 1; j < n_terms; j++) {
			a[i] -= ATA[i][j] * a[j];
		}
	}

	return a;
}

void ChebyshevFitter::FinalizeSolution(vector<double>& cvSolution)
{
	m_cvSolution = cvSolution;
	double rSum = 0;
	for (int index = 0; index < m_iNrPoints; ++index) {
		double rYFit = 0;
		for (int term = 0; term < m_iNrTerms; ++term)
			rYFit += cvSolution[term] * mtBasis[index][term];

		rSum += (rYFit - m_y[index]) * (rYFit - m_y[index]);
	}
	if (m_iNrPoints == m_iNrTerms)
		m_rStdDev = 0;
	else
		m_rStdDev = sqrt(rSum / (m_iNrPoints - m_iNrTerms));

	TransformSolution();
}

void ChebyshevFitter::TransformSolution()
{
	int index, term;
	double rSum;
	vector<double>cvOldConversion(m_iNrTerms);
	vector<double>cvConversion(m_iNrTerms);
	vector<double>cvOldSolution(m_cvSolution);

	/*CVector cvOldConversion(iNrTerms), cvConversion(iNrTerms);
	CVector cvOldSolution(cvSolution);*/
	for (int i = 0; i < m_iNrTerms; ++i)
		cvOldConversion[i] = .0;

	for (index = 0; index < m_iNrTerms; ++index) 
	{
		rSum = 0;
		if (index > 0)
			cvConversion[index - 1] = 0;
		for (term = index; term < m_iNrTerms; ++term) {
			if (term == 0)
				cvConversion[term] = 1.0;
			else if (term == 1)
				if (index == 0)
					cvConversion[term] = m_rConstant;
				else
					cvConversion[term] = m_rMultiplier;
			else
				cvConversion[term] = 2 * m_rMultiplier * cvOldConversion[term - 1]
				+ 2 * m_rConstant * cvConversion[term - 1]
				- cvConversion[term - 2];
			rSum += cvConversion[term] * cvOldSolution[term];
		}
		m_cvSolution[index] = rSum;
		cvOldConversion = cvConversion;
	}

}

double ChebyshevFitter::rFuncPoly(int iTerm, double rX) const
{
	if (iTerm + 1 == m_iNrTerms)
		return m_cvSolution[iTerm];
	else
		return m_cvSolution[iTerm] + rX * rFuncPoly(iTerm + 1, rX);
	return 0;
}

double ChebyshevFitter::rCalc(double rX) const
{
	return rFuncPoly(0, rX);
}


vector<double> ChebyshevFitter::LeastSquaresFitPolynomial(const vector<double>& x,	const vector<double>& y, int n_terms)
{
	vector<double> coef = Chebyshev_fit(x, y, n_terms);

	vector<double> out;

	FinalizeSolution(coef);

	for (int i = 0; i < x.size(); i++)
		out.push_back(rCalc(x[i]));
	
	return out;
}


TableHortonPlots::TableHortonPlots()
{
}


TableHortonPlots::TableHortonPlots(quint64 metaid, const Ilwis::OperationExpression& expr) : OperationImplementation(metaid, expr)
{

}

bool TableHortonPlots::execute(ExecutionContext* ctx, SymbolTable& symTable)
{
	if (_prepState == sNOTPREPARED)
		if ((_prepState = prepare(ctx, symTable)) != sPREPARED)
			return false;

	bool resource = executeTableHortonPlots();

	if (resource && ctx != 0)
	{
		if (m_hortonTable.isValid())
		{
			Table* tbl = static_cast<Table*>(m_hortonTable->clone());
			QVariant value;
			value.setValue<ITable>(tbl);
			ctx->addOutput(symTable, value, tbl->name(), itTABLE, tbl->resource());
		}

		if (m_ratioTable.isValid())
		{
			Table* tbl = static_cast<Table*>(m_ratioTable->clone());
			QVariant value;
			value.setValue<ITable>(tbl);
			ctx->addOutput(symTable, value, tbl->name(), itTABLE, tbl->resource());
		}
	}

	return resource;
}


Ilwis::OperationImplementation* TableHortonPlots::create(quint64 metaid, const Ilwis::OperationExpression& expr)
{
	return new TableHortonPlots(metaid, expr);
}


quint64 TableHortonPlots::createMetadata()
{
	OperationResource operation({ "ilwis://operations/TableHortonPlots" });
	operation.setSyntax("TableHortonPlots(CatchmentMap,MergedCatchmentMap,DrainageNetworkOrderingMap");
	operation.setDescription(TR("Calculate for each (Strahler) stream order number and for each merged catchment: the number of streams, the average stream length, the average area of catchments, and expected values by means of least squares fit"));
	operation.setInParameterCount({ 3 });
	operation.addInParameter(0, itRASTER, TR("catchment Map"), TR("input raster that is the output of a Catchment extraction operation")); 
	operation.addInParameter(1, itRASTER, TR("merged catchment map"), TR("input raster that is the output of a Catchment merge operation"));
	operation.addInParameter(2, itRASTER, TR("drainage networ ordering map"), TR("input raster that is the result of a Drainage network ordering operation"));

	operation.parameterNeedsQuotes(1);
	operation.setOutParameterCount({ 2 });
	operation.addOutParameter(0, itTABLE, TR("output table"), TR("output table which can be used to construct so-called Horton plots"));
	operation.addOutParameter(1, itTABLE, TR("output table"), TR("output thble which contains the Bifurcation (Rb), Length (Rl), and Area (Ra) ratios"));
	operation.setKeywords("catchment,catchment merge,dranage ordering,horton table");
		
	operation.checkAlternateDefinition();
	mastercatalog()->addItems({ operation });
	return operation.id();

}


Ilwis::OperationImplementation::State TableHortonPlots::prepare(ExecutionContext* ctx, const SymbolTable& st)
{
	OperationImplementation::prepare(ctx, st);
	QString incatchment = _expression.parm(0).value();
	QString incatchmentmerge = _expression.parm(1).value();
	QString indrainage = _expression.parm(2).value();

	QString outputName1 = _expression.parm(0, false).value();
	QString outputName2 = _expression.parm(1, false).value();


	if (!m_catchmentMap.prepare(incatchment, itRASTER)) {
		ERROR2(ERR_COULD_NOT_LOAD_2, incatchment, "");
		return sPREPAREFAILED;
	}

	if (!m_drainageMap.prepare(indrainage, itRASTER)) {
		ERROR2(ERR_COULD_NOT_LOAD_2, indrainage, "");
		return sPREPAREFAILED;
	}

	if (!m_catchMergeMap.prepare(incatchmentmerge, itRASTER)) {
		ERROR2(ERR_COULD_NOT_LOAD_2, incatchmentmerge, "");
		return sPREPAREFAILED;
	}

	if (!m_ratioTable.prepare(outputName1, itFLATTABLE)) {
		ERROR1(ERR_NO_INITIALIZED_1, "output table");
		return sPREPAREFAILED;
	}

	if (!m_hortonTable.prepare(outputName2, itFLATTABLE)) {
		ERROR1(ERR_NO_INITIALIZED_1, "output table");
		return sPREPAREFAILED;
	}


	_xsize = m_catchmentMap->size().xsize();
	_ysize = m_catchmentMap->size().ysize();

	return sPREPARED;
}


bool TableHortonPlots::executeTableHortonPlots()
{
	PixelIterator iterDEM = PixelIterator(m_catchmentMap, BoundingBox(), PixelIterator::fXYZ);

	//Create columns in the output table
	ITable tblAtt = m_catchMergeMap->attributeTable();
	std::vector<QVariant> colDrainageID = tblAtt->column("DrainageID");

	int iSize = colDrainageID.size();

	CreateTable(iSize);
	long iMaxOrderNumber = 0;

	//check if column DrainageLen is valid in the input attribute table of the catchment
	//this column was calculared during catchemnt merging process using outlet locations
	//it is used to store the drainage length from/starting at the outlet coordinate!  
	bool fDrainageLen = false;

	std::vector<QVariant> colDrainageLen = tblAtt->column("DrainageLen");
	 if (colDrainageLen.size())
		 fDrainageLen = true;

	 int catchnum;
	 for (long i = 0; i < iSize; i++) 
	 {
		 std::vector<long> vDrainageIDs;
		 vDrainageIDs.clear();
		 SplitString(colDrainageID[i].toString(), QString(","), vDrainageIDs);

		 std::vector<double> vDrainageLen;
		 vDrainageLen.clear();
		 SplitString(colDrainageLen[i].toString(), QString(","), vDrainageLen);

		 std::vector<long> vOrder = GetOrderNumber(vDrainageIDs);
		 //returns the position of the first element with max. in vOrder
		 std::vector<long>::iterator pos = max_element(vOrder.begin(), vOrder.end());
		 long iMax = *pos;
		 catchnum = i + 1;
		 ComputeHortonPlots(catchnum, iMax, vDrainageIDs, vOrder, fDrainageLen, vDrainageLen);
		 if (iMax > iMaxOrderNumber)
			 iMaxOrderNumber = iMax;
	
	 }


	return true;
}


////Compute parameters for Horton Plots, and
////also store values in the ouput table
void TableHortonPlots::ComputeHortonPlots(long CatchmentNr, long iMaxOrder, std::vector<long> vDrainageIDs, std::vector<long> vOrder, bool fDrainageLen, std::vector<double> vDrainageLen)
{
	////Retrieve drainage attributes from related table 
	ITable tblAtt = m_drainageMap->attributeTable();
	std::vector<QVariant> colLength = tblAtt->column("Length");

	////Retrieve catchment attributes from the related table 
	ITable tblCam = m_catchmentMap->attributeTable();
	std::vector<QVariant> colArea = tblCam->column("CatchmentArea");

	////Compute the number of channels of order starting from 1 to iMaxOrder
	QString sCol1 = QString(TR("C%1_N")).arg(CatchmentNr);
	QString sCol2 = QString(TR("C%1_L")).arg(CatchmentNr);
	QString sCol3 = QString(TR("C%1_A")).arg(CatchmentNr);
	QString sCol4 = QString(TR("Order"));
	QString sCol5 = QString(TR("C%1_N_LSq")).arg(CatchmentNr);
	QString sCol6 = QString(TR("C%1_L_Sq")).arg(CatchmentNr);
	QString sCol7 = QString(TR("C%1_A_LSq")).arg(CatchmentNr);

	QString sCol8 = QString(TR("C%1_RB")).arg(CatchmentNr);
	QString sCol9 = QString(TR("C%1_RL")).arg(CatchmentNr);
	QString sCol10 = QString(TR("C%1_RA")).arg(CatchmentNr);

	QString sCol5b = QString(TR("C%1_N_LSq_b")).arg(CatchmentNr);
	QString sCol6b = QString(TR("C%1_L_LSq_b")).arg(CatchmentNr);
	QString sCol7b = QString(TR("C%1_A_LSq_b")).arg(CatchmentNr);
	QString sCol8b = QString(TR("C%1_RB_b")).arg(CatchmentNr);
	QString sCol9b = QString(TR("C%1_RL_b")).arg(CatchmentNr);
	QString sCol10b = QString(TR("C%1_RA_b")).arg(CatchmentNr);

	std::vector<QVariant> col5;
	std::vector<QVariant> col8;

	std::vector<QVariant> col5b;
	std::vector<QVariant> col8b;

	std::vector<QVariant> col6;
	std::vector<QVariant> col9;

	std::vector<QVariant> col6b;
	std::vector<QVariant> col9b;

	std::vector<QVariant> col7;
	std::vector<QVariant> col7b;

	std::vector<QVariant> col10;
	std::vector<QVariant> col10b;

	m_hortonTable->addColumn(sCol5, IlwisObject::create<IDomain>("value"), true);
	ColumnDefinition& coldef5 = m_hortonTable->columndefinitionRef(sCol5);
	coldef5.datadef().range(new NumericRange(0, 1.0e300, 0.01));

	m_hortonTable->addColumn(sCol6, IlwisObject::create<IDomain>("value"), true);
	ColumnDefinition& coldef6 = m_hortonTable->columndefinitionRef(sCol6);
	coldef6.datadef().range(new NumericRange(0, 1.0e300, 0.01));


	m_hortonTable->addColumn(sCol7, IlwisObject::create<IDomain>("value"), true);
	ColumnDefinition& coldef7 = m_hortonTable->columndefinitionRef(sCol7);
	coldef7.datadef().range(new NumericRange(0, 1.0e300, 0.01));



	////total upstream length and area 
	double rTotalLength = 0;
	double rTotalArea = 0;
	double rPreArea = 0;

	for (long i = 1; i <= iMaxOrder; i++)
	{
		int rec = i - 1;
		m_hortonTable->setCell(sCol4, rec, QVariant(i));
		long iCout = std::count(vOrder.begin(), vOrder.end(), i);
		long iNrOrder = 0;

		if (iCout > 0)
		{
			//computer length of streams of each order
			//find the first element with given value
			double rLength = 0;
			double rArea = 0;
			std::vector<long>::iterator pos = find(vOrder.begin(), vOrder.end(), i);
			int iIndex = pos - vOrder.begin();

			

			while (pos != vOrder.end())
			{
				long id = vDrainageIDs[iIndex]-1;
				if (i == 1)
					iNrOrder = iCout;
				else
				{
					if (!fUpstreamSameOrder(id, i, vDrainageIDs))
						iNrOrder++;
				}
				if (fDrainageLen)
					if (vDrainageLen[iIndex] > 0.0)
					{
						rLength += vDrainageLen[iIndex];
						rArea += colArea[id].toDouble();
					}
					else
					{
						rLength += colLength[id].toDouble();
						rArea += colArea[id].toDouble();
					}
				else
				{
					rLength += colLength[id].toDouble();
					rArea += colArea[id].toDouble();
				}
				pos = find(++pos, vOrder.end(), i);
				iIndex = pos - vOrder.begin();
			}
			m_hortonTable->setCell(sCol1, rec, QVariant(iNrOrder));
			
			rLength = (rLength / iNrOrder) / 1000;
			rArea = (rArea / iNrOrder) / 1000000;

			rTotalLength = rTotalLength + rLength;
			rTotalArea = (rPreArea + rArea * iNrOrder) / iNrOrder;
			rPreArea = rPreArea + rArea * iNrOrder;

			m_hortonTable->setCell(sCol2, rec, QVariant(rTotalLength));
			m_hortonTable->setCell(sCol3, rec, QVariant(rTotalArea));
		}
		else
		{
			iMaxOrder = i - 1;
			break;
		}
	}
	//Compute regression coefficients
	std::vector<QVariant> col1 = m_hortonTable->column(sCol1);
	std::vector<QVariant> col2 = m_hortonTable->column(sCol2);
	std::vector<QVariant> col3 = m_hortonTable->column(sCol3);
	std::vector<QVariant> col4 = m_hortonTable->column("Order");


	if (iMaxOrder > 2)
	{
		//Calculate Rb
		Regress(CatchmentNr, col1, col4, col5, iMaxOrder, false);
		CalculateRB(col5, col8, true);
		Regress(CatchmentNr, col1, col4, col5b, iMaxOrder, true); //using all orders  
		CalculateRB(col5b, col8b, true);

		//Calculate Rl
		Regress(CatchmentNr, col2, col4, col6, iMaxOrder, false);
		CalculateRB(col6, col9, false);
		Regress(CatchmentNr, col2, col4, col6b, iMaxOrder, true); //using all orders
		CalculateRB(col6b, col9b, false);

		//Calculate Ra
		Regress(CatchmentNr, col3, col4, col7, iMaxOrder, false);
		CalculateRB(col7, col10, false);
		Regress(CatchmentNr, col3, col4, col7b, iMaxOrder, true); //using all orders
		CalculateRB(col7b, col10b, false);


		for (long i = 1; i <= iMaxOrder; i++)
		{
			int rec = i - 1;
			m_hortonTable->setCell(sCol5, rec, QVariant(col5[rec]));
			m_hortonTable->setCell(sCol6, rec, QVariant(col6[rec]));
			m_hortonTable->setCell(sCol7, rec, QVariant(col7[rec]));
		}

		QString sCol_a = QString(TR("C%1_a")).arg(CatchmentNr);
		QString sCol_b = QString(TR("C%1_b")).arg(CatchmentNr);

		m_ratioTable->addColumn(sCol_a, IlwisObject::create<IDomain>("value"), true);
		ColumnDefinition& coldefa = m_ratioTable->columndefinitionRef(sCol_a);
		coldefa.datadef().range(new NumericRange(0, 1.0e300, 0.01));

		m_ratioTable->addColumn(sCol_b, IlwisObject::create<IDomain>("value"), true);
		ColumnDefinition& coldefb = m_ratioTable->columndefinitionRef(sCol_b);
		coldefb.datadef().range(new NumericRange(0, 1.0e300, 0.01));


		m_ratioTable->setCell(sCol_a, 0, QVariant(col8[0]));
		m_ratioTable->setCell(sCol_a, 1, QVariant(col9[0]));
		m_ratioTable->setCell(sCol_a, 2, QVariant(col10[0]));

		m_ratioTable->setCell(sCol_b, 0, QVariant(col8b[0]));
		m_ratioTable->setCell(sCol_b, 1, QVariant(col9b[0]));
		m_ratioTable->setCell(sCol_b, 2, QVariant(col10b[0]));



		////	//Extrat ratios to a seperate table
		//QString sCol_a = QString(TR("C%li_a")).arg(CatchmentNr);
		//QString sCol_b = QString(TR("C%li_b")).arg(CatchmentNr);
		//std::vector<QVariant> cC_a = m_ratioTable->column(sCol_a);
		//std::vector<QVariant> cC_b = m_ratioTable->column(sCol_b);

		//double val;
		//long iRaw;

		////Extract Ra 
		//double col8val = col8[1].toDouble();
		//double col8bval = col8b[1].toDouble();
		//m_ratioTable->setCell(sCol_a, 1, QVariant(col8val));
		//m_ratioTable->setCell(sCol_b, 1, QVariant(col8bval));

		////Extract Rl 
		//double col9val = col9[1].toDouble();
		//double col9bval = col9b[1].toDouble();
		//m_ratioTable->setCell(sCol_a, 2, QVariant(col9val));
		//m_ratioTable->setCell(sCol_b, 2, QVariant(col9bval));

		////Extract Ra

		//double col10val = col10[1].toDouble();
		//double col10bval = col10b[1].toDouble();
		//m_ratioTable->setCell(sCol_a, 3, QVariant(col10val));
		//m_ratioTable->setCell(sCol_b, 3, QVariant(col10bval));
	}

}

bool TableHortonPlots::fUpstreamSameOrder(long id, long iOrder, std::vector<long> vDrainageIDs)
{
	ITable tblAtt = m_drainageMap->attributeTable();
	std::vector<QVariant> cUpstreamLink = tblAtt->column("UpstreamLinkID");
	std::vector<QVariant> cOrder = tblAtt->column("Strahler");

	QString sLinks = cUpstreamLink[id].toString();
	std::vector<long> vLinks;
	SplitString(sLinks, QString(","), vLinks);


	for ( int i = 0; i< vLinks.size(); i++ )
	{
		int idx = vLinks[i];
		if (cOrder[idx-1].toInt() == iOrder)
		{
			long drainageID = idx;
			std::vector<long>::iterator pos = find(vDrainageIDs.begin(), vDrainageIDs.end(), drainageID);
			if (pos != vDrainageIDs.end())
				return true;
		}
	}


	return false;
}


void TableHortonPlots::Regress(long CatchmentNr, std::vector<QVariant> colY, std::vector<QVariant> colX, std::vector<QVariant>& outCol, int iOrders, bool usAllOrders)
{
	int sz = colY.size();
	std::vector<double> valY;
	std::vector<double> valX;
	std::vector<double> valNewY;
	std::vector<double> valLogNewY;
	std::vector<double> valfitNewY;

	for (int i = 0; i < sz; i++)
	{
		valY.push_back(colY[i].toDouble());
		valX.push_back(colX[i].toDouble());
		valNewY.push_back(colY[i].toDouble());
	}

	if (usAllOrders != true) {

		if (iOrders >= 5)
		{
			valNewY[0] = rUNDEF;
			valNewY[valNewY.size() - 1] = rUNDEF;
		}
		else if (iOrders == 4)
		{
			valNewY[valNewY.size() - 1] = rUNDEF;
		}
	}

	long iNrValidPnts = 0;
	for (int i = 0; i < sz; ++i) {
		if (valX[i] == rUNDEF || valNewY[i] == rUNDEF)
			continue;
		iNrValidPnts++;     // count all valid X,Y points
	}

	// log math

	for (int i = 0; i < valY.size(); i++)
	{
		if (valNewY[i] == rUNDEF)
			valLogNewY.push_back(valNewY[i]);
		else
			valLogNewY.push_back(log10(valNewY[i]));
	}

	///// below code can be replaced by the commentout ilwis3 code below
	//ChebyshevFitter fitter;
	//valfitNewY = fitter.LeastSquaresFitPolynomial(valX, valLogNewY, 2);

	//for (int i = 0; i < sz; ++i) {
	//	double rX = valX[i];
	//	if (rX == rUNDEF)
	//		outCol.push_back(rUNDEF); // no Y-value computed if X-value undefined
	//	else
	//		outCol.push_back(pow(10, valfitNewY[i]));
	//}

	// ilwis algroithm
	
	CVector cvX(iNrValidPnts), cvY(iNrValidPnts);
	long j = 0;
	for (int i = 0; i < sz; ++i) {
		if (valX[i] == rUNDEF || valLogNewY[i] == rUNDEF )
			continue;         // skip non_valid points
		cvX(j) = valX[i];
		cvY(j) = valLogNewY[i];
		j++;
	}

	LeastSquaresFitPolynomial lsf (cvX, cvY, 2);

	for (int i = 0; i < sz; ++i) {
		double rX = valX[i];
		if (rX == rUNDEF)
			outCol.push_back(rUNDEF); // no Y-value computed if X-value undefined
		else
			outCol.push_back(pow(10,lsf.rCalc(rX)));
	}
	
}


std::vector<long> TableHortonPlots::GetOrderNumber(std::vector<long> vIDs)
{
	std::vector<long> vOrder;
	ITable tblAtt = m_drainageMap->attributeTable();
	std::vector<QVariant> colStralhaer = tblAtt->column("Strahler");

	for (int pos = 0; pos < vIDs.size(); ++pos)
	{
		vOrder.push_back(colStralhaer[vIDs[pos]-1].toInt());
	}
	return vOrder;
}

void TableHortonPlots::SplitString(QString s, QString mid, std::vector<double>& results)
{
	s.replace("{", "");
	s.replace("}", "");
	QStringList strlst = s.split(mid);

	results.clear();
	for (unsigned int i = 0; i < strlst.size(); i++)
	{
		double res = strlst[i].toDouble();
		if (res != iUNDEF)
			results.push_back(res);
	}
}

void TableHortonPlots::SplitString(QString s, QString mid, std::vector<long>& results)
{
	s.replace("{", "");
	s.replace("}", "");
	QStringList strlst = s.split(mid);

	results.clear();
	for (unsigned int i = 0; i < strlst.size(); i++)
	{
		long res = strlst[i].toLong();
		if (res != iUNDEF)
			results.push_back(res);
	}
}


void TableHortonPlots::CreateTable(long iSize)
{
	IFlatTable newTable;
	newTable.prepare();
	newTable->addColumn("Order", IlwisObject::create<IDomain>("value"), true);
	ColumnDefinition& coldef0 = newTable->columndefinitionRef("Order");
	coldef0.datadef().range(new NumericRange(0, 32767, 1));

	for (long i = 1; i <= iSize; i++)
	{
		QString sCol1 = QString(TR("C%1_N")).arg(i);
		QString sCol2 = QString(TR("C%1_L")).arg(i);
		QString sCol3 = QString(TR("C%1_A")).arg(i);

		newTable->addColumn(sCol1, IlwisObject::create<IDomain>("value"), true);
		ColumnDefinition& coldef0 = newTable->columndefinitionRef(sCol1);
		coldef0.datadef().range(new NumericRange(0, 32767, 1));

		newTable->addColumn(sCol2, IlwisObject::create<IDomain>("value"), true);
		ColumnDefinition& coldef1 = newTable->columndefinitionRef(sCol2);
		coldef1.datadef().range(new NumericRange(0, 1.0e300, 0.01));

		newTable->addColumn(sCol3, IlwisObject::create<IDomain>("value"), true);
		ColumnDefinition& coldef2 = newTable->columndefinitionRef(sCol3);
		coldef2.datadef().range(new NumericRange(0, 1.0e300, 0.01));
	}

	m_hortonTable = newTable;

	IFlatTable ratioTable;
	ratioTable.prepare();
	QString itemname("ItemName");
	ratioTable->addColumn(itemname, IlwisObject::create<IDomain>("text"), true);
	
	for (long i = 1; i <= iSize; i++)
	{
		QString sCol1 = QString(TR("C%1_a")).arg(i);
		QString sCol2 = QString(TR("C%1_b")).arg(i);

		ratioTable->addColumn(sCol1, IlwisObject::create<IDomain>("value"), true);
		ColumnDefinition& coldef10 = ratioTable->columndefinitionRef(sCol1);
		coldef10.datadef().range(new NumericRange(0, 1.0e300, 0.01));

		ratioTable->addColumn(sCol2, IlwisObject::create<IDomain>("value"), true);
		ColumnDefinition& coldef11 = ratioTable->columndefinitionRef(sCol2);
		coldef11.datadef().range(new NumericRange(0, 1.0e300, 0.01));
	}

	m_ratioTable = ratioTable;

	m_ratioTable->setCell(itemname, 0, QVariant(QString("Rb")));
	m_ratioTable->setCell(itemname, 1, QVariant(QString("Rl")));
	m_ratioTable->setCell(itemname, 2, QVariant(QString("Ra")));


}


double TableHortonPlots::rDistance(Coordinate cd1, Coordinate cd2)
{
	double rDist;
	if (fLatLonCoords())
	{
		double rRadi = rDefaultEarthRadius;
		IConventionalCoordinateSystem projectedCsy = m_catchmentMap->coordinateSystem().as<ConventionalCoordinateSystem>();
		if (projectedCsy->isValid())
			rRadi = projectedCsy->ellipsoid()->majorAxis();

		LatLon llStart = LatLon(cd1.y, cd1.x);
		LatLon llEnd = LatLon(cd2.y, cd2.x);
		if (projectedCsy->ellipsoid()->isSpherical())
		{
			if (llStart.Lat() == llEnd.Lat() && llStart.Lon() == llEnd.Lon())
				rDist = 0; //seems a bug in rEllipsoidalDistance, always get some value, even when llStart and llEnd the same? 
			else
				rDist = projectedCsy->ellipsoid()->distance(llStart, llEnd);
			if (rDist < 8000)
				return rDist;
		}

		rDist = rSphericalDistance(rRadi, llStart, llEnd);
	}
	else
	{
		double dx = (cd1.x - cd2.x);
		double dy = (cd1.y - cd2.y);
		rDist = sqrt(dx * dx + dy * dy);
	}
	return rDist;
}


bool TableHortonPlots::fLatLonCoords()
{
	return m_catchmentMap->coordinateSystem()->isLatLon();
}


void TableHortonPlots::CalculateRB(std::vector<QVariant> col, std::vector<QVariant>& sCol, bool isRB)
{
	for (int i = 0; i < col.size(); i++)
	{
		if ( i == col.size()-1 || col[i].toDouble() == rUNDEF || col[i+1].toDouble() == rUNDEF )
		{
			sCol.push_back(rUNDEF);
			continue;
		}

		if (isRB)
			sCol.push_back(col[i].toDouble() / col[i + 1].toDouble());
		else
			sCol.push_back(col[i + 1].toDouble() / col[i].toDouble());
	}

}
