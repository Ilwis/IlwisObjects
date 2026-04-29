
/* IlwisObjects is a framework for analysis, processingand visualization of remote sensingand gis data
Copyright(C) 2018  52n North

This program is free software : you can redistribute it and /or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see < http://www.gnu.org/licenses/>.*/

#include <functional>
#include <future>
#include "kernel.h"
#include "raster.h"
#include "symboltable.h"
#include "ilwisoperation.h"
#include "blockiterator.h"
#include "pixeliterator.h"
#include "MapFlowDirection.h"
#include "itemdomain.h"
#include "thematicitem.h"



using namespace Ilwis;
using namespace Hydroflow;

REGISTER_OPERATION(MapFlowDirection)

MapFlowDirection::MapFlowDirection()
{
}

MapFlowDirection::MapFlowDirection(quint64 metaid, const Ilwis::OperationExpression& expr) : OperationImplementation(metaid, expr)
{
}

bool MapFlowDirection::execute(ExecutionContext* ctx, SymbolTable& symTable)
{
    if (_prepState == sNOTPREPARED)
        if ((_prepState = prepare(ctx, symTable)) != sPREPARED)
            return false;

	executeFlowDirection();

    bool resource = true;
    if (resource && ctx != 0) {
        QVariant value;
        value.setValue<IRasterCoverage>(_outRaster);
        ctx->setOutput(symTable, value, _outRaster->name(), itRASTER, _outRaster->resource());
    }
    return resource;
}

Ilwis::OperationImplementation* MapFlowDirection::create(quint64 metaid, const Ilwis::OperationExpression& expr)
{
    return new MapFlowDirection(metaid, expr);
}

Ilwis::OperationImplementation::State MapFlowDirection::prepare(ExecutionContext* ctx, const SymbolTable& st)
{
	m_fParallel = false;
	m_fmMethods = fmSlope;

	OperationImplementation::prepare(ctx, st);
    QString inraster = _expression.parm(0).value();
    QString outputName = _expression.parm(0, false).value();
    QString methodstr = _expression.parm(1).value().toLower();
	QString parallelstr = _expression.parm(2).value().toLower();

    std::map<QString, FlowMethod> methods = { {"slope",fmSlope},{"height",fmHeight} };

    auto iter = methods.find(methodstr);
    if (iter == methods.end())
    {
        ERROR2(ERR_NOT_FOUND2, methodstr, TR("in method for flow dirction"));
        return sPREPAREFAILED;
    }
    m_fmMethods = iter->second;

	if (parallelstr == "yes" || parallelstr == "1")
		m_fParallel = true;
	else
		m_fParallel = false;


    if (!_inRaster.prepare(inraster, itRASTER)) {
        ERROR2(ERR_COULD_NOT_LOAD_2, inraster, "");
        return sPREPAREFAILED;
    }

    int copylist = itRASTERSIZE | itENVELOPE | itCOORDSYSTEM | itGEOREF;
    _outRaster = OperationHelperRaster::initialize(_inRaster.as<IlwisObject>(), itRASTER, copylist);
    if (!_outRaster.isValid()) {
        ERROR1(ERR_NO_INITIALIZED_1, "output rastercoverage");
        return sPREPAREFAILED;
    }

	IThematicDomain flowdirectionDom;
	flowdirectionDom.prepare();
	flowdirectionDom->addItem(new ThematicItem({ "E", "", "To the east", }, 1));
	flowdirectionDom->addItem(new ThematicItem({ "SE", "", "To the south east", }, 2));
	flowdirectionDom->addItem(new ThematicItem({ "S", "", "To the south", }, 3));
	flowdirectionDom->addItem(new ThematicItem({ "SW", "", "To the south west", }, 4));
	flowdirectionDom->addItem(new ThematicItem({ "W", "", "To the west", }, 5));
	flowdirectionDom->addItem(new ThematicItem({ "NW", "", "To the north west", }, 6));
	flowdirectionDom->addItem(new ThematicItem({ "N", "", "To the north", }, 7));
	flowdirectionDom->addItem(new ThematicItem({ "NE", "", "To the north east", }, 8));
	flowdirectionDom->name("FlowDirection");
	_outRaster->datadefRef() = DataDefinition(flowdirectionDom);

    for (quint32 i = 0; i < _outRaster->size().zsize(); ++i) {
        QString index = _outRaster->stackDefinition().index(i);
        _outRaster->setBandDefinition(index, DataDefinition(flowdirectionDom));
    }

    if (outputName != sUNDEF)
        _outRaster->name(outputName);

    PixelIterator iterDEM = PixelIterator(_inRaster, BoundingBox(), PixelIterator::fXYZ);
    PixelIterator iterFlow = PixelIterator(_outRaster, BoundingBox(), PixelIterator::fXYZ);
    PixelIterator inEnd = iterDEM.end();

    
	std::fill(iterFlow, iterFlow.end(), 0);

    _xsize = _outRaster->size().xsize();
    _ysize = _outRaster->size().ysize();

    // initialize tranquilizer
    initialize(_xsize * _ysize);

    return sPREPARED;
}

quint64 MapFlowDirection::createMetadata()
{
	OperationResource operation({ "ilwis://operations/MapFlowDirection" });
	operation.setSyntax("MapFlowDirection(inputraster,method=slope|height,useparalleldrainagecorrection=yes|no)");
	operation.setDescription(TR("generates a new raster containing flow directions"));
	operation.setInParameterCount({ 3 });
	operation.addInParameter(0, itRASTER, TR("rastercoverage"), TR("input a sink-free DEM with numeric domain"));
	operation.addInParameter(1, itSTRING, TR("method definition"), TR("Flow direction should be calculated according to the steepst slope or the smallest height"), OperationResource::ueCOMBO);
	operation.addInParameter(2, itSTRING, TR("Parallel drainage correction algorithm"), TR("Option of flow direction algorithm."), OperationResource::ueNONE);
	operation.parameterNeedsQuotes(1);
	operation.setOutParameterCount({ 1 });
	operation.addOutParameter(0, itRASTER, TR("output raster"), TR("output raster with a numeric domain"));
	operation.setKeywords("flow direction,raster,image processing, numeric");

	operation.checkAlternateDefinition();
	mastercatalog()->addItems({ operation });
	return operation.id();

}


void MapFlowDirection::executeFlowDirection()
{
	InitPars();

	m_inDem = new double[_xsize * _ysize];
	m_flowDem = new double[_xsize * _ysize];
	m_flag = new int[_xsize * _ysize];

	for (int y = 0; y < _ysize; y++)
	{
		for (int x = 0; x < _xsize; x++) 
		{
			double val = _inRaster->pix2value(Pixel(x, y, 0)); // replace with access
			m_inDem[idx(x, y)] = val;
			m_flowDem[idx(x, y)] = val;
			m_flag[idx(x, y)] = 0;
		}
	}

	if (m_fParallel) 
	{
		FlowDirectionAlgorithm fda(m_inDem, m_flowDem, _xsize,_ysize);
		fda.calculate(m_fmMethods == fmSlope ? "slope" : "height");
	}
	else 
	{
		for (int y = 0; y < _ysize; ++y) 
		{
			for (int x = 0; x < _xsize; ++x)
			{
				int offset = idx(x,y);

				std::vector<double> vValue;
				FillArray(x, y, vValue);

				std::vector<int> vPos;
				int iCout = 0;
				double rMax = rFindMaxLocation(vValue, vPos, iCout);

				m_flowDem[offset] = iLookUp(rMax, iCout, vPos);
			}
		}
		TreatFlatAreas();
	}

	PixelIterator iterFlow = PixelIterator(_outRaster, BoundingBox(), PixelIterator::fXYZ);
	// cleanup edges and rewrite raster pixel
	for (int y = 0; y < _ysize; y++) 
	{
		for (int x = 0; x < _xsize; x++)
		{

			if (onEdge(x, y) || m_flowDem[idx(x, y)] > 8)
				m_flowDem[idx(x, y)] = 0;

			Pixel pxl(x, y, 0);
			*iterFlow(pxl) = (m_flowDem[idx(x, y)] == 0) ? rUNDEF : (m_flowDem[idx(x,y)] - 1);

		}
	}
	
	delete[]m_inDem;
	delete[]m_flowDem;
	delete[]m_flag;

}


double MapFlowDirection::rFindMaxLocation(std::vector<double>& vValue, std::vector<int>& vPos, int& iCout)
{
	//finds the maximum value in the input vector
		//returns the maximum value, number of elements with max.
		//returns posision(s) for the element(s) with max. in a vector

	std::vector<double>::iterator pos;

	//returns the position of the first element with max. in vValue
	pos = max_element(vValue.begin(), vValue.end());
	double rMax = *pos;

	//count number of elements with max
	iCout = std::count(vValue.begin(), vValue.end(), rMax);

	//find the first element with max value
	pos = find(vValue.begin(), vValue.end(), rMax);
	int iIndex = pos - vValue.begin();

	while (pos != vValue.end())
	{
		vPos.push_back(iIndex);  //push it into a vector
		pos = find(++pos, vValue.end(), rMax);
		iIndex = pos - vValue.begin();
	}
	return rMax;
}

bool MapFlowDirection::onEdge(int x, int y)
{
	return x == 0 || x == _xsize - 1 ||
		y == 0 || y == _ysize - 1;
}


void MapFlowDirection::FillArray(int x, int y, std::vector<double>& vValue)
{
	double rCurH = m_inDem[idx(x, y)];

	int dx[8] = { 1,1,0,-1,-1,-1,0,1 };
	int dy[8] = { 0,1,1,1,0,-1,-1,-1 };

	for (int i = 0; i < 8; i++) 
	{

		int nx = x + dx[i];
		int ny = y + dy[i];

		double rNbH = m_inDem[idx(nx, ny)];

		if (rNbH == rUNDEF) 
		{
			vValue.push_back(rUNDEF);
		}
		else 
		{
			double rValue = (m_fmMethods == fmSlope)
				? rComputeSlope(rCurH, rNbH, i + 1)
				: rComputeHeightDifference(rCurH, rNbH);

			vValue.push_back(rValue);
		}
	}
}


void MapFlowDirection::SetFlowsInFlatArea(std::vector<int>& vOutlets)
{
	std::vector<int> vNbs; // Neighbors without flow direction

	do
	{
		vNbs.clear();

		// Use a standard index loop or iterator for the offsets
		for (int offset : vOutlets)
		{
			m_flag[offset] = 0;
		}

		for (int offset : vOutlets)
		{
			int currX = offset % _xsize;
			int currY = offset / _xsize;

			int iNb = 0;
			// Neighborhood iteration (3x3 grid)
			for (int i = -1; i <= 1; ++i)
			{
				int nbY = currY + i;
				for (int j = -1; j <= 1; ++j)
				{
					int nbX = currX + j;
					int nbOffset = nbY * _xsize + nbX;

					// Bounds check
					if (nbX >= 0 && nbX < _xsize && nbY >= 0 && nbY < _ysize)
					{
						// Logic for flat area flow assignment
						if (m_flowDem[nbOffset] == 9 &&
							m_vFlowSelf[iNb] != m_flowDem[offset] &&
							m_flag[nbOffset] == m_ContFlat)
						{
							m_flowDem[nbOffset] = m_vDirection[iNb];
							vNbs.push_back(nbOffset);
						}

						if ((m_flag[nbOffset] == m_ContFlat) &&
							!(isEven(iNb)) &&
							m_vFlowSelf[iNb] != m_flowDem[offset])
						{
							m_flowDem[nbOffset] = m_vDirection[iNb];
						}
					}
					iNb++;
				}
			}
		}
		vOutlets.swap(vNbs);
	} while (!vOutlets.empty());
}


void MapFlowDirection::TreatFlatAreas()
{
	size_t totalPixels = (size_t)_xsize * _ysize;
	std::fill(m_flag, m_flag + totalPixels, 0);
	m_ContFlat = 0;

	for (int y = 0; y < _ysize; ++y)
	{
		for (int x = 0; x < _xsize; ++x) 
		{

			if (!onEdge(x, y) && m_flowDem[idx(x,y)] == 9) 
			{
				std::vector<int> vOutlets;
				m_vFlatIndices.clear(); // Store offsets instead of Pixels

				LocateOutlets(x, y, vOutlets);

				if (vOutlets.size() > 0) 
				{
					SetFlowsInFlatArea(vOutlets);
				}
				else 
				{
					for (int flatIdx : m_vFlatIndices)
					{
						m_flowDem[flatIdx] = 0;
					}
				}
			}
		}
	}
}

void MapFlowDirection::LocateOutlets(int x, int y, std::vector<int>& vOutlets)
{
	std::vector<int> vStarters;
	int startOffset = y * _xsize + x;
	vStarters.push_back(startOffset);

	double rHeight = m_inDem[startOffset];
	m_ContFlat++;

	while (!vStarters.empty())
	{
		std::vector<int> vNbs;
		for (int currentOffset : vStarters) 
		{
			m_flag[currentOffset] = m_ContFlat;
			m_vFlatIndices.push_back(currentOffset);

			int currX = currentOffset % _xsize;
			int currY = currentOffset / _xsize;

			for (int i = -1; i <= 1; ++i) 
			{
				for (int j = -1; j <= 1; ++j) 
				{
					if (i == 0 && j == 0) continue;

					int nbX = currX + j;
					int nbY = currY + i;
					int nbOffset = nbY * _xsize + nbX;

					if (nbX < 0 || nbX >= _xsize || nbY < 0 || nbY >= _ysize) continue;

					if (m_flag[nbOffset] != m_ContFlat) 
					{
						if (m_flowDem[nbOffset] == 9) 
						{
							m_flag[nbOffset] = m_ContFlat;
							m_vFlatIndices.push_back(nbOffset);
							vNbs.push_back(nbOffset);
						}
						else if (m_inDem[nbOffset] == rHeight && m_flowDem[nbOffset] != 0) 
						{
							m_flag[nbOffset] = m_ContFlat;
							vOutlets.push_back(nbOffset);
						}
					}
				}
			}
		}
		vStarters = std::move(vNbs);
	}
}


double MapFlowDirection::rComputeSlope(double rCurH, double rNbH, int iPos)
{
		//The slope is calculated by subscribing the neighbor's value from the center
		//distance 1.14 is concerned for diagonal cells
		double rVal;
		if (isEven(iPos))
			rVal = (rCurH - rNbH)/1.41;
		else
			rVal = rCurH - rNbH;
		return rVal;
}

double MapFlowDirection::rComputeHeightDifference(double rCurH, double rNbH)
{
		return rCurH - rNbH;
}

bool MapFlowDirection::isEven(int elem)
{
		return elem % 2 == 0;
}


bool MapFlowDirection::isInOneEdge(int iPos1, int iPos2, int iPos3, std::vector<int>& vPos)
{
		bool fCondition1 = find(vPos.begin(),vPos.end(),iPos1) != vPos.end();
		bool fCondition2 = find(vPos.begin(),vPos.end(),iPos2) != vPos.end();
		bool fCondition3 = find(vPos.begin(),vPos.end(),iPos3) != vPos.end();
		return fCondition1 && fCondition2 && fCondition3;
}


void MapFlowDirection::InitPars()
{
		//	Location number				Order in m_vDirection
		//	-------								-------	 looping order of the neighbors
	  //	|6|7|8|								|0|1|2|
		//	-------								-------
		//	|5| |1|								|3|4|5|
		//	-------								-------
		//	|4|3|2|								|6|7|8|
		//	-------								-------
		//
		m_vDirection.resize(9);

		m_vDirection[0] = 2;
		m_vDirection[1] = 3;
		m_vDirection[2] = 4;
		m_vDirection[3] = 1;
		m_vDirection[4] = 0;
		m_vDirection[5] = 5;
		m_vDirection[6] = 8;
		m_vDirection[7] = 7;
		m_vDirection[8] = 6;

		m_vFlowSelf.resize(9);
		m_vFlowSelf[0] = 6;
		m_vFlowSelf[1] = 7;
		m_vFlowSelf[2] = 8;
		m_vFlowSelf[3] = 5;
		m_vFlowSelf[4] = 0;
		m_vFlowSelf[5] = 1;
		m_vFlowSelf[6] = 4;
		m_vFlowSelf[7] = 3;
		m_vFlowSelf[8] = 2;
}

FlowDirectionAlgorithm::Method FlowDirectionAlgorithm::methodValueOf(QString val) {
	if (val== QString("slope")) {
		return FlowDirectionAlgorithm::slope;
	}
	return FlowDirectionAlgorithm::height;
}


FlowDirectionAlgorithm::FlowDirectionAlgorithm(double* dem, double* flow, long xsz, long ysz)
	: flatcell(9)
	, flag(10)
	, increment(1)
	, _xsize(xsz)
	, _ysize(ysz)
	, m_inDem(dem)
	, m_flowDem(flow)
{
	//	Location number				Order in m_vDirection
	//	-------						-------	 looping order of the neighbors 	
	//	|6|7|8|						|0|1|2|
	//	-------						-------
	//	|5| |1|						|3|4|5|
	//	-------						-------
	//	|4|3|2|						|6|7|8|
	//	-------						-------
	//

	Location[0] = 6;
	Location[1] = 7;
	Location[2] = 8;
	Location[3] = 5;
	Location[4] = 1;
	Location[5] = 4;
	Location[6] = 3;
	Location[7] = 2;
}


void FlowDirectionAlgorithm::calculate(QString methodInput)
{
	QString sl = methodInput.toLower();
	method = methodValueOf(sl);

	// -------------------------------
	// STEP 1: Initial flow direction
	// -------------------------------
	for (int y = 0; y < _ysize; y++) {
		for (int x = 0; x < _xsize; x++) {

			int idx = y * _xsize + x;
			m_flowDem[idx] = noflow;

			if (onEdge(x, y) || m_inDem[idx] == rUNDEF)
				continue;

			double listVal[8];
			double max = maxAdj(x, y, listVal);

			if (max > 0) {
				std::vector<FlowDirection> listPos;
				findDirection(listVal, max, listPos);
				FlowDirection fd = getFlowDirection(listPos);
				m_flowDem[idx] = Location[(unsigned char)fd];
			}
			else if (max == rUNDEF) {
				m_flowDem[idx] = noflow;
			}
			else {
				m_flowDem[idx] = flatcell;
			}
		}
	}

	// -------------------------------
	// STEP 2: allocate gradients
	// -------------------------------
	double* grad1 = new double[_xsize * _ysize];
	double* grad2 = new double[_xsize * _ysize];

	std::fill(grad1, grad1 + _xsize * _ysize, 0.0);
	std::fill(grad2, grad2 + _xsize * _ysize, 0.0);

	// -------------------------------
	// STEP 3: Flat area processing
	// -------------------------------
	for (int y = 0; y < _ysize; y++) {
		for (int x = 0; x < _xsize; x++) {

			//int idx = y * _xsize + x;

			if (m_flowDem[idx(x,y)] == flatcell) {

				std::vector<Cell> flatList;
				std::vector<Cell> outletList;

				locateOutlet(x, y, flatList, outletList);

				if (!outletList.empty()) {

					imposeGradient2LowerElevation(outletList, flatList, grad1);
					imposeGradientFromElevation(flatList, grad2);

					combineGradient(grad1, grad2, flatList);

					assignFlowInFlat(flatList, grad2);

					// assign outlet flows
					for (auto& c : outletList) {
						m_flowDem[idx(c.x,c.y)] = c.val;
					}

					// second pass
					assignFlowInFlat(flatList, grad1);
					iniGradient(grad1, grad2, flatList);

				}
			}
		}
	}

	// -------------------------------
	// STEP 4: cleanup edges
	// -------------------------------
	for (int y = 0; y < _ysize; y++) {
		for (int x = 0; x < _xsize; x++) {

			int idx = y * _xsize + x;

			if (onEdge(x, y) || m_flowDem[idx] > 8)
				m_flowDem[idx] = 0;
		}
	}

	// -------------------------------
	// STEP 5: cleanup memory
	// -------------------------------
	delete[] grad1;
	delete[] grad2;
}


bool FlowDirectionAlgorithm::hasFlow(unsigned char flowdirection) 
{
	return flowdirection >= 1 && flowdirection <= 8;
}

bool FlowDirectionAlgorithm::onEdge(int x, int y) 
{
	return x == 0 || x == _xsize - 1 ||
		y == 0 || y == _ysize - 1;
}

void FlowDirectionAlgorithm::findDirection(
	double listA[], double val,
	std::vector<FlowDirection>& listPos)
{
	for (int i = 0; i < 8; i++) {
		if (listA[i] == val)
			listPos.push_back(mapFlowLocation(i));
	}
}


double FlowDirectionAlgorithm::maxAdj(int x, int y, double* gradient, double listVal[])
{
	double center = gradient[idx(x,y)];
	double max = -1;
	int index = 0;
	int pos = 0;

	for (int dy = -1; dy <= 1; dy++) {
		for (int dx = -1; dx <= 1; dx++) {

			pos++;
			if (pos == 5) continue;

			int nx = x + dx;
			int ny = y + dy;

			double val;

			if (m_inDem[idx(nx,ny)] > m_inDem[idx(x,y)]) 
				val = rUNDEF;
			else 
				val = computeHeightDifference(center,gradient[idx(nx,ny)]);

			listVal[index++] = val;

			if (val > max)
				max = val;
		}
	}
	return max;
}

double FlowDirectionAlgorithm::maxAdj(int x, int y, double listVal[])
{
	double center = m_inDem[idx(x,y)];
	double max = -1;
	int index = 0;
	int pos = 0;

	for (int dy = -1; dy <= 1; dy++) {
		for (int dx = -1; dx <= 1; dx++) {

			pos++;
			if (pos == 5) continue;

			int nx = x + dx;
			int ny = y + dy;

			double nheight = m_inDem[idx(nx,ny)];
			if (nheight == rUNDEF)
				return rUNDEF;

			double val = (method == slope)
				? computeSlope(center, nheight, pos)
				: computeHeightDifference(center, nheight);

			listVal[index++] = val;

			if (val > max)
				max = val;
		}
	}
	return max;
}


//The slope is calculated by subscribing the neighbor's elevation value from the center
//distance 1.14 is concerned for diagonal cells
//Parameters
//h1 - the elevation of the center cell
//h2 - the elevation of the neighboring cell
//pos - the location number according to flow direction definition in ILWIS 	
//Returns - slope value of the center cell	
double FlowDirectionAlgorithm::computeSlope(double h1, double h2, int pos)
{
	return isEven(pos) ? (h1 - h2) : (h1 - h2) / 1.41;
}

bool FlowDirectionAlgorithm::isEven(int elem)
{
	return elem % 2 == 0;
}

double FlowDirectionAlgorithm::computeHeightDifference(double h1, double h2)
{
	return h1 - h2;
}

//Examine the flow direction location to perform one of the following:
//---1. if there are only two elements in the given list (maxium slop occurs in two neiboring cells):	
//---2. if there are only two elements in the given list (maxium slop occurs in two neiboring cells):
//make flow direction with S or W or E or N, if such a element exists in the given list
//otherwise, make the flow direction with the first element in the list.
//---3. if there are more than two cells with maximum slop, perform one of the following:
//if three cells located in one edge, make the flow direction value with the middle cell of the edge, otherwise,
//make flow direction with S or W or E or N, if such a element exists in the given list, otherwise
//make the flow direction with the first element in the list 	
//Parameters
//	listPos - the list filled with flow directions
//Returns - flow direction 	
FlowDirectionAlgorithm::FlowDirection FlowDirectionAlgorithm::getFlowDirection(
	const std::vector<FlowDirection>& listPos)
{
	if (listPos.size() == 1)
		return listPos[0];

	if (listPos.size() == 2) {
		for (auto fd : listPos) {
			if (fd == E || fd == S || fd == W || fd == N)
				return fd;
		}
	}

	if (isInOneEdge(listPos, SW, S, SE)) return S;
	if (isInOneEdge(listPos, NW, W, SW)) return W;
	if (isInOneEdge(listPos, NW, N, NE)) return N;
	if (isInOneEdge(listPos, NE, E, SE)) return E;

	for (auto fd : listPos) {
		if (fd == E || fd == S || fd == W || fd == N)
			return fd;
	}

	return listPos[0];
}


long MapFlowDirection::iLookUp(double rMax, int iCout, std::vector<int>& vPos)
{
	//rMax - max. value in slope / height diff.
	//iCout - number of positions with max. elements in vPos
	//vPos - positions with elements having max. value
	//iPos - location which the target cell flows to

	//This function examines the positons for elements with max value
	//then assigns flow direction accordinngly
	long iPos;
	if (rMax <= 0)
		iPos = 9; //a sink or flat cell
	else if (iCout < 3)
		iPos = vPos[0] + 1;
	else {
		//examine the positions of elements with max.
		//if three cells are adjacent, assign it to flow to the center
		//else flow to the first elem. in vector vPos for positions
		if (isInOneEdge(0, 1, 7, vPos))
			iPos = 1;
		else if (isInOneEdge(1, 2, 3, vPos))
			iPos = 3;
		else if (isInOneEdge(3, 4, 5, vPos))
			iPos = 5;
		else if (isInOneEdge(5, 6, 7, vPos))
			iPos = 7;
		else
			iPos = vPos[0] + 1;
	}
	return iPos;
}


int FlowDirectionAlgorithm::noflow = 0;

bool FlowDirectionAlgorithm::isInOneEdge(const std::vector<FlowDirection>& listPos, FlowDirection fd1, FlowDirection fd2, FlowDirection fd3) {

	bool fCondition1 = find(listPos.begin(), listPos.end(), fd1) != listPos.end();
	bool fCondition2 = find(listPos.begin(), listPos.end(), fd2) != listPos.end();
	bool fCondition3 = find(listPos.begin(), listPos.end(), fd3) != listPos.end();

	return fCondition1 && fCondition2 && fCondition3;
}


FlowDirectionAlgorithm::FlowDirection FlowDirectionAlgorithm::mapFlowLocation(int pos)
{
	switch (pos) {
	case 0: return NW;
	case 1: return N;
	case 2: return NE;
	case 3: return W;
	case 4: return E;
	case 5: return SW;
	case 6: return S;
	case 7: return SE;
	}
	return E;
}


void FlowDirectionAlgorithm::locateOutlet( int startX, int startY,std::vector<Cell>& flatList,std::vector<Cell>& outList)
{
	std::vector<Cell> srcList;
	std::vector<Cell> desList;

	double flatElevation = m_inDem[idx(startX,startY)];

	m_flowDem[idx(startX, startY)] = flag;

	srcList.push_back(Cell(startX, startY,10));
	flatList.push_back(Cell(startX, startY,10));

	while (!srcList.empty()) {

		desList.clear();

		for (auto& c : srcList) {

			for (int dy = -1; dy <= 1; dy++) {
				for (int dx = -1; dx <= 1; dx++) {

					int nx = c.x + dx;
					int ny = c.y + dy;

					if (nx < 0 || ny < 0 || nx >= _xsize || ny >= _ysize)
						continue;

					int nidx = idx(nx,ny);

					if (onEdge(nx, ny) || m_flowDem[nidx] == flag)
						continue;

					if (m_inDem[nidx] == flatElevation && hasFlow(m_flowDem[nidx])) {
						outList.push_back(Cell(nx, ny, m_flowDem[nidx]));
						desList.push_back(Cell(nx, ny,10));
						m_flowDem[nidx] = flag;
					}
					else if (m_inDem[nidx] == flatElevation) {
						m_flowDem[nidx] = flag;
						flatList.push_back(Cell(nx, ny,10));
						desList.push_back(Cell(nx, ny,10));
					}
				}
			}
		}
		srcList.swap(desList);
	}
}


void FlowDirectionAlgorithm::imposeGradient2LowerElevation( std::vector<Cell>& outletList,
	std::vector<Cell>& flatList,
	double* gradient)
{
	std::vector<Cell> srcList = outletList;
	std::vector<Cell> desList;

	while (!srcList.empty()) {

		// mark processed outlets
		for (auto& c : srcList) {
			m_flowDem[idx(c.x,c.y)] = flatcell;
		}

		// increment gradient
		for (auto& c : flatList) {
			if (m_flowDem[idx(c.x, c.y)] == flag)
				gradient[idx(c.x, c.y)] += increment;
		}

		// expand
		desList.clear();

		for (auto& c : srcList) {

			for (int dy = -1; dy <= 1; dy++) {
				for (int dx = -1; dx <= 1; dx++) {

					int nx = c.x + dx;
					int ny = c.y + dy;

					if (nx < 0 || ny < 0 || nx >= _xsize || ny >= _ysize)
						continue;

					if (m_flowDem[idx(nx,ny)] == flag) 
					{
						m_flowDem[idx(nx,ny)] = flatcell;
						desList.emplace_back(nx, ny);
					}
				}
			}
		}

		srcList.swap(desList);
	}
}


void FlowDirectionAlgorithm::imposeGradientFromElevation(
	std::vector<Cell>& flatList,
	double* gradient)
{
	bool changed;

	do {
		changed = false;
		std::vector<Cell> flagList;

		for (auto& c : flatList) {

			if (m_flowDem[idx(c.x,c.y)] == flag)
				continue;

			bool surrounded = true;

			for (int dy = -1; dy <= 1 && surrounded; dy++) {
				for (int dx = -1; dx <= 1; dx++) {

					int nx = c.x + dx;
					int ny = c.y + dy;

					if (nx < 0 || ny < 0 || nx >= _xsize || ny >= _ysize)
						continue;

					if (m_flowDem[idx(nx,ny)] != flatcell || m_flowDem[idx(nx, ny)] == flag) 
					{
						flagList.push_back(c);
						surrounded = false;
						break;
					}
				}
				if (!surrounded)
					break;
			}
		}

		if (!flagList.empty()) {
			changed = true;

			for (auto& c : flagList) {
				m_flowDem[idx(c.x, c.y)] = flag;
			}

			for (auto& c : flatList) {
				if (m_flowDem[idx(c.x, c.y)] == flag)
					gradient[idx(c.x, c.y)] += increment;
			}
		}

	} while (changed);
}


void FlowDirectionAlgorithm::combineGradient(
	double* grd1,
	double* grd2,
	std::vector<Cell>& flatList)
{
	for (auto& c : flatList) 
	{
		grd2[idx(c.x,c.y)] += grd1[idx(c.x, c.y)];
	}
}


void FlowDirectionAlgorithm::assignFlowInFlat(
	std::vector<Cell>& flatList,
	double* gradient)
{
	for (auto& c : flatList) {

		if (!hasFlow(m_flowDem[idx(c.x,c.y)])) 
		{

			double listVal[8];
			double max = maxAdj(c.x, c.y, gradient, listVal);

			if (max > 0) {
				std::vector<FlowDirection> listPos;
				findDirection(listVal, max, listPos);
				FlowDirection fd = getFlowDirection(listPos);
				m_flowDem[idx(c.x, c.y)] = Location[(byte)fd];
			}
		}
	}
}


void FlowDirectionAlgorithm::iniGradient(
	double* grd1,
	double* grd2,
	std::vector<Cell>& flatList)
{
	for (auto& c : flatList) {
		grd1[idx(c.x, c.y)] = 0;
		grd2[idx(c.x, c.y)] = 0;
	}
}


