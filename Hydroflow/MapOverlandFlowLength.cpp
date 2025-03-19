
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
#include "coordinatedomain.h"
#include "ellipsoid.h"
#include "conventionalcoordinatesystem.h"
#include "table.h"
#include "MapOverlandFlowLength.h"

using namespace Ilwis;
using namespace Hydroflow;
const double rDefaultEarthRadius = 6371007.0;

REGISTER_OPERATION(MapOverlandFlowLength)

#define	sDownstreamCoord  "DownstreamCoord"
#define	sUpstreamCoord  "UpstreamCoord"
#define sUpstreamID "UpstreamLinkID"

MapOverlandFlowLength::MapOverlandFlowLength()
{
}


MapOverlandFlowLength::MapOverlandFlowLength(quint64 metaid, const Ilwis::OperationExpression& expr) : OperationImplementation(metaid, expr)
{

}

bool MapOverlandFlowLength::execute(ExecutionContext* ctx, SymbolTable& symTable)
{
    if (_prepState == sNOTPREPARED)
        if ((_prepState = prepare(ctx, symTable)) != sPREPARED)
            return false;

    executeLandFlowLength();

    bool resource = true;
    if (resource && ctx != 0) {
        QVariant value;
        value.setValue<IRasterCoverage>(_outRaster);
        ctx->setOutput(symTable, value, _outRaster->name(), itRASTER, _outRaster->resource());
    }
    return resource;
}

Ilwis::OperationImplementation* MapOverlandFlowLength::create(quint64 metaid, const Ilwis::OperationExpression& expr)
{
    return new MapOverlandFlowLength(metaid, expr);
}

Ilwis::OperationImplementation::State MapOverlandFlowLength::prepare(ExecutionContext* ctx, const SymbolTable& st)
{
	OperationImplementation::prepare(ctx, st);
	QString drnstreamStr = _expression.parm(0).value();
	QString flowrasterStr = _expression.parm(1).value();


	if (!_inDrainRaster.prepare(drnstreamStr, itRASTER)) {
		ERROR2(ERR_COULD_NOT_LOAD_2, drnstreamStr, "");
		return sPREPAREFAILED;
	}

	if (!_inFlowRaster.prepare(flowrasterStr, itRASTER)) {
		ERROR2(ERR_COULD_NOT_LOAD_2, flowrasterStr, "");
		return sPREPAREFAILED;
	}

	// Check if we are dealing with FlowDirection.dom; if yes, then recalc raw vals of flow direction map. We also allow Value maps with values 0 til 8 (0 = flat / undef, 1 til 8 are the directions).
	IDomain itemdom = _inFlowRaster->datadefRef().domain();
	if (itemdom.isValid() && hasType(itemdom->valueType(), itTHEMATICITEM | itNUMERICITEM | itTIMEITEM | itNAMEDITEM))
	{
		_inFlowRaster.set(_inFlowRaster->clone());
		m_iterFlow = PixelIterator(_inFlowRaster, BoundingBox(), PixelIterator::fXYZ);
		PixelIterator iterFldEnd = m_iterFlow.end();
		while (m_iterFlow != iterFldEnd)
		{
			*m_iterFlow = (*m_iterFlow != rUNDEF) ? (*m_iterFlow + 1) : 0; // shift the values, to make them the same as ilwis3 (0 is undef, 1..8 are the direction)
			++m_iterFlow;
		}
	}

	IDomain itemdomDno = _inDrainRaster->datadefRef().domain();
	if (itemdomDno.isValid() && hasType(itemdomDno->valueType(), itTHEMATICITEM | itNUMERICITEM | itTIMEITEM | itNAMEDITEM))
	{
		_inDrainRaster.set(_inDrainRaster->clone());
		m_iterDrain = PixelIterator(_inDrainRaster, BoundingBox(), PixelIterator::fXYZ);
		PixelIterator iterDnoEnd = m_iterDrain.end();
		while (m_iterDrain != iterDnoEnd)
		{
			*m_iterDrain = (*m_iterDrain != rUNDEF) ? (*m_iterDrain + 1) : 0; // shift the values, to make them the same as ilwis3 (0 is undef, 1..8 are the direction)
			++m_iterDrain;
		}
	}
	
	int copylist = itRASTERSIZE | itENVELOPE | itCOORDSYSTEM | itGEOREF | itNUMERICDOMAIN;
	_outRaster = OperationHelperRaster::initialize(_inDrainRaster.as<IlwisObject>(), itRASTER, copylist);
	if (!_outRaster.isValid()) {
		ERROR1(ERR_NO_INITIALIZED_1, "output rastercoverage");
		return sPREPAREFAILED;
	}
	else
	{
		IDomain dom("code=domain:value");
		_outRaster->datadefRef() = DataDefinition(dom);

		for (quint32 i = 0; i < _outRaster->size().zsize(); ++i) {
			QString index = _outRaster->stackDefinition().index(i);
			_outRaster->setBandDefinition(index, DataDefinition(dom));
		}

	}

	_xsize = _inDrainRaster->size().xsize();
	_ysize = _inDrainRaster->size().ysize();

    return sPREPARED;
}

quint64 MapOverlandFlowLength::createMetadata()
{
	OperationResource operation({ "ilwis://operations/MapOverlandFlowLength" });
    operation.setSyntax("MapOverlandFlowLength(DrainageNetworkOrderMap,FlowDiractionMap)");

    operation.setDescription(TR("calculates for each pixel the overland distance towards the 'nearest' drainage"));
    operation.setInParameterCount({ 2 });
    operation.addInParameter(0, itRASTER, TR("Drainage NetWork Ordering Map"), TR("input raster that is the output of the Drainage Network Ordering operation"));
    operation.addInParameter(1, itRASTER, TR("Flow Direction Map"), TR("input raster that is the output of the Flow direction operation"));

    operation.parameterNeedsQuotes(1);

    operation.setOutParameterCount({ 1 });

    operation.addOutParameter(0, itRASTER, TR("Output Raster Map"), TR("output raster map that will contain flow length"));
    operation.setKeywords("raster,table");
	operation.checkAlternateDefinition();
	mastercatalog()->addItems({ operation });
    return operation.id();


}

void InitFlowNums(std::vector<byte>& vReceiveNum)
{
	//	Flow number				Receive number			Loop ordering
	//	-------								-------					-------	  	
  //	|6|7|8|								|2|3|4|					|0|1|2|
	//	-------								-------					-------
	//	|5| |1|								|1| |5|					|3|4|5|
	//	-------								-------					-------
	//	|4|3|2|								|8|7|6|					|6|7|8|
	//	-------								-------					-------
	//
	vReceiveNum.resize(9);

	vReceiveNum[0] = 2;
	vReceiveNum[1] = 3;
	vReceiveNum[2] = 4;
	vReceiveNum[3] = 1;
	vReceiveNum[4] = 0;
	vReceiveNum[5] = 5;
	vReceiveNum[6] = 8;
	vReceiveNum[7] = 7;
	vReceiveNum[8] = 6;
}

bool MapOverlandFlowLength::executeLandFlowLength()
{

	//m_iterDrain = PixelIterator(_inDrainRaster, BoundingBox(), PixelIterator::fXYZ);
	//m_iterFlow = PixelIterator(_inFlowRaster, BoundingBox(), PixelIterator::fXYZ);
	m_iterOut = PixelIterator(_outRaster, BoundingBox(), PixelIterator::fXYZ);

	
	PixelIterator inEnd = m_iterDrain.end();

	// init output
	while (m_iterOut != inEnd)
	{
		*m_iterOut = rUNDEF;
		++m_iterOut;
	}

	ITable tblDrnAtt = _inDrainRaster->attributeTable();

	std::vector<QVariant> colDownstreamCoord;
	colDownstreamCoord = tblDrnAtt->column(sDownstreamCoord);

	std::vector<QVariant> colUpstreamCoord;
	colUpstreamCoord = tblDrnAtt->column(sUpstreamCoord);

	std::vector<QVariant> colUpstreamID;
	colUpstreamID = tblDrnAtt->column(sUpstreamID);

	std::vector<QVariant> colStreamID = tblDrnAtt->column(_inDrainRaster->primaryKey());
	long iSize = colStreamID.size();

	InitFlowNums(m_vReceiveNum);

	for (long i = 0; i < iSize; i++) {
		long iDrainageID = colStreamID[i].toInt()+1;
		if (iDrainageID == iUNDEF)
			continue;

		Pixel downCoord = CoordinateStringToPixel(colDownstreamCoord[i].toString());
		downCoord.z = 0;
		m_rcUpstream = CoordinateStringToPixel(colUpstreamCoord[i].toString());

		//---here also upstream links needed to be able to process drainage links seperately
		std::vector<long> vUpstreamLinks;
		SplitString(colUpstreamID[i].toString(), QString(","), vUpstreamLinks);

		bool fUpstreamlins;
		if (vUpstreamLinks.size() > 0)
			fUpstreamlins = true;
		else
			fUpstreamlins = false;

		Lengths2Stream(iDrainageID, downCoord, fUpstreamlins);
	
	}
	
	return true;
}

void MapOverlandFlowLength::SplitString(QString s, QString mid, std::vector<long>& results)
{
	results.clear();

	s.replace("{", "");
	s.replace("}", "");
	QStringList strlst = s.split(mid);

	for (unsigned int i = 0; i < strlst.size(); i++)
	{
		long res = strlst[i].toLong();
		if (res != iUNDEF && res>0 )
			results.push_back(res);
	}
}



Pixel MapOverlandFlowLength::CoordinateStringToPixel(QString coordStr)
{
	coordStr.replace("{", "");
	coordStr.replace("}", "");

	QStringList coodrs = coordStr.split(" ");

	QString xstr = coodrs[0];
	QString ystr = coodrs[1];

	Coordinate crd;
	crd.x = xstr.toDouble();
	crd.y = ystr.toDouble();
	crd.z = 0;
	return _inDrainRaster->georeference()->coord2Pixel(crd);
}


void MapOverlandFlowLength::Lengths2Stream(long iStreamID, Pixel rc, bool fUpstreamlins)
{
	//To match the index from 0 in the input matrix
	PixelIterator iterOut = PixelIterator(_outRaster, BoundingBox(), PixelIterator::fXYZ);

	*iterOut(rc) = 0;  //It should be 0 to outlet cell of the segment
	std::vector<Pixel> vStartCells; 
	vStartCells.push_back(rc);

	std::vector<Pixel>::iterator pos;
	std::vector<Pixel> vFlow2Cells; //*stores cells that drain to start cells 
	m_vCellsOnDivide.resize(0);
	do
	{
		vFlow2Cells.resize(0);
		for (pos = vStartCells.begin(); pos < vStartCells.end(); ++pos)
		{
			Pixel rcStart(*pos);
			int index = 0;
			bool isDivide = true;
			for (int i = -1; i < 2; ++i)
			{
				for (int j = -1; j < 2; ++j)
				{
					Pixel rcCur;
					rcCur.z = 0;
					rcCur.y = rcStart.y + i;
					rcCur.x = rcStart.x + j;
					if ( !IsEdgeCell(rcCur) )
					{
						bool isFlowNum =( *m_iterFlow(rcCur) == m_vReceiveNum[index]);
						bool isUpStream = ((rcCur == m_rcUpstream) && (fUpstreamlins == true));
						bool isUndef = (*iterOut(rcCur) == rUNDEF);
						bool isOverlandPixel = ((*m_iterDrain(rcCur) < 1) || (*m_iterDrain(rcCur) == iStreamID));
						bool isFlow = ((isFlowNum) && (isUpStream != true) && (isUndef));
						if (isFlow && isOverlandPixel)
						{
							if ((*m_iterDrain(rcCur) == iStreamID) )
								*iterOut(rcCur) = 0; //set dist value to 0, it should be
							else
							{
								Coordinate c1 = _inDrainRaster->georeference()->pixel2Coord(rcStart);
								Coordinate c2 = _inDrainRaster->georeference()->pixel2Coord(rcCur);
								double  rDist = rDistance(c1, c2);
								*iterOut(rcCur) = *iterOut(rcStart) + rDist;
							}
							vFlow2Cells.push_back(rcCur);
							isDivide = false;
						}
					}
					index++;
				}
			}
			if (isDivide)
				m_vCellsOnDivide.push_back(rcStart);
		}
		vStartCells.swap(vFlow2Cells);

	} while (vFlow2Cells.size() != 0); //no more element can be located 
}


bool MapOverlandFlowLength::IsEdgeCell(Pixel pxl)
{
	if (pxl.y == 0 || pxl.y == _ysize - 1 ||
		pxl.x == 0 || pxl.x == _xsize - 1)
		return true;
	else
		return false;
}


void MapOverlandFlowLength::Lengths2Divide(long iStreamID, Pixel rcDownstreamCoord)
{
	//std::vector<Pixel>::iterator pos;
	//for (pos = m_vCellsOnDivide.begin(); pos < m_vCellsOnDivide.end(); ++pos)
	//{
	//	Pixel rcStart = (*pos);
	//	m_vOutput_d[rcStart.Row][rcStart.Col] = 0;
	//	bool fCheck = true;
	//	do
	//	{
	//		int iFlowNum = m_vFlowDir[rcStart.Row][rcStart.Col];
	//		RowCol rcEnd;
	//		switch (iFlowNum)
	//		{
	//		case 1: {	//location 1/E
	//			rcEnd.Row = rcStart.Row;
	//			rcEnd.Col = rcStart.Col + 1;
	//		}
	//			  break;
	//		case 2: { //location 2/SE 
	//			rcEnd.Row = rcStart.Row + 1;
	//			rcEnd.Col = rcStart.Col + 1;
	//		}
	//			  break;
	//		case 3: {	//location 3/S
	//			rcEnd.Row = rcStart.Row + 1;
	//			rcEnd.Col = rcStart.Col;
	//		}
	//			  break;
	//		case 4: { //Location 4/SW
	//			rcEnd.Row = rcStart.Row + 1;
	//			rcEnd.Col = rcStart.Col - 1;
	//		}
	//			  break;
	//		case 5: {	//location 5/W
	//			rcEnd.Row = rcStart.Row;
	//			rcEnd.Col = rcStart.Col - 1;
	//		}
	//			  break;
	//		case 6: {	//location 6/NW 
	//			rcEnd.Row = rcStart.Row - 1;
	//			rcEnd.Col = rcStart.Col - 1;
	//		}
	//			  break;
	//		case 7: {	//location 7/N
	//			rcEnd.Row = rcStart.Row - 1;
	//			rcEnd.Col = rcStart.Col;
	//		}
	//			  break;
	//		case 8: {	//location 8/NE
	//			rcEnd.Row = rcStart.Row - 1;
	//			rcEnd.Col = rcStart.Col + 1;
	//		}
	//			  break;
	//		default: {
	//			fCheck = false;
	//		}
	//			   break;
	//		}
	//		Coord c1 = mp->gr()->cConv(rcStart);
	//		Coord c2 = mp->gr()->cConv(rcEnd);
	//		double  rDist = rDistance(c1, c2) + m_vOutput_d[rcStart.Row][rcStart.Col];

	//		if ((m_vOutput_d[rcEnd.Row][rcEnd.Col] == rUNDEF) || (m_vOutput_d[rcEnd.Row][rcEnd.Col] > rDist))
	//		{
	//			m_vOutput_d[rcEnd.Row][rcEnd.Col] = rDist;
	//		}
	//		if (IsEdgeCell(rcEnd.Row, rcEnd.Col))
	//			fCheck = false;
	//		else if ((rcDownstreamCoord.Row - 1 == rcEnd.Row) && (rcDownstreamCoord.Col - 1 == rcEnd.Col))
	//			fCheck = false;

	//		rcStart = rcEnd;
	//	} while (fCheck);
	//}
}

bool MapOverlandFlowLength::fLatLonCoords()
{
	return _inDrainRaster->coordinateSystem()->isLatLon();
}


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


double MapOverlandFlowLength::rDistance(Coordinate cd1, Coordinate cd2)
{
	double rDist;
	if (fLatLonCoords())
	{
		double rRadi = rDefaultEarthRadius;
		IConventionalCoordinateSystem projectedCsy = _inDrainRaster->coordinateSystem().as<ConventionalCoordinateSystem>();
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


