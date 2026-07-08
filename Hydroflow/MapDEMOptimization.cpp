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
#include "table.h"
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

#include "domainitem.h"
#include "itemdomain.h"
#include "identifieritem.h"
#include "identifierrange.h"


#include "MapDEMOptimization.h"


using namespace Ilwis;
using namespace Hydroflow;


REGISTER_OPERATION(MapDEMOptimization)

MapDEMOptimization::MapDEMOptimization()
{
}


MapDEMOptimization::MapDEMOptimization(quint64 metaid, const Ilwis::OperationExpression& expr) : OperationImplementation(metaid, expr)
{

}

bool MapDEMOptimization::execute(ExecutionContext* ctx, SymbolTable& symTable)
{
	if (_prepState == sNOTPREPARED)
		if ((_prepState = prepare(ctx, symTable)) != sPREPARED)
			return false;

	bool resource = executeDemOptimization();

	if (resource && ctx != 0) {
		QVariant value;
		value.setValue<IRasterCoverage>(_outRaster);
		ctx->setOutput(symTable, value, _outRaster->name(), itRASTER, _outRaster->resource());
	}
	return resource;
}

Ilwis::OperationImplementation* MapDEMOptimization::create(quint64 metaid, const Ilwis::OperationExpression& expr)
{
	return new MapDEMOptimization(metaid, expr);
}


quint64 MapDEMOptimization::createMetadata()
{
	OperationResource operation({ "ilwis://operations/MapDEMOptimization" });
	operation.setSyntax("MapDEMOptimization(Input DEM Map,Drainage segment Map, Buffer distance (m),Smooth drop (m),Sharp drop (m))");

	operation.setDescription(TR("calculates 'corrected' height values of output DEM" ));
	operation.setInParameterCount({ 5 });
	operation.addInParameter(0, itRASTER, TR("Input DEM Map"), TR("input raster that is a un optimization DEM"));
	operation.addInParameter(1, itLINE, TR("Drainage segment Map"), TR("input drinage segment map "));
	operation.addInParameter(2, itFLOAT, TR("input value"), TR("Buffer distance indicates the 'area' at either side of the segments that will be affected by the gradual drop or raise"));
	operation.addInParameter(3, itFLOAT, TR("input value"), TR("Smooth drop  indicate the height with which segments should additionally be dropped or be raised in the output DEM"));
	operation.addInParameter(4, itFLOAT, TR("input value"), TR("Sharp drop indicate the height with which segments should additionally be dropped or be raised in the output DEM"));

	operation.parameterNeedsQuotes(1);
	operation.setOutParameterCount({ 1 });

	operation.addOutParameter(0, itRASTER, TR("Output Raster Map"), TR("output raster map that will contain DEM data"));
	operation.setKeywords("raster,drainage, segment");
	operation.checkAlternateDefinition();
	mastercatalog()->addItems({ operation });
	return operation.id();

}


bool MapDEMOptimization::executeDemOptimization()
{

	//PixelIterator iterPos = PixelIterator(_iterEmptyRaster, BoundingBox(), PixelIterator::fXYZ);
	//PixelIterator iterDEM = PixelIterator(_inRaster, BoundingBox(), PixelIterator::fXYZ);
	//PixelIterator iterOut = PixelIterator(_outRaster, BoundingBox(), PixelIterator::fXYZ);

	//PixelIterator inEnd = iterDEM.end();

	//PixelIterator iterTshd;
	//PixelIterator iterFlow;

	//std::fill(iterOut, iterOut.end(), 0);

	//if (m_tmMethods == tmMap)
	//{
	//	iterTshd = PixelIterator(_inThreshdRaster, BoundingBox(), PixelIterator::fXYZ);
	//	iterFlow = PixelIterator(_inFldRaster, BoundingBox(), PixelIterator::fXYZ);
	//	if (!_inThreshdRaster.isValid() || !_inFldRaster.isValid())
	//		return false;
	//}

	//while (iterPos != inEnd)
	//{
	//	Pixel pxl = iterPos.position();
	//	pxl.z = 0;
	//	if (*iterDEM(pxl) == iUNDEF || (m_tmMethods == tmMap && *iterTshd(pxl) == iUNDEF))
	//	{
	//		*iterOut(pxl) = 0;
	//	}
	//	else if (m_tmMethods == tmMap && *iterDEM(pxl) > *iterTshd(pxl))
	//	{
	//		*iterOut(pxl) = 2;
	//	}
	//	else if (m_tmMethods == tmValue && *iterDEM(pxl) > m_iThreshold)
	//	{
	//		*iterOut(pxl) = 2;
	//	}
	//	else
	//		*iterOut(pxl) = 1;
	//	iterPos++;
	//}


	////Construct dead-end drainage lines, if stream threshold map is used
	//if (m_tmMethods == tmMap)
	//{
	//	Pixel pospxl;
	//	PixelIterator iterPos = PixelIterator(_iterEmptyRaster, BoundingBox(), PixelIterator::fXYZ);
	//	while (iterPos != inEnd)
	//	{
	//		Pixel pxl = iterPos.position();
	//		if (!IsEdgeCell(pxl) && *iterOut(pxl) == 2 && *iterFlow(pxl) != 0)
	//		{

	//			pospxl = pxl;
	//			pospxl.z = 0;
	//			while (IsDeadendLine(pospxl))
	//			{
	//				*iterOut(pospxl) = 2;
	//			}
	//		}

	//		iterPos++;
	//	}
	//}

	//iterOut = PixelIterator(_outRaster, BoundingBox(), PixelIterator::fXYZ);
	//while (iterOut != inEnd)
	//{
	//	*iterOut = (*iterOut != 0) ? (*iterOut - 1) : rUNDEF;
	//	++iterOut;
	//}

	return true;
}


bool MapDEMOptimization::IsEdgeCell(Pixel pxl)
{
	if (pxl.y == 0 || pxl.y == _ysize - 1 ||
		pxl.x == 0 || pxl.x == _xsize - 1)
		return true;
	else
		return false;
}

bool MapDEMOptimization::IsDeadendLine(Pixel& pxl)
{
	bool fCondition = true;
	/*PixelIterator iterOut = PixelIterator(_outRaster, BoundingBox(), PixelIterator::fXYZ);
	PixelIterator iterFlow = PixelIterator(_inFldRaster, BoundingBox(), PixelIterator::fXYZ);

	int iFlow = *iterFlow(pxl);
	long in, jn;
	bool fCondition;
	Pixel pospxl;
	pospxl.z = 0;
	switch (iFlow)
	{
	case 1:
	{
		pospxl.y = pxl.y;
		pospxl.x = pxl.x + 1;
		fCondition = (*iterOut(pospxl) != 2);
		break;
	}
	case 2:
	{
		pospxl.y = pxl.y + 1;
		pospxl.x = pxl.x + 1;
		fCondition = (*iterOut(pospxl) != 2);
		break;
	}
	case 3:
	{
		pospxl.y = pxl.y + 1;
		pospxl.x = pxl.x;
		fCondition = (*iterOut(pospxl) != 2);
		break;
	}
	case 4:
	{
		pospxl.y = pxl.y + 1;
		pospxl.x = pxl.x - 1;
		fCondition = (*iterOut(pospxl) != 2);
		break;
	}
	case 5:
	{
		pospxl.y = pxl.y;
		pospxl.x = pxl.x - 1;
		fCondition = (*iterOut(pospxl) != 2);
		break;
	}
	case 6:
	{
		pospxl.y = pxl.y - 1;
		pospxl.x = pxl.x - 1;
		fCondition = (*iterOut(pospxl) != 2);
		break;
	}
	case 7:
	{
		pospxl.y = pxl.y - 1;
		pospxl.x = pxl.x;
		fCondition = (*iterOut(pospxl) != 2);
		break;
	}
	case 8:
	{
		pospxl.y = pxl.y - 1;
		pospxl.x = pxl.x + 1;
		fCondition = (*iterOut(pospxl) != 2);
		break;
	}
	default:
		pospxl.y = pxl.y;
		pospxl.x = pxl.x;
		fCondition = false;
	}

	pxl.y = pospxl.y;
	pxl.x = pospxl.x;

	if ((*iterOut(pospxl) == 0) || (IsEdgeCell(pospxl)))
		fCondition = false;*/

	return fCondition;
}


Ilwis::OperationImplementation::State MapDEMOptimization::prepare(ExecutionContext* ctx, const SymbolTable& st)
{
	OperationImplementation::prepare(ctx, st);
	QString demRasStr = _expression.parm(0).value();
	QString drnSegStr = _expression.parm(1).value();
	QString bufDisStr = _expression.parm(2).value();
	QString smthDropStr = _expression.parm(3).value();
	QString shpDropStr = _expression.parm(4).value();


	if (!_inDemRaster.prepare(demRasStr, itRASTER)) {
		ERROR2(ERR_COULD_NOT_LOAD_2, demRasStr, "");
		return sPREPAREFAILED;
	}

	QString drnStreamSegStr = drnSegStr;
	drnStreamSegStr.replace(".mpr", ".mps");
	if (!_inDranageSegmentMap.prepare(drnStreamSegStr, itFEATURE)) {
		ERROR2(ERR_COULD_NOT_LOAD_2, drnSegStr, "");
		return sPREPAREFAILED;
	}

	m_rBufdist = bufDisStr.toDouble();
	m_rSmooth = smthDropStr.toDouble();
	m_rSharp = shpDropStr.toDouble();
	
	QString outputName = _expression.parm(0, false).value();

	int copylist = itRASTERSIZE | itENVELOPE | itCOORDSYSTEM | itGEOREF | itNUMERICDOMAIN;
	_outRaster = OperationHelperRaster::initialize(_inDemRaster.as<IlwisObject>(), itRASTER, copylist);
	if (_outRaster.isValid())
	{
		INamedIdDomain outDomain;
		QString catchName = "CatchID";
		outDomain.prepare();
		NamedIdentifierRange* idrange = new NamedIdentifierRange();
		outDomain->range(idrange);

		DataDefinition def(outDomain);
		_outRaster->datadefRef() = def;
		for (int band = 0; band < _outRaster->size().zsize(); ++band) {
			_outRaster->datadefRef(band) = def;
		}

	}
	else
	{
		ERROR1(ERR_NO_INITIALIZED_1, "output rastercoverage");
		return sPREPAREFAILED;
	}

	_xsize = _inDemRaster->size().xsize();
	_ysize = _inDemRaster->size().ysize();

	return sPREPARED;
}



