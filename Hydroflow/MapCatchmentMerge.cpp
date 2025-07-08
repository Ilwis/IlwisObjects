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
#include "ilwiscoordinate.h"
#include "blockiterator.h"
#include "pixeliterator.h"
#include "conventionalcoordinatesystem.h"
#include "ellipsoid.h"
#include "table.h"
#include "symboltable.h"
#include "columndefinition.h"
#include "basetable.h"
#include "flattable.h"
#include "domainitem.h"
#include "itemdomain.h"
#include "identifieritem.h"
#include "identifierrange.h"

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

#include "coordinatedomain.h"
#include "ellipsoid.h"
#include "itemdomain.h"
#include "thematicitem.h"

#include "georefimplementation.h"
#include "simpelgeoreference.h"
#include "cornersgeoreference.h"

#include "MapCatchmentMerge.h"
#include <QRegExp>


using namespace Ilwis;
using namespace Hydroflow;

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

MapCatchmentMerge::MapCatchmentMerge()
{
}


MapCatchmentMerge::MapCatchmentMerge(quint64 metaid, const Ilwis::OperationExpression& expr) : OperationImplementation(metaid, expr)
{

}


bool MapCatchmentMerge::execute(ExecutionContext* ctx, SymbolTable& symTable)
{
	if (_prepState == sNOTPREPARED)
		if ((_prepState = prepare(ctx, symTable)) != sPREPARED)
			return false;

	bool resource = executeCatchmentMerging();
	if (resource && ctx != 0)
	{
		//CreatePolygonTableElements();
		/////////////////////////////////////////////
		if (_outMergeRaster.isValid())
		{
			_outMergeRaster->setAttributes(_outputTable);
			std::string connectivity = "8";
			std::string output_name = "polygon_object_" + QString::number(Ilwis::Identity::newAnonymousId()).toStdString();
			QString expr = QString::fromStdString(output_name + "=raster2polygon(" + _outMergeRaster->name().toStdString() + "," + connectivity + ",true)");
			Ilwis::commandhandler()->execute(expr, ctx, symTable);
			Ilwis::Symbol result = symTable.getSymbol(ctx->_results[0]);
			
			if (result._type == itFEATURE && result._var.canConvert<Ilwis::IFeatureCoverage>())
			{

				_outputPolygonMap = Ilwis::IFeatureCoverage(result._var.value<Ilwis::IFeatureCoverage>());
				CreatePolygonTableElements();
				//_outputPolygonMap->setAttributes(_outputTable);

				quint32 count = _outputPolygonMap->featureCount();
				quint32 recordcount = _outputTable->recordCount();
				std::vector<QVariant> colUpstreamLinkCatchment = _outputTable->column("UpStreamLinkCatchment");
				std::vector<QVariant> colTotalUpstreamArea = _outputTable->column("TotalUpstreamArea");

					for (int rec = 0; rec < count; ++rec)
					{
						const SPFeatureI& feature = _outputPolygonMap->feature(rec);
						quint64 id = feature->featureid();
						geos::geom::Geometry* polygon = dynamic_cast<geos::geom::Geometry*>(feature->geometry().get());
						if (polygon)
						{
							double area = polygon->getArea();
							double length = polygon->getLength();
							geos::geom::Point* centroid = polygon->getInteriorPoint();
							Coordinate crd;
							crd.x = centroid->getX();
							crd.y = centroid->getY();
							crd.z = 0;
							
							QVariant c1a;
							c1a.setValue(crd);
							_outputTable->setCell("CenterCatchment", rec, c1a);
							_outputTable->setCell("Perimeter", rec, length);
							_outputTable->setCell("CatchmentArea", rec, area);
							_outputTable->setCell("TotalUpstreamArea", rec, area);

						/*	double rarea = area;
							std::vector<long> vLinks;
							vLinks.clear();
							if (colUpstreamLinkCatchment.size()>0)
								SplitString(colUpstreamLinkCatchment[rec].toString(), QString(","),vLinks);
							if ((vLinks.size() == 1) && (vLinks[0] == 0))
								_outputTable->setCell("TotalUpstreamArea", rec, area);*/

						}

					}
			}
		}
		//////////////////////////////////////////////////////////////////

		if (_outMergeRaster.isValid()) {
			_outMergeRaster->setAttributes(_outputTable);
			QVariant outraster;
			outraster.setValue<IRasterCoverage>(_outMergeRaster);
			logOperation(_outMergeRaster, _expression, { _inDrngOrderRaster, _inFldRaster });
			ctx->addOutput(symTable, outraster, _outMergeRaster->name(), itRASTER, _outMergeRaster->resource());
		}


		if ( _longestPathFeature.isValid())
		{
			_longestPathFeature->setAttributes(_outputPathTable);
			QVariant value;
			value.setValue<IFeatureCoverage>(_longestPathFeature);
			logOperation(_longestPathFeature, _expression, { _inDrngOrderRaster });
			ctx->addOutput(symTable, value, "extractedlongpath", itFEATURE, _longestPathFeature->resource());
		}


		if (_outputExtraxtedSegmentMap.isValid())
		{
			_outputExtraxtedSegmentMap->setAttributes(_outputExtractSegTable);
			QVariant value;
			value.setValue<IFeatureCoverage>(_outputExtraxtedSegmentMap);
			logOperation(_outputExtraxtedSegmentMap, _expression, { _inDrngOrderRaster });
			ctx->addOutput(symTable, value, "extractedmatchsegments", itFEATURE, _outputExtraxtedSegmentMap->resource());
		}


		if (_outputPolygonMap.isValid())
		{
			_outputPolygonMap->setAttributes(_outputTable);
			QVariant value;
			value.setValue<IFeatureCoverage>(_outputPolygonMap);
			logOperation(_outputPolygonMap, _expression, { _inDrngOrderRaster });
			ctx->addOutput(symTable, value, _outputPolygonMap->name(), itFEATURE, _outputPolygonMap->resource());
		}

	}

	return resource;
}


Ilwis::OperationImplementation* MapCatchmentMerge::create(quint64 metaid, const Ilwis::OperationExpression& expr)
{
	return new MapCatchmentMerge(metaid, expr);
}

Ilwis::OperationImplementation::State MapCatchmentMerge::prepare(ExecutionContext* ctx, const SymbolTable& st)
{
	m_UseOutlets = false;
	m_useExtraOrder = false;
	m_includeUndefine = false;

	OperationImplementation::prepare(ctx, st);
	QString drnstreamStr = _expression.parm(0).value();
	QString flowrasterStr = _expression.parm(1).value();
	QString accrasterStr = _expression.parm(2).value();
	QString demrasterStr = _expression.parm(3).value();


	QString drnStreamSegStr = drnstreamStr;
	drnStreamSegStr.replace(".mpr", ".mps");
	if (!_internelDranageSegmentMap.prepare(drnStreamSegStr, itFEATURE)) {
		ERROR2(ERR_COULD_NOT_LOAD_2, drnStreamSegStr, "");
		return sPREPAREFAILED;
	}

	QString outputName = _expression.parm(0, false).value();

	if (!_inDrngOrderRaster.prepare(drnstreamStr, itRASTER)) {
		ERROR2(ERR_COULD_NOT_LOAD_2, drnstreamStr, "");
		return sPREPAREFAILED;
	}

	if (!_inFldRaster.prepare(flowrasterStr, itRASTER)) {
		ERROR2(ERR_COULD_NOT_LOAD_2, flowrasterStr, "");
		return sPREPAREFAILED;
	}

	if (!_inAccRaster.prepare(accrasterStr, itRASTER)) {
		ERROR2(ERR_COULD_NOT_LOAD_2, accrasterStr, "");
		return sPREPAREFAILED;
	}

	if (!_inDemRaster.prepare(demrasterStr, itRASTER)) {
		ERROR2(ERR_COULD_NOT_LOAD_2, demrasterStr, "");
		return sPREPAREFAILED;
	}

	// Check if we are dealing with FlowDirection.dom; if yes, then recalc raw vals of flow direction map. We also allow Value maps with values 0 til 8 (0 = flat / undef, 1 til 8 are the directions).
	IDomain itemdom = _inFldRaster->datadefRef().domain();
	if (itemdom.isValid() && hasType(itemdom->valueType(), itTHEMATICITEM | itNUMERICITEM | itTIMEITEM | itNAMEDITEM))
	{
		_inFldRaster.set(_inFldRaster->clone());
		PixelIterator iterFld = PixelIterator(_inFldRaster, BoundingBox(), PixelIterator::fXYZ);
		PixelIterator iterFldEnd = iterFld.end();
		while (iterFld != iterFldEnd)
		{
			*iterFld = (*iterFld != rUNDEF) ? (*iterFld + 1) : 0; // shift the values, to make them the same as ilwis3 (0 is undef, 1..8 are the direction)
			++iterFld;
		}
	}

	IDomain itemdomDno = _inDrngOrderRaster->datadefRef().domain();
	if (itemdomDno.isValid() && hasType(itemdomDno->valueType(), itTHEMATICITEM | itNUMERICITEM | itTIMEITEM | itNAMEDITEM))
	{
		_inDrngOrderRaster.set(_inDrngOrderRaster->clone());
		PixelIterator iterDno = PixelIterator(_inDrngOrderRaster, BoundingBox(), PixelIterator::fXYZ);
		PixelIterator iterDnoEnd = iterDno.end();
		while (iterDno != iterDnoEnd)
		{
			*iterDno = (*iterDno != rUNDEF) ? (*iterDno + 1) : 0; // shift the values, to make them the same as ilwis3 (0 is undef, 1..8 are the direction)
			++iterDno;
		}
	}

	int copylist = itRASTERSIZE | itENVELOPE | itCOORDSYSTEM | itGEOREF | itNUMERICDOMAIN;
	_outMergeRaster = OperationHelperRaster::initialize(_inDrngOrderRaster.as<IlwisObject>(), itRASTER, copylist);
	if (_outMergeRaster.isValid())
	{
		QString catchName = "CatchID";
		_outDomain.prepare();
		_idrange = new NamedIdentifierRange();
		_outDomain->range(_idrange);

		DataDefinition def(_outDomain);
		_outMergeRaster->datadefRef() = def;
		for (int band = 0; band < _outMergeRaster->size().zsize(); ++band) {
			_outMergeRaster->datadefRef(band) = def;
		}
	
	}
	else
	{
		ERROR1(ERR_NO_INITIALIZED_1, "output rastercoverage");
		return sPREPAREFAILED;
	}

	_xsize = _inDrngOrderRaster->size().xsize();
	_ysize = _inDrngOrderRaster->size().ysize();



	// initialize the georeference, with name, bounding box, coordinatesystem and envelope
	IGeoReference grf(outputName);
	grf->create("corners");
	grf->size(Size<>(_xsize, _ysize, 1)); // sets the bounding box
	grf->coordinateSystem(_inDrngOrderRaster->coordinateSystem());
	QSharedPointer< CornersGeoReference> spGrf = grf->as< CornersGeoReference>();
	spGrf->internalEnvelope(_inDrngOrderRaster->envelope());
	grf->compute(); // all members are set, now the initialization can take place
	_inputgrf = grf;
	_csy = _inputgrf->coordinateSystem();


	//////////////////////////////////////////////////////////////////////////////
	// prepare path feature

	_longestPathFeature.prepare(QString(INTERNAL_CATALOG + "/%1").arg(outputName));
	_longestPathFeature->coordinateSystem(_inDrngOrderRaster->georeference()->coordinateSystem());
	_longestPathFeature->envelope(_inDrngOrderRaster->georeference()->envelope());

	/////////////////////////////////////////////////////////////////
	if (_outMergeRaster.isValid())
	{
		QString catchName = "CatchID";
		_outDomain.prepare();
		_idrange = new NamedIdentifierRange();
		_outDomain->range(_idrange);

		DataDefinition def(_outDomain);
		_outMergeRaster->datadefRef() = def;
		for (int band = 0; band < _outMergeRaster->size().zsize(); ++band) {
			_outMergeRaster->datadefRef(band) = def;
		}

		return sPREPARED;
	}

	return sPREPARED;
}


bool MapCatchmentMerge::isDigitStr(QString str)
{
	QByteArray barr = str.toLatin1();//QString to char
	const char* s = barr.data();

	while (*s && *s >= '0' && *s <= '9')
		s++;

	if (*s)
		return false;

	return true;

}

quint64 MapCatchmentMerge::createMetadata(Ilwis::OperationResource& operation)
{
	operation.checkAlternateDefinition();
	mastercatalog()->addItems({ operation });
	return operation.id();

}


class SortOutletsLessClass //compare two elements for sort algorithm
{
public:
	SortOutletsLessClass(std::vector<OutletLocation>& vol) :
		m_vOutlet(vol)
	{
	}
	bool operator()(OutletLocation ol1, OutletLocation ol2)
	{
		return ol1.StreveOrder < ol2.StreveOrder;
	}
private:
	std::vector<OutletLocation>& m_vOutlet;
};


bool MapCatchmentMerge::executeCatchmentMerging()
{
	//---Create attibute table associated with output map
	CreateTable();

	AddLink2StreamSegments();

	long id = 0;
	if (m_UseOutlets)
	{
		InitOutletVector();
		//---Sort outlets by using operator(elem1, elem2) sorting algorithms by Streve number
		if (m_vOutlet.size() > 0)
		{
			SortOutletsLessClass streve(m_vOutlet);
			sort(m_vOutlet.begin(), m_vOutlet.end(), streve);
		}
	}
	else //---Use stream orders option
	{
		//---Evaluate joint down-flow coords/outlet locations by stream order
		EvaluateJointOutletsbyOrder();
		InitJointOutletVector();
	}

	m_iterInDrng = PixelIterator(_inDrngOrderRaster, BoundingBox(), PixelIterator::fXYZ);
	m_iterFld = PixelIterator(_inFldRaster, BoundingBox(), PixelIterator::fXYZ);
	m_iterOut = PixelIterator(_outMergeRaster, BoundingBox(), PixelIterator::fXYZ);

	PixelIterator iterDEM = PixelIterator(_inDemRaster, BoundingBox(), PixelIterator::fXYZ);
	PixelIterator inEnd = iterDEM.end();

	// init output
	while (m_iterOut != inEnd)
	{
		*m_iterOut = iUNDEF;
		++m_iterOut;
	}



	//Merging sub-catchments here
	m_totalcount = m_vOutlet.size();
	for (std::vector<OutletLocation>::iterator pos = m_vOutlet.begin(); pos < m_vOutlet.end(); ++pos)
	{
		Pixel pxl = pos->pxl;
		OutletLocation ol = (*pos);
		id++;
		m_vUpCatchmentIDs.resize(0);
		m_vStreamsInCatchment.resize(0);
		m_sUpLinkCatchment = QString("");
		m_iOutletVal = *m_iterInDrng(pxl);
		//AddDomainItem(m_dm, id);

		//---put stream id, defined by the joint outlet location in the catchment, if it is a no-flow stream
		if ((m_UseOutlets && ol.isOnNode != true) || (!m_UseOutlets))
			m_vStreamsInCatchment.push_back(*m_iterInDrng(pxl));

		std::vector<long> vUpstreamID;
		MergeCatchment(pxl, id, false);
		UpdateUpLinkCatchment(id);
		UpdateDownLinkCatchment(id);
		UpdateLink2StreamSegments(id, ol);

	}
	if (!m_UseOutlets && m_useExtraOrder) //&& m_iStreamOrders < m_iMaxOrderNumber)
	{
		ExtractOriginalOrder(id);
	}


	//Clean up
	m_vUpCatchmentIDs.resize(0);

	ExtractSegments();
	CreateTableSegmentsExtracted();

	ComputeOtherAttributes();

	m_iterOut = PixelIterator(_outMergeRaster, BoundingBox(), PixelIterator::fXYZ);

	// init output
	while (m_iterOut != inEnd)
	{
		if ((*m_iterOut > 255) || (*m_iterOut == 0))
			*m_iterOut = iUNDEF;
		else
			*m_iterOut = *m_iterOut - 1;
		++m_iterOut;
	}

	m_vOutlet.resize(0);
	m_vStreamCoord.resize(0);
	m_vDrainageAtt.resize(0);
	m_vStreamsInCatchment.resize(0);

	return true;
}


static double rComputeSlope(double rDrop, double rLength, bool fDegree)
{
	if (rLength > 0)
	{
		if (fDegree)
			return (180 / M_PI) * atanl(abs(rDrop / rLength)); //degree
		else
			return (rDrop / rLength) * 100; //percent
	}
	else
		return 0;
}


void MapCatchmentMerge::CreatePolygonTableElements()
{
	ICoordinateDomain crddom;
	crddom.prepare();
	crddom->setCoordinateSystem(_csy);
	_outputTable->addColumn("CenterCatchment", crddom, true);

	_outputTable->addColumn("Perimeter", IlwisObject::create<IDomain>("value"), true);
	ColumnDefinition& coldef0 = _outputTable->columndefinitionRef("Perimeter");
	coldef0.datadef().range(new NumericRange(0, 1.0e300, 0.01));

	_outputTable->addColumn("CatchmentArea", IlwisObject::create<IDomain>("value"), true);
	ColumnDefinition& coldef1 = _outputTable->columndefinitionRef("CatchmentArea");
	coldef1.datadef().range(new NumericRange(0, 1.0e300, 0.01));

	_outputTable->addColumn("TotalUpstreamArea", IlwisObject::create<IDomain>("value"), true);
	ColumnDefinition& coldef2 = _outputTable->columndefinitionRef("TotalUpstreamArea");
	coldef2.datadef().range(new NumericRange(0, 1.0e300, 0.01));



}


void MapCatchmentMerge::CreateTableSegmentsExtracted()
{
	ITable tblAtt = _inDrngOrderRaster->attributeTable();

	std::vector<QVariant> cUpstreamLinkSrc = tblAtt->column("UpstreamLinkID");
	std::vector<QVariant> cUpstreamCoordSrc = tblAtt->column("UpstreamCoord");
	std::vector<QVariant> cUpstreamHeightSrc = tblAtt->column("UpstreamElevation");
	std::vector<QVariant> cDownstreamLinkSrc = tblAtt->column("DownstreamLinkID");
	std::vector<QVariant> cDownstreamCoordSrc = tblAtt->column("DownstreamCoord");
	std::vector<QVariant> cDownstreamHeightSrc = tblAtt->column("DownstreamElevation");
	std::vector<QVariant> cDropSrc = tblAtt->column("ElevationDifference");
	std::vector<QVariant> cStrahlerSrc = tblAtt->column("Strahler");
	std::vector<QVariant> cStrahlerClassSrc = tblAtt->column("StrahlerClass");
	std::vector<QVariant> cStreveSrc = tblAtt->column("Shreve");
	std::vector<QVariant> cLengthSrc = tblAtt->column("Length");
	std::vector<QVariant> cStraightLengthSrc = tblAtt->column("StraightLength");
	std::vector<QVariant> cSlopAlongDrainageSrc = tblAtt->column("SlopeAlongDrainagePerc");
	std::vector<QVariant> cSlopAlongDrainageDegreeSrc = tblAtt->column("SlopeAlongDrainageDegree");
	std::vector<QVariant> cSlopDrainageStraightSrc = tblAtt->column("SlopeDrainageStraightPerc");
	std::vector<QVariant> cSlopDrainageStraightDegreeSrc = tblAtt->column("SlopeDrainageStraightDegree");
	std::vector<QVariant> cSinuositySrc = tblAtt->column("Sinuosity");
	std::vector<QVariant> cTotalUpstreamLengthSrc = tblAtt->column("TotalUpstreamAlongDrainageLength");


	IFlatTable newSegTable;
	newSegTable.prepare();

	newSegTable->addColumn("Catchment", IlwisObject::create<IDomain>("value"), true);
	ColumnDefinition& colcatchment = newSegTable->columndefinitionRef("Catchment");
	colcatchment.datadef().range(new NumericRange(1, 32767, 1));

	newSegTable->addColumn("UpstreamLinkID", IlwisObject::create<IDomain>("text"), true);
	newSegTable->addColumn("UpstreamElevation", IlwisObject::create<IDomain>("value"), true);
	newSegTable->addColumn("DownstreamLinkID", IlwisObject::create<IDomain>("value"), true);
	newSegTable->addColumn("DownstreamElevation", IlwisObject::create<IDomain>("value"), true);
	newSegTable->addColumn("ElevationDifference", IlwisObject::create<IDomain>("value"), true);
	newSegTable->addColumn("Strahler", IlwisObject::create<IDomain>("value"), true);
	newSegTable->addColumn("Shreve", IlwisObject::create<IDomain>("value"), true);
	newSegTable->addColumn("Length", IlwisObject::create<IDomain>("value"), true);
	newSegTable->addColumn("StraightLength", IlwisObject::create<IDomain>("value"), true);
	newSegTable->addColumn("SlopeAlongDrainagePerc", IlwisObject::create<IDomain>("value"), true);
	newSegTable->addColumn("SlopeAlongDrainageDegree", IlwisObject::create<IDomain>("value"), true);
	newSegTable->addColumn("SlopeDrainageStraightPerc", IlwisObject::create<IDomain>("value"), true);
	newSegTable->addColumn("SlopeDrainageStraightDegree", IlwisObject::create<IDomain>("value"), true);
	newSegTable->addColumn("Sinuosity", IlwisObject::create<IDomain>("value"), true);
	newSegTable->addColumn("TotalUpstreamAlongDrainageLength", IlwisObject::create<IDomain>("value"), true);

	IDomain strahlerDom = tblAtt->columndefinitionRef("StrahlerClass").datadef().domain();
	newSegTable->addColumn("StrahlerClass", strahlerDom );

	_outputExtractSegTable = newSegTable;

	ICoordinateDomain crddom0;
	crddom0.prepare();
	crddom0->setCoordinateSystem(_csy);
	_outputExtractSegTable->addColumn("UpstreamCoord", crddom0, true);

	ICoordinateDomain crddom1;
	crddom1.prepare();
	crddom1->setCoordinateSystem(_csy);
	_outputExtractSegTable->addColumn("DownstreamCoord", crddom1, true);

	_outputExtraxtedSegmentMap->setAttributes(_outputExtractSegTable);

	PixelIterator iterDem = PixelIterator(_inDemRaster, BoundingBox(), PixelIterator::fXYZ);
	bool IsExists(false);
	long iRecs = -1;
    for (int i = 0; i < m_extractedsegids.size(); i++)
	{
		long segnum = m_extractedsegids[i].size();
		for (int j = 0; j<segnum; j++)
		{
            long segid = m_extractedsegids[i][j]-1;
			iRecs++;
            _outputExtractSegTable->setCell("Catchment", iRecs, QVariant(i + 1));
			_outputExtractSegTable->setCell("UpstreamLinkID", iRecs, QVariant(cUpstreamLinkSrc[segid].toString()));

			QVariant c1a;
			Coordinate crd = cUpstreamCoordSrc[segid].value<Coordinate>();
			crd.z = 0;
			c1a.setValue(crd);
			_outputExtractSegTable->setCell("UpstreamCoord", iRecs, c1a);

			_outputExtractSegTable->setCell("UpstreamElevation", iRecs, QVariant(cUpstreamHeightSrc[segid].toDouble()));

			_outputExtractSegTable->setCell("DownstreamLinkID", iRecs, QVariant(cDownstreamLinkSrc[segid].toInt()));

			QVariant c1a1;
			Coordinate crd1 = cDownstreamCoordSrc[segid].value<Coordinate>();
			crd1.z = 0;
			c1a1.setValue(crd1);
			_outputExtractSegTable->setCell("DownstreamCoord", iRecs, c1a1);

			_outputExtractSegTable->setCell("DownstreamElevation", iRecs, QVariant(cDownstreamHeightSrc[segid].toDouble()));

			_outputExtractSegTable->setCell("ElevationDifference", iRecs, QVariant(cDropSrc[segid].toDouble()));
			_outputExtractSegTable->setCell("Strahler", iRecs, QVariant(cStrahlerSrc[segid].toInt()));

			int a = cStrahlerClassSrc[segid].toInt();
			QString stra = QString("%1").arg(cStrahlerClassSrc[segid].toInt());

			_outputExtractSegTable->setCell("StrahlerClass", iRecs, QVariant(QString("%1").arg(1 + cStrahlerClassSrc[segid].toInt())));

			_outputExtractSegTable->setCell("Shreve", iRecs, QVariant(cStreveSrc[segid].toDouble()));
			_outputExtractSegTable->setCell("Length", iRecs, QVariant(cLengthSrc[segid].toDouble()));

			_outputExtractSegTable->setCell("StraightLength", iRecs, QVariant(cStraightLengthSrc[segid].toDouble()));
			_outputExtractSegTable->setCell("SlopeAlongDrainagePerc", iRecs, QVariant(cSlopAlongDrainageSrc[segid].toDouble()));
			_outputExtractSegTable->setCell("SlopeAlongDrainageDegree", iRecs, QVariant(cSlopAlongDrainageDegreeSrc[segid].toDouble()));
			_outputExtractSegTable->setCell("SlopeDrainageStraightPerc", iRecs, QVariant(cSlopDrainageStraightSrc[segid].toDouble()));
			_outputExtractSegTable->setCell("SlopeDrainageStraightDegree", iRecs, QVariant(cSlopDrainageStraightDegreeSrc[segid].toDouble()));
			_outputExtractSegTable->setCell("Sinuosity", iRecs, QVariant(cSinuositySrc[segid].toDouble()));
			_outputExtractSegTable->setCell("TotalUpstreamAlongDrainageLength", iRecs, QVariant(cTotalUpstreamLengthSrc[segid].toDouble()));

			if (m_UseOutlets)
			{
				for (std::vector<OutletLocation>::iterator pos = m_vOutlet.begin(); pos < m_vOutlet.end(); ++pos)
				{
					Pixel rc = pos->pxl;
					OutletLocation ol = (*pos);
					if (ol.StreamID == m_extractedsegids[i][j] && ol.isOnNode != true)
					{
						Coordinate cdDownstream = _inDrngOrderRaster->georeference()->pixel2Coord(rc);
						QVariant c1a;
						c1a.setValue(cdDownstream);
						_outputExtractSegTable->setCell("DownstreamCoord", iRecs, c1a);

						double rDownstreamHeight = *(iterDem)(rc);

						_outputExtractSegTable->setCell("DownstreamElevation", iRecs, QVariant(rDownstreamHeight));
						_outputExtractSegTable->setCell("Length", iRecs, QVariant(ol.rLen1));
						Coordinate upcrdstr = cUpstreamCoordSrc[segid].value<Coordinate>();
						double rStraightLength = rDistance(upcrdstr, cdDownstream);
						_outputExtractSegTable->setCell("StraightLength", iRecs, QVariant(rStraightLength));

						double rDrop = cUpstreamHeightSrc[segid].toDouble() - rDownstreamHeight;
						_outputExtractSegTable->setCell("ElevationDifference", iRecs, QVariant(rDrop));

						double rSlop = rComputeSlope(rDrop, ol.rLen1, false);
						_outputExtractSegTable->setCell("SlopeAlongDrainagePerc", iRecs, QVariant(rSlop));

						rSlop = rComputeSlope(rDrop, ol.rLen1, true);
						_outputExtractSegTable->setCell("SlopeAlongDrainageDegree", iRecs, QVariant(rSlop));

						rSlop = rComputeSlope(rDrop, rStraightLength, false);
						_outputExtractSegTable->setCell("SlopeDrainageStraightPerc", iRecs, QVariant(rSlop));

						rSlop = rComputeSlope(rDrop, rStraightLength, true);
						_outputExtractSegTable->setCell("SlopeDrainageStraightDegree", iRecs, QVariant(rSlop));

						double rSinuosity = rComputeSinuosity(ol.rLen1, rStraightLength);
						_outputExtractSegTable->setCell("Sinuosity", iRecs, QVariant(rSinuosity));

					}
				}
			}

		}
	}

}



std::vector<long> MapCatchmentMerge::GetSegmentIDsExtracted( int iCatch )
{
	std::vector<long> vSegIDs;
	vSegIDs.clear();
	std::vector<QVariant> drgIDs = _outputTable->column("DrainageID");
	long iSize = drgIDs.size();

	if (iCatch >= iSize)
		return vSegIDs;

		if (drgIDs[iCatch].isValid())
		SplitString(drgIDs[iCatch].toString(), QString(","), vSegIDs);

	return vSegIDs;
}


void  MapCatchmentMerge::ExtractSegments()
{
	//Remove the part of segments splited by the outlets, if needed
	CleanSegment(_outputExtraxtedSegmentMap, _internelDranageSegmentMap);
}


void MapCatchmentMerge::CleanSegment(IFeatureCoverage smpTo, IFeatureCoverage smpFrom)
{
	m_extractedsegids.clear();
	std::vector<Pixel> vOutlet;
	std::vector<long> vOutletID;
	//check if the outlet has a downflow to catchment merged
	//if no, the rest of the segment at the outlet should be removed from
	//otherwise, just copy the complete segment
	for (std::vector<OutletLocation>::iterator pos = m_vOutlet.begin(); pos < m_vOutlet.end(); ++pos)
	{
		Pixel rc = pos->pxl;
		int iFlow = GetDownStreamCell(rc);
		if (*m_iterOut[rc] == iUNDEF)
		{
			vOutlet.push_back(pos->pxl);
			vOutletID.push_back(pos->StreamID);
		}
	}

	std::vector<QVariant> drgIDs = _outputTable->column("DrainageID");
	long iSize = drgIDs.size();

	m_extractedsegids.resize(iSize);
	for (int i = 0; i < iSize; i++)
	{
		std::vector<long> vSegIDs = GetSegmentIDsExtracted(i);
		std::map<long, long> mapSrcDstIDs;
		////////////////////
		for (auto feature : smpFrom)
		{
			const geos::geom::LineString* ls_geom =
				dynamic_cast<const geos::geom::LineString*>(feature->geometry().get());

			Record rec = feature->record();
			QVariant val = rec.cell(rec.columnCount() - 1);
			int ftid = val.toInt() + 1;

			bool IsExists = find(vSegIDs.begin(), vSegIDs.end(), ftid) != vSegIDs.end();

			if (IsExists && ls_geom && !ls_geom->isEmpty())
			{
				//			long iSegVal = ls_geom->getSRID(); // some iSegVals will be duplicates; only add an item to the domain if its not a duplicate
				long iSegVal = ftid; // some iSegVals will be duplicates; only add an item to the domain if its not a duplicate
				long iDestSegVal;
				std::map<long, long>::iterator item = mapSrcDstIDs.find(iSegVal);
				if (item == mapSrcDstIDs.end()) {
					iDestSegVal = mapSrcDstIDs.size() + 1;
					mapSrcDstIDs[iSegVal] = iDestSegVal;
				}
				else
				{
					iDestSegVal = item->second;
				}

				long iNumC2 = 1;
				geos::geom::CoordinateSequence* crdbufFrom = ls_geom->getCoordinates();
				geos::geom::CoordinateSequence* crdbufnew = crdbufFrom->clone();

				std::vector<long>::iterator pos = find(vOutletID.begin(), vOutletID.end(), iSegVal);
				if (pos != vOutletID.end())
				{
					int iIndex = (long)(pos - vOutletID.begin());
					Pixel rcFrom = m_vDrainageAtt[iSegVal - 1].UpstreamCoord;
					Pixel rcOutlet = vOutlet[iIndex];
					for (int i = 0; i < crdbufFrom->size(); ++i)
					{
						int iFlow = GetDownStreamCell(rcFrom);
						if (rcFrom != rcOutlet)
							iNumC2++;
						else
						{
							iNumC2++;  //fix to include the outlet point
							break;
						}
					}
				}
				else
					iNumC2 = crdbufFrom->size();

				std::vector<geos::geom::Coordinate>* coords = new std::vector<geos::geom::Coordinate>();

				for (int j = 0; j < iNumC2; j++)
					coords->push_back(crdbufnew->getAt(j));

				geos::geom::CoordinateArraySequence* points = new  geos::geom::CoordinateArraySequence(coords);
				geos::geom::Geometry* geometry = smpTo->geomfactory()->createLineString(points);

				if (geometry->isValid())
				{
					SPFeatureI ft = smpTo->newFeature(geometry);
					m_extractedsegids[i].push_back(iSegVal);
				}
			}
		}


	}


}


bool sortID(const int a, const int b)
{
	return a < b;
}

void MapCatchmentMerge::ComputeOtherAttributes()
{
	//Add fields to the table for the output catchments

	_outputTable->addColumn("TotalDrainageLength", IlwisObject::create<IDomain>("value"), true);
	ColumnDefinition& coldef0 = _outputTable->columndefinitionRef("TotalDrainageLength");
	coldef0.datadef().range(new NumericRange(0, 1.0e300, 0.01));

	_outputTable->addColumn("DrainageDensity", IlwisObject::create<IDomain>("value"), true);
	ColumnDefinition& coldef1 = _outputTable->columndefinitionRef("DrainageDensity");
	coldef1.datadef().range(new NumericRange(0, 1.0e300, 0.01));

	_outputTable->addColumn("LongestFlowPathLength", IlwisObject::create<IDomain>("value"), true);
	ColumnDefinition& coldef2 = _outputTable->columndefinitionRef("LongestFlowPathLength");
	coldef2.datadef().range(new NumericRange(0, 1.0e300, 0.01));

	_outputTable->addColumn("LongestDrainageLength", IlwisObject::create<IDomain>("value"), true);
	ColumnDefinition& coldef3 = _outputTable->columndefinitionRef("LongestDrainageLength");
	coldef3.datadef().range(new NumericRange(0, 1.0e300, 0.01));

	ICoordinateDomain crddom1;
	crddom1.prepare();
	crddom1->setCoordinateSystem(_csy);
	_outputTable->addColumn("CenterDrainage", crddom1, true);

	ICoordinateDomain crddom2;
	crddom2.prepare();
	crddom2->setCoordinateSystem(_csy);
	_outputTable->addColumn("OutletCoord", crddom2, true);


	_outputTable->addColumn("OutletElevation", IlwisObject::create<IDomain>("value"), true);
	ColumnDefinition& coldef4 = _outputTable->columndefinitionRef("OutletElevation");
	coldef4.datadef().range(new NumericRange(0, 1.0e300, 0.01));

	ICoordinateDomain crddom3;
	crddom3.prepare();
	crddom3->setCoordinateSystem(_csy);
	_outputTable->addColumn("LFPUpstreamCoord", crddom3, true);

	_outputTable->addColumn("LFPUpstreamElevation", IlwisObject::create<IDomain>("value"), true);
	ColumnDefinition& coldef5 = _outputTable->columndefinitionRef("LFPUpstreamElevation");
	coldef5.datadef().range(new NumericRange(0, 1.0e300, 0.01));

	ICoordinateDomain crddom4;
	crddom4.prepare();
	crddom4->setCoordinateSystem(_csy);
	_outputTable->addColumn("LDPUpstreamCoord", crddom4, true);

	_outputTable->addColumn("LDPUpstreamElevation", IlwisObject::create<IDomain>("value"), true);
	ColumnDefinition& coldef6 = _outputTable->columndefinitionRef("LDPUpstreamElevation");
	coldef6.datadef().range(new NumericRange(0, 1.0e300, 0.01));


	//Retrieve attributes from catchment table

	std::vector<QVariant> colDrainageId;
	colDrainageId = _outputTable->column("DrainageID");

	std::vector<QVariant> colCatchArea;
	colCatchArea = _outputTable->column("CatchmentArea");

	std::vector<QVariant> colDrainageLen;
	colDrainageLen = _outputTable->column("DrainageLen");

	//Retrieve attributes from drainage table 
	ITable tblAtt = _inDrngOrderRaster->attributeTable();

	long iSize = tblAtt->recordCount();

	std::vector<long> pdsrt;

	for (int j = 0; j < iSize; j++)
		pdsrt.push_back(j);

	std::vector<QVariant> colFlowLength;
	colFlowLength = tblAtt->column("Length");

	std::vector<QVariant> colDownLinkID;
	colDownLinkID = tblAtt->column("DownstreamLinkID");

	std::vector<QVariant> colUpCoords;
	colUpCoords = tblAtt->column("UpstreamCoord");

	std::vector<QVariant> colDownCoords;
	colDownCoords = tblAtt->column("DownstreamCoord");

	std::vector<QVariant> colToCoords;
	colToCoords = tblAtt->column("TostreamCoord");


	long iCatchments = colDrainageId.size();

	PixelIterator iterDem = PixelIterator(_inDemRaster, BoundingBox(), PixelIterator::fXYZ);

	AttLongestPath FlowPathAtt;
	std::vector<AttLongestPath> vFlowPathAtt;

	double rHeight1;
	double rHeight2;
	int rec = 0;
	for (long i = 1; i <= iCatchments; i++)
	{
		geos::geom::CoordinateSequence* longestcrdbuf = 0;
		longestcrdbuf = _longestPathFeature->geomfactory()->getCoordinateSequenceFactory()->create();

		rec = i - 1;
		double rTotlaDrainageLen = 0;
		double rLenDownDrainages;
		std::vector<double> vLenDownDrainages;
		vLenDownDrainages.resize(0);

		//	//Get drainages per catchment
		std::vector<long> vDrainageIDs;
		vDrainageIDs.clear();
		SplitString(colDrainageId[rec].toString(), QString(","), vDrainageIDs);
		if (vDrainageIDs.size() <= 0)
			continue;

		//	//Retrieve lengths per drainage per catchment  
		std::vector<double> vLenPerDrainage;
		SplitString(colDrainageLen[rec].toString(), QString(","), vLenPerDrainage);
		if (vLenPerDrainage.size() <= 0)
			continue;

		long iRaw;
		//	//Compute the total drainage length per catchment
		//	//Compute down stream drainage length per drainage 
		for (std::vector<long>::iterator pos = vDrainageIDs.begin(); pos < vDrainageIDs.end(); ++pos)
		{
			int dnopos = *pos - 1;  // vDrainageIDs in iwlis 3.8 starts from 1, here starts from 0
			iRaw = pdsrt[dnopos];
			if (colFlowLength[iRaw].toDouble() == rUNDEF)
				break;

			//The sum of the drainage lengths for each drainage line in the catchment rTotlaDrainageLen
			//Compute the total lengths of the up drainages till to the outlet rLenUpDrainages 
			long iIndex = (long)(pos - vDrainageIDs.begin());
			if ( vLenPerDrainage[iIndex] > 0.0 )
			{
				rTotlaDrainageLen += vLenPerDrainage[iIndex];
				rLenDownDrainages = vLenPerDrainage[iIndex];
			}
			else
			{
				rTotlaDrainageLen += colFlowLength[iRaw].toDouble();
				rLenDownDrainages = colFlowLength[iRaw].toDouble();
			}
			long iDownFlowID = colDownLinkID[iRaw].toInt();
			bool IsExist;
			if (iDownFlowID != iUNDEF)
			{
				do
				{
					std::vector<long>::iterator posID = find(vDrainageIDs.begin(), vDrainageIDs.end(), iDownFlowID);
					IsExist = (posID != vDrainageIDs.end());
					if (IsExist)
					{
						/*	sLbl = String("%li", iDownFlowID);
							iRaw = pdsrtDrainage->iOrd(sLbl);*/
						iRaw = pdsrt[iDownFlowID - 1];
						iIndex = (long)(posID - vDrainageIDs.begin());
						if ( vLenPerDrainage[iIndex] > 0.0 )
								rLenDownDrainages += vLenPerDrainage[iIndex];
						else
							rLenDownDrainages += colFlowLength[iRaw].toDouble();
						iDownFlowID = colDownLinkID[iRaw].toInt();
					}

				} while (IsExist && iDownFlowID != iUNDEF);
			}
			vLenDownDrainages.push_back(rLenDownDrainages);
		}

		std::vector<double>::iterator posMax = max_element(vLenDownDrainages.begin(), vLenDownDrainages.end());
		double rLongestDrainageLength = (*posMax);

		//	//Update longest drainage length, total drainage length, drainage density attributes
		// 
		_outputTable->setCell("LongestDrainageLength", rec, QVariant(rLongestDrainageLength));
		_outputTable->setCell("TotalDrainageLength", rec, QVariant(rTotlaDrainageLen));

		//////////////////////////////
		///  wait for polygon function
		//double rDensity = (rTotlaDrainageLen / colCatchArea[iRaw].toDouble()) * 1000000;
		double rDensity(.0);
		_outputTable->setCell("DrainageDensity", rec, QVariant(rDensity));

		////////////////////////

		//	//Find the headwater drainage for the longest drainage path
		long iIndex = (long)(posMax - vLenDownDrainages.begin());
		long iSourceID = vDrainageIDs[iIndex];

		//	String sLbl("%li", iSourceID);
		//	iRaw = pdsrtDrainage->iOrd(sLbl);

		iRaw = pdsrt[iSourceID - 1];
		Coordinate crddown = colDownCoords[iRaw].value<Coordinate>();
		crddown.z = 0;

		Coordinate crdup = colUpCoords[iRaw].value<Coordinate>();
		crdup.z = 0;

		Pixel rcDownstream = _inDrngOrderRaster->georeference()->coord2Pixel(crddown);
		Pixel rcUpstream = _inDrngOrderRaster->georeference()->coord2Pixel(crdup);

		QVariant c1aup;
		c1aup.setValue(crdup);
		_outputTable->setCell("LDPUpstreamCoord", rec, c1aup);

		double rHeight = *iterDem(rcUpstream);
		_outputTable->setCell("LDPUpstreamElevation", rec, QVariant(rHeight));

		if (*m_iterOut(rcDownstream) != i && m_UseOutlets)
		{
			Pixel rc = rcUpstream;
			Pixel rcDownCell;
			while (*m_iterOut(rc) == i)
			{
				GetDownStreamCell(rc);
			}
			rcDownstream = rc;
		}
		//	//Construct headwater drainage segment per catchment at downstream cell
		InitPars();
		m_vStream.resize(0);
		m_rSourceWaterFlowPathLen = 0;

		m_vStream.push_back(_inDrngOrderRaster->georeference()->pixel2Coord(rcDownstream));

		lastrc = Pixel(0, 0, 0);

		if (m_vStream.size() == 1)
			rHeight1 = *iterDem(rcDownstream);

		if (m_vOutlet[i - 1].fExtractOverlandFlowPath)
			ExtractUpstreamFlowPath(rcDownstream, i);
		else
		{
			m_rSourceWaterFlowPathLen = rLongestDrainageLength;
			long streamid = m_vOutlet[i - 1].StreamID;
			m_vStream.push_back(_inDrngOrderRaster->georeference()->pixel2Coord(m_vDrainageAtt[streamid - 1].DownStreamCoord));
			m_vStream.push_back(_inDrngOrderRaster->georeference()->pixel2Coord(m_vDrainageAtt[streamid - 1].UpstreamCoord));
			iRaw = pdsrt[streamid - 1];
			Coordinate cdDownstreamCoord = StoreSegment(_internelDranageSegmentMap, longestcrdbuf, iRaw, i);
		}

		rHeight2 = *iterDem(lastrc);

		if (rHeight2 > rHeight1)
			std::reverse(m_vStream.begin(), m_vStream.end());


		//Compute the longest flow path length
		double rLongestFlowPathLength;

		if ( vLenPerDrainage[iIndex] > 0.0 )
				rLongestFlowPathLength = rLongestDrainageLength - vLenPerDrainage[iIndex] + m_rSourceWaterFlowPathLen;
		else
			rLongestFlowPathLength = rLongestDrainageLength - colFlowLength[iRaw].toDouble() + m_rSourceWaterFlowPathLen;

		_outputTable->setCell("LongestFlowPathLength", rec, QVariant(rLongestFlowPathLength));

		Pixel vstream = _inDrngOrderRaster->georeference()->coord2Pixel(m_vStream[m_vStream.size() - 1]);
		vstream.z = 0;

		Coordinate crdLFP = _inDrngOrderRaster->georeference()->pixel2Coord(vstream);
		QVariant c1aLFPCoord;
		c1aLFPCoord.setValue(crdLFP);
		_outputTable->setCell("LFPUpstreamCoord", rec, c1aLFPCoord);

		rHeight = *iterDem(vstream);
		_outputTable->setCell("LFPUpstreamElevation", rec, QVariant(rHeight));

		Coordinate crdCenterDrainage = ComputeCenterDrainage(iSourceID, rLongestFlowPathLength / 2, _internelDranageSegmentMap);

		QVariant c1a;
		c1a.setValue(crdCenterDrainage);
		_outputTable->setCell("CenterDrainage", rec, c1a);

		QVariant c1aoutlet;
		c1aoutlet.setValue(_inDrngOrderRaster->georeference()->pixel2Coord(m_vOutlet[i - 1].pxl));
		_outputTable->setCell("OutletCoord", rec, c1aoutlet);

		Pixel pxlpos = m_vOutlet[i - 1].pxl;
		pxlpos.z = 0;
		rHeight = *iterDem(pxlpos);
		_outputTable->setCell("OutletElevation", rec, QVariant(rHeight));

		// generate longest path
		FlowPathAtt.UpstreamCoord = m_vStream[m_vStream.size() - 1];
		FlowPathAtt.DownstreamCoord = m_vStream[0];

		if (m_vOutlet[i - 1].fExtractOverlandFlowPath)
			StoreSourceSegment(longestcrdbuf, i);
		double rFlowPathLen = m_rSourceWaterFlowPathLen;
		long iDownFlowID = colDownLinkID[iRaw].toInt();
		bool IsExist;
		if (iDownFlowID != iUNDEF)
		{
			do
			{
				std::vector<long>::iterator posID = find(vDrainageIDs.begin(), vDrainageIDs.end(), iDownFlowID);
				IsExist = (posID != vDrainageIDs.end());
				if (IsExist)
				{
					/*sLbl = String("%li", iDownFlowID);
					iRaw = pdsrtDrainage->iOrd(sLbl);*/
					iRaw = pdsrt[iDownFlowID - 1];
					iIndex = (long)(posID - vDrainageIDs.begin());
					if ((m_UseOutlets == true) && (vLenPerDrainage[iIndex] > 0.0))
					{
						rFlowPathLen += vLenPerDrainage[iIndex];
						Coordinate crdup = colUpCoords[iRaw].value<Coordinate>();
						crdup.z = 0;
						Pixel rcUpstream = _inDrngOrderRaster->georeference()->coord2Pixel(crdup);
						FlowPathAtt.DownstreamCoord = SplitSegment(_internelDranageSegmentMap, longestcrdbuf, iRaw, rLongestFlowPathLength, i, rcUpstream);
					}
					else
					{
						FlowPathAtt.DownstreamCoord = StoreSegment(_internelDranageSegmentMap, longestcrdbuf, iRaw, i);
						rFlowPathLen += colFlowLength[iRaw].toDouble();
					}
					iDownFlowID = colDownLinkID[iRaw].toInt();
				}
			} while (IsExist && iDownFlowID != iUNDEF);
		}
		FlowPathAtt.rLength = rLongestFlowPathLength;
		vFlowPathAtt.push_back(FlowPathAtt);

		if (longestcrdbuf && longestcrdbuf->getSize() > 0)
		{
			bool dirdown(false);
			Coordinate crdolt = _inDrngOrderRaster->georeference()->pixel2Coord(m_vOutlet[i - 1].pxl);
			crdolt.z = 0;

			double mindis = 1e307;
			int minpos(0);
			for (int j = longestcrdbuf->getSize() - 1; j >= 0; j--)
			{
				Coordinate crd = longestcrdbuf->getAt(j);
				crd.z = 0;
				double dis = rDistance(crdolt, crd);
				if (dis < mindis)
				{
					mindis = dis;
					minpos = j;
				}
			}

			Coordinate mindiscrd = longestcrdbuf->getAt(minpos);
			double otletht = *iterDem[m_vOutlet[i - 1].pxl];
			std::vector<double> hts;
			for (int j = longestcrdbuf->getSize()-1; j >=0; j--)
			{
				if ( j > minpos )
					longestcrdbuf->deleteAt(longestcrdbuf->getSize() - 1);
			}

			Pixel minpospxl = _inDrngOrderRaster->georeference()->coord2Pixel(mindiscrd);
			double minposht = *iterDem(minpospxl);

			if (minposht <= otletht)
				longestcrdbuf->add(crdolt);

			geos::geom::Geometry* geometry = _longestPathFeature->geomfactory()->createLineString(longestcrdbuf);

			if (geometry->isValid())
			{
				SPFeatureI ft = _longestPathFeature->newFeature(geometry);
			}
		}
	}

	CreateTableLongestFlowPath(vFlowPathAtt);

}


void MapCatchmentMerge::CreateTableLongestFlowPath(std::vector<AttLongestPath> vAtt)
{

	IFlatTable newPathTable;
	newPathTable.prepare();

	ICoordinateDomain crddom1;
	crddom1.prepare();
	crddom1->setCoordinateSystem(_csy);
	newPathTable->addColumn("UpstreamCoord", crddom1, true);

	ICoordinateDomain crddom2;
	crddom2.prepare();
	crddom2->setCoordinateSystem(_csy);
	newPathTable->addColumn("DownstreamCoord", crddom2, true);
	newPathTable->addColumn("Length", IlwisObject::create<IDomain>("value"), true);
	newPathTable->addColumn("StraightLength", IlwisObject::create<IDomain>("value"), true);

	newPathTable->addColumn("Sinuosity", IlwisObject::create<IDomain>("value"), true);
	ColumnDefinition& coldef0 = newPathTable->columndefinitionRef("Sinuosity");
	coldef0.datadef().range(new NumericRange(0, 1e10, 0.001));

	_outputPathTable = newPathTable;

	long i = 0;
	for (std::vector<AttLongestPath>::iterator pos = vAtt.begin(); pos < vAtt.end(); ++pos)
	{
		//i++;
		AttLongestPath atts = (*pos);
		QVariant ucrd;
		ucrd.setValue(atts.UpstreamCoord);
		QVariant dcrd;
		dcrd.setValue(atts.DownstreamCoord);

		_outputPathTable->setCell("UpstreamCoord", i, ucrd);
		_outputPathTable->setCell("DownstreamCoord", i, dcrd);
		_outputPathTable->setCell("Length", i, QVariant(atts.rLength));
		double rStraightLength = rDistance(atts.UpstreamCoord, atts.DownstreamCoord);
		_outputPathTable->setCell("StraightLength", i, QVariant(rStraightLength));

		double rSinuousity = rComputeSinuosity(atts.rLength, rStraightLength);
		_outputPathTable->setCell("Sinuosity", i, QVariant(rSinuousity));

		i++;
	}

	_longestPathFeature->setAttributes(_outputPathTable);
}


double MapCatchmentMerge::rComputeSinuosity(double rLength, double rStraightLenght)
{
	if (rStraightLenght > 0)
		return rLength / rStraightLenght;
	else
		return 0;
}

Coordinate MapCatchmentMerge::SplitSegment(IFeatureCoverage smpFrom, long iRaw, double disval, long id, Pixel rc)
{
	Coordinate cd;

	std::vector< const geos::geom::Geometry*> geomarr;

	for (auto feature : smpFrom)
	{
		const geos::geom::LineString* ls_geom =
			dynamic_cast<const geos::geom::LineString*>(feature->geometry().get());

		Record rec = feature->record();
		QVariant val = rec.cell(rec.columnCount() - 1);
		int ftid = val.toInt();

		if (ls_geom && !ls_geom->isEmpty() && ftid == iRaw)
		{
			geos::geom::CoordinateSequence* crdbufFrom = ls_geom->getCoordinates();
			geos::geom::CoordinateSequence* crdbufnew = crdbufFrom->clone();


			int iFlow = 1;
			bool fCalculate = true;
			long iCount = 0;
			while ((iFlow != 0) && fCalculate)
			{
				fCalculate = IsEdgeCell(rc) != true;

				iFlow = GetDownStreamCell(rc);
				if (*m_iterOut[rc] == id && fCalculate && iFlow != 0)   //valid flow
					iCount++;
				else
					fCalculate = false;
			}

			geos::geom::Geometry* geometry = _longestPathFeature->geomfactory()->createLineString(crdbufnew);
			//lines.push_back(geometry);


			if (geometry->isValid())
			{
				geometry->setSRID(id);
				SPFeatureI ft = _longestPathFeature->newFeature(geometry);
				rec.cell(rec.columnCount() - 1) = (int)id;
				ft->record(rec);
			}
			cd = crdbufFrom->getAt(iCount - 1);
			//delete crdbufFrom;
		}

	}

	return cd;
}


void MapCatchmentMerge::StoreSourceSegment(long val)
{
	if (m_vStream.size() == 0)
		return;

	geos::geom::CoordinateSequence* crdbuf = _longestPathFeature->geomfactory()->getCoordinateSequenceFactory()->create();

	for (unsigned long index = 0; index < m_vStream.size(); ++index)
	{
		crdbuf->add(m_vStream[index]);
	}

	geos::geom::Geometry* geometry = _longestPathFeature->geomfactory()->createLineString(crdbuf);

	//lines.push_back(geometry);

	if (geometry->isValid())
	{
		//geometry->setSRID(val);
		SPFeatureI ft = _longestPathFeature->newFeature(geometry);
		Record rec = ft->record();
		rec.cell(rec.columnCount() - 1) = (int)val;
	}

}


Coordinate MapCatchmentMerge::StoreSegment(IFeatureCoverage smpFrom, long id, long val)
{
	Coordinate cd;

	std::vector< const geos::geom::Geometry*> geomarr;

	for (auto feature : smpFrom)
	{
		const geos::geom::LineString* ls_geom =
			dynamic_cast<const geos::geom::LineString*>(feature->geometry().get());

		Record rec = feature->record();
		QVariant idval = rec.cell(rec.columnCount() - 1);

		int ftid = idval.toInt();

		if (ls_geom && !ls_geom->isEmpty())
		{
			if (ftid == id)
			{
				geos::geom::CoordinateSequence* crdbufFrom = ls_geom->getCoordinates();
				geos::geom::CoordinateSequence* crdbufnew = crdbufFrom->clone();

				geos::geom::Geometry* geometry = _longestPathFeature->geomfactory()->createLineString(crdbufnew);

				if (geometry->isValid())
				{
					geometry->setSRID(val);

					SPFeatureI ft = _longestPathFeature->newFeature(geometry);
					Record rec = ft->record();
					rec.cell(rec.columnCount() - 1) = (int)val;
					Record rec1 = ft->record();

					QVariant idval1 = rec.cell(rec.columnCount() - 1);

				}
				cd = crdbufFrom->getAt(crdbufFrom->size() - 1);
			}
			//delete crdbufFrom;
		}

	}

	return cd;

}



//************************************************************************

Coordinate MapCatchmentMerge::SplitSegment(IFeatureCoverage smpFrom, geos::geom::CoordinateSequence* crdbuf, long iRaw, double disval, long id, Pixel rc)
{
	Coordinate cd;

	std::vector< const geos::geom::Geometry*> geomarr;

	for (auto feature : smpFrom)
	{
		const geos::geom::LineString* ls_geom =
			dynamic_cast<const geos::geom::LineString*>(feature->geometry().get());

		Record rec = feature->record();
		QVariant val = rec.cell(rec.columnCount() - 1);
		int ftid = val.toInt();

		if (ls_geom && !ls_geom->isEmpty() && ftid == iRaw)
		{
			const geos::geom::CoordinateSequence* crdbufFrom = ls_geom->getCoordinates();

			int iFlow = 1;
			bool fCalculate = true;
			long iCount = 0;
			while ((iFlow != 0) && fCalculate)
			{
				fCalculate = IsEdgeCell(rc) != true;

				iFlow = GetDownStreamCell(rc);
				if (*m_iterOut[rc] == id && fCalculate && iFlow != 0)   //valid flow
					iCount++;
				else
					fCalculate = false;
			}


			for (long k = 0; k < crdbufFrom->size() - 1; k++)
				crdbuf->add(crdbufFrom->getAt(k));

			cd = crdbufFrom->getAt(iCount - 1);
			//delete crdbufFrom;
		}

	}

	return cd;
}

void MapCatchmentMerge::StoreSourceSegment(geos::geom::CoordinateSequence* crdbuf, long val)
{
	if (m_vStream.size() == 0 || !crdbuf)
		return;

	for (unsigned long index = 0; index < m_vStream.size(); ++index)
	{
		crdbuf->add(m_vStream[index]);
	}

}


Coordinate MapCatchmentMerge::StoreSegment(IFeatureCoverage smpFrom, geos::geom::CoordinateSequence* crdbuf, long id, long val)
{
	Coordinate cd;

	std::vector< const geos::geom::Geometry*> geomarr;

	for (auto feature : smpFrom)
	{
		const geos::geom::LineString* ls_geom =
			dynamic_cast<const geos::geom::LineString*>(feature->geometry().get());

		Record rec = feature->record();
		QVariant idval = rec.cell(rec.columnCount() - 1);

		int ftid = idval.toInt();

		if (ls_geom && !ls_geom->isEmpty())
		{
			if (ftid == id)
			{
				geos::geom::CoordinateSequence* crdbufFrom = ls_geom->getCoordinates();
				geos::geom::CoordinateSequence* crdbufnew = crdbufFrom->clone();

				for (long k = 0; k < crdbufFrom->size() - 1; k++)
					crdbuf->add(crdbufFrom->getAt(k));

				cd = crdbufFrom->getAt(crdbufFrom->size() - 1);
			}
		}

	}

	return cd;

}


//********************************************************************************************

Coordinate MapCatchmentMerge::ComputeCenterDrainage(long iDrainageID, double rLength, IFeatureCoverage sm)
{
	//Retrieve attributes from drainage table 
	ITable tblAtt = _inDrngOrderRaster->attributeTable();
	std::vector<QVariant> colDownstreamLinkID = tblAtt->column(QString("DownstreamLinkID"));
	std::vector<QVariant> colFlowLength = tblAtt->column(QString("Length"));

	double rLenDrainage = 0;
	Coordinate c1;

	if (m_rSourceWaterFlowPathLen > rLength)
	{
		c1 = m_vStream[0];
		for (unsigned long index = 0; index < m_vStream.size(); ++index)
		{
			Coordinate c2 = m_vStream[index];
			rLenDrainage += rDistance(c1, c2);
			if (rLenDrainage >= rLength)
				return c1;
			else
				c1 = c2;
		}
		return c1;
	}
	bool fCondition = true;
	while (fCondition && iDrainageID != iUNDEF)
	{
		long iRaw = iDrainageID - 1;

		rLenDrainage += colFlowLength[iRaw].toDouble();
		if (rLenDrainage > rLength)
			fCondition = false;
		else
			iDrainageID = colDownstreamLinkID[iRaw].toInt();
	}

	std::vector< const geos::geom::Geometry*> geomarr;

	for (auto feature : sm)
	{
		const geos::geom::LineString* ls_geom =
			dynamic_cast<const geos::geom::LineString*>(feature->geometry().get());

		Record rec = feature->record();
		QVariant val = rec.cell(rec.columnCount() - 1);
		int ftid = val.toInt();

		if (ls_geom && !ls_geom->isEmpty() && ftid == iDrainageID)
		{
			geos::geom::CoordinateSequence* crdbufFrom = ls_geom->getCoordinates();
			geos::geom::CoordinateSequence* crdbufnew = crdbufFrom->clone();

			c1 = crdbufFrom->getAt(0);
			for (long i = 0; i < crdbufnew->size(); ++i)
			{
				Coordinate c2 = crdbufnew->getAt(i);
				rLenDrainage += rDistance(c1, c2);
				if (rLenDrainage >= rLength)
					return c1;
				else
					c1 = c2;
			}
			//delete crdbufnew;
		}
	}

	return c1;
}


void MapCatchmentMerge::ExtractUpstreamFlowPath(Pixel rc, long id)
{
	Pixel pos;
	Pixel rcUpstream;
	rcUpstream.z = 0;
	int iNb = 0;
	long iAcc = 0;

	//	//Get drainages per catchment
	
	std::vector<QVariant> colDrainageId;
	colDrainageId = _outputTable->column("DrainageID");

	std::vector<long> vDrainageIDs;
	vDrainageIDs.clear();
	SplitString(colDrainageId[id-1].toString(), QString(","), vDrainageIDs);

	PixelIterator iterDem = PixelIterator(_inDemRaster, BoundingBox(), PixelIterator::fXYZ);
	PixelIterator iterDrainage = PixelIterator(_inDrngOrderRaster, BoundingBox(), PixelIterator::fXYZ);
	
	PixelIterator iterAcc = PixelIterator(_inAccRaster, BoundingBox(), PixelIterator::fXYZ);

	for (int i = -1; i <= 1; ++i) 			//Evaluate neighbors for cell at rc
	{
		pos.y = rc.y + i;
		for (int j = -1; j <= 1; ++j)
		{
			pos.x = rc.x + j;
			pos.z = 0;
			int iFlow = *m_iterFld(pos);
			int iFlowTo = m_vFlowSelf[iNb];
			if (m_vFlowSelf[iNb] == *m_iterFld(pos))
			{
				if (*iterAcc(pos) > iAcc && *m_iterOut(pos) == id)
				{
					iAcc = *iterAcc(pos);
					rcUpstream.y = pos.y;
					rcUpstream.x = pos.x;
				}
			}
			iNb++;
		}
	}

	if ((rcUpstream.isValid()) && (*m_iterOut(rcUpstream) == id))
	{
		
		long iDrainageID = *iterDrainage(rcUpstream);
		std::vector<long>::iterator posID = find(vDrainageIDs.begin(), vDrainageIDs.end(), iDrainageID);
		bool IsExist = (posID != vDrainageIDs.end());
		if (IsExist)
		{
			Coordinate c1 = _inDrngOrderRaster->georeference()->pixel2Coord(rc);
			Coordinate c2 = _inDrngOrderRaster->georeference()->pixel2Coord(rcUpstream);
			m_rSourceWaterFlowPathLen += rDistance(c1, c2);
			m_vStream.push_back(c2);
			lastrc = rcUpstream;
			ExtractUpstreamFlowPath(rcUpstream, id);
		}
	}

}



void MapCatchmentMerge::InitPars()
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

	m_vFlowSelf.resize(9);
	m_vFlowSelf[0] = 2;
	m_vFlowSelf[1] = 3;
	m_vFlowSelf[2] = 4;
	m_vFlowSelf[3] = 1;
	m_vFlowSelf[4] = 0;
	m_vFlowSelf[5] = 5;
	m_vFlowSelf[6] = 8;
	m_vFlowSelf[7] = 7;
	m_vFlowSelf[8] = 6;
}


void MapCatchmentMerge::SplitString(QString s, QString mid, std::vector<long>& results)
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


void MapCatchmentMerge::SplitString(QString s, QString mid, std::vector<double>& results)
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


void MapCatchmentMerge::ExtractOriginalOrder(long& id)
{
	//PixelIterator iterInDrng = PixelIterator(_inDrngOrderRaster, BoundingBox(), PixelIterator::fXYZ);

	bool fContinue = true;
	while (fContinue)
	{
		fContinue = false;
		long iCont = 0;
		for (std::vector<DrainageAtt>::iterator pos = m_vDrainageAtt.begin(); pos < m_vDrainageAtt.end(); ++pos)
		{
			DrainageAtt datt = (*pos);  //This is the current drainage to be evaluated
			std::vector<long>  vUpStreamLink = datt.UpstreamID;
			if (datt.CatchmentLink == iUNDEF && IsUpstreamsMerged(vUpStreamLink))
			{
				fContinue = true;
				id++;
				//AddDomainItem(m_dm, id);
				Merge(id, datt, true);
				Pixel rc = datt.DownStreamCoord;
				bool fExtractOverlandFlowPath = true;
				if (datt.UpstreamID[0] != 0)
					fExtractOverlandFlowPath = false;
				PutOutlet(datt.ID, rc, iUNDEF, fExtractOverlandFlowPath);
			}
			iCont++;
		}
	}
}


void MapCatchmentMerge::Merge(long id, DrainageAtt datt, bool fExtractOriginalOrder)
{
	m_vStreamsInCatchment.resize(0);
	m_vUpCatchmentIDs.resize(0);
	m_sUpLinkCatchment = QString("");

	Pixel rc2 = datt.TostreamCoord;

	m_vStreamsInCatchment.push_back(datt.ID);
	m_iOutletVal = datt.ID; //m_vStreamMap[rc2.Row][rc2.Col];
	MergeCatchment(rc2, id, fExtractOriginalOrder);
	UpdateUpLinkCatchment(id);
	UpdateDownLinkCatchment(id);
	OutletLocation ol;
	UpdateLink2StreamSegments(id, ol);
}

bool MapCatchmentMerge::IsUpstreamsMerged(std::vector<long> vUpstreams)
{
	//---Return true, if it's upstreams have been merged  
	bool fCheck = true;
	for (std::vector<long>::iterator pos = vUpstreams.begin(); pos < vUpstreams.end(); ++pos)
	{
		if ((*pos) >= 1 && (m_vDrainageAtt[(*pos) - 1].iOrder > m_iStreamOrders) && (m_vDrainageAtt[(*pos) - 1].CatchmentLink == iUNDEF))
		{
			fCheck = false;
			break;
		}
	}
	return fCheck;
}


void MapCatchmentMerge::UpdateUpLinkCatchment(long id)
{
	if (m_sUpLinkCatchment.length() == 0)
		m_sUpLinkCatchment = "0";
	/*else
		m_sUpLinkCatchment = QString(",%S", m_sUpLinkCatchment);*/
	quint32 record = id - 1;
	_outputTable->setCell("UpstreamLinkCatchment", record, QVariant(m_sUpLinkCatchment));
}

void MapCatchmentMerge::UpdateDownLinkCatchment(long id)
{
	std::vector<long>::iterator pos;
	for (pos = m_vUpCatchmentIDs.begin(); pos < m_vUpCatchmentIDs.end(); ++pos)
	{
		long iUpLinkCatchmentID = *pos;
		quint32 record = iUpLinkCatchmentID - 1;
		_outputTable->setCell("DownstreamLinkCatchment", record, QVariant((int)id));
	}
}


void MapCatchmentMerge::UpdateLink2StreamSegments(long iCatchmentID, OutletLocation ol)
{
	//Update catchment link for each drainage
	std::vector<long>::iterator pos;
	for (pos = m_vStreamsInCatchment.begin(); pos < m_vStreamsInCatchment.end(); ++pos)
	{
		m_vDrainageAtt[(*pos) - 1].CatchmentLink = iCatchmentID;
	}

	//---Catchment table
	QString sStreamInCatchment;
	QString sStreamsInCatchment;
	for (pos = m_vStreamsInCatchment.begin(); pos < m_vStreamsInCatchment.end(); ++pos)
	{
		sStreamInCatchment = sStreamsInCatchment.length() != 0 ?
			QString(",%1").arg((*pos), 0, 10) : QString("%1").arg((*pos), 0, 10);

		if (*pos > 0)
			sStreamsInCatchment = sStreamsInCatchment + sStreamInCatchment;
	}

	quint32 record = iCatchmentID;

	QString id = QString::number(record);

	if (id != "")
		if (_idrange->contains(id))
		{
			++m_totalcount;
			id = QString::number(m_totalcount);
		}
	*_idrange << id;

	record = record - 1;

	_outputTable->setCell("DrainageID", record, QVariant(sStreamsInCatchment));

	//Also update column for length, this column will be used in Horton Plot function
	std::vector<long> vStreamID;
	QString sLength;
	vStreamID.clear();
	SplitString(sStreamsInCatchment, QString(","), vStreamID);

	for (pos = vStreamID.begin(); pos < vStreamID.end(); ++pos)
	{
		QString sLen = "0,";
		long iStreamID = (*pos);
		if (iStreamID == ol.StreamID)
		{
			sLen = QString::number(ol.rLen1);
		}
		else
		{
			double rLen2 = GetSplitSegmentLength(iStreamID);
			if (rLen2 > 0.001)
				sLen = QString::number(ol.rLen1);
		}
		sLength = sLength + sLen;
	}

	_outputTable->setCell("DrainageLen", record, QVariant(sLength));
}


double MapCatchmentMerge::GetSplitSegmentLength(long iStreamID)
{
	for (std::vector<OutletLocation>::iterator pos = m_vOutlet.begin(); pos < m_vOutlet.end(); ++pos)
	{
		OutletLocation ol = (*pos);
		if (iStreamID == ol.StreamID)
		{
			return ol.rLen2;
		}
	}
	return 0;
}


long MapCatchmentMerge::MergeCatchment(Pixel pxl, long iFlag, bool fExtractOriginalOrder)
{
	//****Merge sub-catchments and assigned id iFlag to the merged catchment 
	//For the specified outlet cell in loaction rc, 
	//check whether its neighboring cells flow to it,
	//If true, flag the cells with iFlag in m_vOutput, 
	//the function is called recursively to all the neighboring cells that flow into it. 
	//The recursion stops when it reaches a cell that has no flow to it, or a cell that 
	//has been evaluated.  

	//****Build uplink catchments topology, which are stored in m_sUpLinkCatchment
	//****Store streams in a merged catchment in m_vStreamsInCatchment  

	//location number
	//	-------
	//	|6|7|8|
	//	-------
	//	|5| |1|
	//	-------
	//	|4|3|2|
	//	-------
	// 

	long iFlow = 1;
	pxl.z = 0;
	if (IsEdgeCell(pxl)) return iFlow;
	bool isFlow; //determine if the neighboring cell flows to the cell in location rc 

	/*PixelIterator iterFld = PixelIterator(_inFldRaster, BoundingBox(), PixelIterator::fXYZ);
	PixelIterator iterOut = PixelIterator(_outMergeRaster, BoundingBox(), PixelIterator::fXYZ);
	PixelIterator iterInDrng = PixelIterator(_inDrngOrderRaster, BoundingBox(), PixelIterator::fXYZ);*/

	Pixel pospxl;
	pospxl.z = 0;

	for (int iNr = 1; iNr < 9; iNr++)
	{
		isFlow = false;
		switch (iNr)
		{
		case 1: {	//East
			if (pxl.x != _xsize - 1)
			{
				pospxl.x = pxl.x + 1;
				pospxl.y = pxl.y;

				isFlow = ((*m_iterFld(pospxl) == 5) && (*m_iterOut(pospxl) == iUNDEF));
				BuildUpLinkCatchment(pospxl, 5, iFlag);
			}
		}
			  break;
		case 2: { //South East 
			if (pxl.x != _xsize - 1 && pxl.y != _ysize - 1)
			{
				pospxl.y = pxl.y + 1;
				pospxl.x = pxl.x + 1;
				isFlow = ((*m_iterFld(pospxl) == 6) && (*m_iterOut(pospxl) == iUNDEF));
				BuildUpLinkCatchment(pospxl, 6, iFlag);
			}
		}
			  break;
		case 3: {	//South
			if (pxl.y != _ysize - 1)
			{
				pospxl.y = pxl.y + 1;
				pospxl.x = pxl.x;
				isFlow = ((*m_iterFld(pospxl) == 7) && (*m_iterOut(pospxl) == iUNDEF));
				BuildUpLinkCatchment(pospxl, 7, iFlag);
			}
		}
			  break;
		case 4: { //South West
			if (pxl.x != 0 && pxl.y != _ysize - 1)
			{
				pospxl.y = pxl.y + 1;
				pospxl.x = pxl.x - 1;
				isFlow = ((*m_iterFld(pospxl) == 8) && (*m_iterOut(pospxl) == iUNDEF));
				BuildUpLinkCatchment(pospxl, 8, iFlag);
			}
		}
			  break;
		case 5: {	//West
			if (pxl.x != 0)
			{
				pospxl.x = pxl.x - 1;
				pospxl.y = pxl.y;
				isFlow = ((*m_iterFld(pospxl) == 1) && (*m_iterOut(pospxl) == iUNDEF));
				BuildUpLinkCatchment(pospxl, 1, iFlag);
			}
		}
			  break;
		case 6: {	//North West 
			if (pxl.x != 0 && pxl.y != 0)
			{
				pospxl.y = pxl.y - 1;
				pospxl.x = pxl.x - 1;
				isFlow = ((*m_iterFld(pospxl) == 2) && (*m_iterOut(pospxl) == iUNDEF));
				BuildUpLinkCatchment(pospxl, 2, iFlag);
			}
		}
			  break;
		case 7: {	//North
			if (pxl.y != 0)
			{
				pospxl.y = pxl.y - 1;
				pospxl.x = pxl.x;
				isFlow = ((*m_iterFld(pospxl) == 3) && (*m_iterOut(pospxl) == iUNDEF));
				BuildUpLinkCatchment(pospxl, 3, iFlag);
			}
		}
			  break;
		case 8: {	//North East
			if (pxl.x != _xsize - 1 && pxl.y != 0)
			{
				pospxl.y = pxl.y - 1;
				pospxl.x = pxl.x + 1;
				isFlow = ((*m_iterFld(pospxl) == 4) && (*m_iterOut(pospxl) == iUNDEF));
				BuildUpLinkCatchment(pospxl, 4, iFlag);
			}
		}
			  break;
		}

		if (isFlow)
		{
			if (!fExtractOriginalOrder)
				iFlow += MergeCatchment(pospxl, iFlag, fExtractOriginalOrder);
			else if ((*m_iterInDrng(pospxl) == -2 || *m_iterInDrng(pospxl) == m_iOutletVal))
			{
				iFlow += MergeCatchment(pospxl, iFlag, fExtractOriginalOrder);
			}
		}
		*m_iterOut(pxl) = iFlag;
		//Identify streams in the merged catchment and put stream IDs into m_vStreamsInCatchment
		if (!fExtractOriginalOrder)
			IdentifyStreamsInCatchment(pxl);
	}
	return iFlow;
}


//after merging, a catchment will have more than one streams corresponding with      
//this function will identify the streams in a catchment, put stream's ID into m_vStreamsInCatchment
void MapCatchmentMerge::IdentifyStreamsInCatchment(Pixel pxl)
{
	long iVal = *m_iterInDrng(pxl);
	if (iVal != iUNDEF && iVal != m_iOutletVal && iVal > 0)
	{
		bool IsExists = find(m_vStreamsInCatchment.begin(), m_vStreamsInCatchment.end(), iVal) != m_vStreamsInCatchment.end();
		if (IsExists != true)
			m_vStreamsInCatchment.push_back(iVal);
	}
}


//find the up-link-catchment IDs, 
//store the uplink catchment ID in m_sUpLinkCatchment   
//A catchment has more than one catchment(s) that flow into it.
void MapCatchmentMerge::BuildUpLinkCatchment(Pixel pxl, int iFlow, long iFlag)
{

	long iVal = (long)*(m_iterOut[pxl]);

	if ((*m_iterFld(pxl) == iFlow) && (iVal != iUNDEF) && (iVal != iFlag))
	{
		bool IsExists = find(m_vUpCatchmentIDs.begin(), m_vUpCatchmentIDs.end(), iVal) != m_vUpCatchmentIDs.end();
		if (IsExists != true)
		{
			m_vUpCatchmentIDs.push_back(iVal);
			QString sUpstreamLink = m_sUpLinkCatchment.length() != 0 ?
				QString(",%1").arg(iVal, 0, 10) : QString("%1").arg(iVal, 0, 10);
			m_sUpLinkCatchment = m_sUpLinkCatchment + sUpstreamLink;
		}
	}

}



//---Scan the stream topology for joint outlet locations by order
//---Put the joint outlet locations in m_vOutlet for merging process  
void MapCatchmentMerge::EvaluateJointOutletsbyOrder()
{
	//Catch based on stream connectivity topology info.
	m_vOutlet.resize(0);
	for (std::vector<DrainageAtt>::iterator pos = m_vDrainageAtt.begin(); pos < m_vDrainageAtt.end(); ++pos)
	{
		DrainageAtt datt = (*pos);

		if (datt.iOrder == m_iStreamOrders)
		{
			//---Determine the joint outlet 
			long downstreamID = datt.iDownStreamID;
			//if (downstreamID==339)
			//	long id = m_vDrainageAtt[downstreamID-1].iOrder;
			if (datt.iDownStreamID == iUNDEF || datt.iDownStreamID == 0 || datt.iDownStreamID == shUNDEF)
			{
				PutOutlet(datt.ID, datt.DownStreamCoord, iUNDEF, true);
			}
			else if (m_vDrainageAtt[downstreamID - 1].iOrder > m_iStreamOrders)
			{
				PutOutlet(datt.ID, datt.TostreamCoord, iUNDEF, true);
			}
		}
	}
}



void MapCatchmentMerge::PutOutlet(long id, Pixel rc, long iFlag, bool fExtractOverlandFlowPath)
{
	OutletLocation outlet;
	outlet.StreamID = id;
	outlet.pxl.y = rc.y;
	outlet.pxl.x = rc.x;
	outlet.pxl.z = 0;
	outlet.StreveOrder = iFlag;
	outlet.fExtractOverlandFlowPath = fExtractOverlandFlowPath;
	outlet.isOnNode = false;
	outlet.rLen1 = rUNDEF;
	outlet.rLen2 = rUNDEF;
	//if (IsOutletExists(m_vOutlet, outlet.rc) == false)
	m_vOutlet.push_back(outlet);
}

void MapCatchmentMerge::CreateTable()
{
	//Create a table associated with the merged catchment map
	//A catchment has a unique ID, which will act as a link to relate with
	//its corresponding drainage. 
	IFlatTable newTable;
	newTable.prepare();

	newTable->addColumn("DrainageID", IlwisObject::create<IDomain>("text"), true);
	newTable->addColumn("UpstreamLinkCatchment", IlwisObject::create<IDomain>("text"), true);
	newTable->addColumn("DownstreamLinkCatchment", IlwisObject::create<IDomain>("value"), true);
	ColumnDefinition& coldef0 = newTable->columndefinitionRef("DownstreamLinkCatchment");
	coldef0.datadef().range(new NumericRange(1, 32767, 1));


	//this column will be used to store the length of segments 
	newTable->addColumn("DrainageLen", IlwisObject::create<IDomain>("text"), true);
	_outputTable = newTable;
}


//Put the outlet locations and ite attributes e.g. Streve in to a vector
void MapCatchmentMerge::InitOutletVector()
{
	PixelIterator iterInDrng = PixelIterator(_inDrngOrderRaster, BoundingBox(), PixelIterator::fXYZ);

	ITable tblAtt = _inDrngOrderRaster->attributeTable();

	std::vector<QVariant> colStreve;
	colStreve = tblAtt->column("Shreve");
	OutletLocation ol;

	// we create a pixeliterator to move around on the output raster coverage
	PixelIterator pixiter(_outMergeRaster);
	// create an empty raster (all pixels nodata)
	PixelIterator end = _outMergeRaster->end();
	while (pixiter != end) {
		*pixiter = iUNDEF;
		++pixiter;
	}

	m_iterOut = PixelIterator(_outMergeRaster, BoundingBox(), PixelIterator::fXYZ);

	// we need a coordinate transformation when the two coordinatesystem dont match
	bool needCoordinateTransformation = _inputgrf->coordinateSystem() != _inPointMap->coordinateSystem();
	// loop over all features
	// the 'feature' object is actually something of the SPFeatureI class which is a simple wrapper for a featureinterface
	for (auto feature : _inPointMap)
	{
		OutletLocation ol;
		// we are only interested in points, other geometry types are skipped.
		if (feature->geometryType() != itPOINT)
			continue;
		// get the coordinate of the feature. fortunately points are simple, only one coordinate. Note that we might have multipoints but for
		// the simplicity of the example we ignore that for the moment
		Coordinate crd(*feature->geometry()->getCoordinate());
		// if the coordinate system of input and output dont match we need to transform the input coordinate to something that
		// makes sense in the system of the output raster

		if (needCoordinateTransformation)
			crd = _outMergeRaster->coordinateSystem()->coord2coord(_inPointMap->coordinateSystem(), crd);
		// a raster contains pixels, each pixel represents a location(world coordinate) in the coordinate system of the raster
		// the georefence translate (both directions) between pixel locations and world locations
		// now we move the pixel iterator to the right location. normally we use operators like ++, --. +=, -= for this job
		// but as a pointmap is representing a very sparsely filled raster we can use a slightly less performant operator (assignment)
		// to keep the code simple. There are not that many pixels to be filled
		if (crd.isValid())
		{
			Pixel pix = _outMergeRaster->georeference()->coord2Pixel(crd);
			pix.z = 0;
			if (!pix.isValid() || IsEdgeCell(pix) || !m_iterOut.contains(pix))
				continue;
			m_iterOut = pix;
			// and we assign a value to the pixel
			quint64 id = feature->featureid();
			*m_iterOut = id + 1;
			ol.pxl = pix;
			ol.StreamID = *iterInDrng(pix);

			if (ol.StreamID == iUNDEF || ol.StreamID <= 0)
			{
				if (fRelocatOutlet(ol.pxl, 1) == false)  // based on 3 * 3 window
					if (fRelocatOutlet(ol.pxl, 2) == false) // based on 5 * 5 window
						continue;
				ol.StreamID = *iterInDrng(ol.pxl);
			}
			ol.StreveOrder = colStreve[ol.StreamID - 1].toInt();

			Pixel pxlFrom = m_vDrainageAtt[ol.StreamID - 1].UpstreamCoord;

			if (pxlFrom.y != ol.pxl.y || pxlFrom.x != ol.pxl.x)
			{
				ol.rLen1 = CalculateLength(pxlFrom, ol.pxl, ol.StreamID);
				ol.rLen2 = m_vDrainageAtt[ol.StreamID - 1].rLenght - ol.rLen1;
				ol.isOnNode = false;
			}
			else
			{
				ol.rLen1 = rUNDEF;
				ol.rLen2 = rUNDEF;
				ol.isOnNode = true;
			}
			ol.fExtractOverlandFlowPath = true;
			//ol.StreveOrder = iUNDEF;  //outlet is located on the junction! 
			m_vOutlet.push_back(ol);
		}
	}

}


void MapCatchmentMerge::InitJointOutletVector()
{
	PixelIterator iterInDrng = PixelIterator(_inDrngOrderRaster, BoundingBox(), PixelIterator::fXYZ);

	ITable tblAtt = _inDrngOrderRaster->attributeTable();

	std::vector<QVariant> colStreve;
	colStreve = tblAtt->column("Shreve");
	OutletLocation ol;

	// we create a pixeliterator to move around on the output raster coverage
	PixelIterator pixiter(_outMergeRaster);
	// create an empty raster (all pixels nodata)
	PixelIterator end = _outMergeRaster->end();
	while (pixiter != end) {
		*pixiter = iUNDEF;
		++pixiter;
	}

	m_iterOut = PixelIterator(_outMergeRaster, BoundingBox(), PixelIterator::fXYZ);

	for (std::vector<OutletLocation>::iterator pos = m_vOutlet.begin(); pos < m_vOutlet.end(); ++pos)
	{
		Pixel pxl = pos->pxl;
		
		pos->StreamID = *iterInDrng(pxl);

		if (pos->StreamID == iUNDEF || pos->StreamID <= 0)
		{
			if (fRelocatOutlet(ol.pxl, 1) == false)  // based on 3 * 3 window
				if (fRelocatOutlet(ol.pxl, 2) == false) // based on 5 * 5 window
					continue;
			pos->StreamID = *iterInDrng(ol.pxl);
		}
		pos->StreveOrder = colStreve[pos->StreamID - 1].toInt();

		Pixel pxlFrom = m_vDrainageAtt[pos->StreamID - 1].UpstreamCoord;

		if (pxlFrom.y != pos->pxl.y || pxlFrom.x != pos->pxl.x)
		{
			pos->rLen1 = CalculateLength(pxlFrom, pos->pxl, pos->StreamID);
			pos->rLen2 = m_vDrainageAtt[pos->StreamID - 1].rLenght - pos->rLen1;
			pos->isOnNode = false;
		}
		else
		{
			pos->rLen1 = rUNDEF;
			pos->rLen2 = rUNDEF;
			pos->isOnNode = true;
		}
		pos->fExtractOverlandFlowPath = true;
	}
}

double MapCatchmentMerge::CalculateLength(Pixel pxl1, Pixel pxl2, long iDrainage)
{
	double rLength = 0;
	Pixel rcDownCell = pxl1;
	int iFlow = 1; //anyway should do onece
	bool fCalculate = true;

	PixelIterator iterInDrng = PixelIterator(_inDrngOrderRaster, BoundingBox(), PixelIterator::fXYZ);
	PixelIterator iterFld = PixelIterator(_inFldRaster, BoundingBox(), PixelIterator::fXYZ);

	Coordinate cd = _outMergeRaster->georeference()->pixel2Coord(pxl1);
	m_vStreamCoord.push_back(cd);

	while ((iFlow != 0) && fCalculate)
	{
		fCalculate = ((rcDownCell != pxl2) && (IsEdgeCell(rcDownCell) != true));
		iFlow = GetDownStreamCell(rcDownCell);
		if ((*iterInDrng(rcDownCell) == iDrainage) && (*iterFld(rcDownCell) != 0))   //valid flow
		{
			Coordinate c1 = _outMergeRaster->georeference()->pixel2Coord(pxl1);
			Coordinate c2 = _outMergeRaster->georeference()->pixel2Coord(rcDownCell);
			rLength = rLength + rDistance(c1, c2);
			pxl1 = rcDownCell;
			m_vStreamCoord.push_back(c2);
		}
	}
	return rLength;
}


double MapCatchmentMerge::rDistance(Coordinate cd1, Coordinate cd2)
{
	double rDist;
	if (fLatLonCoords())
	{
		double rRadi = rDefaultEarthRadius;
		IConventionalCoordinateSystem projectedCsy = _inDrngOrderRaster->coordinateSystem().as<ConventionalCoordinateSystem>();
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


bool MapCatchmentMerge::fLatLonCoords()
{
	return _inDrngOrderRaster->coordinateSystem()->isLatLon();
}


int MapCatchmentMerge::GetDownStreamCell(Pixel& pxl)
{
	//Return a cell that the current given cell rc flows to,
	//otherwise, return the same cell as the given cell, This means
	//that the given cell doesn't flow to any other cell.   

	PixelIterator iterFld = PixelIterator(_inFldRaster, BoundingBox(), PixelIterator::fXYZ);

	int iPos = *iterFld(pxl);
	switch (iPos)
	{
	case 1: 	//East
		pxl.x = pxl.x + 1;
		break;
	case 2:  //South East 
		pxl.y = pxl.y + 1;
		pxl.x = pxl.x + 1;
		break;
	case 3: 	//South
		pxl.y = pxl.y + 1;
		break;
	case 4: //South West
		pxl.y = pxl.y + 1;
		pxl.x = pxl.x - 1;
		break;
	case 5:	//West
		pxl.x = pxl.x - 1;
		break;
	case 6:	//North West 
		pxl.y = pxl.y - 1;
		pxl.x = pxl.x - 1;
		break;
	case 7:	//North
		pxl.y = pxl.y - 1;
		break;
	case 8:	//North East
		pxl.y = pxl.y - 1;
		pxl.x = pxl.x + 1;
		break;
	default:
		iPos = 0;
	}
	return iPos;
}


//if the outlet is off the drainage, relocate it to the neast drainage line
//in defined iSize by iSize pixel window
bool MapCatchmentMerge::fRelocatOutlet(Pixel& rc, int iSize)
{
	PixelIterator iterInDrng = PixelIterator(_inDrngOrderRaster, BoundingBox(), PixelIterator::fXYZ);

	for (int i = -iSize; i <= iSize; ++i)
	{
		for (int j = -iSize; j <= iSize; ++j)
		{
			Pixel pxl;
			pxl.x = rc.x + j;
			pxl.y = rc.y + i;
			pxl.z = 0;

			if (*iterInDrng(pxl) > 0)
			{
				rc.y = rc.y + i;
				rc.x = rc.x + j;
				return true;
			}
		}
	}
	return false;
}


bool MapCatchmentMerge::IsEdgeCell(Pixel pxl)
{
	if (pxl.y == 0 || pxl.y == _ysize - 1 ||
		pxl.x == 0 || pxl.x == _xsize - 1)
		return true;
	else
		return false;
}

void MapCatchmentMerge::AddLink2StreamSegments()
{
	//Add a catchment link for drainage attribute
	ITable tblAtt = _inDrngOrderRaster->attributeTable();

	std::vector<QVariant> colOrder;
	if (m_sOrderSystem == "strahler")
		colOrder = tblAtt->column("Strahler");
	else
		colOrder = tblAtt->column("Shreve");

	std::vector<QVariant> colDownStreamLink = tblAtt->column("DownstreamLinkID");
	std::vector<QVariant> colUpstreamLink = tblAtt->column("UpstreamLinkID");
	std::vector<QVariant> colDownstreamCoord = tblAtt->column("DownstreamCoord");

	std::vector<QVariant> colTostreamCoord = tblAtt->column("TostreamCoord");
	std::vector<QVariant> colUpstreamCoord = tblAtt->column("UpstreamCoord");

	std::vector<QVariant> colLength = tblAtt->column("Length");

	long iSize = tblAtt->recordCount();

	m_iMaxOrderNumber = 0;
	DrainageAtt da;

	for (long i = 0; i < iSize; i++)
	{
		da.ID = i + 1;
		if (da.ID == iUNDEF)
			continue;
		da.UpstreamID.clear();
		SplitString(colUpstreamLink[i].toString(), QString(","), da.UpstreamID);
		da.iOrder = colOrder[i].toInt();
		if (da.iOrder > m_iMaxOrderNumber)
			m_iMaxOrderNumber = da.iOrder;

		da.iDownStreamID = colDownStreamLink[i].toInt();
		if (colDownstreamCoord[i].isNull())
			continue;

		da.DownStreamCoord = _inDrngOrderRaster->georeference()->coord2Pixel(colDownstreamCoord[i].value<Coordinate>());
		da.DownStreamCoord.z = 0;
		da.TostreamCoord = _inDrngOrderRaster->georeference()->coord2Pixel(colTostreamCoord[i].value<Coordinate>());
		da.TostreamCoord.z = 0;
		da.UpstreamCoord = _inDrngOrderRaster->georeference()->coord2Pixel(colUpstreamCoord[i].value<Coordinate>());
		da.UpstreamCoord.z = 0;
		da.CatchmentLink = iUNDEF;
		da.rLenght = colLength[i].toDouble();
		m_vDrainageAtt.push_back(da);
	}

}


Pixel MapCatchmentMerge::CoordinateStringToPixel(QString coordStr)
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
	return _inDrngOrderRaster->georeference()->coord2Pixel(crd);
}

void MapCatchmentMerge::InitInOutMaps()
{

}

REGISTER_OPERATION(MapCatchmentMergeWithOutlet)

MapCatchmentMergeWithOutlet::MapCatchmentMergeWithOutlet(quint64 metaid, const Ilwis::OperationExpression& expr) : MapCatchmentMerge(metaid, expr)
{}

Ilwis::OperationImplementation* MapCatchmentMergeWithOutlet::create(quint64 metaid, const Ilwis::OperationExpression& expr)
{
	return new MapCatchmentMergeWithOutlet(metaid, expr);
}


Ilwis::OperationImplementation::State MapCatchmentMergeWithOutlet::prepare(ExecutionContext* ctx, const SymbolTable& st)
{
	if (MapCatchmentMerge::prepare(ctx, st) == sPREPAREFAILED)
		return sPREPAREFAILED;

	QString outletStr = _expression.parm(4).value();
	QString includeStr = _expression.parm(5).value().toLower();

	m_includeUndefine = (includeStr == "yes" || includeStr == '1') ? true : false;

	if (!_inPointMap.prepare(outletStr, itFEATURE)) {
		ERROR2(ERR_COULD_NOT_LOAD_2, outletStr, "");
		return sPREPAREFAILED;
	}

	m_UseOutlets = true;

	QString outputName = _expression.parm(0, false).value();

	// prepare extracted segment feature
	_outputExtraxtedSegmentMap.prepare(QString(INTERNAL_CATALOG + "/%1").arg(outputName));
	_outputExtraxtedSegmentMap->coordinateSystem(_inDrngOrderRaster->georeference()->coordinateSystem());
	_outputExtraxtedSegmentMap->envelope(_inDrngOrderRaster->georeference()->envelope());

	return sPREPARED;
}

quint64 MapCatchmentMergeWithOutlet::createMetadata()
{
	OperationResource operation({ "ilwis://operations/MapCatchmentMergeWithOutlet" });
	operation.setSyntax(
		"MapCatchmentMerge(DrainageNetworkOrderMap,FlowDiractionMap,FlowaccumulationMap,DEM,OutletPointMap,IncludeUndefinedPixels)");

	operation.setDescription(TR("New merged catchments will be created based on the stream outlets within a catchment; all adjacent catchments that drain into such outlets will be merged"));
	operation.setInParameterCount({ 6 });
	operation.addInParameter(0, itRASTER, TR("Drainage NetWork Ordering Map"), TR("input raster that is the output of the Drainage network ordering operation"));
	operation.addInParameter(1, itRASTER, TR("Flow Direction Map"), TR("input raster that is the output of the Flow direction operation"));
	operation.addInParameter(2, itRASTER, TR("Flow accumulation Map"), TR("input raster that is the output of the Flow accumulation operation"));
	operation.addInParameter(3, itRASTER, TR("DEM"), TR("input raster that is supposed to be a Digital Elevation Model (DEM)"));
	operation.addInParameter(4, itPOINT, TR("Outlet location Map"), TR("input point map that contains the outlet locations; all sub-catchments draining to these outlets will be merged into new catchments."), OperationResource::ueNONE);
	operation.addInParameter(5, itSTRING, TR("Includ undefined pixels=yes|no"), TR("This option should be used when your study area contains lakes"), OperationResource::ueNONE);

	operation.parameterNeedsQuotes(1);

	operation.setOutParameterCount({ 4 });

	operation.addOutParameter(0, itRASTER, TR("Output Raster Map"), TR("output raster map that will contain the merged catchments"));
	operation.addOutParameter(1, itPOLYGON, TR("output polygon map"), TR("output polygons of the extracted catchments"));
	operation.addOutParameter(2, itLINE, TR("Output Segment Map"), TR("output segment map that will contain the longest flow path in each merged catchment."));
	operation.addOutParameter(3, itLINE, TR("Output segment Map"), TR("output segment that will contain only those segments, that fall within the new catchments."));

	operation.setKeywords("raster,table,segment,catchment, merge");

	return MapCatchmentMerge::createMetadata(operation);
}

REGISTER_OPERATION(MapCatchmentMergeWithStreamOrder)

MapCatchmentMergeWithStreamOrder::MapCatchmentMergeWithStreamOrder(quint64 metaid, const Ilwis::OperationExpression& expr) : MapCatchmentMerge(metaid, expr)
{}


Ilwis::OperationImplementation* MapCatchmentMergeWithStreamOrder::create(quint64 metaid, const Ilwis::OperationExpression& expr)
{
	return new MapCatchmentMergeWithStreamOrder(metaid, expr);
}


Ilwis::OperationImplementation::State MapCatchmentMergeWithStreamOrder::prepare(ExecutionContext* ctx, const SymbolTable& st)
{
	if (MapCatchmentMerge::prepare(ctx, st) == sPREPAREFAILED)
		return sPREPAREFAILED;

	QString streamordervalueStr = _expression.parm(5).value();
	m_iStreamOrders = streamordervalueStr.toInt();

	//QString extractOrderStr = _expression.parm(6).value().toLower();
	//m_useExtraOrder = (extractOrderStr == "yes" || extractOrderStr == '1') ? true : false;

	m_useExtraOrder = false;
	m_sOrderSystem = _expression.parm(4).value().toLower();

	QString outputName = _expression.parm(0, false).value();
	// prepare extracted segment feature
	_outputExtraxtedSegmentMap.prepare(QString(INTERNAL_CATALOG + "/%1").arg(outputName));
	_outputExtraxtedSegmentMap->coordinateSystem(_inDrngOrderRaster->georeference()->coordinateSystem());
	_outputExtraxtedSegmentMap->envelope(_inDrngOrderRaster->georeference()->envelope());

	m_UseOutlets = false;

	return sPREPARED;

}


quint64 MapCatchmentMergeWithStreamOrder::createMetadata()
{
	OperationResource operation({ "ilwis://operations/MapCatchmentMergeWithStreamOrder" });

	operation.setSyntax(
		"MapCatchmentMerge(DrainageNetworkOrderMap,FlowDiractionMap,FlowaccumulationMap,DEM,StreamOrderSystem,StreamOrderValue,ExtractOriginalOrder)");
	operation.setDescription(TR("New merged catchments will be created according to the specified stream order system and the order number; All contiguous catchments which drainages have the specified order number will be merged"));

	operation.setInParameterCount({ 6 });
	operation.addInParameter(0, itRASTER, TR("Drainage Net Work Ordering Map"), TR("input raster map that is the output of the Drainage Network Ordering operation"));
	operation.addInParameter(1, itRASTER, TR("Flow Direction Map"), TR("input raster map that is the output of the Flow Direction operation"));
	operation.addInParameter(2, itRASTER, TR("Flow accumulation Map"), TR("input raster that is the output of the Flow accumulation operation"));
	operation.addInParameter(3, itRASTER, TR("DEM"), TR("input raster that is supposed to be a Digital Elevation Model (DEM)"));
	operation.addInParameter(4, itSTRING, TR("Stream Ordering System = straher|shreve"), TR("Specify Straher or Shreve stream ordering method"));
	operation.addInParameter(5, itSTRING, TR("Stream Order number"), TR("Straher or Shreve stream order number"), OperationResource::ueNONE);

	operation.parameterNeedsQuotes(3);

	operation.setOutParameterCount({ 4 });

	operation.addOutParameter(0, itRASTER, TR("Output Raster Map"), TR("output raster map that will contain the merged catchments"));
	operation.addOutParameter(1, itPOLYGON, TR("output polygon map"), TR("output polygons of the extracted catchments"));
	operation.addOutParameter(2, itLINE, TR("Output Segment Map"), TR("output segment map that will contain the longest flow path in each merged catchment."));
	operation.addOutParameter(3, itLINE, TR("Output segment Map"), TR("output segment that will contain only those segments, that fall within the new catchments."));

	operation.setKeywords("raster,table,segment,catchment, merge");

	return MapCatchmentMerge::createMetadata(operation);

}

