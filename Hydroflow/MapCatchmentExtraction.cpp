
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
#include "geos/geom/CoordinateSequence.h"
#include "geos/geom/CoordinateSequenceFactory.h"
#include "geos/geom/GeometryFactory.h"
#include "geos/geom/LineString.h"
#include "coordinatedomain.h"
#include "ellipsoid.h"
#include "itemdomain.h"
#include "thematicitem.h"
#include "geos/geom/Point.h"
#include "geos/geom/Polygon.h"

//#include "raster2polygon.h"


#include "MapCatchmentExtraction.h"

using namespace Ilwis;
using namespace Hydroflow;

REGISTER_OPERATION(MapCatchmentExtraction)

const double rDefaultEarthRadius = 6371007.0;

using namespace std;
//
//static CoordSystem csyLamCyl(FileName fn)
//{
//    CoordSystem csy;
//    FileName fnCoordSyst(fn, ".csy", true);
//
//    CoordSystemProjection* cspr = new CoordSystemProjection(fnCoordSyst, 1);
//    csy.SetPointer(cspr);
//    cspr->datum = new MolodenskyDatum("WGS 1984", "");
//    cspr->ell = cspr->datum->ell;
//
//    Projection prj = cspr->prj;
//    String sPrj("Lambert Cylind EqualArea");
//    prj = Projection(sPrj, cspr->ell);
//    cspr->prj = prj;
//
//    return csy;
//}
//
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

//static void VerifyColumns(Map mpMap)
//{
//
//    Table tblAtt = mpMap->tblAtt();
//    if (!tblAtt.fValid())
//        throw ErrorNoAttTable(mpMap->fnObj);
//
//    if (!tblAtt[sUpstreamLink].fValid())
//        ColumnNotFoundError(tblAtt->fnObj, sUpstreamLink);
//
//    if (!tblAtt[sDownstreamLink].fValid())
//        ColumnNotFoundError(tblAtt->fnObj, sDownstreamLink);
//
//    if (!tblAtt[sDownstreamCoord].fValid())
//        ColumnNotFoundError(tblAtt->fnObj, String("TostreamCoord")); //sDownstreamCoord );
//}
//

MapCatchmentExtraction::MapCatchmentExtraction()
{
}


MapCatchmentExtraction::MapCatchmentExtraction(quint64 metaid, const Ilwis::OperationExpression& expr) : OperationImplementation(metaid, expr)
{

}

bool MapCatchmentExtraction::execute(ExecutionContext* ctx, SymbolTable& symTable)
{
    if (_prepState == sNOTPREPARED)
        if ((_prepState = prepare(ctx, symTable)) != sPREPARED)
            return false;

    bool resource =  executeCatchmentExtraction();

	if (resource && ctx != 0) 
	{
		if (_outCatchmentRaster.isValid()) {
			_outCatchmentRaster->setAttributes(_outputTable);
			std::string connectivity = "8";
			std::string output_name = "polygon_object_" + QString::number(Ilwis::Identity::newAnonymousId()).toStdString();
			QString expr = QString::fromStdString(output_name + "=raster2polygon(" + _outCatchmentRaster->name().toStdString() + "," + connectivity + ",true)");
			Ilwis::commandhandler()->execute(expr, ctx, symTable);
			Ilwis::Symbol result = symTable.getSymbol(ctx->_results[0]);
			if (result._type == itFEATURE && result._var.canConvert<Ilwis::IFeatureCoverage>()) {

				_outputfeatures = Ilwis::IFeatureCoverage(result._var.value<Ilwis::IFeatureCoverage>());
				quint32 count = _outputfeatures->featureCount();
				quint32 recordcount = _outputTable->recordCount();
				for (int rec = 0; rec < count; ++rec) {
					SPFeatureI& feature = _outputfeatures->feature(rec);
					quint64 id = feature->featureid();
					geos::geom::Polygon* polygon = dynamic_cast<geos::geom::Polygon*>(feature->geometry().get());
					if (polygon) {
						double area = polygon->getArea();
						double length = polygon->getLength();
						geos::geom::Point* centroid = polygon->getInteriorPoint();
						Coordinate crd;
						crd.x = centroid->getX();
						crd.y = centroid->getY();
						crd.z = 0;
						QString crdstr = CoordinateFormatString(crd.toString());
						_outputTable->setCell("CenterCatchment", rec, QVariant(crdstr));
						_outputTable->setCell("Perimeter", rec, length);
						_outputTable->setCell("CatchmentArea", rec, area);

						AttUpstreamLink vULs(m_vvUpstreamLinks[rec]);

						if (vULs.UpstreamLink.size() == 1 && vULs.UpstreamLink[0] == 0)
						{
							_outputTable->setCell("TotalUpstreamArea", rec, area);
							_outputfeatures->attributeTable()->setCell("TotalUpstreamArea", rec, area);
						}
						else
						{
							_outputTable->setCell("TotalUpstreamArea", rec, rUNDEF);
							_outputfeatures->attributeTable()->setCell("TotalUpstreamArea", rec, rUNDEF);
						}

						_outputfeatures->attributeTable()->setCell("CenterCatchment", rec, QVariant(crdstr));
						_outputfeatures->attributeTable()->setCell("Perimeter", rec, length);
						_outputfeatures->attributeTable()->setCell("CatchmentArea", rec, area);

					}

					/*		Ilwis::IFeatureCoverage polygons(result._var.value<Ilwis::IFeatureCoverage>());
							quint32 count = polygons->featureCount();
							quint32 recordcount = _outputTable->recordCount();
							for (int rec = 0; rec < count; ++rec) {
								SPFeatureI& feature = polygons->feature(rec);
								quint64 id = feature->featureid();
								geos::geom::Polygon* polygon = dynamic_cast<geos::geom::Polygon*>(feature->geometry().get());
								if (polygon) {
									double area = polygon->getArea();
									double length = polygon->getLength();
									geos::geom::Point * centroid = polygon->getInteriorPoint();
									Coordinate crd;
									crd.x = centroid->getX();
									crd.y = centroid->getY();
									crd.z = 0;
									QString crdstr = CoordinateFormatString(crd.toString());
									_outputTable->setCell("CenterCatchment", rec, QVariant(crdstr));
									_outputTable->setCell("Perimeter", rec, length);
									_outputTable->setCell("CatchmentArea", rec, area);

									AttUpstreamLink vULs(m_vvUpstreamLinks[rec]);

									if (vULs.UpstreamLink.size() == 1 && vULs.UpstreamLink[0] == 0)
										_outputTable->setCell("TotalUpstreamArea", rec, area);
									else
										_outputTable->setCell("TotalUpstreamArea", rec, rUNDEF);
								}*/
				}
			}

			if (_outCatchmentRaster.isValid())
			{
				QVariant outraster;
				outraster.setValue<IRasterCoverage>(_outCatchmentRaster);
				logOperation(_outCatchmentRaster, _expression, { _inDrngOrderRaster, _inFldRaster });
				ctx->addOutput(symTable, outraster, _outCatchmentRaster->name(), itRASTER, _outCatchmentRaster->resource());
			}

			if (_outputfeatures.isValid())
			{
				_outputfeatures->setAttributes(_outputTable);
				QVariant value;
				value.setValue<IFeatureCoverage>(_outputfeatures);
				logOperation(_outputfeatures, _expression, { _inDrngOrderRaster });
				ctx->addOutput(symTable, value, _outputfeatures->name(), itFEATURE, _outputfeatures->resource());
			}
		}
		m_vvUpstreamLinks.resize(0);
	}

	return resource;
}

bool MapCatchmentExtraction::executeCatchmentExtraction()
{
    bool ret(true);

    quint64 copylist = itRASTERSIZE | itENVELOPE | itINTEGER;
    _flagRaster = OperationHelperRaster::initialize(_inDrngOrderRaster.as<IlwisObject>(), itRASTER, copylist);

    PixelIterator iterFlag = PixelIterator(_flagRaster, BoundingBox(), PixelIterator::fXYZ);
    PixelIterator iterInDrng = PixelIterator(_inDrngOrderRaster, BoundingBox(), PixelIterator::fXYZ);
    PixelIterator iterFld = PixelIterator(_inFldRaster, BoundingBox(), PixelIterator::fXYZ);


	NumericStatistics stats;
	stats.calculate(iterInDrng, iterInDrng.end(), NumericStatistics::pBASIC);
	double maxDrainage = stats[NumericStatistics::pMAX];
	iterInDrng = PixelIterator(_inDrngOrderRaster, BoundingBox(), PixelIterator::fXYZ); //rewind

    PixelIterator iterOut = PixelIterator(_outCatchmentRaster, BoundingBox(), PixelIterator::fXYZ);
    PixelIterator inEnd = iterInDrng.end();

    while (iterInDrng != inEnd)
    {
        Pixel pxl = iterInDrng.position();
        pxl.z = 0;

		*iterFlag[pxl] = iUNDEF;
		if (IsEdgeCell(pxl) || *iterFld(pxl) == rUNDEF)
			*iterFld(pxl) = 0;

		iterInDrng++;
    }

    GetAttributes();

	while (m_vRecords.size() != 0)
	{
		for (unsigned long index = 0; index < m_vRecords.size(); ++index)
		{
			AttCols attFields(m_vRecords[index]);

			//should start with a no upstream channel 
			//or with a channel that upstream(s) have been processed.  
			if (attFields.UpstreamLink.size() == 1 && attFields.UpstreamLink[0] == 0)
			{
				//Extract catchment, and flag the cells with DrainageID that flow to Downstream Coord 
				DelineateCatchment(attFields.DownstreamCoord, attFields.DrainageID);

				//Remove the upstreamlink for the downstream that it drains to	
				UpdateUpstreamLinkID(attFields.DownstreamLink, attFields.DrainageID);

				EraseDrainage(attFields.DrainageID);
			}
		}
	}

	iterFlag = PixelIterator(_flagRaster, BoundingBox(), PixelIterator::fXYZ);
	iterOut = PixelIterator(_outCatchmentRaster, BoundingBox(), PixelIterator::fXYZ);

	while (iterOut != inEnd)
	{
		*iterOut = (int)*iterFlag;
		*iterOut = (*iterOut != 0) ? (*iterOut - 1) : rUNDEF;
		*iterFlag++;
		*iterOut++;
	}
	m_vRecords.resize(0);
	SetAttributeTable();



	//bool canConnect;
	//NetworkPoint previousPoint;
	//std::vector<PrevLinePoint> previousLine(_outCatchmentRaster->size().xsize());
	//BlockIterator iterBlock(_outCatchmentRaster, Size<>(3, 3, 1), BoundingBox(), { 1,1,1 }, true);
	//BlockIterator end = iterBlock.end();
	//while (iterBlock != end) {
	//	auto block = (*iterBlock).to3DVector();
	//	if (block.size() == 1)
	//		continue;
	//	GridBlock& blockOut = *iterBlock;

	//	auto currentValue = block[1][1][0];
	//	auto pCenter = (*iterBlock).position();
	//	if (currentValue != block[0][1][0] && currentValue != block[0][0][0]) {
	//		// temporary points are points that do not end up in the network but are in _previousLine. As they represent a
	//		// intermediate point on a straight vertical line segment they have no place in the final network but are needed
	//		// to give the _previousline the correct linkage.
	//		bool  isTempPoint = block[0][0][0] != currentValue &&
	//			block[0][0][0] == block[0][2][0] &&
	//			currentValue == block[0][1][0];
	//		if (canConnect) {
	//			previousPoint = makeConnection(pCenter, isTempPoint, previousLine);
	//		}
	//		canConnect = true;
	//	}
	//	// if the topmid value equals the current value we are looking at a point in the middle of a filled area. so no lines we be drawn here
	//	if (block[1][1][0] == currentValue) {
	//		canConnect = false;
	//		previousPoint = NetworkPoint();
	//	}
	//}




	//Check the type of coordinate system  
	//if (fLatLonCoords())
	//{
	//	//Transform the map to Lambert Cylind EqualArea projection coordinate system needed to be able
	//	//to calculate the perimeter and area parameters
	//	//LAMCYL is a pre-defined coordinate system, should be placed in the ILWIS system directory
	//	FileName fnTmpTFPol(fnObj, ".mpa");
	//	fnTmpTFPol = FileName::fnUnique(fnTmpTFPol);
	//	//String sExprPMT("PolygonMapTransform(%S, %S, %g)", fnPol.sFullPathQuoted(false), String("lamcyl"), 0.000000); 
	//	CoordSystem csy = csyLamCyl(fnTmpTFPol);
	//	csy->Store(); // explicit Store(), due to different behavior between Debug and Release build!! (csyLamCyl will auto-store when the internal csy object is destructed (as it is supposed to do) in the debug build, but not in the release build)
	//	csy->fErase = true;
	//	String sExprPMT("PolygonMapTransform(%S, %S, %g)", fnPol.sFullPathQuoted(false), csy->sName(), 0.000000);
	//	PolygonMap polTmpTFMap;
	//	polTmpTFMap = PolygonMap(fnTmpTFPol, sExprPMT);
	//	polTmpTFMap->Calc();
	//	polTmpTFMap->fErase = true;
	//	//***a temporary histogram file needed to retrieve attributes about Area and Perimeter   
	//	fnTmpHsa = FileName(fnTmpTFPol, ".hsa");
	//	sExprTbl = String("TableHistogramPol(%S)", fnTmpTFPol.sFullPathQuoted());
	//	tblHsa = Table(fnTmpHsa, sExprTbl);
	//}
	//else
	//{
	//	fnTmpHsa = FileName(fnObj, ".hsa");
	//	sExprTbl = String("TableHistogramPol(%S)", fnPol.sFullPathQuoted());
	//	tblHsa = Table(fnTmpHsa, sExprTbl);
	//}

	ComputeCatchmentAttributes();
	ComputeCenterDrainage();


    return ret;
}


Ilwis::OperationImplementation* MapCatchmentExtraction::create(quint64 metaid, const Ilwis::OperationExpression& expr)
{
    return new MapCatchmentExtraction(metaid, expr);
}

Ilwis::OperationImplementation::State MapCatchmentExtraction::prepare(ExecutionContext* ctx, const SymbolTable& st)
{
    OperationImplementation::prepare(ctx, st);
    QString drnetwkraster = _expression.parm(0).value();
	QString outputName = _expression.parm(0, false).value();

    QString flowraster = _expression.parm(1).value();
    QString outcatchment = _expression.parm(2).value();

    if (!_inDrngOrderRaster.prepare(drnetwkraster, itRASTER)) {
        ERROR2(ERR_COULD_NOT_LOAD_2, drnetwkraster, "");
        return sPREPAREFAILED;
    }

    if (!_inFldRaster.prepare(flowraster, itRASTER)) {
        ERROR2(ERR_COULD_NOT_LOAD_2, flowraster, "");
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

	int copylist = itRASTERSIZE | itENVELOPE | itCOORDSYSTEM | itGEOREF;
    _outCatchmentRaster = OperationHelperRaster::initialize(_inDrngOrderRaster.as<IlwisObject>(), itRASTER, copylist);
    if (!_outCatchmentRaster.isValid()) {
        ERROR1(ERR_NO_INITIALIZED_1, "output rastercoverage");
        return sPREPAREFAILED;
    }

    _xsize = _inDrngOrderRaster->size().xsize();
    _ysize = _inDrngOrderRaster->size().ysize();

  
	if (_outCatchmentRaster.isValid())
	{
		QString catchName = "CatchID";
		_outDomain.prepare();
		_idrange = new NamedIdentifierRange();
		_outDomain->range(_idrange);

		DataDefinition def(_outDomain);
		_outCatchmentRaster->datadefRef() = def;
		for (int band = 0; band < _outCatchmentRaster->size().zsize(); ++band) {
			_outCatchmentRaster->datadefRef(band) = def;
		}
		return sPREPARED;
	}

	_outputfeatures.prepare(QString(INTERNAL_CATALOG + "/%1").arg(outputName));

	_inputgrf = _inDrngOrderRaster->georeference();
	_csy = _inputgrf->coordinateSystem();
	_outputfeatures->coordinateSystem(_csy);
	Envelope env = _inDrngOrderRaster->georeference()->envelope();
	_outputfeatures->envelope(env);

    return sPREPARED;
}

quint64 MapCatchmentExtraction::createMetadata()
{
    OperationResource operation({ "ilwis://operations/MapCatchmentExtraction" });
    operation.setSyntax("MapCatchmentExtraction(DrainageNetworkOrderingMap,FlowDiractionMap)");
    operation.setDescription(TR("Constructs catchments;a catchment will be calculated for each stream found in the Drainage Network Ordering map"));
    operation.setInParameterCount({ 2 });
    operation.addInParameter(0, itRASTER, TR("Drainage NetWork Ordering Map"), TR("input raster that is the output of a previous Drainage Network Ordering operation"));
    operation.addInParameter(1, itRASTER, TR("Flow Direction Map"), TR("input raster that is the output of a previous Flow Direction operation"));
    operation.parameterNeedsQuotes(1);

	operation.setOutParameterCount({ 1 });

	operation.addOutParameter(0, itRASTER, TR("output raster"), TR("output raster with the results of the catchment extraction"));
	operation.addOutParameter(1, itPOLYGON, TR("output polygon"), TR("output polygon with the results of the catchment extraction"));
	//operation.addOutParameter(2, itTABLE, TR("output table"), TR("output table with the results of the catchment extraction"));

	operation.setKeywords("raster,table,polygon,catchment, extraction");

    mastercatalog()->addItems({ operation });
    return operation.id();
}



void MapCatchmentExtraction::SetAttributeTable()
{
	// prepare table

	IFlatTable newTable;
	newTable.prepare();
	newTable->addColumn(_outCatchmentRaster->primaryKey(), _outDomain);

	_outputTable = newTable;

}


bool MapCatchmentExtraction::IsEdgeCell(Pixel pxl)
{
    if (pxl.y == 0 || pxl.y == _ysize - 1 ||
        pxl.x == 0 || pxl.x == _xsize - 1)
        return true;
    else
        return false;
}


bool MapCatchmentExtraction::fLatLonCoords()
{
	return _inDrngOrderRaster->coordinateSystem()->isLatLon();
}

//void MapCatchmentExtraction::CompitableGeorefs(FileName fn, Map mp1, Map mp2)
//{
//	bool fIncompGeoRef = false;
//	if (mp1->gr()->fGeoRefNone() && mp2->gr()->fGeoRefNone())
//		fIncompGeoRef = mp1->rcSize() != mp2->rcSize();
//	else
//		fIncompGeoRef = mp1->gr() != mp2->gr();
//	if (fIncompGeoRef)
//		throw ErrorIncompatibleGeorefs(mp1->gr()->sName(true, fn.sPath()),
//			mp2->gr()->sName(true, fn.sPath()), fn, errMapCatchmentExtraction);
//
//}

void MapCatchmentExtraction::GetAttributes()
{
	//Retrieve the network link attributes needed to be able 
	//to extract the catchment area

	ITable tblAtt = _inDrngOrderRaster->attributeTable();
	std::vector<QVariant> colUpstreamLink = tblAtt->column(sUpstreamLink);
	std::vector<QVariant> colDownstreamLink = tblAtt->column(sDownstreamLink);
	std::vector<QVariant> colDownstreamCoord = tblAtt->column(QString("TostreamCoord"));

	std::vector<QVariant> colPrimaryKey = tblAtt->column(_inDrngOrderRaster->primaryKey());

	long iSize = tblAtt->recordCount();
	AttCols ac;
	for (long i = 0; i < iSize; i++) 
	{
		QVariant val = colPrimaryKey[i];
		ac.DrainageID = val.toInt()+1;
		if (ac.DrainageID == iUNDEF )
			continue;
		ac.DownstreamLink = colDownstreamLink[i].toInt();
		if (ac.DownstreamLink == 0)
			ac.DownstreamLink = iUNDEF;

		QString coordsStr = colDownstreamCoord[i].toString();

		if (coordsStr.isEmpty())
			continue;

		coordsStr.replace("{", "");
		coordsStr.replace("}", "");

		QStringList coodrs = coordsStr.split(" ");

		QString xstr = coodrs[0];
		QString ystr = coodrs[1];

		Coordinate crd;
		crd.x = xstr.toDouble();
		crd.y = ystr.toDouble();
		crd.z = 0;
		ac.DownstreamCoord = _inDrngOrderRaster->georeference()->coord2Pixel(crd);

		SplitString(colUpstreamLink[i].toString(), QString(","), ac.UpstreamLink);
		m_vRecords.push_back(ac);
		m_vDrnIDs.push_back(ac.DrainageID);
	}

}


void MapCatchmentExtraction::SplitString(QString s, QString mid, std::vector<long>& results)
{
	s.replace("{", "");
	s.replace("}", "");
	QStringList strlst = s.split(mid);

	results.clear();
	for (unsigned int i = 0; i < strlst.size(); i++)
	{
		long res = strlst[i].toInt();
		if (res != iUNDEF)
			results.push_back(res);
	}
}

void MapCatchmentExtraction::ComputeCatchmentAttributes()
{
	////Retrieve attributes from drainage network attribute table needed 
	////for updating the catchment attributes  

	ITable tblAtt = _inDrngOrderRaster->attributeTable();

	std::vector<QVariant> colStreamID = tblAtt->column(_inDrngOrderRaster->primaryKey());

	std::vector<QVariant> colUpstreamLink = tblAtt->column(sUpstreamLink);
	IDomain colUplinkDomain = tblAtt->columndefinitionRef(sUpstreamLink).datadef().domain();

	if ( !colUplinkDomain.isValid() )
		throw ErrorObject(TR("Source map must have domain class or id"));

	std::vector<QVariant> colDownstreamLink = tblAtt->column(sDownstreamLink);
	std::vector<QVariant> colFlowLength = tblAtt->column(QString("Length"));

	long iSize = colStreamID.size();

	//Define catchment attributes
	_outputTable->addColumn("DrainageID", IlwisObject::create<IDomain>("value"),true);
	ColumnDefinition& coldef0 = _outputTable->columndefinitionRef("DrainageID");
	coldef0.datadef().range(new NumericRange(1, 32767, 1)); 

	_outputTable->addColumn("UpstreamLinkCatchment", IlwisObject::create<IDomain>("text"), true);

	_outputTable->addColumn("DownstreamLinkCatchment", IlwisObject::create<IDomain>("value"), true);
	ColumnDefinition& coldef1 = _outputTable->columndefinitionRef("DownstreamLinkCatchment");
	coldef1.datadef().range(new NumericRange(1, 32767, 1));

	_outputTable->addColumn("Perimeter", IlwisObject::create<IDomain>("value"));
	ColumnDefinition& coldef2 = _outputTable->columndefinitionRef("Perimeter");
	coldef2.datadef().range(new NumericRange(1, 1.0e300, 0.01));

	_outputTable->addColumn("CatchmentArea", IlwisObject::create<IDomain>("value"));
	ColumnDefinition& coldef3 = _outputTable->columndefinitionRef("CatchmentArea");
	coldef3.datadef().range(new NumericRange(1, 1.0e300, 0.01));

	_outputTable->addColumn("CenterCatchment", IlwisObject::create<IDomain>("text"), true);

	_outputTable->addColumn("TotalUpstreamArea", IlwisObject::create<IDomain>("value"));
	ColumnDefinition& coldef4 = _outputTable->columndefinitionRef("TotalUpstreamArea");
	coldef4.datadef().range(new NumericRange(1, 1.0e300, 0.01));

	_outputTable->addColumn("LongestFlowLength", IlwisObject::create<IDomain>("value"));
	ColumnDefinition& coldef5 = _outputTable->columndefinitionRef("LongestFlowLength");
	coldef5.datadef().range(new NumericRange(1, 1.0e300, 0.01));

	AttUpstreamLink vUpstreamLinks;
	for (long i = 0; i <iSize; i++)
	{
		//long iid = colUplinkDomain->da
		long iDrainageID = colStreamID[i].toInt()+1;
		vUpstreamLinks.DrainageID = iDrainageID;
		SplitString(colUpstreamLink[i].toString(), QString(","),vUpstreamLinks.UpstreamLink);
		m_vvUpstreamLinks.push_back(vUpstreamLinks);
	}

	//Retrieve upstream link attributes needed to compute the total area
	quint32 record = 0;
	quint32 totalcount = m_vDrnIDs.size();
	for (std::vector<long>::iterator pos = m_vDrnIDs.begin(); pos < m_vDrnIDs.end(); ++pos)
	{
		long rec(*pos);
		record = rec;

		QString id = QString::number(record);

		if (id != "")
		{
			if (_idrange->contains(id))
			{
				++totalcount;
				id = QString::number(totalcount);
			}
			*_idrange << id;
		}

		record = record - 1;
		long iDrainageID = colStreamID[record].toInt()+1;
		if (iDrainageID == iUNDEF)
			continue;
		_outputTable->setCell("DrainageID", record, QVariant(iDrainageID));
		_outputTable->setCell("UpstreamLinkCatchment", record, QVariant(colUpstreamLink[record].toString()));
		_outputTable->setCell("DownstreamLinkCatchment", record, QVariant(colDownstreamLink[record].toInt()));
		_outputTable->setCell("LongestFlowLength", record, QVariant(colFlowLength[record].toDouble()));
	
		////Initialize the totalarea  to the catchment area itsself, if it is a source link, it should be
		////Otherwise, assign an no data value, this will be computed later

	/*	AttUpstreamLink vULs(m_vvUpstreamLinks[i - 1]);

		if (vULs.UpstreamLink.size() == 1 && vULs.UpstreamLink[0] == 0)
			cTotalUpstreamArea->PutVal(iDrainageID, m_cArea->rValue(i));
		else
			cTotalUpstreamArea->PutVal(iDrainageID, rUNDEF);*/

		//cLongestFlowLength->PutVal(iDrainageID, colFlowLength->rValue(i));

	//}


		_outputTable->setCell(_outCatchmentRaster->primaryKey(), record, QVariant(record));
	}

//	m_vvUpstreamLinks.resize(0);
}



//void MapCatchmentExtraction::ComputeTotalUpstreamArea(DomainSort* pdsrt, Column cArea, Column cTotalArea)
//{
//	bool fComputeTotalArea = true;
//	long iSize = pdsrt->iSize();
//	trq.SetText(TR("Calculate the total upstream catchment area"));
//	while (fComputeTotalArea)
//	{
//		fComputeTotalArea = false;
//		for (long i = 1; i <= iSize; i++)
//		{
//			long iDrainageID = pdsrt->iOrd(i);
//			if (iDrainageID == iUNDEF)
//				continue;
//
//			if (cTotalArea->rValue(i) == rUNDEF)
//			{
//				fComputeTotalArea = true;
//				AttUpstreamLink vULs(m_vvUpstreamLinks[i - 1]);
//				vector<long> vLinks = vULs.UpstreamLink;
//				double rArea = 0;
//				for (vector<long>::iterator pos = vLinks.begin();
//					pos < vLinks.end(); ++pos)
//				{
//					//search ID in domain, return index, if the ID is found 
//					String sLbl("%li", (*pos));
//					long iRaw = pdsrt->iOrd(sLbl);
//					if (cTotalArea->rValue(iRaw) == rUNDEF)
//						break;
//					rArea += cTotalArea->rValue(iRaw);
//					cTotalArea->PutVal(iDrainageID, rArea);
//				}
//			}
//			if (trq.fUpdate(i, iSize)) return;
//		}
//	} //while()
//}

//void MapCatchmentExtraction::ComputerCenterPolygon(FileName fn)
//{
	////First lable point map of the polygon
	////then, put the point coordinates to the catchment attribute table
	//trq.SetText(TR("Compute the center of polygon"));
	//FileName fnTmpPoint(fnObj, ".mpp");
	//fnTmpPoint = FileName::fnUnique(fnTmpPoint);
	//String sExpr("PointMapPolLabels(%S)", fn.sFullPathQuoted(false));
	//PointMap ptTmpMap;
	//ptTmpMap = PointMap(fnTmpPoint, sExpr);
	//ptTmpMap->Calc();

	//Domain dmcrd;
	//dmcrd.SetPointer(new DomainCoord(ptTmpMap->cs()->fnObj));
	//Column cCenterCatchment = m_tbl->colNew(String("CenterCatchment"), dmcrd);

	//Coord crd;
	//long iPoint = ptTmpMap->iFeatures();
	//for (long i = 0; i < iPoint; ++i) {
	//	long iRaw = ptTmpMap->iRaw(i);
	//	if (iRaw == iUNDEF)
	//		continue;
	//	Coord crd = ptTmpMap->cValue(i);
	//	cCenterCatchment->PutVal(iRaw, crd);
	//	if (trq.fUpdate(i, iPoint)) return;
	//}
	//ptTmpMap->fErase = true;
//}

double MapCatchmentExtraction::GetDistance(Pixel& rc)
{
	double dist(0);
	//Return a cell that the current given cell rc flows to,
	//otherwise, return the same cell as the given cell, This means
	//that the given cell doesn't flow to any other cell.  

	PixelIterator iterFld = PixelIterator(_inFldRaster, BoundingBox(), PixelIterator::fXYZ);

	Pixel pospxl;
	pospxl.x = rc.x - 1;
	pospxl.y = rc.y - 1;
	pospxl.z = 0;

	int iPos = *iterFld[pospxl];
	Pixel rc2;
	rc2.z = 0;

	switch (iPos)
	{
	case 1: 	//East
		rc2.y = rc.y;
		rc2.x = rc.x + 1;
		break;
	case 2:  //South East 
		rc2.y = rc.y + 1;
		rc2.x = rc.x + 1;
		break;
	case 3: 	//South
		rc2.y = rc.y + 1;
		rc2.x = rc.x;
		break;
	case 4: //South West
		rc2.y = rc.y + 1;
		rc2.x = rc.x - 1;
		break;
	case 5:	//West
		rc2.y = rc.y;
		rc2.x = rc.x - 1;
		break;
	case 6:	//North West 
		rc2.y = rc.y - 1;
		rc2.x = rc.x - 1;
		break;
	case 7:	//North
		rc2.y = rc.y - 1;
		rc2.x = rc.x;
		break;
	case 8:	//North East
		rc2.y = rc.y - 1;
		rc2.x = rc.x + 1;
		break;
	default:
		rc2.y = rc.y;
		rc2.x = rc.x;
		break;
	}
	Coordinate c1 = _inDrngOrderRaster->georeference()->pixel2Coord(rc);
	Coordinate c2 = _inDrngOrderRaster->georeference()->pixel2Coord(rc2);

	dist = rDistance(c1, c2);
	rc.y = rc2.y;
	rc.x = rc2.x;
	return dist;
}

double MapCatchmentExtraction::rDistance(Coordinate cd1, Coordinate cd2)
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

void MapCatchmentExtraction::ComputeCenterDrainage()
{
	//Retrieve UpstreamCoord attributes from drainage network attribute table needed
  //to be able to trace the flow path and get the center of drainage.
	ITable tblAtt = _inDrngOrderRaster->attributeTable();
	std::vector<QVariant> colUpstreamCoord = tblAtt->column(QString("UpstreamCoord"));

	std::vector<QVariant> colUpstreamLink = tblAtt->column(sUpstreamLink);
	std::vector<QVariant> colLength = tblAtt->column(QString("Length"));

	std::vector<QVariant> colPrimaryKey = tblAtt->column(_inDrngOrderRaster->primaryKey());


	_outputTable->addColumn("CenterDrainage", IlwisObject::create<IDomain>("text"));

	long iSize = colPrimaryKey.size();
	for (long i = 0; i < iSize; i++)
	{
		long iDrainageID = colPrimaryKey[i].toInt();
		if (iDrainageID == iUNDEF )
			continue;

		QString coordsStr = colUpstreamCoord[i].toString();

		if (coordsStr.isEmpty())
			continue;

		coordsStr.replace("(", "");
		coordsStr.replace(")", "");

		coordsStr.replace("{", "");
		coordsStr.replace("}", "");

		Coordinate crd = Coordinate(coordsStr);

		crd.z = 0;
		Pixel pxl = _inDrngOrderRaster->georeference()->coord2Pixel(crd);
		double rLength = colLength[i].toDouble() / 2;
		double rDistance = 0;
		while (rDistance < rLength)
		{
			double rDist = GetDistance(pxl);
			if (rDist == 0)
				break;
			rDistance = rDistance + rDist;
		}

		pxl.x -= 1;
		pxl.y -= 1;

		crd = _inDrngOrderRaster->georeference()->pixel2Coord(pxl);

		QString crdstr = CoordinateFormatString(crd.toString());
		_outputTable->setCell("CenterDrainage", iDrainageID,QVariant(crdstr));

	}
}


QString MapCatchmentExtraction::CoordinateFormatString(QString crd)
{
	QString crdstring;
	crdstring = QString("{") + crd;
	crdstring = crdstring + QString("}");
	return crdstring;
}


long MapCatchmentExtraction::DelineateCatchment(Pixel pxl, long iFlag)
{
	/*For the specified downstream cell in loaction rc, 
	check whether its neighboring cells flow to it,
	If true, flag the cells with iFlag in m_vDrainage, then call the function recursively   
	The recursion stops when it reaches a cell that has no flow to it 
	location number
		-------
		|6|7|8|
		-------
		|5| |1|
		-------
		|4|3|2|
		-------*/
	 
	long iFlow = 1;
	bool isFlow; //determine if the neighboring cell flows to the cell in location rc 
	long in, jn;
	Pixel pospxl;
	pospxl.z = 0;

	PixelIterator iterFld = PixelIterator(_inFldRaster, BoundingBox(), PixelIterator::fXYZ);
	PixelIterator iterFlag = PixelIterator(_flagRaster, BoundingBox(), PixelIterator::fXYZ);


	for (int iNr = 1; iNr < 9; iNr++)
	{
		isFlow = false;
		switch (iNr)
		{
		case 1: 
		{	//East
			if (pxl.x != _xsize - 1)
			{
				pospxl.y = pxl.y;
				pospxl.x = pxl.x + 1;
				isFlow = (*iterFld(pospxl) == 5 && *iterFlag(pospxl) == iUNDEF);
			}
		}
			  break;
		case 2: 
		{ //South East 
			if (pxl.x != _xsize - 1 && pxl.y != _ysize - 1)
			{
				pospxl.y = pxl.y + 1;
				pospxl.x = pxl.x + 1;
				isFlow = (*iterFld(pospxl) == 6 && *iterFlag(pospxl) == iUNDEF);
			}
		}
			  break;
		case 3: {	//South
			if (pxl.y != _ysize - 1)
			{
				pospxl.y = pxl.y + 1;
				pospxl.x = pxl.x;
				isFlow = (*iterFld(pospxl) == 7 && *iterFlag(pospxl) == iUNDEF);
			}
		}
			  break;
		case 4: 
		{ //South West
			if ( pxl.x != 0 && pxl.y != _ysize - 1)
			{
				pospxl.y = pxl.y + 1;
				pospxl.x = pxl.x - 1;
				isFlow = (*iterFld(pospxl) == 8 && *iterFlag(pospxl) == iUNDEF);
			}
		}
			  break;
		case 5:
		{	//West
			if (pxl.x != 0)
			{
				pospxl.y = pxl.y;
				pospxl.x = pxl.x - 1;
				isFlow = (*iterFld(pospxl) == 1 && *iterFlag(pospxl) == iUNDEF);
			}
		}
			  break;
		case 6: 
		{	//North West 
			if (pxl.x != 0 && pxl.y != 0)
			{
				isFlow = false;
				pospxl.y = pxl.y - 1;
				pospxl.x = pxl.x - 1;
				isFlow = (*iterFld(pospxl) == 2 && *iterFlag(pospxl) == iUNDEF);
			}
		}
			  break;
		case 7: {	//North
			if (pxl.y != 0)
			{
				pospxl.y = pxl.y - 1;
				pospxl.x = pxl.x;
				isFlow = (*iterFld(pospxl) == 3 && *iterFlag(pospxl) == iUNDEF);
			}
		}
			  break;
		case 8: 
		{	//North East
			if (pxl.x != _xsize - 1 && pxl.y != 0)
			{
				pospxl.y = pxl.y - 1;
				pospxl.x = pxl.x + 1;
				isFlow = (*iterFld(pospxl) == 4 && *iterFlag(pospxl) == iUNDEF);
			}
		}
			  break;
		}
		if (isFlow)
		{
			iFlow += DelineateCatchment( pospxl, iFlag);
		}
		*iterFlag(pxl) = iFlag;
	}
	return iFlow;
}


long MapCatchmentExtraction::FindDownstreamIndex(long DowmstreamID)
{

	for (vector<AttCols>::iterator pos = m_vRecords.begin(); pos < m_vRecords.end(); ++pos)
	{
		AttCols ac(*pos);
		if (ac.DrainageID == DowmstreamID)
			return (long)(pos - m_vRecords.begin());
	}
	return -1;
}

void MapCatchmentExtraction::UpdateUpstreamLinkID(long DrainageID, long UpstreamID)
{

	if (DrainageID != iUNDEF )
	{
		//Find the downstream index in m_vRecords
		//returns position for the downstream in m_vRecords
		long iIndex = FindDownstreamIndex(DrainageID);

		if (iIndex > 0)
		{
			AttCols attFields(m_vRecords[iIndex]);

			if ((attFields.UpstreamLink.size() == 1) && (attFields.UpstreamLink[0] == UpstreamID))
				m_vRecords[iIndex].UpstreamLink[0] = 0;
			else if (attFields.UpstreamLink.size() > 1)
			{
				vector<long>::iterator pos;
				pos = find(attFields.UpstreamLink.begin(), attFields.UpstreamLink.end(), UpstreamID);
				if (pos != attFields.UpstreamLink.end())
				{
					attFields.UpstreamLink.erase(pos);
					m_vRecords[iIndex].UpstreamLink = attFields.UpstreamLink;
				}
			}
		}
	}
}

void MapCatchmentExtraction::EraseDrainage(long DrainageID)
{
	for (vector<AttCols>::iterator pos = m_vRecords.begin(); pos < m_vRecords.end(); ++pos)
	{
		AttCols ac(*pos);
		if (ac.DrainageID == DrainageID)
		{
			m_vRecords.erase(pos);
			break;
		}
	}
}


bool MapCatchmentExtraction::fEllipsoidalCoords()
{
	//CoordSystemViaLatLon* csviall = mp->cs()->pcsViaLatLon();
	//bool fSpheric = true;
	//if (csviall)
	//	fSpheric = (csviall->ell.fSpherical());
	//return (0 != csviall && 0 == fSpheric);
	return false;
}

//////////////////////
MapCatchmentExtraction::NetworkPoint MapCatchmentExtraction::makeConnection(const Pixel& pCenter, bool isTemp, std::vector<MapCatchmentExtraction::PrevLinePoint>& previousLine) {

	if (previousLine[pCenter.x].isValid()) {
		MapCatchmentExtraction::NetworkPoint currentPoint;
		MapCatchmentExtraction::PrevLinePoint newPPoint;
		currentPoint._x = pCenter.x;
		currentPoint._y = pCenter.y;
		auto topPoint = previousLine[pCenter.x];
		if (isTemp) {
			if (topPoint._shadowLink != iUNDEF) {
				newPPoint._shadowLink = topPoint._shadowLink;
			}
			else {
			}

		}
	}
	return MapCatchmentExtraction::NetworkPoint();
}