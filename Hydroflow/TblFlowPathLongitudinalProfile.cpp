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
#include "coordinatesystem.h"
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


#include "TblFlowPathLongitudinalProfile.h"


using namespace Ilwis;
using namespace Hydroflow;

REGISTER_OPERATION(TblFlowPathLongitudinalProfile)

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

TblFlowPathLongitudinalProfile::TblFlowPathLongitudinalProfile()
{
}


TblFlowPathLongitudinalProfile::TblFlowPathLongitudinalProfile(quint64 metaid, const Ilwis::OperationExpression& expr) : OperationImplementation(metaid, expr)
{

}

bool TblFlowPathLongitudinalProfile::execute(ExecutionContext* ctx, SymbolTable& symTable)
{
	if (_prepState == sNOTPREPARED)
		if ((_prepState = prepare(ctx, symTable)) != sPREPARED)
			return false;

	bool resource = executeTblFlowPathLongitudinalProfile();

	if (resource && ctx != 0) {
		QVariant value;
		value.setValue<ITable>(m_profileTable);
		ctx->setOutput(symTable, value, m_profileTable->name(), itTABLE, m_profileTable->resource());
	}
	return resource;
}

Ilwis::OperationImplementation* TblFlowPathLongitudinalProfile::create(quint64 metaid, const Ilwis::OperationExpression& expr)
{
	return new TblFlowPathLongitudinalProfile(metaid, expr);
}


quint64 TblFlowPathLongitudinalProfile::createMetadata()
{
	OperationResource operation({ "ilwis://operations/TblFlowPathLongitudinalProfile" });
	operation.setSyntax("TblFlowPathLongitudinalProfile(inputattributemap,inputflowpathsegmentmap,inputsegmentID,inputdistance");
	operation.setDescription(TR("generates a flow path longitudinal profile table"));
	operation.setInParameterCount({ 4 });
	operation.addInParameter(0, itRASTER, TR("attribute map"), TR("input a DEM raster map"));
	operation.addInParameter(1, itLINE, TR("segment map"), TR("input a longest flow path segment map"));
	operation.addInParameter(2, itINT32, TR("Segment ID"), TR("a segment id for produce profile table."), OperationResource::ueNONE);
	operation.addInParameter(3, itDOUBLE, TR("Distance value"), TR("input a distance value which is used a interval for produce attribute points"));

	operation.parameterNeedsQuotes(1);
	operation.setOutParameterCount({ 1 });
	operation.addOutParameter(0, itTABLE, TR("output table"), TR("output an attribute table with a numeric domain"));
	operation.setKeywords("flow path,raster, segment,table");

	operation.checkAlternateDefinition();
	mastercatalog()->addItems({ operation });
	return operation.id();

}


Ilwis::OperationImplementation::State TblFlowPathLongitudinalProfile::prepare(ExecutionContext* ctx, const SymbolTable& st)
{
	OperationImplementation::prepare(ctx, st);
	QString inraster = _expression.parm(0).value();
	QString insegmentmap = _expression.parm(1).value();
	QString insegmentid = _expression.parm(2).value();
	QString indistance = _expression.parm(3).value();


	QString outputName = _expression.parm(0, false).value();

	if (!m_DemMap.prepare(inraster, itRASTER)) {
		ERROR2(ERR_COULD_NOT_LOAD_2, inraster, "");
		return sPREPAREFAILED;
	}

	if (!m_segMap.prepare(insegmentmap, itLINE)) {
		ERROR2(ERR_COULD_NOT_LOAD_2, insegmentmap, "");
		return sPREPAREFAILED;
	}

	m_iSegmentID = insegmentid.toInt();
	m_rDist = indistance.toDouble();

	if (!m_profileTable.prepare(outputName, itFLATTABLE)) {
		ERROR1(ERR_NO_INITIALIZED_1, "output table");
		return sPREPAREFAILED;
	}

	_xsize = m_DemMap->size().xsize();
	_ysize = m_DemMap->size().ysize();

	return sPREPARED;
}


bool TblFlowPathLongitudinalProfile::executeTblFlowPathLongitudinalProfile()
{
	PixelIterator iterDEM = PixelIterator(m_DemMap, BoundingBox(), PixelIterator::fXYZ);

	std::vector<LongitudinalProfile> vlp;
	bool fTransformCoords = m_DemMap->coordinateSystem() != m_segMap->coordinateSystem();

	long count = 0;
	initialize(m_segMap->featureCount(itLINE));

	long iLines = _ysize;
	long iCols  = _xsize;

	LongitudinalProfile lp;
	std::vector<Coordinate> vpoints;
	Coordinate cFrom, cTo;

	Pixel rcFrom, rcTo;
	long i;
	int cvrgIdx = m_segMap->attributeTable()->columnIndex(COVERAGEKEYCOLUMN);

	for (auto feature : m_segMap)
	{
		if (feature->geometryType() == itLINE  && feature->isValid() ) 
		{
			const geos::geom::LineString* ls_geom =
				dynamic_cast<const geos::geom::LineString*>(feature->geometry().get());

			Record rec = feature->record();
			QVariant idval = rec.cell(cvrgIdx);
			int ftid = idval.toInt()+1;

			if (ftid == m_iSegmentID)
			{
				geos::geom::CoordinateSequence* crdSegbuf = ls_geom->getCoordinates();
				cFrom = crdSegbuf->getAt(0);
				cTo = crdSegbuf->getAt(crdSegbuf->size() - 1);

				if (fTransformCoords)
				{
					cFrom = m_DemMap->coordinateSystem()->coord2coord(m_segMap->coordinateSystem(), cFrom);
					cTo = m_DemMap->coordinateSystem()->coord2coord(m_segMap->coordinateSystem(), cTo);
				}

				rcFrom = m_DemMap->georeference()->coord2Pixel(cFrom);
				rcTo = m_DemMap->georeference()->coord2Pixel(cTo);

				rcFrom.z = 0;
				rcTo.z = 0;

				if (*(iterDEM(rcFrom)) >= *(iterDEM(rcTo)))
				{
					for (i = 0; i < crdSegbuf->size() - 1; ++i)
					{
						vpoints.push_back(crdSegbuf->getAt(i));
					}
				}
				else
				{
					for (i = crdSegbuf->size() - 1; i > 0; --i)
					{
						vpoints.push_back(crdSegbuf->getAt(i));
					}
				}
				//delete crdSegbuf;
			}

		}
	}

	long iPointID = 1;
	std::vector<Coordinate>::iterator pos;
	cFrom = *vpoints.begin();
	if (fTransformCoords)
	{
		cFrom = m_DemMap->coordinateSystem()->coord2coord(m_segMap->coordinateSystem(), cFrom);
		cFrom.z = .0;
	}
	rcFrom = m_DemMap->georeference()->coord2Pixel(cFrom);
	rcFrom.z = .0;
	lp.X = iPointID;
	lp.Y = *(iterDEM(rcFrom));
	lp.cd = cFrom;
	vlp.push_back(lp);
	double rDist = 0;

	for (pos = vpoints.begin(); pos < vpoints.end(); ++pos)
	{
		Coordinate cTo = *pos;
		if (fTransformCoords)
		{
			cTo = m_DemMap->coordinateSystem()->coord2coord(m_segMap->coordinateSystem(), cTo);
			cTo.z = .0;	
		}
		rcTo = m_DemMap->georeference()->coord2Pixel(cTo);
		rcTo.z = .0;	
		rDist = rDist + rDistance(cFrom, cTo);
		if (rDist >= (m_rDist - 0.5))
		{
			iPointID++;
			lp.X = iPointID;
			lp.Y = *(iterDEM(rcTo));
			lp.cd = cTo;
			vlp.push_back(lp);
			rDist = 0;
		}
		cFrom = cTo;
	}
	
	CreateTable(vlp);
	vpoints.resize(0);
	vlp.resize(0);

	return true;
}

void TblFlowPathLongitudinalProfile::CreateTable(std::vector<LongitudinalProfile> vlp)
{
	IFlatTable newTable;
	newTable.prepare();
	newTable->addColumn("X", IlwisObject::create<IDomain>("value"), true);
	ColumnDefinition& coldef0 = newTable->columndefinitionRef("X");
	coldef0.datadef().range(new NumericRange(0, 32767, 1));

	newTable->addColumn("Attribute", IlwisObject::create<IDomain>("value"), true);
	ColumnDefinition& coldef1 = newTable->columndefinitionRef("Attribute");
	coldef1.datadef().range(new NumericRange(0, 1.0e300, 0.01));

	ICoordinateDomain crddom0;
	crddom0.prepare();
	crddom0->setCoordinateSystem(m_DemMap->coordinateSystem());
	newTable->addColumn("Coordinates", crddom0, true);

	m_profileTable = newTable;

	std::vector<LongitudinalProfile>::iterator pos;
	long rec = 0;
	for (pos = vlp.begin(); pos < vlp.end(); ++pos)
	{
		LongitudinalProfile lp;
		lp = *pos;
		m_profileTable->setCell("X", rec, QVariant(lp.X));
		m_profileTable->setCell("Attribute", rec, QVariant(lp.Y));

		QVariant crd;
		crd.setValue(lp.cd);
		m_profileTable->setCell("Coordinates", rec, crd);
		rec++;
	}
}

double TblFlowPathLongitudinalProfile::rDistance(Coordinate cd1, Coordinate cd2)
{
	double rDist;
	if (fLatLonCoords())
	{
		double rRadi = rDefaultEarthRadius;
		IConventionalCoordinateSystem projectedCsy = m_DemMap->coordinateSystem().as<ConventionalCoordinateSystem>();
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


bool TblFlowPathLongitudinalProfile::fLatLonCoords()
{
	return m_DemMap->coordinateSystem()->isLatLon();
}



