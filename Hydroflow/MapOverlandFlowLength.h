/***************************************************************
 ILWIS integrates image, vector and thematic data in one unique 
 and powerful package on the desktop. ILWIS delivers a wide 
 range of feautures including import/export, digitizing, editing, 
 analysis and display of data as well as production of 
 quality mapsinformation about the sensor mounting platform
 
 Exclusive rights of use by 52°North Initiative for Geospatial 
 Open Source Software GmbH 2007, Germany

 Copyright (C) 2007 by 52°North Initiative for Geospatial
 Open Source Software GmbH

 Author: Jan Hendrikse, Willem Nieuwenhuis,Wim Koolhoven 
 Bas Restsios, Martin Schouwenburg, Lichun Wang, Jelle Wind 

 Contact: Martin Schouwenburg; schouwenburg@itc.nl; 
 tel +31-534874371

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 version 2 as published by the Free Software Foundation.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program (see gnu-gpl v2.txt); if not, write to
 the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 Boston, MA 02111-1307, USA or visit the web page of the Free
 Software Foundation, http://www.fsf.org.

 Created on: 2007-02-8
//////////////////////////////////////////////////////////////////////
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

#ifndef ILWMAPOVERLANDFLOWLENGTH_H
#define ILWMAPOVERLANDFLOWLENGTH_H

namespace Ilwis {
	namespace Hydroflow {

		class MapOverlandFlowLength : public OperationImplementation
		{
		public:

			MapOverlandFlowLength();

			MapOverlandFlowLength(quint64 metaid, const Ilwis::OperationExpression& expr);

			bool execute(ExecutionContext* ctx, SymbolTable& symTable);
			static Ilwis::OperationImplementation* create(quint64 metaid, const Ilwis::OperationExpression& expr);
			Ilwis::OperationImplementation::State prepare(ExecutionContext* ctx, const SymbolTable&);

			static quint64 createMetadata();



		private:
			IRasterCoverage _inDrainRaster;
			IRasterCoverage _inFlowRaster;
			IRasterCoverage _outRaster;

			PixelIterator m_iterDrain;
			PixelIterator m_iterFlow;
			PixelIterator m_iterOut;
				
			IRasterCoverage _iterEmptyRaster;
			PixelIterator iterPos;


			long _xsize, _ysize;
			std::vector<Pixel> m_vCellsOnDivide;
			Pixel m_rcUpstream;
			std::vector<byte> m_vFlowNum;
			std::vector<byte> m_vReceiveNum;
			std::vector<Pixel> m_vFlowLocation;


		private:
			bool executeLandFlowLength();
			Pixel CoordinateStringToPixel(QString coordStr);
			void SplitString(QString s, QString mid, std::vector<long>& results);
			void Lengths2Stream(long iStreamID, Pixel rc, bool);
			void Lengths2Divide(long iStreamID, Pixel rc);
			bool fLatLonCoords();
			double rDistance(Coordinate cd1, Coordinate cd2);
			bool IsEdgeCell(Pixel pxl);

			NEW_OPERATION(MapOverlandFlowLength);
		};
	
	}
}

#endif // MAPOVERLANDFLOWLENGTH_H





// ***************************************************************/
//// MapSlopeLengths.h: interface for the MapSlopeLengths class.
////
////////////////////////////////////////////////////////////////////////
//
//#ifndef ILWMAPOVERLANDFLOWLENGTH_H
//#define ILWMAPOVERLANDFLOWLENGTH_H
//
//#include "Engine\Applications\MAPFMAP.H"
//#include "Engine\Map\Segment\Seg.h"
//#include "Engine\Map\Point\PNT.H"
//#include "Engine\SpatialReference\Coordsys.h"
//#include "Engine\SpatialReference\Ellips.h"
//#include "Engine\SpatialReference\csviall.h"
//#include "Engine\Map\Polygon\POL.H"
//#include "Engine\Table\Col.h"
//#include "Engine\Map\Point\PNT.H"
//#include "LargeVector.h"
//
//#define	sDownstreamCoord  "DownstreamCoord"
//#define	sUpstreamCoord  "UpstreamCoord"
//#define sUpstreamID "UpstreamLinkID"
//
//IlwisObjectPtr * createMapOverlandFlowLength(const FileName& fn, IlwisObjectPtr& ptr, const String& sExpr, vector<void *> parms=vector<void*>() );
//
//class MapOverlandFlowLength : public MapFromMap  
//{
//	friend MapFromMap;
//public:
//	static const char* sSyntax();
//	virtual bool fFreezing();
//	virtual String sExpression() const;
//	static MapOverlandFlowLength* create(const FileName& fn, MapPtr& p, const String& sExpr);
//	MapOverlandFlowLength(const FileName& fn, MapPtr& p);
//	virtual bool fDomainChangeable() const;
//	virtual bool fGeoRefChangeable() const;
//protected:
//	virtual void Store();
//	MapOverlandFlowLength(const FileName& fn, MapPtr& p,
//										const Map& mpdrainageNetwork,
//										const Map& mpFlow); 
//	 ~MapOverlandFlowLength();
//	
//private:
//	LargeVector<LongBuf>   m_vDrainageMap;   
//	LargeVector<ByteBuf>   m_vFlowDir;      //vector for input flow direction 
//	LargeVector<RealBuf>		m_vOutput_s;    //for lengths to stream output
//	LargeVector<RealBuf>		m_vOutput_d;    //for lengths to divide output
//	vector<RowCol> m_vCellsOnDivide;  //Hold the rowcol locations of cells on divide
//	RowCol m_rcUpstream;
//	vector<byte> m_vFlowNum;
//	vector<byte> m_vReceiveNum;
//	vector<RowCol> m_vFlowLocation;
//	Map m_mpFlow;
//	void init();
//	bool IsEdgeCell(long iRow, long iCol);
//	void CompitableGeorefs(FileName fn, Map mp1, Map mp2);
//	void Lengths2Stream(long iStreamID, RowCol rc, bool);
//	void Lengths2Divide(long iStreamID, RowCol rc);
//	bool fLatLonCoords();
//	bool fEllipsoidalCoords();
//	double rDistance(Coord cd1, Coord cd2);
//};
//
//#endif // ILWMapOverlandFlowLength_H
