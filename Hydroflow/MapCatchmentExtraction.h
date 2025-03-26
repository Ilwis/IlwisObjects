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

#ifndef MAPCATCHMENTEXTRACTION_H
#define MAPCATCHMENTEXTRACTION_H

namespace geos {
	namespace geom {
		class Geometery;
	}
}

namespace Ilwis {
	namespace Hydroflow {

#define sUpstreamLink  "UpstreamLinkID"
#define	sDownstreamLink  "DownstreamLinkID"
#define	sDownstreamCoord  "DownstreamCoord"

		struct AttCols
		{
			AttCols() {}
			long         DrainageID;
			std::vector<long> UpstreamLink;
			long	 DownstreamLink;
			Pixel DownstreamCoord;
		};

		struct AttUpstreamLink
		{
			AttUpstreamLink() {}
			long         DrainageID;
			std::vector<long> UpstreamLink;
		};

		class MapCatchmentExtraction : public OperationImplementation
		{
		public:
			MapCatchmentExtraction();
			MapCatchmentExtraction(quint64 metaid, const Ilwis::OperationExpression& expr);

			bool execute(ExecutionContext* ctx, SymbolTable& symTable);
			static Ilwis::OperationImplementation* create(quint64 metaid, const Ilwis::OperationExpression& expr);
			Ilwis::OperationImplementation::State prepare(ExecutionContext* ctx, const SymbolTable&);

			static quint64 createMetadata();

		private:
			bool executeCatchmentExtraction();
			bool IsEdgeCell(Pixel pxl);
			//void CompitableGeorefs(FileName fn, Map mp1, Map mp2);
			void GetAttributes();
			long DelineateCatchment(PixelIterator iterFld, PixelIterator iterFlag,Pixel pxl, long iFlag);
			long DelineateCatchment(Pixel pxl, long iFlag);
			long FindDownstreamIndex(long DownstreamID);
			void UpdateUpstreamLinkID(long DrainageID, long UpstreamID);
			void EraseDrainage(long DrainageID);
			void ComputeCatchmentAttributes();
			bool fLatLonCoords();
			void ComputeCenterDrainage();
			double GetDistance(Pixel& rc);
			double rDistance(Coordinate cd1, Coordinate cd2);
			bool fEllipsoidalCoords();
			void SetAttributeTable();

			void SplitString(QString s,QString mid,std::vector<long>& results);

		private:
			IRasterCoverage _inDrngOrderRaster;
			IRasterCoverage _inFldRaster;

			IRasterCoverage _outCatchmentRaster;
			IRasterCoverage _flagRaster;
			
			IFeatureCoverage _outputfeatures;
			ITable _outputTable;

			INamedIdDomain _outDomain;
			NamedIdentifierRange* _idrange;

			ICoordinateSystem _csy;
			IGeoReference _inputgrf;

			long _xsize, _ysize;
			std::vector<AttCols> m_vRecords;
			std::vector<long> m_vDrnIDs;

			std::vector<AttUpstreamLink> m_vvUpstreamLinks;

		private:
			struct PrevLinePoint {
				int _networkLink = iUNDEF;
				int _shadowLink = iUNDEF;

				bool isValid() const { return _networkLink != iUNDEF || _shadowLink != iUNDEF; }

			};

			struct NetworkPoint {
				int _x = iUNDEF;
				int _y = iUNDEF;
				int _links[4]; // 0=left. 1=up; 2=right; 3=down

				bool isValid() const { return _x != iUNDEF && _y != iUNDEF; }
			};

			std::vector<NetworkPoint> _points;

			NetworkPoint makeConnection(const Pixel& pCenter, bool isTemp, std::vector<MapCatchmentExtraction::PrevLinePoint>& previousLine);




			NEW_OPERATION(MapCatchmentExtraction);
		};

	}

}

#endif //MAPCATCHMENTEXTRACTION_H


