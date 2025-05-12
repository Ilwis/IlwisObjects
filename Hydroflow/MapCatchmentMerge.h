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

			struct OutletLocation
			{
				OutletLocation() {}
				long   StreamID;
				Pixel pxl;
				long  StreveOrder;
				double rLen1;
				double rLen2;
				bool   isOnNode;
				bool   fExtractOverlandFlowPath;
			};

			struct DrainageAtt
			{
				DrainageAtt() {}
				long ID;
				std::vector<long> UpstreamID;
				Pixel DownStreamCoord;
				Pixel UpstreamCoord;
				Pixel TostreamCoord;
				long iDownStreamID;
				long iOrder;
				long CatchmentLink;
				double rLenght;
			};

			struct AttLongestPath
			{
				Coordinate UpstreamCoord;
				Coordinate DownstreamCoord;
				double rLength;
			};

		class MapCatchmentMerge : public OperationImplementation
		{
		public:
			MapCatchmentMerge();
			MapCatchmentMerge(quint64 metaid, const Ilwis::OperationExpression& expr);

			bool execute(ExecutionContext* ctx, SymbolTable& symTable);
			static Ilwis::OperationImplementation* create(quint64 metaid, const Ilwis::OperationExpression& expr);
			Ilwis::OperationImplementation::State prepare(ExecutionContext* ctx, const SymbolTable&);

		static quint64 createMetadata(Ilwis::OperationResource& operation);
		protected:
			bool executeCatchmentMerging();
			void CreateTable();
			void AddLink2StreamSegments();
			void InitOutletVector();
			Pixel CoordinateStringToPixel(QString coordStr);
			bool IsEdgeCell(Pixel pxl);
			bool fRelocatOutlet(Pixel& pxl, int);
			double CalculateLength(Pixel pxl1, Pixel pxl2, long);
			int GetDownStreamCell(Pixel& pxl);
			bool fLatLonCoords();
			double rDistance(Coordinate cd1, Coordinate cd2);
			void EvaluateJointOutletsbyOrder();
			void PutOutlet(long id, Pixel rc, long iFlag, bool fExtractOverlandFlowPath);
			long MergeCatchment(Pixel pxl, long iFlag, bool fExtractOriginalOrder);
			void BuildUpLinkCatchment(Pixel pxl, int iFlow, long iFlag);
			void IdentifyStreamsInCatchment(Pixel pxl);
			void UpdateUpLinkCatchment(long id);
			void UpdateDownLinkCatchment(long id);
			void UpdateLink2StreamSegments(long iCatchmentID, OutletLocation ol);
			double GetSplitSegmentLength(long iStreamID);
			void ExtractOriginalOrder(long& id);
			bool IsUpstreamsMerged(std::vector<long> vUpstreams);
			void Merge(long id, DrainageAtt datt, bool fExtractOriginalOrder);
			bool isDigitStr(QString str);
			void ComputeOtherAttributes();
			void ExtractUpstreamFlowPath(Pixel rc, long id);
			Coordinate ComputeCenterDrainage(long iDrainageID, double rLength, IFeatureCoverage sm);

			void SplitString(QString s, QString mid, std::vector<long>& results);
			void SplitString(QString s, QString mid, std::vector<double>& results);
			Coordinate StringToCoordinate(QString s, QString rpls);
			QString CoordinateToString(Coordinate crd);
			void InitPars();

			Coordinate StoreSegment(IFeatureCoverage smpFrom, long id, long val);
			Coordinate StoreSegment(IFeatureCoverage smpFrom, geos::geom::CoordinateSequence* crdbuf, long id, long val);

			void StoreSourceSegment(long val);
			void StoreSourceSegment(geos::geom::CoordinateSequence* crdbuf, long val);

			Coordinate SplitSegment(IFeatureCoverage smpFrom, long iRaw, double disval, long id, Pixel rc);
			Coordinate SplitSegment(IFeatureCoverage smpFrom, geos::geom::CoordinateSequence* crdbuf, long iRaw, double disval, long id, Pixel rc);

			void CreateTableLongestFlowPath(std::vector<AttLongestPath> vAtt);
			double rComputeSinuosity(double rLength, double rStraightLenght);
			void ExtractSegments();
			void CleanSegment(IFeatureCoverage smpTo, IFeatureCoverage smpFrom);
			void CreateTableSegmentsExtracted();
			std::vector<long> GetSegmentIDsExtracted(int iCatch);
			void CreatePolygonTableElements();
		protected:
			void InitInOutMaps();

		protected:

			Pixel lastrc;
			long _xsize, _ysize;
			INamedIdDomain _outDomain;
			NamedIdentifierRange* _idrange;

			bool m_UseOutlets;
			bool m_UseStrahler;
			bool m_useExtraOrder;

			bool m_includeUndefine;

			long m_iStreamOrders;
			long m_iMaxOrderNumber;
			double m_rSourceWaterFlowPathLen;
			quint32 m_totalcount;


			QString m_sOrderSystem;

			long m_iOutletVal;
			QString m_sUpLinkCatchment;
			std::vector<long> m_vUpCatchmentIDs;
			std::vector<long> m_vStreamsInCatchment;
			std::vector<DrainageAtt> m_vDrainageAtt;
			std::vector<Coordinate> m_vStreamCoord;	//store coordinate of cells in a linksplited   
			std::vector<Coordinate> m_vStream;
			std::vector<OutletLocation> m_vOutlet;
			std::vector<byte> m_vFlowSelf;

			std::map<const geos::geom::Geometry*, int> lsidmap;

			ITable _outputTable;
			IRasterCoverage  _inDrngOrderRaster;
			IRasterCoverage  _inFldRaster;
			IRasterCoverage  _inAccRaster;
			IRasterCoverage  _inDemRaster;
			IRasterCoverage  _outMergeRaster;

			PixelIterator m_iterFld;
			PixelIterator m_iterOut;
			PixelIterator m_iterInDrng;

			IFeatureCoverage _inPointMap;

			IFeatureCoverage _internelDranageSegmentMap;

			IFeatureCoverage _outputExtraxtedSegmentMap;
			IFeatureCoverage _outputPolygonMap;
			IFeatureCoverage _longestPathFeature;

			ITable _outputPathTable;
			ITable _outputExtractSegTable;
			ICoordinateSystem _csy;
			IGeoReference _inputgrf;

			std::vector<geos::geom::Geometry*> lines;

			NEW_OPERATION(MapCatchmentMerge);
		};


		class MapCatchmentMergeWithOutlet : public MapCatchmentMerge
		{

		public:
			MapCatchmentMergeWithOutlet(quint64 metaid, const Ilwis::OperationExpression& expr);

			static Ilwis::OperationImplementation* create(quint64 metaid, const Ilwis::OperationExpression& expr);
			Ilwis::OperationImplementation::State prepare(ExecutionContext* ctx, const SymbolTable&);

			static quint64 createMetadata();

			NEW_OPERATION(MapCatchmentMergeWithOutlet);
		};

		class MapCatchmentMergeWithStreamOrder : public MapCatchmentMerge
		{

		public:

			MapCatchmentMergeWithStreamOrder(quint64 metaid, const Ilwis::OperationExpression& expr);

			static Ilwis::OperationImplementation* create(quint64 metaid, const Ilwis::OperationExpression& expr);
			Ilwis::OperationImplementation::State prepare(ExecutionContext* ctx, const SymbolTable&);
			static quint64 createMetadata();

			NEW_OPERATION(MapCatchmentMergeWithStreamOrder);
		};


	}
}

#endif // ILWMapCatchmentMerge_H
