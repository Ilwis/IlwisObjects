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
			void SplitString(QString s, QString mid, std::vector<long>& results);
			void Lengths2Stream(long iStreamID, Pixel rc, bool);
			void Lengths2Divide(long iStreamID, Pixel rc);
			bool fLatLonCoords();
			double rDistance(Coordinate cd1, Coordinate cd2);
			bool IsEdgeCell(Pixel pxl);
			void InitFlowNums(std::vector<byte>& vReceiveNum);

			NEW_OPERATION(MapOverlandFlowLength);
		};
	
	}
}

#endif // MAPOVERLANDFLOWLENGTH_H



