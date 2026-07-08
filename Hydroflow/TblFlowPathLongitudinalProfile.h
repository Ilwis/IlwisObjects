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

#ifndef TBLFLOWPATHLONGITUDINALPROFILE_H
#define TBLFLOWPATHLONGITUDINALPROFILE_H

namespace Ilwis {
	namespace Hydroflow {

		struct LongitudinalProfile
		{
			long   X;
			double Y;
			Coordinate  cd;
		};

		class TblFlowPathLongitudinalProfile : public OperationImplementation
		{
		public:
			enum ThresholdMethod { tmValue, tmMap };
			TblFlowPathLongitudinalProfile();

			TblFlowPathLongitudinalProfile(quint64 metaid, const Ilwis::OperationExpression& expr);

			bool execute(ExecutionContext* ctx, SymbolTable& symTable);
			static Ilwis::OperationImplementation* create(quint64 metaid, const Ilwis::OperationExpression& expr);
			Ilwis::OperationImplementation::State prepare(ExecutionContext* ctx, const SymbolTable&);

			static quint64 createMetadata();

		protected:
			IFeatureCoverage _inFlowPathSegmentMap;

			long _xsize, _ysize;

		private:
			void   CreateTable(std::vector<LongitudinalProfile> vlp);
			double rDistance(Coordinate cd1, Coordinate cd2);
			bool fLatLonCoords();
			IFeatureCoverage          m_segMap;
			IRasterCoverage           m_DemMap;
			ITable					  m_profileTable;

			long                      m_iSegmentID;
			double                    m_rDist;

		protected:
			bool executeTblFlowPathLongitudinalProfile();

			NEW_OPERATION(TblFlowPathLongitudinalProfile);
		};


	}
}

#endif //TBLFLOWPATHLONGITUDINALPROFILE_H

