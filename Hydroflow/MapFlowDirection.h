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

#ifndef MAPFLOWDIRECTION_H
#define MAPFLOWDIRECTION_H

namespace Ilwis {
    namespace Hydroflow {

	class MapFlowDirection : public OperationImplementation
		{
		public:

			enum FlowMethod { fmSlope, fmHeight };
			
			MapFlowDirection();

			MapFlowDirection(quint64 metaid, const Ilwis::OperationExpression& expr);

			bool execute(ExecutionContext* ctx, SymbolTable& symTable);
			static Ilwis::OperationImplementation* create(quint64 metaid, const Ilwis::OperationExpression& expr);
			Ilwis::OperationImplementation::State prepare(ExecutionContext* ctx, const SymbolTable&);

			static quint64 createMetadata();


		FlowMethod m_fmMethods;
		bool       m_fParallel;
		

		private:
			IRasterCoverage _inRaster;
			IRasterCoverage _outRaster;

			long _xsize, _ysize;
			std::vector<int> m_vFlatIndices;	//store a continuous flat area  	
			std::vector<byte> m_vDirection;
			std::vector<byte> m_vFlowSelf;

			double* m_inDem;
			double* m_flowDem;
			int* m_flag;

		private:

			inline int idx(int x, int y) const
			{
				return y * _xsize + x;
			}
	
			void InitPars();
			long m_ContFlat;
			long iLookUp(double rMax, int iCout, std::vector<int>&);
			bool isInOneEdge(int iPos1, int iPos2, int iPos3, std::vector<int>& vPos);
			void executeFlowDirection();
			bool onEdge(int x, int y);
			void FillArray(int x, int y, std::vector<double>& vValue);
			void TreatFlatAreas();
			void LocateOutlets(int sx, int sy, std::vector<int>& outlets);
			void SetFlowsInFlatArea(std::vector<int>& vOutlets);
			bool isEven(int elem);
			double rComputeSlope(double rCurH, double rNbH, int iPos);
			double rComputeHeightDifference(double rCurH, double rNbH);
			double rFindMaxLocation(std::vector<double>& vValue, std::vector<int>& vPos, int& iCout);


            NEW_OPERATION(MapFlowDirection);
        };

	
	};

	struct Cell
	{
		int x;
		int y;
		int val;

		Cell(int _x = 0, int _y = 0, int _val = 0)
			: x(_x), y(_y), val(_val) {
		}

		void set(int _x, int _y) {
			x = _x;
			y = _y;
		}
	};

	class FlowDirectionAlgorithm 
	{
	public:
		enum Method { slope, height };
		enum FlowDirection { NW, N, NE, W, E, SW, S, SE };

	private:
		double* m_inDem;
		double* m_flowDem;
		int* m_flag;

	private:
		FlowDirectionAlgorithm::Method method;
		byte Location[8];
		long lines, columns;
		byte increment;
		byte flatcell;
		byte flag;
		long _xsize, _ysize;

		inline int idx(int x, int y) const
		{
			return y * _xsize + x;
		}

		Method methodValueOf(QString val);
		bool onEdge(int x, int y);
		bool isEven(int elem);
		bool isInOneEdge(const std::vector<FlowDirection>& listPos,FlowDirection fd1, FlowDirection fd2, FlowDirection fd3);
		bool hasFlow(unsigned char flowdirection);

		double maxAdj(int x, int y, double listVal[]);
		double maxAdj(int x, int y, double* gradient, double listVal[]);
		void findDirection(double listA[], double val,std::vector<FlowDirection>& listPos);
		FlowDirection getFlowDirection(const std::vector<FlowDirection>& listPos);
		void locateOutlet(int x, int y,std::vector<Cell>& flatList,std::vector<Cell>& outList);
		void imposeGradient2LowerElevation(std::vector<Cell>& outletList,std::vector<Cell>& flatList,double* gradient);
		void imposeGradientFromElevation(std::vector<Cell>& flatList,double* gradient);
		void combineGradient(double* grd1, double* grd2,std::vector<Cell>& flatList);
		void assignFlowInFlat(std::vector<Cell>& flatList,double* gradient);
		void iniGradient(double* grd1, double* grd2,std::vector<Cell>& flatList);
		double computeSlope(double h1, double h2, int pos);
		double computeHeightDifference(double h1, double h2);
		FlowDirection mapFlowLocation(int pos);

	public:
		static int noflow;
			FlowDirectionAlgorithm(double* dem,double* flow, long xsz,long ysz);

			void calculate(QString method);
	};

}

#endif // MAPFLOWDIRECTION_H

