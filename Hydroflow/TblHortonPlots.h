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

#ifndef ILWTABLEHORTONPLOTS_H
#define ILWTABLEHORTONPLOTS_H

#include "realmat.h"

using namespace std;

namespace Ilwis {
    namespace Hydroflow {

        class LeastSquaresFit
        {
        protected:
            LeastSquaresFit(const CVector& cvXData, const CVector& cvYData, int iTerms);
        public:
            virtual ~LeastSquaresFit();
            double rStdDev() { return _rStdDev; }
            virtual QString sFormula() const = 0;
            virtual double rCalc(double X) const = 0;
        protected:
            void Compute();
            virtual void Transform();
            virtual void CreateBasisFunctions() = 0;
            virtual void InverseTransform(double rYFit);
            virtual void TransformSolution();
            const CVector cvX, cvY;
            CVector cvSolution, cvW, cvZ;
            RealMatrix mtBasis;
            double rConstant, rMultiplier;
            const int iNrPoints, iNrTerms;
        private:
            CVector LeastSquares(const RealMatrix& A, const CVector& b);
            void InitializeAndFormBasisVectors();
            void CreateAndSolveEquations();
            void FinalizeSolution();
            double _rStdDev;
        };

        class LeastSquaresFitPolynomial : public LeastSquaresFit
        { // Y = a0 + a1 * X + a2 * X^2 + ...
        public:
            LeastSquaresFitPolynomial(const CVector& cvXData, const CVector& cvYData, int iOrder);
            virtual QString sFormula() const;
            virtual double rCalc(double X) const;
        protected:
            virtual void Transform();
            virtual void CreateBasisFunctions();
            virtual void TransformSolution();
        private:
            void WritePoly(int iTerm, ostream& os) const;
            double rFuncPoly(int iTerm, double rX) const;
        };

        class ChebyshevFitter 
        {
        public:
            ChebyshevFitter() {};
            vector<double> LeastSquaresFitPolynomial(const vector<double>& x, const vector<double>& y, int n_terms);
        private:
            vector<double> Chebyshev_fit(const vector<double>& x, const vector<double>& y, int n_terms);
            void FinalizeSolution(vector<double>& cvSolution);
            double rFuncPoly(int iTerm, double rX) const;
            double rCalc(double rX) const;
            double Map_to_chebyshev_domain(double x, double xmin, double xmax);
            void TransformSolution();
            double Chebyshev(int k, double t);

            vector<vector<double>> mtBasis;
            vector<double> m_cvSolution;
            vector<double> m_x;
            vector<double> m_y;
            vector<double> m_coef;
            double m_rConstant;
            double m_rMultiplier;
            long m_iNrPoints;
            int m_iNrTerms;
            double m_rStdDev;

        };

        struct LongitudinalProfile
        {
            long   X;
            double Y;
            Coordinate  cd;
        };

        class TableHortonPlots : public OperationImplementation
        {
        public:
            enum ThresholdMethod { tmValue, tmMap };
            TableHortonPlots();

            TableHortonPlots(quint64 metaid, const Ilwis::OperationExpression& expr);

            bool execute(ExecutionContext* ctx, SymbolTable& symTable);
            static Ilwis::OperationImplementation* create(quint64 metaid, const Ilwis::OperationExpression& expr);
            Ilwis::OperationImplementation::State prepare(ExecutionContext* ctx, const SymbolTable&);

            static quint64 createMetadata();

        protected:
            IFeatureCoverage _inFlowPathSegmentMap;

            long _xsize, _ysize;

        private:
            void   CreateTable(long iSize);
            double rDistance(Coordinate cd1, Coordinate cd2);
            bool fLatLonCoords();
            void SplitString(QString s, QString mid, std::vector<long>& results);
            void SplitString(QString s, QString mid, std::vector<double>& results);
            std::vector<long> GetOrderNumber(std::vector<long> vIDs);
            void ComputeHortonPlots(long CatchmentNr, long iMaxOrder, std::vector<long> vDrainageIDs, std::vector<long> vOrder, bool fDrainageLen,
                std::vector<double> vDrainageLen);
            bool fUpstreamSameOrder(long id, long iOrder, std::vector<long> vDrainageIDs);
            void Regress(long CatchmentNr, std::vector<QVariant>  colY, std::vector<QVariant>  colX, std::vector<QVariant>& outCol, int iOrders, bool usAllOrders);
            void CalculateRB(std::vector<QVariant> col, std::vector<QVariant>& sCol, bool isRB);

            IRasterCoverage           m_catchmentMap;
            IRasterCoverage           m_drainageMap;
            IRasterCoverage           m_catchMergeMap;

            ITable					  m_hortonTable;
            ITable					  m_ratioTable;

        protected:
            bool executeTableHortonPlots();

            NEW_OPERATION(TableHortonPlots);
        };


    }
}

#endif //ILWTABLEHORTONPLOTS_H



//
//////////////////
//#ifndef ILWTABLEHORTONPLOTS_H
//#define ILWTABLEHORTONPLOTS_H
//
//#include "Engine\Applications\TBLVIRT.H"
//#include "Engine\Map\Raster\Map.h"
//#include "LargeVector.h"
//
//IlwisObjectPtr * createTableHortonPlots(const FileName& fn, IlwisObjectPtr& ptr, const String& sExpr, vector<void *> parms );
//
//struct coefficients
//{
//    double a;
//    double b;
//};
//
//struct PlotPoints
//{
//    double x;
//    double y;
//};
//
//class DATEXPORT TableHortonPlots: public TableVirtual
//{
//    friend class TableVirtual;
//
//    public:
//        static const char*      sSyntax();
//        Column                  colMap1;
//        Column                  colMap2;
//        static TableHortonPlots*  create(const FileName& fn, TablePtr& p, const String& sExpression);
//        TableHortonPlots(const FileName& fn, TablePtr& p);
//
//    protected:
//        TableHortonPlots(const FileName& fn, TablePtr& p, 
//                         const Map& mpCatchment,
//                         const Map& mpMergedCatchment,
//                         const Map& mpDrainage);
//                         
//                   
//        ~TableHortonPlots();
//
//        virtual void     Store();
//        virtual bool     fFreezing();
//        virtual String   sExpression() const;
//    
//    private:
//        void                      Init();
//        void                      CreateColumns(long iSize);
//        vector<long>              GetOrderNumber(vector<long> vIDs);
//        void                      AddDomainItem(Domain dm,long item);
//        bool                      fUpstreamSameOrder(long,long,vector<long>);
//        void                      ComputeHortonPlots(long, long, vector<long>, vector<long>,bool, vector<double>);
//        //coefficients              Regress(int,vector<PlotPoints>);
//        void                      Regress(long CatchmentNr, Column, Column,String,int, bool);
//        void                      CalculateRB(Column, String, bool);
//        Map                       m_mpCatchment;
//        Map                       m_mpDrainage;
//        Map                       m_mpMergedCatchment;
//        LargeVector<LongBuf>     m_vCatchment;   
//	    LargeVector<LongBuf>     m_vDrainage;        
//        LargeVector<LongBuf>     m_vMergedCatchment;
//        Table                     m_tbl_ratio;
//};
//
//#endif 
//
//
//
//
