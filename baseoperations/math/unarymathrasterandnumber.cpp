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
#include "columndefinition.h"
#include "table.h"
#include "basetable.h"
#include "flattable.h"
#include "symboltable.h"
#include "ilwisoperation.h"
#include "unarymath.h"
#include "unarymathrasterandnumber.h"


using namespace Ilwis;
using namespace BaseOperations;

UnaryMathRasterAndNumber::UnaryMathRasterAndNumber() {

}

UnaryMathRasterAndNumber::UnaryMathRasterAndNumber(quint64 metaid, const Ilwis::OperationExpression& expr, const QString& outputdom,UnaryFunction func) :
    UnaryMath(metaid, expr,outputdom, func),
    _number(rUNDEF)
{

}

bool UnaryMathRasterAndNumber::execute(ExecutionContext *ctx, SymbolTable& symTable)
{
    if (_prepState == sNOTPREPARED)
        if((_prepState = prepare(ctx, symTable)) != sPREPARED)
            return false;

    if ( _case == otSPATIAL) {
        BoxedAsyncFunc unaryFun = [&](const ProcessingBoundingBoxes& box, int threadIdx) -> bool {
            PixelIterator iterIn(_inputGC, threadIdx, _box);
            PixelIterator iterOut(_outputGC, threadIdx, _box);

            double v_in = 0;
            std::for_each(iterOut, iterOut.end(), [&](double& v){
                if ( (v_in = *iterIn) != rUNDEF) {
                    v = _unaryFun(v_in);
                }
                ++iterIn;
            });
            return true;
        };


        bool resource = OperationHelperRaster::execute(ctx, unaryFun, { _outputGC });
        if ( resource && ctx != 0) {
            QVariant value;
            value.setValue<IRasterCoverage>(_outputGC);
             _outputGC->addDescription(_expression.toString());
			 logOperation(_outputGC, _expression, { _inputGC });
            ctx->setOutput(symTable,value,_outputGC->name(), itRASTER,_outputGC->resource() );
        }
    } else if (_case == otNUMBER) {
        double v = _unaryFun(_number);
        ctx->setOutput(symTable, QVariant(v), sUNDEF, itDOUBLE, Resource());

    }


    return true;
}

NumericRange *UnaryMathRasterAndNumber::constructRangeFrom(const SPNumericRange &range)
{
    return static_cast<NumericRange *>(range->clone());
}

OperationImplementation::State UnaryMathRasterAndNumber::prepare(ExecutionContext *ctx,const SymbolTable &st)
{
    OperationImplementation::prepare(ctx,st);
    IlwisTypes ptype = _expression.parm(0).valuetype();


    if ( hasType(ptype,itNUMBER) ) {
        _case = otNUMBER;
        bool ok;
        _number = _expression.parm(0).value().toDouble(&ok);
        if (!ok) {
            ERROR2(ERR_NO_OBJECT_TYPE_FOR_2,"Numerical value", "UnaryMathRasterAndNumber operation");
            _number = rUNDEF;
            return sPREPAREFAILED;
        }
        return sPREPARED;

    } else if ( hasType(ptype,itRASTER)) {
        QString raster = _expression.parm(0).value();

        if (!_inputGC.prepare(raster)) {
            ERROR2(ERR_COULD_NOT_LOAD_2,raster,"");
            return sPREPAREFAILED;
        }
        OperationHelperRaster helper;
        _box = helper.initialize(_inputGC, _outputGC, itRASTERSIZE | itENVELOPE | itCOORDSYSTEM | itGEOREF);
        if ( !_outputGC.isValid()) {
            ERROR1(ERR_NO_INITIALIZED_1, "output rastercoverage");
            return sPREPAREFAILED;
        }
        QString outputName = _expression.parm(0,false).value();
        if ( outputName != sUNDEF)
            _outputGC->name(outputName);

        auto nrange = _inputGC->datadef().range<NumericRange>();
        if (nrange.isNull())
            return sPREPAREFAILED;

        NumericRange *newRange = constructRangeFrom(nrange);
        newRange->resolution(0); // all unary operation are floating point

        IDomain dom;
         if(!dom.prepare(_outputDomain))
             return sPREPAREFAILED;

         _outputGC->datadefRef().domain(dom);
         _outputGC->datadefRef().range(newRange);
        for(quint32 i=0; i<_outputGC->size().zsize(); ++i){
                QString index = _outputGC->stackDefinition().index(i);
             _outputGC->setBandDefinition(index,{dom, newRange->clone()});
         }
        _case = otSPATIAL;
        return sPREPARED;
    }
    return sNOTPREPARED;
}

Resource UnaryMathRasterAndNumber::populateMetadata(const QString& item, const QString& longname) {
    OperationResource operation(item);
    QString shortname = item.split("/").last();
    operation.setSyntax(QString("%1(raster|number)").arg(shortname));
    operation.setDescription(TR("generates a new numerical rastercoverage/featurecoverage based on the operation, applied to all the pixels"));
    operation.setLongName(longname);
    operation.setInParameterCount({1});
    operation.addInParameter(0,itRASTER|itNUMBER,"input value");
    operation.setOutParameterCount({1});
    operation.addOutParameter(0,itRASTER|itNUMBER, "output name",sUNDEF, itINTEGER | itFLOAT | itDOUBLE);
    operation.setKeywords("raster,numeric,math");

    return operation;

}



