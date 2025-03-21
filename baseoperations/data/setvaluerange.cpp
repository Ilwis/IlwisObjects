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

#include "raster.h"
#include "table.h"
#include "featurecoverage.h"
#include "feature.h"
#include "symboltable.h"
#include "ilwisoperation.h"
#include "pixeliterator.h"
#include "mastercatalog.h"
#include "setvaluerange.h"

using namespace Ilwis;
using namespace BaseOperations;

REGISTER_OPERATION(SetValueRange)

SetValueRange::SetValueRange()
{
}

SetValueRange::SetValueRange(quint64 metaid, const Ilwis::OperationExpression &expr) : OperationImplementation(metaid, expr)
{

}

bool SetValueRange::execute(ExecutionContext *ctx, SymbolTable &symTable)
{
    if (_prepState == sNOTPREPARED)
        if((_prepState = prepare(ctx, symTable)) != sPREPARED)
            return false;

    DataDefinition& datadef = _outputRaster->datadefRef();
    const SPNumericRange numrange = datadef.range<NumericRange>();
    const SPNumericRange rngDomain = datadef.domain()->range<NumericRange>();
    double lmin = _min == rUNDEF ? numrange->min() : std::max(_min, rngDomain->min());
    double lmax = _max == rUNDEF ? numrange->max() : std::min(_max, rngDomain->max());
    double lstep = _step == rUNDEF ? numrange->resolution() : std::max(_step, rngDomain->resolution());
    NumericRange *rng = new NumericRange(lmin, lmax,lstep);
    datadef.range(rng);


    std::function<bool(const ProcessingBoundingBoxes&, int)> SetValrange = [&](const ProcessingBoundingBoxes& box, int threadIdx) -> bool {
        PixelIterator iter(_raster, threadIdx, box);
        PixelIterator iterOut(_outputRaster, threadIdx, box);
        auto end =  iter.end();
        while(iter != end){
            double val = *iter;
            if ( val != rUNDEF){
                double v = rng->ensureExt(val).value<double>();
                val = v;
            }
            *iterOut = val;
            updateTranquilizer(iterOut.linearPosition(), 1000);
            ++iter;
            ++iterOut;
        };
        return true;
    };
	bool ok = OperationHelperRaster::execute(ctx, SetValrange, { _raster, _outputRaster });
    _outputRaster->datadefRef().range<NumericRange>()->resolution(lstep);
    for(int i=0; i < _outputRaster->size().zsize(); ++i){
        _outputRaster->datadefRef(i).range<NumericRange>()->resolution(lstep);
    }

    QVariant value;
    value.setValue<IRasterCoverage>(_outputRaster);
	logOperation(_outputRaster, _expression, { _raster });
    _outputRaster->addDescription(_expression.toString());
    ctx->setOutput(symTable,value,_outputRaster->name(), itRASTER,_outputRaster->resource() );
    return ok;
}

Ilwis::OperationImplementation *SetValueRange::create(quint64 metaid, const Ilwis::OperationExpression &expr)
{
    return new SetValueRange(metaid, expr);
}

Ilwis::OperationImplementation::State SetValueRange::prepare(ExecutionContext *ctx, const SymbolTable &st)
{
    OperationImplementation::prepare(ctx,st);
    QString objectName = _expression.parm(0).value();
    if ( !_raster.prepare(objectName)) {
        ERROR2(ERR_COULD_NOT_LOAD_2, "raster", objectName);
        return sPREPAREFAILED;
    }
    OperationHelperRaster helper;
    _box = helper.initialize(_raster, _outputRaster, itRASTERSIZE | itENVELOPE | itCOORDSYSTEM | itGEOREF|itDOMAIN);

    QString minTxt = _expression.parm(1).value();
    QString maxTxt = _expression.parm(2).value();

    bool ok = true;
    _step = _raster->datadef().domain()->range<NumericRange>()->resolution();
    if ( minTxt.trimmed() != ""){
        _min = minTxt.toDouble(&ok);
    }
    if ( maxTxt.trimmed() != ""){
        _max = maxTxt.toDouble(&ok);
    }
    if ( _expression.parameterCount() == 4) {
        QString stepTxt = _expression.parm(3).value();
        if ( stepTxt.trimmed() != ""){
            _step = stepTxt.toDouble(&ok);
        }
    }
    auto indexes = _raster->stackDefinition().indexes();
    _raster->setDataDefintions(_raster->datadef().domain(),indexes);
    if ( !ok) {
        ERROR2(ERR_ILLEGAL_VALUE_2, TR("parameter"),TR("expression"));
        return sPREPAREFAILED;
    }
    initialize(_outputRaster->size().linearSize());
    return sPREPARED;
}

quint64 SetValueRange::createMetadata()
{

    OperationResource operation({"ilwis://operations/setvaluerange"});
    operation.setSyntax("setvaluerange(objectname,min,max,resolution)");
    operation.setLongName("Set value range raster");
    operation.setDescription(TR("sets the value range of a coverage or column to a new value; all value outside the range will become undefined"));
    operation.setInParameterCount({4});
    operation.addInParameter(0,itRASTER , TR("input object"),TR("input object. If the syntax uses the [] notation it points to a column of a table"));
    operation.addInParameter(1,itNUMBER, TR("minimum"), TR("Color in which the grid lines are drawn, a question mark if this parameter is not used"));
    operation.addInParameter(2,itNUMBER , TR("maximum"), TR("maximum of the new value range. If the value is undefined it will be ignored"));
    operation.addInParameter(3,itNUMBER , TR("resolution"), TR("the minimum distance between to elements of the value range"));
    operation.setOutParameterCount({1});
    operation.addValidation(0,0,"domain=numericdomain");
    operation.addOutParameter(0,itRASTER, TR("raster coverage"));
    operation.setKeywords("numeric,valuerange");

    operation.checkAlternateDefinition();
    mastercatalog()->addItems({operation});
    return operation.id();
}
