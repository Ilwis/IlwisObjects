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
#include "domainitem.h"
#include "itemdomain.h"
#include "identifieritem.h"
#include "identifierrange.h"
#include "symboltable.h"
#include "ilwisoperation.h"
#include "basetable.h"
#include "flattable.h"
#include "combinationmatrix.h"
#include "comboclassification.h"

using namespace Ilwis;
using namespace RasterOperations;

REGISTER_OPERATION(ComboClassification)

ComboClassification::ComboClassification()
{

}

ComboClassification::ComboClassification(quint64 metaid, const Ilwis::OperationExpression &expr) : OperationImplementation(metaid, expr){

}

std::vector<int> ComboClassification::raw2index(CombinationMatrix::Axis axis){
    std::vector<int> rawindex;
    for(int x=0; x < _comboMatrix->axisValueCount(axis); ++x){
        QString value = _comboMatrix->axisValue(axis,x);
        Raw raw = _comboMatrix->axisDefinition(axis).domain()->impliedValue(value).toDouble();
        if ( rawindex.size() <= raw)
            rawindex.resize(raw+1, -1);
        rawindex[(qint32)raw] = x;
    }
    return rawindex;
}

bool ComboClassification::execute(ExecutionContext *ctx,SymbolTable& symTable){

    if (_prepState == sNOTPREPARED)
        if((_prepState = prepare(ctx, symTable)) != sPREPARED)
            return false;

    if ( !_inputRasterX->georeference()->isCompatible(_inputRasterY->georeference())) {
        if (!OperationHelperRaster::resample(_inputRasterX, _inputRasterY, ctx)) {
            return ERROR2(ERR_COULD_NOT_CONVERT_2, TR("georeferences"), TR("common base"));
        }
    }

    quint64 currentCount = 0;
    BoxedAsyncFunc binaryMath = [&](const ProcessingBoundingBoxes box, int threadIdx) -> bool {
        PixelIterator iterInX(_inputRasterX, threadIdx, box);
        PixelIterator iterInY(_inputRasterY, threadIdx, box);
        PixelIterator iterOut(_outputRaster, threadIdx, box);

        std::vector<int> raw2indexX = raw2index(CombinationMatrix::aXAXIS);
        std::vector<int> raw2indexY = raw2index(CombinationMatrix::aYAXIS);

		std::vector<bool> allRaws = OperationHelper::collectRawsArray(_outputRaster->datadef().domain());
        PixelIterator iterEnd = end(iterInX);
        while(iterOut != iterEnd) {
            Raw rx = *iterInX;
            Raw ry = *iterInY;
			double outraw = (rx != rUNDEF && ry != rUNDEF) ? _comboMatrix->combo(raw2indexX[rx], raw2indexY[ry]) : rUNDEF;
			if (outraw < allRaws.size())
				allRaws[(quint32)outraw] = true;
			*iterOut = outraw;
            ++iterInX;
            ++iterInY;
            ++iterOut;
            updateTranquilizer(currentCount++, 1000);
        };
		auto itemdom = _outputRaster->datadef().domain().as<ItemDomain<DomainItem>>();
		for (quint32 i = 0; i < allRaws.size(); ++i) {
			if (allRaws[i]) {
				auto item = itemdom->item(i);
				_outputRaster->datadef().range()->as<ItemRange>()->add(item->clone());
			}
		}
        return true;
    };

  


	if (OperationHelperRaster::execute(ctx, binaryMath, { _inputRasterX,_inputRasterY, _outputRaster })) {
		IFlatTable tbl;
		tbl.prepare();
		tbl->addColumn(COVERAGEKEYCOLUMN, _outputRaster->datadef().domain());
		int rec = 0;
		ItemRangeIterator iter(_outputRaster->datadef().range<>().data());
		while (iter.isValid()) {
			SPDomainItem item = (*iter);
			tbl->setCell(0, rec++, item->raw());
			++iter;
		}
		_outputRaster->setAttributes(tbl);

        QVariant value;
        value.setValue<IRasterCoverage>(_outputRaster);
		logOperation(_outputRaster, _expression, {_inputRasterX, _inputRasterY, _comboMatrix});
        ctx->setOutput(symTable,value,_outputRaster->name(), itRASTER,_outputRaster->resource() );
        return true;
    }

    return false;
}

Ilwis::OperationImplementation *ComboClassification::create(quint64 metaid,const Ilwis::OperationExpression& expr){
    return new ComboClassification(metaid, expr);
}

Ilwis::OperationImplementation::State ComboClassification::prepare(ExecutionContext *ctx, const SymbolTable &st){
    OperationImplementation::prepare(ctx,st);
    OperationHelper::check([&] ()->bool { return _inputRasterX.prepare(_expression.input<QString>(0), itRASTER); },
    {ERR_COULD_NOT_LOAD_2,_expression.input<QString>(0), "" } );

    OperationHelper::check([&] ()->bool { return _inputRasterY.prepare(_expression.input<QString>(1), itRASTER); },
    {ERR_COULD_NOT_LOAD_2,_expression.input<QString>(1), "" } );

    if ( _inputRasterX->datadef().domain()->ilwisType() != itITEMDOMAIN || _inputRasterY->datadef().domain()->ilwisType() != itITEMDOMAIN){
        kernel()->issues()->log(TR("Both input rasters must have an item domain as domain"));
        return sPREPAREFAILED;
    }
    OperationHelper::check([&] ()->bool { return _comboMatrix.prepare(_expression.input<QString>(2), itCOMBINATIONMATRIX); },
    {ERR_COULD_NOT_LOAD_2,_expression.input<QString>(2), "" } );

    if ( !_comboMatrix->axisDefinition(CombinationMatrix::aXAXIS).domain()->isCompatibleWith(_inputRasterX->datadef().domain().ptr())){
        kernel()->issues()->log(TR("Domain of the X axis input map doesnt match the combination matrix"));
        return sPREPAREFAILED;
    }

    if (!_comboMatrix->axisDefinition(CombinationMatrix::aYAXIS).domain()->isCompatibleWith(_inputRasterY->datadef().domain().ptr())){
        kernel()->issues()->log(TR("Domain of the Y axis input map doesnt match the combination matrix"));
        return sPREPAREFAILED;
    }

    OperationHelperRaster helper;
    _box = helper.initialize(_inputRasterX, _outputRaster, itRASTERSIZE | itENVELOPE | itCOORDSYSTEM | itGEOREF);

    IDomain dom = _comboMatrix->combinationDef().domain();
    _outputRaster->datadefRef().domain(dom);

    for(quint32 i = 0; i < _outputRaster->size().zsize(); ++i){
     QString index = _outputRaster->stackDefinition().index(i);
        _outputRaster->setBandDefinition(index,DataDefinition(dom));
    }

    initialize(_outputRaster->size().linearSize());

    return sPREPARED;
}


quint64 ComboClassification::createMetadata()
{
    OperationResource operation({"ilwis://operations/crossrasterclassification"});
    operation.setSyntax("crossrasterclassification(inputraster1,inputraster2, combinationmatrix)");
    operation.setDescription(TR("unsupervised classification the value of the output cells is determined by the combination of the input rrasters values in the cells of the combinationmatrix"));
    operation.setInParameterCount({3});
    operation.addInParameter(0,itRASTER,  TR("x - input raster"),TR("Raster which has the values of the x-axis in the combination matrix; it must have an item domain as domain"));
    operation.addInParameter(1,itRASTER,  TR("y - input raster"),TR("Raster which has the values of the y-axis in the combination matrix; it must have an item domain as domain"));
    operation.addInParameter(2,itCOMBINATIONMATRIX, TR("Combo matrix"), TR("Matrix which contains all relevant combinations of the combination of value of the domains of the x and y map"));
    operation.setOutParameterCount({1});
    operation.addOutParameter(0,itRASTER, TR("classified raster"), TR("output raster were all the pixels have been classified according to the combinations"));
    operation.setKeywords("raster, classification");

    operation.checkAlternateDefinition();
    mastercatalog()->addItems({operation});
    return operation.id();
}
