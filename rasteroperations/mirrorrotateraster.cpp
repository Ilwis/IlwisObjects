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
#include "basetable.h"
#include "flattable.h"
#include "domainitem.h"
#include "itemdomain.h"
#include "identifieritem.h"
#include "identifierrange.h"
#include "symboltable.h"
#include "ilwisoperation.h"
#include "operationhelpergrid.h"
#include "geometryhelper.h"
#include "mirrorrotateraster.h"

using namespace Ilwis;
using namespace RasterOperations;

REGISTER_OPERATION(MirrorRotateRaster)

MirrorRotateRaster::MirrorRotateRaster()
{
}

MirrorRotateRaster::MirrorRotateRaster(quint64 metaid, const Ilwis::OperationExpression &expr) : OperationImplementation(metaid, expr)
{

}
Pixel MirrorRotateRaster::outputPixel(const Pixel& pos, int width, int height){
    switch(_method)    {

    case tmRotate90:
        return Pixel(std::round((pos.y - (height - 1) / 2.0) + (width - 1) / 2.0),
                     std::round(- (pos.x - (width - 1) / 2.0) + (height - 1) / 2.0),
                     pos.z);
    case tmRotate180:
        return Pixel((width - 1) - pos.x,
                     (height - 1) - pos.y,
                     pos.z);
    case tmRotate270:
         return  Pixel(
                     static_cast<int>(std::round((pos.y - (height - 1) / 2.0) + (width - 1) / 2.0)),
                     static_cast<int>(std::round(-(pos.x - (width - 1) / 2.0) + (height - 1) / 2.0))
                 );
    case tmMirrorDiagonal:
        return Pixel(pos.y, pos.x, pos.z);
    case tmTranspose:
        return Pixel(width - 1 - pos.y, height - pos.x, pos.z);
    case tmMirrorHorizontal:
        return Pixel( pos.x,  height - 1 - pos.y, pos.z);
    case tmMirrorVertical:
        return Pixel(width - 1 - pos.x, pos.y, pos.z);
    default:
        break;
    }
}


bool MirrorRotateRaster::execute(ExecutionContext *ctx, SymbolTable &symTable)
{
    if (_prepState == sNOTPREPARED)
        if((_prepState = prepare(ctx,symTable)) != sPREPARED)
            return false;

    PixelIterator itertt = PixelIterator(_inputRaster);


    BoxedAsyncFunc Transform = [&](const ProcessingBoundingBoxes box, int threadIdx) -> bool {
        PixelIterator iterIn(_inputRaster);

        PixelIterator iterOut(_outputRaster);
        auto end = iterIn.end();
        auto w = _outputRaster->size().xsize();
        auto h = _outputRaster->size().ysize();
        for(; iterIn != end; ++iterIn){
            double v = *iterIn;
            auto pos = iterIn.position();
            auto pix = outputPixel(pos, w, h);
            iterOut.position(pix);
            *iterOut = v;
        }
        return true;

    };

	bool ok = OperationHelperRaster::execute(ctx, Transform, { _inputRaster, _outputRaster });

    if ( ok && ctx != 0) {
        QVariant value;
        value.setValue<IRasterCoverage>(_outputRaster);
		logOperation(_outputRaster, _expression, {_inputRaster});
        ctx->setOutput(symTable,value,_outputRaster->name(), itRASTER, _outputRaster->resource() );
    }
    return ok;
}

Ilwis::OperationImplementation *MirrorRotateRaster::create(quint64 metaid, const Ilwis::OperationExpression &expr)
{
    return new MirrorRotateRaster(metaid, expr);
}

Ilwis::OperationImplementation::State MirrorRotateRaster::prepare(ExecutionContext *ctx, const SymbolTable &st)
{
    OperationImplementation::prepare(ctx,st);
    QString raster = _expression.parm(0).value();
    QString outputName = _expression.parm(0,false).value();
    QString method = _expression.parm(1).value().toLower();

    if (!_inputRaster.prepare(raster, itRASTER)) {
        ERROR2(ERR_COULD_NOT_LOAD_2,raster,"");
        return sPREPAREFAILED;
    }
    std::map<QString, TransPoseMethod> methods={{"mirrhor",tmMirrorHorizontal},{"mirrvert",tmMirrorVertical},
                                                {"mirrdiag",tmMirrorDiagonal},{"transpose",tmTranspose},{"rotate90",tmRotate90},
                                                {"rotate180",tmRotate180},{"rotate270",tmRotate270}};
    auto iter = methods.find(method);
    if ( iter == methods.end()){
        ERROR2(ERR_NOT_FOUND2,method, TR("in method for mirrorrotate"));
        return sPREPAREFAILED;
    }
    _method = iter->second;
    Size<> sz = _inputRaster->size();
    Envelope outputenv = _inputRaster->envelope();

    if ( _method == tmTranspose || _method == tmRotate90 || _method == tmRotate270 || _method == tmMirrorDiagonal){
        sz = Size<>(sz.ysize(), sz.xsize(), sz.zsize());
        Coordinate center = (outputenv.max_corner() + outputenv.min_corner()) / 2.0;
        auto rotated = GeometryHelper::rotate2d(center,90,(std::vector<Coordinate>)outputenv);
        outputenv = {rotated};
    }

    _outputRaster = OperationHelperRaster::initialize(_inputRaster,itRASTER,itCOORDSYSTEM | itDOMAIN);
    _outputRaster->gridRef()->prepare(_outputRaster->id(),sz);

    QString grfs = QString("code=georef:type=corners,csy=%1,envelope=%2,gridsize=%3,cornerofcorners=yes")
            .arg(_outputRaster->coordinateSystem()->id())
            .arg(outputenv.toString())
            .arg(sz.toString());
    _outputRaster->georeference(grfs);
    if (outputName != sUNDEF)
        _outputRaster->name(outputName);


    return sPREPARED;
}

quint64 MirrorRotateRaster::createMetadata()
{
    OperationResource operation({"ilwis://operations/mirrorrotateraster"});
    operation.setSyntax("mirrorrotateraster(inputraster,mirrortype=mirrhor | mirrvert | mirrdiag | transpose | rotate90 | rotate180 | rotate270)");
    operation.setDescription(TR("transpose the raster according to the method indicated by the second parameter"));
    operation.setInParameterCount({2});
    operation.addInParameter(0,itRASTER,  TR("input raster"),TR("ratser to be transposed"));
    operation.addInParameter(1,itSTRING, TR("transpose method"),TR("rotation or mirror of the input map"));
    operation.setOutParameterCount({1});
    operation.addOutParameter(0,itRASTER, TR("output raster"), TR("output raster with a new georef"));
    operation.setKeywords("raster, geometry");

    operation.checkAlternateDefinition();
    mastercatalog()->addItems({operation});
    return operation.id();
}
