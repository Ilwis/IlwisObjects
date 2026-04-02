
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
#include "symboltable.h"
#include "ilwisoperation.h"
#include "blockiterator.h"
#include "pixeliterator.h"
#include "MapFillSinks.h"

using namespace Ilwis;
using namespace Hydroflow;

REGISTER_OPERATION(MapFillSinks)

MapFillSinks::MapFillSinks()
{
}


MapFillSinks::MapFillSinks(quint64 metaid, const Ilwis::OperationExpression& expr) : OperationImplementation(metaid, expr)
{

}

bool MapFillSinks::execute(ExecutionContext* ctx, SymbolTable& symTable)
{
   if (_prepState == sNOTPREPARED)
        if ((_prepState = prepare(ctx, symTable)) != sPREPARED)
            return false;

    executeFillSink();

    bool resource = true;
    if (resource && ctx != 0) {
        QVariant value;
        value.setValue<IRasterCoverage>(_outRaster);
        ctx->setOutput(symTable, value, _outRaster->name(), itRASTER, _outRaster->resource());
    }
    return resource;
}

Ilwis::OperationImplementation* MapFillSinks::create(quint64 metaid, const Ilwis::OperationExpression& expr)
{
    return new MapFillSinks(metaid, expr);
}

Ilwis::OperationImplementation::State MapFillSinks::prepare(ExecutionContext* ctx, const SymbolTable& st)
{
	m_sinkPixels = 0;
    m_sinkHeight = -99999;
    m_sinkPixelsThreshold = 2;
	
    OperationImplementation::prepare(ctx, st);
    QString inraster = _expression.parm(0).value();
    QString outputName = _expression.parm(0, false).value();
    QString methodstr = _expression.parm(1).value().toLower();

    std::map<QString, FillMethod> methods = { {"fill",fmFill},{"cut",fmCut} };


    auto iter = methods.find(methodstr);
    if (iter == methods.end()) 
    {
        ERROR2(ERR_NOT_FOUND2, methodstr, TR("in method for fillsinks"));
        return sPREPAREFAILED;
    }
    method = iter->second;


    if (!_inRaster.prepare(inraster, itRASTER)) {
        ERROR2(ERR_COULD_NOT_LOAD_2, inraster, "");
        return sPREPAREFAILED;
    }
  
    
    int copylist = itRASTERSIZE | itENVELOPE | itCOORDSYSTEM | itGEOREF;
    _outRaster = OperationHelperRaster::initialize(_inRaster.as<IlwisObject>(), itRASTER, copylist);
    if (!_outRaster.isValid()) {
        ERROR1(ERR_NO_INITIALIZED_1, "output rastercoverage");
        return sPREPAREFAILED;
    }

    IDomain dom("code=domain:value");
    _outRaster->datadefRef() = DataDefinition(dom);

    for (quint32 i = 0; i < _outRaster->size().zsize(); ++i) {
        QString index = _outRaster->stackDefinition().index(i);
        _outRaster->setBandDefinition(index, DataDefinition(dom));
    }

    if (outputName != sUNDEF)
        _outRaster->name(outputName);

  
    _xsize = _outRaster->size().xsize();
    _ysize = _outRaster->size().ysize();

	
    return sPREPARED;
}

quint64 MapFillSinks::createMetadata()
{
    OperationResource operation({ "ilwis://operations/MapFillSinks" });
    operation.setSyntax("MapFillSinks(inputraster,method=fill|cut)");
    operation.setDescription(TR("Remove local depressions or cut terrain from a DEM raster according to the method indicated by the second parameter"));
    operation.setInParameterCount({ 2 });
    operation.addInParameter(0, itRASTER, TR("input raster DEM"), TR("input raster DEM with numeric domain"));
    operation.addInParameter(1, itSTRING, TR("fill method"), TR("fill sinks or cut terrain"), OperationResource::ueCOMBO);
    operation.parameterNeedsQuotes(1);
    operation.setOutParameterCount({ 1 });
    //operation.addValidation(0, 1, "values with select code from filters where type='linear' ");
    operation.addOutParameter(0, itRASTER, TR("output raster"), TR("output raster DEM that will contain height values without local depressions"));
    operation.setKeywords("fill sink,raster,image processing, numeric");

    operation.checkAlternateDefinition();
    mastercatalog()->addItems({ operation });
    return operation.id();
}


void MapFillSinks::executeFillSink()
{
    m_sinkPixelsThreshold = 2; //TO initiallize!!!

    /*iterFlag 1. Each element is initialled to zero;
                2. Elements in depres.contributing area are flagged,
                   when the the depres. area are defined;
                3. Elements in flat area flagged to -1*/
    m_iFlag = 0;

    m_inDem = new double[_xsize * _ysize];
    m_outDem = new double[_xsize * _ysize];
    m_flag = new int[_xsize * _ysize];

    PixelIterator iterDEM = PixelIterator(_inRaster, BoundingBox(), PixelIterator::fXYZ);
    PixelIterator iterOut = PixelIterator(_outRaster, BoundingBox(), PixelIterator::fXYZ);
    PixelIterator inEnd = iterDEM.end();

    while (iterDEM != inEnd)
    {
        Pixel px = iterDEM.position();
        m_inDem[idx(px.x, px.y)] = *iterDEM;
        m_outDem[idx(px.x, px.y)] = *iterDEM;
        m_flag[idx(px.x, px.y)] = 0;

        *iterOut = *iterDEM;
        *iterOut++;
        *iterDEM++;
    }
    SingleSinkFill();

    //Fill sinks based on the user specified threshold
    if (method == fmCut)
        GroupSinksFill();

    iterDEM = PixelIterator(_inRaster, BoundingBox(), PixelIterator::fXYZ);

    while (iterDEM != inEnd)
    {
        Pixel pxl = iterDEM.position();
        if (!onEdge(pxl) && fLocateInitialSink(pxl))
        {
            /*increment flag variable, locate the extent of the
                     contributing sink area*/
            m_iFlag++; //increment by the number of init. sink
            FindSinkContArea(pxl);

            //Identify outlet cell
            Pixel rcOutlet;
            if (fIdentifyOutletCell(pxl, rcOutlet))
            {
                /*if the elevation of outlet is greater than that of the
                 *initial sink then fill the depressions, otherwise
                 *if they are equal to, no filling is needed*/
                if (method == fmFill) 
                {
                    if (m_outDem[idx(rcOutlet.x, rcOutlet.y)] >
                        m_outDem[idx(pxl.x,pxl.y)] )
                    {
                        DepresFill(rcOutlet);
                    }
                    else
                        FlatAreaFlag(rcOutlet);
                }
                else
                    CutTerrain(rcOutlet);
            }
            else
                FlatAreaFlag(rcOutlet);
            _vPxlSinks.clear();

        }

        iterDEM++;
    }


    iterOut = PixelIterator(_outRaster, BoundingBox(), PixelIterator::fXYZ);

    while (iterOut != inEnd )
    {
        Pixel px = iterOut.position();
        *iterOut = m_outDem[idx(px.x, px.y)];
        ++iterOut;
    }

    _vPxlSinks.resize(0);
    delete []m_inDem;
	delete []m_outDem;
	delete []m_flag;

}


void MapFillSinks::CutTerrain(Pixel rcOutlet)
{
    // cell with lowest height to the outlet is selected for breaching
    double cutValue = getCutValue(rcOutlet);
    double rHeight = m_outDem[idx(rcOutlet.x, rcOutlet.y)];

    for (const Pixel& pxl : _vPxlSinks)
    {
        int i = idx(pxl.x, pxl.y);

        if (m_outDem[i] <= rHeight)
        {
            m_outDem[i] = cutValue;   // setPixelValue
            m_flag[i] = -1;        // flag the cell within the sink area
        }
    }
}


double MapFillSinks::getCutValue(Pixel rcOutlet)
{
    double rHeight = m_outDem[idx(rcOutlet.x, rcOutlet.y)];

    for (int ny = -1; ny <= 1; ny++)
    {
        for (int nx = -1; nx <= 1; nx++)
        {
            int xi = rcOutlet.x + nx;
            int yi = rcOutlet.y + ny;

            int i = idx(xi, yi);

            if (m_outDem[i] < rHeight)
                rHeight = m_outDem[i];
        }
    }

    return rHeight;
}


void MapFillSinks::FlatAreaFlag(Pixel rcOutlet)
{
    // flag the cells in an existing flat area
    for (const Pixel& pxl : _vPxlSinks)
    {
        m_flag[idx(pxl.x, pxl.y)] = -1;
    }
}


void MapFillSinks::DepresFill(Pixel rcOutlet)
{
    // elevation of outlet
    double rHeight = m_outDem[idx(rcOutlet.x, rcOutlet.y)];

    for (const Pixel& pxl : _vPxlSinks)
    {
        int i = idx(pxl.x, pxl.y);

        if (m_outDem[i] < rHeight)
            m_outDem[i] = rHeight;

        // flag the cell within the sink area
        m_flag[i] = -1;
    }
}


void MapFillSinks::SingleSinkFill()
{
    for (int y = 0; y < _ysize; y++)
    {
        for (int x = 0; x < _xsize; x++)
        {
            int i = idx(x, y);

            if (onEdge(Pixel(x, y, 0)))
            {
                m_flag[i] = -2;
            }
            else if (m_inDem[i] == rUNDEF)
            {
                FlagNeighbors(x, y);
            }
            else
            {
                double rMin = DBL_MAX;

                for (int ny = -1; ny <= 1; ny++)
                {
                    for (int nx = -1; nx <= 1; nx++)
                    {
                        if (nx == 0 && ny == 0)
                            continue;

                        int ni = idx(x + nx, y + ny);
                        rMin = std::min(rMin, m_outDem[ni]);
                    }
                }

                if (m_outDem[i] < rMin)
                    m_outDem[i] = rMin;
            }
        }
    }

}


void MapFillSinks::FlagNeighbors(int x, int y)
{
    for (int ny = -1; ny <= 1; ny++)
    {
        for (int nx = -1; nx <= 1; nx++)
        {
            int xi = x + nx;
            int yi = y + ny;

            if (xi >= 0 && xi < _xsize &&
                yi >= 0 && yi < _ysize)
            {
                m_flag[idx(xi, yi)] = -3;
            }
        }
    }
}
 

bool MapFillSinks::fLocateInitialSink(Pixel pxl)
{
    int x = pxl.x;
    int y = pxl.y;

    int i = idx(x, y);

    // skip if already flagged
    if (m_flag[i] != 0)
        return false;

    double center = m_outDem[i];
    double rMin = DBL_MAX;

    // check 8 neighbors
    for (int ny = -1; ny <= 1; ny++)
    {
        for (int nx = -1; nx <= 1; nx++)
        {
            if (nx == 0 && ny == 0)
                continue;

            int xi = x + nx;
            int yi = y + ny;

            int ni = idx(xi, yi);

            double val = m_outDem[ni];

            if (val < rMin)
                rMin = val;
        }
    }

    return (center <= rMin);
}


double MapFillSinks::getPixelValue(Pixel pxl)
{
    return *(iterOut(pxl));

}

void MapFillSinks::FindSinkContArea(Pixel rcInitSink)
{
    std::vector<Pixel> vStartCells;
    std::vector<Pixel> vAdjaCells;

    // initial expansion
    FlagAdjaCell(rcInitSink, vStartCells);

    m_flag[idx(rcInitSink.x, rcInitSink.y)] = m_iFlag;

    _vPxlSinks.clear();
    _vPxlSinks.push_back(rcInitSink);

    while (!vStartCells.empty())
    {
        vAdjaCells.clear();

        for (const Pixel& pxl : vStartCells)
        {
            _vPxlSinks.push_back(pxl);

            // expand neighbors using same logic
            FlagAdjaCell(pxl, vAdjaCells);
        }

        vStartCells.swap(vAdjaCells);
    }
}


void MapFillSinks::FindSinkContArea2(Pixel rcInitSink)
{
    std::vector<Pixel> vStartCells;
    std::vector<Pixel> vAdjaCells;

    m_sinkPixels = 1;

    // initial expansion
    FlagAdjaCell(rcInitSink, vStartCells);

    m_flag[idx(rcInitSink.x, rcInitSink.y)] = m_iFlag;

    _vPxlSinks.clear();
    _vPxlSinks.push_back(rcInitSink);

    while (!vStartCells.empty())
    {
        vAdjaCells.clear();

        for (const Pixel& pxl : vStartCells)
        {
            _vPxlSinks.push_back(pxl);

            // expand neighbors using SAME logic
            FlagAdjaCell(pxl, vAdjaCells);
        }

        vStartCells.swap(vAdjaCells);

        if (m_sinkPixels > m_sinkPixelsThreshold)
            break;
    }
}


void MapFillSinks::FlagAdjaCell(Pixel rcStartCell, std::vector<Pixel>& vAdj)
{
    int sx = rcStartCell.x;
    int sy = rcStartCell.y;
    int si = idx(sx, sy);

    for (int i = -1; i <= 1; i++)
    {
        for (int j = -1; j <= 1; j++)
        {
            int x = sx + j;
            int y = sy + i;

            if (x < 0 || x >= _xsize ||
                y < 0 || y >= _ysize)
                continue;

            int ai = idx(x, y);

            if (m_flag[ai] != m_iFlag &&
                m_flag[ai] > -2 &&
                m_outDem[ai] >= m_outDem[si])
            {
                Pixel adj(x, y, 0);

                vAdj.push_back(adj);
                m_flag[ai] = m_iFlag;

                if (m_outDem[ai] == m_sinkHeight)
                    m_sinkPixels++;
            }
        }
    }
}


bool MapFillSinks::fIdentifyOutletCell(Pixel rcSink, Pixel& rcOutlet)
{
    // Find outlet cell in the rim of the sink contributing area
    std::vector<Pixel> vOutlets;
    vOutlets.clear();

    for (const Pixel& pxl : _vPxlSinks)
    {
        if (IsPotentialOutlet(pxl))
            vOutlets.push_back(pxl);
    }

    if (!vOutlets.empty())
    {
        // select outlet with minimum elevation
        Pixel best = vOutlets[0];
        double minVal = m_outDem[idx(best.x, best.y)];

        for (size_t i = 1; i < vOutlets.size(); ++i)
        {
            Pixel p = vOutlets[i];
            double val = m_outDem[idx(p.x, p.y)];

            if (val < minVal)
            {
                minVal = val;
                best = p;
            }
        }

        rcOutlet = best;
        return true;
    }

    return false;
}


bool MapFillSinks::IsPotentialOutlet(Pixel pxl)
{
    double rHeight = m_outDem[idx(pxl.x, pxl.y)];

    for (int i = -1; i <= 1; i++)
    {
        for (int j = -1; j <= 1; j++)
        {
            int x = pxl.x + j;
            int y = pxl.y + i;

            if (x < 0 || x >= _xsize || y < 0 || y >= _ysize)
                continue;

            int ni = idx(x, y);

            if (m_flag[ni] != m_iFlag)
            {
                if (m_flag[ni] == -2)   // flat area at DEM edge
                {
                    if (rHeight >= m_outDem[ni])
                        return true;
                }
                else
                {
                    if (rHeight > m_outDem[ni])
                        return true;
                }
            }
        }
    }

    return false;
}


bool MapFillSinks::onEdge(Pixel pix) {
    return pix.x == 0 || pix.x == _xsize - 1 ||
        pix.y == 0 || pix.y == _ysize - 1;
}


void MapFillSinks::GroupSinksFill()
{
    // scan DEM seeking an initial sink

    for (int y = 0; y < _ysize; y++)
    {
        for (int x = 0; x < _xsize; x++)
        {
            Pixel rcSink(x, y, 0);

            if (!onEdge(rcSink) && fLocateInitialSink(rcSink))
            {
                m_sinkHeight = m_outDem[idx(x, y)];

                // increment flag variable
                m_iFlag++;

                // locate contributing area
                FindSinkContArea2(rcSink);

                if (m_sinkPixels <= m_sinkPixelsThreshold)
                {
                    // Identify outlet cell
                    Pixel rcOutlet;

                    if (fIdentifyOutletCell(rcSink, rcOutlet))
                    {
                        int sinkIndex = idx(rcSink.x, rcSink.y);
                        int outletIndex = idx(rcOutlet.x, rcOutlet.y);

                        if (m_outDem[outletIndex] > m_outDem[sinkIndex])
                            DepresFill(rcOutlet);
                        else
                            FlatAreaFlag(rcOutlet);
                    }
                    else
                    {
                        FlatAreaFlag(rcOutlet);
                    }

                    _vPxlSinks.clear();
                }
            }
        }
    }

    // reset positive flags to zero
    for (int y = 0; y < _ysize; y++)
    {
        for (int x = 0; x < _xsize; x++)
        {
            int i = idx(x, y);

            if (m_flag[i] > 0)
                m_flag[i] = 0;
        }
    }
}

