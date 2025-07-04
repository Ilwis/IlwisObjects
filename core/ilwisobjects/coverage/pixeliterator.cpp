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
#include "geos/geom/PrecisionModel.h"
#include "geos/algorithm/locate/SimplePointInAreaLocator.h"
#include "geos/geom/Point.h"
#ifdef Q_OS_WIN
#include "geos/geom/Envelope.h"
//#include "geos/geom/PrecisionModel.h" // already included in featurecoverage.cpp
#endif
#include "geos/geom/GeometryFactory.h"
#include "tranquilizer.h"
#include "table.h"
#include "feature.h"
#include "pixeliterator.h"
#include "bresenham.h"
#include "vertexiterator.h"
#include "operationhelpergrid.h"

using namespace Ilwis;

PixelIterator::PixelIterator() :
    _grid(0),
    _x(0),
    _y(0),
    _z(0),
    _flow(fXYZ),
    _isValid(true),
    _endx(0),
    _endy(0),
    _endz(0),
    _linearposition(0),
    _endposition(0),
    _xChanged(false),
    _yChanged(false),
    _zChanged(false)
{

}

PixelIterator::PixelIterator(const IRasterCoverage& raster, geos::geom::Geometry* selection) :
    _raster(raster),
    _flow(fXYZ),
    _isValid(false)
{
    //geos::geom::GeometryTypeId type = selection->getGeometryTypeId();
    const geos::geom::Envelope *env = selection->getEnvelopeInternal();
    Coordinate crdMin = Coordinate(env->getMinX(), env->getMinY());
    Coordinate crdMax = Coordinate(env->getMaxX(), env->getMaxY());
    _box = raster->georeference()->coord2Pixel(Box<Coordinate>(crdMin, crdMax));
    _box.min_corner().x = std::max(_box.min_corner().x, 0);
    _box.min_corner().y = std::max(_box.min_corner().y, 0);
    _box.max_corner().x = std::min(_box.max_corner().x, (int)(raster->size().xsize() - 1));
    _box.max_corner().y = std::min(_box.max_corner().y, (int)(raster->size().ysize() - 1));
    init();
    Bresenham algo(_raster->georeference());
    VertexIterator iter(selection);

    std::vector<Pixel> selectionPix = algo.rasterize(::begin(iter), ::end(iter));
    if ( selectionPix.size() == 0){
        _isValid = false;
        return;
    }
    _selectionPixels.resize(_raster->size().ysize());
    // we use  a vector of sets to store the final cleaned selection as sets are both ordered and have no duplicates,
    if ( selection->getGeometryTypeId() == geos::geom::GEOS_POLYGON || selection->getGeometryTypeId() == geos::geom::GEOS_MULTIPOLYGON)
        cleanUp4PolyBoundaries(selectionPix, selection);
    else{
        for(qint64 i = 0; i < selectionPix.size(); ++i){
            _selectionPixels[selectionPix[i].y].push_back(selectionPix[i].x) ;
        }
        _x = selectionPix[0].x;
        _y = selectionPix[0].y;
    }
    initPosition();
    _selectionIndex = 1;
    _insideSelection = true;
}

PixelIterator::PixelIterator(const IRasterCoverage &raster, const BoundingBox& box, Flow flow) :
    _raster(raster),
    _box(box),
    _flow(flow),
    _isValid(false)
{
    init();
}

PixelIterator::PixelIterator(const IRasterCoverage& raster, int threadIdx, const BoundingBox& box, Flow flow)  :
	_raster(raster),
	_box(box),
	_flow(flow),
	_isValid(false),
	_threadIndex(threadIdx)
{
	init();
}

PixelIterator::PixelIterator(const IRasterCoverage& raster, int threadIdx, const ProcessingBoundingBoxes& boxes, Flow flow)  :
    _raster(raster),
    _box(boxes.box(raster, threadIdx)),
    _flow(flow),
    _isValid(false),
    _threadIndex(threadIdx)
{
    init();
}

PixelIterator::PixelIterator(const IRasterCoverage &raster, PixelIterator::Flow flow) :
    _raster(raster),
    _box(BoundingBox()),
    _flow(flow),
    _isValid(false)
{
    init();
}

PixelIterator::PixelIterator(PixelIterator&& iter) :
    _raster(std::move(iter._raster)),
    _grid(std::move(iter._grid)),
    _box(std::move(iter._box)),
    _x(iter._x),
    _y(iter._y),
    _z(iter._z),
    _flow(iter._flow),
    _isValid(iter._isValid),
    _endx(iter._endx),
    _endy(iter._endy),
    _endz(iter._endz),
    _linearposition(iter._linearposition),
    _endposition(iter._endposition),
    _xChanged(iter._xChanged),
    _yChanged(iter._yChanged),
    _zChanged(iter._zChanged),
    _selectionPixels(iter._selectionPixels),
    _selectionIndex(iter._selectionIndex),
    _insideSelection (iter._insideSelection)
{
}

PixelIterator::PixelIterator(const PixelIterator& iter)  {
    copy(iter);
}

PixelIterator::~PixelIterator(){
    if ( _raster.isValid()){
        _raster->gridRef()->closure();
    }
}

void PixelIterator::copy(const PixelIterator &iter) {
    _raster = iter._raster;
    if ( _raster.isValid())
        _grid = _raster->_grid.get();
    _box = iter._box;
    _isValid = iter._isValid;
    _flow = iter._flow;
    _x = iter._x;
    _y = iter._y;
    _z = iter._z;
    _endx = iter._endx;
    _endy = iter._endy;
    _endz = iter._endz;
    _linearposition = iter._linearposition;
    _endposition = iter._endposition;
    _selectionPixels  = iter._selectionPixels;
    _selectionIndex = iter._selectionIndex;
    _insideSelection = iter._insideSelection;

}

void PixelIterator::init() {
    const Size<>& sz = _raster->size();
    if ( !_box.isValid()) {
        _box = BoundingBox(sz);
    }
    _box.ensure(sz);

    _x = _box.min_corner().x;
    _y = _box.min_corner().y;
    _z = _box.min_corner().z;
    if ( isNumericalUndef(_z))
        _z = 0;

    _endx = _box.max_corner().x;
    _endy = _box.max_corner().y;
    _endz = _box.max_corner().z;
    if ( isNumericalUndef(_endz))
        _endz = 0;

    _grid = _raster->gridRef().get();
    if ( _grid == 0) {
        _isValid = false;
        if (!isValid())
            throw ErrorObject(TR("Using invalid pixeliterator, are all data sources accessible?"));
    }

    initPosition();
    bool inside = contains(Pixel(_x,_y, _z));
    _isValid = inside;
    _xChanged = _yChanged = _zChanged = false;
    _grid->setOrientation(_flow == Flow::fZXY ? Grid::Orientation::oZXY : Grid::Orientation::oXYZ);
}

bool PixelIterator::moveXY(qint64 delta){
    _zChanged = (_z - delta) %  (qint64)_box.zlength() != 0;
    qint64 tempx = _x + (_z - _box.min_corner().z) / _box.zlength();
    _z = _box.min_corner().z + (_z - _box.min_corner().z) % (qint64)_box.zlength();
    _xChanged = tempx != _x;
    std::swap(_x,tempx);
    if ( _x > _endx) {
        quint64 newy = _y + (_x - _box.min_corner().x) / _box.xlength();
        _yChanged = newy != _y;
        _y = newy;
        _x = _box.min_corner().x + (_x - _box.min_corner().x) % (qint64)_box.xlength();
        _xChanged = _x != tempx;
        if ( _y > _endy) { // done with this iteration block
            _linearposition = _endposition;
            return false;
        }
    }
    _linearposition = _x + _y * _grid->size().xsize() + _z * _grid->size().xsize() * _grid->size().ysize();
    return true;
}

bool PixelIterator::moveXZ(qint64 delta)
{
    qint64 tempy = _y;
    qint64 tempx = _x;
    _y = _box.min_corner().y + ( _y - _box.min_corner().y) % (qint64)_box.ylength();
    _x += (tempy - _y) / _box.ylength();
    _xChanged = tempx != _x;
    _linearposition = _x + _y * _grid->size().xsize() + _z * _grid->size().xsize() * _grid->size().ysize();

    if ( _x > _endx){
        quint64 newz = _z + (_x - _box.min_corner().x) / _box.xlength();
        _zChanged = newz != _z;
        _z = newz;
        tempx = _x;
        _x = _box.min_corner().x + (_x - _box.min_corner().x) % (qint64)_box.xlength();
        _xChanged = _x != tempx;
        if ( _z > _endz) { // done with this iteration block
            _linearposition = _endposition;
            return false;
        }
    }
    return true;
}


bool PixelIterator::moveYZ(qint64 delta){
    qint64 tempy;
    if ( _x < _box.min_corner().x){

        int xdelta = _box.min_corner().x - _x;
        _x += _box.xlength() * ceil((_box.min_corner().x - _x) / _box.xlength());
        if ( _endposition == _linearposition - delta) { // special case, we were at one beyond the end;
            tempy = _box.max_corner().y;
            --_z;
            _linearposition = _z * _box.size().linearSize() + _box.max_corner().y * _grid->size().xsize() + _box.max_corner().x + delta + 1; // +1 because we were at 1 beyound the end.
        } else
            tempy = _y - ceil(xdelta / _box.xlength()) + (_x - _box.min_corner().x) / _box.xlength();

        _xChanged = delta %  (qint64)_box.xlength() != 0;
    }
    else{
        tempy = _y + (_x - _box.min_corner().x) / _box.xlength(); // with a big delta, y might jump more number of lines == xdistance / box xsize
        _xChanged = (_x - delta) %  (qint64)_box.xlength() != 0;
        _x = _box.min_corner().x + (_x - _box.min_corner().x) % (qint64)_box.xlength(); //  the modulo xlength gives the jump in one line
    }

    _yChanged = tempy != _y;
    std::swap(_y,tempy);

    if ( _y > _endy) {
        quint64 newz = _z + (_y - _box.min_corner().y) / _box.ylength();
        _zChanged = newz != _z;
        _z = newz;
        _y = _box.min_corner().y + (_y - _box.min_corner().y) % (qint64)_box.ylength();
        _yChanged = _y != tempy;
        if ( _z > _endz) { // done with this iteration block
            _linearposition = _endposition;
            return false;
        }
    }else if ( _y < _box.min_corner().y){
        _y = _box.max_corner().y;
        --_z;
        _zChanged = true;
        if ( _z < 0) {
            return false;
        }

    }
    return true;
}

Pixel PixelIterator::position() const
{
    return Pixel(_x, _y, _z);
}

void PixelIterator::position(const Pixel& pos) {
    _xChanged = pos.x != _x;
    _yChanged = pos.y != _y;
    _zChanged = pos.z != _z;
    if ( _xChanged || _yChanged || _zChanged){
        _x = pos.x;
        _z = pos.z;
        _y = pos.y;
        _linearposition = (_endx + 1) * (_endy + 1) * _z + _y * (_endx + 1) + _x;
        if (_linearposition >= _endposition || _linearposition < 0){
            _isValid = false;
            _x = _endx;
            _y = _endy;
            _z = _endz;
        }
    }
}

const BoundingBox &PixelIterator::box() const
{
    return _box;
}

void PixelIterator::box(const BoundingBox& box) {
	_box = box;
	init();
}

quint64 PixelIterator::linearPosition() const
{
   // return sz.xsize() * sz.ysize() * _z + sz.xsize() * _y + _x;
    return _linearposition;
}

void PixelIterator::setRaster(const IRasterCoverage &raster) {
    const BoundingBox box;

    _raster = raster;
    _box = box;
    _flow = fXYZ;
    _isValid = false;

    init();
}

PixelIterator& PixelIterator::operator=(const PixelIterator& iter) {
    copy(iter);
    return *this;
}

PixelIterator& PixelIterator::operator=(const PixelIterator&& iter) {
    copy(iter);
    return *this;
}

QVariant PixelIterator::operator()(const QString &column)
{
    quint32 raw = operator *();
    return _raster->attributeTable()->cell(column, raw);
}

inline PixelIterator PixelIterator::operator++(int) {
    PixelIterator temp(*this);
    if(!move(1))
        return end();
    return temp;
}
inline PixelIterator PixelIterator::operator--(int) {
    PixelIterator temp(*this);
    move(-1);
    return temp;
}

bool PixelIterator::operator==(const PixelIterator& iter) const{
    return _linearposition == iter._linearposition;
}

bool PixelIterator::operator!=(const PixelIterator& iter) const{
    return ! operator ==(iter);
}

bool PixelIterator::operator<(const PixelIterator &iter) const
{
    return _linearposition < iter._linearposition;
}

bool PixelIterator::operator>(const PixelIterator &iter) const
{
    return _linearposition > iter._linearposition;
}

bool PixelIterator::operator<=(const PixelIterator &iter) const
{
    return _linearposition <= iter._linearposition;
}

bool PixelIterator::operator>=(const PixelIterator &iter) const
{
    return _linearposition >= iter._linearposition;
}

PixelIterator PixelIterator::end() const {
    PixelIterator iter(*this);
	iter += _box.size().linearSize();
    return iter;
}

void PixelIterator::toEnd()
{
    _x = _endx;
    _z = _endz;
    _y = _endy;
    _linearposition = (_endx + 1) * (_endy + 1) * (_endz + 1)  -1;
}

void PixelIterator::setFlow(Flow flw) {
    _flow = flw;
}

bool PixelIterator::contains(const Pixel& pix) {
    bool ok =  pix.x >= _box.min_corner().x &&
            pix.x < _box.max_corner().x  &&
            pix.y >= _box.min_corner().y &&
            pix.y < _box.max_corner().y;
    if (isNumericalUndef(pix.z)){
        if( _box.min_corner().z == 0 && _box.max_corner().z == 0)
            return ok;
        if ( isNumericalUndef(_box.min_corner().z) && isNumericalUndef(_box.max_corner().z))
             return ok;
        return false;
    }
    if ( _box.is3D())
        ok = ok && pix.z  >= _box.min_corner().z &&  pix.z <= _box.max_corner().z;
    else
        ok = ok && ( pix.z == 0 || isNumericalUndef(pix.z));
    return ok;
}

bool PixelIterator::xchanged() const {
    return _xChanged;
}

bool PixelIterator::ychanged() const {
    return _yChanged;
}

bool PixelIterator::zchanged() const {
    return _zChanged;
}

void PixelIterator::initPosition() {
    const Size<>& sz = _raster->size();
    quint64 linpos = _y * sz.xsize() + _x;
    quint64 endpos = _endy * sz.xsize() + _endx;
    _linearposition = sz.xsize() * sz.ysize() * _z + linpos;
    _endposition = sz.xsize() * sz.ysize() * _endz + endpos + 1; // one past the last valid position
}

qint64 PixelIterator:: operator-(const PixelIterator& iter) {
    return linearPosition() - iter.linearPosition();
}

const IRasterCoverage& PixelIterator::raster() const
{
    return _raster;
}


bool PixelIterator::move2NextSelection(qint64 delta)
{
    if ( _selectionIndex  >= _selectionPixels[_y].size() - 1){ // this was the last boundary on this row
        _x = _endx + 1; // put x beyond the edge of the box so moveXY will trigger a y shift
        if(!moveYZ(delta))
            return false;
        while (_y < _selectionPixels.size() && _selectionPixels[_y].size() == 0) { // search for the next non-empty row
            _x = _endx + 1; // put x beyond the edge of the box so moveXY will trigger a y shift
            if(!moveYZ(delta))
                return false;
        }
        if ( _y >= _selectionPixels.size())
            return false;
        qint64 xnew = _selectionPixels[_y][0];
        _linearposition += xnew - _box.min_corner().x;
        _selectionIndex = 0;
        _insideSelection = false; //  we are not yet in a selection
        if ( _selectionPixels[_y].size() > 0){
            _x = _selectionPixels[_y][0] - delta; // -delta because the next iteration will add delta again to it setting it exact at the edge again

        }
    }else{
        qint64 xnew = _selectionPixels[_y][++_selectionIndex];
        _linearposition += xnew - _box.min_corner().x;
        _x = xnew;
        moveYZ(delta);
        _x -= 1;
    }
    return true;
}

void PixelIterator::cleanUp4PolyBoundaries(const std::vector<Pixel>& selectionPix, geos::geom::Geometry* selection)
{
    std::vector<std::set<qint32>> boundaries(_selectionPixels.size());
    for(const Pixel& pixel : selectionPix){
        if ( (pixel.y >= 0) && (pixel.y < _selectionPixels.size())){
            boundaries[pixel.y].insert(pixel.x);
        }
    }
    int y = 0;
    int ystart = iUNDEF;
    // getting rid of all unneeded pixels. adjacent pixels are not needed as you can't rasterize between them
    geos::geom::PrecisionModel *precisionModel = new geos::geom::PrecisionModel(geos::geom::PrecisionModel::FLOATING);
    geos::geom::GeometryFactory *geometryFactory = new geos::geom::GeometryFactory(precisionModel,-1);

    for(std::set<qint32> &borders : boundaries){
        borders.insert(_box.min_corner().x);
        borders.insert(_box.max_corner().x);
        std::vector<qint32> cleaned;
        std::vector<std::pair<qint32,qint32>> pairs;
        qint32 prev = iUNDEF;
        qint32 run = iUNDEF;
        for (qint32 xpos : borders) {
            if (run == iUNDEF)
                run = xpos;
            if (prev != iUNDEF) {
                if (xpos - prev != 1) {
                    if (run != prev) {
                        // if a run was detected save it
                        std::pair<qint32,qint32> p1(run, prev);
                        pairs.push_back(p1);
                    }
                    // regular case: save pair
                    std::pair<qint32,qint32> p(prev, xpos);
                    pairs.push_back(p);
                    run = xpos;
                } else if (xpos == _box.max_corner().x) {   // special case
                    // happens when we have a run up until and including the right edge
                    std::pair<qint32,qint32> p(run, xpos);
                    pairs.push_back(p);
                }
            }
            prev = xpos;
        }
        for (std::pair<qint32,qint32> p : pairs) {
            PIXVALUETYPE middle = (p.first + p.second - 1) / 2.0;
            Pixeld position (middle + 0.5, y + 0.5); // inspect relationship with the given geometry in the middle between two borders
            Coordinate crd = _raster->georeference()->pixel2Coord(position);
            geos::geom::Point *pnt = geometryFactory->createPoint(crd);
            if (selection->contains(pnt)) {
                if ((cleaned.size() > 0) && (cleaned[cleaned.size() - 1]) == p.first) {
                    cleaned[cleaned.size() - 1] = p.second;
                } else {
                    cleaned.push_back(p.first);
                    cleaned.push_back(p.second);
                }
            }
            delete pnt;
        }

        _selectionPixels[y].resize(cleaned.size());
        std::copy(cleaned.begin(), cleaned.end(), _selectionPixels[y].begin() ); //copy pixels to output vector
        if (ystart == iUNDEF && cleaned.size() != 0) {
            ystart = y;
        }
        ++y;
    }
    delete geometryFactory;
    delete precisionModel;
    if ( ystart != iUNDEF){
        _y = ystart;
        _x =  _selectionPixels[_y][0];
    }
}
