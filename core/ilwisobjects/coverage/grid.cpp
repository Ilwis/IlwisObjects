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
#include "ilwiscontext.h"
#include "connectorinterface.h"
#include "geometries.h"
#include "grid.h"

using namespace Ilwis;

/*
 * The planes member, a 3D nested set of vectors, contains the 'real' data. It is structured differently for the XYZ flow and the ZXY flow
 * the XYZ flow is structured planes[zsize][ysize][xsize], the ZXY flow planes[ysize][zsize][xsize]. The amount of data actually loaded in
 * these vectors is limited. For XYZ a block like planes[ZIndex][YLower ... YUpper][0..XMAX]. The moment y passes YUpper,
 * planes[ZIndex][YLower ... YUpper][0..XMAX] will be written as block to the cache file and the memory occuppied by that block is freed.
 * Whenever it is needed it can be very quickly reloaded as it is a simple dump of memory
 */
//----------------------------------------------------------------------

Grid::Grid()   {
    //setOrientation(ori);
}


Grid::~Grid() {
    clear();

}

Size<> Grid::size() const {
    if ( _igrid)
        return _igrid->size();
    return Size<>();
}

Grid *Grid::clone(quint64 newRasterId, int dLimit1, int dLimit2)
{
    Locker<> lock(_mutex);

    Grid *grid = new Grid();
    grid->setOrientation(_igrid->orientation());
    grid->useCache(_igrid->useCache());
    _igrid->clone(newRasterId, dLimit1, dLimit2, grid->igrid().get());

    grid->closure();
    return grid;

}

void Grid::clear() {
    Locker<> lock(_mutex);
    _igrid->clear();
}

void Grid::dump(unsigned int y1, unsigned int z1, unsigned int z2, int threadIndex) {
    if ( _igrid)
        _igrid->dump(y1,z1,z2,threadIndex);
}

void Grid::load(unsigned int y1, unsigned int z1, unsigned int z2, unsigned int x, int threadIndex){
    if ( _igrid) 
        _igrid->load(y1,z1,z2,x,threadIndex);
}

bool inside(int x, int y, int z, const BoundingBox& box){

    const auto& pmin = box.min_corner();
    const auto& pmax = box.max_corner();
    return x >= pmin.x && x <= pmax.x && y >= pmin.y && y <= pmax.y && z >= pmin.z && z <= pmax.z;
}
PIXVALUETYPE& Grid::valueRef(const Pixel& pix, int threadIndex) {
    quint32 z = pix.z == iUNDEF ? 0 : pix.z;
#ifdef QT_DEBUG
    if ( inside(pix.x, pix.y, pix.z, _igrid->boundingBox())){
        return (*_igrid.*(_igrid->_valueFuncPtr))(z, pix.y, pix.x, threadIndex);
    }
          
    throw ErrorObject("Pixel location out of bounds");
#else
    return (*_igrid.*(_igrid->_valueFuncPtr))(z, pix.y, pix.x, threadIndex);
#endif
}

PIXVALUETYPE Grid::value(const Pixel& pix, int threadIndex) {
   quint32 z = pix.z == iUNDEF ? 0 : pix.z;
#ifdef QT_DEBUG
    if ( inside(pix.x, pix.y, pix.z, _igrid->size()))
        return (*_igrid.*(_igrid->_valueFuncPtr))(z, pix.y, pix.x, threadIndex);
    else
       return PIXVALUEUNDEF;
#else
    return (*_igrid.*(_igrid->_valueFuncPtr))(z, pix.y, pix.x, threadIndex);
#endif
}

PIXVALUETYPE &Grid::value(int z, int y, int x, int threadIndex)  {
    z = z == iUNDEF ? 0 : z;
#ifdef QT_DEBUG
    if ( inside(x, y, z, _igrid->boundingBox())){
        return (*_igrid.*(_igrid->_valueFuncPtr))(z, y, x, threadIndex);
    }
    throw ErrorObject("Pixel location out of bounds");
#else
    return (*_igrid.*(_igrid->_valueFuncPtr))(z, y, x, threadIndex);
#endif
}

void Grid::setBlockData(int z, int y, const std::vector<PIXVALUETYPE>& data){
    if ( _igrid)
        _igrid->setBlockData(z,y,data);
}

void Grid::setValue(int z, int y, int x, PIXVALUETYPE v, int threadIndex ){

    if ( _igrid)
        _igrid->setValue(z,y,x,v,threadIndex);
}

bool Grid::prepare(quint64 rasterid, const Size<> &sz,Grid::Orientation ori, const BoundingBox& box) {
    Locker<> lock(_mutex);
    if ( !_igrid || _igrid->orientation() != ori){
        if ( ori == oXYZ){
            _igrid.reset(new GridXYZ());
        }else
            _igrid.reset(new GridZXY());
    }
    if ( !_igrid->isDefined())
        _igrid->prepare(rasterid, sz, ori, box);
    return true;
}


void Grid::setOrientation(Grid::Orientation ori){
    if ( !_igrid || _igrid->orientation() != ori){
        if ( ori == oXYZ){
            _igrid.reset(new GridXYZ());
        }else
            _igrid.reset(new GridZXY());

    }
    _igrid->setOrientation(ori);
}

void Grid::unload(int threadIndex) {
    Locker<> lock(_mutex);
    if ( _igrid)
        _igrid->unload(threadIndex);
}

bool Grid::isValid() const
{
    if (_igrid)
        return !(_igrid->size().isNull() || _igrid->size().isValid());
    return false;
}

void Grid::prepare4Operation(int nThreads) {
    if ( _igrid)
        _igrid->prepare4Operation(nThreads);
}

void Grid::unprepare4Operation() {

}

quint32 Grid::blockCacheLimit() const
{
    if ( _igrid){
        return _igrid->orientation() == Grid::oXYZ ? 100 : 8;
    }
}

quint32 Grid::blocksCount() const{
    if ( _igrid)
        return _igrid->blocksCount();

    return iUNDEF;
}

quint64 Grid::blockSize(int d) const {
    if ( _igrid)
        return _igrid->blockSize(d);

    return i64UNDEF;
}


void Grid::setd3Size(int z, int y, int xzs){
    if ( _igrid)
        _igrid->set3DSize(z,y,xzs);
}


quint32 Grid::linesPerBlock(int d) const
{
    if ( _igrid)
        return _igrid->linesPerBlock(d);

    return iUNDEF;
}

quint64 Grid::seekPosition(int z, int y, int x) const{
    if ( _igrid)
        return _igrid->seekPosition(z,y,x);
    return i64UNDEF;
}

void Grid::closure()
{
    if ( _igrid)
        _igrid->closure();
}


void Grid::useCache(bool yesno){
   if ( _igrid)
       _igrid->useCache(yesno);
}

Grid::Orientation Grid::flow() const{    void dump(unsigned int y1, unsigned int y2, unsigned int z1, unsigned int z2, int threadIndex);
     if ( _igrid)
         return _igrid->flow();
     return oUNKNOWN;
}


//------------------------------------------------------------------------------------
InternalGrid::InternalGrid()
{
    
}

InternalGrid::InternalGrid(Grid::Orientation ori) : _orientation(ori), _rasterid(i64UNDEF)  {
}


void InternalGrid::sizeData(unsigned int y2, unsigned int z2, unsigned int& ysize, unsigned int& zsize, std::vector<char>& data)
{
    quint64 dataSize = _size.xsize() * sizeof(PIXVALUETYPE);
    zsize = std::min(_size.zsize(), z2);
    ysize = std::min(_size.ysize(),y2);
    data.resize(dataSize);
}

bool InternalGrid::prepare(quint64 rasterid, const Size<> &sz,Grid::Orientation ori, const BoundingBox& box) {
    clear();
    _rasterid = rasterid;
    _size = sz;

    if (box.isValid()){
        _box =  box;
    } else
        _box = BoundingBox(sz);
    setOrientation(ori);
    return true;
}

void InternalGrid::dumpBlock(int z, int yBase, int threadIndex)
{
    if (_diskImages.size() <= (size_t)threadIndex && _useCache){
        _diskImages.resize(threadIndex+1);
        _imageNames.resize(threadIndex + 1);
        createCacheFile(threadIndex);
    }
    if ( _orientation == Grid::oXYZ && blockStatus(z, yBase)._valid){
        dump(yBase , z, z + 1,threadIndex);
    }else if ( blockStatus(0, yBase)._valid){
        dump(yBase, 0, _size.zsize()-1, threadIndex);
    }
}

void InternalGrid::loadFromCache(unsigned int yBase, unsigned int z1, unsigned int z2, unsigned int x, int threadIndex)
{
    quint64 dataSize = _size.xsize() * sizeof(PIXVALUETYPE);
    std::vector<char> data(dataSize);
    for (unsigned int z = z1; z <= std::min(_size.zsize()-1,z2); ++z){
        for(unsigned int y = yBase; y <= std::min(_size.ysize()-1,yBase + XYZYBLOCKS); ++y){
            quint64 seekposition = calcSeekPosition(z, y);
            _diskImages[threadIndex]->seek(seekposition);
            _diskImages[threadIndex]->read(&data[0], dataSize);
            loadFromCache(y,z, data);
            setStartIndexes();

        }
    }
    setBlockStatus(z1, yBase, true);
}

quint64 InternalGrid::calcSeekPosition(unsigned int z, unsigned int y ) const{
   quint64 pos = (z * _size.xsize() * _size.ysize() + y * _size.xsize()) * sizeof(PIXVALUETYPE) + 3 * sizeof(int);
   return pos;
}

void InternalGrid::loadFromSource(int z1, int x, int y1){
    IIlwisObject obj = mastercatalog()->get(_rasterid);
    if ( obj.isValid() ){
        IRasterCoverage raster = obj.as<RasterCoverage>();
        IOOptions options;
        quint64 blockFloor = y1; //blockNr * blockCacheLimit() ;
        options.addOption({"z", z1});
        options.addOption({"bands", _size.zsize()});
        options.addOption({"y", blockFloor});
        options.addOption({"size", blockCacheLimit() * _size.xsize()});
        options.addOption({"orientation", _orientation == Grid::oZXY ? "ZXY" : "XYZ"});
        raster->getData(options);
    }
}

bool InternalGrid::pastHorizon(quint32 v) const {
    return v % blockCacheLimit() == 0 && v != 0;
}

void InternalGrid::load(unsigned int y1, unsigned int z1, unsigned int z2, unsigned int x, int threadIndex)
{
    quint64 ybase = y1;
    IIlwisObject obj = mastercatalog()->get(_rasterid);
    if (obj.isValid() || _imageNames.size() > 0 ) { //.is used in the clone function; an emptry grid is created without(yet) a raster belonging to it
        if (blockStatus(z1,ybase)._loadedFromSource && _useCache){
                loadFromCache(y1,z1, z2,  x, threadIndex);
            }else{
                loadFromSource(z1,x,y1);
            }
    }else {
       setBlockStatus(z1, ybase, true);
    }


    if ( pastHorizon(ybase)){
        dumpBlock(z1,ybase,threadIndex);
    }
}

void InternalGrid::clear()
{
    _validStripe = RelativeVector<RelativeVector<BlockStatus>>();
    _planes = RelativeVector<Plane>();
    if ( _diskImages.size() > 0){
        for(auto *file : _diskImages){
            file->close();
            file->remove();
        }
        _diskImages = std::vector<QFile*>();
    }
    _size = Size<>();
    _rasterid = iUNDEF;
    _orientation = Grid::oXYZ;
}

bool InternalGrid::useCache() const
{
    return _useCache;
}

void InternalGrid::useCache(bool newUseCache)
{
    _useCache = newUseCache;
}

bool InternalGrid::createCacheFile(int threadIndex){

    QString name = QString("gridblocks_%1_%2.temp").arg(threadIndex).arg(_rasterid);
    QDir localDir(context()->cacheLocation().toLocalFile());
    if ( !localDir.exists()) {
        localDir.mkpath(localDir.absolutePath());
    }
    QString filepath = localDir.absolutePath() + "/" + name;
    _imageNames[threadIndex] = filepath;
    _diskImages[threadIndex] = new QFile(filepath);
    if (_diskImages[threadIndex]->open(QIODevice::ReadWrite)) {
        if (!_diskImages[threadIndex]->resize(_size.linearSize())){
            return false;
        }
    }else
        return false;

    return true;
}

void InternalGrid::prepare4Operation(int nThreads) {
    _diskImages.resize(nThreads + 1);
    _imageNames.resize(nThreads + 1);
    for(quint32 threadIndex = 1 ; threadIndex <= nThreads; ++threadIndex)
        createCacheFile(threadIndex);
}

Grid::Orientation InternalGrid::flow() const
{
    return _orientation;
}

bool InternalGrid::isDefined() const{
    return _validStripe.size() > 0 && _planes.size() > 0 && _rasterid != i64UNDEF;
}
//------------------------------------------------------------------
GridXYZ::GridXYZ()
{
   _valueFuncPtr = reinterpret_cast<ValueFunc>(&GridXYZ::value);
}


quint32 GridXYZ::blocksCount() const{
    return int(_size.ysize()/ XYZYBLOCKS) * _size.zsize();
}

unsigned int GridXYZ::blockFloor(unsigned int y) const{
    quint32 blockNr = quint32(rY(y) / XYZYBLOCKS);
    return  blockNr * XYZYBLOCKS;
}

void GridXYZ::clone(quint64 newRasterId, int dLimit1, int dLimit2, InternalGrid *grid)
{
    auto igrid = static_cast<GridXYZ *>(grid);
    igrid->_planes.resize(_size.zsize());
    quint32 start = dLimit1 == iUNDEF ? 0 : dLimit1;
    quint32 end = dLimit2 == iUNDEF ? _size.zsize() : dLimit2 + 1;
    igrid->prepare(newRasterId,Size<>(_size.xsize(), _size.ysize(), end - start), _orientation, _box);
    for(quint32 z = start;  z < end; ++z) {
        igrid->_planes[z] = Plane(_size.ysize(),0);
        for(quint32 y = 0; y < _size.ysize(); ++y){
            igrid->_planes[z][y] = RelativeVector<double>(_size.xsize(),0);
            for(quint32 x=0; x < _size.xsize(); ++x){
                igrid->value(z,y,x,0) = value(z,y,x,0);
            }
        }
    }
}

void GridXYZ::setStartIndexes() {
    _planes.startIndex(_box.min_corner().z);
    for(auto& plane1 : _planes){
        plane1.startIndex(_box.min_corner().y);
        for(auto& plane2 : plane1)
            plane2.startIndex(_box.min_corner().x);
    }
}

PIXVALUETYPE &GridXYZ::value(int z, int y, int x, int threadIndex)
{
    if (!_validStripe[z][y]._valid){
        load(y, z, z, x, threadIndex);
    }
    return _planes[z][y][x];
}

void GridXYZ::loadFromCache(int y, int z, std::vector<char>& data)
{
    if (_planes[z][y].size() != _size.xsize()){
        _planes[z][y] = RelativeVector<double>(_size.xsize(), _planes[0][0].startIndex());
    }
    std::memcpy(_planes[z][y].data(), &data[0], _size.xsize() * sizeof(double));
}


void GridXYZ::setBlockData(int z, int y, const std::vector<PIXVALUETYPE> &data)
{
    unsigned int yBase = y;
    quint32 maxStripeSize = linesPerBlock(y); // normalizedY + blockCacheLimit() < _size.ysize() ? blockCacheLimit() : _size.ysize() % blockCacheLimit();
    auto inpIter = data.begin();
    for(unsigned int localy = 0; localy < maxStripeSize; ++localy){
        auto yAbs = yBase + localy;
        if (!_validStripe[z][yBase]._valid){
            // initialzing the whole block as it is a new block
            _planes[z].resize(_size.ysize());
            _planes[z].startIndex(_box.min_corner().y);
            for(quint32 stripeD = 0 ; stripeD < maxStripeSize; ++stripeD){
                _planes[z][yBase + stripeD].resize(_size.xsize());
                _planes[z][yBase + stripeD].startIndex(_box.min_corner().x);
             }

            _validStripe[z][yAbs]._valid = true;
            _planes[z].startIndex(_box.min_corner().y);
        }
        auto it = _planes[z][yAbs].begin();
        std::copy(inpIter, inpIter + _size.xsize(), it);
        inpIter +=  _size.xsize();
        _validStripe[z][yAbs]._loadedFromSource = true;
    }
}

void GridXYZ::setValue(int z, int y, int x, PIXVALUETYPE v, int threadIndex)
{std::vector<char> data;
    if (!_validStripe[z][y]._valid)
        load(y, 0, z, x, threadIndex);
    _planes[z][x][y] = v;
}

void GridXYZ::setOrientation(Grid::Orientation ori)
{
    bool unloadMap = _validStripe.size() == 0 || _imageNames.size() == 0;
    _orientation = ori;
    if ( _size.isValid()){
        _validStripe.resize(_size.zsize());
        _validStripe.startIndex(_box.min_corner().z);
        for(int z = _box.min_corner().z; z <= _box.max_corner().z; ++z){
            _validStripe[z].resize(_size.ysize());
            for(int y=_box.min_corner().y; y <= _box.max_corner().y; ++y)
                _validStripe[z][y] = InternalGrid::BlockStatus(unloadMap ? false : true, false);
        }
        _planes = RelativeVector<Plane>();
        _planes.resize(rZ(_size.zsize()));
        _planes.startIndex(_box.min_corner().z);
    }
}

// same as ZXY
void GridXYZ::unload(int threadIndex)
{
    auto yBase = blockCacheLimit();
    for( int z = _box.min_corner().z; z < _box.max_corner().z; ++z){
        for(int y = _box.min_corner().y; y < _box.max_corner().y; y += yBase){
            dumpBlock(z,y, threadIndex);
        }
    }
}

void GridXYZ::setBlockStatus( unsigned int z1, quint64 blockNr, bool status)
{
    _validStripe[z1][blockNr]._valid = status;
}

quint32 GridXYZ::blockCacheLimit() const
{
    return XYZYBLOCKS;
}

void GridXYZ::set3DSize(int z, int y, int xzs)
{
    if ( rZ(z) < _planes.size()){
        auto& plane = _planes[z];
        if ( plane.size() != _size.ysize()){
            plane.resize( _size.ysize());
            plane.startIndex(_box.min_corner().y);
        }
        quint64 yBase = blockFloor(y);
        quint32 maxd2 = linesPerBlock(y);
        for(quint32 yl=0; yl < maxd2; ++yl){
            plane[yBase + yl].resize(xzs, rUNDEF);
            plane[yBase + yl].startIndex(_box.min_corner().y);
            _validStripe[z][yBase + yl]._valid = true;
            _validStripe[z][yBase + yl]._loadedFromSource = true;
        }

    }
}

quint32 GridXYZ::linesPerBlock(int y) const
{
    quint64 blockNr = int((rY(y) / XYZYBLOCKS));
    unsigned int normalizedd = blockNr * XYZYBLOCKS;
    return normalizedd + XYZYBLOCKS <= _size.ysize() ? XYZYBLOCKS : _size.ysize() % XYZYBLOCKS;
}

void GridXYZ::dump(unsigned int yBase, unsigned int z1, unsigned int rz2, int threadIndex) {
    unsigned int zsize, ysize;
    std::vector<char> data;
    sizeData(yBase, rz2, ysize, zsize, data);
    for (unsigned int z = z1; z < zsize; ++z){
        for(unsigned int y = yBase - linesPerBlock(yBase); y < ysize; ++y){
            quint64 seekposition = calcSeekPosition(z,y);
            memcpy(&data[0], _planes[z][y].data(), _size.xsize() * sizeof(PIXVALUETYPE));
            _planes[z][y] = RelativeVector<double>();
            _diskImages[threadIndex]->seek(seekposition);
            _diskImages[threadIndex]->write(&data[0], _size.xsize()*sizeof(double));
            setBlockStatus(z, y, false);
        }
    }

}

void GridXYZ::closure()
{
    auto clim = blockCacheLimit();
    for( unsigned int z = _box.min_corner().z; z < _box.max_corner().z; ++z){
        for(quint32 yb = _box.min_corner().y; yb < _box.max_corner().y; yb += clim){
            const auto& st = _validStripe[z][yb];
            if ( st._valid){
                dumpBlock(z,yb,0);
            }
        }
    }
}

quint64 GridXYZ::blockSize(int d) const {
    auto lines = linesPerBlock(d);
    return lines * _size.xsize();
}

const InternalGrid::BlockStatus &GridXYZ::blockStatus(int z, int y) const
{
    return _validStripe[z][y];
}

quint64 GridXYZ::seekPosition(int rz, int ry, int rx) const
{
    auto yBase = rY(blockFloor(ry));
    return (rz * _size.xsize() * _size.ysize() + yBase * _size.xsize()) * sizeof(PIXVALUETYPE);
}

//----------------------------------------------------------------------------------

GridZXY::GridZXY()
{
    _valueFuncPtr = reinterpret_cast<ValueFunc>(&GridZXY::value);
}

void GridZXY::clone(quint64 newRasterId, int dLimit1, int dLimit2, InternalGrid *grid)
{
    auto igrid = static_cast<GridZXY *>(grid);
    igrid->_planes.resize(_size.ysize());
    quint32 start = dLimit1 == iUNDEF ? 0 : dLimit1;
    quint32 end = dLimit2 == iUNDEF ? _size.ysize() : dLimit2 + 1;
    igrid->prepare(newRasterId,Size<>(_size.xsize(), _size.ysize(), _size.zsize()), _orientation, _box);
    for(quint32 y = start;  y < end; ++y) {
         igrid->_planes[y] = RelativeVector<RelativeVector<double>>(_size.zsize());
         for(quint32 z = 0; z < _size.zsize(); ++z){
             igrid->_planes[y][z] = RelativeVector<double>(_size.xsize(),0);
             for(quint32 x=0; x < _size.xsize(); ++x)
                 igrid->value(z,y,x,0) = value(z,y,x,0);
         }
    }
}


PIXVALUETYPE &GridZXY::value(int z, int y, int x, int threadIndex)
{
    z = z == iUNDEF ? 0 : z;
    if ( !_validStripe[y][z]._valid){
        load(y, z, _size.zsize()-1, x, threadIndex);
    }
    return _planes[y][z][x];
}

void GridZXY::setStartIndexes() {
    _planes.startIndex(_box.min_corner().y);
    for(auto& plane1 : _planes){
        plane1.startIndex(_box.min_corner().z);
        for(auto& plane2 : plane1)
            plane2.startIndex(_box.min_corner().x);
    }
}
void GridZXY::setBlockData(int z, int y, const std::vector<PIXVALUETYPE> &data)
{
    int yBase = y;
    int lpb = linesPerBlock(y);
    if (!_validStripe[yBase][z]._valid){

        for(int ly = yBase; ly < yBase + lpb; ++ly){
            _planes[ly].resize(_size.zsize());
            _planes[ly].startIndex(_box.min_corner().z);
            for(quint32 zl = 0 ; zl < _size.zsize(); ++zl){
                _planes[ly][zl].resize(_size.xsize());
                _planes[ly][zl].startIndex(_box.min_corner().x);
            }
        }
    }
    auto inpIter = data.begin();
    for(int lz = _box.min_corner().z; lz <= _box.max_corner().z; ++lz){
        for(int ly = yBase; ly < yBase + lpb; ++ly){
            std::copy(inpIter, inpIter + _size.xsize(), _planes[ly][lz].begin() );
            inpIter += _size.xsize();
            _validStripe[ly][lz]._loadedFromSource = true;
            _validStripe[ly][lz]._valid = true;
        }
    }
}

void GridZXY::setValue(int z, int y, int x, PIXVALUETYPE v, int threadIndex)
{
    if (!_validStripe[y][z]._valid)
         load(y, z, _size.zsize()-1, x, threadIndex);
     _planes[z][x][y] = v;
}

void GridZXY::setOrientation(Grid::Orientation ori)
{
    bool unloadMap = _validStripe.size() == 0 || _imageNames.size() == 0;
    _orientation = ori;
    _validStripe.resize(_size.ysize());
    _validStripe.startIndex(blockFloor(_box.min_corner().y));
    for(unsigned int y = _box.min_corner().y; y <= _box.max_corner().y; ++y){
        _validStripe[y].resize(_size.zsize());
        for(unsigned int z = _box.min_corner().z; z <= _box.max_corner().z; ++z){
            _validStripe[y][z] = InternalGrid::BlockStatus(unloadMap ? false : true, false);
        }
    }
    _planes = RelativeVector<Plane>();
    _planes.resize(rY(_size.ysize()));
}

void GridZXY::unload(int threadIndex)
{
    auto clim = blockCacheLimit();
    for( unsigned int z = _box.min_corner().z; z < _box.max_corner().z; ++z){
        for(unsigned int yb = _box.min_corner().y; yb < _box.max_corner().y; yb += clim){
            dumpBlock(z,yb, threadIndex);
        }
    }
}

void GridZXY::set3DSize(int z, int y, int xzs)
{
    if ( z < _planes.size()){

        quint32 maxd2 = size().zsize();
        quint64 blockNr = int((y / blockCacheLimit()));
        quint64 ylines = linesPerBlock(y);
        quint64 yBase = blockFloor(y);
        for ( quint32 y1 = yBase; y1 < yBase + ylines; ++y1){
            auto& plane = _planes[y1];
            if ( plane.size() != _size.zsize())
                plane.resize( _size.zsize());
            for(quint32 zl=0; zl < maxd2; ++zl)
                plane[zl].resize(xzs, rUNDEF);
        }

        _validStripe[blockNr][0]._valid = true;
        _validStripe[blockNr][0]._loadedFromSource = true;

    }
}

quint32 GridZXY::linesPerBlock(int d) const
{
    quint64 blockNr = int((d / ZXYYBLOCKS));
    unsigned int normalizedd = blockNr * ZXYYBLOCKS;
    return normalizedd + ZXYYBLOCKS <= _size.ysize() ? ZXYYBLOCKS: _size.ysize() % ZXYYBLOCKS;
}

void GridZXY::loadFromCache(int y, int z, std::vector<char>& data)
{
    if ( _planes[y].size() != _size.zsize())
        _planes[y].resize(_size.zsize());
    if ( _planes[y][z].size() != _size.xsize())
        _planes[y][z] = RelativeVector<double>(_size.xsize(),0);
    std::memcpy(_planes[y][z].data(), &data[0], _size.xsize() * sizeof(double));
}
void GridZXY::closure()
{
    auto clim = blockCacheLimit();
    for(quint32 yb = _box.min_corner().y; yb < _box.max_corner().y; yb += clim){
        for(quint32 z = _box.min_corner().z; z < _box.max_corner().z; ++z){
            const auto& st = _validStripe[yb][z];
            if ( st._valid){
                dumpBlock(z,yb,0);
            }
        }
    }
}

void GridZXY::dump(unsigned int yBase, unsigned int z1, unsigned int z2, int threadIndex) {
    unsigned int zsize, ysize;
    std::vector<char> data;
    sizeData(yBase, z2, ysize, zsize, data);
    for (unsigned int z = z1; z < zsize; ++z){
        for(unsigned int y = yBase - linesPerBlock(yBase); y < ysize; ++y){
            quint64 seekposition = calcSeekPosition(z,y);
            _diskImages[threadIndex]->seek(seekposition);
            memcpy(&data[0], _planes[y][z].data(), _size.xsize() * sizeof(PIXVALUETYPE));
            _planes[y][z] = RelativeVector<double>();
            _diskImages[threadIndex]->write(&data[0], _size.xsize()*sizeof(double));
            setBlockStatus(z, y, false);
        }
    }

}
quint32 GridZXY::blocksCount() const{
   return int(_size.xsize()/ ZXYYBLOCKS) * _size.ysize();
}

quint64 GridZXY::blockSize(int d) const {
    auto lines = linesPerBlock(d);
    return lines * _size.ysize();
}

const InternalGrid::BlockStatus& GridZXY::blockStatus(int z, int blockIndex) const{
    return _validStripe[blockIndex][z];
}

quint64 GridZXY::seekPosition(int z, int y, int x) const
{
    auto normalizedY = blockFloor(y);
    return (normalizedY * _size.xsize()* _size.zsize() + x * _size.zsize()) * sizeof(PIXVALUETYPE);
}

void GridZXY::setBlockStatus( unsigned int z, quint64 y, bool status)
{
    _validStripe[y][z]._valid = status;
}

quint32 GridZXY::blockCacheLimit() const
{
    return ZXYYBLOCKS;
}

unsigned int GridZXY::blockFloor(unsigned int y) const{
    quint64 blockNr = int(y / ZXYYBLOCKS);
    return  blockNr * ZXYYBLOCKS;
}
