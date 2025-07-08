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

Grid::Grid(Orientation ori) : _orientation(ori), _rasterid(i64UNDEF)  {
}


Grid::~Grid() {
    clear();

}

Size<> Grid::size() const {
    return _size;
}


Grid *Grid::clone(quint64 newRasterId, int dLimit1, int dLimit2)
{
    Locker<> lock(_mutex);
    quint32 start = dLimit1 == iUNDEF ? 0 : dLimit1;
    Grid *grid = new Grid(_orientation);
    grid->_useCache = _useCache;
    if ( _orientation == oXYZ){
        quint32 end = dLimit2 == iUNDEF ? _size.zsize() : dLimit2 + 1;
        grid->prepare(newRasterId,Size<>(_size.xsize(), _size.ysize(), end - start), _orientation);
        grid->_planes.resize(_size.zsize());
        for(quint32 z = start;  z < end; ++z) {
            grid->_planes[z] = std::vector<std::vector<double>>(_size.ysize());
            for(quint32 y = 0; y < _size.ysize(); ++y){
                grid->_planes[z][y] = std::vector<double>(_size.xsize());
                for(quint32 x=0; x < _size.xsize(); ++x)
                    grid->value(z,y,x,0) = value(z,y,x,0);
            }
        }
    }else {
       grid->_planes.resize(_size.ysize());
       quint32 end = dLimit2 == iUNDEF ? _size.ysize() : dLimit2 + 1;
       grid->prepare(newRasterId,Size<>(_size.xsize(), _size.ysize(), _size.zsize()), _orientation);
       for(quint32 y = start;  y < end; ++y) {
            grid->_planes[y] = std::vector<std::vector<double>>(_size.zsize());
            for(quint32 z = 0; z < _size.zsize(); ++z){
                grid->_planes[y][z] = std::vector<double>(_size.xsize());
                for(quint32 x=0; x < _size.xsize(); ++x)
                    grid->value(z,y,x,0) = value(z,y,x,0);
            }
       }
    }
    grid->closure();
    return grid;

}

void Grid::clear() {
    Locker<> lock(_mutex);
    _validStripe = std::vector<std::vector<BlockStatus>>();
    _planes = std::vector<Plane>();
    if ( _diskImages.size() > 0){
        for(auto *file : _diskImages){
            file->close();
            file->remove();
        }
        _diskImages = std::vector<QFile*>();
    }
    _size = Size<>();
    _rasterid = iUNDEF;
    _orientation = oXYZ;

}

quint64 Grid::calcSeekPosition(unsigned int z, unsigned int y ) const{
   quint64 pos = (z * _size.xsize() * _size.ysize() + y * _size.xsize()) * sizeof(PIXVALUETYPE);
   return pos;
}

void Grid::dump(unsigned int y1, unsigned int y2, unsigned int z1, unsigned int z2, int threadIndex) {
    quint64 dataSize = _size.xsize() * sizeof(PIXVALUETYPE);
    auto zsize = std::min(_size.zsize()-1, z2);
    auto ysize = std::min(_size.ysize()-1,y2);
    std::vector<char> data(dataSize);
    for (unsigned int z = z1; z <= zsize; ++z){
        for(unsigned int y = y1; y <= ysize; ++y){
            quint64 seekposition = calcSeekPosition(z,y);
            _diskImages[threadIndex]->seek(seekposition);
            if ( _orientation == Grid::oXYZ){
                memcpy(&data[0], _planes[z][y].data(), _size.xsize() * sizeof(PIXVALUETYPE));
                _planes[z][y] = std::vector<double>();
            }else {
                memcpy(&data[0], _planes[y][z].data(), _size.xsize() * sizeof(PIXVALUETYPE));
                 _planes[y][z] = std::vector<double>();
            }
            _diskImages[threadIndex]->write(&data[0], _size.xsize()*sizeof(double));
        }
    }
    quint64 blockNr = int((y1 / blockCacheLimit()));
    setBlockStatus(z1, blockNr, false);
}

void Ilwis::Grid::loadFromCache(unsigned int y1, unsigned int y2, unsigned int z1, unsigned int z2, unsigned int x, int threadIndex)
{
    quint64 blockNr = int((y1 / blockCacheLimit()));

    quint64 dataSize = _size.xsize() * sizeof(PIXVALUETYPE);
    std::vector<char> data(dataSize);
    for (unsigned int z = z1; z <= std::min(_size.zsize()-1,z2); ++z){
        for(unsigned int y = y1; y <= std::min(_size.ysize()-1,y2); ++y){
            quint64 seekposition = calcSeekPosition(z, y);
            _diskImages[threadIndex]->seek(seekposition);
            _diskImages[threadIndex]->read(&data[0], dataSize);
            if ( _orientation == Grid::oXYZ){
                if (_planes[z][y].size() != _size.xsize())
                    _planes[z][y] = std::vector<double>(_size.xsize());
                std::memcpy(_planes[z][y].data(), &data[0], _size.xsize() * sizeof(double));
            }else{
                if ( _planes[y].size() != _size.zsize())
                    _planes[y].resize(_size.zsize());
                if ( _planes[y][z].size() != _size.xsize())
                    _planes[y][z] = std::vector<double>(_size.xsize());
                std::memcpy(_planes[y][z].data(), &data[0], _size.xsize() * sizeof(double));
            }
        }
    }
    setBlockStatus(z1, blockNr, true);
}

void Ilwis::Grid::setBlockStatus( unsigned int z1, quint64 blockNr, bool status)
{
    if ( _orientation == oXYZ)
        _validStripe[z1][blockNr]._valid = status;
    else
        _validStripe[blockNr][0]._valid = status;
}

void Grid::load(unsigned int y1, unsigned int y2, unsigned int z1, unsigned int z2, unsigned int x, int threadIndex){
    quint64 blockNr = int((y1 / blockCacheLimit()));
    IIlwisObject obj = mastercatalog()->get(_rasterid);
    if (obj.isValid() || _imageNames.size() > 0 ) { //.is used in the clone function; an emptry grid is created without(yet) a raster belonging to it
        if (blockStatus(z1,blockNr)._loadedFromSource && _useCache){
                loadFromCache(y1,y2,z1, z2,  x, threadIndex);
            }else{
                loadFromSource(z1,x,y1);
            }
    }else {
       setBlockStatus(z1, blockNr, true);
    }


    if ( pastHorizon(y1)){
        dumpBlock(z1,y1,x, threadIndex);
    }

}

void Ilwis::Grid::dumpBlock(int z, int y, int x, int threadIndex)
{
    if (_diskImages.size() <= (size_t)threadIndex && _useCache){
        _diskImages.resize(threadIndex+1);
        _imageNames.resize(threadIndex + 1);
        createCacheFile(threadIndex);
    }
    quint32 blockNr = int(y/blockCacheLimit()) - 1;
        auto yb = blockNr  * blockCacheLimit();
        auto yl = ( blockNr + 1) * blockCacheLimit() - 1;
        if ( _orientation == Grid::oXYZ && blockStatus(z, blockNr)._valid){
            dump(yb , yl, z, z,threadIndex);
        }else if ( blockStatus(0, blockNr)._valid){
            dump(yb , yl, 0, _size.zsize()-1, threadIndex);
        }

}

void Grid::loadFromSource(int z1, int x, int y1){
    IIlwisObject obj = mastercatalog()->get(_rasterid);
    if ( obj.isValid() ){
        quint64 blockNr = int((y1 / blockCacheLimit()));
        IRasterCoverage raster = obj.as<RasterCoverage>();
        IOOptions options;
        quint64 normalizedY = blockNr * blockCacheLimit() ;
        options.addOption({"z", z1});
        options.addOption({"y", normalizedY});
        options.addOption({"size", blockCacheLimit() * _size.xsize()});
        options.addOption({"orientation", _orientation == oZXY ? "ZXY" : "XYZ"});
        raster->getData(options);
    }
}

PIXVALUETYPE& Grid::valueRef(const Pixel& pix, int threadIndex) {
    quint32 z = pix.z == iUNDEF ? 0 : pix.z;
    if (z >= 0 && z < _size.zsize() && pix.y >= 0 && pix.y < _size.ysize() && pix.x >= 0 && pix.x < _size.xsize())
        return value(z, pix.y, pix.x, threadIndex);
    throw ErrorObject("Pixel location out of bounds");
}

PIXVALUETYPE Grid::value(const Pixel& pix, int threadIndex) {
    quint32 z = pix.z == iUNDEF ? 0 : pix.z;
    if (z >= 0 && z < _size.zsize() && pix.y >= 0 && pix.y < _size.ysize() && pix.x >= 0 && pix.x < _size.xsize())
        return value(z, pix.y, pix.x, threadIndex);
    else
        return PIXVALUEUNDEF;
}

PIXVALUETYPE &Grid::value(int z, int y, int x, int threadIndex)  {
     quint64 blockNr = int((y / blockCacheLimit()));
     auto yb = blockNr  * blockCacheLimit();
     auto yl = ( blockNr + 1) * blockCacheLimit() - 1;
     switch(_orientation){
     case oZXY:{

             if (!_validStripe[blockNr][0]._valid){
                 load(yb, yl, 0, _size.zsize()-1, x, threadIndex);
             }
             return _planes[y][z][x];
        }
    default:{

         if (!_validStripe[z][blockNr]._valid)
             load(yb,yl, z, z, x, threadIndex);
         return _planes[z][y][x];
         }
     }
}

void Grid::setBlockDataXYZ(int z, int y, const std::vector<PIXVALUETYPE>& data){
    quint64 blockNr = int((y / blockCacheLimit()));
    int normalizedY = blockNr * blockCacheLimit();
    quint32 maxStripeSize = linesPerBlock(y); // normalizedY + blockCacheLimit() < _size.ysize() ? blockCacheLimit() : _size.ysize() % blockCacheLimit();
    auto inpIter = data.begin();
    for(unsigned int localy = 0; localy < maxStripeSize; ++localy){
        if (!_validStripe[z][blockNr]._valid){
            _planes[z].resize(_size.ysize());
            for(quint32 stripeD = 0 ; stripeD < maxStripeSize; ++stripeD)
                _planes[z][normalizedY + stripeD].resize(_size.xsize());
            _validStripe[z][blockNr]._valid = true;
        }
        std::copy(inpIter, inpIter + _size.xsize(), _planes[z][normalizedY + localy].begin());
        inpIter +=  _size.xsize();
    }
    _validStripe[z][blockNr]._loadedFromSource = true;
}

void Grid::setBlockDataZXY(int z, int y, const std::vector<PIXVALUETYPE>& data){
    auto climit = blockCacheLimit();
    quint64 blockNr = int((y / climit));
    int normalizedY = blockNr * blockCacheLimit();
    quint32 maxStripeSize = _size.ysize() - y < climit ?_size.ysize() % climit : climit;

    if (!_validStripe[blockNr][0]._valid){

        for(unsigned int localy = 0; localy < maxStripeSize; ++localy){
            _planes[normalizedY +localy].resize(_size.zsize());
            for(quint32 zl = 0 ; zl < _size.zsize(); ++zl)
                _planes[normalizedY +localy][zl].resize(_size.xsize());
        }
    }
    auto inpIter = data.begin();
    for(unsigned int lz = 0; lz < _size.zsize(); ++lz){
        for(unsigned int localy = 0; localy < maxStripeSize; ++localy){
            std::copy(inpIter, inpIter + _size.xsize(), _planes[normalizedY + localy][lz].begin() );
            inpIter += _size.xsize();
        }
    }
    _validStripe[blockNr][0]._loadedFromSource = true;
    _validStripe[blockNr][0]._valid = true;
}
void Grid::setBlockData(int z, int y, const std::vector<PIXVALUETYPE>& data){
    switch(_orientation){
    case oZXY:
        setBlockDataZXY(z,y,data);
        break;
    default:
        setBlockDataXYZ(z,y,data);
        break;
    }
}

void Grid::setValue(int z, int y, int x, PIXVALUETYPE v, int threadIndex ){
    quint64 blockNr = int((y / blockCacheLimit()));
    switch(_orientation){
    case oZXY:{
        if (!_validStripe[y][blockNr]._valid)
             load(y, y + blockNr * blockCacheLimit(), z, _size.zsize()-1, x, threadIndex);
         _planes[z][x][y] = v;
        }
        break;
    default:{
            if (!_validStripe[y][blockNr]._valid)
                load(y, y + blockNr * blockCacheLimit(), 0, z, x, threadIndex);
             _planes[z][x][y] = v;
        }
        break;
    }
}

bool Grid::prepare(quint64 rasterid, const Size<> &sz,Grid::Orientation ori) {
    Locker<> lock(_mutex);
    clear();
    _rasterid = rasterid;
    _size = sz;

    int maxBlock = sz.ysize() / blockCacheLimit() + 1;
    for(auto& bs : _validStripe){
        bs.resize(maxBlock);
    }
    setOrientation(ori);
    return true;
}

void Grid::setOrientation(Grid::Orientation ori){
    auto unloadMap = _validStripe.size() == 0 || _imageNames.size() ==0;
    if ( ori != _orientation || unloadMap){
        unloadMap = true;
        closure();
        _orientation = ori;
        _validStripe = std::vector<std::vector<BlockStatus>>();

        switch(_orientation){
            case oZXY:
            _validStripe.resize(_size.ysize()/blockCacheLimit() + 1);
            for(unsigned int y = 0; y < _validStripe.size(); ++y){
                _validStripe[y] = std::vector<BlockStatus>(_size.zsize(), {false,unloadMap ? false : true});
            }
            _planes = std::vector<Plane>();
             _planes.resize(_size.ysize());
                break;
            default:{
                _validStripe.resize(_size.zsize());
                for(unsigned int z = 0; z < _validStripe.size(); ++z){
                    unsigned int numBlocks = _size.ysize() / blockCacheLimit();
                    _validStripe[z] = std::vector<BlockStatus>(numBlocks + 1, {unloadMap ? false : true, false});
                }
                _planes = std::vector<Plane>();
                _planes.resize(_size.zsize());
            }
        }
    }

}
bool Grid::createCacheFile(int threadIndex){

    QString name = QString("gridblocks_%1_%2.temp").arg(threadIndex).arg(_rasterid);
    QDir localDir(context()->cacheLocation().toLocalFile());
    if ( !localDir.exists()) {
        localDir.mkpath(localDir.absolutePath());
    }
    QString filepath = localDir.absolutePath() + "/" + name;
    _imageNames[threadIndex] = filepath;
    _diskImages[threadIndex] = new QFile(filepath);
    if (_diskImages[threadIndex]->open(QIODevice::ReadWrite)) {
        if (!_diskImages[threadIndex]->resize(_size.linearSize() * sizeof(PIXVALUETYPE))){
            return false;
        }
    }else
        return false;
    return true;
}


void Grid::unload(int threadIndex) {
    Locker<> lock(_mutex);
    auto clim = blockCacheLimit();
    if ( _orientation == oXYZ){
        for( unsigned int z = 0; z < _size.zsize(); ++z){
            for(unsigned int block = 0; block < _validStripe[z].size(); ++block){
                auto y = block * clim;
                dumpBlock(z,y,iUNDEF, threadIndex);
            }
        }
    }else{
        for( unsigned int y = 0; y < _size.ysize(); ++y){
            for(unsigned int block = 0; block < _validStripe[y].size(); ++block){
                auto y = block * clim;
                dumpBlock(0,y,iUNDEF, threadIndex);
            }
        }
    }
}

bool Grid::isValid() const
{
    return !(_size.isNull() || _size.isValid());
}

void Grid::prepare4Operation(int nThreads) {
    _diskImages.resize(nThreads + 1);
    _imageNames.resize(nThreads + 1);
    for(quint32 threadIndex = 1 ; threadIndex <= nThreads; ++threadIndex)
        createCacheFile(threadIndex);
}

void Grid::unprepare4Operation() {

}

quint32 Grid::blocksCount() const{
    if ( _orientation == oZXY)
        return int(_size.xsize()/ blockCacheLimit()) * _size.ysize();
    return int(_size.ysize()/ blockCacheLimit()) * _size.zsize();
}

quint64 Grid::blockSize(int d) const {
    auto lines = linesPerBlock(d);
    if (_orientation == oZXY){
        return lines * _size.ysize();
    }
    return lines * _size.xsize();
}

void Grid::setd3Size(int z, int y, int xzs){
    if ( _orientation == oXYZ){
        if ( z < _planes.size()){
            auto& plane = _planes[z];
            if ( plane.size() != _size.ysize())
                plane.resize( _size.ysize());
            quint64 blockNr = int((y / blockCacheLimit()));
            quint32 maxd2 = linesPerBlock(y);
            for(quint32 yl=0; yl < maxd2; ++yl)
                plane[blockNr * blockCacheLimit() + yl].resize(xzs, rUNDEF);

            _validStripe[z][blockNr]._valid = true;
            _validStripe[z][blockNr]._loadedFromSource = true;

        }
    }else {
        if ( z < _planes.size()){

            quint32 maxd2 = size().zsize();
            quint64 blockNr = int((y / blockCacheLimit()));
            quint64 ylines = linesPerBlock(y);
            quint64 ymin = blockNr * blockCacheLimit();
            for ( quint32 y1 = ymin; y1 < ymin + ylines; ++y1){
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
}
quint32 Grid::linesPerBlock(int d) const
{
    if ( _orientation == oXYZ){
        quint64 blockNr = int((d / blockCacheLimit()));
        unsigned int normalizedd = blockNr * blockCacheLimit();
        quint32 maxLinesPerBlock = normalizedd + blockCacheLimit() <= _size.ysize() ? blockCacheLimit() : _size.ysize() % blockCacheLimit();
        return maxLinesPerBlock;
    }else{
        quint64 blockNr = int((d / blockCacheLimit()));
        unsigned int normalizedd = blockNr * blockCacheLimit();
        quint32 maxLinesPerBlock = normalizedd + blockCacheLimit() <= _size.ysize() ? blockCacheLimit() : _size.ysize() % blockCacheLimit();
        return maxLinesPerBlock;
    }
}

quint32 Grid::blockCacheLimit() const{
    if ( _orientation == oZXY)
        return 8;
    return 100;
}

const Grid::BlockStatus& Grid::blockStatus(int z, int blockIndex) const{
    if ( _orientation == oZXY)
        return _validStripe[blockIndex][0];
    return _validStripe[z][blockIndex];
}

quint64 Grid::seekPosition(int z, int y, int x) const{
    auto normalizedY = normY(y);
    if ( _orientation == oZXY)    {
        return (normalizedY * _size.xsize()* _size.zsize() + x * _size.zsize()) * sizeof(PIXVALUETYPE);
    }else
        return (z * _size.xsize() * _size.ysize() + normalizedY * _size.xsize()) * sizeof(PIXVALUETYPE);
}

quint32 Grid::normY(int y) const{
    quint64 blockNr = int((y / blockCacheLimit()));
    return blockNr * blockCacheLimit() ;
}

bool Grid::pastHorizon(quint32 v) const {
    return (v + 1) % blockCacheLimit() == 0 && v != 0;
}

void Grid::closure()
{
    if ( _validStripe.size() > 0){
        if ( _orientation == oXYZ){
            for(quint32 z=0; z < _size.zsize(); ++z){
                for(quint32 yb=0; yb < _validStripe[z].size(); ++yb){
                    const auto& st = _validStripe[z][yb];
                    if ( st._valid){
                        auto y = (yb + 1) * blockCacheLimit();
                        dumpBlock(z,y,0,0);
                    }
                }
            }
        }else{
            for(quint32 yb=0; yb < _validStripe.size(); ++yb){
                for(quint32 z=0; z < _size.zsize(); ++z){
                    const auto& st = _validStripe[yb][z];
                    if ( st._valid){
                        auto y = (yb + 1) * blockCacheLimit();
                        dumpBlock(z,y,0,0);
                    }
                }
            }
        }
    }
}

void Grid::useCache(bool yesno){
    _useCache = yesno;
}

Grid::Orientation Grid::flow() const{
    return _orientation;
}
