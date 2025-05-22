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
    quint32 end = dLimit2 == iUNDEF ? (_orientation == oXYZ ? _size.zsize() : _size.xsize()) : dLimit2 + 1;
    Grid *grid = new Grid(_orientation);
    grid->prepare(newRasterId,Size<>(_size.xsize(), _size.ysize(), end - start));

    if ( _orientation == oXYZ){
        grid->_planes.resize(_size.zsize());
        for(quint32 z = start;  z < end; ++z) {
            grid->_planes[z] = std::vector<std::vector<double>>(_size.ysize());
            for(quint32 y = 0; y < _size.ysize(); ++y){
                grid->_planes[z][y] = std::vector<double>(_size.xsize());
                for(quint32 x=0; x < _size.xsize(); ++x)
                    grid->value(z,y,x,0) = value(z,y,x,0);
            }
        }
    }
    grid->_validStripe = _validStripe;
    grid->_useCache = _useCache;
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

void Grid::dumpXYZ(int z, int blockNr, int threadIndex){

    quint64 normalizedY = blockNr * blockCacheLimit() ;
    quint32 maxStripeSize = linesPerBlock(normalizedY);
    if ( _planes[z][normalizedY].size() == 0) // already dumped
        return;
    if (!_useCache){
        for(unsigned int localy = 0; localy < maxStripeSize; ++localy){
            _planes[z][normalizedY + localy] = std::vector<double>();
        }
    } else {
        quint64 seekposition = (z * _size.xsize() * _size.ysize() + normalizedY * _size.xsize()) * sizeof(PIXVALUETYPE);
        quint64 dataSize = _size.xsize() * maxStripeSize * sizeof(PIXVALUETYPE);
        char data[dataSize];
        quint64 p = 0;
        for(unsigned int localy = 0; localy < maxStripeSize; ++localy){
            memcpy(&data[p], _planes[z][normalizedY + localy].data(), _size.xsize() * sizeof(PIXVALUETYPE));
            _planes[z][normalizedY + localy] = std::vector<double>();
            p += _size.xsize() * sizeof(PIXVALUETYPE);

        }

        _diskImages[threadIndex]->seek(seekposition);
        _diskImages[threadIndex]->write(data, maxStripeSize * _size.xsize() * sizeof(PIXVALUETYPE));

    }
    _validStripe[z][blockNr]._valid = false;
}

void Grid::dumpZXY(int z, int blockNr, int threadIndex){
    quint64 normalizedY = blockNr * blockCacheLimit() ;
    quint32 maxStripeSize = linesPerBlock(normalizedY);
    if ( _planes[normalizedY][z].size() == 0) // already dumped
        return;
    if (!_useCache){
        for(unsigned int localy = 0; localy < maxStripeSize; ++localy){
            for(quint32 zl=0 ; zl < _size.zsize(); ++zl){
                quint32 ly = normalizedY + localy;
                _planes[ly][zl] = std::vector<double>();
            }
        }
    }else {
        quint64 seekposition = (blockNr * blockCacheLimit() * _size.xsize() * _size.zsize()) * sizeof(PIXVALUETYPE);
        _diskImages[threadIndex]->seek(seekposition);

        quint64 dataSize = _size.xsize() * sizeof(PIXVALUETYPE);
        char data[dataSize];

        for(unsigned int localy = 0; localy < maxStripeSize; ++localy){
            for(quint32 zl=0 ; zl < _size.zsize(); ++zl){
                quint32 ly = normalizedY + localy;
                memcpy(&data[0], _planes[ly][zl].data(), _size.xsize() * sizeof(PIXVALUETYPE));
                _diskImages[threadIndex]->write(data, _size.xsize()*sizeof(double));

                _planes[ly][zl] = std::vector<double>();

                seekposition += _size.xsize() * sizeof(double);
                _diskImages[threadIndex]->seek(seekposition);
            }
        }
    }
     _validStripe[blockNr][0]._valid = false;
 }

void Ilwis::Grid::dumpBlock(int z, int y, int x, int threadIndex)
{
    if (_diskImages.size() <= (size_t)threadIndex && _useCache){
        _diskImages.resize(threadIndex+1);
        _imageNames.resize(threadIndex + 1);
        createCacheFile(threadIndex);
    }
    quint32 blockNr = int(y/blockCacheLimit()) - 1;
    if ( _orientation == Grid::oXYZ){
        dumpXYZ(z, blockNr, threadIndex);
    }else{
        dumpZXY(z,blockNr, threadIndex);
    }
}
void Ilwis::Grid::loadBlockXYZ(int z, int y, int x, int threadIndex){
    quint64 blockNr = int((y / blockCacheLimit()));
    if (blockStatus(z,blockNr)._loadedFromSource && _useCache){
        quint64 normalizedY = blockNr * blockCacheLimit() ;
        quint64 seekposition = (z * _size.xsize() * _size.ysize() + normalizedY * _size.xsize()) * sizeof(PIXVALUETYPE);
        _diskImages[threadIndex]->seek(seekposition);
        quint32 maxStripeSize = linesPerBlock(normalizedY);
        quint64 dataSize = _size.xsize() * maxStripeSize * sizeof(PIXVALUETYPE);
        char data[dataSize];

        _diskImages[threadIndex]->read(data, dataSize);
        quint64 p = 0;
        for(unsigned int localy = 0; localy < maxStripeSize; ++localy){
            _planes[z][normalizedY + localy].resize(_size.xsize());
             std::memcpy(_planes[z][normalizedY + localy].data(), &data[p], _size.xsize() * sizeof(double));
             p += _size.xsize()*sizeof(PIXVALUETYPE);

        }
        _validStripe[z][blockNr]._valid = true;
    }else
         loadFromSource(z,x,y);

    if ( pastHorizon(y)){
        dumpBlock(z,y,x, threadIndex);
    }
}

void Ilwis::Grid::loadBlockZXY(int z, int y, int x, int threadIndex){
    quint64 blockNr = int((y / blockCacheLimit()));
    if (blockStatus(z,blockNr)._loadedFromSource && _useCache){
        quint64 normalizedY = blockNr * blockCacheLimit() ;
        quint64 seekposition = seekPosition(z,y,x);
        _diskImages[threadIndex]->seek(seekposition);
        quint32 maxStripeSize = linesPerBlock(normalizedY);
        quint64 dataSize = _size.xsize() * sizeof(PIXVALUETYPE);
        char data[dataSize];
        for(unsigned int localy = 0; localy < maxStripeSize; ++localy){
            for(quint32 zl=0 ; zl < _size.zsize(); ++zl){
                _diskImages[threadIndex]->read(data, dataSize);
                _planes[normalizedY + localy][zl].resize(_size.xsize());
                 std::memcpy(_planes[normalizedY + localy][zl].data(), &data, _size.xsize() * sizeof(double));
                 seekposition += _size.xsize() * sizeof(double);
                 _diskImages[threadIndex]->seek(seekposition);
            }
        }
        _validStripe[blockNr][0]._valid = true;
    }else{
        loadFromSource(z,x,y);
    }

    if ( pastHorizon(blockNr)){
        dumpBlock(z,y,x, threadIndex);
    }
}

void Grid::loadFromSource(int z, int x, int y){
    IIlwisObject obj = mastercatalog()->get(_rasterid);
    if ( obj.isValid() ){
        quint64 blockNr = int((y / blockCacheLimit()));
        IRasterCoverage raster = obj.as<RasterCoverage>();
        IOOptions options;
        quint64 normalizedY = blockNr * blockCacheLimit() ;
        options.addOption({"z", z});
        options.addOption({"y", normalizedY});
        options.addOption({"size", blockCacheLimit() * _size.xsize()});
        options.addOption({"orientation", _orientation == oZXY ? "ZXY" : "XYZ"});
        raster->getData(options);
    }
}

PIXVALUETYPE& Grid::value(const Pixel &pix, int threadIndex) {
    quint32 z = pix.z == iUNDEF ? 0 : pix.z;
    if ( z>= 0 && z < _size.zsize() && pix.y >=0 && pix.y < _size.ysize() && pix.x>=0 && _size.xsize())
        return value(z, pix.y, pix.x, threadIndex);
    throw ErrorObject("Pixel location out of bounds");
}

PIXVALUETYPE &Grid::value(int z, int y, int x, int threadIndex)  {

     switch(_orientation){
     case oZXY:{
         quint64 blockNr = int((y / blockCacheLimit()));
         if (!_validStripe[blockNr][0]._valid)
             loadBlockZXY(z,y,x, threadIndex);
         return _planes[y][z][x];
        }
    default:{
         quint64 blockNr = int((y / blockCacheLimit()));
         if (!_validStripe[z][blockNr]._valid)
             loadBlockXYZ(z,y, x, threadIndex);
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
    switch(_orientation){
    case oZXY:{
        quint64 blockNr = int((y / blockCacheLimit()));
        if (!_validStripe[y][blockNr]._valid)
            loadBlockZXY(z,y,x, threadIndex);
         _planes[z][x][y] = v;
        }
        break;
    default:{
            quint64 blockNr = int((y / blockCacheLimit()));
            if (!_validStripe[y][blockNr]._valid)
                loadBlockXYZ(z,y,x, threadIndex);
             _planes[z][x][y] = v;
        }
        break;
    }
}

bool Grid::prepare(quint64 rasterid, const Size<> &sz) {
    Locker<> lock(_mutex);
    clear();
    _rasterid = rasterid;
    _size = sz;

    int maxBlock = sz.ysize() / blockCacheLimit() + 1;
    for(auto& bs : _validStripe){
        bs.resize(maxBlock);
    }
    setOrientation(_orientation);
    return true;
}

void Grid::setOrientation(Grid::Orientation ori){
    if ( ori != _orientation || _validStripe.size() == 0){
        _orientation = ori;
        _validStripe = std::vector<std::vector<BlockStatus>>();

        switch(_orientation){
            case oZXY:
            _validStripe.resize(_size.ysize()/blockCacheLimit() + 1);
            for(unsigned int y = 0; y < _validStripe.size(); ++y){
                _validStripe[y].resize(1);
            }
            _planes = std::vector<Plane>();
             _planes.resize(_size.ysize());
                break;
            default:{
                _validStripe.resize(_size.zsize());
                for(unsigned int z = 0; z < _validStripe.size(); ++z){
                    unsigned int numBlocks = _size.ysize() / blockCacheLimit();
                    _validStripe[z].resize(numBlocks + 1);
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
            auto& plane = _planes[y];
            if ( plane.size() != _size.zsize())
                plane.resize( _size.zsize());
            quint32 maxd2 = std::min(blockCacheLimit(), size().ysize());
            for(quint32 yl=0; yl < maxd2; ++yl)
                plane[y + yl].resize(xzs, rUNDEF);

            _validStripe[y][0]._valid = true;
            _validStripe[y][0]._loadedFromSource = true;

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
    return v % blockCacheLimit() == 0 && v != 0;
}

void Grid::closure()
{
    if ( _validStripe.size() > 0){
        for(quint32 z=0; z < _size.zsize(); ++z){
            for(quint32 yb=0; yb < _validStripe[z].size(); ++yb){
                const auto& st = _validStripe[z][yb];
                if ( st._valid){
                    auto y = (yb + 1) * blockCacheLimit(); // +1 because dump block , dumps the block before the current y, while we need the current y
                    dumpBlock(z,y,0,0);
                }
            }
        }
    }
}

void Grid::useCache(bool yesno){
    _useCache = yesno;
}

