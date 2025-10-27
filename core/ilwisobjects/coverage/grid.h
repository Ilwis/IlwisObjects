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

#ifndef Grid_H
#define Grid_H

#include <list>
#include <mutex>
#include <QDir>
#include <QTemporaryFile>
#include <iostream>
#include <deque>
#include "kernel.h"
#include "errorobject.h"
#include "size.h"
#include "location.h"
#include "box.h"
#include "relativevector.h"

namespace Ilwis {

class RasterCoverage;
class IOOptions;
class InternalGrid;

typedef std::unique_ptr<InternalGrid> IGrid;

typedef RelativeVector<RelativeVector<double>> Plane;

class KERNELSHARED_EXPORT Grid

{
public:
    enum Orientation {oXYZ, oZXY, oUNKNOWN};
    friend class PixelIterator;

    Grid();
    virtual ~Grid();

    void clear();


    PIXVALUETYPE &valueRef(const Pixel& pix, int threadIndex = 0);
    PIXVALUETYPE value(const Pixel& pix, int threadIndex = 0);
    PIXVALUETYPE& value(int z, int y, int x, int threadIndex = 0);
    void setValue(int z, int y, int x, PIXVALUETYPE v , int threadIndex=0);
    void setOrientation(Grid::Orientation ori);

    void setBlockData(int d1, int d2, const std::vector<PIXVALUETYPE>& data);
    bool prepare(quint64 rasterid, const Size<> &sz, Orientation ori = oXYZ, const BoundingBox &box = BoundingBox()) ;
    Size<> size() const;
    Grid * clone(quint64 newRasterId, int index1=iUNDEF, int index2=iUNDEF) ;
    void unload(int threadIndex=0);
    bool isValid() const;
    quint32 linesPerBlock(int y) const;
    quint32 blocksCount() const;
    quint64 blockSize(int d) const;
    void setd3Size(int d1, int d2, int d3sz);

    void prepare4Operation(int nThreads);
    void unprepare4Operation();
    quint32 blockCacheLimit() const;
    void useCache(bool yesno);
    Orientation flow() const;
    IGrid& igrid() { return _igrid;}
    void dump(unsigned int y1, unsigned int z1, unsigned int z2, int threadIndex);
    void load(unsigned int y1, unsigned int z1, unsigned int z2, unsigned int x, int threadIndex);
    quint64 seekPosition(int z, int y, int x) const;



private:
    void closure();

    std::recursive_mutex _mutex;
    IGrid _igrid;
};

class InternalGrid {
public:
    friend class GridXYZ;
    friend class GridZXY;

    struct BlockStatus{
        BlockStatus(bool valid=false, bool loaded=false) : _valid(valid), _loadedFromSource(loaded){}
        bool _valid = false;
        bool _loadedFromSource = false;
    };

    typedef PIXVALUETYPE& (InternalGrid::*ValueFunc)(int z, int y, int x, int threadIndex);

    InternalGrid();
    InternalGrid(Grid::Orientation ori);
    void dumpBlock(int z, int yBase, int threadIndex);
    void load(unsigned int y1, unsigned int z1, unsigned int z2, unsigned int x, int threadIndex);
    bool prepare(quint64 rasterid, const Size<> &sz, Grid::Orientation ori, const BoundingBox &box = BoundingBox());
    Size<> size() const { return _size; };
    BoundingBox boundingBox() const { return _box;}
    Grid::Orientation orientation() const { return _orientation;}
    bool unloadMap() const {return (_validStripe.size() == 0 || _imageNames.size() ==0);}
    void clear();
    bool isDefined() const;


    virtual void setOrientation(Grid::Orientation ori) = 0;
    virtual void clone(quint64 newRasterId, int dLimit1, int dLimit2, InternalGrid *grid) = 0;
    virtual void closure() = 0;
    virtual void dump(unsigned int y1, unsigned int z1, unsigned int z2, int threadIndex) =0;
    virtual void setBlockData(int z, int y, const std::vector<PIXVALUETYPE>& data) = 0;
    virtual void setValue(int z, int y, int x, PIXVALUETYPE v, int threadIndex) = 0;
    virtual void unload(int threadIndex)= 0;
    virtual quint32 blocksCount() const = 0;
    virtual quint64 blockSize(int d) const = 0;
    virtual void set3DSize(int z, int y, int xzs) = 0;
    virtual quint32 linesPerBlock(int d) const = 0;
    virtual quint32 blockCacheLimit() const = 0;
    virtual const BlockStatus& blockStatus(int z, int blockIndex) const = 0;
    virtual quint64 seekPosition(int z, int y, int x) const = 0;
    virtual void loadFromCache(int y, int z, std::vector<char>& data) = 0;
    virtual void setBlockStatus( unsigned int z1, quint64 blockNr, bool status) = 0;
    virtual unsigned int blockFloor(unsigned int y) const = 0;
    virtual void setStartIndexes() = 0;

    bool useCache() const;
    void useCache(bool newUseCache);
    bool createCacheFile(int i);
    void prepare4Operation(int nThreads);
    Grid::Orientation flow() const;
    unsigned int rZ(unsigned int z) const { return z - _box.min_corner().z;}
    unsigned int rY(unsigned int y) const { return y - _box.min_corner().y;}
    unsigned int rX(unsigned int x) const { return x - _box.min_corner().x;}

    ValueFunc _valueFuncPtr = nullptr;

protected:
    const unsigned int BLOCKYHORIZON = 4;
    const unsigned int BLOCKXHORIZON = 8;
    RelativeVector<Plane> _planes;
    RelativeVector<RelativeVector<BlockStatus>> _validStripe;
    Size<> _size;
    Grid::Orientation _orientation;
    quint64 _rasterid = i64UNDEF;
    std::vector<QFile*> _diskImages;
    std::vector<QString> _imageNames;
    const quint32 XYZYBLOCKS = 100;
    const quint32 ZXYYBLOCKS = 8;
    bool _useCache = true;
    BoundingBox _box;

    void loadFromCache(unsigned int y1, unsigned int z1, unsigned int z2, unsigned int x, int threadIndex);
    quint64 calcSeekPosition(unsigned int z, unsigned int y ) const;

    void loadFromSource(int z1, int x, int y1);
    bool pastHorizon(quint32 v) const;
    void sizeData(unsigned int y2, unsigned int z2, unsigned int& ysize, unsigned int& zsize, std::vector<char>& data);
};

class GridXYZ : public InternalGrid {
public:
    GridXYZ();
    void clone(quint64 newRasterId, int dLimit1, int dLimit2, InternalGrid *grid);
    PIXVALUETYPE &value(int z, int y, int x, int threadIndex);
    void setBlockData(int z, int y, const std::vector<PIXVALUETYPE>& data);
    void setValue(int z, int y, int x, PIXVALUETYPE v, int threadIndex );
    void setOrientation(Grid::Orientation ori);
    void unload(int threadIndex);
    void set3DSize(int z, int y, int xzs);
    quint32 linesPerBlock(int d) const;
    void closure();
    void dump(unsigned int yBase, unsigned int z1, unsigned int z2, int threadIndex) ;
    quint32 blocksCount() const;
    quint64 blockSize(int d) const;
    const InternalGrid::BlockStatus& blockStatus(int z, int blockIndex) const;
    quint64 seekPosition(int z, int y, int x) const;
    void loadFromCache(int y, int z, std::vector<char>& data);
    void setBlockStatus(unsigned int z1, quint64 ybase, bool status);
    quint32 blockCacheLimit() const;
    unsigned int blockFloor(unsigned int y) const;
    void setStartIndexes();
};

class GridZXY : public InternalGrid {
public:
    GridZXY();
    void clone(quint64 newRasterId, int dLimit1, int dLimit2, InternalGrid *grid);
    PIXVALUETYPE &value(int z, int y, int x, int threadIndex);
    void setBlockData(int z, int y, const std::vector<PIXVALUETYPE>& data);
    void setValue(int z, int y, int x, PIXVALUETYPE v, int threadIndex);
    void setOrientation(Grid::Orientation ori);
    void unload(int threadIndex);
    void set3DSize(int z, int y, int xzs);
    quint32 linesPerBlock(int d) const;
    void closure();
    void dump(unsigned int yBase, unsigned int z1, unsigned int z2, int threadIndex);
    quint32 blocksCount() const;
    quint64 blockSize(int d) const;
    const InternalGrid::BlockStatus& blockStatus(int z, int blockIndex) const;
    quint64 seekPosition(int z, int y, int x) const;
    void loadFromCache(int y, int z, std::vector<char>& data);
    void setBlockStatus( unsigned int z1, quint64 blockNr, bool status);
    quint32 blockCacheLimit() const;
    unsigned int blockFloor(unsigned int y) const;
    void setStartIndexes();
};

typedef std::unique_ptr<Grid> UPGrid;
}



#endif // Grid_H
