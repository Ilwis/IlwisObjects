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

namespace Ilwis {

class RasterCoverage;
class IOOptions;


typedef std::vector<std::vector<double>> Plane;

class KERNELSHARED_EXPORT Grid

{
public:
    typedef std::vector<std::vector<double>> Plane;
    enum Orientation {oXYZ, oZXY};
    friend class PixelIterator;
    friend class GridBlock;

    Grid(Orientation ori=oXYZ);
    virtual ~Grid();

    void clear();


    PIXVALUETYPE &valueRef(const Pixel& pix, int threadIndex = 0);
    PIXVALUETYPE value(const Pixel& pix, int threadIndex = 0);
    void setValue(int z, int y, int x, PIXVALUETYPE v , int threadIndex=0);
    void setOrientation(Grid::Orientation ori);

    void setBlockData(int d1, int d2, const std::vector<PIXVALUETYPE>& data);
    bool prepare(quint64 rasterid, const Size<> &sz) ;
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


protected:

private:
    struct BlockStatus{
        bool _valid = false;
        bool _loadedFromSource = false;
    };
    typedef std::vector<std::deque<quint64>> PlaneBlockCache;

    void dumpBlock(int z, int y, int x, int threadIndex);
    bool createCacheFile(int i);
    //void getBlockData(int d1, int d2, std::vector<PIXVALUETYPE>& data);
    PIXVALUETYPE& value(int d1, int d2, int d3, int threadIndex = 0); //d1=z,d2=y,d3=x (XYZ) or d1=y, d2=x, d3=z (ZXY)
    void setBlockDataZXY(int z, int y, const std::vector<PIXVALUETYPE>& data);
    void setBlockDataXYZ(int z, int y, const std::vector<PIXVALUETYPE>& data);
    const BlockStatus& blockStatus(int z, int y) const;
    quint64 seekPosition(int z, int y, int x) const;
    quint32 normY(int y) const;
    void loadFromSource(int z, int x, int y)    ;
    bool pastHorizon(quint32 v) const ;
    void closure();
    void dump(unsigned int y1, unsigned int y2, unsigned int z1, unsigned int z2, int threadIndex) ;
    void load(unsigned int y1, unsigned int y2, unsigned int z1, unsigned int z2, unsigned int x, int threadIndex);

    const unsigned int BLOCKYHORIZON = 4;
    const unsigned int BLOCKXHORIZON = 8;
    std::vector<Plane> _planes;
    std::vector<std::vector<BlockStatus>> _validStripe;
    std::recursive_mutex _mutex;
    Size<> _size;
    Orientation _orientation;
    quint64 _rasterid;
    std::vector<QFile*> _diskImages;
    std::vector<QString> _imageNames;
    const quint32 XYZYBLOCKS = 100;
    const quint32 ZXYYBLOCKS = 8;
    bool _useCache = true;

    void loadFromCache(unsigned int y1, unsigned int y2, unsigned int z1, unsigned int z2, unsigned int x, int threadIndex);
    void setBlockStatus(unsigned int z1, quint64 blockNr, bool status);
    quint64 calcSeekPosition(unsigned int z, unsigned int y ) const;
};

typedef std::unique_ptr<Grid> UPGrid;
}



#endif // Grid_H
