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

#include <QSqlQuery>
#include <QSqlError>
#include <QFile>
#include <QDir>
#include <fstream>
#include <iterator>
#include <future>

#include "kernel.h"
#include "raster.h"
#include "columndefinition.h"
#include "table.h"
#include "module.h"
#include "numericrange.h"
#include "connectorinterface.h"
#include "mastercatalog.h"
#include "ilwisobjectconnector.h"
#include "catalogexplorer.h"
#include "catalogconnector.h"
#include "inifile.h"
#include "numericrange.h"
#include "numericdomain.h"
#include "representation.h"
#include "catalog.h"
#include "ilwiscontext.h"
#include "pixeliterator.h"
#include "ilwis3connector.h"
#include "rawconverter.h"
#include "coverageconnector.h"
#include "tableconnector.h"
#include "operationhelpergrid.h"
#include "ilwisrastercoverageconnector.h"

using namespace Ilwis;
using namespace Ilwis3;

ConnectorInterface *RasterCoverageConnector::create(const Resource &resource, bool load, const IOOptions &options) {
    return new RasterCoverageConnector(resource, load, options);

}



RasterCoverageConnector::RasterCoverageConnector(const Resource &resource, bool load, const IOOptions &options) : CoverageConnector(resource, load, options),_storesize(1)
{
}

bool RasterCoverageConnector::loadMapList(IlwisObject *data,const IOOptions& options) {
    Ilwis3Connector::loadMetaData(data, options);

    RasterCoverage *gcoverage = static_cast<RasterCoverage *>(data);

    QString file = filename2FullPath(_odf->value("MapList","Map0"));
    if ( file == sUNDEF)
        return ERROR2(ERR_COULD_NOT_LOAD_2,"RasterCoverage",_odf->url());

    IRasterCoverage mp;
    if (!mp.prepare(file))
        return ERROR2(ERR_COULD_NOT_LOAD_2,"RasterCoverage",file);

    bool ok;
    qint32 z = _odf->value("MapList","Maps").toInt(&ok);
    Size<> sz = mp->size();
    sz.zsize(z);

    if (!ok || z < 0)
        return ERROR2(ERR_INVALID_PROPERTY_FOR_2,"Number of maps", gcoverage->name());

    gcoverage->gridRef()->prepare(gcoverage->id(), sz);
    std::vector<double> bands(sz.zsize());
    for(int i =0; i < sz.zsize(); ++i)
        bands[i] = i;
    gcoverage->stackDefinitionRef().setSubDefinition(IDomain("count"),bands);
    double gcvmin = rUNDEF, gcvmax = rUNDEF;
    for(int i = 0; i < z; ++i) {
        QString file = _odf->value("MapList",QString("Map%1").arg(i));
        //file = filename2FullPath(file);
        file = _resource.container(true).toLocalFile()+ "/" + file;
        if ( file != sUNDEF) {
            ODF odf(new IniFile(file));
            //QString dataFile = filename2FullPath(odf.value("MapStore","Data"));
            QUrl url (QUrl::fromLocalFile(_resource.container(true).toLocalFile() + "/" + odf->value("MapStore","Data")));
            _dataFiles.push_back(url);

            DataDefinition def = determineDataDefinition(odf, options);
            if ( !def.isValid()) {
                return false;
            }
            gcoverage->setBandDefinition(i, def);
            NumericRange* defnr = def.range()->as<NumericRange>();
            if (defnr) {
                gcvmin = (gcvmin != rUNDEF) ? std::min(gcvmin, defnr->min()) : defnr->min();
                gcvmax = (gcvmax != rUNDEF) ? std::max(gcvmax, defnr->max()) : defnr->max();
            }

        } else {
            ERROR2(ERR_COULD_NOT_LOAD_2,"files","maplist");
            --z;
        }
    }

    IniFile ini;
    if (!ini.setIniFile(QUrl(file).toLocalFile()))
        return ERROR2(ERR_COULD_NOT_LOAD_2,"files","maplist");

    gcoverage->datadefRef().domain(mp->datadef().domain<>());

    double vmax, vmin, scale, offset;
    //QString range = _odf->value("MapList", "Range"); // Do not use the range from the MapList ODF. It is different and cannot be used to compute the storetype's scale/offset
    QString range = ini.value("BaseMap", "Range"); // Get the range from the first map. This is for computing the storetype's scale/offset.
    if (range != sUNDEF) {
        if (getRawInfo(range, vmin, vmax, scale, offset)) {
            if (scale == 1.0) {
                gcoverage->datadefRef().range(new NumericRange(gcvmin, gcvmax, 1));
            }
            else {
                gcoverage->datadefRef().range(new NumericRange(gcvmin, gcvmax));
            }
        }
    }
    else if (gcvmin != rUNDEF && gcvmax != rUNDEF) {
        gcoverage->datadefRef().range(new NumericRange(gcvmin, gcvmax));
    }

    QString storeType = ini.value("MapStore","Type");
    setStoreType(storeType);
    const DataDefinition & def = gcoverage->datadef(0);
    if ( def.domain<>()->valueType() & itNUMBER)
        updateConverter(ini);

    gcoverage->georeference(mp->georeference());
    gcoverage->size(sz);
    gcoverage->coordinateSystem(mp->coordinateSystem());
    gcoverage->envelope(mp->envelope());
    _dataType = mp->datadef().range()->valueType();

    return true;

}

void RasterCoverageConnector::setStoreType(const QString& storeType) {
    _storetype = itUINT8;
    if ( storeType == "Int") {
        _storesize = 2;
        _storetype = itINT16;
    }
    else if ( storeType == "Long"){
        _storesize = 4;
        _storetype = itINT32;
    }
    else if ( storeType == "Float"){
        _storesize = 4;
        _storetype = itFLOAT;
    }
    else if ( storeType == "Real") {
        _storesize = 8;
        _storetype = itDOUBLE;
    }
    _converter.storeType(_storetype);
}

void RasterCoverageConnector::updateConverter(const IniFile & odf)
{
    QString dminfo = odf.value("BaseMap","DomainInfo");
    if ( dminfo != sUNDEF) {
        int index = dminfo.indexOf("class;");
        if ( index != -1) {
            _converter = RawConverter("class");
        }else {
            index = dminfo.indexOf("group;");
            if ( index != -1) {
                _converter = RawConverter("group");
            }else {
                index = dminfo.indexOf("id;");
                if ( index != -1) {
                    _converter = RawConverter("id");
                } else {
                    index = dminfo.indexOf("UniqueID;");
                    if ( index != -1) {
                        _converter = RawConverter("UniqueID");
                    }
                    index = dminfo.indexOf("color;");
                    if ( index != -1) {
                        _converter = RawConverter("color");
                    }
                }
            }
        }
    }
}

bool RasterCoverageConnector::setDataType(IlwisObject *data, const IOOptions &options) {

    RasterCoverage *raster = static_cast<RasterCoverage *>(data);

    DataDefinition def = determineDataDefinition(_odf, options);
    if ( !def.isValid()) {
        return false;
    }

    raster->datadefRef() = def;
    return true;
}

bool RasterCoverageConnector::loadMetaData(IlwisObject *data, const IOOptions &options)
{
    Locker<> lock(_mutex);

    QFileInfo inf(_resource.toLocalFile());

    bool isMapList  = inf.suffix().toLower() == "mpl";

    if (isMapList ){
        return loadMapList(data, options);
    }
    else {
        if (!setDataType(data, options))
            return false;
        _dataFiles.clear();

        if (!CoverageConnector::loadMetaData(data, options))
            return false;
    }

    RasterCoverage *gcoverage = static_cast<RasterCoverage *>(data);

    QString grfName = _odf->value("Map","GeoRef");
    grfName = filename2FullPath(grfName, _resource);
    IGeoReference grf;
    if (!grf.prepare(grfName, itGEOREF, options)) {
        kernel()->issues()->log(TR(ERR_COULDNT_CREATE_OBJECT_FOR_2).arg("Georeference",grfName));
        return false;
    }
	if (gcoverage->resourceRef().hasProperty("georeference")) {
		QString altUrl = gcoverage->resourceRef()["georeference"].toString();
		if (altUrl != grf->resourceRef().url(true).toString()) {
			grf.prepare(altUrl);
		}
	}
    if ( grf->code() == "undetermined"){ // grf none comes from the internal factory; it doesnt read some properties
        grf->name(data->name().left(data->name().indexOf(".")));
        QString szstr = _odf->value("Map","Size");
        QStringList parts = szstr.split(" ");
        if ( parts.size() == 2){
            Size<> sz(parts[1].toLong(), parts[0].toLong(),1);
            grf->size(sz);
        }
        grf->coordinateSystem(gcoverage->coordinateSystem());
    } else
        gcoverage->envelope(grf->envelope()); // always take over the envelope from the georeference (see ILWIS3, Map.cpp, line 834)
    QString dataFile = filename2FullPath(_odf->value("MapStore","Data"), _resource);
    if ( dataFile != sUNDEF)
         _dataFiles.push_back(dataFile);

    QString storeType = _odf->value("MapStore","Type");
    std::vector<double> v={0};
    gcoverage->stackDefinitionRef().setSubDefinition(IDomain("count"),v);
    gcoverage->setBandDefinition(0, gcoverage->datadef());

    setStoreType(storeType);
   // const DataDefinition & def = gcoverage->datadef();
   // if ( def.domain<>()->valueType() & itNUMBER)
   updateConverter(*_odf);

    gcoverage->gridRef()->prepare(gcoverage->id(), grf->size());

    gcoverage->georeference(grf);
    if ( !gcoverage->datadefRef().range().isNull() )
        _dataType = gcoverage->datadefRef().range()->valueType();
    else
        _dataType = itUNKNOWN;

    return true;

}

IlwisObject *RasterCoverageConnector::create() const
{
    return new RasterCoverage(_resource);
}

inline double RasterCoverageConnector::value(const char *block, int index) const{
    double v = rUNDEF;
    const char *c = &(block[index * _storesize]);
    switch (_storetype) {
    case itUINT8:
        v = *(quint8 *)c; break;
    case itINT16:
        v =  *(qint16 *)c; break;
    case itINT32:
        v = *(qint32 *)c; break;
    case itFLOAT:
        v = *(float *)c; break;
    case itINT64:
        return *(qint64 *)c;
    case itDOUBLE:
        v = *(double *)c; break;
    }
    return v;
}

/*qint64  RasterCoverageConnector::conversion(QFile& file, Grid *grid, int& count) {
    qint64 blockSizeBytes = grid->blockSize(0) * _storesize;
    qint64 szLeft = grid->size().xsize() * grid->size().ysize() * _storesize;
    qint64 result = 0;
    qint64 totalRead =0;
    char *block = new char[blockSizeBytes];
    bool noconversionneeded = _converter.isNeutral();
    while(szLeft > 0) {
        if ( szLeft >= blockSizeBytes)
            result = file.read((char *)block,blockSizeBytes);
        else {
            result = file.read((char *)block,szLeft);
        }
        if ( result == -1){
            kernel()->issues()->log(TR("Reading past the end of file %1").arg(file.fileName()));
            break;
        }
        quint32 noItems = grid->blockSize(count);
         if ( noItems == iUNDEF)
            return 0;
        vector<double> values(noItems);
        for(quint32 i=0; i < noItems; ++i) {
            double v = value(block, i);
            values[i] = noconversionneeded ? v :_converter.raw2real(v);
        }
        grid->setBlockData(count, values);
        totalRead += result;
        ++count;
        szLeft -= blockSizeBytes;

    }
    delete [] block;

    return totalRead;
}*/
vector<double> RasterCoverageConnector::loadBlock(UPGrid& grid,std::unique_ptr<QFile> &file, int y, quint32 z) {
    auto lines = grid->linesPerBlock(y);
    //quint64 blockNr = int((y / grid->blockCacheLimit()));
    auto blockSizeBytes = lines * grid->size().xsize() * _storesize;
    qint64 seekPos = y * grid->size().xsize() * _storesize ;

    if (file->seek(seekPos)) {
        QByteArray bytes = file->read(blockSizeBytes);
        quint32 noItems = lines * grid->size().xsize();
        if (noItems != iUNDEF) {
            vector<double> values(noItems);
            for (quint32 i = 0; i < noItems; ++i) {
                double v = value(bytes.constData(), i);
                if (_converter.isNeutral()) {

                    if (v != iILW3UNDEF && v != shILW3UNDEF)
                        values[i] = v;
                    else
                        values[i] = rUNDEF;
                }
                else
                    values[i] = _converter.raw2real(v);
            }
            return values;
        } // else we are trying to read beyond the available data, perhaps because a new band was added; just return
    }else
        ERROR2(ERR_COULD_NOT_OPEN_READING_2,file->fileName(),TR("seek failed"));
}

bool RasterCoverageConnector::openDataFile(int z, std::unique_ptr<QFile>& file){

    QString  datafile = _dataFiles[z].toLocalFile();
    if ( datafile.right(1) != "#") { // can happen, # is a special token in urls
        datafile += "#";
    }
    QFileInfo localfile(datafile);
    file.reset(new QFile(localfile.absoluteFilePath()));
    if ( !file->exists()){
        return ERROR1(ERR_MISSING_DATA_FILE_1,datafile);
    }
    if (!file->open(QIODevice::ReadOnly )) {
        return ERROR1(ERR_COULD_NOT_OPEN_READING_1,datafile);
    }
    return true;
}

bool RasterCoverageConnector::loadData(IlwisObject* data, const IOOptions &options)
{
    Locker<> lock(_mutex);
    IOOptions iooptions = options.isEmpty() ? ioOptions() : options;

    if ( _dataFiles.size() == 0) {
        return ERROR1(ERR_MISSING_DATA_FILE_1,_resource.name());
    }

    RasterCoverage *raster = static_cast<RasterCoverage *>(data);

    UPGrid& grid = raster->gridRef();
    if ( _openFiles.size() == 0){
        _openFiles.resize(grid->size().zsize());
        for(quint32 lz=0 ; lz < _openFiles.size(); ++lz){
            if (!openDataFile(lz, _openFiles[lz])){
                return false;
            }
        }
    }

    quint32 z = options.contains("z") ? options["z"].toUInt(): 0;
    quint32 y = options.contains("y") ? options["y"].toUInt(): 0;
    QString orientation = options.contains("orientation") ? options["orientation"].toString(): "XYZ";
    if (orientation == "ZXY"){

        std::vector<double> block;
        for(quint32 lz = 0; lz < _openFiles.size(); ++lz){
               auto values = loadBlock(grid, _openFiles[lz], y , lz);
               block.insert(block.end(), std::make_move_iterator(values.begin()), std::make_move_iterator(values.end()));

        }
        grid->setBlockData(z, y, block);
    }else {
           auto values = loadBlock(grid, _openFiles[z], y , z);
           grid->setBlockData(z, y, values);
    }
    if ( raster->attributeTable().isValid()) {
        ITable tbl = raster->attributeTable();
        QString primkey = tbl->primaryKey();
        if (primkey == sUNDEF)
            primkey = COVERAGEKEYCOLUMN; // fallback, should not happen?
        for(quint32 i=0; i < tbl->recordCount() ; ++i) {
            tbl->setCell(primkey,i, QVariant(i));
        }
        raster->primaryKey(primkey);
    }
    _binaryIsLoaded = true;
    return true;

}

bool RasterCoverageConnector::storeBinaryData(IlwisObject *obj)
{
    Locker<> lock(_mutex);

    if ( obj == nullptr)
        return false;
    IRasterCoverage raster = mastercatalog()->get(obj->id());
    if ( !raster.isValid())
        return false;

    if (!raster->georeference().isValid())
        return false;
    if ( raster->size().zsize() > 1) // mpl doesnt have binary data
        return true;

    const IDomain dom = raster->datadef().domain<>();
    if (!dom.isValid())
        return ERROR2(ERR_NO_INITIALIZED_2, "Domain", raster->name());

    QFileInfo inf(obj->resource(IlwisObject::cmOUTPUT).toLocalFile());
    QString filename;

    filename = inf.absolutePath() + "/" + QString(inf.baseName()).replace(QRegExp("[/ .'\"]"),"_") + ".mp#";

    Size<> sz = raster->size();
    bool ok = false;
    if ( dom->ilwisType() == itNUMERICDOMAIN) {
        calcStatistics(obj, NumericStatistics::pBASIC);
        const NumericStatistics& stats = raster->statistics(PIXELVALUE);
        double resolution = raster->datadef().range()->as<NumericRange>()->resolution();
		double precision = resolution;
        if (precision < 1e-06)
            precision = 0.0;
        bool hasUndefs = stats[NumericStatistics::pCOUNT] != stats[NumericStatistics::pNETTOCOUNT];
        if (_domainName == "image") // image does not support undef
            hasUndefs = false; // nothing else needed, as long as conv.real2raw() converts rUNDEF to 0
        RawConverter conv(stats[NumericStatistics::pMIN], stats[NumericStatistics::pMAX], precision, hasUndefs);

        std::ofstream output_file(filename.toLatin1(),ios_base::out | ios_base::binary | ios_base::trunc);
        if ( !output_file.is_open())
            return ERROR1(ERR_COULD_NOT_OPEN_WRITING_1,filename);

        qint32 delta = stats[NumericStatistics::pDELTA];
        if ( delta >= 0 && delta < (hasUndefs ? 255 : 256) && resolution == 1 && stats[NumericStatistics::pMIN] >= 0){ // if there is an undef, restrict to one less than 255, to make space for the undef
            ok = save<quint8>(output_file,(conv.scale() == 1 && conv.offset() == 0) ? RawConverter() : conv, raster,sz);
        } else if ( conv.storeType() == itUINT8) {
            ok = save<quint8>(output_file,(conv.scale() == 1 && conv.offset() == 0) ? RawConverter() : conv, raster,sz);
        } else if ( conv.storeType() == itINT16) {
            ok = save<qint16>(output_file,conv, raster,sz);
        } else if ( conv.storeType() == itINT32) {
            ok = save<qint32>(output_file,conv, raster,sz);
        } else { // itDOUBLE
            ok = save(output_file, raster, sz);
        }
        output_file.close();

    } else if ( dom->ilwisType() == itITEMDOMAIN ){
        if ( hasType(dom->valueType(), itTHEMATICITEM | itNAMEDITEM | itNUMERICITEM | itINDEXEDITEM)) {
            std::ofstream output_file(filename.toLatin1(),ios_base::out | ios_base::binary | ios_base::trunc);
            if ( !output_file.is_open()){
                return ERROR1(ERR_COULD_NOT_OPEN_WRITING_1,filename);
            }

            if( hasType(dom->valueType(), itTHEMATICITEM | itNUMERICITEM)){
                RawConverter conv(dom->valueType() == itTHEMATICITEM ? "class" : "group");
                ok = save<quint8>(output_file,conv, raster,sz);
            }
            else{
                RawConverter conv("ident");
                ok = save<quint16>(output_file,conv, raster,sz);
            }
        }
    }
    ITable attTable = raster->attributeTable();
    if ( attTable.isValid() && attTable->isValid()) {
         attTable->store({"storemode",IlwisObject::smBINARYDATA});
    }
    return ok;

}

QString RasterCoverageConnector::format() const
{
    return "map";
}

void RasterCoverageConnector::calcStatistics(const IlwisObject *obj, NumericStatistics::PropertySets set) const {
    IRasterCoverage raster = mastercatalog()->get(obj->id());
    if ( !raster->statistics(PIXELVALUE).isValid()) {
        PixelIterator iter(raster,BoundingBox(raster->size()));
        raster->statisticsRef(PIXELVALUE).calculate(iter, iter.end(),set);
    }
}

QString RasterCoverageConnector::createBandName(const IRasterCoverage& raster, quint32 i){
    QString base = sourceRef().name();
    QString band = raster->stackDefinition().index(i);
    return base + "_band_" + band;
}

bool RasterCoverageConnector::storeMetaDataMapList(IlwisObject *obj) {
    bool ok = Ilwis3Connector::storeMetaData(obj, itRASTER | itCOLLECTION);
    if ( !ok)
        return false;


    IRasterCoverage raster = mastercatalog()->get(obj->id());

    QString localName = getGrfName(raster);
    if ( localName == sUNDEF)
        return false;
    _odf->setKeyValue("Ilwis","Type","MapList");
    _odf->setKeyValue("MapList","GeoRef",QFileInfo(localName).fileName());
    Size<> sz = raster->size();
    _odf->setKeyValue("MapList","Size",QString("%1 %2").arg(sz.ysize()).arg(sz.xsize()));
    _odf->setKeyValue("MapList","Maps",IniFile::FormatElement(sz.zsize()));

    QString mpl_name = obj->name();
    for(quint32 i = 0; i < sz.zsize(); ++i) {
        QString mapName = createBandName(raster, i);
        //QString mapName = QString("%1_band_%2").arg(obj->name()).arg(i);
        mapName = mapName.replace(QRegExp("[/ .'\"]"),"_");
        _odf->setKeyValue("MapList",QString("Map%1").arg(i),mapName + ".mpr");

        IRasterCoverage gcMap;
        OperationHelperRaster::initialize(raster, gcMap, itENVELOPE | itGEOREF | itCOORDSYSTEM | itDOMAIN);
        gcMap->size(Size<>(sz.xsize(), sz.ysize(),1));
        gcMap->name(mapName);

        int index = _odf->url().lastIndexOf("/");
        gcMap->copyBinary(raster, i);
        QString path = _odf->url().left(index);
        QUrl url =  path + "/" + mapName;
        gcMap->connectTo(url, "map", "ilwis3", Ilwis::IlwisObject::cmOUTPUT);
        IOOptions options;
        options.addOption({"storemode",IlwisObject::smBINARYDATA | IlwisObject::smMETADATA});
        options.addOption({"collection", mpl_name});
        gcMap->store(options);
    }

    _odf->store("mpl",sourceRef().toLocalFile());
    return true;
}

QString RasterCoverageConnector::getGrfName(const IRasterCoverage& raster) {
    const IGeoReference grf = raster->georeference();
    if (!grf.isValid()) {
        ERROR2(ERR_NO_INITIALIZED_2, "Georeference", raster->name());
        return sUNDEF;
    }
    if ( grf->code() == "undetermined"){
        return "none.grf";
    }
    QString name = grf->resource(IlwisObject::cmOUTPUT).url().toString();
    if ( grf->isAnonymous()) { // get a suitable output name
        name = raster->resource(IlwisObject::cmOUTPUT).url().toString();
    }
    QString localName = Resource::toLocalFile(QUrl(name),false, "grf");
    if ( localName == sUNDEF){
        QFileInfo localGrfPath(_resource.url(true).toLocalFile());
        QString namePart = name.mid(name.lastIndexOf("/"));
        localName = localGrfPath.absolutePath() + namePart;

    }
    QString destinationPath = QFileInfo(sourceRef().toLocalFile()).absolutePath();
    QFileInfo sourceGrf(localName);
    QString fileName (sourceGrf.fileName());
    QFileInfo destGrf(destinationPath + '/' + fileName);

    if ( !destGrf.exists()) { // if it is not an existing ilwis3 grf, we create one from scratch
        QUrl url = QUrl::fromLocalFile(destGrf.absoluteFilePath());
        grf->connectTo(url, "georef", "ilwis3", Ilwis::IlwisObject::cmOUTPUT);
        grf->store({"storemode", IlwisObject::smMETADATA});
    }

    return destGrf.absoluteFilePath();
}

bool RasterCoverageConnector::storeMetaData( IlwisObject *obj, const IOOptions& options)  {
    Locker<> lock(_mutex);

    IRasterCoverage raster = mastercatalog()->get(obj->id());
    if (!raster.isValid())
        return false;
    if (!raster->georeference().isValid())
        return false;

    if ( raster->size().zsize() > 1)
        return storeMetaDataMapList(obj);

    bool ok = CoverageConnector::storeMetaData(obj, itRASTER, raster->datadef());
    if ( !ok)
        return false;

    _odf->setKeyValue("BaseMap","Type","Map");


    if ( !raster.isValid())
        return ERROR2(ERR_COULD_NOT_LOAD_2,"RasterCoverage", obj->name());

    QString localName = getGrfName(raster);
    if ( localName == sUNDEF)
        return false;
    if ( raster->georeference()->code() == "undetermined"){
        Envelope bounds = raster->envelope();
        if ( bounds.isNull() || !bounds.isValid())
            bounds = raster->coordinateSystem()->envelope();
        _odf->setKeyValue("BaseMap","CoordBounds",QString("%1 %2 %3 %4").
                          arg(bounds.min_corner().x,0,'f',10).
                          arg(-bounds.max_corner().y,0,'f',10).
                          arg(bounds.max_corner().x,0,'f',10).
                          arg(-bounds.min_corner().y,0,'f',10));
    }

    QFileInfo grfInf(localName);

    _odf->setKeyValue("Map","GeoRef",grfInf.fileName());
    Size<> sz = raster->size();
    _odf->setKeyValue("Map","Size",QString("%1 %2").arg(sz.ysize()).arg(sz.xsize()));
    _odf->setKeyValue("Map","Type","MapStore");

    const IDomain dom = raster->datadef().domain<>();
    if ( dom->ilwisType() == itNUMERICDOMAIN) {
        const NumericStatistics& stats = raster->statistics(PIXELVALUE);
        double resolution = raster->datadef().range()->as<NumericRange>()->resolution();
        //double precision = (resolution == 0.0) ? resolution : pow(10, -stats.significantDigits());
		double precision = resolution;
        if (precision < 1e-06)
            precision = 0.0;
        bool hasUndefs = stats[NumericStatistics::pCOUNT] != stats[NumericStatistics::pNETTOCOUNT];
        if (_domainName == "image") // image does not support undef
            hasUndefs = false;
        RawConverter conv(stats[NumericStatistics::pMIN], stats[NumericStatistics::pMAX], precision, hasUndefs);
        qint32 delta = stats[NumericStatistics::pDELTA];
        QString minmax;
        if (_domainName == "image")
            minmax = QString("%1:%2").arg(dom->range<NumericRange>()->min()).arg(stats[NumericStatistics::pMAX]);
        else
            minmax = QString("%1:%2").arg(stats[NumericStatistics::pMIN]).arg(stats[NumericStatistics::pMAX]);
        _odf->setKeyValue("BaseMap","MinMax",minmax);

        // Add the BaseMap:Range section to match the default valuerange with the Type
        QString range;
        if ( delta >= 0 && delta < (hasUndefs ? 255 : 256) &&  resolution == 1 && stats[NumericStatistics::pMIN] >= 0){ // if there is an undef, restrict to one less than 255, to make space for the undef
           _odf->setKeyValue("MapStore","Type","Byte");
           range = "0:255:offset=0";
        } else if ( conv.storeType() == itUINT8){
           _odf->setKeyValue("MapStore","Type","Byte");
           range = "0:255:offset=0";
        } else if ( conv.storeType() == itINT16){
            _odf->setKeyValue("MapStore","Type","Int");
            range = "-32768:32767:offset=0";
        } else if ( conv.storeType() == itINT32){
            _odf->setKeyValue("MapStore","Type","Long");
            range = "-2147483648:2147483647:offset=0";
        } else if ( conv.storeType() == itDOUBLE){
            _odf->setKeyValue("MapStore","Type","Real");
            range = "-1e300:1e300:0:offset=0";
        }
        // if stats where not calculated (properly) at least specify the default range for the data type
        if (!stats.isValid())
            _odf->setKeyValue("BaseMap", "Range", range);

        if(!isSystemObject(_domainName) && _domainName != "value.dom") { // check!!
            QString filename = context()->workingCatalog()->resolve(_domainName);
            if ( filename == sUNDEF){
                int index = _odf->url().lastIndexOf("/");
                filename = _odf->url().left(index) + "/" + _domainName;
            }
            dom->connectTo(filename,"domain","ilwis3", Ilwis::IlwisObject::cmOUTPUT);
            dom->store();
        }
    } if ( hasType(dom->ilwisType(),itITEMDOMAIN)) {
        if ( hasType(dom->valueType(), itTHEMATICITEM | itNUMERICITEM)  )
            _odf->setKeyValue("MapStore","Type","Byte");
        else if ( hasType(dom->valueType(), itNAMEDITEM | itINDEXEDITEM)) {
            _odf->setKeyValue("MapStore","Type","Int");
        }
        if ( _domainName.indexOf(".dom") != -1 && !isSystemObject(_domainName)){
            QString filename = context()->workingCatalog()->resolve(_domainName);
            if ( filename == sUNDEF){
                int index = _odf->url().lastIndexOf("/");
                filename = _odf->url().left(index) + "/" + _domainName;
            }
            dom->connectTo(filename,"domain","ilwis3", Ilwis::IlwisObject::cmOUTPUT);
            dom->store();

            if (hasType(dom->valueType(), itTHEMATICITEM | itNUMERICITEM)) {
                const IRepresentation rpr = raster->datadef().representation();
                if (rpr.isValid()) {
                    QString rprName = !rpr->isAnonymous() ? rpr->name() : (QFileInfo(_domainName).baseName() + ".rpr");
                    QString rprFilename = context()->workingCatalog()->resolve(rprName);
                    if (rprFilename == sUNDEF) {
                        int index = _odf->url().lastIndexOf("/");
                        rprFilename = _odf->url().left(index) + "/" + rprName;
                    }
                    rpr->connectTo(rprFilename, "representation", "ilwis3", Ilwis::IlwisObject::cmOUTPUT);
                    rpr->store();
                }
            }
        }
    }

    ITable attTable = raster->attributeTable();
    if ( attTable.isValid() && attTable->columnCount() > 1) {
        QFileInfo basename(QUrl(_odf->url()).toLocalFile());
        //QScopedPointer<TableConnector> conn(createTableStoreConnector(attTable, raster.ptr(), itRASTER,basename.baseName()));
        //conn->storeMetaData(attTable.ptr(),options);
        createTableStoreConnector(attTable, raster.ptr(), itRASTER,basename.baseName());
        attTable->store({"storemode",IlwisObject::smMETADATA});
    }

    QFileInfo inf(_resource.toLocalFile());
    QString file = inf.baseName() + ".mp#";
    QString exts = "mprmpamppmpsdomtbtgrfcsympl";
    if ( exts.indexOf(inf.suffix()) == -1)
        file = inf.fileName() + ".mp#";

    _odf->setKeyValue("MapStore","Data",file);
    _odf->setKeyValue("MapStore","Structure","Line");
    _odf->setKeyValue("MapStore","StartOffset","0");
    _odf->setKeyValue("MapStore","RowLength",IniFile::FormatElement(sz.xsize()));
    _odf->setKeyValue("MapStore","PixelInterLeaved","No");
    _odf->setKeyValue("MapStore","SwapBytes","No");
    _odf->setKeyValue("MapStore","UseAs","No");

    QString coll_name = sUNDEF;
    if (options.contains("collection"))
        coll_name = options["collection"].toString();
    if (coll_name != sUNDEF){
        _odf->setKeyValue("Collection", "NrOfItems", "1");
        _odf->setKeyValue("Collection", "Item0", coll_name);
    }

    _odf->store("mpr", sourceRef().toLocalFile());

    return true;
}

