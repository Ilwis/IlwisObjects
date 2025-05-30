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
#include "version.h"
#include "connectorinterface.h"
#include "versionedserializer.h"
#include "domain.h"
#include "table.h"
#include "basetable.h"
#include "flattable.h"
#include "pixeliterator.h"
#include "grid.h"
#include "factory.h"
#include "abstractfactory.h"
#include "versioneddatastreamfactory.h"
#include "ilwisobjectconnector.h"
#include "streamconnector.h"
#include "coverageserializerv1.h"
#include "rawconverter.h"
#include "rasterserializerv1.h"

using namespace Ilwis;
using namespace Stream;

RasterSerializerV1::RasterSerializerV1(QDataStream& stream, const QString &version) : CoverageSerializerV1(stream, version)
{
}

template<typename T> void storeBulk(const RawConverter& converter, QDataStream& stream, StreamConnector *streamconnector, const BoundingBox& box, const IRasterCoverage& raster){
    quint64 count = streamconnector->position();
    if ( streamconnector->isFileBased()){
        const UPGrid& grid = raster->grid();
        std::vector<T> rawData(1000 * grid->size().xsize());
        for(quint32 z =0; z < grid->size().zsize(); ++z){
            quint32 yrel = 0;
            quint32 pos = 0;
            for(quint32 y =0; y < grid->size().ysize(); ++y){
                for(quint32 x =0; x < grid->size().xsize(); ++x){
                    pos = x + grid->size().xsize() * yrel;
                    rawData[pos ] = converter.real2raw(grid->value(Pixel(x,y,z)));
                }
                if ( pos == rawData.size() - 1) {
                    stream.writeRawData((const char *)rawData.data(), rawData.size() * sizeof(T));
                    yrel = 0;
                }else
                    yrel = y;
            }

        }
        streamconnector->flush(true);
    }else {
        PixelIterator iter(raster, box);
        while(iter != iter.end()){
            if ( count >= STREAMBLOCKSIZE - 9 ) {
                streamconnector->flush(false);
                count = 0;

            }
            count += sizeof(T);
            stream << (T)converter.real2raw(*iter);
            ++iter;
        }
    }
}

template<typename T> void loadBulk(std::vector<T>& rawdata, std::vector<PIXVALUETYPE>& realdata,const RawConverter& converter, QDataStream& stream, StreamConnector *streamconnector, const BoundingBox& box, IRasterCoverage& raster){

    if ( streamconnector->isFileBased()){

        UPGrid &grid = raster->gridRef();
        std::vector<T> rawdata(1000 * grid->size().xsize());
        int maxblocks = grid->blocksCount();
        for ( int b=0; b < maxblocks; ++b){
            int count = stream.readRawData((char *)&rawdata[0],rawdata.size() * sizeof(T) );
            for(int p=0; p < count; ++p){
                realdata[p] = converter.raw2real(rawdata[p]);
            }
        }



        streamconnector->flush(true);
    } else {
     // TODO network restore, roughly code below though that was based on a older version
        //    UPGrid &grid = raster->gridRef();
        //    qint64 blockSizeBytes = grid->blockSize(block) * sizeof(T);
        //    qint64 szLeft = data.size();
        //    T value;
        //    QBuffer buf(&data);
        //    buf.open(QIODevice::ReadWrite);
        //    QDataStream stream(&buf);
        //    std::vector<PIXVALUETYPE> values(grid->blockSize(0));
        //    quint32 noItems = grid->blockSize(block);
        //    if ( noItems == iUNDEF)
        //        return 0;

        //    values.resize(noItems);
        //    for(quint32 i=0; i < noItems; ++i) {
        //        stream >> value;

        //        values[i] = converter.raw2real(value);
        //    }

        //    grid->setBlockData(block, values);
        //    szLeft -= blockSizeBytes;
    }
}

bool RasterSerializerV1::store(IlwisObject *obj, const IOOptions &opt)
{
	auto options = addParent(obj, opt);
    if (!CoverageSerializerV1::store(obj, options))
        return false;
    RasterCoverage *raster = static_cast<RasterCoverage *>(obj);


    VersionedDataStreamFactory *factory = kernel()->factory<VersionedDataStreamFactory>("ilwis::VersionedDataStreamFactory");
    if (!factory)
        return false;

    _stream << raster->size().xsize() << raster->size().ysize() << raster->size().zsize();


    if(!storeDataDefintion(raster->datadef(), _stream, options))
        return false;

    for(int index = 0; index < raster->size().zsize(); ++index)   {
        const DataDefinition& def = raster->datadef(index);
        storeDataDefintion(def,_stream, options);

    }

    std::unique_ptr<DataInterface> domainStreamer(factory->create(Version::interfaceVersion41, itDOMAIN,_stream));
    if ( !domainStreamer)
        return false;
    auto vtype = raster->stackDefinition().domain()->valueType();
    _stream << vtype;
    storeSystemPath(raster->stackDefinition().domain()->resource());
    domainStreamer->store( raster->stackDefinition().domain().ptr(), options);

    std::vector<QString> indexes = raster->stackDefinition().indexes();
    _stream << (quint32)indexes.size();
    for(auto index : indexes)
        _stream << index;

    std::unique_ptr<DataInterface> grfstreamer(factory->create(Version::interfaceVersion41, itGEOREF,_stream));
    if ( !grfstreamer)
        return false;

    storeSystemPath(raster->georeference()->resource());
	IOOptions newOptions;
	newOptions.addOption("storename", raster->name());
    if(!grfstreamer->store(raster->georeference().ptr(), newOptions))
        return false;
    _stream << raster->hasAttributes();
    if ( raster->hasAttributes()){
        std::unique_ptr<DataInterface> tblstreamer(factory->create(Version::interfaceVersion41, itTABLE,_stream));
        if ( !tblstreamer)
            return false;

        if(!tblstreamer->store(raster->attributeTable().ptr(), options))
            return false;
		if (!tblstreamer->storeData(raster->attributeTable().ptr(), options))
			return false;
        _stream << raster->primaryKey();
    }

    return true;

}

bool RasterSerializerV1::storeData(IlwisObject *obj, const IOOptions &options )
{
    qint64 pos = _stream.device()->pos();
    _stream << pos + sizeof(qint64);
    _stream << itRASTER;
    _stream << Version::interfaceVersion40;
    RasterCoverage *raster = static_cast<RasterCoverage *>(obj);
    RawConverter converter;
    if ( hasType(raster->datadef().domain()->ilwisType() , itNUMERICDOMAIN)){
        NumericStatistics& stats = raster->statistics(PIXELVALUE, ContainerStatistics<PIXVALUETYPE>::pBASIC);
		PIXVALUETYPE scale = raster->datadef().range()->as<NumericRange>()->resolution();
		bool hasUndefs = stats[ContainerStatistics<PIXVALUETYPE>::pCOUNT] != stats[ContainerStatistics<PIXVALUETYPE>::pNETTOCOUNT];
        converter = RawConverter(stats[ContainerStatistics<PIXVALUETYPE>::pMIN], stats[ContainerStatistics<PIXVALUETYPE>::pMAX],scale, hasUndefs);

        _stream << stats[ContainerStatistics<PIXVALUETYPE>::pMIN] << stats[ContainerStatistics<PIXVALUETYPE>::pMAX] << scale;
		_stream << (quint32) hasUndefs;
    }else{
        if(hasType(raster->datadef().domain()->ilwisType() ,itITEMDOMAIN) )
            converter = RawConverter("ident");
        if(hasType(raster->datadef().domain()->ilwisType() ,itCOLORDOMAIN) )
            converter = RawConverter("color");

    }
    if ( !converter.isValid()){
        kernel()->issues()->log(QString(TR("Couldnt find a correct converter for raster data of %1")).arg(obj->name()));
        return false;
    }
    BoundingBox box;
    if (options.contains("lines")) {
        QStringList parts = options["lines"].toString().split(" ");
        quint32 layer = parts[0].toUInt();
        quint32 minlines = parts[1].toUInt();
        quint32 maxlines = parts[2].toUInt();
        _stream <<  layer <<minlines << maxlines;
        box = BoundingBox(Pixel(0, minlines,layer), Pixel(raster->size().xsize(),maxlines,layer));

    }else {
        quint32 undef = iUNDEF;
        _stream << undef << undef << undef;
    }
    IRasterCoverage rcoverage(raster);
    switch (converter.storeType()){
    case itUINT8:
        storeBulk<quint8>(converter, _stream, _streamconnector, box,rcoverage); break;
    case itINT16:
        storeBulk<qint16>(converter, _stream, _streamconnector, box, rcoverage); break;
    case itUINT16:
        storeBulk<quint16>(converter, _stream, _streamconnector, box, rcoverage); break;
    case itINT32:
        storeBulk<qint32>(converter, _stream, _streamconnector, box, rcoverage); break;
    case itUINT32:
        storeBulk<quint32>(converter, _stream, _streamconnector, box, rcoverage); break;
    case itDOUBLE:
        storeBulk<double>(converter, _stream, _streamconnector, box, rcoverage); break;
	case itFLOAT:
		storeBulk<float>(converter, _stream, _streamconnector, box, rcoverage); break;
        break;
    case itINT64:
    default:
        for(PIXVALUETYPE v : rcoverage)
            _stream << (qint64)v;
        break;
    }

    return true;
}

bool RasterSerializerV1::loadMetaData(IlwisObject *obj, const IOOptions &options)
{
    if (!CoverageSerializerV1::loadMetaData(obj, options))
        return false;
    VersionedDataStreamFactory *factory = kernel()->factory<VersionedDataStreamFactory>("ilwis::VersionedDataStreamFactory");
    if (!factory)
        return false;

    RasterCoverage *raster = static_cast<RasterCoverage *>(obj);


    quint32 xsize, ysize, zsize;
    _stream >> xsize >> ysize >> zsize;
    raster->size(Size<>(xsize, ysize, zsize));


    loadDataDefinition(raster->datadefRef(),_stream, options);
    for(int band = 0; band < raster->size().zsize(); ++band) {
        loadDataDefinition(raster->datadefRef(band), _stream, options)    ;
    }
    IlwisTypes valueType;
    _stream >> valueType;
    quint64 type;
    QString version, url;
    _stream >> url;
    _stream >> type;
    _stream >> version;

    std::unique_ptr<DataInterface> domainStreamer(factory->create(version, itDOMAIN,_stream));
    if ( !domainStreamer)
        return false;

    IDomain systemDomain = makeSystemObject<IDomain>(url);
    IDomain dom(type|valueType);
    domainStreamer->loadMetaData( dom.ptr(), options);
    quint32 nrOfBands;
    _stream >> nrOfBands;
    std::vector<QString> variants(nrOfBands);
    for(int i =0; i < nrOfBands; ++i){
        _stream >> variants[i];
    }
    raster->stackDefinitionRef().setSubDefinition(systemDomain.isValid() ? systemDomain : dom, variants);


    _stream >> url;
    _stream >> type;
    _stream >> version;

    std::unique_ptr<DataInterface> grfstreamer(factory->create(version, itGEOREF,_stream));
    if ( !grfstreamer)
        return false;
    IGeoReference systemGrf = makeSystemObject<IGeoReference>(url);
    IGeoReference georeference (type);
    grfstreamer->loadMetaData(georeference.ptr(), options)    ;
	georeference->resourceRef().addContainer(raster->resourceRef().container(), true);
	georeference->resourceRef().addContainer(raster->resourceRef().container());
    raster->georeference(systemGrf.isValid() ? systemGrf : georeference);


    bool hasAttr;
    _stream >> hasAttr;
    if ( hasAttr){
        _stream >> type;
        _stream >> version;
        VersionedDataStreamFactory *factory = kernel()->factory<VersionedDataStreamFactory>("ilwis::VersionedDataStreamFactory");
        std::unique_ptr<DataInterface> tableStreamer(factory->create(version, itTABLE,_stream));
        if ( !tableStreamer)
            return false;
		static_cast<VersionedSerializer *>(tableStreamer.get())->connector(_streamconnector);
        ITable tbl;
        tbl.prepare();

        tableStreamer->loadMetaData(tbl.ptr(),options);
        _stream >> type;
        _stream >> version;
        tableStreamer->loadData(tbl.ptr(),options);
        QString primkey;
        _stream >> primkey;
          raster->primaryKey(primkey);
        tbl->resourceRef().setExtendedType(itRASTER);
        tbl->resourceRef().addProperty("rasterid", raster->id());
        raster->setAttributes(tbl);
    }
    qint64 beginData;
    _stream >> beginData;
    _streamconnector->beginDataSection(itRASTER, beginData);

    return true;
}


bool RasterSerializerV1::loadData(IlwisObject *data, const IOOptions &options)
{
    RasterCoverage *raster = static_cast<RasterCoverage *>(data);
    BoundingBox box;
    RawConverter converter;
    if ( hasType(raster->datadef().domain()->ilwisType(), itNUMERICDOMAIN)){
		quint32 hasUndefs;
        PIXVALUETYPE mmin, mmax, mscale;
        _stream >> mmin >> mmax >> mscale;
		_stream >> hasUndefs;
        converter = RawConverter(mmin, mmax, mscale, (bool)hasUndefs);
    }else {
        if (  hasType(raster->datadef().domain()->ilwisType(), itITEMDOMAIN))
            converter = RawConverter("ident");
        if ( hasType(raster->datadef().domain()->ilwisType(), itCOLORDOMAIN))
            converter = RawConverter("color");
    }
    if ( !converter.isValid()){
        kernel()->issues()->log(QString(TR("Couldnt find a correct converter for raster data of %1")).arg(data->name()));
        return false;
    }

    quint32 layerIndex, minLines, maxLines; //only defined in some cases, if not defined it assumed that the whole coverage is there
    _stream >> layerIndex >> minLines >> maxLines;
    Resource resource = data->resource();
    if ( resource.code().indexOf("band=") == 0){
        bool ok;
        QString part = resource.code().mid(5);
        int band = part.toInt(&ok);
        if ( ok){
            Size<> sz = raster->size();
            box = BoundingBox(Pixel(0,0,band),Pixel(sz.xsize(),sz.ysize(),band));
        }
	}
	/*else if (options.contains("blockindex")) {
		double currentBlock = options["blockindex"].toInt();
		int band = currentBlock / raster->grid()->blocksPerBand();
		int relativeBlock = currentBlock - band * raster->grid()->blocksPerBand();
		unsigned int minLine = raster->grid()->maxLines() * relativeBlock;
		unsigned int maxLine = std::min(minLine + raster->grid()->maxLines(), raster->size().ysize());
		box = BoundingBox(Pixel(0, minLine, band), Pixel(raster->size().xsize(), maxLine, band));
	}*/
	std::vector<PIXVALUETYPE> realdata;
    IRasterCoverage rcoverage(raster);
    switch (converter.storeType()){
	case itUINT8: {
		std::vector<quint8> rawdata;
		loadBulk<quint8>(rawdata, realdata, converter, _stream, _streamconnector, box, rcoverage); break;
	}
	case itINT16: {
		std::vector<qint16> rawdata;
		loadBulk<qint16>(rawdata, realdata,converter, _stream, _streamconnector, box, rcoverage); break;
	}
	case itUINT16: {
		std::vector<quint16> rawdata;
		loadBulk<quint16>(rawdata, realdata, converter, _stream, _streamconnector, box, rcoverage); break;
	}
	case itINT32: {
		std::vector<qint32> rawdata;
		loadBulk<qint32>(rawdata, realdata, converter, _stream, _streamconnector, box, rcoverage); break;
	}
	case itUINT32: {
		std::vector<qint32> rawdata;
		loadBulk<qint32>(rawdata, realdata, converter, _stream, _streamconnector, box, rcoverage); break;
	}
	case itDOUBLE: {
		std::vector<double> rawdata;
		loadBulk<double>(rawdata, realdata, converter, _stream, _streamconnector, box, rcoverage); break;
	}
	case itFLOAT: {
		std::vector<float> rawdata;
		loadBulk<float>(rawdata, realdata, converter, _stream, _streamconnector, box, rcoverage); break;
	}
    case itINT64:
	default: {
		std::vector<qint64> rawdata;
		loadBulk<qint64>(rawdata, realdata, converter, _stream, _streamconnector, box, rcoverage); break;
	}
    }
    _dataLoaded = true;
    return true;
}


VersionedSerializer *RasterSerializerV1::create(QDataStream &stream, const QString &version)
{
    return new RasterSerializerV1(stream, version);
}
