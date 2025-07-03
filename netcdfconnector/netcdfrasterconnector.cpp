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

#include "kernel.h"
#include "version.h"
#include "ilwisdata.h"
#include "connectorinterface.h"
#include "domain.h"
#include "datadefinition.h"
#include "columndefinition.h"
#include "factory.h"
#include "abstractfactory.h"
#include "ilwisobjectfactory.h"
#include "connectorfactory.h"
#include "table.h"
#include "raster.h"
#include "ilwisobjectconnector.h"
#include "netcdfrasterconnector.h"
#include "catalogexplorer.h"
#include "netcdfcatalogconnector.h"
#include "ilwisobjectconnector.h"
#include "ncFile.h"


using namespace Ilwis;
using namespace NetCdf;

ConnectorInterface *NetCdfRasterConnector::create(const Ilwis::Resource &resource, bool load, const IOOptions& options)
{
    return new NetCdfRasterConnector(resource, load, options);
}

NetCdfRasterConnector::NetCdfRasterConnector(const Ilwis::Resource &resource, bool load, const IOOptions& options) : IlwisObjectConnector(resource, load, options)
{
    _version = 1;
}

IlwisObject *NetCdfRasterConnector::create() const {
    return new RasterCoverage(_resource);
}

bool NetCdfRasterConnector::loadMetaData(IlwisObject *obj, const IOOptions &)
{
    auto url = obj->resource().url();
    auto idx = url.toString().lastIndexOf("/");
    QString urls = url.toString();
    QString ncdfFile = urls.left(idx);
    QString product = url.toString().right(url.toString().size() - idx - 1);
    auto items = NetCdfCatalogExplorer::createResources(ncdfFile, product);
    if ( items.size() == 0)
        return false;

    const Resource& res = items[0];
    QString csyCode = "code=proj4:" + res["coordinatesystem"].toString();
    ICoordinateSystem csy;
    if (!csy.prepare(csyCode))
        return false;
    RasterCoverage *raster = static_cast<RasterCoverage *>(obj);
    raster->coordinateSystem(csy);


    Size<> sz(res.dimensions());

    IGeoReference grf;
    QString grfP = res["georeference"].toString();
    if (!grf.prepare(grfP))
        return false;
    raster->georeference(grf);
    raster->size(sz);
    std::vector<double> bands(sz.zsize());
    for(unsigned int z=0 ; z < bands.size(); ++z)
        bands[z] = z;
    QString dom = res["domain"].toString();
    raster->setDataDefintions({ "code=domain:value" }, bands);

    auto keys = res["dimension vars"].toString();
    _dimVars = keys.split("|");
    auto &resR = raster->resourceRef();
    if ( res.hasProperty("undefined_value"))
        resR.addProperty("undefined_value", res["undefined_value"]);
    if ( res.hasProperty("used var"))
        resR.addProperty("used var", res["used var"]);
    auto props = res.properties();
    for (auto iter = props.begin(); iter != props.end(); ++iter){
        auto key = iter.key();
        if ( key.indexOf("metadata") != -1)
            resR.addProperty(key, res[key]);
    }


    return true;
}
QString getAttrValue(const netCDF::NcVar& v, const QString& label) {
    try{
        auto att2 = v.getAtt(label.toStdString());
        std::string s;
        att2.getValues(s);
        return QString::fromStdString(s)    ;

    } catch(const std::exception& ex){
        return sUNDEF;
    }
    return sUNDEF;



}

template<class T> void NetCdfRasterConnector::getValues(const std::vector<size_t>& index, const std::vector<size_t>& count, int maxX, int maxY,  bool xy, netCDF::NcVar& var,std::vector<T>&data, double undef, std::vector<double>& values){
    try
    {

    var.getVar(index,count, data.data());


    for(int y =0; y < maxY; ++y){
        for(int x = 0; x < maxX; ++x){
            quint64 localLinLoc = y * maxX + x;
            quint64 linLoc = xy ? localLinLoc: x * maxY + y;
            double v = data[linLoc];
            values[localLinLoc] = v == undef ? rUNDEF : v;
        }
    }
    }catch(netCDF::exceptions::NcException& e)
    {
      e.what();
      std::cout<<"FAILURE*************************************"<<endl;
    }
}
void NetCdfRasterConnector::getValue(const std::vector<size_t>& index, const std::vector<size_t>& count, int maxX, int maxY,  bool xy, netCDF::NcVar& var, std::vector<double>& values, double undef){


    netCDF::NcType tp = var.getType();
    auto dataType = tp.getTypeClass();

    if ( dataType == netCDF::NcType::nc_FLOAT){
        std::vector<float> data(maxX * maxY, flUNDEF);
        getValues(index, count, maxX, maxY, xy, var, data, undef, values);
    }
    if ( dataType == netCDF::NcType::nc_DOUBLE){
        std::vector<double> data(maxX * maxY, flUNDEF);
        getValues(index, count, maxX, maxY, xy, var, data, undef, values);
    }
    if ( dataType == netCDF::NcType::nc_INT){
        std::vector<quint32> data(maxX * maxY, flUNDEF);
        getValues(index, count, maxX, maxY, xy, var, data, undef, values);
    }
    if ( dataType == netCDF::NcType::nc_UINT){
        std::vector<quint32> data(maxX * maxY, flUNDEF);
        getValues(index, count, maxX, maxY, xy, var, data, undef, values);
    }
    if ( dataType == netCDF::NcType::nc_INT64){
        std::vector<qint64> data(maxX * maxY, flUNDEF);
        getValues(index, count, maxX, maxY, xy, var, data, undef, values);
    }
    if ( dataType == netCDF::NcType::nc_UINT64){
        std::vector<quint64> data(maxX * maxY, flUNDEF);
        getValues(index, count, maxX, maxY, xy, var, data, undef, values);
    }
    if ( dataType == netCDF::NcType::nc_SHORT){
        std::vector<qint16> data(maxX * maxY, flUNDEF);
        getValues(index, count, maxX, maxY, xy, var, data, undef, values);
    }
    if ( dataType == netCDF::NcType::nc_USHORT){
        std::vector<quint16> data(maxX * maxY, flUNDEF);
        getValues(index, count, maxX, maxY, xy, var, data, undef, values);
    }
    if ( dataType == netCDF::NcType::nc_BYTE){
        std::vector<qint8> data(maxX * maxY, flUNDEF);
        getValues(index, count, maxX, maxY, xy, var, data, undef, values);
    }
    if ( dataType == netCDF::NcType::nc_UBYTE){
        std::vector<quint8> data(maxX * maxY, flUNDEF);
        getValues(index, count, maxX, maxY, xy, var, data, undef, values);
    }
}

bool NetCdfRasterConnector::loadData(IlwisObject* obj, const IOOptions& options) {

    const Resource& res = obj->resourceRef();
    auto urls = res.url().toString();
    auto idx = urls.lastIndexOf("/");
    QString ncdfFile = urls.left(idx);
    QString product = urls.right(urls.size() - idx - 1);
    auto filename = QUrl(ncdfFile).toLocalFile();

    netCDF::NcFile file(filename.toStdString(), netCDF::NcFile::read);
    std::multimap<std::string,netCDF::NcVar> vars = file.getVars();
    std::multimap<std::string,netCDF::NcVar>::const_iterator iter = vars.find(product.toStdString());

    if ( iter == vars.end())
        return false;
    netCDF::NcVar var = iter->second;
    RasterCoverage *raster = static_cast<RasterCoverage *>(obj);
    UPGrid& grid = raster->gridRef();

    quint32 y = options.contains("y") ? options["y"].toUInt(): 0;
    quint32 z = options.contains("z") ? options["z"].toUInt(): 0;
    Size<> sz = raster->size();

    int n = var.getDimCount();
    std::vector< netCDF::NcDim > dims = var.getDims();
    std::vector<int> axisLocations(n);
    int acount = 0;
    for(netCDF::NcDim dim : dims){
        auto axisName = dim.getName();
        std::multimap< std::string, netCDF::NcVar >::const_iterator iterAxis= vars.find(axisName);
        auto axisType = NetCdfCatalogExplorer::getAxisType(iterAxis);
        if ( axisType == NetCdfCatalogExplorer::atX )
           axisLocations[0] = acount;
        if ( axisType == NetCdfCatalogExplorer::atY)
           axisLocations[1] = acount;
        if ( axisType == NetCdfCatalogExplorer::atZ)
           axisLocations[2] = acount;
        ++acount;
    }
    std::vector<size_t> index = {0,0,0,0};
    index[axisLocations[0]] = 0;
    index[axisLocations[1]] = y;
    index[axisLocations[2]] =  z;
    std::vector<size_t> count = {1,0,0,0};
    int maxX = sz.xsize();
    int maxY = grid->blockSize(options["blockindex"].toInt()) / maxX;
    count[axisLocations[2]] = 1;
    count[axisLocations[0]] = maxX;
    count[axisLocations[1]] = maxY;
    std::vector<double> values(maxX * maxY, rUNDEF);
    bool xy = axisLocations[0] < axisLocations[1];
    QString s;
    if ( _resource.hasProperty("undefined_value")){
       s = _resource["undefined_value"].toString();
    }else
        s = getAttrValue(var, "_FillValue");
    double undef = s == sUNDEF ? rUNDEF : s.toDouble();
    getValue(index,count,maxX,maxY,xy,var,values, undef);
    grid->setBlockData(z, y, values);
     _binaryIsLoaded = true;


    return true;
}

QString NetCdfRasterConnector::provider() const
{
    return "Ilwis4";
}
