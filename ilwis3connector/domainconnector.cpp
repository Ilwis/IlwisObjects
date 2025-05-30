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

#include <QColor>
#include "kernel.h"
#include "geos/geom/Coordinate.h"
#include "connectorinterface.h"
#include "mastercatalog.h"
#include "ilwisobjectconnector.h"
#include "catalogexplorer.h"
#include "catalogconnector.h"
#include "inifile.h"
#include "ilwisdata.h"
#include "domainitem.h"
#include "domain.h"
#include "itemdomain.h"
#include "identifieritem.h"
#include "thematicitem.h"
#include "numericrange.h"
#include "datadefinition.h"
#include "colorrange.h"
#include "colordomain.h"
#include "itemrange.h"
#include "identifierrange.h"
#include "rawconverter.h"
#include "ilwis3connector.h"
#include "numericdomain.h"
#include "binaryilwis3table.h"
#include "ilwis3range.h"
#include "textdomain.h"
#include "interval.h"
#include "intervalrange.h"
#include "domainconnector.h"
//#include "coordinatesystem.h"
//#include "coverage.h"
#include "coordinatedomain.h"

#include <QUrlQuery>

using namespace Ilwis;
using namespace Ilwis3;

ConnectorInterface *DomainConnector::create(const Resource& resource, bool load, const IOOptions &options) {
    return new DomainConnector(resource, load, options);

}

DomainConnector::DomainConnector(const Resource& resource, bool load, const IOOptions &options) : Ilwis3Connector(resource, load, options)
{
}

bool DomainConnector::loadMetaData(IlwisObject* data, const IOOptions& options)
{
    Ilwis3Connector::loadMetaData(data, options);

    if (type() == itUNKNOWN) {
        kernel()->issues()->log(TR(ERR_INVALID_PROPERTY_FOR_2).arg("Domain type",_resource.name()));
        return false;
    }
    if ( type() == itNUMERICDOMAIN) {
        return handleValueDomains(data);
    } else if ( type() == itITEMDOMAIN || type() == itDOMAIN) {
        return handleItemDomains(data);
    }

    return false;
}

bool DomainConnector::handleIdDomain(IlwisObject* data) {
    ItemDomain<IndexedIdentifier> *iddomain = static_cast<ItemDomain<IndexedIdentifier> *>(data);
    bool ok;
    quint32 nritems = _odf->value("DomainIdentifier","Nr").toLong(&ok);
    if ( !ok) {
        return ERROR2(ERR_INVALID_PROPERTY_FOR_2,"domain items", data->name());
    }
    QString prefix = _odf->value("DomainSort","Prefix");
    iddomain->addItem(new IndexedIdentifier(prefix,0, nritems));

    return true;
}

bool DomainConnector::handleItemDomains(IlwisObject* data) {
    QString domtype = _odf->value("Domain","Type");
    if (domtype == "DomainCoord") {
        return true;
    }
    bool hasDataFile = _odf->value("TableStore","Col1") == "Ord";
    if ( (domtype == "DomainUniqueID" || domtype ==  "DomainIdentifier") && !hasDataFile)
    { // no table found? internal domain
        return handleIdDomain(data);
    }
    Ilwis3::BinaryIlwis3Table tbl ;
    tbl.load(_odf);
    quint32 indexName = tbl.index("Name");
    if (indexName == iUNDEF) { // no name column in the table ?
        kernel()->issues()->log(TR(ERR_COLUMN_MISSING_2).arg("Name",_odf->url()));
        return false;
    }
    quint32 indexCode = tbl.index("Code"); // not mandatory
    if ( domtype == "DomainGroup"){
        ItemDomain<Interval> *intervaldomain = static_cast<ItemDomain<Interval> *>(data);
        QString itemName, itemCode = sUNDEF;
        double bound = 0, oldBound = -1e307; // -1e307 to represent minus-infinity; with rUNDEF/-1e308 the validity-check of the NumericRange fails.
        quint32 indexBounds = tbl.index("Bounds");
        std::vector<std::pair<double, std::pair<QString, QString>>> domainItems; // first: load the data from the domain-table
        for (quint32 i = 0; i < tbl.rows(); ++i) {
            tbl.get(i, indexName, itemName);
            if (indexCode != iUNDEF)
                tbl.get(i, indexCode, itemCode);
            tbl.get(i, indexBounds, bound);
            std::pair<double, std::pair<QString, QString>> domainItem = std::pair<double, std::pair<QString, QString>>(bound, std::pair<QString, QString>(itemName, itemCode));
            domainItems.push_back(domainItem);
        }
        std::sort(domainItems.begin(), domainItems.end()); // second: sort the table data on "Bounds", and generate ilwis4-style min/max bounds
        std::map<QString, Interval*> itemIntervalMap;
        for (std::vector<std::pair<double, std::pair<QString, QString>>>::iterator it = domainItems.begin(); it != domainItems.end(); ++it) {
            bound = it->first;
            itemName = it->second.first;
            itemCode = it->second.second;
            Interval* item = new Interval({ itemName,NumericRange(oldBound, bound) });
            itemIntervalMap[itemName] = item;
            oldBound = bound;
        }
        for (quint32 i = 0; i < tbl.rows(); ++i) { // third: add the Interval items to the intervaldomain in the original order
            tbl.get(i, indexName, itemName);
            Interval* item = itemIntervalMap[itemName];
            intervaldomain->addItem(item);
        }
    }else if  (domtype == "DomainIdentifier"){
        ItemDomain<NamedIdentifier> *tdomain = static_cast<ItemDomain<NamedIdentifier> *>(data);
        QString itemName;
        for(quint32 i = 0; i < tbl.rows(); ++i) {
            tbl.get(i,indexName,itemName);
            NamedIdentifier *item = new NamedIdentifier(itemName);
            tdomain->addItem(item);
        }
    }else {
        ItemDomain<ThematicItem> *tdomain = static_cast<ItemDomain<ThematicItem> *>(data);
        QString itemName, itemCode = sUNDEF;
        for(quint32 i = 0; i < tbl.rows(); ++i) {
            tbl.get(i,indexName,itemName);
            if ( indexCode != iUNDEF)
                tbl.get(i,indexCode,itemCode);
            ThematicItem *item = new ThematicItem({itemName,itemCode});
            tdomain->addItem(item);
        }
    }

    return true;
}


bool DomainConnector::handleValueDomains(IlwisObject* data) {
    NumericDomain *vdata = static_cast<NumericDomain*>(data);
    IlwisTypes odfType = _resource.ilwisType();
    Range *range = 0;
    if ( (odfType & itDOMAIN) != 0){ // domain objects
        range = handleValueDomainObjects();
    } else {
        QString section;
        if ( odfType <= itRASTER) { //  the 'basemaps' objects
            section = "BaseMap"    ;
        } else if (odfType == itTABLE) {
            QUrlQuery queryItem(_resource.url());
            QString colName = queryItem.queryItemValue("column");
            section = QString("Col:%1").arg(colName);
        }
        if ( section != "") {
            range = Ilwis3Range::findRange(_odf, section);

        } else {
            kernel()->issues()->log(TR("Illegal type %1 for domain").arg(type()));
            return false;
        }
    }
    // TODO: columns domain handling, no seperate odf, so special type(??)

    if (!range) {
        kernel()->issues()->log(TR(ERR_NO_INITIALIZED_1).arg(data->name()));
        return false;
    }
    vdata->range(range);
    IDomain dom("value");
    vdata->setParent(dom);


    return true;
}

Range * DomainConnector::handleValueDomainObjects() const {
    bool isOk = true;
    double vmin, vmax;
    QString type = _odf->value("DomainValue", "Type");
    vmin = _odf->value(type, "Min").toDouble(&isOk);
    if (!isOk) {
        kernel()->issues()->log(TR(ERR_INVALID_PROPERTY_FOR_2).arg("Minimum", _resource.name()));
        return 0;
    }

    vmax = _odf->value(type ,"Max").toDouble(&isOk);
    if (!isOk) {
        kernel()->issues()->log(TR(ERR_INVALID_PROPERTY_FOR_2).arg("Maximum", _resource.name()));
        return 0;
    }
    double resolution = _odf->value(type , "Step").toDouble(&isOk);
    if (!isOk || fmod(resolution , 1.0) == 0)
        return new NumericRange(vmin, vmax,1);
    return new NumericRange(vmin, vmax, resolution);

}


QString DomainConnector::parseDomainInfo(const QString& inf) const{
    if ( inf != sUNDEF) {
        QStringList parts = inf.split(";");
        if ( parts.size() > 3) {
            QString dmtype = parts[2];
            return dmtype.toLower();
        }
    }
    return sUNDEF;
}

bool DomainConnector::storeMetaDataSortDomain(Domain *dom, IlwisTypes valueType) {

    ItemDomain<NamedIdentifier> *piddomain = static_cast< ItemDomain<NamedIdentifier> *>(dom);
    INamedIdDomain iddomain;
    iddomain.set(piddomain);

    auto writeColumnFunc = [&] (const QString& name, const QString& domName, const QString& domInfo, const QString& rng, const QString& storeType) -> void {
        auto tm = IniFile::FormatElement((quint32)Time::now().toTime_t());
        _odf->setKeyValue(name, "Time", tm);
        _odf->setKeyValue(name, "Version", "3.1");
        _odf->setKeyValue(name, "Class", "Column");
        _odf->setKeyValue(name, "Domain", domName);
        _odf->setKeyValue(name, "DomainInfo", domInfo);
        if ( rng != sUNDEF)
            _odf->setKeyValue(name, "Range", rng);
        _odf->setKeyValue(name, "ReadOnly", "No");
        _odf->setKeyValue(name, "OwnedByTable", "No");
        _odf->setKeyValue(name, "Type", "ColumnStore");
        _odf->setKeyValue(name, "StoreType", storeType);

    };

    auto tm = IniFile::FormatElement((quint32)Time::now().toTime_t());
    _odf->setKeyValue("Table", "Time", tm);
    _odf->setKeyValue("Table","Version","3.1");
    _odf->setKeyValue("Table","Class","Table");
    _odf->setKeyValue("Table","Domain","String.dom");
    _odf->setKeyValue("Table","Type","TableStore");
    _odf->setKeyValue("Table","DomainInfo", "String.dom;String;string;0;;");
    int nrOfColumns = 3;
    if ( valueType == itTHEMATICITEM)
        nrOfColumns = 5;
    if ( valueType == itNUMERICITEM)
        nrOfColumns = 6;
    _odf->setKeyValue("Table","Columns",IniFile::FormatElement(nrOfColumns));
    _odf->setKeyValue("Table","Records", IniFile::FormatElement(iddomain->count()));
    _odf->setKeyValue("DomainSort","Prefix", "");
    _odf->setKeyValue("DomainSort","Sorting","Alphabetical");
    QString typeName = "DomainIdentifier";
    if ( valueType == itTHEMATICITEM)
        typeName = "DomainClass";
    if ( valueType == itNUMERICITEM)
        typeName = "DomainGroup";

    QFileInfo inf(QUrl(_odf->url()).toLocalFile());
    _odf->setKeyValue("Domain", "Type", typeName);
    if (hasType(valueType, itTHEMATICITEM | itNUMERICITEM))
        _odf->setKeyValue("Domain", "Representation", inf.baseName() + ".rpr"); // a gamble that the rpr will have the same name as the domain; at this stage we don't have access to tbe representation (we would have to find the datadefinition object)
    _odf->setKeyValue( hasType(valueType,itTHEMATICITEM | itNUMERICITEM) ? "DomainClass" : "DomainIdentifier", "Nr", IniFile::FormatElement(iddomain->count()));

    QString dataName  = inf.baseName() + ".dm#";
    _odf->setKeyValue("TableStore", "Data", dataName);
    _odf->setKeyValue("TableStore", "Col0", "Name");
    _odf->setKeyValue("TableStore", "Col1", "Ord");
    _odf->setKeyValue("TableStore", "Col2", "Ind");
    _odf->setKeyValue("TableStore", "Type", "TableBinary");
    if (hasType(valueType,itTHEMATICITEM | itNUMERICITEM)) {
        _odf->setKeyValue("TableStore", "Col3", "Code");
        _odf->setKeyValue("TableStore", "Col4", "Description");
        if ( valueType == itNUMERICITEM)
             _odf->setKeyValue("TableStore", "Col5", "Bounds");
    }

    writeColumnFunc("Col:Name","String.dom","String.dom;String;string;0;;", sUNDEF, "String");
    writeColumnFunc("Col:Ord","value.dom","value.dom;Long;value;0;-9999999.9:9999999.9:0.1:offset=0;", "-32766:32767:offset=0","Int");
    writeColumnFunc("Col:Ind","value.dom","value.dom;Long;value;0;-9999999.9:9999999.9:0.1:offset=0;", "-32766:32767:offset=0","Int");
    if ( hasType(valueType,itTHEMATICITEM | itNUMERICITEM)) {
        writeColumnFunc("Col:Code","String.dom","String.dom;String;string;0;;", sUNDEF, "String");
        writeColumnFunc("Col:Description","String.dom","String.dom;String;string;0;;", sUNDEF, "String");
        if ( valueType == itNUMERICITEM)
             writeColumnFunc("Col:Bounds","value.dom","value.dom;Long;value;0;-9999999.9:9999999.9:0.1:offset=0", "-1e+100:1e+100:0.000000:offset=0", "Real");
    }


    BinaryIlwis3Table ilw3tbl;
    std::ofstream output_file;
    if(!ilw3tbl.openOutput(inf.absolutePath() + "/"+ dataName, output_file))
        return false;
    ITextDomain txtdom;
    txtdom.prepare();
    DataDefinition deftxt(txtdom);
    INumericDomain numdom;
    numdom.prepare("integer");
    DataDefinition deford(numdom);
    deford.range(new NumericRange(-32766,32767,1));
    DataDefinition defind(numdom);
    defind.range(new NumericRange(-32766,32767,1));

    ilw3tbl.addStoreDefinition(deftxt);
    ilw3tbl.addStoreDefinition(deford);
    ilw3tbl.addStoreDefinition(defind);
    if ( valueType == itNUMERICITEM){
        ilw3tbl.addStoreDefinition(deftxt);
        ilw3tbl.addStoreDefinition(deftxt);
        ilw3tbl.addStoreDefinition(DataDefinition(IDomain("value")));
    }
    if ( valueType == itTHEMATICITEM) {
        ilw3tbl.addStoreDefinition(deftxt);
        ilw3tbl.addStoreDefinition(deftxt);
    }

    std::map<quint32, std::vector<QVariant>> orderedRecords;
    for(DomainItem *item : iddomain){
        int columnCount=  3;
        if ( valueType == itTHEMATICITEM)
            columnCount = 5;
        if ( valueType == itNUMERICITEM)
            columnCount = 6;
        std::vector<QVariant> record(columnCount);
        record[0] = item->name();
        record[1] = item->raw() + 1;
        record[2] = item->raw() + 1;
        if ( valueType == itNUMERICITEM){
            Interval *intervalitem = static_cast<Interval *>(item);
            record[3] = intervalitem->code();
            record[4] = intervalitem->description();
            record[5] = intervalitem->range().max();
        }
        else if ( valueType == itTHEMATICITEM ) {
            ThematicItem *thematicItem = static_cast<ThematicItem *>(item);
            record[3] = thematicItem->code();
            record[4] = thematicItem->description();
        }
        orderedRecords[item->raw()] = record;
    }
    for(const auto& rec : orderedRecords){

        ilw3tbl.storeRecord(output_file, rec.second);
    }

    output_file.close();

    return true;

}

bool DomainConnector::storeMetaData(IlwisObject *data, const IOOptions &options)
{
    Domain *dom = static_cast<Domain *>(data);
    QString dmName = dom->name();
    QString alias = kernel()->database()->findAlias(dmName,"domain","ilwis3");
    if ( alias != sUNDEF)
        return true; // nothing to be done, already exists as a system domain
    if(!Ilwis3Connector::storeMetaData(data, itDOMAIN))
        return false;

    _odf->setKeyValue("Ilwis", "Type", "Domain");
    if ( dom->ilwisType() == itNUMERICDOMAIN) {
        SPNumericRange numRange = dom->range<NumericRange>();
        int width=12;
        QString type = "DomainValueInt";
        if ( dom->valueType() & (itINT8 | itUINT8)){
            width=3;
        } else if ( dom->valueType() & (itINT16 | itUINT16) ){
            width = 8;
        } else if ( dom->valueType() & (itINT32 | itUINT32) ){
            width=10;
        } else
            type = "DomainValueReal";

        _odf->setKeyValue("Domain", "Type", "DomainValue");
        _odf->setKeyValue("Domain", "Width", IniFile::FormatElement(width));
        _odf->setKeyValue("DomainValue", "Type", type);
        _odf->setKeyValue(type, "Min", (type == "DomainValueReal") ? IniFile::FormatElement(numRange->min()) : IniFile::FormatElement((long)(numRange->min())));
        _odf->setKeyValue(type, "Max", (type == "DomainValueReal") ? IniFile::FormatElement(numRange->max()) : IniFile::FormatElement((long)(numRange->max())));
        _odf->setKeyValue("Domain","Representation","gray.rpr");
        if ( numRange->resolution() != 1) {
            _odf->setKeyValue(type, "Step", IniFile::FormatElement(numRange->resolution()));
        }
    } else if ( dom->valueType() == itTHEMATICITEM) {
        storeMetaDataSortDomain(dom, itTHEMATICITEM);
    } else if ( dom->valueType() == itNUMERICITEM) {
        storeMetaDataSortDomain(dom, itNUMERICITEM);
    } else if ( dom->valueType() & itIDENTIFIERITEM) {
        storeMetaDataSortDomain(dom, itIDENTIFIERITEM);
    } else if ( dom->ilwisType() == itTEXTDOMAIN) {
    } else if ( dom->ilwisType() == itCOLORDOMAIN) {
    } else if ( dom->ilwisType() == itCOORDDOMAIN) {
    }

    _odf->store("dom", sourceRef().toLocalFile());
    return true;
}


IlwisObject *DomainConnector::fromValueRange() const{
    QString range = _odf->value("BaseMap","Range");
    if ( range != sUNDEF) {
        QStringList parts = range.split(":");
        if ( parts.size() > 2){
            double vmin = parts[0].toDouble();
            double vmax = parts[1].toDouble();
            double resolution = 1;
            if ( parts.size() > 3 ){
                resolution = parts[2].toDouble();
            }
            return new NumericDomain( new NumericRange(vmin, vmax, resolution));
        }
    }
    return 0;
}
IlwisObject *DomainConnector::create() const
{
    //TODO: other domain types time, coordinatesystem
    QString subtype = sUNDEF;
    if ( type() & itCOVERAGE) {
        subtype = parseDomainInfo( _odf->value("BaseMap","DomainInfo"));
        if ( subtype == "image.dom" || subtype == "value.dom"){
            return fromValueRange();
        }

    } else if( type() & itTABLE) {
        QUrlQuery queryItem(_resource.url());
        QString colName = queryItem.queryItemValue("column");
        if ( colName != sUNDEF) {
            subtype = parseDomainInfo( _odf->value(QString("Col:%1").arg(colName),"DomainInfo"));
        }

    }

    if ( type() == itNUMERICDOMAIN)
        return new NumericDomain(_resource);
    else if (type() == itITEMDOMAIN || type() == itDOMAIN) { // second case is for internal domains
        subtype =_odf->value("Domain", "Type");
        if (subtype == "DomainCoord") {
            CoordinateDomain * crddom = new CoordinateDomain(_resource);
            QString csyName = _resource.url().toString();
            if (csyName.right(4) == ".csy") {
                ICoordinateSystem csy;
                csy.prepare(csyName);
                crddom->setCoordinateSystem(csy);
            }
            return crddom;
        }
        bool internal = _odf->value("TableStore","Col1") != "Ord"; // if not exists it is no table so no datafile
        if ( subtype == "DomainUniqueID" || internal)
            return new ItemDomain<IndexedIdentifier>(_resource);
        if ( subtype == "DomainIdentifier")
            return new ItemDomain<NamedIdentifier>(_resource);
        if ( subtype == "DomainClass" || subtype == "DomainSort")
            return new ItemDomain<ThematicItem>(_resource);
        if ( subtype == "DomainGroup")
            return new ItemDomain<Interval>(_resource);
        subtype = parseDomainInfo( _odf->value("BaseMap","DomainInfo"));
        if ( subtype.left(5) == "image" || subtype.left(5) == "value"){
            return fromValueRange();
        }
        if ( subtype.left(5) == "color" ){
            ColorDomain *dm = new ColorDomain(_resource);
            ContinuousColorRange *rng = new ContinuousColorRange(QColor("#000000"), QColor("(#ffffff"));
            dm->range(rng);
            return dm;
        }
    }
    return 0;
}

QString DomainConnector::format() const
{
    return "domain";
}
