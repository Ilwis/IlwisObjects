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

#include <QCoreApplication>
#include <QDir>
#include <QDesktopServices>
#include <QJsonDocument>
#include <QLibrary>
#include "kernel.h"
#include "ilwisdata.h"
#include "errorobject.h"
#include "abstractfactory.h"
#include "connectorinterface.h"
#include "ilwisobjectconnector.h"
#include "catalogexplorer.h"
#include "catalogconnector.h"
#include "catalog.h"
#include "ilwiscontext.h"
#include "mastercatalog.h"
#include "mastercatalogcache.h"

Ilwis::IlwisContext *Ilwis::IlwisContext::_context = 0;

using namespace Ilwis;

IlwisContext* Ilwis::context(const QString & ilwisDir, int runMode) {
    if (Ilwis::IlwisContext::_context == 0) {
        Ilwis::IlwisContext::_context = new Ilwis::IlwisContext(runMode);
        Ilwis::IlwisContext::_context->init(ilwisDir);

    }
    return Ilwis::IlwisContext::_context;
}



IlwisContext::IlwisContext(int runMode) :  _memoryLimit(9e8), _memoryLeft(_memoryLimit), _runMode(runMode)
{
    // _workingCatalog = new Catalog(); // empty catalog>

}

IlwisContext::~IlwisContext()
{
   ICatalog wcatalog = workingCatalog();
    if ( wcatalog.isValid()){
        _configuration.putValue("users/" + currentUser() + "/workingcatalog",wcatalog->resource().url().toString());
        _configuration.store();
    }
  
}

void IlwisContext::addSystemLocation(const QUrl &resource)
{
    if ( std::find(_systemlocations.begin(), _systemlocations.end(), resource) == _systemlocations.end() ) {
        _systemlocations.push_back(resource);
    }

}

void IlwisContext::removeSystemLocation(const QUrl &)
{
    //TODO:
}

QFileInfo IlwisContext::ilwisFolder() const {
    return this->_ilwisDir;
}

QString IlwisContext::setCacheLocation(const QString loc){
    QString locat = loc;
    if ( loc == "" || loc == sUNDEF){
        QString location = ilwisconfig("users/" + Ilwis::context()->currentUser() + "/cache-location",QString(sUNDEF));
        if ( location == sUNDEF){
            QDir localDir(QStandardPaths::writableLocation(QStandardPaths::CacheLocation));
            locat = localDir.absolutePath();
        }else
            locat = location;
    }
    _cacheLocation = QUrl::fromLocalFile(locat);
    QDir cacheDir(locat);

    if (!cacheDir.exists()){
        cacheDir.mkpath(locat);
    }
    QStringList files = cacheDir.entryList(QStringList() << "gridblock*.*" << "osm*.png", QDir::Files);
    for(QString file : files)
        cacheDir.remove(file);

    return locat;
}
QString Ilwis::IlwisContext::setInternalCatalog(const QString& loc)
{
    QString datalocation = loc;
    if ( loc == "" || loc == sUNDEF){

        datalocation = ilwisconfig("users/" + Ilwis::context()->currentUser() + "/internalcatalog-location",QString(sUNDEF));
        if ( datalocation == sUNDEF){
            datalocation = QStandardPaths::writableLocation(QStandardPaths::DataLocation) + "/internalcatalog";
        }else
            datalocation = QUrl(datalocation).toLocalFile();
        datalocation = OSHelper::neutralizeFileName(datalocation);
    }

    QDir localDir(datalocation);

    if (!localDir.exists()){
        localDir.mkpath(datalocation);
    }
    _persistentInternalCatalog = QUrl::fromLocalFile(datalocation);
    auto files = localDir.entryList(QStringList() << "*", QDir::Files);
    for(QString file : files)
        localDir.remove(file);
    return datalocation;
}

void IlwisContext::init(const QString &ilwisDir)
{
    if (ilwisDir.length() > 0) {
        this->_ilwisDir = QFileInfo(ilwisDir);
        if (!this->_ilwisDir.isDir()) {
            printf("User-supplied Ilwis directory '%s' not found\n",this->_ilwisDir.filePath().toStdString().c_str());
            this->_ilwisDir = QFileInfo(qApp->applicationDirPath());
        } else
            qApp->addLibraryPath(this->_ilwisDir.absolutePath() + "/plugins"); // also inform Qt where its "plugins" folder is installed, so that it doesn't use the hardcoded qt_plugpath inside Qt5Core.dll
    } else
        this->_ilwisDir = QFileInfo(qApp->applicationDirPath());

    QString loc = QStandardPaths::writableLocation(QStandardPaths::ConfigLocation);

    QString configfile = loc + "/" + "ilwis.config";
    QFileInfo file;
    file.setFile(configfile);
    if ( !file.exists()){
       QFileInfo resourceInf(resourcesLocation());
        configfile = resourcesLocation() + "/ilwis.config";
        file.setFile(configfile);
    }
    OSHelper::loadExtraLibs(_ilwisDir.absoluteFilePath());
 
    _configuration.prepare(file.absoluteFilePath());

    setCacheLocation();
    setInternalCatalog();

    mastercatalog()->addContainer(INTERNAL_CATALOG_URL);
    Resource res = mastercatalog()->name2Resource(INTERNAL_CATALOG_URL.toString(),itCATALOG);
    res.name("temporary catalog",false,true);
    mastercatalog()->addContainer(persistentInternalCatalog());
    mastercatalog()->addContainer(QUrl("ilwis://operations"));

    _systemCatalog.prepare("ilwis://system");
    mastercatalog()->addContainer(QUrl("ilwis://system/domains"));
    mastercatalog()->addContainer(QUrl("ilwis://system/coordinatesystems"));
    mastercatalog()->addContainer(QUrl("ilwis://system/representations"));
	mastercatalog()->addContainer(QUrl("ilwis://system/representations/item"));
	mastercatalog()->addContainer(QUrl("ilwis://system/representations/value"));
    mastercatalog()->addContainer(QUrl("ilwis://system/ellipsoids"));
    mastercatalog()->addContainer(QUrl("ilwis://system/projections"));
    mastercatalog()->addContainer(QUrl("ilwis://system/datums"));
    mastercatalog()->addContainer(QUrl("ilwis://system/coverages"));
    mastercatalog()->addContainer(QUrl("ilwis://system/scripts"));
	mastercatalog()->addContainer(QUrl("ilwis://system/tables"));

    if (!hasType(_runMode, rmDESKTOP)){
        initializationFinished(true);
    }

}

ICatalog IlwisContext::workingCatalog() const{
//    if ( _workingCatalog.hasLocalData())
//        return static_cast<Catalog *>(_workingCatalog.localData());
    Locker<std::mutex> lock(_lock);
    const QVariant *var = kernel()->getFromTLS("workingcatalog");
    if ( var && var->isValid()){
        ICatalog cat = var->value<ICatalog>();
        return cat;
    }
    return ICatalog();
}

const ICatalog &IlwisContext::systemCatalog() const
{
    return _systemCatalog;
}

void IlwisContext::setWorkingCatalog(const ICatalog &cat)
{
    // the ilwis default workspace is just is a placeholder for everything goes; so we don't assign it
    if ( !cat.isValid() || cat->resource().url().toString() == Catalog::DEFAULT_WORKSPACE)
        return;

    Locker<std::mutex> lock(_lock);
    QVariant *var = new QVariant();
    var->setValue(cat);
    kernel()->setTLS("workingcatalog", var);
    context()->configurationRef().putValue("users/" + currentUser() + "/workingcatalog",cat->resource().url().toString());
    //QFileInfo inf(cat->resource().url().toLocalFile());
}

QUrl IlwisContext::cacheLocation() const
{
    return _cacheLocation;
}

QUrl IlwisContext::persistentInternalCatalog() const
{
    return _persistentInternalCatalog;
}

quint64 IlwisContext::memoryLeft() const
{
    return _memoryLeft;
}
quint64 IlwisContext::changeMemoryLeft(qint64 amount)
{
    if ( (_memoryLeft + amount) > 0) {
        _memoryLeft += amount;
    }
    else
        _memoryLeft = 0;

    return _memoryLeft;
}

IlwisConfiguration &IlwisContext::configurationRef()
{
    return _configuration;
}

const IlwisConfiguration &IlwisContext::configuration() const
{
    return _configuration;
}

QFileInfo IlwisContext::resourceRoot() const
{
    QString root = ilwisconfig("system-settings/resource-root",QString("app-base"));
    if ( root == "app-base"){
        return QFileInfo(ilwisFolder().absoluteFilePath() + "/resources");
    }
    QFileInfo inf(root + "/resources");
    if ( inf.exists())
        return inf;

    return QFileInfo(ilwisFolder().absoluteFilePath() + "/resources");
}

QString IlwisContext::ipv4() const
{
    if ( _ipv4 == sUNDEF){
        const_cast<IlwisContext *>(this)->_ipv4 = const_cast<IlwisContext *>(this)->configurationRef()("server-settings/ipv4-address", QString(""));
    }
    return _ipv4;
}

QString IlwisContext::currentUser() const
{
    return "user-0"; // by default atm,.
}

int IlwisContext::runMode() const
{
    return _runMode;
}

void IlwisContext::runMode(int mode)
{
    _runMode = mode;
}

bool IlwisContext::initializationFinished() const
{
   // Locker<stdV::mutex> locker(_lock);
    return _initializationFinished;
}

void IlwisContext::initializationFinished(bool yesno)
{
   //Locker<std::mutex> locker(_lock);
    _initializationFinished = yesno;
}

QString IlwisContext::resourcesLocation(const QString &internalName) const{
    QString loc;
    if ( internalName == "")
        loc = _ilwisDir.absoluteFilePath() + "/resources";
    else
        loc = _ilwisDir.absoluteFilePath() + "/extensions/" + internalName + "/resources";
    return loc;
}




