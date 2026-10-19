#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include "../../core/kernel.h"
#include "../../core/ilwiscontext.h"
#include "../../core/catalog/catalog.h"
#include "../../core/version.h"

#include "../../core/ilwisobjects/ilwisdata.h"
#include "../../core/ilwisobjects/operation/operationmetadata.h"
#include "../../core/ilwisobjects/operation/symboltable.h"
#include "../../core/ilwisobjects/operation/commandhandler.h"
#include "../../core/ilwisobjects/operation/operationExpression.h"
#include "../../core/ilwisobjects/operation/operation.h"

#include "../../core/ilwisobjects/ilwisobject.h"

#include "../../core/ilwisobjects/domain/domain.h"
#include "../../core/ilwisobjects/domain/datadefinition.h"
#include "../../core/ilwisobjects/table/columndefinition.h"
#include "../../core/ilwisobjects/table/table.h"
#include "../../core/ilwisobjects/table/attributedefinition.h"

#include "../../core/ilwisobjects/coverage/raster.h"
#include "../../core/ilwisobjects/coverage/coverage.h"

#include "../../core/ilwisobjects/coverage/rastercoverage.h"

#include "../../core/util/box.h"

#include "../../core/ilwisobjects/coverage/featurecoverage.h"
#include "../../core/ilwisobjects/coverage/feature.h"
#include "../../core/kernel.h"
#include "workflownode.h"
#include "conditionNode.h"
#include "operationExpression.h"
#include "workflownode.h"
#include "operationnode.h"
#include "workflow.h"

#include "pythonapi_object.h"
#include "pythonapi_engine.h"
#include "pythonapi_collection.h"
#include "pythonapi_rastercoverage.h"
#include "pythonapi_featurecoverage.h"
#include "pythonapi_pyobject.h"
#include "pythonapi_catalog.h"
#include "pythonapi_table.h"
#include "pythonapi_booleanobject.h"

using namespace pythonapi;

std::map<quint64, Ilwis::IWorkflow> Engine::_workflows;

Engine::Engine(){
}

qint64 Engine::_do2(std::string output_name, std::string operation, std::string c3, std::string c4, std::string c5,std::string c6, std::string c7, std::string c8, std::string c9, std::string c10, std::string c11){
    Ilwis::SymbolTable symtbl;
    Ilwis::ExecutionContext ctx;
    ctx.clear();
    //is no internal result name is given it will look like operation_id
    //but the id is to be added afterwards
    bool rename = false;
    if (output_name.empty()){
        output_name = operation;
        rename = true;
    }
    QString command;
    if (!c3.empty()){
        if(!c4.empty()){
            if(!c5.empty()){
                if(!c6.empty()){
                    if(!c7.empty()){
                        if(!c8.empty()){
                            if(!c9.empty()){
                                if (!c10.empty()){
                                    if (!c11.empty()){
                                        command = QString("script %1=%2(%3,%4,%5,%6,%7,%8,%9,").arg(output_name.c_str(),operation.c_str(),c3.c_str(),c4.c_str(),c5.c_str(),c6.c_str(),c7.c_str(),c8.c_str(),c9.c_str()) + QString("%1,%2)").arg(c10.c_str(),c11.c_str());
                                    }else{
                                        command = QString("script %1=%2(%3,%4,%5,%6,%7,%8,%9,").arg(output_name.c_str(),operation.c_str(),c3.c_str(),c4.c_str(),c5.c_str(),c6.c_str(),c7.c_str(),c8.c_str(),c9.c_str()) + QString("%1)").arg(c10.c_str());
                                    }
                                }else{
                                    command = QString("script %1=%2(%3,%4,%5,%6,%7,%8,%9)").arg(output_name.c_str(),operation.c_str(),c3.c_str(),c4.c_str(),c5.c_str(),c6.c_str(),c7.c_str(),c8.c_str(),c9.c_str());
                                }
                            }else{
                                command = QString("script %1=%2(%3,%4,%5,%6,%7,%8)").arg(output_name.c_str(),operation.c_str(),c3.c_str(),c4.c_str(),c5.c_str(),c6.c_str(),c7.c_str(),c8.c_str());
                            }
                        }else{
                            command = QString("script %1=%2(%3,%4,%5,%6,%7)").arg(output_name.c_str(),operation.c_str(),c3.c_str(),c4.c_str(),c5.c_str(),c6.c_str(),c7.c_str());
                        }
                    }else{
                        command = QString("script %1=%2(%3,%4,%5,%6)").arg(output_name.c_str(),operation.c_str(),c3.c_str(),c4.c_str(),c5.c_str(),c6.c_str());
                    }
                }else{
                    command = QString("script %1=%2(%3,%4,%5)").arg(output_name.c_str(),operation.c_str(),c3.c_str(),c4.c_str(),c5.c_str());
                }
            }else{
                command = QString("script %1=%2(%3,%4)").arg(output_name.c_str(),operation.c_str(),c3.c_str(),c4.c_str());
            }
        }else{
            command = QString("script %1=%2(%3)").arg(output_name.c_str(),operation.c_str(),c3.c_str());
        }
    }else{
        command = QString("script %1=%2").arg(output_name.c_str(),operation.c_str());
    }
    if (Ilwis::commandhandler()->execute(command,&ctx, symtbl) && !ctx._results.empty()){
        //std::vector<Object*> results;
        for (int i = 0; i < ctx._results.size(); ++i) {
            Ilwis::Symbol result = symtbl.getSymbol(ctx._results[i]);
            if (result._type == itRASTER){
                if (result._var.canConvert<Ilwis::IRasterCoverage>()){
                    Ilwis::IRasterCoverage obj (result._var.value<Ilwis::IRasterCoverage>());
                    return obj->id();
                }
            }else if (result._type == itFEATURE){
                if (result._var.canConvert<Ilwis::IFeatureCoverage>()){
                    Ilwis::IFeatureCoverage obj (result._var.value<Ilwis::IFeatureCoverage>());
                    return obj->id();
                }
            }else if (result._type == itCOORDSYSTEM){
                if (result._var.canConvert<Ilwis::ICoordinateSystem>()){
                    Ilwis::ICoordinateSystem obj (result._var.value<Ilwis::ICoordinateSystem>());
                    return obj->id();
                }
            }else if (result._type == itGEOREF){
                if (result._var.canConvert<Ilwis::IGeoReference>()){
                    Ilwis::IGeoReference obj (result._var.value<Ilwis::IGeoReference>());
                    return obj->id();
                }
            }else if (result._type & itTABLE){
                if (result._var.canConvert<Ilwis::ITable>()){
                    Ilwis::ITable obj (result._var.value<Ilwis::ITable>());
                    return obj->id();
                }
            } else {
                auto list = result._var.toList();
                if ( list.size() > 0)
                    return list[0].toInt();
                return iUNDEF;
            }
        }
        /*
        if (results.size() == 0)
            throw Ilwis::ErrorObject(QString("couldn't handle return type of \"%1\"").arg(command.mid(8 + output_name.size())));
        else if (results.size() == 1)
            return results[0];
        else {
            return new Collection(results);
        }
        */
    }else{
        QString filter = QString("(type=%1 or type=%2)").arg(itSINGLEOPERATION).arg(itWORKFLOW);
        std::vector<Ilwis::Resource> ops = Ilwis::mastercatalog()->select(filter);
        bool found = false;
        for(auto it = ops.begin(); it != ops.end(); it++){
            if (it->name().toStdString() == operation) {
                found = true;
                break;
            }
        }
        if (found)
            throw Ilwis::ErrorObject(QString("Failed to execute command \"%1\"; Please check the parameters provided.").arg(command.mid(8 + output_name.size())));
        else
            throw Ilwis::ErrorObject(QString("Command \"%1\" does not exist; See ilwis.Engine.operations() for the full list.").arg(operation.c_str()));
    }
    return iUNDEF;
}

Object* Engine::_do(std::string output_name, std::string operation, std::string c3, std::string c4, std::string c5,std::string c6, std::string c7, std::string c8, std::string c9, std::string c10, std::string c11){
    Ilwis::SymbolTable symtbl;
    Ilwis::ExecutionContext ctx;
    ctx.clear();
    //is no internal result name is given it will look like operation_id
    //but the id is to be added afterwards
    bool rename = false;
    if (output_name.empty()){
        output_name = operation + "_object_" + QString::number(Ilwis::Identity::newAnonymousId()).toStdString();
        rename = true;
    }
    QString command;
    if (!c3.empty()){
        c3 = addQuotesIfNeeded(c3);
        if(!c4.empty()){
            c4 = addQuotesIfNeeded(c4);
            if(!c5.empty()){
                c5 = addQuotesIfNeeded(c5);
                if(!c6.empty()){
                    c6 = addQuotesIfNeeded(c6);
                    if(!c7.empty()){
                        c7 = addQuotesIfNeeded(c7);
                        if(!c8.empty()){
                            c8 = addQuotesIfNeeded(c8);
                            if(!c9.empty()){
                                c9 = addQuotesIfNeeded(c9);
                                if (!c10.empty()){
                                    c10 = addQuotesIfNeeded(c10);
                                    if (!c11.empty()){
                                        c11 = addQuotesIfNeeded(c11);
                                        command = QString("script %1=%2(%3,%4,%5,%6,%7,%8,%9,").arg(output_name.c_str(),operation.c_str(),c3.c_str(),c4.c_str(),c5.c_str(),c6.c_str(),c7.c_str(),c8.c_str(),c9.c_str()) + QString("%1,%2)").arg(c10.c_str(),c11.c_str());
                                    }else{
                                        command = QString("script %1=%2(%3,%4,%5,%6,%7,%8,%9,").arg(output_name.c_str(),operation.c_str(),c3.c_str(),c4.c_str(),c5.c_str(),c6.c_str(),c7.c_str(),c8.c_str(),c9.c_str()) + QString("%1)").arg(c10.c_str());
                                    }
                                }else{
                                    command = QString("script %1=%2(%3,%4,%5,%6,%7,%8,%9)").arg(output_name.c_str(),operation.c_str(),c3.c_str(),c4.c_str(),c5.c_str(),c6.c_str(),c7.c_str(),c8.c_str(),c9.c_str());
                                }
                            }else{
                                command = QString("script %1=%2(%3,%4,%5,%6,%7,%8)").arg(output_name.c_str(),operation.c_str(),c3.c_str(),c4.c_str(),c5.c_str(),c6.c_str(),c7.c_str(),c8.c_str());
                            }
                        }else{
                            command = QString("script %1=%2(%3,%4,%5,%6,%7)").arg(output_name.c_str(),operation.c_str(),c3.c_str(),c4.c_str(),c5.c_str(),c6.c_str(),c7.c_str());
                        }
                    }else{
                        command = QString("script %1=%2(%3,%4,%5,%6)").arg(output_name.c_str(),operation.c_str(),c3.c_str(),c4.c_str(),c5.c_str(),c6.c_str());
                    }
                }else{
                    command = QString("script %1=%2(%3,%4,%5)").arg(output_name.c_str(),operation.c_str(),c3.c_str(),c4.c_str(),c5.c_str());
                }
            }else{
                command = QString("script %1=%2(%3,%4)").arg(output_name.c_str(),operation.c_str(),c3.c_str(),c4.c_str());
            }
        }else{
            command = QString("script %1=%2(%3)").arg(output_name.c_str(),operation.c_str(),c3.c_str());
        }
    }else{
        command = QString("script %1=%2").arg(output_name.c_str(),operation.c_str());
    }
    if (Ilwis::commandhandler()->execute(command,&ctx, symtbl)){
        std::vector<Object*> results;
        for (int i = 0; i < ctx._results.size(); ++i) {
            Ilwis::Symbol result = symtbl.getSymbol(ctx._results[i]);
            if (result._type == itRASTER){
                if (result._var.canConvert<Ilwis::IRasterCoverage>()){
                    Ilwis::IRasterCoverage obj (result._var.value<Ilwis::IRasterCoverage>());
                    if (rename)
                        obj->name(QString("%1_%2").arg(operation.c_str()).arg(obj->id()));
                    results.push_back(new RasterCoverage(obj));
                }
            }else if (result._type == itFEATURE){
                if (result._var.canConvert<Ilwis::IFeatureCoverage>()){
                    Ilwis::IFeatureCoverage obj (result._var.value<Ilwis::IFeatureCoverage>());
                    if (rename)
                        obj->name(QString("%1_%2").arg(operation.c_str()).arg(obj->id()));
                    results.push_back(new FeatureCoverage(obj));
                }
            }else if (result._type == itCOORDSYSTEM){
                if (result._var.canConvert<Ilwis::ICoordinateSystem>()){
                    Ilwis::ICoordinateSystem obj (result._var.value<Ilwis::ICoordinateSystem>());
                    if (rename)
                        obj->name(QString("%1_%2").arg(operation.c_str()).arg(obj->id()));
                    results.push_back(new CoordinateSystem(obj));
                }
            }else if (result._type == itGEOREF){
                if (result._var.canConvert<Ilwis::IGeoReference>()){
                    Ilwis::IGeoReference obj (result._var.value<Ilwis::IGeoReference>());
                    if (rename)
                        obj->name(QString("%1_%2").arg(operation.c_str()).arg(obj->id()));
                    results.push_back(new GeoReference(obj));
                }
            }else if (result._type & itTABLE){
                if (result._var.canConvert<Ilwis::ITable>()){
                    Ilwis::ITable obj (result._var.value<Ilwis::ITable>());
                    if (rename)
                        obj->name(QString("%1_%2").arg(operation.c_str()).arg(obj->id()));
                    results.push_back(new Table(obj));
                }
            } else if (result._type == itBOOL){

            }
        }
        if (results.size() == 0)
            return new BooleanObject();
        else if (results.size() == 1)
            return results[0];
        else {
            return new Collection(results);
        }
    }else{
        QString filter = QString("(type=%1 or type=%2)").arg(itSINGLEOPERATION).arg(itWORKFLOW);
        std::vector<Ilwis::Resource> ops = Ilwis::mastercatalog()->select(filter);
        bool found = false;
        for(auto it = ops.begin(); it != ops.end(); it++){
            if (it->name().toStdString() == operation) {
                found = true;
                break;
            }
        }
        if (found) {
            QString lastIssue = Ilwis::kernel()->issues()->lastIssue();
            if (lastIssue.length() > 0)
                lastIssue = "\nAdditional info: " + lastIssue;
            throw Ilwis::ErrorObject(QString("Failed to execute command \"%1\"; Please check the parameters provided.%2").arg(command.mid(8 + output_name.size())).arg(lastIssue));
        } else
            throw Ilwis::ErrorObject(QString("Command \"%1\" does not exist; See ilwis.Engine.operations() for the full list.").arg(operation.c_str()));
    }
}

void Engine::_setWorkingCatalog(const std::string& location) {
    setWorkingCatalog(location);
}

void Engine::setWorkingCatalog(const std::string& location) {
    QString loc (QString::fromStdString(location));
    loc.replace('\\','/');
    // if it is file:// (or http:// etc) leave it untouched; if not, append file:// and the working catalog path if it is missing
    if (loc.indexOf("://") < 0) {
        int pos = loc.indexOf('/');
        if (pos > 0) { // full path starting with drive-letter (MS-DOS-style)
            loc = "file:///" + loc;
            if (loc.endsWith('/')) // workaround an IlwisObjects problem that scans the folder twice if it ends with a slash
                loc = loc.left(loc.length() - 1);
        } else if (pos == 0) { // full path starting with path-separator (UNIX-style)
            loc = "file://" + loc;
            if (loc.endsWith('/'))
                loc = loc.left(loc.length() - 1);
        }
    }

    Ilwis::ICatalog cat;
    cat.prepare(loc);
    if(cat.isValid()){
        Ilwis::context()->setWorkingCatalog(cat);
        Ilwis::mastercatalog()->addContainer(QUrl(loc));
    }else
        throw Ilwis::ErrorObject(QString("invalid container location: '%1'").arg(location.c_str()));
}

std::string Engine::getLocation(){
    Ilwis::ICatalog cat = Ilwis::context()->workingCatalog();
    QUrl location = cat->filesystemLocation();
    return location.toString().toStdString();
}

PyObject* Engine::operations(const std::string& mode){
    QString filter = QString("(type=%1 or type=%2)").arg(itSINGLEOPERATION).arg(itWORKFLOW);
    std::vector<Ilwis::Resource> ops = Ilwis::mastercatalog()->select(filter);
    PyObject* list = newPyTuple(ops.size());
    int i = 0;
    for(auto it = ops.begin(); it != ops.end(); it++){
        if ( mode == ""){
            if (!setTupleItem(list, i++, PyUnicodeFromString(it->name().toStdString().data()))){
                throw Ilwis::ErrorObject(QString("internal conversion error while trying to add '%1' to list of attributes").arg( it->name()));
            }
        }else if ( mode == "unique"){
            QString nameidtup = it->name() + "|" + QString::number(it->id());
            if (!setTupleItem(list, i++, PyUnicodeFromString(nameidtup.toStdString().data()))){
                throw Ilwis::ErrorObject(QString("internal conversion error while trying to add '%1' to list of attributes").arg( it->name()));
            }
        }
    }
    return list;
}

void pythonapi::Engine::addParameter(QString &parms, int currIndex, std::map<QString, QString>& parm, bool input)
{
    parms += "\"parameter_" + QString((input ? "i" : "o")) + "_"+ QString::number(currIndex) + "\" : {";
    for(auto it = parm.begin() ; it != parm.end(); ++it){
        QString value = it->first != "type" ? it->second : Ilwis::TypeHelper::type2names(it->second.toLongLong());
        parms += "\"" + it->first + "\":\"" + value + "\",";
    }
    parms += QString("\"index\":") + "\"" + QString::number(currIndex) + "\",";
    parms += QString("\"input\":") + (input ? "\"yes\"" : "\"no\"");
    parm.clear();

    parms += "},";
}
std::map<QString, QStringList> extractFixedValues(const QString& input) {
    std::map<QString, QStringList> results;

    // 1. Identify the boundaries of the parameters
    int firstParen = input.indexOf('(');
    int lastParen = input.lastIndexOf(')');

    if (firstParen == -1 || lastParen == -1 || lastParen <= firstParen) {
        return results;
    }

    // Extract everything between the parentheses
    QString content = input.mid(firstParen + 1, lastParen - firstParen - 1);
    if (content.isEmpty()) return results;

    // 2. Split into individual definitions (e.g., "name=a|b")
    QStringList paramDefs = content.split(',', QString::SkipEmptyParts);

    for (const QString& def : paramDefs) {
        QString trimmedDef = def.trimmed();

        if (trimmedDef.contains('=')) {
            // 3. Handle key-value pairs: name=val1|val2
            QStringList parts = trimmedDef.split('=');
            QString key = parts.value(0).trimmed();
            QStringList values;

            if (parts.size() > 1) {
                values = parts.at(1).split('|', QString::SkipEmptyParts);
                // Clean up whitespace for each allowed value
                for(QString& val : values) val = val.trimmed();
            }
            results[key] = values;
        } else {
            // 4. Handle standalone parameter names
            results[trimmedDef] = QStringList();
        }
    }

    return results;
}

std::string Engine::operationMetaData(const std::string &name, const std::string &element){
    QString filter = QString("(type=%1 or type=%2)").arg(itSINGLEOPERATION).arg(itWORKFLOW);
    std::vector<Ilwis::Resource> ops = Ilwis::mastercatalog()->select(filter);
    bool isid;
    QString::fromStdString(name).toDouble(&isid);
    bool all = name == "*";
    int count = 0;
    QString ret = all ? "[" : "";
    for(auto it = ops.begin();it != ops.end(); it++){
        if (isid || all){
            if (!all){
                if (it->id() != QString::fromStdString(name).toULongLong())
                    continue;
            }
            if ( count > 0) ret += ",";
            if ( it->name() == "atan"){
                qDebug() << "stop";
            }
            ret += "{" ;
            ret += "\"id\" : \""+ QString::number(it->id()) + "\",";
            ret += "\"name\": \""+ it->name() + "\",";
            ret += "\"description\": \""+ it->description() + "\",";
            ret += "\"url\" :\""+ it->url().toString() + "\",";
            ret += "\"rawurl\": \""+ it->url(true).toString() + "\",";
            ret += "\"container\" : \""+ it->container().toString() + "\",";
            ret += "\"rawcontainer\":\""+ it->container(true).toString() + "\",";
            ret += "\"size\" :\"" + QString::number(it->size()) + "\",";
            ret += "\"dimensions\": \"" + it->dimensions() + "\",";
            ret += "\"type\": \"" + Ilwis::TypeHelper::type2name( it->ilwisType()) + "\",";

            QString parms = "\"parameters\" : ";
            std::map<int, std::map<QString, QString>> iparms;
            auto props = it->properties();
            auto keys = props.keys();
            std::sort(keys.begin(), keys.end(), [](const auto& a, const auto& b) {
                    return a < b;
                });

            QJsonObject root;
            std::map<QString, QStringList> fixeds = extractFixedValues(props["syntax"].toString());

            for (const QString& key : keys) {
                // 1. Split the string by underscores
                // parameter_X_Y -> ["parameter", "X", "Y"]
                QStringList parts = key.split('_');

                // Ensure we have at least 3 parts (parameter, X, and Y)
                if (parts.size() >= 3) {
                    // Reconstruct the unique parameter group (e.g., "parameter_1")
                    QString paramGroup = parts[0] == "pin" ? "parameter_i_" + parts[1] : "parameter_o_" + parts[1];

                    // The rest is the name (Y). Join back in case Y contains underscores.
                    QString attrName = parts.mid(2).join('_');

                    // 2. Retrieve existing group object or create a new one
                    QJsonObject paramObject = root.value(paramGroup).toObject();

                    // 3. Insert the Y name and the value from prop()
                    QString valueTypes;
                    if ( parts[2] != "type")
                        paramObject.insert(attrName, props[key].toString().trimmed());
                    else{
                         valueTypes = Ilwis::TypeHelper::type2names(props[key].toLongLong());
                         paramObject.insert("type", valueTypes);
                    }
                    paramObject.insert("index", parts[1]);
                    paramObject.insert("input", parts[0] == "pin" ? "yes" : "no");
                    if ( parts[0] == "pout" && parts[2] == "type"){
                        root.insert("boolean_operation", valueTypes.indexOf("bool") >= 0 ? "yes": "no");
                    }
                    if ( key.indexOf("parm_name") > 0){
                        QString pname = props[key].toString();
                        paramObject.insert("valueset", fixeds[pname].join('|'));
                    }

                    // 4. Put the updated group back into the root object
                    root.insert(paramGroup, paramObject);
                }else {
                    root.insert(key, props[key].toString().trimmed());
                }
            }

            QJsonDocument doc(root);
            parms +=  QString::fromLocal8Bit( doc.toJson(QJsonDocument::Compact)) + "}";

            ret += parms;

            ++count;
        }else if (QString::fromStdString(name).compare(it->name(),Qt::CaseInsensitive) == 0){
            if(!ret.isEmpty())
                ret.append("; ");
            if (element.compare("description") == 0)
                ret.append(it->description());
            else if ( element.compare("id") == 0)
                ret.append(QString::number(it->id()));
            else
                ret.append((*it)[element.c_str()].toString());
        }
    }
    ret += all ? "]" : "";
    return ret.toStdString();
}

std::string Engine::_operationMetaData(const std::string &name, const std::string &element1, int ordinal, const std::string &element2)
{
    std::string element;
    if ( element1 == "input"){
        element = "pin_" + std::to_string(ordinal);
    }else if ( element1 == "output"){
        element = "pout_" + std::to_string(ordinal);

    }
    if ( element != "" ){
        if ( element2 == "description"){
           element += "_desc";
        }else
            element += "_" + element2;
    }else
        element = element1;


    auto retValue =  operationMetaData(name, element);
    if (element1 == "type" || element2 == "type"){
         QString typeNames = QString::fromStdString(retValue);
         QString tpNames = Ilwis::TypeHelper::type2names(typeNames.toULongLong()," or ");
         retValue = tpNames.toStdString();
     }
     return retValue;
}

std::string Engine::_updateWorkflow(const std::string& json){
    QJsonDocument doc = QJsonDocument::fromJson(QString::fromStdString(json).toUtf8());
    if (!doc.isNull()) {
        Ilwis::IWorkflow wf;
        QJsonValue v = doc["workflow_id"];
        auto v2 = v.toInt();
        if ( v2 == -1){
            wf.prepare(Ilwis::Resource(QString("wf%1").arg((Ilwis::Time::now().toString(true, itUNKNOWN, true))), itWORKFLOW));
            Ilwis::mastercatalog()->addItems({wf->resource()});
            _workflows[v2] = wf;
        }else {
           wf.prepare(v2);
        }
        if ( wf.isValid()){
            auto p = wf->id();
            QString action= doc["action"].toString();
            QString s = QString("{\"workflow_id\" :%1").arg(p);
            if ( action == "addnode"){
                QString opid = doc["id"].toString();
                Ilwis::SPOperationNode opNode;
                auto opNodePtr = new Ilwis::OperationNode(opid.toLongLong());
                opNode.reset(opNodePtr);
                wf->addNode(opNode);
                wf->createMetadata();
                s += ",\"node_id\":" + QString::number(opNodePtr->id());


            }else if ( action == "addconnection"){
                int parmInIndex = doc["paramInIndex"].toInt();
                if ( parmInIndex != 9999){ // 9999 is magic number for indicating the test part of condition
                    int sourceNodeId= doc["source"].toString().toInt();
                    int targetNodeId= doc["target"].toString().toInt();
                    int parmOutIndex = doc["paramOutIndex"].toInt();
                    wf->addFlow(sourceNodeId, targetNodeId, parmInIndex,parmOutIndex);
                    wf->createMetadata();
                }else {
                    int sourceNodeId= doc["source"].toString().toInt();
                    int targetNodeId= doc["target"].toString().toInt();
                    auto testNode = wf->nodeById(sourceNodeId);
                    auto condNode = std::static_pointer_cast<Ilwis::WorkFlowCondition>(wf->nodeById(targetNodeId));
                    condNode->addTest(testNode);
                }
                s += QString(",\"node_id\":") + doc["source"].toString();
            }else if ( action == "fixed"){
                int nodeid = doc["node_id"].toString().toInt();
                auto parms = doc["parameters"].toObject();
                int count = 0;
                for (auto it = parms.begin(); it != parms.end(); ++it) {
                    QString key = it.key();
                    auto value = it.value().toObject();
                    QString fixed =value["fixed_value"].toString();
                    if ( fixed != ""){
                        wf->setFixedParameter(fixed, nodeid, count);
                    }
                    ++count;
                }
            }else if ( action == "addcondition"){
                Ilwis::SPWorkFlowNode conNode;
                auto conNodePtr = new Ilwis::WorkFlowCondition();
                conNodePtr->nodeId(wf->generateId());
                conNode.reset(conNodePtr);
                wf->addNode(conNode);
                s += ",\"node_id\":" + QString::number(conNode->id());
            }
            s += "}";
            return s.toStdString();
        }

    }
    return "{\'workflow_id\': -1}";
}

PyObject* Engine::_catalogItems(quint64 filter){
    Ilwis::ICatalog cat = Ilwis::context()->workingCatalog();
    std::vector<Ilwis::Resource> resVec = cat->items();
    std::vector<Ilwis::Resource> result;
    if ( filter != itUNKNOWN){
        for(Ilwis::Resource& res : resVec){
            if ( hasType(res.ilwisType(), filter)){
                result.push_back(res);
            }
        }
    }else
        result = resVec;
    PyObject* tup = newPyTuple(result.size());
    int i = 0;
    for(auto it = result.begin();it != result.end(); it++){
        if (!setTupleItem(tup, i++, PyUnicodeFromString(it->name().toStdString().data()))){
            throw Ilwis::ErrorObject(QString("internal conversion error while trying to add '%1' to list of files").arg( it->name()));
        }
    }
    return tup;
}

std::string Engine::_version()
{
    return Ilwis::kernel()->version()->verionNumber().toStdString();
}

PyObject *Engine::_operations(const std::string &mode)
{
    return operations(mode);
}

std::string Engine::addQuotesIfNeeded(std::string parameter) {
    if (parameter.front() != '\'' && parameter.back() != '\'') { // if it does not already have quotes
        double d;
        if (sscanf(parameter.c_str(), "%lf", &d) <= 0) { // if it is not a number (int,float). Do we also expect hex/oct numbers here?
            parameter = "'" + parameter + "'";
        }
    }
    return parameter;
}

void Engine::_showProgress(bool show) {
    if (show)
        Ilwis::context()->runMode(Ilwis::rmCOMMANDLINE);
    else
        Ilwis::context()->runMode(Ilwis::rmCOMMANDLINE | Ilwis::rmNOUI);
}

void Engine::_removeObject(quint64 id) {
    Ilwis::Resource res = Ilwis::mastercatalog()->id2Resource(id);
    if ( res.isValid())
        Ilwis::mastercatalog()->removeItems({res});
}

std::string Engine::contextProperty(const std::string& propName){
    QString pn = QString::fromStdString(propName);
    pn = pn.toLower();
    QString result = sUNDEF;
    if ( propName == "cachelocation"){
        result = Ilwis::context()->cacheLocation().toLocalFile();
    }else if (propName == "ilwisfolder"){
        result = Ilwis::context()->ilwisFolder().filePath();
    }else if ( propName == "workingcatalog"){
       Ilwis::ICatalog cat = Ilwis::context()->workingCatalog();
       QUrl location = cat->filesystemLocation();
       result = location.toString();
    }else if ( propName == "systemcatalog"){
        const Ilwis::ICatalog cat = Ilwis::context()->systemCatalog();
        QUrl location = cat->filesystemLocation();
        result = location.toString();
    }else if ( propName == "loglocation"){
      result = Ilwis::kernel()->issues()->location();
    }
    return result.toStdString();
}

void Engine::_setContextProperty(const std::string& propName, const std::string& value){
    QString pn = QString::fromStdString(propName);
    pn = pn.toLower();
    QString result = sUNDEF;
    if ( propName == "cachelocation"){
        Ilwis::context()->setCacheLocation(QString::fromStdString(value));
    }else if (propName == "internalcatalog"){
        Ilwis::context()->setInternalCatalog(QString::fromStdString(value));
    } else if ( propName == "initlogger"){
        Ilwis::kernel()->issues()->initLogger(QString::fromStdString(value));
    }
}
