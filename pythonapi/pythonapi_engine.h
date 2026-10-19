#ifndef PYTHONAPI_ENGINE_H
#define PYTHONAPI_ENGINE_H

#include "pythonapi_object.h"

typedef struct _object PyObject;
typedef std::tuple<QString, QStringList> FParameters;

namespace Ilwis{
    template<class T> class IlwisData;
    class Workflow;
    typedef IlwisData<Workflow> IWorkflow;
}
namespace pythonapi {
    class Catalog;
    class Engine{
    public:
        Engine();
        static qint64 _do2(std::string output_name, std::string operation,std::string c3 = "",std::string c4 = "",std::string c5 = "",std::string c6 = "",std::string c7="", std::string c8="", std::string c9="", std::string c10="", std::string c11="");
        static Object* _do(std::string output_name, std::string operation,std::string c3 = "",std::string c4 = "",std::string c5 = "",std::string c6 = "",std::string c7="", std::string c8="", std::string c9="", std::string c10="", std::string c11="");
        static void setWorkingCatalog(const std::string& location);
        static void _setWorkingCatalog(const std::string& location);
        static std::string getLocation();
        static PyObject* operations(const std::string &mode);
        static std::string operationMetaData(const std::string& name, const std::string &element = "syntax");
        static std::string _operationMetaData(const std::string& name, const std::string &element1 = "syntax", int ordinal=-1, const std::string &element2 = "");
        static std::string _updateWorkflow(const std::string& json);
        static PyObject* _catalogItems(quint64 filter);
        static std::string _version();
        static PyObject* _operations(const std::string& mode="");
        static void _showProgress(bool show);
        static void _removeObject(quint64 objectid);
        static std::string contextProperty(const std::string& propName);
        static void _setContextProperty(const std::string& propName, const std::string& value);
    private:
        static std::string addQuotesIfNeeded(std::string parameter);
        static std::map<quint64, Ilwis::IWorkflow> _workflows;
        static void addParameter(QString &parms, int currIndex, std::map<QString, QString>& parm, bool input=true);
    };

}
#endif // PYTHONAPI_ENGINE_H
