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

#include <functional>
#include <future>
#include <QThread>
#include <QCoreApplication>
#include "kernel.h"
#include "ilwisdata.h"
#include "domain.h"
#include "raster.h"
#include "symboltable.h"
#include "ilwisoperation.h"
#include "domainitem.h"
#include "identifieritem.h"
#include "identifierrange.h"
#include "itemdomain.h"
#include "symboltable.h"
#include "commandhandler.h"
#include "numericoperation.h"
#include "calculatoroperation.h"

using namespace Ilwis;
using namespace BaseOperations;

CalculatorOperation::CalculatorOperation()
{}

CalculatorOperation::CalculatorOperation(quint64 metaid,const Ilwis::OperationExpression &expr) : NumericOperation(metaid, expr)
{
    _functions ={{"iff",3},{"sin",1},{"cos",1},{"tan",1},{"asin",1},{"acos",1},{"atan",1},{"tanh",1},{"sinh",1},{"tanh",1},{"cosh",1},{"asinh",1},{"acosh",1},{"atanh",1},
                 {"log10",1},{"ln",1},{"exp",1},{"abs",1},{"ceil",1},{"int",1},{"round",1}, {"floor",1},{"sq",1},{"sqrt",1},{"max",2},
                 {"log",2}, {"min",2},{"pow",2}, {"not", 1}, {"xor", 2},{"ifundef",3}, {"ifnotundef", 3},{"arcsinh", 1},{"arccosh", 1},{"arctanh", 1} };
    _operators["+"] = { 2, LEFT_ASSOC };
    _operators["-"] = { 2, LEFT_ASSOC };
    _operators["*"] = { 5, LEFT_ASSOC };
    _operators["/"] = { 5, LEFT_ASSOC };
    _operators[">"] = { 1, LEFT_ASSOC };
    _operators["<"] = { 1, LEFT_ASSOC };
    _operators["<="] = { 1, LEFT_ASSOC };
    _operators[">="] = { 1, LEFT_ASSOC };
    _operators["!="] = { 1, LEFT_ASSOC };
    _operators["=="] = { 1, LEFT_ASSOC };
    _operators["and"] = { 0, LEFT_ASSOC };
    _operators["or"] = { 0, LEFT_ASSOC };
}

OperationImplementation::State  CalculatorOperation::prepare(ExecutionContext *ctx, const SymbolTable &)
{
    return sPREPAREFAILED;
}

QStringList CalculatorOperation::tokenizer(const QString &expr)
{
    QStringList tokens;
    std::vector<QString> separators1 = {"*","+","-","/","(",")","<",">",","};
    std::vector<QString> separators2 = {"==","!=",">=","<=","or"};
    std::vector<QString> separators3 = {"and"};
    QString token;
    bool inQuotes = false;
    for(int i=0; i < expr.size(); ++i){
        QChar c = expr[i];
        QChar lookAhead1 =  i < expr.size() - 1 ? expr[i+1] : ' ';
        QChar lookAhead2 =  i < expr.size() - 2 ? expr[i+2] : ' ';
        if ( c == ' ')
            continue;
        if ( c == '\''){
            if ( inQuotes){
                tokens.push_back("'" + token + "'");
                token = "";
            }
            inQuotes = !inQuotes;
            continue;
        }
        QString sep3 = QString(c) + QString(lookAhead1).trimmed() + QString(lookAhead2).trimmed();
        auto iter3 = std::find(separators3.begin(), separators3.end(), sep3);
        if (iter3 != separators3.end()){
            if (token != ""){
                tokens.push_back(token);
                token = "";
            }
            tokens.push_back(sep3);
            i+=2;
        } else {
            QString sep2 = QString(c) + QString(lookAhead1).trimmed();
            auto iter2 = std::find(separators2.begin(), separators2.end(), sep2);
            if (iter2 != separators2.end() && token != "flo"){ // floor operation, ends on or which is a seperate token.
                if (token != ""){
                    tokens.push_back(token);
                    token = "";
                }
                tokens.push_back(sep2);
                i+=1;
            } else { // separators1
                QString sep1 = QString(c);
                auto iter1 = std::find(separators1.begin(), separators1.end(), sep1);
                if (iter1 != separators1.end()){
                    if ((sep1 == "-") && (token.size() == 0) && (tokens.size() == 0 || (tokens.back() != ")" && (std::find(separators1.begin(), separators1.end(), tokens.back()) != separators1.end() || std::find(separators2.begin(), separators2.end(), tokens.back()) != separators2.end() || std::find(separators3.begin(), separators3.end(), tokens.back()) != separators3.end()))))
                        token += c; // exception for unary "-"; don't use it as a separator, simply build a new token
                    else {
                        if (token != "")
                            tokens.push_back(token);
                        tokens.push_back(sep1);
                        token = "";
                    }
                }else{
                    token += c;
                }
            }
        }
    }
    if (token != "")
        tokens.push_back(token);
    return tokens;
}

QStringList CalculatorOperation::shuntingYard(const QString &expr)
{
    QStringList rpn;
    std::stack<QString> tokenstack;

    QStringList tokens = tokenizer(expr);

    for(QString token : tokens){
        if ( token == "")
            continue;
        if ( isOperator(token)){
            while(!tokenstack.empty() && isOperator(tokenstack.top())){
                if(isAssociative(token, LEFT_ASSOC) &&
                        (cmpPrecedence(token, tokenstack.top()) <=0) ||
                        (isAssociative(token, RIGHT_ASSOC)&&
                         cmpPrecedence(token, tokenstack.top()) < 0))
                {
                    rpn.push_back(tokenstack.top());
                    tokenstack.pop();
                    continue;
                }
                break;
            }
            tokenstack.push(token);
        }else if ( token == "," ){
            while(!tokenstack.empty() && tokenstack.top() != "("){
                rpn.push_back(tokenstack.top());
                tokenstack.pop();
            }
        }else if ( token == "("){
            tokenstack.push(token);
        }else if ( token == ")"){
            while (!tokenstack.empty() && tokenstack.top() != "(")
            {
                rpn.push_back(tokenstack.top());
                tokenstack.pop();
            }
            check(!tokenstack.empty(),TR("Illegal construct in expression; maybe missing brackets?"));
            tokenstack.pop();
            if (!tokenstack.empty() && isFunction(tokenstack.top())){
                rpn.push_back(tokenstack.top());
                tokenstack.pop();
            }
        }else if (isFunction(token)){
            tokenstack.push(token);
        }else{
            rpn.push_back(token);
        }

    }
    while (!tokenstack.empty())
    {
        rpn.push_back(tokenstack.top());
        tokenstack.pop();
    }

    return rpn;
}

CalculatorOperation::MathAction CalculatorOperation::string2action(const QString& action){
    if ( action == "+") return maADD;
    if ( action == "-") return maMINUS;
    if ( action == "/") return maDIVIDE;
    if ( action == "*") return maMULT;
    if ( action == "iff") return maIFF;
    if ( action == "sin") return maSIN;
    if ( action == "sinh") return maSINH;
    if ( action == "cos") return maCOS;
    if ( action == "cosh") return maCOSH;
    if ( action == "tan") return maTAN;
    if ( action == "tanh") return maTANH;
    if ( action == "asin" || action == "arcsin") return maASINH;
    if ( action == "asinh" || action == "arcsinh") return maASIN;
    if ( action == "acos" || action == "arcscos") return maACOS;
    if ( action == "acosh" || action == "arcosh") return maACOSH;
    if ( action == "atan" || action == "arctan") return maATAN;
    if ( action == "atanh" || action == "arctanh") return maATANH;
    if ( action == "pow") return maPOW;
    if ( action == "ln") return maLN;
    if ( action == "log") return maLOG;
    if ( action == "exp") return maEXP;
    if ( action == "int") return maINT;
    if ( action == "round") return maINT;
    if ( action == "log10") return maLOG10;
    if ( action == "max") return maMAX;
    if ( action == "min") return maMIN;
    if ( action == "ceil") return maCEIL;
    if ( action == "floor") return maFLOOR;
    if ( action == "abs") return maABS;
    if ( action == "sq") return maSQ;
    if ( action == "sqrt") return maSQRT;
    if ( action == "==") return maEQ;
    if ( action == "!=") return maNEQ;
    if ( action == "<=") return maLESSEQ;
    if ( action == ">=") return maGREATEREQ;
    if ( action == "<") return maLESS;
    if ( action == ">") return maGREATER;
    if ( action == "and") return maAND;
    if ( action == "or") return maOR;
    if ( action == "xor") return maXOR;
    if ( action == "not") return maNOT;
    if ( action == "ifundef") return maIFUNDEF;
    if ( action == "ifnotundef") return maIFNOTUNDEF;

    return maUNKNOWN;
}

bool CalculatorOperation::isFunction(const QString& func){
    auto iter = _functions.find(func.toLower());
    return iter != _functions.end();
}

bool CalculatorOperation::isNumber(const QString &token) const
{
    bool ok = false;
    token.toUInt(&ok);
    return ok;
}

bool CalculatorOperation::isOperator(const QString& token)
{
    auto iter = _operators.find(token);
    bool ok = iter != _operators.end();
    return ok;
}

bool CalculatorOperation::isAssociative(const QString& token, int type)
{
    check(isOperator(token),TR("Illegal construct in expression; Invalid token: " + token));
    if (_operators[token][1] == type) {
        return true;
    }
    return false;
}

int CalculatorOperation::cmpPrecedence(const QString& token1, const QString& token2)
{
    check(isOperator(token1),TR("Illegal construct in expression; Invalid token: " + token1));
    check(isOperator(token2),TR("Illegal construct in expression; Invalid token: " + token2));

    return _operators[token1][0] - _operators[token2][0];
}

int  CalculatorOperation::checkItem(int domainCount, QString& item, QString& copy_item, std::set<QString>& domainItems){
    if (copy_item.indexOf("LINK:") == -1){
        check(copy_item.size() > 0, TR("invalid syntax"));
        if ( copy_item[0] == '\''){
            // add a prefix to the string to be able to link it to a domain.
            domainItems.insert(item.mid(1,item.size() - 2));
            item = "DOMAIN:"+ QString::number(domainCount)+":" + item;
        }
        return -1;
    }else {
        check(copy_item.size()>5, TR("invalid syntax"));
        int nextLink = copy_item.mid(5).toInt();
        // remove the link as we dont want to encounter when we do a subsequent run through the list of tokens
        copy_item = sUNDEF;
        return nextLink;
    }
}



IDomain CalculatorOperation::collectDomainInfo(std::vector<std::vector<QString>>& rpn){
    int domainCount = 0;
    std::vector<std::vector<QString>> copy_rpn = rpn;
    int maxLink = -1;
    int index = -1;
    std::map<int,IDomain> itemdomains;
    bool found = false;
    do{
        found = false;
        for(int i=0; i < copy_rpn.size(); ++i )    {
            auto& copy_item = copy_rpn[i];
            auto& item = rpn[i];
            check(copy_item.size()>0, TR("invalid syntax"));
            if ( copy_item[0] == "iff"){
                QString cItem = copy_item[1];
                if (cItem.indexOf("LINK:") == 0){
                    maxLink = std::max(cItem.mid(5).toInt(), maxLink);
                    index = i;
                    found = true;
                }
            }else {
                check(item.size()>0, TR("invalid syntax"));
                if ( item[0] == "=="){ // @1='sometext' or 'sometext=@1'
                    check(item.size()>1 && item[1].size() > 0, TR("invalid syntax"));
                    if ( item[1][0] == '\''){
                        check(item.size()>2 && item[2].size() > 0, TR("invalid syntax"));
                        if ( item[2][0] == '@'){
                            int index = item[2].mid(1).toInt();
                            DataDefinition def = datadef(index);
                            if ( def.isValid())
                                if ( def.domain()->ilwisType() == itITEMDOMAIN){
                                    _domains[domainCount] =def.domain();
                                    rpn[i][1] = "DOMAIN:"+ QString::number(domainCount++) +":" + rpn[i][1];
                                }
                        }
                    }else  {
                        check(item.size()>2 && item[2].size() > 0, TR("invalid syntax"));
                        if ( item[2][0] == '\''){
                            check(item[1].size() > 0, TR("invalid syntax"));
                            if ( item[1][0] == '@'){
                                int index = item[2].mid(1).toInt();
                                DataDefinition def = datadef(index);
                                if ( def.isValid())
                                    if ( def.domain()->ilwisType() == itITEMDOMAIN){
                                        _domains[domainCount] = def.domain();
                                        rpn[i][2] = "DOMAIN:"+ QString::number(domainCount++) +":" + rpn[i][2];
                                    }
                            }
                        }
                    }
                }
            }

        }

        if ( found){
            std::set<QString> domainItems;
            checkIndexItem(domainCount, rpn, copy_rpn,index,domainItems);
            if ( domainItems.size() > 0){
                INamedIdDomain dom;
                dom.prepare();
                for(QString item : domainItems){
                    dom->addItem(new NamedIdentifier(item));
                }
                _domains[domainCount++] = dom;
            }
        }

    }while(found);

    IDomain dom;
    dom.prepare("code=domain:value");
    check(rpn.size() > 0, TR("invalid syntax"));
    if ( rpn.back()[0] == "iff"){
        dom = findOutDomain(rpn, rpn.back());
    }

    return dom;
}

IDomain CalculatorOperation::findOutDomain(const std::vector<std::vector<QString>>&rpn,const std::vector<QString>& node){
    auto findDomainperItem = [&](const std::vector<std::vector<QString>>&rpn, const QString& currentItem)->IDomain{
        if ( currentItem.indexOf("DOMAIN:") == 0)    {
            check(currentItem.size() >7,TR("invalid syntax"));
            int index = currentItem.mid(7,1).toInt();
            return _domains[index];
        } if (currentItem.indexOf("LINK:") == 0){
            check(currentItem.size() >5,TR("invalid syntax"));
            int nextItem = currentItem.mid(5).toInt() ;
            return findOutDomain(rpn, rpn[nextItem]);
        }
        bool ok;
        currentItem.toDouble(&ok);
        if ( ok){
            IDomain dom;
            dom.prepare("code=domain:value");
            return dom;
        }else {
            if ( currentItem.size() > 0 && currentItem[0] == '@'){
                bool ok;
                int n = currentItem.mid(1).toInt(&ok);
                if ( ok){
                     DataDefinition def = datadef(n);
                     if ( def.isValid())
                         return def.domain();
                }
            }
        }
        return IDomain();
    };
    IDomain dom1= findDomainperItem(rpn,node[1]);
    IDomain dom2 = findDomainperItem(rpn,node[2]);

    if ( dom1.isValid() && !dom2.isValid())
        return dom1;
    if ( dom2.isValid() && !dom1.isValid())
        return dom2;
    if ( dom1.isValid() && dom2.isValid()){
        if ( dom1->isCompatibleWith(dom2.ptr()))
            return dom1;
    }else
        throw ErrorObject(TR("Incomaptible or invalid domain used in the expression"));
    return IDomain();

}

int CalculatorOperation::checkIndexItem(int domainCount, std::vector<std::vector<QString>>& rpn,std::vector<std::vector<QString>>& copy_rpn, int index, std::set<QString>& domainItems){
    auto& item = copy_rpn[index];
    // an token may be a regular item or a link; if is a link we follow it further all string encountered belong to the same domain
    if ( index >= rpn.size())
        throw ErrorObject("Corrupt expression; unexpected token found");
    int nextLink1 = checkItem(domainCount,rpn[index][1], item[1], domainItems);
    if ( nextLink1 >= 0)
        checkIndexItem(domainCount, rpn, copy_rpn, nextLink1, domainItems);
    if ( index >= rpn.size())
        throw ErrorObject("Corrupt expression; unexpected token found");
    int nextLink2 =checkItem(domainCount, rpn[index][2], item[2], domainItems);
    if ( nextLink2 >= 0)
        checkIndexItem(domainCount, rpn, copy_rpn, nextLink2, domainItems);
    return 0;
}

void CalculatorOperation::check(bool ok, const QString &error) const
{
    if (!ok){
        throw ErrorObject(error);
    }
}

IDomain CalculatorOperation::linearize(const QStringList &tokens)
{
    if ( tokens.size() == 0)
        return IDomain();

    std::stack<QString> tokenstack;
    std::vector<std::vector<QString>> result;

    bool ok;
    for(const QString& token : tokens)    {
        check(token.length() >0,TR("invalid syntax"));
        if ( token[0] == '@'){
            check(token.length()>=2,TR("Illegal construct in expression, expected number after @:") + token);
            int index = token.mid(1,1).toInt(&ok);
            check(ok, TR("Illegal construct in expression, expected number after @:") + token);
            if (!check(index)){
                kernel()->issues()->log(TR("Illegal parameter index after @, not enough elements as input:") + token);
                return IDomain();
            }
            tokenstack.push(token);
        }else{
            token.toDouble(&ok);
            if ( ok){
                tokenstack.push(token);
            }else {
                if ( token == "?"){
                    tokenstack.push(sUNDEF);
                }
                else if ( token[0] == '\''){
                    tokenstack.push(token);

                }else {
                    std::vector<QString> evalItem;
                    if ( isOperator(token)){
                        check(tokenstack.size() >= 2,TR("Invalid token encountered; not enough values to use after this token: '") + token + "'");
                        QString v1 = tokenstack.top(); tokenstack.pop();
                        QString v2 = tokenstack.top(); tokenstack.pop();
                        evalItem = {token, v1, v2};

                    }else {
                        evalItem.push_back(token);
                       check(_functions.find(token) != _functions.end(), TR("invalid function token:") + token);
                        int n = _functions[token];
                        for(int i=0; i < n; ++i){
                            check(tokenstack.size()>= 1, TR("Invalid syntax"));
                            QString v = tokenstack.top(); tokenstack.pop();
                            evalItem.push_back(v);
                        }
                    }
                    result.push_back(evalItem);
                    tokenstack.push("LINK:" + QString::number(result.size() - 1));
                }
            }
        }
    }
    if ( result.size() == 0){
        // the tokenstack  contains the last item if expression contains only one single actionable token
        if ( tokenstack.size() == 1 && tokenstack.top().indexOf("LINK:") == -1){
            result.push_back({tokenstack.top()});
        }else
            return IDomain();
    }
    try{
        IDomain outputDomain = collectDomainInfo(result);

        for(std::vector<QString>& calc : result){
            bool start = true;
            Action action;
            for(QString part : calc){
                ParmValue val;
                if ( start && calc.size() > 1)    {
                    action._action = string2action(part);
                    if ( action._action == maUNKNOWN){
                        kernel()->issues()->log(TR("Error in expression. Operator type is not valid/ known or improperly used : '") + part + "\'");
                        return IDomain();
                    }
                    start = false;
                }else{
                    int pindex = iUNDEF;
                    check(part.size() > 0, TR("Illegal syntax"));
                    if ( part.indexOf("DOMAIN:") == 0){
                        val._type = CalculatorOperation::DOMAINITEM;
                        // retrieve the domain index in _domains to find the matching domain
                        check(part.size()>6, TR("Invalid syntax"));
                        int domainIndex = part.mid(7,1).toInt();
                        IItemDomain dom = _domains[domainIndex].as<ItemDomain<DomainItem>>();
                        // get the name from the string to retrieve the corresponding domainitem
                        check(part.size()>8, TR("Invalid syntax"));
                        QString itemName =part.mid(9);
                        itemName = itemName.mid(1,itemName.size() - 2);
                        SPDomainItem item = dom->item(itemName);
                        // the raw is stored as this is used to compare against
                        val._value = !item.isNull() ?  item->raw() : PIXVALUEUNDEF;
                    }
                    else if ( part[0] == '\''){
                        val._type = CalculatorOperation::STRING;
                        val._string =part.mid(1,part.size() - 2);
                    }
                    else if ( part[0] == '@'){
                        pindex = part.mid(1,1).toInt(&ok);
                        fillValues(pindex, part, val, action._action);
                    }else {
                        PIXVALUETYPE number = part.toDouble(&ok);
                        if (ok){
                            val._type = CalculatorOperation::NUMERIC;
                            val._value = number;
                        }else if (part.indexOf("LINK:") == 0){
                            int link = part.mid(5).toInt();
                            val._type = CalculatorOperation::LINK;
                            val._link = link;
                        }else if ( part == sUNDEF){
                            val._type = CalculatorOperation::NUMERIC;
                            val._value = PIXVALUEUNDEF;
                        }else {
                            kernel()->issues()->log(TR("Error in expression. Index value is not valid :") + part);
                            return IDomain();
                        }
                    }
                    action._values.push_back(val);
                }
            }
            std::reverse(action._values.begin(), action._values.end());
            _actions.push_back(action);
        }
        return outputDomain;
    } catch (ErrorObject& err){
        return IDomain();
    }
}

PIXVALUETYPE CalculatorOperation::calc(const std::vector<Action>& localActions) {
    auto GetValue = [&](const ParmValue& parm,const std::vector<PIXVALUETYPE>& result, bool *isNumeric=0)->PIXVALUETYPE{
        switch(parm._type){
        case ParmType::LINK:
            return result[parm._link];break;
        case ParmType::ITERATOR:
            return *(*(parm._source));break;
        case ParmType::NUMERIC:{
            if ( isNumeric ) *isNumeric=true;
            return parm._value;break;
        }
        case ParmType::DOMAINITEM:
            return parm._value;break;
        case ParmType::COLUMN:
            return parm._columnValues[_record].toDouble();break;
        }
        return PIXVALUEUNDEF;
    };

    auto CalcBinary = [](MathAction act, PIXVALUETYPE v1, PIXVALUETYPE v2) ->PIXVALUETYPE{
        if (isNumericalUndef(v1)  || isNumericalUndef(v2))
            return PIXVALUEUNDEF;
        switch(act){
        case maADD:
            return v1+v2;
        case maDIVIDE:
            return v2 == 0 ? PIXVALUEUNDEF : v1/v2;
        case maMULT:
            return v1 * v2;
        case maMINUS:
            return v1 - v2;
        case maPOW:
            return std::pow(v1,v2);
        case maMAX:
            return std::max(v1,v2);
        case maMIN:
            return std::min(v1,v2);


        default:
            return PIXVALUEUNDEF;
        }
        return PIXVALUEUNDEF;
    };

    PIXVALUETYPE calcResult;
    std::vector<PIXVALUETYPE> result(localActions.size(),PIXVALUEUNDEF);
    for(int i=0; i < localActions.size(); ++i){
        const Action& action = localActions[i];
        switch(action._action){
            case maADD:
            case maMULT:
            case maDIVIDE:
            case maMAX:
            case maMIN:
            case maMINUS:
            case maPOW:
            {
                PIXVALUETYPE v1 = GetValue(action._values[0],result);
                PIXVALUETYPE v2 = GetValue(action._values[1],result);
                calcResult = CalcBinary(action._action,v1,v2);
                break;
            }
            case maSIN:
            {
                PIXVALUETYPE v = GetValue(action._values[0],result);
                calcResult =  isNumericalUndef(v) ? PIXVALUEUNDEF : std::sin(v);
                break;
            }
            case maCOS:
            {
                PIXVALUETYPE v = GetValue(action._values[0],result);
                calcResult =  isNumericalUndef(v) ? PIXVALUEUNDEF : std::cos(v);
                break;
            }
            case maACOSH:
            {
                PIXVALUETYPE v = GetValue(action._values[0],result);
                calcResult =  isNumericalUndef(v) ? PIXVALUEUNDEF : std::acosh(v);
                break;
            }
            case maASINH:
            {
                PIXVALUETYPE v = GetValue(action._values[0],result);
                calcResult =  isNumericalUndef(v) ? PIXVALUEUNDEF : std::asinh(v);
                break;
            }
            case maTAN:
            {
                PIXVALUETYPE v = GetValue(action._values[0],result);
                calcResult =  (std::abs(v) == M_PI / 2 || isNumericalUndef(v)) ? PIXVALUEUNDEF : std::tan(v);
                break;
            }
            case maTANH:
            {
                PIXVALUETYPE v = GetValue(action._values[0],result);
                calcResult =  (isNumericalUndef(v)) ? PIXVALUEUNDEF : std::tanh(v);
                break;
            }
            case maACOS:
            {
                PIXVALUETYPE v = GetValue(action._values[0],result);
                calcResult =  ( v < -1 || v > 1 || isNumericalUndef(v)) ? PIXVALUEUNDEF : std::acos(v);
                break;
            }
            case maCOSH:
            {
                PIXVALUETYPE v = GetValue(action._values[0],result);
                calcResult =  (isNumericalUndef(v)) ? PIXVALUEUNDEF : std::cosh(v);
                break;
            }
            case maSINH:
            {
                PIXVALUETYPE v = GetValue(action._values[0],result);
                calcResult =  ( isNumericalUndef(v)) ? PIXVALUEUNDEF : std::sinh(v);
                break;
            }
            case maASIN:
            {
                PIXVALUETYPE v = GetValue(action._values[0],result);
                calcResult =  ( v < -1 || v > 1 || isNumericalUndef(v)) ? PIXVALUEUNDEF : std::asin(v);
                break;
            }
            case maATAN:
            {
                PIXVALUETYPE v = GetValue(action._values[0],result);
                calcResult =  ( v < -M_PI/2 || v > M_PI/2 || isNumericalUndef(v)) ? PIXVALUEUNDEF : std::acos(v);
                break;
            }
            case maATANH:
            {
                PIXVALUETYPE v = GetValue(action._values[0],result);
                calcResult =  (isNumericalUndef(v)) ? PIXVALUEUNDEF : std::atanh(v);
                break;
            }
            case maLOG10:
            {
                PIXVALUETYPE v = GetValue(action._values[0],result);
                calcResult =  ( v <= 0 || isNumericalUndef(v)) ? PIXVALUEUNDEF : std::log10(v);
                break;
            }
            case maLOG:
            {
                PIXVALUETYPE v1 = GetValue(action._values[0],result);
                PIXVALUETYPE v2 = GetValue(action._values[1],result);
                if ( v2 > 0){
                    v2 = log(v2);
                    calcResult =  ( v1 <= 0 || isNumericalUndef(v1) || v2 == 0) ? PIXVALUEUNDEF : std::log2(v1) / v2;
                } else
                    calcResult = PIXVALUEUNDEF;
                break;
            }
            case maLN:
            {
                PIXVALUETYPE v = GetValue(action._values[0],result);
                calcResult =  ( v <=0 ||  isNumericalUndef(v)) ? PIXVALUEUNDEF : std::log(v);
                break;
            }
            case maEXP:
            {
                PIXVALUETYPE v = GetValue(action._values[0],result);
                calcResult =  (isNumericalUndef(v)) ? PIXVALUEUNDEF : std::exp(v);
                break;
            }
            case maINT:
            {
                PIXVALUETYPE v = GetValue(action._values[0],result);
                calcResult =  (isNumericalUndef(v)) ? PIXVALUEUNDEF : (qlonglong)(v);
                break;
            }
            case maROUND:
            {
                PIXVALUETYPE v = GetValue(action._values[0],result);
                calcResult =  (isNumericalUndef(v)) ? PIXVALUEUNDEF : std::round(v);
                break;
            }
            case maABS:
            {
                PIXVALUETYPE v = GetValue(action._values[0],result);
                calcResult =  (isNumericalUndef(v)) ? PIXVALUEUNDEF : std::abs(v);
                break;
            }
            case maSQ:
            {
                PIXVALUETYPE v = GetValue(action._values[0],result);
                calcResult =  (isNumericalUndef(v)) ? PIXVALUEUNDEF : v * v;
                break;
            }
            case maSQRT:
            {
                PIXVALUETYPE v = GetValue(action._values[0],result);
                calcResult =  (v < 0 || isNumericalUndef(v)) ? PIXVALUEUNDEF : std::sqrt(v);
                break;
            }
            case maFLOOR:
            {
                PIXVALUETYPE v = GetValue(action._values[0],result);
                calcResult =  isNumericalUndef(v) ? PIXVALUEUNDEF : std::floor(v);
                break;
            }
            case maCEIL:
            {
                PIXVALUETYPE v = GetValue(action._values[0],result);
                calcResult =  isNumericalUndef(v) ? PIXVALUEUNDEF : std::ceil(v);
                break;
            }
            case maEQ:
            {
                bool isNumeric=false;
                PIXVALUETYPE v1 = GetValue(action._values[0],result,&isNumeric);
                PIXVALUETYPE v2 = GetValue(action._values[1],result,&isNumeric);
                calcResult = ( !isNumeric && (isNumericalUndef(v1) || isNumericalUndef(v2))) ? PIXVALUEUNDEF : v1 == v2;
                break;
            }
            case maNEQ:
            {
                bool isNumeric=false;
                PIXVALUETYPE v1 = GetValue(action._values[0],result,&isNumeric);
                PIXVALUETYPE v2 = GetValue(action._values[1],result,&isNumeric);
                calcResult = ( !isNumeric && (isNumericalUndef(v1) || isNumericalUndef(v2))) ? PIXVALUEUNDEF : v1 != v2;
                break;
            }
            case maLESSEQ:
            {
                PIXVALUETYPE v1 = GetValue(action._values[0],result);
                PIXVALUETYPE v2 = GetValue(action._values[1],result);
                calcResult =  isNumericalUndef(v1) || isNumericalUndef(v2) ? PIXVALUEUNDEF : ( v1 <= v2);
                break;
            }
            case maLESS:
            {
                PIXVALUETYPE v1 = GetValue(action._values[0],result);
                PIXVALUETYPE v2 = GetValue(action._values[1],result);
                calcResult = isNumericalUndef(v1) || isNumericalUndef(v2) ? PIXVALUEUNDEF : ( v1 < v2);
                break;
            }
            case maGREATEREQ:
            {
                PIXVALUETYPE v1 = GetValue(action._values[0],result);
                PIXVALUETYPE v2 = GetValue(action._values[1],result);
                calcResult =  isNumericalUndef(v1) ||isNumericalUndef(v2) ? PIXVALUEUNDEF : ( v1 >= v2);
                break;
            }
            case maGREATER:
            {
                PIXVALUETYPE v1 = GetValue(action._values[0],result);
                PIXVALUETYPE v2 = GetValue(action._values[1],result);
                calcResult =  isNumericalUndef(v1) || isNumericalUndef(v2) ? PIXVALUEUNDEF : ( v1 > v2);
                break;
            }
            case maAND:
            {
                PIXVALUETYPE v1 = GetValue(action._values[0],result);
                PIXVALUETYPE v2 = GetValue(action._values[1],result);
                if (!(bool)v1 || !(bool)v2)
                    calcResult = false;
                else if (isNumericalUndef(v1) || isNumericalUndef(v2))
                    calcResult = PIXVALUEUNDEF;
                else
                    calcResult = true;
                break;
            }
            case maOR:
            {
                PIXVALUETYPE v1 = GetValue(action._values[0],result);
                PIXVALUETYPE v2 = GetValue(action._values[1],result);
                if (isNumericalUndef(v1)) {
                    if (isNumericalUndef(v2))
                        calcResult = PIXVALUEUNDEF;
                    else
                        calcResult = (bool)v2 ? true : PIXVALUEUNDEF;
                } else if (isNumericalUndef(v2))
                    calcResult = (bool)v1 ? true : PIXVALUEUNDEF;
                else
                    calcResult = (bool)v1 || (bool)v2;
                break;
            }
            case maIFF:
            {
                // if the action is a iterator we can directly get its value from the iterator else if it is a
                // comparisson it will be calculated previously and its result will be in calcresult
                PIXVALUETYPE test = GetValue(action._values[0],result);
                PIXVALUETYPE v2 = GetValue(action._values[1],result);
                PIXVALUETYPE v3 = GetValue(action._values[2],result);
                calcResult = test == PIXVALUEUNDEF ? PIXVALUEUNDEF : ((bool)test ? v2 : v3);
                break;
            }
            case maATTRIBUTE:
            {
                PIXVALUETYPE v1 = GetValue(action._values[0], result);
                calcResult = action._values[0]._keyMapping.at((quint32)v1);

            }
            case maNOT:
            {
                PIXVALUETYPE v = GetValue(action._values[0], result);
                calcResult = isNumericalUndef(v) ? PIXVALUEUNDEF : ~(qint64)v;
                break;

            }
            case maXOR:
            {
                PIXVALUETYPE v1 = GetValue(action._values[0], result);
                PIXVALUETYPE v2 = GetValue(action._values[1], result);
                calcResult = (isNumericalUndef(v1) || isNumericalUndef(v2))? PIXVALUEUNDEF : ((qint64)v1) ^ ((qint64)v2);
                break;

            }
            case maIFUNDEF:
            {
                PIXVALUETYPE v1 = GetValue(action._values[0], result);
                PIXVALUETYPE v2 = GetValue(action._values[1], result);
                PIXVALUETYPE v3 = GetValue(action._values[2], result);
                calcResult = isNumericalUndef(v1) ? v2 : v3;
                break;
            }
            case maIFNOTUNDEF:
            {
                PIXVALUETYPE v1 = GetValue(action._values[0], result);
                PIXVALUETYPE v2 = GetValue(action._values[1], result);
                PIXVALUETYPE v3 = GetValue(action._values[2], result);
                calcResult = (!isNumericalUndef(v1)) ? v2 : v3;
                break;
            }
        }
        result[i] = calcResult;
    }
    return result.back();
}
