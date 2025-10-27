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
#include <iostream>
#include "kernel.h"
#include "symboltable.h"
#include "commandhandler.h"
#include "ilwiscontext.h"
#include "catalog.h"
#include "raster.h"
#include "pixeliterator.h"
#include "errorobject.h"
#include "rastercoverage.h"

int main(int argc, char *argv[])
{
    try{
        QCoreApplication a(argc, argv);
        double sum = 0, sum2 = 0;
        Ilwis::initIlwis(Ilwis::rmCOMMANDLINE | Ilwis::rmNOUI);
        Ilwis::ExecutionContext ctx;
        Ilwis::SymbolTable syms;
        Ilwis::IRasterCoverage raster;
        Ilwis::context()->setWorkingCatalog(Ilwis::ICatalog("file:///home/mschouwen/data/ilwisdata"));
        QString expr = QString("aaas3{format(ilwis3,map)}=aggregaterasterstatistics(%1,max)").arg("file:///home/mschouwen/temp/Sentinel2_CA_2018-11-18.dat");
        //QString expr = QString("aaas4{format(ilwis3,map)}=mirrorrotateraster(%1,mirrdiag)").arg("file:///home/mschouwen/data/ilwisdata/small9_6.mpr");
        //QString expr = QString("aaas4{format(ilwis3,map)}=selection(%1,boundingbox(2 3, 5 6))").arg("file:///home/mschouwen/data/ilwisdata/small9_6.mpr");
        expr = "script " + expr;
        Ilwis::kernel()->startClock();
        Ilwis::commandhandler()->execute(expr, &ctx, syms) ;
        //Ilwis::kernel()->endClock();
        //raster.prepare("file:///home/mschouwen/data/ilwisdata/netcdf/urkdn2019.nc/NDVI");
        auto v = raster->pix2value(Ilwis::Pixel(3,3,2));


        Ilwis::exitIlwis();
        return 0;

    } catch(const Ilwis::ErrorObject& err){
        std::cerr << err.message().toStdString();

    }
    catch(std::exception& ex){
        std::cerr << ex.what();
    }
    catch(...){
        std::cerr << "unknown error";
    }
    return 1;
}
