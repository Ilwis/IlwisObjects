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

#include <QUrl>
//#include "module.h"
#include "kernel.h"
#include "ilwisdata.h"
#include "geometries.h"
#include "domain.h"
#include "coordinatesystem.h"
#include "coverage.h"
#include "coordinatedomain.h"

using namespace Ilwis;

CoordinateDomain::CoordinateDomain()
{
}

CoordinateDomain::CoordinateDomain(const Resource &resource) : Domain(resource)
{
    if ( hasType(resource.ilwisType(),itCOVERAGE)) {
        ICoverage cov;
        if( cov.prepare(resource)) {
            _csy = cov->coordinateSystem();
            _envelope.reset(cov->envelope().clone());
        }
    } else if ( hasType(resource.ilwisType(),itCOORDSYSTEM)) {
        _csy.prepare(resource);
    }
}

void CoordinateDomain::setRange(Range *vr)
{
    _envelope = QSharedPointer<Range>(vr);
}

Domain::Containement CoordinateDomain::contains(const QVariant &value) const
{
    Coordinate crd = value.value<Coordinate>();
    if ( _envelope->valueType() == itCOORDINATE) {
        QSharedPointer<Envelope> box = _envelope.dynamicCast<Envelope>();
        if ( box)
            return box->contains(crd) ? Domain::cSELF : Domain::cNONE;
    }
    return Domain::cNONE;
}

bool CoordinateDomain::isCompatibleWith(const IlwisObject *dom, bool strict) const {
    if ( !dom || !dom->isValid())
        return false;
    if(dom->ilwisType() != itCOORDDOMAIN)
        return false;
    ICoordinateDomain crddom;
    crddom.prepare(dom->id());
    if ( crddom.isValid())
        return  coordinateSystem() == crddom->coordinateSystem();
    return false;
}


IlwisTypes CoordinateDomain::ilwisType() const
{
    return itCOORDDOMAIN;
}

ICoordinateSystem CoordinateDomain::coordinateSystem() const
{
    return _csy;
}

void CoordinateDomain::setCoordinateSystem(ICoordinateSystem& csy)
{
    _csy = csy;
}

void CoordinateDomain::range(Range *)
{

}

SPRange CoordinateDomain::getRange() const
{
    return _envelope;
}

bool CoordinateDomain::isOrdered() const
{
    return false;
}

IlwisTypes CoordinateDomain::valueType() const
{
    return itCOORDINATE;
}


QVariant CoordinateDomain::impliedValue(const QVariant &v) const
{
    Coordinate crd = v.value<Coordinate>();
    QString c = QString("%1 %2 %3").arg(crd.x, 0, 'g').arg(crd.y, 0, 'g').arg(crd.z, 0, 'g');
    return c;
}
