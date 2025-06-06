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

#ifndef COORDINATESYSTEM_H
#define COORDINATESYSTEM_H

#include "kernel_global.h"

namespace Ilwis {

class CoordinateSystem;
template<class T> class Box;

typedef IlwisData<CoordinateSystem> ICoordinateSystem;

class KERNELSHARED_EXPORT CoordinateSystem : public IlwisObject
{
public:
    CoordinateSystem();
    CoordinateSystem(const Ilwis::Resource &resource);

    virtual Coordinate coord2coord(const ICoordinateSystem& sourceCs, const Coordinate& crdSource) const =0;
    virtual LatLon coord2latlon(const Coordinate &crdSource) const =0;
    virtual Coordinate latlon2coord(const LatLon& ll) const = 0;
    virtual Ilwis::Envelope convertEnvelope(const ICoordinateSystem& sourceCs, const Envelope& envelope) const;
    virtual bool canConvertToLatLon() const;
    virtual bool canConvertToCoordinate() const;
    virtual Coordinate inverseCoordinateConversion(const CoordinateSystem& cs, const Coordinate& crd) const;
    Ilwis::Envelope envelope(bool tolatlon=false) const;
    void envelope(const Envelope &env);
    virtual bool isLatLon() const = 0;
    virtual bool isUnknown() const = 0;
    virtual QString toWKT(quint32 spaces=0) const=0;
    virtual QString toEpsg() const { return sUNDEF;}
    static ICoordinateSystem fromWKT(const QString& wkt);

	static bool addCsyProperty(const ICoordinateSystem& csy, Resource& resource);
	static Envelope latLonEnvelope(const ICoordinateSystem& cs, const Envelope& env);



protected:
    void copyTo(IlwisObject *obj);
private:
    Ilwis::Envelope _envelope;

};

}

Q_DECLARE_METATYPE(Ilwis::ICoordinateSystem);


#endif // COORDINATESYSTEM_H
