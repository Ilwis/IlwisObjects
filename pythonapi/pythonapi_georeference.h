#ifndef PYTHONAPI_GEOREFERENCE_H
#define PYTHONAPI_GEOREFERENCE_H

#include "pythonapi_ilwisobject.h"
#include "pythonapi_coordinatesystem.h"
#include "pythonapi_util.h"

namespace Ilwis{
    class GeoReference;
    typedef IlwisData<GeoReference> IGeoReference;
}

namespace pythonapi {

    class GeoReference : public IlwisObject{
        friend class Engine;
        friend class RasterCoverage;
        friend class Catalog;
        private:
            GeoReference(const Ilwis::IGeoReference& gr);
        protected:
            virtual const QString getStoreFormat() const;
        public:
            GeoReference(const std::string &resource);
            GeoReference(const std::string& csyCode, const Envelope& env, const Size& psize, const std::string cornerstype = "cornerofcorners"); // corners grf
            GeoReference(const CoordinateSystem& csy,const Envelope &env,const Size &psize, const std::string cornerstype = "cornerofcorners");
            static GeoReference* toGeoReference(Object *obj);

            CoordinateSystem coordinateSystem() const;
            void setCoordinateSystem(const CoordinateSystem& csy);

            Coordinate pixel2Coord(const PixelD& pixel) const;
            Coordinate pixel2Coord(const Pixel& pixel) const;
            PixelD coord2Pixel(const Coordinate& coord) const;
            virtual Envelope box2Envelope(const Box &box ) const;
            virtual Box envelope2Box(const Envelope &box) const;
            double pixelSize() const;

            Size size() const;
            void setSize(const Size& sz);
            bool centerOfPixel() const;
            void setCenterOfPixel(bool yesno);
            bool isCompatible(const GeoReference& other) const;
            bool compute();
            IlwisTypes ilwisType() const;
            const std::string toString() const;
    };

} // namespace pythonapi

#endif // PYTHONAPI_GEOREFERENCE_H
