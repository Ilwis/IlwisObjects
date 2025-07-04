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

#ifndef PixelIterator_H
#define PixelIterator_H

#include <iterator>

namespace geos{
namespace geom {
class Geometery;
}
}

namespace Ilwis {

typedef std::unique_ptr<geos::geom::Geometry> UPGeometry;
class ProcessingBoundingBoxes;

class Tranquilizer;
typedef std::shared_ptr<Tranquilizer> SPTranquilizer;

/*!
 * \brief The PixelIterator class an iterator class that iteratos over all the pixels in an grid (or subsection of it)
 *
 * The pixeliterator is the main access mechanism (together with the blockiterator) to pixels in 2D or 3D gridcoverages. Basically it sees the pixels as one long linear space and moves over it. The movement (flow) can have several directions
 * which resembles directions in the 'real' world. In the default case it moves first over the x dirdction, than the y and finally the z direction. But one could just as well first move in the z direction than x, than y. Internally these movements
 * are translated to offsets in the linear space. The pixel iterator obeys the normal rules for iterators in the STL and thus can be combined with the algorithms in this library.
 *
 * Here are basically two flavors of the  pixeliterator : the pixel iterator itself and the blockiterator. The pixeliterator moves over single pixels, the blockiterator  moves over an [n x m x p] size block (some algorithms need this).
 * So how does this work? Suppose you want to count all the pixels in a raster-coverage with a value greater than 100. First the crude way (there are other more elegant ways to do this)
 *
 *    1 PixelIterator iterator(someraster);
 *    2 PixelIterator endlocation = end(someraster);
 *    3 int count = 0;
 *    4 while ( iterator != endlocation){
 *    5     if (*iterator > 100)
 *    6           ++count;
 *    7 }
 *
 * In line 1 a new pixeliterator is created. The simplest constructor only needs the raster-coverage as parameter. In line 2 we create the end iterator.  By definition (in c++) the end iterator is located at one container cell beyond the last valid location within the container. With the pixel iterator this means +1 layer beyond the last pixel in the last layer. The end method in line 2 is a convenience method that generates an iterator that satisfies this.
 * Line 4 is end-condition for the traversing the raster-coverage. In line 5 the data of a pixel is accessed by using the ‘*’ operator on the iterator. This operator returns a reference to a PIXVALUETYPE (64-bit) numerical value that is the value of the pixels.
 *
 * A more modern way (using lambda’s)
 *
 *    1 int count = 0;
 *    2 std::for_each (begin(someraster), end(someraster) , [&](PIXVALUETYPE& v)){
 *    3     if v > 100) ++count;
 *    4 });
 *
 * It becomes more interesting  when we are combining iterators. Suppose in the example above we want to create a new raster-coverage that only contains the pixels with a numerical value greater than 100
 *
 *    1 PixelIterator iterInput(someraster);
 *    2 std::for_each(begin(outputRaster), end(outputRaster) , [&](PIXVALUETYPE& v)){
 *    3     if v > 100) v = *iterInput;
 *    4    ++iterInput;
 *    5 });
 *
 * In line 1 we define a iterator for the input, the output operator(invisible) is hidden behind convenience begin() and end() methods. Remember that the ‘*’operator gives a reference to the value of (in this case output) raster-coverage. ‘for_each’ uses this operator to give access to values of the container and exposes this through the PIXVALUETYPE& (so also a reference). We only need to set this value to actually change values in the output raster-coverage. That’s all. Assuming that input and output have the same geometry.
 * The other interesting part is that the pixeliterator integrates with the existing STL library of C++. Giving access to a large number of existing basic algorithms. Suppose You need to copy all of the values of a raster-coverage to a vector(array) and do some operation on it. Due to memory limitations this might not always be a good idea, but there are certainly enough use cases were this is useful.
 *
 *    std::vector<PIXVALUETYPE> data(someraster.size().totalSize());
 *    std::copy(begin(someraster), end(someraster), data.begin());
 *
 * or swapping of the pixels of two raster-coverages
 *
 *    swap_ranges(begin(someraster),end(someraster), begin(otherraster));
 *
 * Note that because the iterator ‘automatically’ moves in layer-index direction all algorithms also work on stacks of raster layers.
 *
 */
class KERNELSHARED_EXPORT PixelIterator   {
public:
	typedef std::random_access_iterator_tag iterator_category;
    typedef PIXVALUETYPE value_type;
    typedef ptrdiff_t difference_type;
    typedef PIXVALUETYPE* pointer;
    typedef PIXVALUETYPE& reference;
    /*!
     * The possible flows, not all are implemented (yet).<br>
     * atm only xyz works
     */
    enum Flow { fXYZ, fZXY};

    /*!
     * \brief isValid tells if an iterator is in a valid state.
     *
     * An iterator is valid if bounding box and grid are a valid combination.
     * \return true if the bounding box fits in the grid.
     */
    bool isValid() const{
        return _isValid;
    }


    /*!
     * \brief The empty constructor for PixelIterator
     *
     * This constructor creates an empty PixelIterator
     */
    PixelIterator();
    PixelIterator(const IRasterCoverage& raster, geos::geom::Geometry *selection);
    PixelIterator(const IRasterCoverage& raster, const UPGeometry& selection) : PixelIterator(raster, selection.get()){

    }

    /*!::value
     * \brief Constructs a PixelIterator from a raster and a bounding box
     *
     * Constructs a PixelIterator from the raster and the bounding box<br>
     * The bounding box is the area within the raster which should be walked by the PixelIterator<br>
     * If one of the parameters was invalid the PixelIterator will also be invalid
     *
     * \param raster The raster from which this PixelIterator should be created
     * \param box The bounding box which desides what part of the raster should be walked
     */
    PixelIterator(const IRasterCoverage& raster, const BoundingBox& box=BoundingBox(), Flow flow=fXYZ);
	PixelIterator(const IRasterCoverage& raster, int threadIdx, const BoundingBox& box = BoundingBox(), Flow flow = fXYZ);
    PixelIterator(const IRasterCoverage& raster, int threadIdx, const ProcessingBoundingBoxes& boxes, Flow flow = fXYZ);
    PixelIterator(const IRasterCoverage& raster, Flow flow);
	
    /*!
     * Copy's all the values from the existing PixelIterator onto this one
     *
     * \param iter The iterator whose values should be duplicated
     */
    PixelIterator(const PixelIterator& iter);

    /*!
     * \brief Copy constructor
     * \param iter PixelIterator that must be copied
     */
    PixelIterator(PixelIterator &&iter);
    ~PixelIterator();

    /*!
     * override of the operator=<br>
     * copies the values of the supplied iterator onto this one<br>
     *  so: this=iter;<br>
     * results in this getting the same values as iter
     *
     * \param iter pixeliterator that must be copied
     * \return this, with the modified values
     */
    PixelIterator& operator=(const PixelIterator& iter);
    PixelIterator& operator=(const Pixel &pix){
			qint64 dIndex = 0;
			auto xz = pix.z == iUNDEF ? 0 : pix.z;
			switch (_flow) {
			case fXYZ:
				dIndex = (pix.x - _x) + (pix.y - _y) * _box.xlength() + (xz - _z) * _box.xlength() * _box.ylength();
				move(dIndex);
				break;
			case fZXY:
				dIndex = (pix.z - _z) + (pix.x - _x) * _box.ylength() + (pix.y - _y) * _box.xlength() * _box.zlength();
				move(dIndex);
				break;
			default:
				_x = pix.x;
				_y = pix.y;
				_z = pix.z;
				_yChanged = _xChanged = _zChanged = true;
				initPosition();
			}
			return *this;

			/*_x = pix.x;
			_y = pix.y;
			_z = pix.z == iUNDEF ? 0 : pix.z;
			_yChanged = _xChanged = _zChanged = true;
			initPosition();

		return *this;*/
    }

    /*!
     * override of the operator=<br>
     * copies the values of the supplied iterator onto this one<br>
     *  so: this=iter;<br>
     * results in this getting the same values as iter
     *
     * \param iter pixeliterator that must be copied
     * \return this, with the modified values
     */
    PixelIterator& operator=(const PixelIterator&& iter);

    /*!
     * \brief Adds 1 to this PixelIterator
     * the direction moved in the raster depends on the boundingbox in combination with the flow
     * \return This iterator moved 1
     */
    PixelIterator& operator++() {
        move(_step);
        return *this;
    }

    /*!
     * \brief Substracts 1 from this PixelIterator
     * the direction moved in the raster depends on the boundingbox in combination with the flow
     * \return This iterator moved -1
     */
    PixelIterator& operator--() {
        move(-_step);
        return *this;
    }

    /*!
     * \brief Adds n to this PixelIterator
     * the direction moved in the raster depends on the boundingbox in combination with the flow
     * \param n amount to add to this PixelIterator
     * \return this iterator moved n
     */
    PixelIterator& operator+=(qint64 n) {
        move(n);
        return *this;
    }

    /*!
     * \brief Substracts n from this PixelIterator
     * the direction moved in the raster depends on the boundingbox in combination with the flow
     * \param n amount to substract
     * \return this iterator moved -n
     */
    PixelIterator& operator-=(qint64 n) {
        move(-n);
        return *this;
    }

    /*!
     * \brief Random acces operator
     *
     * Moves to the specified index and returns a pointer to the value at this index
     *
     * \param index the target index
     * \return the value at the index
     */
    PIXVALUETYPE& operator[](qint64 index){
        _x = 0;
        _y = 0;
        _z = 0;
        _yChanged = _xChanged = _zChanged = true;
        initPosition();
        move(index);
        return this->operator *();
    }

    QVariant operator()(const QString& column);

    /*!
     * \brief Moves the PixelIterator to the given Pixel
     *
     * Moves the PixelIterator to the given Pixel, also adjusts the lineairposition (index) of this iterator<br>
     * the voxel should have all three coordinates defined.<br>
     * invalid voxel will cause this iterator to become invalid.
     *
     * \param vox the Pixel to move to
     * \return this at the given Pixel
     */
    PixelIterator &operator ()(const Pixel &pix)
    {

			_x = pix.x;
			_y = pix.y;
			_z = pix.z;
			_yChanged = _xChanged = _zChanged = true;
			initPosition();
        return *this;
    }

    /*!
     * \brief Adds 1 to this PixelIterator
     * the direction moved in the raster depends on the boundingbox in combination with the flow
     * \return this iterator moved 1
     */
    PixelIterator operator++(int);

    /*!
     * \brief Substracts 1 from this PixelIterator
     * the direction moved in the raster depends on the boundingbox in combination with the flow
     * \return this iterator moved -1
     */
    PixelIterator operator--(int);

    /*!
     * \brief Random acces operator
     *
     * Moves to the specified voxel, requires a valid voxel
     *
     * \param index the target index
     * \return this iterator at the specified voxel
     */
    PixelIterator& operator[](const Pixel& pix){
        return operator ()(pix);
    }

    /*!
     * \brief Checks if this PixelIterator has the same linearPosition as another PixelIterator
     * \param iter another PixelIterator
     * \return true when linearPosition is equal
     */
    bool operator==(const PixelIterator& iter) const;

    /*!
     * \brief Checks if this PixelIterator has a different linearPosition as another PixelIterator
     * \param iter another PixelIterator
     * \return true when linearPosition is not equal
     */
    bool operator!=(const PixelIterator& iter) const;

    /*!
     * \brief Checks if this PixelIterator has a smaller linearPosition as another PixelIterator
     * \param iter another PixelIterator
     * \return true when linearPosition(this)<linearPosition(iter)
     */
    bool operator<(const PixelIterator& iter) const;

    /*!
     * \brief Checks if this PixelIterator has a smaller or equal linearPosition as another PixelIterator
     * \param iter another PixelIterator
     * \return true when linearPosition(this)<=linearPosition(iter)
     */
    bool operator<=(const PixelIterator& iter) const;

    /*!
     * \brief Checks if this PixelIterator has a bigger linearPosition as another PixelIterator
     * \param iter another PixelIterator
     * \return true when linearPosition(this)>linearPosition(iter)
     */
    bool operator>(const PixelIterator& iter) const;

    /*!
     * \brief Checks if this PixelIterator has a bigger or equal linearPosition as another PixelIterator
     * \param iter another PixelIterator
     * \return true when linearPosition(this)>=linearPosition(iter)
     */
    bool operator>=(const PixelIterator& iter) const;

    /*!
     * \brief Query for a reference to the current value of the PixelIterator
     * \return reference to the currentvalue
     */
    PIXVALUETYPE& operator*() {
        return _grid->value(_z, _y, _x, _threadIndex);
    }

    /*!
     * \brief Query for a reference to the current vallue of the PixelIterator
     * \return reference to the currentvalue
     */
    const PIXVALUETYPE& operator*() const {
        return  _grid->value(_z, _y, _x, _threadIndex);
    }

    /*!
     * \brief Query for the current value of the PixelIterator
     * \return ->value(this(current))
     */
    PIXVALUETYPE* operator->() {
        return &(_grid->value(_z, _y, _x, _threadIndex));
    }

    /*!
     * \brief Returns the end position of this PixelIterator, this is 1 past the actual lastblock of the boundingbox
     * \return the endvalue of the lineairposition
     */
    PixelIterator end() const ;
    void toEnd();

    /*!
     * \brief Changes the flow of this PixelIterator
     *
     * Changes the flow according to the submitted value for flw,<br>
     * if flw =fXYZ<br>
     * this means that the pixeliterator will start at the first row in x direction,<br>
     * at the end of this row, it will go 1 in y direction and start over in x direction,<br>
     * at the end of the y row, it will go 1 in z direction and start over<br>
     *
     * note at this moment only the xyz flow is implemented
     * \param flw the requested flow
     */
    void setFlow(Flow flw);

    /*!
     * \brief Checks if a certain pixel is inside this PixelIterator
     *
     * \param pix the pixel to be checked
     * \return true when the pixel is in this PixelIterator
     */
    bool contains(const Pixel& pix) ;
	qint64 step() const { return _step; }
	void step(qint64 n) { if (n != isNumericalUndef(n)) _step = n; }

    /*!
     * \brief Checks if the x coordinate has changed in the last step taken
     * all movement and/or position change commands automatically set the x y and z flags when appropiate
     * \return true if the x has changed
     */
    bool xchanged() const;

    /*!
     * \brief Checks if the y coordinate has changed in the last step taken
     * all movement and/or position change commands automatically set the x y and z flags when appropiate
     * \return true if the y has changed
     */
    bool ychanged() const;

    /*!
     * \brief Checks if the z coordinate has changed in the last step taken
     * all movement and/or position change commands automatically set the x y and z flags when appropiate
     * \return true if the z has changed
     */
    bool zchanged() const;
    quint64 x() const { return _x;}
    quint64 y() const { return _y;}
    quint64 z() const { return _z;}

    /*!
     * \brief Checks if this PixelIterator is at its endpoint
     * \return true if it is at the end
     */
    bool isAtEnd() const {
        return _x == _endx &&
               _y == _endy &&
               _z == _endz;
    }

    /*!
     * \brief Query for the current non lineair position of this PixelIterator
     * \return a Pixel with the current non lineair position of this PixelIterator
     */
    Pixel position() const;
    void position(const Pixel& pos) ;

    /*!
     * \brief Query for the bounding box of this PixelIterator
     * the bounding box decides which part of the raster this PixelIterator should walk, thus its size can never be bigger than the rastersize
     * \return the bounding box of this PixelIterator
     */
    const BoundingBox& box() const;
	void box(const BoundingBox& box);

    /*!
     * \brief Query for the linearPosition of this PixelIterator
     * \return the lineairposition
     */
    quint64 linearPosition() const;

    /*!
     * \brief Sets the raster of this PixelIterator
     *
     * Sets the raster of this PixelIterator
     *
     * \param raster the given raster
     */
    void setRaster(const IRasterCoverage &raster);
	void threadIndex(int idx) {
		_threadIndex = idx >= 0 ? idx : 0;
	}

    /*!
     * \brief Substracts another pixeliterator fomr this PixelIterator
     * the direction moved in the raster depends on the boundingbox in combination with the flow<br>
     * this might have unexpected results if both iterators have a different flow.
     * \param iter2 another PixelIterator
     * \return this iterator moved -iter2.linearPosition()
     */
    qint64 operator -(const PixelIterator &iter2);

    /*!
     * \brief Adds n to this PixelIterator and returns it
     * the direction moved in the raster depends on the boundingbox in combination with the flow
     * \param n the amount to be added
     * \return this iterator moved n
     */
    PixelIterator operator+(qint64 n) {
        PixelIterator iter(*this);
        iter.move(n);
        return iter;
    }

    const Ilwis::IRasterCoverage &raster() const;

protected:
    PixelIterator(quint64 endpos ) :
        _grid(0),
        _x(0),
        _y(0),
        _z(0),
        _flow(fXYZ),
        _isValid(true),
        _endx(0),
        _endy(0),
        _endz(0),
        _linearposition(endpos),
        _endposition(endpos),
        _xChanged(false),
        _yChanged(false),
        _zChanged(false)
    {
    }

    void init();
    void initPosition();
    void copy(const PixelIterator& iter);

    IRasterCoverage _raster;
    Grid *_grid = 0;
    BoundingBox _box;
    qint64 _x = 0;
    qint64 _y = 0;
    qint64 _z = 0;
    Flow _flow;
    bool _isValid;
    qint64 _endx;
    qint64 _endy;
    qint64 _endz;
    qint64 _linearposition;
    qint64 _endposition;
    bool _xChanged =false;
    bool _yChanged = false;
    bool _zChanged = false;
    std::vector<std::vector<qint64>> _selectionPixels;
    qint64 _selectionIndex = -1;
    bool _insideSelection = false;
    qint64 _step = 1;
	int _threadIndex = 0;


    bool move(qint64 n) {
        if ( n == 0){
            _xChanged = _yChanged = _zChanged = false;
            return true;
        }

        bool ok = false;
        if (isAtEnd() && n >= 0) {
            _linearposition = _endposition;
            _xChanged = _yChanged = _zChanged = true;
            return false;
        }
        if ( _flow == fXYZ) {
            ok = moveXYZ(n);
        }
        else if ( _flow == fZXY){
            ok = moveZXY(n);
        }
        return ok;
    }


private:


    bool moveZXY(qint64 delta){
        _z += delta;
        _linearposition += delta * _box.xlength() * _box.ylength();
        _zChanged = true;
        _xChanged = _yChanged = false;
        if (_selectionIndex < 0){
            if ( _z > _endz || _z < _box.min_corner().z){
                return moveXY(delta);
            }
        }
        return true;
    }

    bool moveYXZ(qint64 delta){
        _y += delta;
        _linearposition += delta * _box.xlength();
        _yChanged = true;
        _xChanged = _zChanged = false;
        if (_selectionIndex < 0){
            if ( _y > _endy || _y < _box.min_corner().y){
                return moveXZ(delta);
            }
        }
        return true;
    }

    bool moveXYZ(qint64 delta) {
        _x += delta;
        _linearposition += delta;
        _xChanged = true;
        _yChanged = _zChanged = false;
        if ( _selectionIndex < 0){
            if ( _x > _endx || _z > _endz || _x < _box.min_corner().x) {
                return moveYZ(delta);
            }
        } else {
            qint64 selectionPix = (qint64)_selectionPixels[_y].size();
            if (  selectionPix  == 0 ){
                 _x = _endx + 1;
                 if(!moveYZ(delta))
                    return false;
                _selectionIndex = 0;
                _insideSelection = false;
            }
            else if ( _x == _selectionPixels[_y][_selectionIndex]){ // passed a boundary on this row
                _insideSelection = !_insideSelection;

                if (!_insideSelection ) {
                    if(!move2NextSelection(delta))
                        return false;
                }else
                    ++_selectionIndex;
            }
        }
        return true;
    }

    bool moveYZ(qint64 delta);
    bool moveXY(qint64 delta);
    bool moveXZ(qint64 delta);
    bool move2NextSelection(qint64 delta);
    void cleanUp4PolyBoundaries(const std::vector<Ilwis::Pixel> &selectionPix, geos::geom::Geometry *selection);
};

inline Ilwis::PixelIterator begin(const Ilwis::IRasterCoverage& raster) {
    return PixelIterator(raster);
}

inline Ilwis::PixelIterator end(const Ilwis::IRasterCoverage& raster) {
    PixelIterator iter(raster);
    return iter.end();
}

inline Ilwis::PixelIterator end(const Ilwis::PixelIterator& iter) {
    return iter.end();
}

}
#endif // PixelIterator_H
