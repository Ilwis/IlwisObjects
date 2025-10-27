#ifndef RELATIVEVECTOR_H
#define RELATIVEVECTOR_H
#include <iostream>
#include <vector>
#include <stdexcept>
#include <iterator>

namespace Ilwis {
template<typename T>
class RelativeVector {
private:
    std::vector<T> dataVector;
    int start_index = 0;

public:
    class iterator {
    public:
        using iterator_category = std::random_access_iterator_tag;
        using value_type = T;
        using difference_type = std::ptrdiff_t;
        using pointer = T*;
        using reference = T&;

        // Constructor
        explicit iterator(typename std::vector<T>::iterator it) : m_it(it) {}

        // Dereference operator
        reference operator*() const {
            return *m_it;
        }

        pointer operator->() const {
            return &(*m_it);
        }

        // Pre-increment
        iterator& operator++() {
            ++m_it;
            return *this;
        }

        // Post-increment
        iterator operator++(int) {
            iterator temp = *this;
            ++(*this);
            return temp;
        }

        // ... (Other operators like --, +, -, ==, !=, <, >, <=, >=) ...
        // For a full implementation, you'd need to overload all of these.
        // Let's just do a few to illustrate.

        // Pre-decrement
        iterator& operator--() {
            --m_it;
            return *this;
        }

        // Equality check
        bool operator==(const iterator& other) const {
            return m_it == other.m_it;
        }

        // Inequality check
        bool operator!=(const iterator& other) const {
            return !(*this == other);
        }

        // Addition operator
        iterator operator+(difference_type n) const {
            return iterator(m_it + n);
        }

        // Subtraction operator
        iterator operator-(difference_type n) const {
            return iterator(m_it - n);
        }

        // Difference between two iterators
        difference_type operator-(const iterator& other) const {
            return m_it - other.m_it;
        }

    private:
        typename std::vector<T>::iterator m_it;
    };

    class const_iterator {
     public:
         using iterator_category = std::random_access_iterator_tag;
         using value_type = T;
         using difference_type = std::ptrdiff_t;
         using pointer = const T*;
         using reference = const T&;

         // Constructor
         explicit const_iterator(typename std::vector<T>::const_iterator it) : m_it(it) {}

         // Allows implicit conversion from a non-const iterator
         const_iterator(const typename RelativeVector<T>::iterator& it) : m_it(it.get_internal_iterator()) {}


         // Dereference operator
         reference operator*() const {
             return *m_it;
         }

         pointer operator->() const {
             return &(*m_it);
         }

         // Pre-increment
         const_iterator& operator++() {
             ++m_it;
             return *this;
         }

         // Post-increment
         const_iterator operator++(int) {
             const_iterator temp = *this;
             ++(*this);
             return temp;
         }

         // Pre-decrement
         const_iterator& operator--() {
             --m_it;
             return *this;
         }

         // Post-decrement
         const_iterator operator--(int) {
             const_iterator temp = *this;
             --(*this);
             return temp;
         }

         // Addition assignment
         const_iterator& operator+=(difference_type n) {
             m_it += n;
             return *this;
         }

         // Subtraction assignment
         const_iterator& operator-=(difference_type n) {
             m_it -= n;
             return *this;
         }

         // Subscript operator
         reference operator[](difference_type n) const {
             return m_it[n];
         }

         // Equality check
         bool operator==(const const_iterator& other) const {
             return m_it == other.m_it;
         }

         // Inequality check
         bool operator!=(const const_iterator& other) const {
             return !(*this == other);
         }

         // Comparison operators
         bool operator<(const const_iterator& other) const {
             return m_it < other.m_it;
         }
         bool operator>(const const_iterator& other) const {
             return m_it > other.m_it;
         }
         bool operator<=(const const_iterator& other) const {
             return m_it <= other.m_it;
         }
         bool operator>=(const const_iterator& other) const {
             return m_it >= other.m_it;
         }

         // Addition operator
         const_iterator operator+(difference_type n) const {
             return const_iterator(m_it + n);
         }

         // Subtraction operator
         const_iterator operator-(difference_type n) const {
             return const_iterator(m_it - n);
         }

         // Difference between two iterators
         difference_type operator-(const const_iterator& other) const {
             return m_it - other.m_it;
         }

     private:
         typename std::vector<T>::const_iterator m_it;

         // Friend class to allow non-const iterator to access m_it
         friend class RelativeVector<T>::iterator;

         // Private getter for the non-const iterator to use for conversion
         typename std::vector<T>::const_iterator get_internal_iterator() const {
             return m_it;
         }
     };
    // Constructor with a starting index
    explicit RelativeVector(int start_idx=0) : start_index(start_idx) {}

    // Constructor with size and starting index
    RelativeVector(size_t size, int start_idx)
       : dataVector(size), start_index(start_idx) {}

    // Constructor with size, initial value, and starting index
    RelativeVector(size_t size, const T& value, int start_idx = 0)
        : dataVector(size, value), start_index(start_idx) {}

    // Copy constructor
    RelativeVector(const RelativeVector& other)
        : dataVector(other.dataVector), start_index(other.start_index) {}

    // Move constructor
    RelativeVector(RelativeVector&& other) noexcept
        : dataVector(std::move(other.dataVector)), start_index(other.start_index) {}

    // Copy assignment operator
    RelativeVector& operator=(const RelativeVector& other) {
        if (this != &other) {
            dataVector = other.dataVector;
            start_index = other.start_index;
        }
        return *this;
    }

    // Move assignment operator
    RelativeVector& operator=(RelativeVector&& other) noexcept {
        if (this != &other) {
            dataVector = std::move(other.dataVector);
            start_index = other.start_index;
        }
        return *this;
    }

    // Access element with bounds checking
    T& at(int index) {
        if (index < start_index || index >= start_index + static_cast<int>(dataVector.size())) {
            throw std::out_of_range("Index out of bounds");
        }
        return dataVector[index - start_index];
    }

    const T& at(int index) const {
        if (index < start_index || index >= start_index + static_cast<int>(dataVector.size())) {
            throw std::out_of_range("Index out of bounds");
        }
        return dataVector[index - start_index];
    }

    // Access element without bounds checking
    T& operator[](int index) {
        return dataVector[index - start_index];
    }

    const T& operator[](int index) const {
        return dataVector[index - start_index];
    }

    // Returns the size of the vector
    size_t size() const {
        return dataVector.size();
    }

    // Returns true if the vector is empty, false otherwise
    bool empty() const {
        return dataVector.empty();
    }

    // Clears the vector
    void clear() {
        dataVector.clear();
    }

    // Resizes the vector
    void resize(size_t new_size) {
        dataVector.resize(new_size);
    }

    // Resizes the vector with an initial value
    void resize(size_t new_size, const T& value) {
        dataVector.resize(new_size, value);
    }

    // Add an element to the end of the vector
    void push_back(const T& value) {
        dataVector.push_back(value);
    }

    // Get the starting index of the vector
    int startIndex() const {
        return start_index;
    }

    void startIndex(int idx)  {
        start_index = idx;
    }

    T* data(){
        return dataVector.data();
    }

    const T* data() const{
        return dataVector.data();
    }

    iterator begin() {
        return iterator(dataVector.begin());
    }

    const_iterator begin() const {
        return const_iterator(dataVector.cbegin());
    }

    iterator end() {
        return iterator(dataVector.end());
    }

    const_iterator end() const {
        return const_iterator(dataVector.cend());
    }

};

}
#endif // RELATIVEVECTOR_H
