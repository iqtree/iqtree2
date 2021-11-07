#pragma once
#ifndef SUBCLASS_POINTER_VECTOR_H
#define SUBCLASS_POINTER_VECTOR_H

template <class T, class S> class SubclassPointerVector: public S {
    //
    //A subclass of vector class S, of pointers,
    //where the entries in S are pointer types of T
    //(where T is a subclass of S::value_type).
    //
    //eg. SubclassPointerVector<PhyloNode, NodeVector>
    //    is a subclass of NodeVector, such that operator[] and
    //    begin() and end() work in terms of PhyloNode* rather
    //    than Node*.  So client code won't need to cast
    //    iterator dereferences or [] entry-lookups to PhyloNode*.
    //
    //
public:
    typedef S super;
    typedef SubclassPointerVector<T, S> this_type;
    typedef T* value_type;
    using typename S::size_type;

    SubclassPointerVector() : super() {}
    SubclassPointerVector(const this_type& rhs) : super(rhs) {}
    this_type& operator= (const this_type& rhs){
        super::operator=(rhs);
        return *this;
    }

    explicit SubclassPointerVector(const S& rhs) : super(rhs) {}
    this_type& operator= (const S& rhs) {
        super::operator=(rhs);
        return *this;
    }

    class reference {
        private:
            S& to_vector;
            size_type at_index;
        public:
            reference(S& vector, size_type index):
                to_vector(vector), at_index(index) {}
            operator T*()   { typename S::value_type s = to_vector[at_index]; return dynamic_cast<T*> ( s ) ; }
            T* operator->() { typename S::value_type s = to_vector[at_index]; return dynamic_cast<T*> ( s ) ; }
            reference& operator= (T* new_value) {
                to_vector[at_index] = new_value;
                return *this;
            }
            reference& operator= (const reference& elsewhere) {
                to_vector[at_index] = elsewhere.to_vector[elsewhere.at_index];
                return *this;
            }
    };

    class iterator : public super::iterator {
        public:
            typedef typename super::iterator::difference_type step;
            iterator() : super::iterator() {}
            explicit iterator( typename super::iterator i) : super::iterator(i) {}
            T* operator*()  { return dynamic_cast<T*> ( super::iterator::operator*() ); }
            iterator& operator+=(step move) { super::iterator::operator+=(move); return *this;   }
            iterator  operator+(step move)  { return iterator(super::iterator::operator+(move)); }
    };
    class const_iterator: public super::const_iterator {
        public:
            typedef typename super::const_iterator::difference_type step;
            const_iterator() : super::const_iterator() {}
            explicit const_iterator( typename super::const_iterator i) : super::const_iterator(i) {};
            T* operator*() { return dynamic_cast<T*>( super::const_iterator::operator*() ); }
            const_iterator& operator+=(step move) { super::const_iterator::operator+=(move); return *this; }
            const_iterator  operator+(step move) { return iterator(super::iterator::operator+(move));  }
    };
    iterator begin() { return iterator( super::begin() ); }
    iterator end()   { return iterator( super::end() ); }
    const_iterator begin() const { return const_iterator ( super::begin() ); }
    const_iterator end() const { return const_iterator ( super::end() ); }
    T* operator[] (typename super::size_type i) const {
        return dynamic_cast<T*>( super::operator[] (i) );
    }
    reference operator[] (size_type i) {
        return reference(*(dynamic_cast<S*>(this)), i );
    }
    this_type& reverseAll() {
        std::reverse(super::begin(), super::end());
        return *this;
    }
};

#endif //SUBCLASS_POINTER_VECTOR_H
