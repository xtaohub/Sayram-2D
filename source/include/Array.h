#ifndef ARRAY_H_
#define ARRAY_H_

/*
  
  We define a very basic multi-dimensional array container

*/

#include <array>
#include <iostream>
#include <iterator>
#include <cassert>

#ifdef DEBUG
#define CHECK(A) assert(A)
#else
#define CHECK(A)
#endif

typedef long size_tp;
typedef std::array<size_tp, 1> Shape1;  // not the definition in array std::array<data_type,array_size>
typedef std::array<size_tp, 2> Shape2;
typedef std::array<size_tp, 3> Shape3;

template<class T, int ndim_>
class Array{
public:

  Array() { 
    for (int i=0; i<ndim_; ++i) shape_[i] = 0;
    setup_storage(); 
  }
  
  explicit Array(size_tp n0) {
    assert(ndim_ == 1); 
    shape_[0] = n0;
    setup_storage();
  }

  Array(size_tp n0, size_tp n1) {
    assert(ndim_ == 2);

    shape_[0] = n0;
    shape_[1] = n1;

    setup_storage();
  }
  
  Array(size_tp n0, size_tp n1, size_tp n2) {
    assert(ndim_ == 3); 

    shape_[0] = n0;
    shape_[1] = n1;
    shape_[2] = n2;

    setup_storage();
  }

  explicit Array(const std::array<size_tp, ndim_>& shape){
    shape_ = shape;
    setup_storage(); 
  }

  Array(const Array<T, ndim_>& b){
    shape_ = b.shape(); 
    setup_storage();
    
    std::copy(b.data(), b.data() + b.size(), data_);
  }

  Array<T,ndim_>& operator=(const Array<T,ndim_>& rhs){
    assert(shape_ == rhs.shape());
    if (this != &rhs) std::copy(rhs.data(), rhs.data() + rhs.size(), data_);

    return *this; 
  }

  Array<T,ndim_>& operator=(const T& rhs){
    std::fill(data_, data_ + size_, rhs);
    return *this;
  }

  void resize(size_tp n0) {
    assert(ndim_ == 1);

    if (n0 != size()) {
      delete [] data_;
      shape_[0] = n0;
      setup_storage();
    }
    
  }

  void resize(size_tp n0, size_tp n1){
    assert(ndim_ == 2);

    if (shape_[0] != n0 || shape_[1] != n1){
      delete [] data_;

      shape_[0] = n0;
      shape_[1] = n1;
    
      setup_storage();
    }
    
  }

  void resize(size_tp n0, size_tp n1, size_tp n2){
    assert(ndim_ == 3);

    if (shape_[0] != n0 || shape_[1] != n1 || shape_[2] != n2){
      delete [] data_;

      shape_[0] = n0;
      shape_[1] = n1;
      shape_[2] = n2;
    
      setup_storage();
    }
    
  }

  void resize(const std::array<size_tp, ndim_>& shape){

    if (shape != shape_) {
      shape_ = shape;
      delete [] data_; 
      setup_storage();
    }
    
  }

  const std::array<size_tp,ndim_>& shape() const { return shape_; }
  size_tp shape(int i) const { return shape_[i]; }  

  constexpr int ndim() const { return ndim_; }
  constexpr size_tp length() const { return size_; }
  constexpr size_tp nrows() const { return shape_[0]; }
  constexpr size_tp ncols() const { return shape_[1]; }
  constexpr size_tp size() const { return size_; }

  T& operator()(size_tp i0) {
    CHECK(ndim_ == 1);
    CHECK(i0 >= 0 && i0 < shape(0)); 

    return data_[i0*stride_[0]];
  }

  const T& operator()(size_tp i0) const {
    CHECK(ndim_ == 1);
    CHECK(i0 >= 0 && i0 < shape(0)); 
    return data_[i0*stride_[0]];
  }

  T& operator()(size_tp i0, size_tp i1) {
    CHECK(ndim_ == 2);
    CHECK(i0 >= 0 && i0 < shape(0) && i1 >= 0 && i1 < shape(1)); 
    return data_[i0*stride_[0] + i1*stride_[1]];
  }

  const T& operator()(size_tp i0, size_tp i1) const {
    CHECK(ndim_ == 2);
    CHECK(i0 >= 0 && i0 < shape(0) && i1 >= 0 && i1 < shape(1)); 
    return data_[i0*stride_[0] + i1*stride_[1]];
    
  }

  T& operator()(size_tp i0, size_tp i1, size_tp i2) {
    CHECK(ndim_ == 3);
    CHECK(i0 >= 0 && i0 < shape(0) && i1 >= 0 && i1 < shape(1) && i2 >= 0 && i2 < shape(2)); 
    return data_[i0*stride_[0] + i1*stride_[1] + i2*stride_[2]];
  }

  const T& operator()(size_tp i0, size_tp i1, size_tp i2) const {
    CHECK(ndim_ == 3);
    CHECK(i0 >= 0 && i0 < shape(0) && i1 >= 0 && i1 < shape(1) && i2 >= 0 && i2 < shape(2)); 
    return data_[i0*stride_[0] + i1*stride_[1] + i2*stride_[2]];
  }

  T* data() { return data_; }
  const T* data() const { return data_; }

  ~Array() { delete [] data_; }


private:

  T* data_;
  std::array<size_tp,ndim_> shape_;
  std::array<size_tp,ndim_> stride_;
  size_tp size_; 

  void setup_storage() {

    size_ = shape_[0]; 
    for (int i=1; i<ndim_; ++i) size_ *= shape_[i]; 

    if (size_ > 0) {
      data_ = new (std::nothrow) T [size_];
      assert(data_);

      for (int i=0; i<ndim_; ++i){
	stride_[i] = 1;
	for (int j=i+1; j<ndim_; ++j)
	  stride_[i] *= shape_[j];
      }
    }
    
    else {
      size_ = 0; 
      data_ = NULL;
      for (int i=0; i<ndim_; ++i) { shape_[i] = 0; stride_[i] = 0; }
    }
  }

};

template <class T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& arr){ 
  std::copy(arr.cbegin(), arr.cend(), std::ostream_iterator<T>(os, " "));
  return os;
}

template <class T, int ndim>
std::ostream& operator<<(std::ostream& os, const Array<T, ndim>& arr){

  switch (ndim){
  case 1:
    for (size_tp i = 0; i < arr.length(); ++i) os << arr(i) << " ";
    break; 

  case 2:
    for (size_tp i = 0; i < arr.nrows(); ++i){
      for (size_tp j = 0; j < arr.ncols(); ++j) os << arr(i, j) << " ";
      os << std::endl;
    }
    break;

  case 3:

    for (size_tp i = 0; i < arr.shape(0); ++i){
      for (size_tp j = 0; j < arr.shape(1); ++j){
	for (size_tp k = 0; k < arr.shape(2); ++k)
	  os << arr(i, j, k) << " ";
	os << std::endl;
      }
      os << std::endl;
    }
    
    break;
  }

  return os;
}


#endif
