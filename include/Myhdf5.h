#ifndef MYHDF5_H_
#define MYHDF5_H_

// define a parallel hdf5 wrapper
// this assumes that mpi has been initiated.

#include <string>
#include <cassert>
#include "hdf5.h"
#include "mpi.h"

using std::string;

class Myhdf5{
public:
  Myhdf5(const string& filename);
  ~Myhdf5(); 

  hid_t fid() const { return fid_; }
  
  void add_group(const string& groupname) const;
  
  const string& filename() const { return h5f_; } 

  hid_t cdouble() const { return H5T_NATIVE_DOUBLE; }
  hid_t cint() const { return H5T_NATIVE_INT; }
  hid_t clong() const { return H5T_NATIVE_LONG; }

  template <class data_tp>
  void write_scalar_sp(data_tp scalar, hid_t h5_datatype, const string& dset_name, int root=0) const;

  template <class data_tp>
  void write_vector_sp(data_tp* ptr, long size, hid_t h5_datatype, const string& dset_name, int root=0) const;

  template <class data_tp>
  void write_matrix_sp(data_tp* ptr, long size_row, long size_col, hid_t h5_datatype, const string& dset_name, int root=0) const;

  template <class data_tp>
  void write_cube_sp(data_tp* ptr, long size0, long size1, long size2, hid_t h5_datatype, const string& dset_name, int root=0) const;  

  template <class data_tp>
  void write_vector(data_tp* ptr, // pointer to the data buff
  		    long size_local, 
  		    long size_global,
  		    long start, 
  		    hid_t h5_datatype,
  		    const string& dset_name,
  		    MPI_Comm comm=MPI_COMM_WORLD) const;

  template <class data_tp>
  void write_matrix(data_tp* ptr, // pointer to the data buff
  		    long size_local[2], 
  		    long size_global[2],
  		    long start[2], 
  		    hid_t h5_datatype,
  		    const string& dset_name,
  		    MPI_Comm comm=MPI_COMM_WORLD) const;

  
private:
  int rank_; // used for single process read & write 
  string h5f_;
  hid_t fid_;
  
  void open(const string& filename); 
};


inline Myhdf5::Myhdf5(const string& filename){
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  open(filename);
}

inline void Myhdf5::open(const string& filename){
  hid_t fapl_id;
  h5f_ = filename + ".h5"; 

  /* 
   * Set up file access property list with parallel I/O access
   */

  fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);

  /*
   * Create a new file collectively and release property list identifier.
   */
  fid_ = H5Fcreate(h5f_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
  assert(fid_ > 0);

  H5Pclose(fapl_id);

}

inline Myhdf5::~Myhdf5(){
  H5Fclose(fid_);
}

inline void Myhdf5::add_group(const string& gname) const{

  // /* 
  //  * Set up file access property list with parallel I/O access
  //  */

  // hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  // H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);

  // hid_t fid = H5Fopen(h5f_.c_str(), H5F_ACC_RDWR, fapl_id);
  // assert(fid > 0);
  // H5Pclose(fapl_id);

  // /*
  //  * Create a new file collectively and release property list identifier.
  //  */
  // hid_t fid = H5Fcreate(h5f_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
  // assert(fid > 0);
  // H5Pclose(fapl_id);

  hid_t gid; 
  
  gid = H5Gcreate(fid_, gname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(gid);

}

template <class data_tp>
void Myhdf5::write_scalar_sp(data_tp scalar, hid_t h5_datatype, const string& dset_name, int root) const{

  hsize_t dim[1];
  dim[0] = 1;

  hid_t dataspace, dataset_id;
  dataspace = H5Screate_simple(0, dim, NULL);

  dataset_id = H5Dcreate(fid_, dset_name.c_str(), h5_datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // set up dataset transfer property list for collective MPI I/O
  hid_t xferplist = H5Pcreate(H5P_DATASET_XFER);

  // H5Pset_dxpl_mpio(xferplist, H5FD_MPIO_COLLECTIVE);
  // H5Dwrite(dataset_id, h5_datatype, H5S_ALL, H5S_ALL, xferplist, &scalar);

  H5Pset_dxpl_mpio(xferplist, H5FD_MPIO_INDEPENDENT);

  if (rank_ == root)
    H5Dwrite(dataset_id, h5_datatype, H5S_ALL, H5S_ALL, xferplist, &scalar);  

  H5Sclose(dataspace);
  H5Dclose(dataset_id);
  
}

template <class data_tp>
void Myhdf5::write_vector_sp(data_tp* ptr, long size, hid_t h5_datatype, const string& dset_name, int root) const{

  hsize_t dim[1];
  dim[0] = size;

  hid_t dataspace, dataset_id;
  dataspace = H5Screate_simple(1, dim, NULL);

  dataset_id = H5Dcreate(fid_, dset_name.c_str(), h5_datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // set up dataset transfer property list for collective MPI I/O
  hid_t xferplist = H5Pcreate(H5P_DATASET_XFER);

  H5Pset_dxpl_mpio(xferplist, H5FD_MPIO_INDEPENDENT);

  if (rank_ == root)
    H5Dwrite(dataset_id, h5_datatype, H5S_ALL, H5S_ALL, xferplist, ptr);  
  
  H5Sclose(dataspace);
  H5Dclose(dataset_id);

}

template <class data_tp>
void Myhdf5::write_matrix_sp(data_tp* ptr, long size_row, long size_col, hid_t h5_datatype, const string& dset_name, int root) const{

  hsize_t dim[2];
  dim[0] = size_row;
  dim[1] = size_col;

  hid_t dataspace, dataset_id;
  dataspace = H5Screate_simple(2, dim, NULL);

  dataset_id = H5Dcreate(fid_, dset_name.c_str(), h5_datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // set up dataset transfer property list for collective MPI I/O
  hid_t xferplist = H5Pcreate(H5P_DATASET_XFER);

  // H5Pset_dxpl_mpio(xferplist, H5FD_MPIO_COLLECTIVE);
  // H5Dwrite(dataset_id, h5_datatype, H5S_ALL, H5S_ALL, xferplist, ptr);

  H5Pset_dxpl_mpio(xferplist, H5FD_MPIO_INDEPENDENT);

  if (rank_ == root)
    H5Dwrite(dataset_id, h5_datatype, H5S_ALL, H5S_ALL, xferplist, ptr);    

  H5Sclose(dataspace);
  H5Dclose(dataset_id);

}

template <class data_tp>
void Myhdf5::write_cube_sp(data_tp* ptr, long size0, long size1, long size2, hid_t h5_datatype, const string& dset_name, int root) const{

  hsize_t dim[3];
  dim[0] = size0;
  dim[1] = size1;
  dim[2] = size2;

  hid_t dataspace, dataset;

  /*
   * Describe the size of the array and create the data space for fixed
   * size dataset.
   */
  int RANK=3; 
  dataspace = H5Screate_simple(RANK, dim, NULL);
  
  /*
   * Create a new dataset within the file using defined dataspace and
   * default dataset creation properties.
   */
  dataset = H5Dcreate(fid_, dset_name.c_str(), h5_datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


  hid_t xferplist;   

  xferplist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(xferplist, H5FD_MPIO_INDEPENDENT);  
  
  /*
   * Write the data to the dataset using default transfer properties if current rank is rank.
   */

  if (rank_ == root)
    H5Dwrite(dataset, h5_datatype, H5S_ALL, H5S_ALL, xferplist, ptr);

  H5Sclose(dataspace);
  H5Dclose(dataset);
  
}

template <class data_tp>
void Myhdf5::write_vector(data_tp* ptr,
			  long size_local,
			  long size_global,
			  long start_loc, 
			  hid_t h5_datatype,
			  const string& dset_name, MPI_Comm comm) const{


  hid_t dset_id;
  hid_t filespace, memspace;
  hid_t	plist_id;

  const int drank = 1; 
  hsize_t dimsf[drank];                 /* global dataset dimensions */
  hsize_t dimsm[drank];                 /* local dataset dimensions */
  hsize_t count[drank];	                /* hyperslab selection parameters */
  hsize_t stride[drank];
  hsize_t block[drank];
  hsize_t start[drank];
  
  /*
     Create the dataspace for the dataset.
  */
  dimsf[0] = size_global; 
  dimsm[0] = size_local;
  
  filespace = H5Screate_simple(drank, dimsf, NULL); 
  memspace  = H5Screate_simple(drank, dimsm, NULL);

  /*
   * Create the dataset with default properties.
   */
  dset_id = H5Dcreate(fid_, dset_name.c_str(), h5_datatype, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  count[0] = 1;
  stride[0] = 1;
  block[0] = dimsm[0];
  start[0] = start_loc;

  /*
   * Select hyperslab in the file.
   */
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, stride, count, block);
  
  /*
   * Create property list for collective dataset write.
   */
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT); 

  H5Dwrite(dset_id, h5_datatype, memspace, filespace, plist_id, ptr);

  H5Pclose(plist_id);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Sclose(memspace);

}

template <class data_tp>
void Myhdf5::write_matrix(data_tp* ptr, // pointer to the data buff
			  long size_local[2], 
			  long size_global[2],
			  long start_loc[2], 
			  hid_t h5_datatype,
			  const string& dset_name,
			  MPI_Comm comm) const{

    hid_t dset_id;
    hid_t filespace, memspace;
    hid_t plist_id;

    const int drank = 2; 
    hsize_t dimsf[drank];                 /* global dataset dimensions */
    hsize_t dimsm[drank];                 /* local dataset dimensions */
    hsize_t count[drank];	          /* hyperslab selection parameters */
    hsize_t stride[drank];
    hsize_t block[drank];
    hsize_t start[drank];
  
    /*
      Create the dataspace for the dataset.
    */
    dimsf[0] = size_global[0];
    dimsf[1] = size_global[1];
    dimsm[0] = size_local[0];
    dimsm[1] = size_local[1];
    
    filespace = H5Screate_simple(drank, dimsf, NULL); 
    memspace  = H5Screate_simple(drank, dimsm, NULL); 

    /*
     * Create the dataset with default properties.
     */
    dset_id = H5Dcreate(fid_, dset_name.c_str(), h5_datatype, filespace,
    			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    count[0] = 1;
    count[1] = 1;
    stride[0] = 1;
    stride[1] = 1;

    block[0] = size_local[0];
    block[1] = size_local[1];

    start[0] = start_loc[0];
    start[1] = start_loc[1]; 
  
    /*
     * Select hyperslab in the file.
     */
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, stride, count, block); 

    /*
     * Create property list for collective dataset write.
     */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);     

    H5Dwrite(dset_id, h5_datatype, memspace, filespace, plist_id, ptr); 

    H5Pclose(plist_id);
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    
}


#endif
