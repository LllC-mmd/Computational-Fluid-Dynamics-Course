#include <iostream>
#include <string>
#include "preproc.hpp"


template <typename T>
void readHDF5(hid_t &file_id, const char* dataset_name, T* store_addr, const string dtype){
    herr_t status;
    hid_t dataset_id, dataspace_id, file_dataspace_id;
    hsize_t* dims, *max_dims;
    int rank;
    
    // open existing dataset
    dataset_id = H5Dopen(file_id, dataset_name, H5P_DEFAULT);
    // determine dataset parameters
    file_dataspace_id = H5Dget_space(dataset_id);
    // get the number of dimensions in the dataspace
    rank = H5Sget_simple_extent_ndims(file_dataspace_id); 
    // get the max size for each dimension in max_dims
    dims = (hsize_t*) malloc(rank*sizeof(hsize_t));
    max_dims = (hsize_t*) malloc(rank*sizeof(hsize_t));
    H5Sget_simple_extent_dims(file_dataspace_id, dims, max_dims);
    // create memory dataspace
    dataspace_id = H5Screate_simple(rank, dims, max_dims);
    // read matrix data from file
    // HDF5 works only with single continuous array in memory
    if(dtype=="int"){
        status = H5Dread(dataset_id, H5T_NATIVE_INT32, dataspace_id, file_dataspace_id, H5P_DEFAULT, store_addr);
    }
    else if(dtype=="float"){
        status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, dataspace_id, file_dataspace_id, H5P_DEFAULT, store_addr);
    }
    else if(dtype=="double"){
        status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, dataspace_id, file_dataspace_id, H5P_DEFAULT, store_addr);
    }

    // release resources and close file
    status = H5Dclose(dataset_id);
    status = H5Sclose(dataspace_id);
    status = H5Sclose(file_dataspace_id);
}

template void readHDF5<int>(hid_t &file_id, const char* dataset_name, int* store_addr, const string dtype);
template void readHDF5<double>(hid_t &file_id, const char* dataset_name, double* store_addr, const string dtype);

// read 1D string array
void readHDF5_str(H5File& file, const char* dataset_name, char** store_addr){
    StrType dt;
    DataSpace dspace;
    int rank, ndims;
    hsize_t* dims, *max_dims;

    // open existing HDF5 file and the dataset in it
    DataSet dset = file.openDataSet(dataset_name);
    // 1) get the data type of the dataset
    dt = dset.getStrType();
    // 2) get the dataspace of the dataset
    dspace = dset.getSpace();
    // -----get the the number of dimensions in the dataspace
    rank = dspace.getSimpleExtentNdims();

    // create memory dataspace
    dims = (hsize_t*) malloc(rank*sizeof(hsize_t));
    max_dims = (hsize_t*) malloc(rank*sizeof(hsize_t));
    ndims = dspace.getSimpleExtentDims(dims, max_dims);
    DataSpace mspace(rank, max_dims);

    // read HDF5 data
    dset.read(store_addr, dt, mspace, dspace);

    dset.close();
    mspace.close();
    dspace.close();
}


/*
int main() {
    // read input
    double* A;
    int num_edge = 12;
    int col_edge = 2;
    A = (double*) malloc(num_edge*col_edge*sizeof(double));
    readHDF5("dataset/node.h5", "/NodeXY", A, "double");
    // do something with A
    for (int i=0; i<num_edge; i++){
        std::cout << "Row: "<< i << std::endl;
        for (int j=0; j<col_edge; j++){
            std::cout<< A[i*col_edge+j] << " ";
        }
        std::cout << std::endl;
    }
    free(A);

    return 0;
}
*/
