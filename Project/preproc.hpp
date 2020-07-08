#include "H5Cpp.h"

using namespace std;
using namespace H5;

template <typename T> extern void readHDF5(hid_t &file_id, const char* dataset_name, T* store_addr, const string dtype);
void readHDF5_str(H5File& file, const char* dataset_name, char** store_addr);