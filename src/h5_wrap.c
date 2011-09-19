/**
 * R/C Interface code for HDF5 file format. 
 */
#include <hdf5.h>
#include <Rinternals.h>    
#include <R.h>

#define DEBUG 1
#define HID(argname) (((h5_holder*) R_ExternalPtrAddr(argname))->id)
#define NM(argname) (CHAR(STRING_ELT(argname, 0)))
#define SUCCESS ScalarLogical(1)
#define FAILURE ScalarLogical(0)

typedef struct h5_holder {
  int is_file;
  hid_t id;
} h5_holder;

void h5R_finalizer(SEXP h5_obj) {
  h5_holder* h = (h5_holder*) R_ExternalPtrAddr(h5_obj);
  if (! h) {
    return;
  }
  if (h->is_file == 1) {
    H5Fflush(HID(h5_obj), H5F_SCOPE_GLOBAL);
    H5Fclose(HID(h5_obj));
  } else {
    switch (H5Iget_type(HID(h5_obj))) {
    case H5I_DATASET:
      H5Dclose(HID(h5_obj));
      break;
    case H5I_ATTR:
      H5Aclose(HID(h5_obj));
      break;
    case H5I_GROUP:
      H5Gclose(HID(h5_obj));
      break;
    default:
      break;
    }
  }
  Free(h);
  R_ClearExternalPtr(h5_obj);
}

SEXP _h5R_make_holder (hid_t id, int is_file) {
  if (id < 0) {
    return R_NilValue;
  } 
  h5_holder* holder = (h5_holder*) Calloc(1, h5_holder);
  holder->id = id;
  holder->is_file = is_file;
  SEXP e_ptr = R_MakeExternalPtr(holder, R_NilValue, R_NilValue); 
  PROTECT(e_ptr);
  R_RegisterCFinalizerEx(e_ptr, h5R_finalizer, TRUE);
  UNPROTECT(1); 
  return e_ptr;
}

SEXP h5R_flush(SEXP h5_file) {
  H5Fflush(HID(h5_file), H5F_SCOPE_GLOBAL);
  return SUCCESS;
}

SEXP h5R_open(SEXP filename, SEXP mode) {
  int _mode_ = (INTEGER(mode)[0] == 1) ? H5F_ACC_RDWR : H5F_ACC_RDONLY;
  return _h5R_make_holder(H5Fopen(NM(filename), _mode_, H5P_DEFAULT), 1);
}

SEXP h5R_create(SEXP filename) {
  return _h5R_make_holder(H5Fcreate(NM(filename), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT), 1);
}

SEXP h5R_get_group(SEXP h5_obj, SEXP group_name) {
  return _h5R_make_holder(H5Gopen2(HID(h5_obj), NM(group_name), H5P_DEFAULT), 0);
}

SEXP h5R_create_group(SEXP h5_obj, SEXP group_name) {
  return _h5R_make_holder(H5Gcreate2(HID(h5_obj), NM(group_name), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), 0);
}

SEXP h5R_get_dataset(SEXP h5_obj, SEXP dataset_name) {
  return _h5R_make_holder(H5Dopen2(HID(h5_obj), NM(dataset_name), H5P_DEFAULT), 0);
}

SEXP h5R_get_attr(SEXP h5_obj, SEXP attr_name) {
  return _h5R_make_holder(H5Aopen(HID(h5_obj), NM(attr_name), H5P_DEFAULT), 0);
}

hid_t _h5R_get_space(SEXP h5_obj) {
  hid_t space = -1;

  switch (H5Iget_type(HID(h5_obj))) {
  case H5I_DATASET:
    space = H5Dget_space(HID(h5_obj));
    break;
  case H5I_ATTR:
    space = H5Aget_space(HID(h5_obj));
    break;
  default:
    error("Unknown object in %s.\n", __func__);
  }
  return space;
}

hid_t _h5R_get_type(SEXP h5_obj) {
  hid_t dtype = -1;

  switch (H5Iget_type(HID(h5_obj))) {
  case H5I_DATASET:
    dtype = H5Dget_type(HID(h5_obj));
    break;
  case H5I_ATTR:
    dtype = H5Aget_type(HID(h5_obj));
    break;
  default:
    error("Unknown object in %s.\n", __func__);
  }
  return dtype;
}

int _h5R_get_size (SEXP h5_obj) {
  int s = -1;

  switch (H5Iget_type(HID(h5_obj))) {
  case H5I_DATASET:
    s = H5Tget_size(H5Dget_type(HID(h5_obj)));
    break;
  case H5I_ATTR:
    s = H5Tget_size(H5Aget_type(HID(h5_obj)));
    break;
  default:
    error("Unknown object in %s.\n", __func__);
  }
  return s;
}

int _h5R_is_vlen (SEXP h5_obj) {
  hid_t hdty = _h5R_get_type(h5_obj);
  int rval = H5Tis_variable_str(hdty);
  H5Tclose(hdty);
  return rval;
}

int _h5R_get_ndims(SEXP h5_obj) {
  hid_t space = _h5R_get_space(h5_obj);
  int ndims = H5Sget_simple_extent_ndims(space);
  H5Sclose(space);
  return ((ndims < 0) ? 1 : ndims);
}

int _h5R_get_nelts(SEXP h5_obj) {
  int v = 1; int i;
  int ndims = _h5R_get_ndims(h5_obj);
  hid_t space = _h5R_get_space(h5_obj);
  hsize_t* dims = (hsize_t* ) Calloc(ndims, hsize_t);
  H5Sget_simple_extent_dims(space, dims, NULL);
  for (i = 0; i < ndims; i++)
    v *= dims[i];
  Free(dims);
  H5Sclose(space);
  return(v);
}

SEXP h5R_get_type(SEXP h5_obj) {
  SEXP dtype = R_NilValue;
  hid_t hdty = _h5R_get_type(h5_obj);
  PROTECT(dtype = ScalarInteger(H5Tget_class(hdty)));
  H5Tclose(hdty);
  UNPROTECT(1);
  return(dtype);
}

SEXP h5R_get_dims(SEXP h5_obj) {
  int i; SEXP res; 
  int ndims = _h5R_get_ndims(h5_obj);
  hid_t space = _h5R_get_space(h5_obj);

  hsize_t* dims = (hsize_t* ) Calloc(ndims, hsize_t);
  H5Sget_simple_extent_dims(space, dims, NULL);
    
  PROTECT(res = allocVector(INTSXP, ndims)); 
  for (i = 0; i < ndims; i++)
    INTEGER(res)[i] = dims[i];
  UNPROTECT(1);
    
  Free(dims);
  H5Sclose(space);
    
  return res;
}

/** ******************************************************************************
 *
 * Read API
 *
 ******************************************************************************* **/
SEXP _h5R_read_integer(SEXP h5_obj, int nelts, hid_t memspace, hid_t filespace) {
  SEXP dta;
  PROTECT(dta = allocVector(INTSXP, nelts));
  if (H5Iget_type(HID(h5_obj)) == H5I_DATASET) {
    H5Dread(HID(h5_obj), H5T_NATIVE_INT, memspace, filespace, H5P_DEFAULT, INTEGER(dta));
  } else {
    H5Aread(HID(h5_obj), H5T_NATIVE_INT, INTEGER(dta)); 
  }
  UNPROTECT(1);
  return dta;
}

SEXP _h5R_read_float(SEXP h5_obj, int nelts, hid_t memspace, hid_t filespace) {
  SEXP dta;
  PROTECT(dta = allocVector(REALSXP, nelts));
  if (H5Iget_type(HID(h5_obj)) == H5I_DATASET) {
    H5Dread(HID(h5_obj), H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, REAL(dta));
  } else {
    H5Aread(HID(h5_obj), H5T_NATIVE_DOUBLE, REAL(dta));
  }
  UNPROTECT(1);
  return dta;
}

SEXP _h5R_read_string(SEXP h5_obj, int nelts, hid_t memspace, hid_t filespace) {
  SEXP dta = R_NilValue;
  void* buf; 
  char** rdata = (char **) Calloc(nelts, char*);
  hid_t memtype = H5Tcopy (H5T_C_S1);

  if (! _h5R_is_vlen(h5_obj)) { 
    int sdim = _h5R_get_size(h5_obj) + 1; 
    rdata[0] = (char *) Calloc(nelts * sdim, char);
    for (int i = 1; i < nelts; i++)
      rdata[i] = rdata[0] + i * sdim;
    H5Tset_size (memtype, sdim);
    buf = rdata[0];
  } 
  else {
    H5Tset_size(memtype, H5T_VARIABLE);
    buf = rdata;
  }

  if (H5Iget_type(HID(h5_obj)) == H5I_DATASET) {
    H5Dread(HID(h5_obj), memtype, memspace, filespace, H5P_DEFAULT, buf);
  } else {
    H5Aread(HID(h5_obj), memtype, buf);
  }

  PROTECT(dta = allocVector(STRSXP, nelts));
  for (int i = 0; i < nelts; i++) {
    if (rdata[i]) 
      SET_STRING_ELT(dta, i, mkChar(rdata[i])); 
  }
  UNPROTECT(1); 

  if (_h5R_is_vlen(h5_obj)) {
    // It doesn't like the H5S_ALL space in the vlen_reclaim.
    if (memspace == H5S_ALL) {
      hid_t space = _h5R_get_space(h5_obj);
      H5Dvlen_reclaim (memtype, space, H5P_DEFAULT, buf);
      H5Sclose(space);
    } else {
      H5Dvlen_reclaim (memtype, memspace, H5P_DEFAULT, buf);
    }
  } else {
    Free(buf);
  }
  Free(rdata);

  return dta;
}

SEXP _h5R_read_compound(SEXP h5_obj, int nelts, hid_t memspace, hid_t filespace) {
  return R_NilValue;
}

SEXP h5R_read_dataset_all(SEXP h5_dataset) {
  SEXP dta = R_NilValue;

  switch (INTEGER(h5R_get_type(h5_dataset))[0]) {
  case H5T_INTEGER: 
    dta = _h5R_read_integer(h5_dataset, _h5R_get_nelts(h5_dataset), H5S_ALL, H5S_ALL);
    break;
  case H5T_FLOAT:
    dta = _h5R_read_float(h5_dataset, _h5R_get_nelts(h5_dataset), H5S_ALL, H5S_ALL);
    break;
  case H5T_STRING:
    dta = _h5R_read_string(h5_dataset, _h5R_get_nelts(h5_dataset), H5S_ALL, H5S_ALL);
    break;
  case H5T_COMPOUND:
    dta = _h5R_read_compound(h5_dataset, _h5R_get_nelts(h5_dataset), H5S_ALL, H5S_ALL);
    break;
  default:
    dta = R_NilValue;
  }
  return dta;
}

SEXP h5R_read_dataset(SEXP h5_dataset, SEXP _offsets, SEXP _counts) {
  SEXP dta = R_NilValue;
  hid_t filespace = -1, memspace = -1;
  int nelts = 1;

  hsize_t* _h_offsets = (hsize_t*) Calloc(length(_counts), hsize_t);
  hsize_t* _h_counts  = (hsize_t*) Calloc(length(_counts), hsize_t);
  
  for (int i = 0; i < length(_counts); i++) {
    nelts *= INTEGER(_counts)[i];
    _h_offsets[i] = INTEGER(_offsets)[i];
    _h_counts[i] = INTEGER(_counts)[i];
  }
  filespace = _h5R_get_space(h5_dataset);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, _h_offsets, NULL, _h_counts, NULL);
  memspace = H5Screate_simple(length(_counts), _h_counts, NULL);
  
  switch (INTEGER(h5R_get_type(h5_dataset))[0]) {
  case H5T_INTEGER: 
    dta = _h5R_read_integer(h5_dataset, nelts, memspace, filespace);
    break;
  case H5T_FLOAT:
    dta = _h5R_read_float(h5_dataset, nelts, memspace, filespace);
    break;
  case H5T_STRING:
    dta = _h5R_read_string(h5_dataset, nelts, memspace, filespace);
    break;
  case H5T_COMPOUND:
    dta = _h5R_read_compound(h5_dataset, nelts, memspace, filespace);
    break;
  default:
    dta = R_NilValue;
  }
  Free(_h_offsets);
  Free(_h_counts);
  H5Sclose(memspace);
  H5Sclose(filespace);
  
  if (dta == R_NilValue) {
    error("Unsupported class in %s\n", __func__);
  }
  return dta;
}

SEXP h5R_read_points(SEXP h5_dataset, SEXP _points, SEXP _nr, SEXP _nc) {
  SEXP dta = R_NilValue;
  hid_t filespace = -1, memspace = -1;

  int nr = INTEGER(_nr)[0];
  int nc = INTEGER(_nc)[0];

  hsize_t* points = (hsize_t*) Calloc(nr*nc, hsize_t);
  for (int i = 0; i < nr*nc; i++)	
    points[i] = INTEGER(_points)[i];
  hsize_t hnr = nr;

  filespace = _h5R_get_space(h5_dataset);
  H5Sselect_elements(filespace, H5S_SELECT_SET, nr, points);
  memspace = H5Screate_simple(nc, &hnr, NULL);

  switch (INTEGER(h5R_get_type(h5_dataset))[0]) {
  case H5T_INTEGER: 
    dta = _h5R_read_integer(h5_dataset, length(_points), memspace, filespace);
    break;
  case H5T_FLOAT:
    dta = _h5R_read_float(h5_dataset, length(_points), memspace, filespace);
    break;
  case H5T_STRING:
    dta = _h5R_read_string(h5_dataset, length(_points), memspace, filespace);
    break;
  case H5T_COMPOUND:
    dta = _h5R_read_compound(h5_dataset, length(_points), memspace, filespace);
    break;
  default:
    dta = R_NilValue;
  }

  /** clean up. **/
  Free(points);
  H5Sclose(memspace);
  H5Sclose(filespace);
    
  if (dta == R_NilValue) {
    error("Unsupported class in %s\n", __func__);
  }
  return dta;
}

SEXP h5R_read_attr(SEXP h5_attr) {
  SEXP dta = R_NilValue;

  switch (INTEGER(h5R_get_type(h5_attr))[0]) {
  case H5T_INTEGER:
    dta = _h5R_read_integer(h5_attr, _h5R_get_nelts(h5_attr), H5S_ALL, H5S_ALL);
    break;
  case H5T_FLOAT:
    dta = _h5R_read_float(h5_attr, _h5R_get_nelts(h5_attr), H5S_ALL, H5S_ALL);
    break;
  case H5T_STRING:
    dta =  _h5R_read_string(h5_attr, _h5R_get_nelts(h5_attr), H5S_ALL, H5S_ALL);
    break;
  default:
    error("Unsupported class in %s.\n", __func__);
  }
  return dta;
}

/** This is an optimization to a very common use case of 1-d slabs. **/
SEXP h5R_read_1d_slabs(SEXP h5_dataset, SEXP _offsets, SEXP _counts) {
  int rlen = length(_counts);
  SEXP r_lst, _SEXP_offsets, _SEXP_counts;
  int i;
  int* counts = INTEGER(_counts);
  int* offsets = INTEGER(_offsets);
    
  PROTECT(r_lst = allocVector(VECSXP, rlen));
  PROTECT(_SEXP_offsets = allocVector(INTSXP, 1));
  PROTECT(_SEXP_counts  = allocVector(INTSXP, 1));
    
  for (i = 0; i < rlen; i++) {
    INTEGER(_SEXP_offsets)[0] = offsets[i];
    INTEGER(_SEXP_counts)[0] = counts[i];
    SET_VECTOR_ELT(r_lst, i, h5R_read_dataset(h5_dataset, _SEXP_offsets, _SEXP_counts));
  }
  UNPROTECT(3);
  return(r_lst);
}


/** ******************************************************************************
 *
 * Write API
 *
 ******************************************************************************* **/
SEXP h5R_create_dataset(SEXP h5_obj, SEXP name, SEXP h5_type, SEXP dims, SEXP chunks) {
  int i; 
    
  hsize_t* current_dims = (hsize_t*) Calloc(length(dims), hsize_t);
  hsize_t* max_dims = (hsize_t*) Calloc(length(dims), hsize_t);
  hsize_t* chunk_lens = (hsize_t*) Calloc(length(chunks), hsize_t);

  for (i = 0; i < length(dims); i++) {
    current_dims[i] = INTEGER(dims)[i];
    max_dims[i]     = H5S_UNLIMITED;
    chunk_lens[i]   = INTEGER(chunks)[i];
  }
    
  hid_t memtype = -1;

  switch (INTEGER(h5_type)[0]) {
  case H5T_INTEGER: 
    memtype = H5T_NATIVE_INT;
    break;
  case H5T_FLOAT:
    memtype = H5T_NATIVE_DOUBLE;
    break;
  case H5T_STRING:
    memtype = H5Tcopy (H5T_C_S1);
    H5Tset_size (memtype, H5T_VARIABLE);   
    break;
  default:
    error("Unsupported class in %s.\n", __func__);
  }
   
  hid_t cparms = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(cparms, length(chunks), chunk_lens);
    
  hid_t dataspace = H5Screate_simple(length(dims), current_dims, max_dims);
  hid_t dataset   = H5Dcreate2(HID(h5_obj), NM(name), memtype, 
			       dataspace, H5P_DEFAULT, cparms, H5P_DEFAULT);
    
  H5Pclose(cparms);
  H5Sclose(dataspace);

  Free(max_dims);
  Free(current_dims);
  Free(chunk_lens);
  
  return _h5R_make_holder(dataset, 0);
}

SEXP h5R_write_slab(SEXP h5_dataset, SEXP _offsets, SEXP _counts, SEXP data) {
  int __ERROR__ = 0;
  hid_t space = -1, memspace = -1, memtype = -1;
  int i; 
  void* _data; 
  char** tmp;

  int* offsets  = INTEGER(_offsets);
  int* counts   = INTEGER(_counts);

  /** I'm surprised I have to do this, but it seems to be necessary. **/
  hsize_t* _h_offsets = (hsize_t*) Calloc(length(_counts), hsize_t);
  hsize_t* _h_counts  = (hsize_t*) Calloc(length(_counts), hsize_t);

  for (i = 0; i < length(_counts); i++) {
    _h_offsets[i] = offsets[i];
    _h_counts[i]  = counts[i];
  }

  space = _h5R_get_space(h5_dataset);
  H5Sselect_hyperslab(space, H5S_SELECT_SET, _h_offsets, NULL, _h_counts, NULL);
  memspace = H5Screate_simple(length(_counts), _h_counts, NULL);

  int memtype_int = INTEGER(h5R_get_type(h5_dataset))[0];
    
  switch (memtype_int) {
  case H5T_INTEGER: 
    memtype = H5T_NATIVE_INT;
    _data = (void*) INTEGER(data);
    break;
  case H5T_FLOAT:
    memtype = H5T_NATIVE_DOUBLE;
    _data = (void*) REAL(data);
    break;
  case H5T_STRING:
    memtype = H5Tcopy (H5T_C_S1);
    H5Tset_size (memtype, H5T_VARIABLE);   
    tmp = (char**) Calloc(length(data), char*);
    for (i = 0; i < length(data); i++)
      tmp[i] = (char*) CHAR(STRING_ELT(data, i));
    _data = (void*) tmp;
    break;
  default:
    __ERROR__ = 1;
  }

  if (__ERROR__ == 0) {
    H5Dwrite(HID(h5_dataset), memtype, memspace, space, H5P_DEFAULT, _data);
    if (memtype_int == H5T_STRING) Free(tmp);
  }

  /** clean up. **/
  Free(_h_offsets);
  Free(_h_counts);
  H5Sclose(memspace);
  H5Sclose(space);
    
  if (__ERROR__ == 1) {
    error("Unsupported class in %s\n", __func__);
    return FAILURE;
  }
  return SUCCESS;
}

SEXP h5R_delete_object(SEXP h5_obj, SEXP name) {
  H5Ldelete(HID(h5_obj), NM(name), H5P_DEFAULT);
  return SUCCESS;
}

SEXP h5R_delete_attribute(SEXP h5_obj, SEXP name) {
  H5Adelete(HID(h5_obj), NM(name));
  return SUCCESS;
}

SEXP h5R_create_attribute(SEXP h5_obj, SEXP name, SEXP h5_type, SEXP dims) {
  int i; 
  hsize_t* current_dims = (hsize_t*) Calloc(length(dims), hsize_t);
  hsize_t* max_dims = (hsize_t*) Calloc(length(dims), hsize_t);

  for (i = 0; i < length(dims); i++) {
    current_dims[i] = INTEGER(dims)[i];
    max_dims[i]     = H5S_UNLIMITED;
  }
    
  hid_t memtype = -1;

  switch (INTEGER(h5_type)[0]) {
  case H5T_INTEGER: 
    memtype = H5T_NATIVE_INT;
    break;
  case H5T_FLOAT:
    memtype = H5T_NATIVE_DOUBLE;
    break;
  case H5T_STRING:
    memtype = H5Tcopy (H5T_C_S1);
    H5Tset_size (memtype, H5T_VARIABLE);   
    break;
  default:
    error("Unsupported class in %s.\n", __func__);
  }
    
  hid_t dataspace = H5Screate_simple(length(dims), current_dims, max_dims);
  hid_t attribute = H5Acreate2(HID(h5_obj), NM(name), memtype, 
			       dataspace, H5P_DEFAULT, H5P_DEFAULT);
    
  H5Sclose(dataspace);
  Free(max_dims);
  Free(current_dims);
  
  return _h5R_make_holder(attribute, 0);
}

SEXP h5R_write_attribute(SEXP h5_attr, SEXP data) {
  hid_t memtype = -1;
  void* buf = NULL; 
  char** tmp = NULL;
  int i;
  int atype = INTEGER(h5R_get_type(h5_attr))[0];
    
  switch (atype) {
  case H5T_INTEGER:
    buf = INTEGER(data);
    memtype = H5T_NATIVE_INT;
    break;
  case H5T_FLOAT:
    buf = REAL(data);
    memtype = H5T_NATIVE_DOUBLE;
    break;
  case H5T_STRING:
    memtype = H5Tcopy (H5T_C_S1);
    H5Tset_size (memtype, H5T_VARIABLE);   
    tmp = (char**) Calloc(length(data), char*);
    for (i = 0; i < length(data); i++)
      tmp[i] = (char*) CHAR(STRING_ELT(data, i));
    buf = (void*) tmp;
    break;
  default:
    error("Unsupported class in %s.\n", __func__);
  }	
  H5Awrite(HID(h5_attr), memtype, buf);
    
  if (atype == H5T_STRING) {
    Free(tmp);
    H5Tclose(memtype);
  }
  H5Aclose(HID(h5_attr));
  return SUCCESS;
}


/** ******************************************************************************
 *
 * File Inspection
 *
 ******************************************************************************* **/
SEXP h5R_attribute_exists(SEXP h5_obj, SEXP name) {
  if (H5Aexists(HID(h5_obj), NM(name)) == 1) {
    return SUCCESS;
  } 
  else {
    return FAILURE;
  }
}

SEXP h5R_dataset_exists(SEXP h5_obj, SEXP name) {
  if (H5Lexists(HID(h5_obj), NM(name), H5P_DEFAULT) == 1) {
    return SUCCESS;
  } 
  else {
    return FAILURE;
  }
}

/** Iteration **/
typedef struct __index_and_SEXP__ {
  int  i;
  SEXP s;
} __index_and_SEXP__;

herr_t _h5R_count_func(hid_t loc_id, const char *name, const H5O_info_t *info,
		       void *operator_data) {
  int* counter = ((int *) operator_data);
  (*counter)++;
    
  return 0;
}

herr_t _h5R_capture_name(hid_t loc_id, const char *name, const H5O_info_t *info,
			 void *operator_data) {
  __index_and_SEXP__* od = (__index_and_SEXP__*) operator_data;
  SET_STRING_ELT(od->s, (od->i)++, mkChar(name));
    
  return 0;
}

herr_t _h5R_capture_name_and_type(hid_t loc_id, const char *name, const H5O_info_t *info,
				  void *operator_data) {
  __index_and_SEXP__* od = (__index_and_SEXP__*) operator_data;
  SEXP lst, str;

  PROTECT(lst = allocVector(VECSXP, 2));
  PROTECT(str = allocVector(STRSXP, 1));
  SET_STRING_ELT(str, 0, mkChar(name));
  SET_VECTOR_ELT(lst, 0, str);
  SET_VECTOR_ELT(lst, 1, ScalarInteger(info->type));
  SET_VECTOR_ELT(od->s, (od->i)++, lst);

  UNPROTECT(2);

  return 0;
}

SEXP h5R_list_contents(SEXP h5_obj) {
  int counter = 0; 
  __index_and_SEXP__* isxp = (__index_and_SEXP__*) Calloc(1, __index_and_SEXP__);
  SEXP dta;

  H5Ovisit (HID(h5_obj), H5_INDEX_NAME, H5_ITER_NATIVE, _h5R_count_func, (void*) &counter);

  PROTECT(dta = allocVector(VECSXP, counter));
  isxp->s = dta;
  isxp->i = 0;
  H5Ovisit (HID(h5_obj), H5_INDEX_NAME, H5_ITER_NATIVE, _h5R_capture_name_and_type, (void*) isxp);

  Free(isxp);
  UNPROTECT(1);

  return(dta);
}

SEXP h5R_list_attributes(SEXP h5_obj) {
  int counter = 0;
  hsize_t n   = 0;
  SEXP dta;

  H5Aiterate2(HID(h5_obj), H5_INDEX_NAME, H5_ITER_NATIVE, &n, (H5A_operator2_t) _h5R_count_func, (void*) &counter);
    
  __index_and_SEXP__* isxp = (__index_and_SEXP__*) Calloc(1, __index_and_SEXP__);
  PROTECT(dta = allocVector(STRSXP, counter));
  isxp->s = dta;
  isxp->i = 0;
    
  n = 0;
  H5Aiterate2(HID(h5_obj), H5_INDEX_NAME, H5_ITER_NATIVE, &n, 
	      (H5A_operator2_t) _h5R_capture_name, (void*) isxp);

  Free(isxp);
  UNPROTECT(1);

  return(dta);
}

herr_t _h5R_name_exists(hid_t loc_id, const char *name, const H5O_info_t *info,
			void *operator_data) {
  const char* probe = (const char*) operator_data;

  if (strcmp(name, probe) == 0) {
    return 1; //short-circuit.
  } else {
    return 0; //continue
  }
}
