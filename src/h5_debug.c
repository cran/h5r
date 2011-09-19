/**
 * Code which is used to debug things. 
 */
#include <hdf5.h>
#include <Rinternals.h>    
#include <R.h>

#define myfree(p) (free(p))
#define myalloc(s) (malloc(s))

void h5R_allocate_finalizer(SEXP eptr) {
    char* vector = R_ExternalPtrAddr(eptr);
    myfree(vector);
    R_ClearExternalPtr(eptr);
}

SEXP h5R_allocate_meg() {
    char* vector = (char*) myalloc(1048576);
    for (int j = 0; j < 1048576; j++) {
    	vector[j] = 'c';
    }
    SEXP e_ptr = R_MakeExternalPtr(vector, R_NilValue, R_NilValue); 
    PROTECT(e_ptr);
    R_RegisterCFinalizerEx(e_ptr, h5R_allocate_finalizer, TRUE);
    UNPROTECT(1);
    return e_ptr;
}

SEXP h5R_allocate_k() {
    char* vector = (char*) myalloc(1024);
    for (int j = 0; j < 1024; j++) {
    	vector[j] = 'c';
    }
    SEXP e_ptr = R_MakeExternalPtr(vector, R_NilValue, R_NilValue); 
    PROTECT(e_ptr);
    R_RegisterCFinalizerEx(e_ptr, h5R_allocate_finalizer, TRUE);
    UNPROTECT(1);
    return e_ptr;
}

SEXP h5R_allocate_gig() {
    char* vector = (char*) myalloc(1073741824);
    for (int j = 0; j < 1073741824; j++) {
    	vector[j] = 'c';
    }
    SEXP e_ptr = R_MakeExternalPtr(vector, R_NilValue, R_NilValue); 
    PROTECT(e_ptr);
    R_RegisterCFinalizerEx(e_ptr, h5R_allocate_finalizer, TRUE);
    UNPROTECT(1);
    return e_ptr;
}
