
/* This file is generated from numpy_distutils/system_info.py */
#ifdef __CPLUSPLUS__
extern "C" {
#endif
#include "Python.h"
static PyMethodDef module_methods[] = { {NULL,NULL} };
PyMODINIT_FUNC initatlas_version(void) {
  void ATL_buildinfo(void);
  ATL_buildinfo();
  Py_InitModule("atlas_version", module_methods);
}
#ifdef __CPLUSCPLUS__
}
#endif
