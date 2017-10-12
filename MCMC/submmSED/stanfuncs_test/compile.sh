PYSTAN=/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/pystan
clang++ -dynamiclib -o blackbody.dylib  -I/Users/jaffe/home/dev/stan/stan/src/ \
-I${PYSTAN}/stan/lib/stan_math/lib/boost_1.64.0 \
-I${PYSTAN}/stan/lib/stan_math/lib/eigen_3.3.3/ \
-I${PYSTAN}/stan/lib/stan_math \
-I${PYSTAN}/stan/lib/stan_math/lib/cvodes_2.9.0/include \
./blackbody.cpp

# -I/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7 \
