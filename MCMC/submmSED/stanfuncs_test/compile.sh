CXX=clang++

PYSTAN=/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/pystan
STAN_MATH=${PYSTAN}/stan/lib/stan_math

${CXX} -dynamiclib -o blackbody.dylib \
-I/Users/jaffe/home/dev/stan/stan/src/ \
-I${STAN_MATH} -I${STAN_MATH}/lib/boost_1.64.0 -I${STAN_MATH}/lib/eigen_3.3.3/ -I${STAN_MATH}/lib/cvodes_2.9.0/include \
./blackbody.cpp

# -I/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7 \
#
