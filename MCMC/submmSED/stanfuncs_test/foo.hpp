// Code generated by Stan version 2.17.0

// [[Rcpp::depends(rstan)]]
#include <stan/model/standalone_functions_header.hpp>

namespace foo_functions {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using namespace stan::math;

// [[Rcpp::export]]
boost::ecuyer1988 __create_rng(int seed) {
  return(boost::ecuyer1988(seed));
}

typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_d;
typedef Eigen::Matrix<double,1,Eigen::Dynamic> row_vector_d;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_d;

stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "unknown file name");
    reader.add_event(3, 3, "end", "unknown file name");
    return reader;
}

template <typename T0__>
Eigen::Matrix<typename boost::math::tools::promote_args<T0__>::type, Eigen::Dynamic,1>
foo(const Eigen::Matrix<T0__, Eigen::Dynamic,1>& x, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

    int current_statement_begin__ = -1;
    try {

        current_statement_begin__ = 2;
        return stan::math::promote_scalar<fun_return_scalar_t__>(multiply(2,x));
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}


struct foo_functor__ {
    template <typename T0__>
        Eigen::Matrix<typename boost::math::tools::promote_args<T0__>::type, Eigen::Dynamic,1>
    operator()(const Eigen::Matrix<T0__, Eigen::Dynamic,1>& x, std::ostream* pstream__) const {
        return foo(x, pstream__);
    }
};

// [[Rcpp::export]]
extern "C"
Eigen::Matrix<double, Eigen::Dynamic,1>
foo(const Eigen::Matrix<double, Eigen::Dynamic,1>& x, std::ostream* pstream__ = 0){
  return foo<double>(x, pstream__);
}

 }
