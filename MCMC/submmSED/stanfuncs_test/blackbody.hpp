// Code generated by Stan version 2.17.0

// [[Rcpp::depends(rstan)]]
#include <stan/model/standalone_functions_header.hpp>

namespace blackbody_functions { 
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
    reader.add_event(13, 13, "end", "unknown file name");
    return reader;
}

template <typename T0__, typename T1__, typename T2__>
typename boost::math::tools::promote_args<T0__, T1__, T2__>::type
greybody(const T0__& beta,
             const T1__& T,
             const T2__& nu, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 5;
        local_scalar_t__ h_over_k;
        (void) h_over_k;  // dummy to suppress unused var warning

        stan::math::initialize(h_over_k, DUMMY_VAR__);
        stan::math::fill(h_over_k,DUMMY_VAR__);
        stan::math::assign(h_over_k,0.047992369999999999);
        current_statement_begin__ = 6;
        local_scalar_t__ nu_bar;
        (void) nu_bar;  // dummy to suppress unused var warning

        stan::math::initialize(nu_bar, DUMMY_VAR__);
        stan::math::fill(nu_bar,DUMMY_VAR__);
        stan::math::assign(nu_bar,1000);
        current_statement_begin__ = 8;
        local_scalar_t__ x;
        (void) x;  // dummy to suppress unused var warning

        stan::math::initialize(x, DUMMY_VAR__);
        stan::math::fill(x,DUMMY_VAR__);
        stan::math::assign(x,((h_over_k * nu) / T));
        current_statement_begin__ = 9;
        local_scalar_t__ x_bar;
        (void) x_bar;  // dummy to suppress unused var warning

        stan::math::initialize(x_bar, DUMMY_VAR__);
        stan::math::fill(x_bar,DUMMY_VAR__);
        stan::math::assign(x_bar,((h_over_k * nu_bar) / T));


        current_statement_begin__ = 11;
        return stan::math::promote_scalar<fun_return_scalar_t__>(((pow((nu / nu_bar),(3 + beta)) * stan::math::expm1(x_bar)) / stan::math::expm1(x)));
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}


struct greybody_functor__ {
    template <typename T0__, typename T1__, typename T2__>
        typename boost::math::tools::promote_args<T0__, T1__, T2__>::type
    operator()(const T0__& beta,
             const T1__& T,
             const T2__& nu, std::ostream* pstream__) const {
        return greybody(beta, T, nu, pstream__);
    }
};

// [[Rcpp::export]]
double
greybody(const double& beta,
             const double& T,
             const double& nu, std::ostream* pstream__ = 0){
  return greybody<double, double, double>(beta, T, nu, pstream__);
}

 } 
