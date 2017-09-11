#include <RcppEigen.h>
#include <boost/math/tools/promotion.hpp>
#include <stan/math/prim/scal.hpp>
#include <stan/math.hpp>
// [[Rcpp::depends(RcppEigen)]]
//[[Rcpp::depends(BH)]]
using namespace Rcpp;
typedef Eigen::Map<Eigen::ArrayXd> mapa;
typedef Eigen::Map<Eigen::MatrixXd> mapmat;
typedef Eigen::Map<Eigen::VectorXd> mapvec;



//' evd_dnorm
//' 
//' This function computes the RSSp negative log likelihood
//' @param par a length 2 numeric vector.  par[1] corresponds to the sigma_u parameter, and par[2] corresponds to the confounding parameter
//' @param dvec. The eigenvalues of the LD matrix
//' @param quh The precomputed matrix vector product Q%*%u_hat (passed as a vector)
//[[Rcpp::export]]
double evd_dnorm(const Eigen::Map<Eigen::ArrayXd> par,const Eigen::Map<Eigen::ArrayXd> dvec, const Eigen::Map<Eigen::ArrayXd> quh){
  const double sigu=par(0);
  const double varu=sigu*sigu;
  const double a=(par.size()>1)?par(1):0;
  double tsum = ((dvec*dvec*varu+dvec+a).log()).sum();
  double tprod = ((quh*(1/(dvec*dvec*varu+dvec+a)))*(quh)).sum();
  return -0.5*(tsum+tprod);
}




//' evd_dnorm
//' 
//' This function computes the RSSp negative log likelihood gradient
//' @param par a length 2 numeric vector.  par[1] corresponds to the sigma_u parameter, and par[2] corresponds to the confounding parameter
//' @param dvec. The eigenvalues of the LD matrix
//' @param quh The precomputed matrix vector product Q%*%u_hat (passed as a vector)
//[[Rcpp::export]]
Eigen::ArrayXd evd_dnorm_grad(const Eigen::Map<Eigen::ArrayXd> par,const Eigen::Map<Eigen::ArrayXd> dvec, const Eigen::Map<Eigen::ArrayXd> quh){
  const double sigu=par(0);
  const double varu=sigu*sigu;
  const double a=(par.size()>1)?par(1):0;
  
  double sgrad = -((dvec.square()*sigu*(a+dvec.square()*varu+dvec-quh.square()))/(a+dvec.square()*varu+dvec).square()).sum();
  double  agrad = -((0.5*(a+dvec.square()*varu+dvec-quh.square()))/(a+dvec.square()*varu+dvec).square()).sum();

  Eigen::ArrayXd retvec(2);
  retvec(0)=sgrad;
  retvec(1)=agrad;
  return(retvec);
}


//Borrowed heavily from https://arxiv.org/pdf/1509.07164.pdf

struct evd_dens {
  const Eigen::Map<Eigen::ArrayXd> dvec;
  const Eigen::Map<Eigen::ArrayXd> quh;
  const size_t p;
  evd_dens(const Eigen::Map<Eigen::ArrayXd> dvec_,const Eigen::Map<Eigen::ArrayXd> quh_): dvec(dvec_),quh(quh_),p(quh.size()){}
  
  template <typename T>
  T operator()(const Eigen::Matrix<T,Eigen::Dynamic,1> &theta) const{
    using std::log;
    T sigu = theta[0];
    // Rcpp::Rcout<<"sigu:"<<sigu<<std::endl;
    // Rcpp::Rcout<<"theta[0]:"<<theta[0]<<std::endl;
    T a=(theta.size()>1)?theta[1]:0;
    // std::vector<T> tlsum(p);
    // std::vector<T> tlprod(p)
    //    Rcpp::Rcout<<"a:"<<a<<std::endl;
    T lp=0;
    for(size_t i=0;i<p;i++){
      lp+=-0.5*log(dvec[i]*dvec[i]*(sigu*sigu)+dvec[i]+a);
      lp+=-0.5*(quh[i]*quh[i])/(dvec[i]*dvec[i]*sigu*sigu+dvec[i]+a);
      
    }

    // Rcpp::Rcout<<"tsum:"<<lp<<std::endl;
//    lp+=-0/5*((quh*(1/(dvec*dvec*(sigu*sigu)+dvec+a)))*(quh)).sum();
//    lp*=-0.5;
    return(lp);
  }
  
};



//' evd_dnorm_grad_stan
//' 
//' This is an attempt to use stan's AD features to calculate a gradient 
//' for the RSSp likelihood
//' 
//' //' @param par a length 2 numeric vector.  par[1] corresponds to the sigma_u parameter, and par[2] corresponds to the confounding parameter
//' @param dvec. The eigenvalues of the LD matrix
//' @param quh The precomputed matrix vector product Q%*%u_hat (passed as a vector)
//[[Rcpp::export]]
Eigen::ArrayXd evd_dnorm_grad_stan(const Eigen::Map<Eigen::ArrayXd> par,const Eigen::Map<Eigen::ArrayXd> dvec, const Eigen::Map<Eigen::ArrayXd> quh){
//  Rcpp::Rcout<<"Function started!"<<std::endl;
  const double sigu=par(0);
//  const double varu=sigu*sigu;
  const double a=(par.size()>1)?par(1):0;
  // Rcpp::Rcout<<"sigu is :"<<sigu<<std::endl;
  double fx=0;
  Eigen::Matrix <double,Eigen::Dynamic,1> grad_fx;
  Eigen::Matrix <double,Eigen::Dynamic,1> param(2);
  param(0)=sigu;
  param(1)=a;
  // Rcpp::Rcout<<"param:"<<param<<std::endl;
  evd_dens f(dvec,quh);
//  Rcpp::Rcout<<"Struct constructed!"<<std::endl;
  stan::math::gradient(f,param,fx,grad_fx);
//  Rcpp::Rcout<<"Gradient computed!"<<std::endl;
  Eigen::ArrayXd retvec(grad_fx.size()+1);
  retvec(0)=fx;
  retvec.segment(1,grad_fx.size())=grad_fx;
  return(retvec);
}
