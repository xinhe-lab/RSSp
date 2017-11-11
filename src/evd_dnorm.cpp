#include <RSSp.h>
#include <boost/math/tools/promotion.hpp>
#include <stan/math/prim/scal.hpp>
#include <stan/math.hpp>
#include <stan/math/mix/mat.hpp>
// [[Rcpp::depends(RcppEigen)]]
//[[Rcpp::depends(BH)]]
using namespace Rcpp;

template<typename T> struct ParamArray{
  typedef typename Eigen::Map<Eigen::Array<double,T::value,1> > ParamType;
};


typedef std::integral_constant<int, 2> two_parameter_type;
typedef std::integral_constant<int, 1> one_parameter_type;

template <typename T> double RSSp_lik(const MapA par, const MapA dvec, const MapA quh);

template <> double RSSp_lik<two_parameter_type>(const MapA par, const MapA dvec, const MapA quh){
  const double varu=par(0);
  const double a=par(1);
  const double tsum = ((dvec*dvec*varu+dvec+a).log()).sum();
  const double tprod = (quh*(1/(dvec*dvec*varu+dvec+a))*quh).sum();
  return -0.5*(tsum+tprod);
}

template <> double RSSp_lik<one_parameter_type>(const MapA par, const MapA dvec, const MapA quh){
  const double varu=par(0);
  const double tsum = ((dvec*dvec*varu+dvec).log()).sum();
  const double tprod = (quh*(1/(dvec*dvec*varu+dvec))*quh).sum();
  return -0.5*(tsum+tprod);
}

//' RSSp marginalized likelihood function
//' 
//' Compute the negative of the marginalized RSSp log-likelihood
//' @template RSSp_stat
//[[Rcpp::export]]
double evd_dnorm(const MapA par,const MapA dvec, const MapA quh){
  const double varu=par(0);
  const double a=(par.size()>1)?par(1):0;
  const double tsum = ((dvec*dvec*varu+dvec+a).log()).sum();
  const double tprod = (quh*(1/(dvec*dvec*varu+dvec+a))*quh).sum();
  return -0.5*(tsum+tprod);
}





//' variance  
//' @template RSSp_stat
Eigen::ArrayXd evd_marg_var(const MapA par, const MapA dvec){
  const double varu=par(0);
  const double a=(par.size()>1)?par(1):0;
  return((varu*(a+dvec))/(a+dvec.square()*varu+dvec));
}

//' Compute the posterior mean of `crossprod(Q,u)` given eigenvalues, `crossprod(Q,u_hat)`,sigma_u^2 and confounding levels.
//' @template RSSp_stat
//[[Rcpp::export]]
Eigen::ArrayXd evd_post_mean(const MapA par,const MapA dvec, const MapA quh){
  const double varu=par(0);
  const double a=(par.size()>1)?par(1):0;
  return((dvec*varu/(a+dvec.square()*varu+dvec))*quh);
}





  


//' Compute the posterior variance of `crossprod(Q,u)` given eigenvalues, `crossprod(Q,u_hat)`,sigma_u^2 and confounding levels.
//' @template RSSp_stat
//[[Rcpp::export]]
Eigen::ArrayXd evd_post_var(const MapA par,const MapA dvec){
  const double varu=par(0);
  const double a=(par.size()>1)?par(1):0;
  return((varu*(a+dvec))/(a+dvec+varu));
}



//' evd_dnorm_grad
//' 
//' This function computes the RSSp negative log likelihood gradient
//' @template RSSp_stat
//[[Rcpp::export]]
Eigen::ArrayXd evd_dnorm_grad(const MapA par,const MapA dvec, const MapA quh){
  const bool useConfound = par.size()>1;
  const double varu=par(0);
  const double a=useConfound?par(1):0;
  
  const double sgrad = (-(dvec.square()*(a + dvec + dvec.square()*varu - quh.square()))/(2*(a + dvec + dvec.square()*varu).square())).sum();
  const double  agrad = (-(a+dvec.square()*varu+dvec-quh.square())/(2*(a+dvec.square()*varu+dvec).square())).sum();
  
  
  Eigen::ArrayXd retvec(useConfound?2:1);
  retvec(0)=sgrad;
  if(useConfound){
    retvec(1)=agrad;
  }
  return(retvec);
}







//Borrowed heavily from https://arxiv.org/pdf/1509.07164.pdf

struct evd_dens {
  const MapA dvec;
  const MapA quh;
  const size_t p;
  evd_dens(const MapA dvec_,const MapA quh_): dvec(dvec_),quh(quh_),p(quh.size()){}
  
  template <typename T>
  T operator()(const Eigen::Matrix<T,Eigen::Dynamic,1> &theta) const{
    using std::log;
    T varu = theta[0];
    // Rcpp::Rcout<<"sigu:"<<sigu<<std::endl;
    // Rcpp::Rcout<<"theta[0]:"<<theta[0]<<std::endl;
    T a=(theta.size()>1)?theta[1]:0;
    // std::vector<T> tlsum(p);
    // std::vector<T> tlprod(p)
    //    Rcpp::Rcout<<"a:"<<a<<std::endl;
    T lp=0;
    for(size_t i=0;i<p;i++){
      lp+=-0.5*log(dvec[i]*dvec[i]*(varu)+dvec[i]+a);
      lp+=-0.5*(quh[i]*quh[i])/(dvec[i]*dvec[i]*varu+dvec[i]+a);
      
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
//' @template RSSp_stat
//[[Rcpp::export]]
Eigen::ArrayXd evd_dnorm_grad_stan(const MapA par,const MapA dvec, const MapA quh){
  //  Rcpp::Rcout<<"Function started!"<<std::endl;
  const double varu=par(0);
  const double a=par(1);
  // Rcpp::Rcout<<"sigu is :"<<sigu<<std::endl;
  double fx=0;
  Eigen::Matrix <double,Eigen::Dynamic,1> grad_fx;
  Eigen::Matrix <double,Eigen::Dynamic,1> param(2);
  param(0)=varu;
  param(1)=a;
  // Rcpp::Rcout<<"param:"<<param<<std::endl;
  evd_dens f(dvec,quh);
  //  Rcpp::Rcout<<"Struct constructed!"<<std::endl;
  stan::math::gradient(f,param,fx,grad_fx);
  //  Rcpp::Rcout<<"Gradient computed!"<<std::endl;
  return(grad_fx.array());
}


//' evd_dnorm_hess_stan
//' 
//' This is an attempt to use stan's AD features to calculate a hessian 
//' for the RSSp likelihood
//' 
//' @template RSSp_stat
//[[Rcpp::export]]
Eigen::MatrixXd evd_dnorm_hess_stan(const MapA par,const MapA dvec, const MapA quh){
  //  Rcpp::Rcout<<"Function started!"<<std::endl;
  const double varu=par(0);
  const double a=par(1);
  // Rcpp::Rcout<<"sigu is :"<<sigu<<std::endl;
  double fx=0;
  Eigen::Matrix <double,Eigen::Dynamic,1> grad_fx;
  Eigen::Matrix <double,Eigen::Dynamic,Eigen::Dynamic> H_fx;
  
  Eigen::Matrix <double,Eigen::Dynamic,1> param(2);
  param(0)=varu;
  param(1)=a;
  // Rcpp::Rcout<<"param:"<<param<<std::endl;
  evd_dens f(dvec,quh);
  //  Rcpp::Rcout<<"Struct constructed!"<<std::endl;
  stan::math::hessian(
    f,
    param,
    fx,
    grad_fx,
    H_fx);
  //  Rcpp::Rcout<<"Gradient computed!"<<std::endl;

  return(H_fx);
}


//[[Rcpp::export]]
Eigen::MatrixXd  posterior_mean_D(const Eigen::ArrayXd sigu,const Eigen::ArrayXd confound,const Eigen::ArrayXd dvec){
  const size_t n_gene_ests=sigu.size();
  const size_t p =dvec.size();
  if(n_gene_ests!=confound.size()){
    Rcpp::stop("sigu and confound must be the same size");
  }
  
  Eigen::ArrayXd varu=sigu.square();
  return( (dvec.matrix()*varu.matrix().transpose()).array()/((((dvec.square().matrix()*varu.matrix().transpose()).array()).rowwise()+confound.transpose()).colwise()+dvec).array());
}
//[[Rcpp::export]]
Eigen::MatrixXd  posterior_mean_U(const MapA sigu,const MapA confound,const MapA dvec,const MapMat quh,const MapMat Q){
  Eigen::MatrixXd pD=posterior_mean_D(sigu,confound,dvec);
  return(Q*(pD.array()*quh.array()).matrix());
}

//[[Rcpp::export]]
Eigen::MatrixXd  posterior_mean_Beta(const MapA sigu,const MapA confound,const MapA dvec,const MapMat quh,const MapMat Q,const MapMat se){
  const size_t g=sigu.size();
  const size_t p=dvec.size();
  if((Q.rows()!=p)||(Q.cols()!=p)){
    Rcpp::stop("Eigenvector shape (Q) does not match eigenvalue shape (dvec)");
  }
  if(quh.rows()!=p ||(quh.cols()!=g)){
    Rcpp::stop("quh shape does not match eigenvalue shape (dvec) or sigu shape");
  }
  if(se.rows()!=p || se.cols()!=g){
    Rcpp::stop("se shape does not match eigenvalue shape (dvec) or sigu shape");
  }

  Eigen::MatrixXd pD=posterior_mean_D(sigu,confound,dvec);
  return(se.array()*(Q*(pD.array()*quh.array()).matrix()).array());
}
  
//[[Rcpp::export]]
Eigen::MatrixXd  posterior_mean_Y(const MapA sigu,const MapA confound,const MapA dvec,const MapMat quh,const MapMat Q,const MapMat se,const MapMat x){
  Eigen::MatrixXd pD=posterior_mean_D(sigu,confound,dvec);
  return(x*(se.array()*(Q*(pD.array()*quh.array()).matrix()).array()).matrix());
}





