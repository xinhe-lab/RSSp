#include <RSSp.h>
#include <boost/math/tools/promotion.hpp>
#include <stan/math/prim/scal.hpp>
#include <stan/math.hpp>
#include <stan/math/mix/mat.hpp>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(unwindProtect)]]
//[[Rcpp::depends(BH)]]
using namespace Rcpp;

template<typename T> struct ParamArray{
  typedef typename Eigen::Map<Eigen::Array<double,T::value,1> > ParamType;
};


template<int N>
double t_evd_dnorm(const Eigen::Array<double,N,1> cvec, const MapA D, const MapA quh);


template<int N>
double t_evd_dnorm(const Eigen::Array<double,N,1> cvec ,const MapA D, const MapA quh){
  const double p=static_cast<double>(D.size());
  double tot_sum=0;
  for(int i=0; i<p; i++){
    double tvar=D[i];
    for (int k = 0; k < N; k++) {
      tvar += cvec[k] * D[i] * std::pow(D[i], 2 - k);
    }
    tot_sum += std::log(tvar) + (quh[i] * quh[i]) / tvar;
  }
  return -(-0.5 * (tot_sum)-0.5 * p * log(2 * M_PI));
}

template<>
double t_evd_dnorm<1>(const Eigen::Array<double,1,1> cvec, const MapA D, const MapA quh){
  //  Rcpp::Rcout<<"1"<<std::endl;
  const double p=static_cast<double>(D.size());
  const double tsum = ((D*D*cvec[0]+D).log()+ (quh*(1/(D*D*cvec[0]+D))*quh)).sum();
  return -(-0.5*(tsum)-0.5*p*log(2*M_PI));
}

template <>
double t_evd_dnorm<2>(const Eigen::Array<double, 2, 1> cvec, const MapA D,
                      const MapA quh) {
  //  Rcpp::Rcout<<"2"<<std::endl;
  const double p = static_cast<double>(D.size());
  const double tsum =
      ((D + cvec[0] * D.square() + cvec[1] * D).log() +
       (quh * (1 / (D + cvec[0] * D.square() + cvec[1] * D)) * quh))
          .sum();
  return -(-0.5 * (tsum)-0.5 * p * log(2 * M_PI));
}

template<>
double t_evd_dnorm<3>(const Eigen::Array<double,3,1> cvec, const MapA D, const MapA quh){
  //  Rcpp::Rcout<<"3"<<std::endl;
  const double p=static_cast<double>(D.size());
  const double tsum =          (D+cvec[0]*D.square()+cvec[1]*D+cvec[2]).log().sum();
  const double tprod = (quh*(1/(D+cvec[0]*D.square()+cvec[1]*D+cvec[2]))*quh).sum();
  return -(-0.5*(tsum+tprod)-0.5*p*log(2*M_PI));
}


template<>
double t_evd_dnorm<4>(const Eigen::Array<double,4,1> cvec, const MapA D, const MapA quh){
  //  Rcpp::Rcout<<"4"<<std::endl;
  const double p=static_cast<double>(D.size());
  const double tsum =          (D+cvec[0]*D.square()+cvec[1]*D+cvec[2]+cvec[3]*D.pow(-1)).log().sum();
  const double tprod = (quh*(1/(D+cvec[0]*D.square()+cvec[1]*D+cvec[2]+cvec[3]*D.pow(-1)))*quh).sum();
  return -(-0.5*(tsum+tprod)-0.5*p*log(2*M_PI));
}





//' RSSp marginalized likelihood function
//'
//' Compute the negative of the marginalized RSSp log-likelihood
//' @template RSSp_stat
//' @export
//[[Rcpp::export]]
double evd_dnorm(const MapA par, const MapA D, const MapA quh) {
  const int num_param=par.size();
  switch(num_param){
  case 1:{
    Eigen::Array<double,1,1> tcvec(par);
    return(t_evd_dnorm(tcvec,D,quh));
  }
  case 2:{
    Eigen::Array<double,2,1> tcvec(par);
    return(t_evd_dnorm(tcvec,D,quh));
    return(t_evd_dnorm(tcvec,D,quh));
  }
  case 3:{
    Eigen::Array<double,3,1> tcvec(par);
    return(t_evd_dnorm(tcvec,D,quh));
  }
  case 4:{
    Eigen::Array<double,4,1> tcvec(par);
    return(t_evd_dnorm(tcvec,D,quh));
  }
  case 5:{
    Eigen::Array<double,5,1> tcvec(par);
    return(t_evd_dnorm(tcvec,D,quh));
  }
  default: { Rcpp::stop("too many confounding terms!"); }
  }
}

//Borrowed heavily from https://arxiv.org/pdf/1509.07164.pdf
struct evd_dens {
  const MapA D;
  const MapA quh;
  const size_t p;
  evd_dens(const MapA D_,const MapA quh_): D(D_),quh(quh_),p(quh.size()){}
  
  template <typename T>
  T operator()(const Eigen::Matrix<T,Eigen::Dynamic,1> &theta) const{
    const double p=static_cast<double>(D.size());
    const size_t N=theta.size();
    T tot_sum=0;
    for(int i=0; i<p; i++){
      T tvar=D[i];
      for (int k = 0; k < N; k++) {
        tvar += theta[k] * D[i] * std::pow(D[i], 2 - k);
      }
      tot_sum += log(tvar) + (quh[i] * quh[i]) / tvar;
    }
    return -0.5 * (tot_sum)-0.5 * p * log(2 * M_PI);
  }
};





//' evd_dnorm_grad_stan
//' 
//' This is an attempt to use stan's AD features to calculate a gradient 
//' for the RSSp likelihood
//' 
//' @template RSSp_stat
//' @export
//[[Rcpp::export]]
Eigen::ArrayXd evd_dnorm_grad_stan(const MapA par,const MapA D, const MapA quh){
  //  Rcpp::Rcout<<"Function started!"<<std::endl;
  // const double varu=par(0);
  // const double a=par(1);
  // Rcpp::Rcout<<"sigu is :"<<sigu<<std::endl;
  double fx=0;
  Eigen::Matrix <double,Eigen::Dynamic,1> grad_fx;
  Eigen::Matrix <double,Eigen::Dynamic,1> param=par;
  // param(0)=varu;
  // param(1)=a;
  // Rcpp::Rcout<<"param:"<<param<<std::endl;
  evd_dens f(D,quh);
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
//' @export
//[[Rcpp::export]]
Eigen::MatrixXd evd_dnorm_hess_stan(const MapA par, const MapA D,
                                    const MapA quh) {
  //  Rcpp::Rcout<<"Function started!"<<std::endl;
  const double varu = par(0);
  const double a = par(1);
  // Rcpp::Rcout<<"sigu is :"<<sigu<<std::endl;
  double fx = 0;
  Eigen::Matrix<double, Eigen::Dynamic, 1> grad_fx;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H_fx;

  Eigen::Matrix<double, Eigen::Dynamic, 1> param = par;
  // param(0)=varu;
  // param(1)=a;
  // Rcpp::Rcout<<"param:"<<param<<std::endl;
  evd_dens f(D, quh);
  //  Rcpp::Rcout<<"Struct constructed!"<<std::endl;
  stan::math::hessian(f, param, fx, grad_fx, H_fx);
  //  Rcpp::Rcout<<"Gradient computed!"<<std::endl;

  return (H_fx);
}




template <int I>
struct pve_D_impl {
  template <int M>
  static auto run(const Eigen::Array<double, M, 1> &short_vec,
                  const Eigen::ArrayXd &long_vec)
      -> decltype(short_vec(I) * long_vec.pow(2 - I) +
                  pve_D_impl<I - 1>::run(short_vec, long_vec)) {
    static_assert(I>0,"can't have negative size!");
    static_assert(M>0,"can't have negative size!");
    return short_vec(I) * long_vec.pow(2 - I) +
           pve_D_impl<I - 1>::run(short_vec, long_vec);
  }
};

template <>
struct pve_D_impl<0> {
  template <int M>
  static auto run(const Eigen::Array<double, M, 1> &short_vec,
                  const Eigen::ArrayXd &long_vec)
      -> decltype(0 * long_vec) {
    static_assert(M>0,"can't have negative size!");
    return 0 * long_vec;
  }
};

template <int N>
auto pve_D(const Eigen::Array<double, N, 1> &short_vec,
           const Eigen::ArrayXd &long_vec)
    -> decltype(pve_D_impl<N - 1>::run(short_vec, long_vec)) {
  return pve_D_impl<N - 1>::run(short_vec, long_vec);
}

// template <>
// auto pve_D<1>(const Eigen::Array<double, 1, 1> &short_vec,
//               const Eigen::ArrayXd &long_vec)
//     -> decltype(Eigen::ArrayXd::Zero(long_vec.size())) {
//   return Eigen::ArrayXd::Zero(long_vec.size());
// }

// template<int N,int IB=0>
// Eigen::ArrayXd confound_D(const Eigen::Array<double,N,1> &cvec, const MapA & D){
//   static_assert(IB<=N,"IB must be less than N in confound_D");
//   Eigen::ArrayXd return_vec = Eigen::ArrayXd::Zero(D.size());
//   for(int i=IB; i<N; i++){
//     return_vec += cvec(i)*D.pow(2-i);
//   }
//   return return_vec;
// }


template<int N>
double t_estimate_pve(const Eigen::Array<double,N,1> &cvec, const MapA D,const MapA quh,const int sample_size){

  const size_t p(D.size());
  Eigen::ArrayXd sigma_hat(p);
  Eigen::ArrayXd mu_hat(p);
  sigma_hat = (D * (D + pve_D<N>(cvec, D)).inverse() * D + 1 / cvec(0))
            .inverse();
  mu_hat=sigma_hat*D*(D +pve_D<N>(cvec,D)).inverse()*quh;
  return((D*sigma_hat+mu_hat.square()*D).sum()/sample_size);
}


//' compute pve with confounding
//' @export
//[[Rcpp::export]]
double estimate_pve(const MapA cvec,  const MapA D,const MapA quh,int sample_size){
  const int num_param=cvec.size();
  switch(num_param){
  case 1:{
    Eigen::Array<double,1,1> tcvec(cvec);
    return(t_estimate_pve(tcvec,D,quh,sample_size));
  }
  case 2:{
    Eigen::Array<double,2,1> tcvec(cvec);
    return(t_estimate_pve(tcvec,D,quh,sample_size));
  }
  case 3:{
    Eigen::Array<double,3,1> tcvec(cvec);
    return(t_estimate_pve(tcvec,D,quh,sample_size));
  }
  case 4:{
    Eigen::Array<double,4,1> tcvec(cvec);
    return(t_estimate_pve(tcvec,D,quh,sample_size));
  }
  default: { Rcpp::stop("too many confounding terms!"); }
  }
}
