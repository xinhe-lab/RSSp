// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(StanHeaders)]]
#include <RSSp.h>

#include <stan/math/prim/scal.hpp>
#include <stan/math.hpp>
#include <stan/math/mix/mat.hpp>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(unwindProtect)]]
//[[Rcpp::depends(BH)]]
using namespace Rcpp;

typedef Eigen::Map<Eigen::ArrayXd>  MapA;
typedef Eigen::Map<Eigen::MatrixXd> MapMat;
typedef Eigen::Map<Eigen::VectorXd> MapVec;
typedef Eigen::Ref<Eigen::VectorXd> RefVec;
typedef Eigen::Ref<Eigen::MatrixXd> RefMat;
typedef Eigen::Ref<Eigen::ArrayXd>  RefA;


// [[Rcpp::plugins(cpp14)]]




template<typename T>
T t_evd_dnorm(const Eigen::Array<T,Eigen::Dynamic,1> cvec ,const MapA D, const MapA quh){
  using std::log;
  const int N=cvec.size();
  const T p=static_cast<T>(D.size());
  T tot_sum=0;
  for(int i=0; i<p; i++){
    T tvar=D[i];
    for (int k = 0; k < N; k++) {
      tvar += cvec[k] * std::pow(D[i], 2 - k);
    }
    tot_sum += log(tvar) + (quh[i] * quh[i]) / tvar;
  }
  return -(-0.5 * (tot_sum) -0.5 * p * log(2 * M_PI));
}

template<typename T>
T t_evd_dnorm(const T cvec, const MapA D, const MapA quh){
  using std::log;
  //  Rcpp::Rcout<<"1"<<std::endl;
  const T p=static_cast<T>(D.size());
  const T tsum = ((cvec*D.square()+D).log()+ (quh.square()*(1/(cvec*D.square()+D)))).sum();
  return -(-0.5*(tsum)-0.5*p*log(2*M_PI));
}

template<typename T>
T t_evd_dnorm(const Eigen::Array<T, 2, 1> cvec, const MapA D,
                      const MapA quh) {
  using std::log;
  //  Rcpp::Rcout<<"2"<<std::endl;
  const T p = static_cast<T>(D.size());
  const T tsum =
      ((D + cvec[0] * D.square() + cvec[1] * D).log() +
       (quh * (1 / (D + cvec[0] * D.square() + cvec[1] * D)) * quh))
          .sum();
  return -(-0.5 * (tsum)-0.5 * p * log(2 * M_PI));
}

template<typename T>
T t_evd_dnorm(const Eigen::Array<T,3,1> cvec, const MapA D, const MapA quh){
  using std::log;
  //  Rcpp::Rcout<<"3"<<std::endl;
  const T p=static_cast<T>(D.size());
  const T tsum =          (D+cvec[0]*D.square()+cvec[1]*D+cvec[2]).log().sum();
  const T tprod = (quh*(1/(D+cvec[0]*D.square()+cvec[1]*D+cvec[2]))*quh).sum();
  return -(-0.5*(tsum+tprod)-0.5*p*log(2*M_PI));
}



template<typename T>
T t_evd_dnorm(const Eigen::Array<T,4,1> cvec, const MapA D, const MapA quh){
  using std::log;
  //  Rcpp::Rcout<<"4"<<std::endl;
  const T p=static_cast<T>(D.size());
  const T tsum =          (D+cvec[0]*D.square()+cvec[1]*D+cvec[2]+cvec[3]*D.pow(-1)).log().sum();
  const T tprod = (quh*(1/(D+cvec[0]*D.square()+cvec[1]*D+cvec[2]+cvec[3]*D.pow(-1)))*quh).sum();
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
    double tcvec=par[0];
    return(t_evd_dnorm(tcvec,D,quh));
  }
  case 2:{
    Eigen::Array<double,2,1> tcvec(par);
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
    Eigen::Array<double,Eigen::Dynamic,1> tcvec(par);
    return(t_evd_dnorm(tcvec,D,quh));
  }
  default: { Rcpp::stop("too many confounding terms!"); }
  }
}

//Borrowed heavily from https://arxiv.org/pdf/1509.07164.pdf




  // int num_terms=par.size();
  // evd_dens f(D,quh,num_terms);
  // //  Rcpp::Rcout<<"Struct constructed!"<<std::endl;
  // Eigen::Matrix<double,Eigen::Dynamic,1> param(par);
  // Eigen::Matrix <double,Eigen::Dynamic,1> grad_fx;
  // stan::math::gradient(f,param,fx,grad_fx);

  // return (grad_fx.array());



struct evd_dens {
  const MapA D;
  const MapA quh;
  const size_t p;
  const size_t N;
  evd_dens(const MapA D_,const MapA quh_,const size_t N_): D(D_),quh(quh_),p(quh.size()),N(N_){}
  template <typename T>
  T operator()(const Eigen::Matrix<T,Eigen::Dynamic,1> &par) const{
      Eigen::Array<T,Eigen::Dynamic,1> tcvec(par);
      return(t_evd_dnorm(tcvec,D,quh));
  }
};




//' Gradient for RSSp likelihood
//' This is an attempt to use stan's AD features to optimize the RSSp likelihood
//' @template RSSp_stat
//' @export
//[[Rcpp::export]]
Eigen::ArrayXd evd_dnorm_grad_stan(const MapA par,const MapA D, const MapA quh){
  //  Rcpp::Rcout<<"Function started!"<<std::endl;
  // const double varu=par(0);
  // const double a=par(1);
  // Rcpp::Rcout<<"sigu is :"<<sigu<<std::endl;
  double fx=0;
  
  int num_terms=par.size();
  evd_dens f(D,quh,num_terms);
  //  Rcpp::Rcout<<"Struct constructed!"<<std::endl;
  Eigen::Matrix<double,Eigen::Dynamic,1> param(par);
  Eigen::Matrix <double,Eigen::Dynamic,1> grad_fx;
  stan::math::gradient(f,param,fx,grad_fx);
  
  return (grad_fx.array());
}



//' evd_dnorm_grad_stan
//'
//' This is an attempt to use stan's AD features to calculate a gradient
//' for the RSSp likelihood
//'
//' @template RSSp_stat
//' @export
//[[Rcpp::export]]
double evd_dnorm_stan(const MapA par,const MapA D, const MapA quh){

    double fx=0;

  int num_terms=par.size();
  evd_dens f(D,quh,num_terms);
  //  Rcpp::Rcout<<"Struct constructed!"<<std::endl;
  Eigen::Matrix<double,Eigen::Dynamic,1> param(par);
  Eigen::Matrix <double,Eigen::Dynamic,1> grad_fx;
  stan::math::gradient(f,param,fx,grad_fx);

  return (fx);

}

  //  Rcpp::Rcout<<"Gradient computed!"<<std::endl;

//};


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

  // Rcpp::Rcout<<"sigu is :"<<sigu<<std::endl;
  double fx = 0;
  Eigen::Matrix<double, Eigen::Dynamic, 1> grad_fx;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> H_fx;

  Eigen::Matrix<double, Eigen::Dynamic, 1> param = par;
  // param(0)=varu;
  // param(1)=a;
  // Rcpp::Rcout<<"param:"<<param<<std::endl;
  evd_dens f(D, quh,par.size());
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
//' @param cvec vector with length equal to the number of terms used to fit the model that contains parameter estimates
//' @param D vector of eigenvalues
//' @param quh vector of transformed summary statistics (must have length equal to quh)
//' @param sample_size integer giving sample size of original GWAS
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



template <int I>
struct do_transformation_impl {
  template<int M>
  static auto run(const Eigen::Array<double,M,1> &short_vec, const Eigen::ArrayXd &long_vec)
    -> decltype(short_vec(I)*long_vec.pow(2-I) + do_transformation_impl<I-1>::run(short_vec,long_vec))
    {
      return short_vec(I)*long_vec.pow(2-I) + do_transformation_impl<I-1>::run(short_vec,long_vec);
    }
};

template <>
struct do_transformation_impl<0> {
  template<int M>
  static auto run(const Eigen::Array<double,M,1> &short_vec, const Eigen::ArrayXd &long_vec)
    -> decltype(short_vec(0)*long_vec.pow(2))
    {
      return short_vec(0)*long_vec.pow(2);
    }
};

template<int N>
auto do_transformation(const Eigen::Array<double,N,1> &short_vec, const Eigen::ArrayXd &long_vec)
  -> decltype(do_transformation_impl<N-1>::run(short_vec,long_vec))
  {
    return do_transformation_impl<N-1>::run(short_vec,long_vec);
  }

