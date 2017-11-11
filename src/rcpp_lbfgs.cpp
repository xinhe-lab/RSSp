// #include <RcppEigen.h>
// #include <boost/math/tools/promotion.hpp>
// #include <stan/math/prim/scal.hpp>
// #include <stan/math.hpp>
// #include <stan/math/mix/mat.hpp>
// 
// using namespace Rcpp;
// typedef Eigen::Map<Eigen::ArrayXd> Mapa;
// typedef Eigen::Map<Eigen::MatrixXd> mapmat;
// typedef Eigen::Map<Eigen::VectorXd> mapvec;
// 
// 
// using namespace Numer;
// typedef const Eigen::Ref<const Eigen::VectorXd> Constvec;
// typedef Eigen::Ref<Eigen::VectorXd> Refvec;
// typedef Eigen::Map<Eigen::MatrixXd> MapMat;
// typedef Eigen::Map<Eigen::VectorXd> MapVec;
// 
// class RSSp_mvd_nc: public MFuncGrad
// {
// private:
//   const Mapa quh;
//   const Mapa d;
//   const int p;
// public:
//   RSSp_mvd_nc(const Mapa quh_, const Mapa d_) : quh(quh_), d(d_) ,p(d.size()){}
//   
//   double f_grad(Constvec& par, Refvec grad)
//   {
//     const double varu=par(0);
//     const double a=0;
//     
//     const double tsum = (-0.5*(((d.square()*varu+d+a).log())+(quh.square()*(1/(d.square()*varu+d+a))))).sum();
//     grad.noalias() = (-(dvec.square()*(a + dvec + dvec.square()*varu - quh.square()))/(2*(a + dvec + dvec.square()*varu).square())).sum();
//     grad(0)=sgrad
//       //    const double  agrad = (-(a+dvec.square()*varu+dvec-quh.square())/(2*(a+dvec.square()*varu+dvec).square())).sum();
//       return(tsum);
//   }
// };
// 
// 
// 
// class RSSp_mvd_c: public MFuncGrad
// {
// private:
//   const Mapa quh;
//   const Mapa d;
//   const int p;
// public:
//   RSSp_mvd_c(const Mapa quh_, const Mapa d_) : quh(quh_), d(d_) ,p(d.size()){}
//   
//   double f_grad(Constvec& par, Refvec grad)
//   {
//     const double varu=par(0);
//     const double a=par(1)*varu;
//     
//     const double tsum = (-0.5*(((d.square()*varu+d+a).log())+(quh.square()*(1/(d.square()*varu+d+a))))).sum();
//     grad(0) = (-(d.square()*(a + d + d.square()*varu - quh.square()))/(2*(a + d + d.square()*varu).square())).sum();
//     grad(1) = (-(a+dvec.square()*varu+dvec-quh.square())/(2*(a+dvec.square()*varu+dvec).square())).sum();
//       return(tsum);
//   }
// };
// // 
// // Rcpp::DataFrame optim_RSSP(const Eigen::Map<Eigen::ArrayXd> par,const Eigen::Map<Eigen::ArrayXd> dvec, const Eigen::Map<Eigen::ArrayXd> quh))
// // 
// 
