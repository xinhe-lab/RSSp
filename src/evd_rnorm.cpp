#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
typedef Eigen::Map<Eigen::ArrayXd> mapa;
typedef Eigen::Map<Eigen::MatrixXd> mapmat;
typedef Eigen::Map<Eigen::VectorXd> mapvec;

// 

//[[Rcpp::export]]
Eigen::MatrixXd simuh_dir_cpp(double sigu,  double bias, int nreps, Eigen::MatrixXd &Q, Eigen::ArrayXd &D, Rcpp::StringVector fgeneids, Eigen::MatrixXd usim){
  size_t p = D.size();
  if(p!=Q.rows()){
    Rcpp::stop("Q and D do not match in size");
  }
  if(usim.rows()!=nreps){
    Rcpp::stop("usim must have nreps cols");
  }
  Eigen::VectorXd nd=(sigu*sigu*D.square()+D+bias).sqrt();
  return((usim*(Q*(nd.asDiagonal())*Q.transpose())).transpose());
}


//[[Rcpp::export]]
Rcpp::NumericMatrix block_mat_mul(Rcpp::ListOf<Rcpp::NumericMatrix> &mat_l, Rcpp::NumericMatrix &ymat,bool transpose_mat_l = false){
  Rcpp::NumericMatrix retmat(ymat.nrow(),ymat.ncol());
  Eigen::Map<Eigen::MatrixXd> tretmat(&retmat[0],retmat.nrow(),retmat.ncol());
  Eigen::Map<Eigen::MatrixXd> tymat(&ymat[0],ymat.nrow(),ymat.ncol());
  const size_t list_size =mat_l.size();
  size_t offput=0;
  Rcpp::NumericMatrix* tptr;
  for(size_t i=0;i<list_size;i++){
    Eigen::Map<Eigen::MatrixXd> ttmat(&((mat_l[i])[0]),(mat_l[i]).nrow(),(mat_l[i]).ncol());
    size_t ttmat_rows=ttmat.rows();
    size_t ttmat_cols=ttmat.cols();
    if(transpose_mat_l){
      tretmat.block(offput,0,ttmat_rows,tretmat.cols())=ttmat.transpose()*tymat.block(offput,0,ttmat_rows,tretmat.cols());
    }else{
        tretmat.block(offput,0,ttmat_rows,tretmat.cols())=ttmat*tymat.block(offput,0,ttmat_rows,tretmat.cols());
    }
    offput+=ttmat_rows;
  }
  return(retmat);
}
