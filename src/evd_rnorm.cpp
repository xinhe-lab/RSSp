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

// simuh_dir <-  function(sigu,bias,Q,D,fgeneids,seed=NULL,chr=NULL,range_id=NULL,gen_quh=F,ret_d=F){
//   p <- nrow(Q)
//   if(!is.null(seed)){
//     stopifnot(is.integer(seed))
//     set.seed(seed)
//   }
//   nreps <- length(fgeneids[[1]])
//     fgeneids <- as.character(fgeneids[[1]])
//     nd <- sigu^2*D^2+D+bias
// #  A <- Q%*%diag(sqrt(nd))
//     A <- Q %*% (t(Q) * sqrt(pmax(nd, 0)))
//     usim <- matrix(rnorm(n = p*nreps),nrow=nreps,byrow = T)
//     uhmat <- t(usim%*%A)
//     if(gen_quh){
//       quhmat <- crossprod(Q,uhmat)
//     }
//     
// # mcov <- sigu^2*Rsq+R+bias*diag(p)
// # uhmat <- t(mvtnorm::rmvnorm(n = nreps,mean = rep(0,p),sigma = mcov))
//     
//     colnames(uhmat) <- fgeneids
//       if(gen_quh){
//         colnames(quhmat) <- fgeneids
//       }
//       retdf <- as_data_frame(uhmat) %>% mutate(snp_index=as.character(1:n())) %>% gather(fgeneid,uh,-snp_index) %>% mutate(tsigu=sigu,tbias=bias)
//         if(gen_quh){  
//           retdf <- as_data_frame(quhmat) %>% mutate(snp_index=as.character(1:n())) %>% gather(fgeneid,quh,-snp_index) %>% mutate(tsigu=sigu,tbias=bias) %>% inner_join(retdf,by=c("snp_index","fgeneid","tsigu","tbias"))
//         }
//         if(!is.null(range_id)){
//           retdf <- mutate(retdf,range_id=range_id)
//         }
//         if(!is.null(chr)){
//           retdf <- mutate(retdf,chrom=chr)
//         }
//         if(ret_d){
//           retdf <- data_frame(rd=D) %>% mutate(snp_index=as.character(1:n())) %>% inner_join(retdf,by = "snp_index")
//         }
//         return(retdf)
// }