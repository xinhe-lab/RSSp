#ifndef RSSP_H
#define RSSP_H
// #pragma GCC diagnostic ignored "-Wignored-attributes" 
// #pragma GCC diagnostic ignored "-Wdeprecated-declarations" 

#include <RcppEigen.h>
typedef Eigen::Map<Eigen::ArrayXd>  MapA;
typedef Eigen::Map<Eigen::MatrixXd> MapMat;
typedef Eigen::Map<Eigen::VectorXd> MapVec;
typedef Eigen::Ref<Eigen::VectorXd> RefVec;
typedef Eigen::Ref<Eigen::MatrixXd> RefMat;
typedef Eigen::Ref<Eigen::ArrayXd>  RefA;


typedef const Eigen::Ref<const Eigen::VectorXd> ConstVec;
typedef const Eigen::Ref<const Eigen::ArrayXd> ConstA;

#endif