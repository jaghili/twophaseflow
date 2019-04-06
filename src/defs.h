#ifndef DEFS_H
#define DEFS_H

#define LONGLINE "-------------------------------------------------"

#include <iostream>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <vector>

typedef double R;
typedef Eigen::VectorXi iVec;
typedef Eigen::MatrixXi iMat;
typedef Eigen::VectorXd Vec;
typedef Eigen::MatrixXd Mat;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> Triplet;

using std::min;
using std::max;

template <typename T>
uint compare(T& A, T& B) { return (A.idx < B.idx); };


// count the nnz of a Eigen::Matrix m
/*
  void countnnz(Mat &m) {
  uint nnz=0;
  // browse columnwise
  for (int j = 0 ; j < m.cols(); j++) {
  for (int i = 0; i < m.rows(); i++) {
  if ( m(i,j) != 0.) { nnz++; }
  }
  }
  std::cout << "nnz=" << nnz << "/" << m.cols()*m.rows()
  << "\t ("<< 100 * ((float)nnz / (float)((m.cols()*m.rows()))) << "%)"
  << std::endl;
  }
*/
#endif
