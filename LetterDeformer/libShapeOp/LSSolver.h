///////////////////////////////////////////////////////////////////////////////
// This file is part of ShapeOp, a lightweight C++ library
// for static and dynamic geometry processing.
//
// Copyright (C) 2014 Sofien Bouaziz <sofien.bouaziz@gmail.com>
// Copyright (C) 2014 LGG EPFL
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
///////////////////////////////////////////////////////////////////////////////
#ifndef LSSOLVER_H
#define LSSOLVER_H
///////////////////////////////////////////////////////////////////////////////
#include "Types.h"
///////////////////////////////////////////////////////////////////////////////
/** \file
This file contains all the linear system solvers of the ShapeOp library.*/
///////////////////////////////////////////////////////////////////////////////
namespace ShapeOp {
///////////////////////////////////////////////////////////////////////////////
/** \brief Base class of any sparse linear system solver. This class defines the main functionalities of the ShapeOp sparse linear system solvers (Ax = b).*/
class SHAPEOP_API LSSolver {
 public:
  virtual ~LSSolver() {};
  /** \brief Initialize the linear system solver using the sparse matrix A.*/
  virtual void initialize(const SparseMatrix &A) = 0;
  /** \brief Solve the linear system Ax = b.*/
  virtual VectorX solve(const VectorX &b) const = 0;
};
///////////////////////////////////////////////////////////////////////////////
/** \brief Sparse linear system solver based on Cholesky. This class implements a sparse linear system solver based on the Cholesky LDL^T algorithm from Eigen.*/
class SHAPEOP_API SimplicialLDLTSolver : public LSSolver {
 public:
  virtual ~SimplicialLDLTSolver() {};
  /** \brief Prefactorize the sparse matrix (A = LDL^T).*/
  virtual void initialize(const SparseMatrix &A) override final;
  /** \brief Solve the linear system by applying twice backsubstitution.*/
  virtual VectorX solve(const VectorX &b) const override final;
 private:
  Eigen::SimplicialLDLT<SparseMatrix> solver_;
};
///////////////////////////////////////////////////////////////////////////////
} // namespace ShapeOp
///////////////////////////////////////////////////////////////////////////////
#ifdef SHAPEOP_HEADER_ONLY
#include "LSSolver.cpp"
#endif
///////////////////////////////////////////////////////////////////////////////
#endif // LSSOLVER_H
///////////////////////////////////////////////////////////////////////////////
