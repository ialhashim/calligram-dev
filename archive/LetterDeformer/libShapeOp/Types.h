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
#ifndef TYPES_H
#define TYPES_H
///////////////////////////////////////////////////////////////////////////////
#include "Common.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
///////////////////////////////////////////////////////////////////////////////
/** \file
This file redefines EIGEN types using the scalar type ::ShapeOpScalar defined in Common.h.*/
///////////////////////////////////////////////////////////////////////////////
//TODO: For windows 32 bit we may need Eigen::DontAlign
#ifdef SHAPEOP_DONT_ALIGN
#define SHAPEOP_ALIGNMENT Eigen::DontAlign
#else
#define SHAPEOP_ALIGNMENT Eigen::AutoAlign
#endif
///////////////////////////////////////////////////////////////////////////////
namespace ShapeOp {
typedef ShapeOpScalar Scalar;
//Dense
template < int Rows, int Cols, int Options = (Eigen::ColMajor | SHAPEOP_ALIGNMENT) >
using MatrixT = Eigen::Matrix<Scalar, Rows, Cols, Options>;
typedef MatrixT<2, 1> Vector2;
typedef MatrixT<2, 2> Matrix22;
typedef MatrixT<2, 3> Matrix23;
typedef MatrixT<3, 1> Vector3;
typedef MatrixT<3, 2> Matrix32;
typedef MatrixT<3, 3> Matrix33;
typedef MatrixT<3, 4> Matrix34;
typedef MatrixT<4, 1> Vector4;
typedef MatrixT<4, 4> Matrix44;
typedef MatrixT<3, Eigen::Dynamic> Matrix3X;
typedef MatrixT<Eigen::Dynamic, 3> MatrixX3;
typedef MatrixT<Eigen::Dynamic, 1> VectorX;
typedef MatrixT<Eigen::Dynamic, Eigen::Dynamic> MatrixXX;
//Sparse
template<int Options = Eigen::ColMajor>
using SparseMatrixT = Eigen::SparseMatrix<Scalar, Options>;
typedef SparseMatrixT<> SparseMatrix;
typedef Eigen::Triplet<Scalar> Triplet;
///////////////////////////////////////////////////////////////////////////////
} // namespace ShapeOp
///////////////////////////////////////////////////////////////////////////////
#endif // TYPES_H
///////////////////////////////////////////////////////////////////////////////
