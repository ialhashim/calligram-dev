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
#include "Constraint.h"
#include <cassert>
///////////////////////////////////////////////////////////////////////////////
#define SHAPEOP_INNER_ITERATIONS 4 //TODO: fix this
///////////////////////////////////////////////////////////////////////////////
namespace ShapeOp {
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE Scalar clamp(Scalar v, Scalar vMin, Scalar vMax) {
  Scalar result = v > vMin ? v : vMin;
  return result > vMax ? vMax : result;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE Constraint::Constraint(const std::vector<int> &idI, Scalar weight) :
  idI_(idI),
  weight_(std::sqrt(weight)) {
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE EdgeStrainConstraint::EdgeStrainConstraint(const std::vector<int> &idI,
                                                          Scalar weight,
                                                          const Matrix3X &positions,
                                                          Scalar rangeMin,
                                                          Scalar rangeMax) :
  Constraint(idI, weight),
  rangeMin_(rangeMin),
  rangeMax_(rangeMax) {
  assert(idI.size() == 2);
  Scalar length = (positions.col(idI_[1]) - positions.col(idI_[0])).norm();
  rest_ = 1.0f / length;
  weight_ *= std::sqrt(length);
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void EdgeStrainConstraint::project(const Matrix3X &positions, Matrix3X &projections) const {
  Vector3 edge = positions.col(idI_[1]) - positions.col(idI_[0]);
  Scalar l = edge.norm();
  edge /= l;
  l = clamp(l * rest_, rangeMin_, rangeMax_);
  projections.col(idO_) = weight_ * l * edge;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void EdgeStrainConstraint::addConstraint(std::vector<Triplet> &triplets, int &idO) const {
  idO_ = idO;
  triplets.push_back(Triplet(idO_, idI_[0], -weight_ * rest_));
  triplets.push_back(Triplet(idO_, idI_[1], weight_ * rest_));
  idO += 1;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void EdgeStrainConstraint::setEdgeLength(Scalar length) {
  rest_ = 1.0f / length;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE TriangleStrainConstraint::TriangleStrainConstraint(const std::vector<int> &idI,
                                                                  Scalar weight,
                                                                  const Matrix3X &positions, bool is2D,
                                                                  Scalar rangeMin,
                                                                  Scalar rangeMax) :
  Constraint(idI, weight),
  rangeMin_(rangeMin),
  rangeMax_(rangeMax), is2D(is2D) {
  assert(idI.size() == 3);
  Matrix32 edges, P;
  edges.col(0) = positions.col(idI_[1]) - positions.col(idI_[0]);
  edges.col(1) = positions.col(idI_[2]) - positions.col(idI_[0]);
  P.col(0) = is2D ? Vector3::UnitX() : edges.col(0).normalized();
  P.col(1) = is2D ? Vector3::UnitY() : (edges.col(1) - edges.col(1).dot(P.col(0)) * P.col(0)).normalized();
  rest_ = (P.transpose() * edges).inverse();
  Scalar A = (P.transpose() * edges).determinant() / 2.0f;
  weight_ *= std::sqrt(std::abs(A));
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void TriangleStrainConstraint::project(const Matrix3X &positions, Matrix3X &projections) const {
  Matrix32 edges, P;
  edges.col(0) = (positions.col(idI_[1]) - positions.col(idI_[0]));
  edges.col(1) = (positions.col(idI_[2]) - positions.col(idI_[0]));
  P.col(0) = is2D ? Vector3::UnitX() : edges.col(0).normalized();
  P.col(1) = is2D ? Vector3::UnitY() : (edges.col(1) - edges.col(1).dot(P.col(0)) * P.col(0)).normalized();
  Matrix22 F = P.transpose() * edges * rest_;
  Eigen::JacobiSVD<Matrix22> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Vector2 S = svd.singularValues();
  S(0) = clamp(S(0), rangeMin_, rangeMax_);
  S(1) = clamp(S(1), rangeMin_, rangeMax_);
  if (is2D && svd.matrixU().determinant()*svd.matrixV().determinant() < 0.0f) S(1) = -S(1); //For 2D
  F = svd.matrixU() * S.asDiagonal() * svd.matrixV().transpose();
  projections.block<3, 2>(0, idO_) = (weight_ * P * F);
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void TriangleStrainConstraint::addConstraint(std::vector<Triplet> &triplets, int &idO) const {
  idO_ = idO;
  int n = 2;
  for (int i = 0; i < n; ++i) {
    triplets.push_back(Triplet(idO_ + i, idI_[0], -weight_ * (rest_(0, i) + rest_(1, i))));
    triplets.push_back(Triplet(idO_ + i, idI_[1], weight_ * rest_(0, i)));
    triplets.push_back(Triplet(idO_ + i, idI_[2], weight_ * rest_(1, i)));
  }
  idO += n;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE TetrahedronStrainConstraint::TetrahedronStrainConstraint(const std::vector<int> &idI,
                                                                        Scalar weight,
                                                                        const Matrix3X &positions,
                                                                        Scalar rangeMin,
                                                                        Scalar rangeMax) :
  Constraint(idI, weight),
  rangeMin_(rangeMin),
  rangeMax_(rangeMax) {
  assert(idI.size() == 4);
  Matrix33 edges;
  for (int i = 0; i < 3; ++i) edges.col(i) = positions.col(idI_[i + 1]) - positions.col(idI_[0]);
  rest_ = edges.inverse();
  Scalar V = (edges).determinant() / 6.0f;
  weight_ *= std::sqrt(std::abs(V));
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void TetrahedronStrainConstraint::project(const Matrix3X &positions, Matrix3X &projections) const {
  Matrix33 edges;
  for (int i = 0; i < 3; ++i) edges.col(i) = positions.col(idI_[i + 1]) - positions.col(idI_[0]);
  Matrix33 F = edges * rest_;
  Eigen::JacobiSVD<Matrix33> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Vector3 S = svd.singularValues();
  S(0) = clamp(S(0), rangeMin_, rangeMax_);
  S(1) = clamp(S(1), rangeMin_, rangeMax_);
  S(2) = clamp(S(2), rangeMin_, rangeMax_);
  if (svd.matrixU().determinant()*svd.matrixV().determinant() < 0.0f) S(2) = -S(2);
  F = svd.matrixU() * S.asDiagonal() * svd.matrixV().transpose();
  projections.block<3, 3>(0, idO_) = weight_ * F;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void TetrahedronStrainConstraint::addConstraint(std::vector<Triplet> &triplets, int &idO) const {
  idO_ = idO;
  int n = 3;
  for (int i = 0; i < n; ++i) {
    triplets.push_back(Triplet(idO_ + i, idI_[0], -weight_ * (rest_(0, i) + rest_(1, i) + rest_(2, i))));
    triplets.push_back(Triplet(idO_ + i, idI_[1], weight_ * rest_(0, i)));
    triplets.push_back(Triplet(idO_ + i, idI_[2], weight_ * rest_(1, i)));
    triplets.push_back(Triplet(idO_ + i, idI_[3], weight_ * rest_(2, i)));
  }
  idO += n;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE AreaConstraint::AreaConstraint(const std::vector<int> &idI,
                                              Scalar weight,
                                              const Matrix3X &positions,
                                              Scalar rangeMin,
                                              Scalar rangeMax) :
  Constraint(idI, weight),
  rangeMin_(rangeMin),
  rangeMax_(rangeMax) {
  assert(idI.size() == 3);
  Matrix32 edges, P;
  edges.col(0) = positions.col(idI_[1]) - positions.col(idI_[0]);
  edges.col(1) = positions.col(idI_[2]) - positions.col(idI_[0]);
  P.col(0) = edges.col(0).normalized();
  P.col(1) = (edges.col(1) - edges.col(1).dot(P.col(0)) * P.col(0)).normalized();
  rest_ = (P.transpose() * edges).inverse();
  Scalar A = (P.transpose() * edges).determinant() / 2.0f;
  weight_ *= std::sqrt(std::abs(A));
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void AreaConstraint::project(const Matrix3X &positions, Matrix3X &projections) const {
  Matrix32 edges, P;
  edges.col(0) = (positions.col(idI_[1]) - positions.col(idI_[0]));
  edges.col(1) = (positions.col(idI_[2]) - positions.col(idI_[0]));
  P.col(0) = edges.col(0).normalized();
  P.col(1) = (edges.col(1) - edges.col(1).dot(P.col(0)) * P.col(0)).normalized();
  Matrix22 F = P.transpose() * edges * rest_;
  Eigen::JacobiSVD<Matrix22> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Vector2 S = svd.singularValues();
  Vector2 d(0.0f, 0.0f);
  for (int i = 0; i < SHAPEOP_INNER_ITERATIONS; ++i) {
    Scalar v = S(0) * S(1);
    Scalar f = v - clamp(v, rangeMin_, rangeMax_);
    Vector2 g(S(1), S(0));
    d = -((f - g.dot(d)) / g.dot(g)) * g;
    S = svd.singularValues() + d;
  }
  F = svd.matrixU() * S.asDiagonal() * svd.matrixV().transpose();
  projections.block<3, 2>(0, idO_) = (weight_ * P * F);
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void AreaConstraint::addConstraint(std::vector<Triplet> &triplets, int &idO) const {
  idO_ = idO;
  int n = 2;
  for (int i = 0; i < n; ++i) {
    triplets.push_back(Triplet(idO_ + i, idI_[0], -weight_ * (rest_(0, i) + rest_(1, i))));
    triplets.push_back(Triplet(idO_ + i, idI_[1], weight_ * rest_(0, i)));
    triplets.push_back(Triplet(idO_ + i, idI_[2], weight_ * rest_(1, i)));
  }
  idO += n;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE VolumeConstraint::VolumeConstraint(const std::vector<int> &idI,
                                                  Scalar weight,
                                                  const Matrix3X &positions,
                                                  Scalar rangeMin,
                                                  Scalar rangeMax) :
  Constraint(idI, weight),
  rangeMin_(rangeMin),
  rangeMax_(rangeMax) {
  assert(idI_.size() == 4);
  Matrix33 edges;
  for (int i = 0; i < 3; ++i) edges.col(i) = positions.col(idI_[i + 1]) - positions.col(idI_[0]);
  rest_ = edges.inverse();
  Scalar V = (edges).determinant() / 6.0f;
  weight_ *= std::sqrt(std::abs(V));
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void VolumeConstraint::project(const Matrix3X &positions, Matrix3X &projections) const {
  Matrix33 edges;
  for (int i = 0; i < 3; ++i) edges.col(i) = positions.col(idI_[i + 1]) - positions.col(idI_[0]);
  Matrix33 F = edges * rest_;
  Eigen::JacobiSVD<Matrix33> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Vector3 S = svd.singularValues();
  Vector3 d(0.0f, 0.0f, 0.0f);
  for (int i = 0; i < SHAPEOP_INNER_ITERATIONS; ++i) {
    Scalar v = S(0) * S(1) * S(2);
    Scalar f = v - clamp(v, rangeMin_, rangeMax_);
    Vector3 g(S(1)*S(2), S(0)*S(2), S(0)*S(1));
    d = -((f - g.dot(d)) / g.dot(g)) * g;
    S = svd.singularValues() + d;
  }
  if (svd.matrixU().determinant()*svd.matrixV().determinant() < 0.0f) S(2) = -S(2);
  F = svd.matrixU() * S.asDiagonal() * svd.matrixV().transpose();
  projections.block<3, 3>(0, idO_) = weight_ * F;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void VolumeConstraint::addConstraint(std::vector<Triplet> &triplets, int &idO) const {
  idO_ = idO;
  int n = 3;
  for (int i = 0; i < n; ++i) {
    triplets.push_back(Triplet(idO_ + i, idI_[0], -weight_ * (rest_(0, i) + rest_(1, i) + rest_(2, i))));
    triplets.push_back(Triplet(idO_ + i, idI_[1], weight_ * rest_(0, i)));
    triplets.push_back(Triplet(idO_ + i, idI_[2], weight_ * rest_(1, i)));
    triplets.push_back(Triplet(idO_ + i, idI_[3], weight_ * rest_(2, i)));
  }
  idO += n;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE BendingConstraint::BendingConstraint(const std::vector<int> &idI,
                                                    Scalar weight,
                                                    const Matrix3X &positions,
                                                    Scalar rangeMin,
                                                    Scalar rangeMax) :
  Constraint(idI, weight),
  rangeMin_(rangeMin),
  rangeMax_(rangeMax) {
  Matrix3X p(3, idI.size());
  for (int i = 0; i < idI_.size(); ++i) p.col(i) = positions.col(idI_[i]);
  Scalar l01 = (p.col(0) - p.col(1)).norm();
  Scalar l02 = (p.col(0) - p.col(2)).norm();
  Scalar l12 = (p.col(1) - p.col(2)).norm();
  Scalar r0 = 0.5 * (l01 + l02 + l12);
  Scalar A0 = std::sqrt(r0 * (r0 - l01) * (r0 - l02) * (r0 - l12));
  Scalar l03 = (p.col(0) - p.col(3)).norm();
  Scalar l13 = (p.col(1) - p.col(3)).norm();
  Scalar r1 = 0.5 * (l01 + l03 + l13);
  Scalar A1 = std::sqrt(r1 * (r1 - l01) * (r1 - l03) * (r1 - l13));
  weight_ *= std::sqrt(3.0 / (A0 + A1));
  Scalar cot02 = ((l01 * l01) - (l02 * l02) + (l12 * l12)) / (4.0 * A0);
  Scalar cot12 = ((l01 * l01) + (l02 * l02) - (l12 * l12)) / (4.0 * A0);
  Scalar cot03 = ((l01 * l01) - (l03 * l03) + (l13 * l13)) / (4.0 * A1);
  Scalar cot13 = ((l01 * l01) + (l03 * l03) - (l13 * l13)) / (4.0 * A1);
  w_ =  Vector4::Zero();
  w_(0) = cot02 + cot03;
  w_(1) = cot12 + cot13;
  w_(2) = -(cot02 + cot12);
  w_(3) = -(cot03 + cot13);
  n_ = (p * w_).norm();
//    int n = static_cast<int>(idI.size());
//    w_ =  VectorX::Zero(n);
//    VectorX A(n-1);
//    for(int i=1; i<n; ++i) {
//        int id1 = i;
//        int id2 = i+1;
//        if(id2 == n) id2 = 1;
//        Scalar l1 = (p.col(0) - p.col(id1)).norm();
//        Scalar l2 = (p.col(0) - p.col(id2)).norm();
//        Scalar l3 = (p.col(id1) - p.col(id2)).norm();
//        Scalar r = 0.5*(l1 + l2 + l3);
//        A[i-1] = std::sqrt(r*(r-l1)*(r-l2)*(r-l3));
//        if(A[i-1] > 1e-6) {
//            Scalar cot1 = (-(l1*l1)+(l2*l2)+(l3*l3))/(4.0*A[i-1]);
//            Scalar cot2 = ((l1*l1)-(l2*l2)+(l3*l3))/(4.0*A[i-1]);
//            w_[id1] += 0.5*cot1;
//            w_[id2] += 0.5*cot2;
//        }
//    }
//    for(int i=1; i<n; ++i) w_[0] -= w_[i];
//    weight_ *= std::sqrt(3.0/A.sum());
//    n_ = (p*w_).norm();
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void BendingConstraint::project(const Matrix3X &positions, Matrix3X &projections) const {
  Vector3 e = Vector3::Zero();
  if (n_ > 1e-6) {
    for (int i = 0; i < idI_.size(); ++i)
      e += w_(i) * positions.col(idI_[i]);
    Scalar l = e.norm();
    if (l > 1e-6) {
      e /= l;
      l = n_ * clamp(l / n_, rangeMin_, rangeMax_);
      e *= l;
    }
  }
  projections.col(idO_) = weight_ * e;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void BendingConstraint::addConstraint(std::vector<Triplet> &triplets, int &idO) const {
  idO_ = idO;
  for (int i = 0; i < idI_.size(); ++i)
    triplets.push_back(Triplet(idO_, idI_[i], weight_ * w_(i)));
  idO += 1;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE ClosenessConstraint::ClosenessConstraint(const std::vector<int> &idI,
                                                        Scalar weight,
                                                        const Matrix3X &positions) :
  Constraint(idI, weight) {
  assert(idI.size() == 1);
  rest_ = positions.col(idI_[0]);
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void ClosenessConstraint::project(const Matrix3X & /*positions*/, Matrix3X &projections) const {
  projections.col(idO_) = rest_ * weight_;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void ClosenessConstraint::addConstraint(std::vector<Triplet> &triplets, int &idO) const {
  idO_ = idO;
  triplets.push_back(Triplet(idO_, idI_[0], weight_));
  idO += 1;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void ClosenessConstraint::setPosition(const Vector3 &position) {
  rest_ = position;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE LineConstraint::LineConstraint(const std::vector<int> &idI,
                                              Scalar weight,
                                              const Matrix3X &) :
  Constraint(idI, weight) {
  assert(idI.size() >= 2);
  input = Matrix3X::Zero(3, idI.size());
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void LineConstraint::project(const Matrix3X &positions, Matrix3X &projections) const {
  for (int i = 0; i < static_cast<int>(idI_.size()); ++i) input.col(i) = positions.col(idI_[i]);
  Vector3 mean_vector = input.rowwise().mean();
  input.colwise() -= mean_vector;
  Eigen::JacobiSVD<Matrix3X> jSVD;
  jSVD.compute(input, Eigen::ComputeFullU);
  Matrix33 basis = jSVD.matrixU();
  input = basis.transpose() * input;
  input.row(1) = Eigen::Matrix<Scalar, 1, Eigen::Dynamic>::Constant(input.cols(), 0.0);
  input.row(2) = Eigen::Matrix<Scalar, 1, Eigen::Dynamic>::Constant(input.cols(), 0.0);
  projections.block(0, idO_, 3, input.cols()) = (basis * input) * weight_;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void LineConstraint::addConstraint(std::vector<Triplet> &triplets, int &idO) const {
  idO_ = idO;
  int n_idx = static_cast<int>(idI_.size());
  double coef1 = (1.0 - 1.0 / n_idx) * weight_;
  double coef2 = - weight_ / n_idx;
  for (int i = 0; i < n_idx; ++i) {
    for (int j = 0; j < n_idx; ++j)
      triplets.push_back(Triplet(idO, idI_[j], (i == j ? coef1 : coef2)));
    idO++;
  }
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE PlaneConstraint::PlaneConstraint(const std::vector<int> &idI,
                                                Scalar weight,
                                                const Matrix3X &) :
  Constraint(idI, weight) {
  assert(idI.size() >= 3);
  input = Matrix3X::Zero(3, idI.size());
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void PlaneConstraint::project(const Matrix3X &positions, Matrix3X &projections) const {
  for (int i = 0; i < static_cast<int>(idI_.size()); ++i) input.col(i) = positions.col(idI_[i]);
  Vector3 mean_vector = input.rowwise().mean();
  input.colwise() -= mean_vector;
  Eigen::JacobiSVD<Matrix3X> jSVD;
  jSVD.compute(input, Eigen::ComputeFullU);
  Matrix33 basis = jSVD.matrixU();
  input = basis.transpose() * input;
  input.row(2) = Eigen::Matrix<Scalar, 1, Eigen::Dynamic>::Constant(input.cols(), 0.0);
  projections.block(0, idO_, 3, input.cols()) = (basis * input) * weight_;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void PlaneConstraint::addConstraint(std::vector<Triplet> &triplets, int &idO) const {
  idO_ = idO;
  int n_idx = static_cast<int>(idI_.size());
  double coef1 = (1.0 - 1.0 / n_idx) * weight_;
  double coef2 = - weight_ / n_idx;
  for (int i = 0; i < n_idx; ++i) {
    for (int j = 0; j < n_idx; ++j)
      triplets.push_back(Triplet(idO, idI_[j], (i == j ? coef1 : coef2)));
    idO++;
  }
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE CircleConstraint::CircleConstraint(const std::vector<int> &idI,
                                                  Scalar weight,
                                                  const Matrix3X &) :
  Constraint(idI, weight) {
  assert(idI.size() >= 3);
  input = Matrix3X::Zero(3, idI.size());
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void CircleConstraint::project(const Matrix3X &positions, Matrix3X &projections) const {
  for (int i = 0; i < static_cast<int>(idI_.size()); ++i) input.col(i) = positions.col(idI_[i]);
  Vector3 mean_vector = input.rowwise().mean();
  input.colwise() -= mean_vector;
  Eigen::JacobiSVD<Matrix3X> jSVD;
  jSVD.compute(input, Eigen::ComputeFullU);
  Matrix33 basis = jSVD.matrixU();
  input = basis.transpose() * input;
  input.row(2) = Eigen::Matrix<Scalar, 1, Eigen::Dynamic>::Constant(input.cols(), 0.0);
  ////////////// 2D Circle fitting
  double Suu = 0.0;
  double Suv = 0.0;
  double Svv = 0.0;
  double Suuu = 0.0;
  double Suvv = 0.0;
  double Svuu = 0.0;
  double Svvv = 0.0;
  for (int j = 0; j < input.cols(); ++j) {
    double uu = input(0, j) * input(0, j);
    double vv = input(1, j) * input(1, j);
    Suu += uu;
    Svv += vv;
    Suv += input(0, j) * input(1, j);
    Suuu += uu * input(0, j);
    Suvv += input(0, j) * vv;
    Svuu += input(1, j) * uu;
    Svvv += vv * input(1, j);
  }
  Matrix22 A;
  A << Suu, Suv,  Suv, Svv;
  if (std::fabs(A.determinant()) > 1e-5) {
    Vector2 b(0.5 * (Suuu + Suvv), 0.5 * (Svvv + Svuu));
    Vector2 center = A.inverse() * b;
    double radius = std::sqrt(center(0) * center(0) + center(1) * center(1) + (Suu + Svv) / static_cast<double>(input.cols()));
    for (int j = 0; j < input.cols(); ++j) {
      Vector2 d = input.block(0, j, 2, 1) - center;
      d.normalize();
      input.block(0, j, 2, 1) = center + d * radius;
    }
  }
  projections.block(0, idO_, 3, input.cols()) = (basis * input) * weight_;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void CircleConstraint::addConstraint(std::vector<Triplet> &triplets, int &idO) const {
  idO_ = idO;
  int n_idx = static_cast<int>(idI_.size());
  double coef1 = (1.0 - 1.0 / n_idx) * weight_;
  double coef2 = - weight_ / n_idx;
  for (int i = 0; i < n_idx; ++i) {
    for (int j = 0; j < n_idx; ++j)
      triplets.push_back(Triplet(idO, idI_[j], (i == j ? coef1 : coef2)));
    idO++;
  }
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE SphereConstraint::SphereConstraint(const std::vector<int> &idI,
                                                  Scalar weight,
                                                  const Matrix3X &) :
  Constraint(idI, weight) {
  assert(idI.size() >= 4);
  input = Matrix3X::Zero(3, idI.size());
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void SphereConstraint::project(const Matrix3X &positions, Matrix3X &projections) const {
  for (int i = 0; i < static_cast<int>(idI_.size()); ++i) input.col(i) = positions.col(idI_[i]);
  Vector3 mean_vector = input.rowwise().mean();
  input.colwise() -= mean_vector;
  ////////////// 3D Sphere fitting
  double Suu = 0.0;
  double Suv = 0.0;
  double Suw = 0.0;
  double Svv = 0.0;
  double Svw = 0.0;
  double Sww = 0.0;
  double Suuu = 0.0;
  double Suvv = 0.0;
  double Suww = 0.0;
  double Svuu = 0.0;
  double Svvv = 0.0;
  double Svww = 0.0;
  double Swuu = 0.0;
  double Swvv = 0.0;
  double Swww = 0.0;
  for (int j = 0; j < input.cols(); ++j) {
    double uu = input(0, j) * input(0, j);
    double vv = input(1, j) * input(1, j);
    double ww = input(2, j) * input(2, j);
    Suu += uu;
    Svv += vv;
    Sww += ww;
    Suv += input(0, j) * input(1, j);
    Suw += input(0, j) * input(2, j);
    Svw += input(1, j) * input(2, j);
    Suuu += input(0, j) * uu;
    Suvv += input(0, j) * vv;
    Suww += input(0, j) * ww;
    Svuu += input(1, j) * uu;
    Svvv += input(1, j) * vv;
    Svww += input(1, j) * ww;
    Swuu += input(2, j) * uu;
    Swvv += input(2, j) * vv;
    Swww += input(2, j) * ww;
  }
  Matrix33 A;
  A << Suu, Suv, Suw,  Suv, Svv, Svw, Suw, Svw, Sww;
  if (std::fabs(A.determinant()) > 1e-5) {
    Vector3 b(0.5 * (Suuu + Suvv + Suww), 0.5 * (Svuu + Svvv + Svww), 0.5 * (Swuu + Swvv + Swww));
    Vector3 center = A.inverse() * b;
    double radius = std::sqrt(center(0) * center(0) + center(1) * center(1) + center(2) * center(2) + (Suu + Svv + Sww) / static_cast<double>(input.cols()));
    for (int j = 0; j < input.cols(); ++j) {
      Vector3 d = input.col(j) - center;
      d.normalize();
      input.col(j) = center + d * radius;
    }
  }
  projections.block(0, idO_, 3, input.cols()) = input * weight_;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void SphereConstraint::addConstraint(std::vector<Triplet> &triplets, int &idO) const {
  idO_ = idO;
  int n_idx = static_cast<int>(idI_.size());
  double coef1 = (1.0 - 1.0 / n_idx) * weight_;
  double coef2 = - weight_ / n_idx;
  for (int i = 0; i < n_idx; ++i) {
    for (int j = 0; j < n_idx; ++j)
      triplets.push_back(Triplet(idO, idI_[j], (i == j ? coef1 : coef2)));
    idO++;
  }
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE UniformLaplacianConstraint::UniformLaplacianConstraint(const std::vector<int> &idI,
                                                                      Scalar weight,
                                                                      const Matrix3X &positions,
                                                                      bool displacement_lap) :
  Constraint(idI, weight) {
  weighted_rest_.setZero();
  if (displacement_lap) {
    int n_idx = static_cast<int>(idI.size());
    for (int i = 1; i < n_idx; ++i) weighted_rest_ += positions.col(idI[i]);
    weighted_rest_ /= double(n_idx - 1);
    weighted_rest_ -= positions.col(idI[0]);
    weighted_rest_ *= weight_;
  }
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void UniformLaplacianConstraint::project(const Matrix3X & /*positions*/, Matrix3X &projections) const {
  projections.col(idO_) = weighted_rest_;
}
///////////////////////////////////////////////////////////////////////////////
SHAPEOP_INLINE void UniformLaplacianConstraint::addConstraint(std::vector<Triplet> &triplets, int &idO) const {
  idO_ = idO;
  triplets.push_back(Triplet(idO_, idI_[0], -weight_));
  double c = weight_ / (idI_.size() - 1);
  for (int i = 1; i < static_cast<int>(idI_.size()); ++i)
    triplets.push_back(Triplet(idO_, idI_[i], c));
  idO ++;
}
///////////////////////////////////////////////////////////////////////////////
} // namespace ShapeOp
///////////////////////////////////////////////////////////////////////////////
