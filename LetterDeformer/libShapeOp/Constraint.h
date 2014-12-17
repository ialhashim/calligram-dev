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
#ifndef CONSTRAINT_H
#define CONSTRAINT_H
///////////////////////////////////////////////////////////////////////////////
#include "Types.h"
///////////////////////////////////////////////////////////////////////////////
/** \file
This file containts all the constraints of the ShapeOp libary.*/
///////////////////////////////////////////////////////////////////////////////
namespace ShapeOp {
///////////////////////////////////////////////////////////////////////////////
/** \brief Base class of any constraints. This class defines the interface of a ShapeOp constraint.*/
class SHAPEOP_API Constraint {
 public:
  /** \brief Constraint constructor.*/
  Constraint(const std::vector<int> &idI, Scalar weight);
  virtual ~Constraint() {}
  /** \brief Find the closest configuration from the input positions that satisfy the constraint.*/
  virtual void project(const Matrix3X &positions, Matrix3X &projections) const = 0;
  /** \brief Add the constraint to the linear system.*/
  virtual void addConstraint(std::vector<Triplet> &triplets, int &idO) const = 0;
  /** \brief Compute the constraint violation for the given input positions.*/
  Scalar error(const Matrix3X &) const { return 0; }
 protected:
  /** \brief ids of the vertices involved in this constraint.*/
  std::vector<int> idI_;
  /** \brief weight for the constraint.*/
  Scalar weight_;
  /** \brief location of this constraint in the linear system.*/
  mutable int idO_;
};
///////////////////////////////////////////////////////////////////////////////
/** \brief Edge strain constraint. See \cite Bouaziz2014 for more details.*/
class SHAPEOP_API EdgeStrainConstraint : public Constraint {
 public:
  /** \brief Constraint constructor.*/
  EdgeStrainConstraint(const std::vector<int> &idI,
                       Scalar weight,
                       const Matrix3X &positions,
                       Scalar rangeMin = 1.0,
                       Scalar rangeMax = 1.0);
  virtual ~EdgeStrainConstraint() {}
  /** \brief Find the closest configuration from the input positions that satisfy the constraint.*/
  virtual void project(const Matrix3X &positions, Matrix3X &projections) const override final;
  /** \brief Add the constraint to the linear system.*/
  virtual void addConstraint(std::vector<Triplet> &triplets, int &idO) const override final;
  /** \brief Set a new edge length.*/
  void setEdgeLength(Scalar length);
 private:
  Scalar rest_;
  Scalar rangeMin_;
  Scalar rangeMax_;
};
///////////////////////////////////////////////////////////////////////////////
/** \brief Triangle strain constraint. See \cite Bouaziz2014 for more details.*/
class SHAPEOP_API TriangleStrainConstraint : public Constraint {
 public:
  /** \brief Constraint constructor.*/
  TriangleStrainConstraint(const std::vector<int> &idI,
                           Scalar weight,
                           const Matrix3X &positions, bool is2D = false,
                           Scalar rangeMin = 1.0,
                           Scalar rangeMax = 1.0);
  virtual ~TriangleStrainConstraint() {}
  /** \brief Find the closest configuration from the input positions that satisfy the constraint.*/
  virtual void project(const Matrix3X &positions, Matrix3X &projections) const override final;
  /** \brief Add the constraint to the linear system.*/
  virtual void addConstraint(std::vector<Triplet> &triplets, int &idO) const override final;
 private:
  Matrix22 rest_;
  Scalar rangeMin_;
  Scalar rangeMax_;
  bool is2D;
};
///////////////////////////////////////////////////////////////////////////////
/** \brief Tetrahedron strain constraint. See \cite Bouaziz2014 for more details.*/
class SHAPEOP_API TetrahedronStrainConstraint : public Constraint {
 public:
  /** \brief Constraint constructor.*/
  TetrahedronStrainConstraint(const std::vector<int> &idI,
                              Scalar weight,
                              const Matrix3X &positions,
                              Scalar rangeMin = 1.0,
                              Scalar rangeMax = 1.0);
  virtual ~TetrahedronStrainConstraint() {}
  /** \brief Find the closest configuration from the input positions that satisfy the constraint.*/
  virtual void project(const Matrix3X &positions, Matrix3X &projections) const override final;
  /** \brief Add the constraint to the linear system.*/
  virtual void addConstraint(std::vector<Triplet> &triplets, int &idO) const override final;
 private:
  Matrix33 rest_;
  Scalar rangeMin_;
  Scalar rangeMax_;
};
///////////////////////////////////////////////////////////////////////////////
/** \brief Area constraint. See \cite Bouaziz2014 for more details.*/
class SHAPEOP_API AreaConstraint : public Constraint {
 public:
  /** \brief Constraint constructor.*/
  AreaConstraint(const std::vector<int> &idI,
                 Scalar weight,
                 const Matrix3X &positions,
                 Scalar rangeMin = 1.0,
                 Scalar rangeMax = 1.0);
  virtual ~AreaConstraint() {}
  /** \brief Find the closest configuration from the input positions that satisfy the constraint.*/
  virtual void project(const Matrix3X &positions, Matrix3X &projections) const override final;
  /** \brief Add the constraint to the linear system.*/
  virtual void addConstraint(std::vector<Triplet> &triplets, int &idO) const override final;
 private:
  Matrix22 rest_;
  Scalar rangeMin_;
  Scalar rangeMax_;
};
///////////////////////////////////////////////////////////////////////////////
/** \brief Volume constraint. See \cite Bouaziz2014 for more details.*/
class SHAPEOP_API VolumeConstraint : public Constraint {
 public:
  /** \brief Constraint constructor.*/
  VolumeConstraint(const std::vector<int> &idI,
                   Scalar weight,
                   const Matrix3X &positions,
                   Scalar rangeMin = 1.0,
                   Scalar rangeMax = 1.0);
  virtual ~VolumeConstraint() {}
  /** \brief Find the closest configuration from the input positions that satisfy the constraint.*/
  virtual void project(const Matrix3X &positions, Matrix3X &projections) const override final;
  /** \brief Add the constraint to the linear system.*/
  virtual void addConstraint(std::vector<Triplet> &triplets, int &idO) const override final;
 private:
  Matrix33 rest_;
  Scalar rangeMin_;
  Scalar rangeMax_;
};
///////////////////////////////////////////////////////////////////////////////
/** \brief Bending constraint. See \cite Bouaziz2014 for more details.*/
class SHAPEOP_API BendingConstraint : public Constraint {
 public:
  /** \brief Constraint constructor.*/
  BendingConstraint(const std::vector<int> &idI,
                    Scalar weight,
                    const Matrix3X &positions,
                    Scalar rangeMin = 1.0,
                    Scalar rangeMax = 1.0);
  virtual ~BendingConstraint() {}
  /** \brief Find the closest configuration from the input positions that satisfy the constraint.*/
  virtual void project(const Matrix3X &positions, Matrix3X &projections) const override final;
  /** \brief Add the constraint to the linear system.*/
  virtual void addConstraint(std::vector<Triplet> &triplets, int &idO) const override final;
 private:
  VectorX w_;
  Scalar n_;
  Scalar rangeMin_;
  Scalar rangeMax_;
};
///////////////////////////////////////////////////////////////////////////////
class SHAPEOP_API ClosenessConstraint : public Constraint {
 public:
  /** \brief Constraint constructor.*/
  ClosenessConstraint(const std::vector<int> &idI,
                      Scalar weight,
                      const Matrix3X &positions);
  virtual ~ClosenessConstraint() {}
  /** \brief Find the closest configuration from the input positions that satisfy the constraint.*/
  virtual void project(const Matrix3X & /*positions*/, Matrix3X &projections) const override final;
  /** \brief Add the constraint to the linear system.*/
  virtual void addConstraint(std::vector<Triplet> &triplets, int &idO) const override final;
  /** \brief Set a new closeness position.*/
  void setPosition(const Vector3 &position);
 private:
  Vector3 rest_;
};
///////////////////////////////////////////////////////////////////////////////
/** \brief Line constraint. See \cite Bouaziz2012 for more details.*/
class SHAPEOP_API LineConstraint : public Constraint {
 public:
  /** \brief Constraint constructor.*/
  LineConstraint(const std::vector<int> &idI,
                 Scalar weight,
                 const Matrix3X &positions);
  virtual ~LineConstraint() {}
  /** \brief Find the closest configuration from the input positions that satisfy the constraint.*/
  virtual void project(const Matrix3X &positions, Matrix3X &projections) const override final;
  /** \brief Add the constraint to the linear system.*/
  virtual void addConstraint(std::vector<Triplet> &triplets, int &idO) const override final;
 private:
  mutable Matrix3X input;
};
///////////////////////////////////////////////////////////////////////////////
/** \brief Plane constraint. See \cite Bouaziz2012 for more details.*/
class SHAPEOP_API PlaneConstraint : public Constraint {
 public:
  /** \brief Constraint constructor.*/
  PlaneConstraint(const std::vector<int> &idI,
                  Scalar weight,
                  const Matrix3X &positions);
  virtual ~PlaneConstraint() {}
  /** \brief Find the closest configuration from the input positions that satisfy the constraint.*/
  virtual void project(const Matrix3X &positions, Matrix3X &projections) const override final;
  /** \brief Add the constraint to the linear system.*/
  virtual void addConstraint(std::vector<Triplet> &triplets, int &idO) const override final;
 private:
  mutable Matrix3X input;
};
///////////////////////////////////////////////////////////////////////////////
/** \brief Circle constraint. See \cite Bouaziz2012 for more details.*/
class SHAPEOP_API CircleConstraint : public Constraint {
 public:
  /** \brief Constraint constructor.*/
  CircleConstraint(const std::vector<int> &idI,
                   Scalar weight,
                   const Matrix3X &positions);
  virtual ~CircleConstraint() {}
  /** \brief Find the closest configuration from the input positions that satisfy the constraint.*/
  virtual void project(const Matrix3X &positions, Matrix3X &projections) const override final;
  /** \brief Add the constraint to the linear system.*/
  virtual void addConstraint(std::vector<Triplet> &triplets, int &idO) const override final;
 private:
  mutable Matrix3X input;
};
///////////////////////////////////////////////////////////////////////////////
/** \brief Sphere constraint. See \cite Bouaziz2012 for more details.*/
class SHAPEOP_API SphereConstraint : public Constraint {
 public:
  /** \brief Constraint constructor.*/
  SphereConstraint(const std::vector<int> &idI,
                   Scalar weight,
                   const Matrix3X &positions);
  virtual ~SphereConstraint() {}
  /** \brief Find the closest configuration from the input positions that satisfy the constraint.*/
  virtual void project(const Matrix3X &positions, Matrix3X &projections) const override final;
  /** \brief Add the constraint to the linear system.*/
  virtual void addConstraint(std::vector<Triplet> &triplets, int &idO) const override final;
 private:
  mutable Matrix3X input;
};
///////////////////////////////////////////////////////////////////////////////
class SHAPEOP_API UniformLaplacianConstraint : public Constraint {
 public:
  /** \brief Constraint constructor.*/
  UniformLaplacianConstraint(const std::vector<int> &idI,
                             Scalar weight,
                             const Matrix3X &positions,
                             bool displacement_lap);
  virtual ~UniformLaplacianConstraint() {}
  /** \brief Find the closest configuration from the input positions that satisfy the constraint.*/
  virtual void project(const Matrix3X & /*positions*/, Matrix3X &projections) const override final;
  /** \brief Add the constraint to the linear system.*/
  virtual void addConstraint(std::vector<Triplet> &triplets, int &idO) const override final;
 private:
  Vector3 weighted_rest_;
};
///////////////////////////////////////////////////////////////////////////////
} // namespace ShapeOp
///////////////////////////////////////////////////////////////////////////////
#ifdef SHAPEOP_HEADER_ONLY
#include "Constraint.cpp"
#endif
///////////////////////////////////////////////////////////////////////////////
#endif // CONSTRAINT_H
///////////////////////////////////////////////////////////////////////////////
