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
#include "API.h"
#include "Solver.h"
#include "Constraint.h"
#include "Force.h"
///////////////////////////////////////////////////////////////////////////////
struct ShapeOpSolver {
  std::shared_ptr<ShapeOp::Solver> s;
};
///////////////////////////////////////////////////////////////////////////////
extern ShapeOpSolver *shapeop_create() {
  ShapeOpSolver *solver = new ShapeOpSolver;
  solver->s = std::make_shared<ShapeOp::Solver>();
  return solver;
}
extern void shapeop_delete(ShapeOpSolver *op) {
  delete op;
}
extern int shapeop_init(ShapeOpSolver *op) {
  return static_cast<int>(op->s->initialize());
}
extern int  shapeop_initDynamic(ShapeOpSolver *op,  ShapeOpScalar masses, ShapeOpScalar damping, ShapeOpScalar timestep) {
  return static_cast<int>(op->s->initialize(true, masses, damping, timestep));
}
extern int shapeop_solve(ShapeOpSolver *op, unsigned int iteration) {
  return static_cast<int>(op->s->solve(iteration));
}
extern void shapeop_setPoints(ShapeOpSolver *op, ShapeOpScalar *points, int nb_points) {
  Eigen::Map<ShapeOp::Matrix3X> p(points, 3, nb_points);
  op->s->setPoints(p);
}
extern void shapeop_getPoints(ShapeOpSolver *op, ShapeOpScalar *points, int nb_points) {
  Eigen::Map<ShapeOp::Matrix3X> p(points, 3, nb_points);
  p = op->s->getPoints();
}
extern void shapeop_setTimeStep(ShapeOpSolver *op, ShapeOpScalar timestep) {
  op->s->setTimeStep(timestep);
}
extern void shapeop_setDamping(ShapeOpSolver *op, ShapeOpScalar damping) {
  op->s->setDamping(damping);
}
///////////////////////////////////////////////////////////////////////////////
extern int shapeop_addEdgeStrainConstraint(ShapeOpSolver *op, int id1, int id2, ShapeOpScalar weight) {
  std::vector<int> ids; ids.push_back(id1); ids.push_back(id2);
  auto c = std::make_shared<ShapeOp::EdgeStrainConstraint>(ids, weight, op->s->getPoints());
  return op->s->addConstraint(c);
}
extern void shapeop_editEdgeStrainConstraint(ShapeOpSolver *op, int constraint_id, ShapeOpScalar length) {
  auto c = std::dynamic_pointer_cast<ShapeOp::EdgeStrainConstraint>(op->s->getConstraint(constraint_id)); //TODO: this will need to be robustify
  c->setEdgeLength(length);
}
extern int shapeop_addTriangleStrainConstraint(ShapeOpSolver *op, int id1, int id2, int id3, ShapeOpScalar weight) {
  std::vector<int> ids; ids.push_back(id1); ids.push_back(id2); ids.push_back(id3);
  auto c = std::make_shared<ShapeOp::TriangleStrainConstraint>(ids, weight, op->s->getPoints());
  return op->s->addConstraint(c);
}
extern int shapeop_addTetrahedronStrainConstraint(ShapeOpSolver *op, int id1, int id2, int id3, int id4, ShapeOpScalar weight) {
  std::vector<int> ids; ids.push_back(id1); ids.push_back(id2); ids.push_back(id3); ids.push_back(id4);
  auto c = std::make_shared<ShapeOp::TetrahedronStrainConstraint>(ids, weight, op->s->getPoints());
  return op->s->addConstraint(c);
}
extern int shapeop_addAreaConstraint(ShapeOpSolver *op, int id1, int id2, int id3, ShapeOpScalar weight) {
  std::vector<int> ids; ids.push_back(id1); ids.push_back(id2); ids.push_back(id3);
  auto c = std::make_shared<ShapeOp::AreaConstraint>(ids, weight, op->s->getPoints());
  return op->s->addConstraint(c);
}
extern int shapeop_addVolumeConstraint(ShapeOpSolver *op, int id1, int id2, int id3, int id4, ShapeOpScalar weight) {
  std::vector<int> ids; ids.push_back(id1); ids.push_back(id2); ids.push_back(id3); ids.push_back(id4);
  auto c = std::make_shared<ShapeOp::VolumeConstraint>(ids, weight, op->s->getPoints());
  return op->s->addConstraint(c);
}
extern int shapeop_addBendingConstraint(ShapeOpSolver *op, int *ids, int nb_ids, ShapeOpScalar weight) {
  std::vector<int> id_vector(ids, ids + nb_ids);
  auto c = std::make_shared<ShapeOp::BendingConstraint>(id_vector, weight, op->s->getPoints());
  return op->s->addConstraint(c);
}
extern int shapeop_addClosenessConstraint(ShapeOpSolver *op, int id, ShapeOpScalar weight) {
  std::vector<int> ids; ids.push_back(id);
  auto c = std::make_shared<ShapeOp::ClosenessConstraint>(ids, weight, op->s->getPoints());
  return op->s->addConstraint(c);
}
extern void shapeop_editClosenessConstraint(ShapeOpSolver *op, int constraint_id, ShapeOpScalar *point) {
  Eigen::Map<ShapeOp::Vector3> p(point, 3, 1);
  auto c = std::dynamic_pointer_cast<ShapeOp::ClosenessConstraint>(op->s->getConstraint(constraint_id)); //TODO: this will need to be robustify
  c->setPosition(p);
}
extern int shapeop_addLineConstraint(ShapeOpSolver *op, int *ids, int nb_ids, ShapeOpScalar weight) {
  std::vector<int> id_vector(ids, ids + nb_ids);
  auto c = std::make_shared<ShapeOp::LineConstraint>(id_vector, weight, op->s->getPoints());
  return op->s->addConstraint(c);
}
extern int shapeop_addPlaneConstraint(ShapeOpSolver *op, int *ids, int nb_ids, ShapeOpScalar weight) {
  std::vector<int> id_vector(ids, ids + nb_ids);
  auto c = std::make_shared<ShapeOp::PlaneConstraint>(id_vector, weight, op->s->getPoints());
  return op->s->addConstraint(c);
}
extern int shapeop_addCircleConstraint(ShapeOpSolver *op, int *ids, int nb_ids, ShapeOpScalar weight) {
  std::vector<int> id_vector(ids, ids + nb_ids);
  auto c = std::make_shared<ShapeOp::CircleConstraint>(id_vector, weight, op->s->getPoints());
  return op->s->addConstraint(c);
}
extern int shapeop_addSphereConstraint(ShapeOpSolver *op, int *ids, int nb_ids, ShapeOpScalar weight) {
  std::vector<int> id_vector(ids, ids + nb_ids);
  auto c = std::make_shared<ShapeOp::SphereConstraint>(id_vector, weight, op->s->getPoints());
  return op->s->addConstraint(c);
}
extern int shapeop_addUniformLaplacianConstraint(ShapeOpSolver *op, int *ids, int nb_ids,
                                                 int displacement_lap, ShapeOpScalar weight) {
  std::vector<int> id_vector(ids, ids + nb_ids);
  auto c = std::make_shared<ShapeOp::UniformLaplacianConstraint>(id_vector, weight, op->s->getPoints(), displacement_lap != 0);
  return op->s->addConstraint(c);
}
extern ShapeOpScalar shapeop_getConstraintError(ShapeOpSolver *op, int constraint_id) {
  return op->s->getError(constraint_id);
}
///////////////////////////////////////////////////////////////////////////////
extern int shapeop_addGravityForce(ShapeOpSolver *op, ShapeOpScalar *force) {
  Eigen::Map<ShapeOp::Vector3> g(force, 3, 1);
  auto f = std::make_shared<ShapeOp::GravityForce>(g);
  return op->s->addForces(f);
}
extern int shapeop_addVertexForce(ShapeOpSolver *op, ShapeOpScalar *force, int id) {
  Eigen::Map<ShapeOp::Vector3> g(force, 3, 1);
  auto f = std::make_shared<ShapeOp::VertexForce>(g, id);
  return op->s->addForces(f);
}
extern void shapeop_editVertexForce(ShapeOpSolver *op, int force_id, ShapeOpScalar *force, int id) {
  Eigen::Map<ShapeOp::Vector3> g(force, 3, 1);
  auto f = std::dynamic_pointer_cast<ShapeOp::VertexForce>(op->s->getForce(force_id)); //TODO: this will need to be robustify
  f->setId(id);
  f->setForce(g);
}
///////////////////////////////////////////////////////////////////////////////
