#include <iostream>
#include <string>
#include <chrono>
#include <math.h>
#include <vector>
using namespace std;

#include "Sai2Model.h"
#include "rbdl/rbdl.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
using namespace Eigen;


const string robot_fname = "../resources/kuka_iiwa.urdf";

Eigen::Matrix3d crossMat(const Eigen::Vector3d& r) {
	Eigen::Matrix3d ret_mat;
	ret_mat << 0, -r(2), r(1),
    r(2), 0, -r(0),
    -r(1), r(0), 0;
    return ret_mat;
}

int main (int argc, char** argv) {
	auto sai2model = new Sai2Model::Sai2Model(robot_fname, false, Affine3d::Identity(), Vector3d::Zero());
	auto rbdlmodel = sai2model->_rbdl_model;

	// set joint angles and update model
	sai2model->_q.setZero();
	// sai2model->_q(1) = 90/180.0*M_PI;
	sai2model->_q << 0/180.0*M_PI,
				90/180.0*M_PI,
				-45/180.0*M_PI,
				-90/180.0*M_PI,
				0/180.0*M_PI,
				60.2/180.0*M_PI,
				187.2/180.0*M_PI;
	sai2model->updateModel();

	uint dof = sai2model->dof();
	uint size = dof*6 + dof*5 + 3;

	// initialize matrix for comparison
	MatrixXd A(size, size); // should this be sparse?
	VectorXd g(size);
	A.setZero();
	g.setZero();

	for(uint i = 0; i < dof; i++) {
		// cout << "------ " << i << " ------" << endl;
		// fill in mass and inertia values
		const auto& body = rbdlmodel->mBodies[i+1];
		double m = body.mMass;
		Matrix3d inertia = body.mInertia;
		A(6*i, 6*i) = m;
		A(6*i+1, 6*i+1) = m;
		A(6*i+2, 6*i+2) = m;
		Matrix3d body_rot = CalcBodyWorldOrientation(*rbdlmodel, sai2model->_q, i+1, false).transpose();
		Vector3d com_world = body_rot*body.mCenterOfMass;
		// cout << "com world " << com_world.transpose() << endl;
		A.block<3,3>(6*i+3, 6*i+3) = body_rot*inertia*body_rot.transpose();

		// fill in gravity
		g.segment<3>(6*i) = m*Vector3d(0,0,-9.8);

		// fill in joint constraint Jacobians
		// Note: we assume all joints are revolute here
		// Note: except for base joint, all are in relative coordinates
		Vector3d jaxis_joint_frame = rbdlmodel->mJoints[i+1].mJointAxes[0].segment<3>(0);
		if(i == 0) {
			MatrixXd Jconstr1(5, 6);
			Jconstr1.setZero();
			Jconstr1.block<3,3>(0,0) = Matrix3d::Identity();
			Jconstr1.block<3,3>(0,3) = crossMat(com_world);
			Vector3d jaxis_world = jaxis_joint_frame;
			Vector3d dir1 = jaxis_world.cross(com_world);
			if(dir1.norm() < 1e-5) {
				dir1 << 1, 0, 0;
			}
			dir1 /= dir1.norm();
			// cout << "J axis " << jaxis_joint_frame.transpose() << endl;
			Vector3d dir2 = jaxis_world.cross(dir1);
			dir2 /= dir2.norm();
			Jconstr1.block<1,3>(3,3) = dir1;
			Jconstr1.block<1,3>(4,3) = dir2;

			A.block<5,6>(6*dof + 5*i, 6*i) = Jconstr1;
			A.block<6,5>(6*i, 6*dof + 5*i) = Jconstr1.transpose();
			// cout << "J constr" << endl;
			// cout << Jconstr1 << endl;
		} else {
			uint p_id = rbdlmodel->GetParentBodyId(i+1);
			// cout << "pid " << p_id << endl;
			Matrix3d pbody_rot = CalcBodyWorldOrientation(*rbdlmodel, sai2model->_q, p_id, false).transpose();
			Matrix3d rot_parent_to_joint = rbdlmodel->X_T[i+1].E;

			MatrixXd Jconstr1(5, 6);
			Jconstr1.setZero();
			Jconstr1.block<3,3>(0,0) = Matrix3d::Identity();
			Jconstr1.block<3,3>(0,3) = crossMat(com_world);
			Vector3d jaxis_world = pbody_rot*rot_parent_to_joint*jaxis_joint_frame;
			// cout <<"jaxis_world " << jaxis_world.transpose() << endl;
			// cout <<"body rot " << endl;
			// cout << body_rot << endl;
			// cout <<"pbody rot " << endl;
			// cout << pbody_rot << endl;
			Vector3d dir1 = jaxis_world.cross(com_world);
			if(dir1.norm() < 1e-5) {
				dir1 << 1, 0, 0;
				if(dir1.dot(jaxis_world) > 0.9999) {
					dir1 << 0, 1, 0;
				}
			}
			dir1 /= dir1.norm();
			Vector3d dir2 = jaxis_world.cross(dir1);
			dir2 /= dir2.norm();
			Jconstr1.block<1,3>(3,3) = dir1;
			Jconstr1.block<1,3>(4,3) = dir2;

			A.block<5,6>(6*dof + 5*i, 6*i) = Jconstr1;
			A.block<6,5>(6*i, 6*dof + 5*i) = Jconstr1.transpose();

			MatrixXd Jconstr2(5, 6);
			Jconstr2.setZero();
			Vector3d pcom_world = pbody_rot*rbdlmodel->mBodies[p_id].mCenterOfMass;
			Vector3d jpos_parent = pbody_rot*rbdlmodel->X_T[i+1].r;
			Vector3d jpos_parent_in_world = jpos_parent - pcom_world;
			// cout << "jpos_parent " << jpos_parent.transpose() << endl;
			Jconstr2.block<3,3>(0,0) = -Matrix3d::Identity();
			Jconstr2.block<3,3>(0,3) = -crossMat(-jpos_parent_in_world);
			Jconstr2.block<1,3>(3,3) = -dir1;
			Jconstr2.block<1,3>(4,3) = -dir2;

			A.block<5,6>(6*dof + 5*i, 6*(p_id-1)) = Jconstr2;
			A.block<6,5>(6*(p_id-1), 6*dof + 5*i) = Jconstr2.transpose();

			// cout << "J constr" << endl;
			// cout << Jconstr1 << endl;

			// cout << "J constr 2" << endl;
			// cout << Jconstr2 << endl;
		}
	}
	// contact at end effector
	A.block(size-3, (dof-1)*6, 3, 3) = Matrix3d::Identity();
	A.block((dof-1)*6, size-3, 3, 3) = Matrix3d::Identity();

	// cout << "A det: " << A.determinant() << endl;
	// cout << A.block(0, 18, size, 6) << endl;
	// cout << "G:" << endl;
	// cout << g.transpose() << endl;

	SparseMatrix<double> S;
	S = A.sparseView();
	// SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;
	// solver.analyzePattern(S);
	// solver.factorize(S);
	// VectorXd sol = solver.solve(g);
	// cout << sol.transpose() << endl;

	// compare with contact space sol
	Matrix3d Lambda_inv;
	uint bodyid = 7;
	string link_name = "link6";
	Vector3d com_link6 = rbdlmodel->mBodies[bodyid].mCenterOfMass;
	MatrixXd Jv(3,dof);
	sai2model->Jv(Jv, link_name, com_link6);
	Lambda_inv = Jv*sai2model->_M_inv*(Jv.transpose());
	VectorXd grav(dof);
	sai2model->gravityVector(grav, Vector3d(0,0,-9.8));
	Vector3d rhs = -Jv*sai2model->_M_inv*grav;
	// cout << "Contact space force " << -Lambda_inv.ldlt().solve(rhs) << endl;

	uint num_tests = 10000;
	auto time_pt1 = std::chrono::high_resolution_clock::now();
	for(uint i = 0; i < num_tests; i++) {
		SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;
		solver.analyzePattern(S);
		solver.factorize(S);
		VectorXd sol = solver.solve(g);
	}
	auto time_pt2 = std::chrono::high_resolution_clock::now();
	double tbody_coord = std::chrono::duration<double>(time_pt2 - time_pt1).count();
	cout << "tbody_coord (s): " << tbody_coord << endl;

	auto time_pt3 = std::chrono::high_resolution_clock::now();
	for(uint i = 0; i < num_tests; i++) {
		VectorXd sol = Lambda_inv.ldlt().solve(rhs);
	}
	auto time_pt4 = std::chrono::high_resolution_clock::now();
	double tcspace_coord = std::chrono::duration<double>(time_pt4 - time_pt3).count();
	cout << "tcspace_coord (s): " << tcspace_coord << endl;

	return 0;
}