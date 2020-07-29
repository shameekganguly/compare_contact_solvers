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


const string robot_fname = "../resources/toro_headless.urdf";

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
	sai2model->_q << 0.0, //0 floating_base_px
				0.0,	//1 floating_base_py
				0.0,	//2 floating_base_pz
				0.0/180.0*M_PI, //3 floating_base_rx
				0.0/180.0*M_PI, //4 floating_base_ry
				0.0/180.0*M_PI,	//5 floating_base_rz
				-30.0/180.0*M_PI,	//6 right thigh adduction
				-45.0/180.0*M_PI,	//7 right thigh pitch
				-15/180.0*M_PI,	//8 right knee roll
				82/180.0*M_PI,	//9 right knee pitch
				15/180.0*M_PI,	//10 right ankle adduction
				-30/180.0*M_PI,	//11 right ankle pitch
				0/180.0*M_PI,	//12 trunk roll
				-15/180.0*M_PI,	//13 right shoulder pitch
				45/180.0*M_PI,	//14 right shoulder adduction
				0/180.0*M_PI,	//15 right shoulder roll
				100/180.0*M_PI,	//16 right elbow pitch
				30/180.0*M_PI,	//17 right elbow roll
				0/180.0*M_PI,	//18 right hand adduction - axis incorrect
				-15/180.0*M_PI,	//19 left shoulder pitch
				45/180.0*M_PI,	//20 left shoulder adduction
				0/180.0*M_PI,	//21 left shoulder roll
				90/180.0*M_PI,	//22 left elbow pitch
				30/180.0*M_PI,	//23 left elbow roll
				0/180.0*M_PI,	//24 left hand adduction - axis incorrect
				0/180.0*M_PI,	//25 neck roll
				// 30/180.0*M_PI,	//26 neck pitch
				30/180.0*M_PI,	//27 left thigh adduction
				-45/180.0*M_PI,	//28 left thigh pitch
				15/180.0*M_PI,	//29 left knee roll
				82/180.0*M_PI,	//30 left knee pitch
				-15/180.0*M_PI,	//31 left ankle adduction
				-30/180.0*M_PI;	//32 left ankle pitch
	sai2model->updateModel();

	uint base_dof = 6;
	uint dof = sai2model->dof();
	uint num_bodies = dof - base_dof + 1;
	uint num_joints = dof - base_dof;
	uint size = num_bodies*6 + num_joints*5 + 6*4;
	cout << "Num bodies " <<  num_bodies << endl;

	// initialize matrix for comparison
	MatrixXd A(size, size); // should this be sparse?
	VectorXd g(size);
	A.setZero();
	g.setZero();

	for(uint i = base_dof; i < dof+1; i++) {
		// cout << "------ " << i << " ------" << endl;
		// fill in mass and inertia values
		const auto& body = rbdlmodel->mBodies[i];
		double m = body.mMass;
		// cout << m << endl;
		Matrix3d inertia = body.mInertia;
		A(6*(i-base_dof), 6*(i-base_dof)) = m;
		A(6*(i-base_dof)+1, 6*(i-base_dof)+1) = m;
		A(6*(i-base_dof)+2, 6*(i-base_dof)+2) = m;
		Matrix3d body_rot = CalcBodyWorldOrientation(*rbdlmodel, sai2model->_q, i, false).transpose();
		Vector3d com_world = body_rot*body.mCenterOfMass;
		// cout << "com world " << com_world.transpose() << endl;
		A.block<3,3>(6*(i-base_dof)+3, 6*(i-base_dof)+3) = body_rot*inertia*body_rot.transpose();

		// fill in gravity
		g.segment<3>(6*(i-base_dof)) = m*Vector3d(0,0,-9.8);

		// fill in joint constraint Jacobians
		// Note: we assume all joints are revolute here
		// Note: except for base joint, all are in relative coordinates
		Vector3d jaxis_joint_frame = rbdlmodel->mJoints[i].mJointAxes[0].segment<3>(0);
		if(i > base_dof) {
			uint p_id = rbdlmodel->GetParentBodyId(i);
			// cout << "pid " << p_id << endl;
			Matrix3d pbody_rot = CalcBodyWorldOrientation(*rbdlmodel, sai2model->_q, p_id, false).transpose();
			Matrix3d rot_parent_to_joint = rbdlmodel->X_T[i].E;

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
				if(dir1.dot(jaxis_world) > 0.9999) {
					dir1 << 0, 0, 1;
				}
			}
			dir1 /= dir1.norm();
			Vector3d dir2 = jaxis_world.cross(dir1);
			dir2 /= dir2.norm();
			Jconstr1.block<1,3>(3,3) = dir1;
			Jconstr1.block<1,3>(4,3) = dir2;

			A.block<5,6>(6*num_bodies + 5*(i-base_dof-1), 6*(i-base_dof)) = Jconstr1;
			A.block<6,5>(6*(i-base_dof), 6*num_bodies + 5*(i-base_dof-1)) = Jconstr1.transpose();

			MatrixXd Jconstr2(5, 6);
			Jconstr2.setZero();
			Vector3d pcom_world = pbody_rot*rbdlmodel->mBodies[p_id].mCenterOfMass;
			Vector3d jpos_parent = pbody_rot*rbdlmodel->X_T[i].r;
			Vector3d jpos_parent_in_world = jpos_parent - pcom_world;
			// cout << "jpos_parent " << jpos_parent.transpose() << endl;
			Jconstr2.block<3,3>(0,0) = -Matrix3d::Identity();
			Jconstr2.block<3,3>(0,3) = -crossMat(-jpos_parent_in_world);
			Jconstr2.block<1,3>(3,3) = -dir1;
			Jconstr2.block<1,3>(4,3) = -dir2;

			// cout << "Row start: " <<  6*num_bodies + 5*(i-6) << endl;
			// cout << "Col start: " <<  6*(p_id-6) << endl;
			A.block<5,6>(6*num_bodies + 5*(i-base_dof-1), 6*(p_id-base_dof)) = Jconstr2;
			A.block<6,5>(6*(p_id-base_dof), 6*num_bodies + 5*(i-base_dof-1)) = Jconstr2.transpose();

			// cout << "J constr" << endl;
			// cout << Jconstr1 << endl;

			// cout << "J constr 2" << endl;
			// cout << Jconstr2 << endl;
		}
	}
	// contact at end effectors
	string left_leg_name = "LL_foot";
	uint left_leg_body_id = sai2model->linkId(left_leg_name);
	cout << "Left leg id " << left_leg_body_id << endl;
	A.block(size-24, (left_leg_body_id-base_dof)*6, 6, 6) = Matrix3d::Identity();
	A.block((left_leg_body_id-base_dof)*6, size-24, 6, 6) = Matrix3d::Identity();

	string right_leg_name = "RL_foot";
	uint right_leg_body_id = sai2model->linkId(right_leg_name);
	cout << "Right leg id " << right_leg_body_id << endl;
	A.block(size-18, (right_leg_body_id-base_dof)*6, 6, 6) = Matrix3d::Identity();
	A.block((right_leg_body_id-base_dof)*6, size-18, 6, 6) = Matrix3d::Identity();

	string left_hand_name = "la_link6";
	uint left_hand_body_id = sai2model->linkId(left_hand_name);
	cout << "Left hand id " << left_hand_body_id << endl;
	A.block(size-12, (left_hand_body_id-base_dof)*6, 6, 6) = Matrix3d::Identity();
	A.block((left_hand_body_id-base_dof)*6, size-12, 6, 6) = Matrix3d::Identity();

	string right_hand_name = "ra_link6";
	uint right_hand_body_id = sai2model->linkId(right_hand_name);
	cout << "Right hand id " << right_hand_body_id << endl;
	A.block(size-6, (right_hand_body_id-base_dof)*6, 6, 6) = Matrix3d::Identity();
	A.block((right_hand_body_id-base_dof)*6, size-6, 6, 6) = Matrix3d::Identity();

	cout << "A det: " << A.determinant() << endl;
	// cout << A.block(num_bodies*6,0,size-num_bodies*6,6) << endl;
	// cout << "G:" << endl;
	// cout << g.transpose() << endl;

	SparseMatrix<double> S;
	S = A.sparseView();
	S.makeCompressed();
	// SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;
	// solver.analyzePattern(S);
	// solver.factorize(S);
	// cout << solver.info() << endl;
	// Eigen::SimplicialCholesky<SparseMatrix<double>> solver(S);  // performs a Cholesky factorization of A
  	// VectorXd sol = solver.solve(g);
  	// VectorXd sol = A.partialPivLu().solve(g);
	// cout << sol.tail(12).transpose() << endl;
	// cout << "Left leg contact force: " << sol.segment<3>(size - 6).transpose() << endl;
	// cout << "Right leg contact force: " << sol.segment<3>(size - 3).transpose() << endl;


	// compare with contact space sol
	MatrixXd Lambda_inv(24,24);
	MatrixXd J_c(24,dof);
	string llink_name = left_leg_name;
	Vector3d com_llink = rbdlmodel->mBodies[left_leg_body_id].mCenterOfMass;
	MatrixXd Jvll(6,dof);
	sai2model->J_0(Jvll, llink_name, com_llink);
	J_c.block(0,0,6,dof) = Jvll;
	string rlink_name = right_leg_name;
	Vector3d com_rlink = rbdlmodel->mBodies[right_leg_body_id].mCenterOfMass;
	MatrixXd Jvrl(6,dof);
	sai2model->J_0(Jvrl, rlink_name, com_rlink);
	J_c.block(6,0,6,dof) = Jvrl;
	string lhlink_name = left_hand_name;
	Vector3d com_lhlink = rbdlmodel->mBodies[left_hand_body_id].mCenterOfMass;
	MatrixXd Jvlh(6,dof);
	sai2model->J_0(Jvlh, lhlink_name, com_lhlink);
	J_c.block(12,0,6,dof) = Jvlh;
	string rhlink_name = right_hand_name;
	Vector3d com_rhlink = rbdlmodel->mBodies[right_hand_body_id].mCenterOfMass;
	MatrixXd Jvrh(6,dof);
	sai2model->J_0(Jvrh, rhlink_name, com_rhlink);
	J_c.block(24,0,6,dof) = Jvrh;

	Lambda_inv = J_c*sai2model->_M_inv*(J_c.transpose());
	VectorXd grav(dof);
	sai2model->gravityVector(grav, Vector3d(0,0,-9.8));
	VectorXd rhs(24);
	rhs = -J_c*sai2model->_M_inv*grav;
	// cout << "Contact space force " << -Lambda_inv.ldlt().solve(rhs).transpose() << endl;

	cout << "--- PSD solver comparison: ---" << endl;
	uint num_tests = 10000;
	auto time_pt1 = std::chrono::high_resolution_clock::now();
	for(uint i = 0; i < num_tests; i++) {
		// SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;
		// solver.analyzePattern(S);
		// solver.factorize(S);
		SimplicialCholesky<SparseMatrix<double>> solver(S);
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

	cout << "--- Non-PSD solver comparison: ---" << endl;
	auto time_pt5 = std::chrono::high_resolution_clock::now();
	for(uint i = 0; i < num_tests; i++) {
		SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;
		solver.analyzePattern(S);
		solver.factorize(S);
		VectorXd sol = solver.solve(g);
	}
	auto time_pt6 = std::chrono::high_resolution_clock::now();
	double tbody_coord2 = std::chrono::duration<double>(time_pt6 - time_pt5).count();
	cout << "tbody_coord (s): " << tbody_coord2 << endl;

	auto time_pt7 = std::chrono::high_resolution_clock::now();
	for(uint i = 0; i < num_tests; i++) {
		VectorXd sol = Lambda_inv.partialPivLu().solve(rhs);
	}
	auto time_pt8 = std::chrono::high_resolution_clock::now();
	double tcspace_coord2 = std::chrono::duration<double>(time_pt8 - time_pt7).count();
	cout << "tcspace_coord (s): " << tcspace_coord2 << endl;

	return 0;
}