# Compare compute time for contact solvers
#### Comparison of computation time for contact solvers using Body and Contact-Space coordinates on a Humanoid 4-contact problem.

### What are Contact Solvers?
Simulation and dynamics computations of robotic manipulators in contact with fixed surfaces, objects or other robots requires real-time computation of contact forces and impulses at the contacting surfaces on the robot's body. These forces are then applied in the numerical integration of the robot's equations of motion to realize physically correct behavior of non-penetration, frictional dissipation and impact dissipation. Contact Solvers are the module in the simulator that compute these contact forces and impulses by solving a set of complementarity equations with additional constraints depending on the frictional nature of the contact [1].

The inputs to a Contact Solver are the points on the robot in contact and the degrees of freedom potentially constrained by the contact (3 for point contact, 6 for surface contact). Also, the Contact Solver requires a kinematic and dynamic model of the robot.

### What solvers are we comparing?
There are many different research-grade and commercial solvers employing different contact models and coordinate representations of contact process. My comparison is focused on solvers that are capable of simultaneous solution to contact forces at multiple points on the robot. In addition, I am most interested in this study to compare the computational time for just one step in the solver, which involves setting up and solving a system of linear equations. Most solvers will typically perform this step several times to iteratively determine constraint satisfaction at each contact as a result of applying a particular set of contact conditions. It is by far the most expensive step in the solver loop.

The size of the linear system depends on the number of contacts, but also the representation of the robot kinematically. Most commonly, free body coordinates (or maximal coordinates) are used to represent the robot which involves 6 degrees-of-freedom (DOF) per rigid link. Another alternative is the generalized (or minimal) coordinates representation which uses _n_-DOFs to represent the robot of _n_ joints (assuming 1-DOF per joint). For the specific case of resolving multi-contacts, a third representation is available -- Contact Space coordinates, which require only the description of the contacts (3-DOF per point contact, 6-DOF per surface contact) [2].

The use of free body coordinates leads to a large, but sparse linear system whereas the use of either generalized coordinates or Contact Space coordinates leads to a smaller, dense linear system. Theoretically, the advantage of sparsity is in very large systems involving a large number of robots and objects. However, simulations of robot tasks typically involves only a few number of objects and robots simultaneously in contact. Therefore, I want to compare for a typical robotic example the computational time to solve the linear system between a free body coordinate representation and a Contact Space coordinates representation. It is expected that the computational time for generalized coordinates representation will lie somewhere in between.

### Method
For this study, I chose to use a scenario of a humanoid robot with 27 links, 26 joints and 33-DOF overall. The robot is modeled in contact at both of its hands and feet simultaneously, leading to 4 surface contacts. [Sai2-Model](https://github.com/manips-sai-org/sai2-model) is used to compute the kinematic and dynamic descriptions of the robot in both coordinate representations [3]. [Eigen](https://eigen.tuxfamily.org/) is used to solve the linear systems in both coordinate representations. The robot falls under the influence of gravity (g = 9.8m/s<sup>2</sup>) and all 4 contacts are assumed to be sticking i.e. no relative motion between the robot and the surface occurs at the contacts. I also assume that robot joint velocities are zero to simplify the computations, but the results are general as the joint velocities simply add additional right hand side terms in the linear system. So the computational time remains the same.

The free body coordinates representation for this problem leads to a sparse system of 316 equations -- 27*6 linear and angular body accelerations, 26*5 joint reaction forces and moments, and 4*6 contact forces and moments must be determined. The left hand matrix is represented in Eigen using its [SparseMatrix](https://eigen.tuxfamily.org/dox/group__Sparse__chapter.html) module.

The Contact Space coordinates representation leads to a dense system of 24 equations -- 4*6 contact forces and moments must be determined. The left hand matrix is represented in Eigen using a MatrixXd type and the rhs vector using a VectorXd type ([dense Eigen algebra](https://eigen.tuxfamily.org/dox/group__DenseLinearSolvers__chapter.html)).

There are two comparisons to be made. In the case of frictionless or perfectly sticking contacts, the left hand matrix is guaranteed to be positive semi-definite (PSD). This allows the use of fast Cholesky decompositions to solve the system [4]. When sliding contacts are involved, the Couloumb friction cone constraint causes the left hand matrix to lose its symmetry and positive semi-definiteness. As a result, more general decompositions for square matrices such as LU-decomposition must be used [5].

I ran the simulations on my 2015 MacBook pro with a 2.7 GHz Dual-Core Intel Core i5 processor in a single threaded application (main2.cpp in source). 

### Results
The total solver time for 10000 iterations are shown in the table below.

Coordinates | Runtime with PSD (s) | Runtime without PSD (s)
--- | --- | ---
**Free Body** | 1.64183 | 6.53227
**Contact Space** | 0.0495608 | 0.0676028

When the matrix is assumed to be positive semi-definite (PSD), the use of Contact Space coordinates results in a 30x speedup computationally over the use of Free Body coordinates. In the non-PSD case, the use of Contact Space coordinates results in a 70x speedup computationally!

### Conclusion
In this simple humanoid benchmark, it is clear that the mimimal representation of the contact dynamics through the use of Contact Space coordinates trumps the sparsity advantage of the free body coordinates. As the problem size is further increased, for instance by considering more robots in simultaneous contact, it is likely that the computational time for both methods become more similar. The Contact Space solver method brings not only computational efficiency, but also physically accurate resolution to multi-contact interactions in articulated bodies [6]. Therefore, for robotics simulation particularly where manipulation and torque control are involved, the Contact Space method is preferable over the more general free body coordinates representation.

### References
[1] Baraff, D. (1994). Fast contact force computation for nonpenetrating rigid bodies. Proceedings of the 21st Annual Conference on Computer Graphics and Interactive Techniques SIGGRAPH 94, 28(May), 23–34. https://doi.org/10.1145/192161.192168

[2] Ruspini, D., & Khatib, O. (2000). A framework for multi-contact multi-body dynamic simulation and haptic display. Proceedings. 2000 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS 2000) (Cat. No.00CH37113), 2, 1322–1327. https://doi.org/10.1109/IROS.2000.893204

[3] Khatib, O., Brock, O., Chang, K., Conti, F., Ruspini, D., Sentis, L., & Nteractive, I. (2002). Robotics and Interactive Simulation. Communications of the ACM, 45(3).

[4] https://en.wikipedia.org/wiki/Cholesky_decomposition

[5] https://en.wikipedia.org/wiki/LU_decomposition

[6] Ganguly, S. and Khatib, O., 2018, November. Experimental studies of contact space model for multi-surface collisions in articulated rigid-body systems. In International Symposium on Experimental Robotics (pp. 425-436). Springer, Cham.

_Acknowledgement: I was part of [Stanford Robotics Lab](http://cs.stanford.edu/groups/manips) when this study was performed._

_Disclaimer: I myself, or Stanford University do not own any rights to the KUKA IIWA or the DLR TORO robots. The URDF files are purely fictitious, based on estimates of kinematic and dynamic parameters. They are developed solely for the purpose of this comparative study._
