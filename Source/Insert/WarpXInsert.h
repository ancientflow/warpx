#include "AMReX_Vector.H"
extern void ParticleInjectionEntrance();
extern void PhiAdjustmentEntrance();
extern void BoundaryPhiSetEntrance();
extern void MyDiag();
extern void MyInit();
extern void MyOutput();
extern void CollisionRecord(amrex::Vector<int>);
extern void DataExamine();
extern void PhiGuardSetEntrance ();
extern void BeforeCollision (int step);
extern void AfterCollision (int step);