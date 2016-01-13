#include "defs.h"
#include "Cell.h"
#include "Compute.h"
#include "physics.h"
#include "ckmulticast.h"
#include <algorithm>
using std::swap;

//compute - Default constructor
Compute::Compute() : stepCount(1) {
  energy[0] = energy[1] = 0;
  usesAtSync = true;
}

Compute::Compute(CkMigrateMessage *msg): CBase_Compute(msg)  { 
  usesAtSync = true;
  delete msg;
}

//interaction within a cell
void Compute::selfInteract(ParticleDataMsg *msg){
  double energyP = 0;
  std::vector<vec2> force1;
//  CkPrintf("entered selfInt (%d,%d,%d,%d)\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z);
  energyP = calcInternalForces(msg, stepCount, force1,(int)thisIndex.w,(int)thisIndex.x);

  //energy assignment only in begining and end
  //if (stepCount == 1) energy[0] = energyP;
  //else if (stepCount == finalStepCount) energy[1] = energyP;

  //contribute to force reduction
  CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  CkGetSectionInfo(mcast1, msg);
  mCastGrp->contribute(sizeof(vec2)*msg->lengthAll, &force1[0], CkReduction::sum_double, mcast1);
 
//  CkPrintf("exiting selfInt (%d,%d,%d,%d)\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z);
  delete msg;
}

//interaction between two cells
void Compute::interact(ParticleDataMsg *msg1, ParticleDataMsg *msg2){
  CkSectionInfo *handleA = &mcast1, *handleB = &mcast2;
//  CkPrintf("entered Interact (%d,%d,%d,%d)\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z);
  if (msg2->x * cellArrayDimY  + msg2->y  <
      msg1->x * cellArrayDimY + msg1->y )
     swap(handleA, handleB);

  std::vector<vec2> force1, force2;
  double energyP = calcPairForces(msg1, msg2, stepCount, force1, force2);

//  CkPrintf("Interact (%d,%d,%d,%d) calculated forces",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z);
  //energy assignment only in begining and end
  //if (stepCount == 1) energy[0] = energyP;
  //else if (stepCount == finalStepCount) energy[1] = energyP;

  //contribute to force reduction
  CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  CkGetSectionInfo(*handleA, msg1);
  mCastGrp->contribute(sizeof(vec2)*msg1->lengthAll, &force1[0], CkReduction::sum_double, *handleA);
  CkGetSectionInfo(*handleB, msg2);
  mCastGrp->contribute(sizeof(vec2)*msg2->lengthAll, &force2[0], CkReduction::sum_double, *handleB);

//  CkPrintf("exiting Interact (%d,%d,%d,%d)\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z);
  delete msg1;
  delete msg2;
}

//pack important information if I am moving
void Compute::pup(PUP::er &p) {
  CBase_Compute::pup(p);
  __sdag_pup(p);
  p | stepCount;
  p | mcast1;
  p | mcast2;
  //PUParray(p, energy, 2);
//  if (p.isUnpacking() && CkInRestarting()) {
//    mcast1.set_redNo(0); 
//	if(mcast2.info.type != 0)
//		mcast2.set_redNo(0);
//  }
}
