#ifndef __PHYSICS_H__
#define __PHYSICS_H__

#include "ckmulticast.h"
#include "defs.h"

//#define BLOCK_SIZE	512

//function to calculate forces among 2 lists of atoms
inline double calcPairForces(ParticleDataMsg* first, ParticleDataMsg* second, int stepCount, std::vector<vec2>& force1, std::vector<vec2> &force2) {
  int i, j;
  int firstLen = first->lengthAll;
  int secondLen = second->lengthAll;
  double  dist_r, rsqd;
  vec2 separation, force;
   double energy = 0;
 // int doEnergy = 0;
 // if(stepCount == 1 || stepCount == finalStepCount)
 //   doEnergy = 1;
  double delta;
  force1.resize(firstLen);
  force2.resize(secondLen);
 
 // int i1, j1;
 //	CkPrintf("firstLen= %d, secondLen= %d\n",firstLen,secondLen);
  for(i = 0; i < firstLen; i=i+1){
    for(j = 0; j < secondLen; j=j+1){
      separation = first->part[i]-second->part[j];
      rsqd = dot(separation, separation);
      bool contact = (rsqd < 4*RADIUS*RADIUS);
          
           if (contact)
           {
		dist_r=sqrt(rsqd);

           //   cos_n   = (separation.x)/sqrt(rsqd);
           //   sin_n   = (separation.y)/sqrt(rsqd);

              delta = (2*RADIUS-dist_r);
		force=(separation * (double)(SPRING_CONSTANT*delta)) / dist_r;
		force1[i]+=  force;
		force2[j]-= force;
	     //  printf("Collision Detected for particle %d with particle %d\n",(p+j)->MyName-1,(p+i)->MyName-1);
 
           }
          }
        }
//	CkPrintf("Exiting force Pair calc");
  return energy;
}

//function to calculate forces among atoms in a single list
inline double calcInternalForces(ParticleDataMsg* first, int stepCount, std::vector<vec2>& force1,int cell_index_x,int cell_index_y) {
  int i, j;
  int firstLen = first->lengthAll;
  double powTwenty, powTen, firstx, firsty, rx, ry,  r, rsqd, fx, fy, f, fr;
  vec2 firstpos, separation, force;
  double w_separation,dist_r;
  double energy = 0;
  double delta;
 // if(stepCount == 1 || stepCount == finalStepCount)
   // doEnergy = 1;
  force1.resize(firstLen);
bool contact;
//  powTen = pow(10.0, -10);
//  powTwenty = pow(10.0, -20);
  for(i = 0; i < firstLen; i++){
	contact=0;
    firstpos = first->part[i];
   for(j = i+1; j < firstLen; j++) {
    
  // computing base values

 //	separation = sqrt(((p+j)->x-(p+i)->x)*((p+j)->x-(p+i)->x) + ((p+j)->y-(p+i)->y)*((p+j)->y-(p+i)->y));
      separation = firstpos - first->part[j];
      rsqd = dot(separation, separation);
      contact = (rsqd < 4*RADIUS*RADIUS);
          
           if (contact)
           {
		dist_r=sqrt(rsqd);
                     //     cos_n   = (separation.x)/dist_r;
         //     sin_n   = (separation.y)/sqrt(rsqd);

              delta = (2*RADIUS-dist_r);
		force=(separation * (double)(SPRING_CONSTANT*delta)) / dist_r;
		force1[i]+=  force;
		force1[j]-= force;
	     //  printf("Collision Detected for particle %d with particle %d\n",(p+j)->MyName-1,(p+i)->MyName-1);
 
           }
    }
  }
  return energy;
}
#endif
