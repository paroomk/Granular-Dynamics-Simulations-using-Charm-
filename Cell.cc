#include <sstream>
#include <fstream>
#include <iostream>

#include "defs.h"
#include "leanmd.decl.h"
#include "Cell.h"
#include "ckmulticast.h"

extern char filename[20];
extern int particles_per_cell_x,particles_per_cell_y;
extern double cell_size_x,cell_size_y;

Cell::Cell() : inbrs(NUM_NEIGHBORS), stepCount(1), updateCount(0), computesList(NUM_NEIGHBORS) {
  //load balancing to be called when AtSync is called
  usesAtSync = true;

  int myid = (thisIndex.y+thisIndex.x*cellArrayDimY); 
  //myNumParts = PARTICLES_PER_CELL_START + (myid*(PARTICLES_PER_CELL_END-PARTICLES_PER_CELL_START))/(cellArrayDimX*cellArrayDimY);
  myNumParts = particles_per_cell_x*particles_per_cell_y;
  // starting random generator
  srand48(myid);

  // Particle initialization
  for(int i = 0; i < myNumParts; i++) {
    particles.push_back(Particle());
    particles[i].mass = SPHERE_MASS;

    //uniformly place particles, avoid close distance among them
    particles[i].pos.x = (GAP/(float)2) + thisIndex.x * cell_size_x + ((i)/(particles_per_cell_y))*GAP;
    particles[i].pos.y = (GAP/(float)2) + thisIndex.y * cell_size_y + ((i)%particles_per_cell_y)*GAP;
    //give random values for velocity
    particles[i].vel.x = (drand48() - 0.5) * 1 * MAX_VELOCITY;
    particles[i].vel.y = (drand48() - 0.5) * 1 * MAX_VELOCITY;
  }

  energy[0] = energy[1] = 0;
  setMigratable(false);
}

//constructor for chare object migration
Cell::Cell(CkMigrateMessage *msg): CBase_Cell(msg) {
  usesAtSync = true;
  setMigratable(false);
  delete msg;
}  

Cell::~Cell() {}

//function to create my computes
void Cell::createComputes() {
  int x = thisIndex.x, y = thisIndex.y;
  int px1, py1, dx, dy, px2, py2;

  /*  The computes X are inserted by a given cell:
   *
   *	^  X  X  X
   *	|  0  X  X
   *	y  0  0  0
   *	   x ---->
   */

  // for round robin insertion
  int currPe = CkMyPe();
  int counter=0;

  if (x==0||x==cellArrayDimX-1){
     if (y==0||y==cellArrayDimY-1)
     {
     inbrs=4;
     computesList.resize(inbrs);
     }
     else 
     {
     inbrs=6;
     computesList.resize(inbrs);
     }
  }
  else {
     if (y==0||y==cellArrayDimY-1)
     {
     inbrs=6;
     computesList.resize(inbrs);
     }
     else 
     {
     inbrs=9;
     computesList.resize(inbrs);
     }
  }

  for (int num = 0; num < NUM_NEIGHBORS; num++) {
    dx = num / (NBRS_Y) - NBRS_X/2;
    dy = num % (NBRS_Y) - NBRS_Y/2;

     px1 = x ;//+ KAWAY_X;
     py1 = y ;//+ KAWAY_Y;
     px2 = px1+dx;
     py2 = py1+dy;

    if(px2<cellArrayDimX&&py2<cellArrayDimY&&px2>=0&&py2>=0)
     {
    if (num >= NUM_NEIGHBORS / 2){
      CkArrayIndex4D index(px1, py1, px2, py2);
      computeArray[index].insert((++currPe) % CkNumPes());
      computesList[counter] = index;
      printf("compute added %d,%d,%d,%d,%d,%d \n",x,y,px1,py1,px2,py2);
    } else {
      // these computes will be created by pairing cells
    /*  px1 = WRAP_X(x + dx) + KAWAY_X;
      py1 = WRAP_Y(y + dy) + KAWAY_Y;
      px2 = px1 - dx;
      py2 = py1 - dy;*/
      CkArrayIndex4D index(px2, py2, px1, py1);
      computesList[counter] = index;
//	printf("paired compute added\n");
    }
    counter++;
   }
  } // end of for loop
  contribute(CkCallback(CkReductionTarget(Main,run),mainProxy));
}

//call multicast section creation
void Cell::createSection() {
  //knit the computes into a section
  mCastSecProxy = CProxySection_Compute::ckNew(computeArray.ckGetArrayID(), &computesList[0], computesList.size());
// printf("computesList size=%d for (%d,%d) with ID=%d\n",computesList.size(),thisIndex.x,thisIndex.y,computeArray.ckGetArrayID());
  //delegate the communication responsibility for this section to multicast library
  CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  mCastSecProxy.ckSectionDelegate(mCastGrp);
  mCastGrp->setReductionClient(mCastSecProxy, new CkCallback(CkReductionTarget(Cell,reduceForces), thisProxy(thisIndex.x, thisIndex.y)));
}

// Function to start interaction among particles in neighboring cells as well as its own particles
void Cell::sendPositions() {
  unsigned int len = particles.size();
  //create the particle and control message to be sent to computes
  ParticleDataMsg* msg = new (len) ParticleDataMsg(thisIndex.x, thisIndex.y, len);

  for(int i = 0; i < len; ++i)
    msg->part[i] = particles[i].pos;

  mCastSecProxy.calculateForces(msg);
//  CkPrintf("Calculate forces msg posted for (%d,%d)\n",thisIndex.x,thisIndex.y);
}

//send the atoms that have moved beyond my cell to neighbors
void Cell::migrateParticles(){
  int x1, y1,x2,y2;
  std::vector<std::vector<Particle> > outgoing;
  outgoing.resize(NUM_NEIGHBORS);

  int size = particles.size();
  for(std::vector<Particle>::reverse_iterator iter = particles.rbegin(); iter != particles.rend(); iter++) {
    migrateToCell(*iter, x1, y1);
    if(x1!=0 || y1!=0 ) {
      outgoing[(x1+KAWAY_X)*NBRS_Y + (y1+KAWAY_Y)].push_back((*iter));
      std::swap(*iter, particles[size - 1]);
      size--;
    }
  }
  particles.resize(size);

  for(int num = 0; num < NUM_NEIGHBORS; num++) {
    x1 = num / (NBRS_Y )            - NBRS_X/2;
    y1 = (num % (NBRS_Y ))  - NBRS_Y/2;
    x2=thisIndex.x+x1;
    y2=thisIndex.y+y1;	
    if(x2<cellArrayDimX&&y2<cellArrayDimY&&x2>=0&&y2>=0)
    cellArray(x2,y2).receiveParticles(outgoing[num]);
  }
}

//check if the particle is to be moved
void Cell::migrateToCell(Particle p, int &px, int &py) {
  double x = thisIndex.x * cell_size_x + CELL_ORIGIN_X;
  double y = thisIndex.y * cell_size_y + CELL_ORIGIN_Y;
  px = py = 0;

  //if (p.pos.x < (x-cell_size_x)) px = -2;
  if (p.pos.x < x) px = -1;
  //else if (p.pos.x > (x+2*cell_size_x)) px = 2;
  else if (p.pos.x > (x+cell_size_x)) px = 1;

  //if (p.pos.y < (y-cell_size_y)) py = -2;
  if (p.pos.y < y) py = -1;
  //else if (p.pos.y > (y+2*cell_size_y)) py = 2;
  else if (p.pos.y > (y+cell_size_y)) py = 1;

}
/*void Cell::writeCell(int stepCount)
{
    int id = thisIndex.x + thisIndex.y*cellArrayDimX ;

    std::stringstream ssParticles;
    ssParticles << "x,";
    ssParticles << "y,";
    ssParticles << "xVelocity,";
    ssParticles << "yVelocity,";
    ssParticles << std::endl;

    std::ofstream fileNameParticles;
    std::stringstream ssFileNameParticles;
    ssFileNameParticles << filename << "/step." << stepCount << ".chare." << id;

    fileNameParticles.open(ssFileNameParticles.str().c_str());
    for(int i1 = 0;i1 < particles.size();i1+=128)
    {
	ssParticles.str(std::string());
	ssParticles.clear();
	for(int i=i1;i<i1+128||i<particles.size();i++)
	{
      	    Particle p = particles[i];
      	    ssParticles << p.pos.x << ',';
            ssParticles << p.pos.y << ',';
            ssParticles << p.vel.x << ',';
      	    ssParticles << p.vel.y << ',';
            ssParticles << std::endl;
	}
    	fileNameParticles<<ssParticles.str();
    }
    fileNameParticles.close();
}*/
void Cell::writeCell(int stepCount)
{
    int id = thisIndex.x + thisIndex.y*cellArrayDimX ;

    std::stringstream ssParticles;
    ssParticles << "x,";
    ssParticles << "y,";
    ssParticles << "xVelocity,";
    ssParticles << "yVelocity,";
    ssParticles << std::endl;

    for(int i = 0;i < particles.size();i++)
    {
      Particle p = particles[i];
      ssParticles << p.pos.x << ',';
      ssParticles << p.pos.y << ',';
      ssParticles << p.vel.x << ',';
      ssParticles << p.vel.y << ',';
      ssParticles << std::endl;
    }

    std::ofstream fileNameParticles;
    std::stringstream ssFileNameParticles;
    ssFileNameParticles << filename << "/step." << stepCount << ".chare." << id;

    fileNameParticles.open(ssFileNameParticles.str().c_str());
    fileNameParticles<<ssParticles.str();
    fileNameParticles.close();
}



// Function to update properties (i.e. acceleration, velocity and position) in particles
void Cell::updateProperties(vec2 *forces) {
  int i;
  double invMassParticle;
  bool contact;
  double w_separation;

  for(i = 0; i < particles.size(); i++) {
    //calculate energy only in begining and end
   if(stepCount == 1) {
      energy[0] += (0.5 * particles[i].mass * dot(particles[i].vel, particles[i].vel) ); // in milliJoules
    } else if(stepCount == finalStepCount) { 
      energy[1] += (0.5 * particles[i].mass * dot(particles[i].vel, particles[i].vel) );
    }
    // applying kinetic equations
    invMassParticle = 1 / particles[i].mass;
    particles[i].acc = forces[i] * invMassParticle; // in m/sec^2
    particles[i].vel += particles[i].acc * DEFAULT_DELTA; // in A/fs

	if(thisIndex.x==0)
	{
		w_separation=particles[i].pos.x-RADIUS;
		contact=w_separation<0;
		if(contact)
			   particles[i].vel.x = -particles[i].vel.x;
	}
	else if(thisIndex.x==cellArrayDimX-1)
	{
		w_separation=particles[i].pos.x+RADIUS-cellArrayDimX*cell_size_x;
		contact=w_separation>0;
		if(contact)
			   particles[i].vel.x = -particles[i].vel.x;
	}
	if(thisIndex.y==0)
	{
		w_separation=particles[i].pos.y-RADIUS;
		contact=w_separation<0;
		if(contact)
			   particles[i].vel.y = -particles[i].vel.y;
	}
	else if(thisIndex.y==cellArrayDimY-1)
	{
		w_separation=particles[i].pos.y+RADIUS-cellArrayDimY*cell_size_y;
		contact=w_separation>0;
		if(contact)
			   particles[i].vel.y = -particles[i].vel.y;
	}
//    limitVelocity(particles[i]);

    particles[i].pos += particles[i].vel * DEFAULT_DELTA; // in A
  }
}
/*
inline double velocityCheck(double inVelocity) {
  if(fabs(inVelocity) > MAX_VELOCITY) {
    if(inVelocity < 0.0 )
      return -MAX_VELOCITY;
    else
      return MAX_VELOCITY;
  } else {
    return inVelocity;
  }
}

void Cell::limitVelocity(Particle &p) {
  p.vel.x = velocityCheck(p.vel.x);
  p.vel.y = velocityCheck(p.vel.y);
}

Particle& Cell::wrapAround(Particle &p) {
  if(p.pos.x < CELL_ORIGIN_X) p.pos.x += cell_size_x*cellArrayDimX;
  if(p.pos.y < CELL_ORIGIN_Y) p.pos.y += cell_size_y*cellArrayDimY;
  if(p.pos.x > CELL_ORIGIN_X + cell_size_x*cellArrayDimX) p.pos.x -= cell_size_x*cellArrayDimX;
  if(p.pos.y > CELL_ORIGIN_Y + cell_size_y*cellArrayDimY) p.pos.y -= cell_size_y*cellArrayDimY;

  return p;
}
*/
//pack important data when I move/checkpoint
void Cell::pup(PUP::er &p) {
  CBase_Cell::pup(p);
  __sdag_pup(p);
  p | particles;
  p | stepCount;
  p | myNumParts;
  p | updateCount;
  p | stepTime;
  p | inbrs;
  p | numReadyCheckpoint;
  PUParray(p, energy, 2);

  p | computesList;

  p | mCastSecProxy;
  //adjust the multicast tree to give best performance after moving
  if (p.isUnpacking()){
    if(CkInRestarting()){
      createSection();
    }
    else{
      CkMulticastMgr *mg = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
      mg->resetSection(mCastSecProxy);
      mg->setReductionClient(mCastSecProxy, new CkCallback(CkReductionTarget(Cell,reduceForces), thisProxy(thisIndex.x, thisIndex.y)));
    }
  }
}

