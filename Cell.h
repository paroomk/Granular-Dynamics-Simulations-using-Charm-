#ifndef __CELL_H__
#define __CELL_H__

#include "ckmulticast.h"
#include <string>

//data message to be sent to computes
struct ParticleDataMsg : public CkMcastBaseMsg, public CMessage_ParticleDataMsg {
  vec2* part; //list of atoms
  int lengthAll;  //length of list
  int x, y;    // (x,y,z) coordinate of cell sending this message

  ParticleDataMsg(const int x_, const int y_, const int numPos)
    : x(x_), y(y_), lengthAll(numPos) { }

  void pup(PUP::er &p){
    CMessage_ParticleDataMsg::pup(p);
    p | lengthAll;
    p | x; p | y;
    if (p.isUnpacking()) part = new vec2[lengthAll];
    PUParray(p, part, lengthAll);
  } 
};

//chare used to represent a cell
class Cell : public CBase_Cell {
private:
  Cell_SDAG_CODE;

  // list of atoms
  std::vector<Particle> particles;
  // my compute locations
  std::vector<CkArrayIndex4D> computesList;
  // to count the number of steps, and decide when to stop
  int stepCount;
  // number of atoms in my cell
  int myNumParts;
  // number of interacting neighbors
  int inbrs;
  double stepTime;
  int updateCount;
  // store kinetic energy - initial and final
  double energy[2];
  int numReadyCheckpoint;
  void migrateToCell(Particle p, int &px, int &py);
  // updates properties after receiving forces from computes
  void updateProperties(vec2 *forces);
  // limit velcities to an upper limit
  void limitVelocity(Particle &p);
  // particles going out of right enters from left
  Particle& wrapAround(Particle &p);
  // handle to section proxy of computes
  CProxySection_Compute mCastSecProxy;

public:
  Cell();
  Cell(CkMigrateMessage *msg);
  ~Cell();
  void createComputes();  //add my computes
  void createSection();   //created multicast section of computes
  void migrateParticles();
  void writeCell(int stepCount);
  void sendPositions();
  void startCheckpoint(int);
  void pup(PUP::er &p);
};

#endif
