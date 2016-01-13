#include <string>
#include "time.h"
#include "ckmulticast.h"
#include "Main.h"
#include "Cell.h"
#include "Compute.h"

/* readonly */ CProxy_Main mainProxy;
/* readonly */ CProxy_Cell cellArray;
/* readonly */ CProxy_Compute computeArray;
/* readonly */ CkGroupID mCastGrpID;

/* readonly */ int cellArrayDimX;
/* readonly */ int cellArrayDimY;
/* readonly */// int cellArrayDimZ;
/* readonly */ int finalStepCount; 
/* readonly */ int firstLdbStep; 
/* readonly */ int ldbPeriod;
/* readonly */ int checkptFreq; 
/* readonly */ int checkptStrategy;
/* readonly */ std::string logs;
/* readonly */ char filename[20];
/* readonly */ int particles_per_cell_x, particles_per_cell_y;
/* readonly */ double cell_size_x,cell_size_y;

// Entry point of Charm++ application
Main::Main(CkArgMsg* m) {
  CkPrintf("\nGRANULAR DYNAMICS START UP ...\n");

  //set variable values to a default set
  cellArrayDimX = CELLARRAY_DIM_X;
  cellArrayDimY = CELLARRAY_DIM_Y;
  particles_per_cell_x=PARTICLES_PER_CELL_X;
  particles_per_cell_y=PARTICLES_PER_CELL_Y;

//  cellArrayDimZ = CELLARRAY_DIM_Z;
  finalStepCount = DEFAULT_FINALSTEPCOUNT;
  firstLdbStep = DEFAULT_FIRST_LDB;
  ldbPeriod = DEFAULT_LDB_PERIOD;
  checkptFreq = DEFAULT_FT_PERIOD;

  mainProxy = thisProxy;

  //branch factor for spanning tree of multicast
  int bFactor = 4;
  //creating the multicast spanning tree
  mCastGrpID = CProxy_CkMulticastMgr::ckNew(bFactor);

  int numPes = CkNumPes();
  int currPe = -1, pe;
  int cur_arg = 1;
  

  CkPrintf("\nNumber of Processors: %d...\n", numPes);
  CkPrintf("\nInput Parameters...\n");

  //read user parameters
  //number of cells in each dimension
  if (m->argc > cur_arg) {
    cellArrayDimX=atoi(m->argv[cur_arg++]);
    cellArrayDimY=atoi(m->argv[cur_arg++]);
 //   cellArrayDimZ=atoi(m->argv[cur_arg++]);
   
  }

if (m->argc > cur_arg) {
    particles_per_cell_x=atoi(m->argv[cur_arg++]);
    particles_per_cell_y=atoi(m->argv[cur_arg++]);
 //   cellArrayDimZ=atoi(m->argv[cur_arg++]);
   
  }

  //number of steps in simulation
  if (m->argc > cur_arg) {
    finalStepCount=atoi(m->argv[cur_arg++]);
    CkPrintf("Final Step Count:%d\n",finalStepCount);
  }

  if (m->argc > cur_arg) {
    sprintf(filename,(m->argv[cur_arg++]));
    CkPrintf("Filename: %s\n",filename);
  }
  //step after which load balancing starts
  if (m->argc > cur_arg) {
    firstLdbStep=atoi(m->argv[cur_arg++]);
    CkPrintf("First LB Step: %d\n",firstLdbStep);
  }

  //periodicity of load balancing
  if (m->argc > cur_arg) {
    ldbPeriod=atoi(m->argv[cur_arg++]);
    CkPrintf("LB Period:%d\n",ldbPeriod);
  }

  //periodicity of checkpointing
  if (m->argc > cur_arg) {
    checkptFreq=atoi(m->argv[cur_arg++]);
    CkPrintf("FT Period:%d\n",checkptFreq);
  }

  checkptStrategy = 1;
  //choose the checkpointing strategy use in disk checkpointing
  if (m->argc > cur_arg) {
  	checkptStrategy = 0;
    logs = m->argv[cur_arg];
  }
  cell_size_x=GAP*particles_per_cell_x;
  cell_size_y=GAP*particles_per_cell_y;
 CkPrintf("Particles Array Dimension in each cell X:%d Y:%d of size %g %g\n",particles_per_cell_x,particles_per_cell_y,cell_size_x,cell_size_y);
 CkPrintf("Cell Array Dimension X:%d Y:%d in domain of size %g %g\n",cellArrayDimX,cellArrayDimY,cellArrayDimX*cell_size_x,cellArrayDimY*cell_size_y); 
 cellArray = CProxy_Cell::ckNew();
  //initializing the 3D Patch array (with a uniform distribution) and 6D compute array
  int patchCount = 0;
  float ratio = ((float)CkNumPes() - 1)/(cellArrayDimX*cellArrayDimY);
  computeArray = CProxy_Compute::ckNew();
  for (int x=0; x<cellArrayDimX; x++)
    for (int y=0; y<cellArrayDimY; y++){
        cellArray(x, y).insert((int)(patchCount++ * ratio));
        cellArray(x, y).createComputes();
      }

  cellArray.doneInserting();
  CkPrintf("\nCells: %d X %d .... created\n", cellArrayDimX, cellArrayDimY);

  delete m;
}

//constructor for chare object migration
Main::Main(CkMigrateMessage* msg): CBase_Main(msg) { 
}

//pup routine incase the main chare moves, pack important information
void Main::pup(PUP::er &p) {
  CBase_Main::pup(p);
  __sdag_pup(p);
}

#include "leanmd.def.h"
