
#ifndef __DEFS__
#define __DEFS__

#include <math.h>
#include "pup.h"

//#define HYDROGEN_MASS           (1.67 * pow( 10.0,-24)) // in g
//#define VDW_A                   (1.1328 * pow(10.0, -133)) // in (g m^2/s^2) m^12
//#define VDW_B                   (2.23224 * pow(10.0, -76)) // in (g m^2/s^2) m^6

#define ENERGY_VAR              (1.0 * pow(10.0,-5))

//average of next two should be what you want as your atom density
//this should comply with the PERDIM parameter; for KAWAY 1 1 1, the maximum number
//of particles can be 10*10*10 = 1000 - 10 comes from PERDIM parameter, which is
//currently set to be 10, using a GAP of 3; as KWAYness increases, the maximum
//number of particles decreases - for 2 1 1, it is 500, for 2 2 1 it is 250; you
//can set them to have lower values but not higher; alternatively a host of
//paramters including PTP_CUT_OFF, PERDIM, GAP can be set to suitable values to
//#define PARTICLES_PER_CELL_START        9 
//#define PARTICLES_PER_CELL_END          9

#define PARTICLES_PER_CELL_X 3
#define PARTICLES_PER_CELL_Y 3



#define DEFAULT_DELTA           1E-3	// in femtoseconds

#define DEFAULT_FIRST_LDB       20
#define DEFAULT_LDB_PERIOD      20
#define DEFAULT_FT_PERIOD       100000

#define KAWAY_X                 1
#define KAWAY_Y                 1
//#define KAWAY_Z                 1
#define NBRS_X	                (2*KAWAY_X+1)
#define NBRS_Y                  (2*KAWAY_Y+1)
//#define NBRS_Z                  (2*KAWAY_Z+1)
#define NUM_NEIGHBORS           (NBRS_X * NBRS_Y )


#define SPHERE_MASS 1E-1
#define SPRING_CONSTANT 1E5
#define DAMPING_CONSTANT 1E-2
#define RADIUS 1E-2
#define G 9.81

#define CELLARRAY_DIM_X         3 // Number chares x coord
#define CELLARRAY_DIM_Y         3 // Number chares y coord
//#define CELLARRAY_DIM_Z         3
//#define CELL_SIZE_Z             (PTP_CUT_OFF + CELL_MARGIN)/KAWAY_Z




//variables to control initial uniform placement of atoms;
//atoms should not be too close at startup for a stable system;  
//PERDIM * GAP should be less than (PTPCUTOFF+CELL_MARGIN);
//max particles per cell should not be greater thatn PERDIM^3 for 1 AWAY;
//#define PERDIM                  64 
#define GAP                     3*RADIUS 

//#define PTP_CUT_OFF            PERDIM * GAP // cut off for atom to atom interactions
//#define CELL_MARGIN             0  // constant diff between cutoff and cell size
//#define CELL_SIZE_X             (PTP_CUT_OFF + CELL_MARGIN)/KAWAY_X
//#define CELL_SIZE_Y             (PTP_CUT_OFF + CELL_MARGIN)/KAWAY_Y
#define CELL_ORIGIN_X           0
#define CELL_ORIGIN_Y	        0
//#define CELL_ORIGIN_Z	        0

#define MIGRATE_STEPCOUNT	        150
#define DEFAULT_FINALSTEPCOUNT	        1001
#define MAX_VELOCITY		        .1  //in m/s

//#define WRAP_X(a)		(((a) + cellArrayDimX) % cellArrayDimX)
//#define WRAP_Y(a)		(((a) + cellArrayDimY) % cellArrayDimY)
//#define WRAP_Z(a)		(((a) + cellArrayDimZ) % cellArrayDimZ)

struct vec2 {
  double x, y;

  vec2(double d = 0.0) : x(d), y(d) { }
  vec2(double x_, double y_) : x(x_), y(y_) { }
  void print()
  {
    printf("x = %f, y = %f \n", x, y);
  }

  inline vec2& operator += (const vec2 &rhs) {
    x += rhs.x; y += rhs.y;
    return *this;
  }
  inline vec2 operator+ (const vec2& rhs) const {
    return vec2(x + rhs.x, y + rhs.y);
  }
  inline vec2 operator+ (const double rhs) const {
    return vec2(x + rhs, y + rhs);
  }
  inline vec2& operator -= (const vec2 &rhs) {
    return *this += (rhs * -1.0);
  }
  inline vec2 operator* (const double d) const {
    return vec2(d*x, d*y);
  }
  inline vec2 operator* (const vec2 a) const {
    return vec2(x * a.x, y * a.y);
  }
  inline vec2 operator- (const vec2& rhs) const {
    return vec2(x - rhs.x, y - rhs.y);
  }
  inline vec2 operator- (const double rhs) const {
    return vec2(x - rhs, y - rhs);
  }
  inline vec2 operator/ (const double a) const {
    return vec2(x/a, y/a);
  }
	
};
inline double dot(const vec2& a, const vec2& b) {
  return a.x*b.x + a.y*b.y;
}
inline double magnitude(const vec2& a){
  return sqrt(dot(a, a));
}
PUPbytes(vec2)

//class for keeping track of the properties for a particle
struct Particle {
  double mass;
  //   Position, acceleration, velocity
  vec2 pos,acc,vel;
};
PUPbytes(Particle);

#include "leanmd.decl.h"

extern /* readonly */ CProxy_Main mainProxy;
extern /* readonly */ CProxy_Cell cellArray;
extern /* readonly */ CProxy_Compute computeArray;
extern /* readonly */ CkGroupID mCastGrpID;

extern /* readonly */ int cellArrayDimX;
extern /* readonly */ int cellArrayDimY;
//extern /* readonly */ int cellArrayDimZ;
extern /* readonly */ int finalStepCount;
extern /* readonly */ int checkptStrategy;
extern /* readonly */ std::string logs;

#endif
