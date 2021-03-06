mainmodule leanmd {
  include "pup_stl.h";

  readonly CProxy_Main mainProxy;       //central controller
  readonly CProxy_Cell cellArray;     //array that houses atoms
  readonly CProxy_Compute computeArray; //computational kernels
  readonly CkGroupID mCastGrpID;      //multicast group handle
  readonly int particles_per_cell_x;
  readonly int particles_per_cell_y; 
  readonly double cell_size_x;
  readonly double cell_size_y;

 	

  readonly int cellArrayDimX;		// X dimension of the Cell array
  readonly int cellArrayDimY;		// Y dimension of the Cell array
//  readonly int cellArrayDimZ;		// Y dimension of the Cell array
  readonly int finalStepCount;		// number of steps in the simulaion
  readonly int firstLdbStep;		// begin load balancing after this many steps
  readonly int ldbPeriod;		// load balancing period
  readonly int checkptFreq;     // do FT every period
  readonly int checkptStrategy; // choose the checkpoint strategy
  readonly char filename[20];
  readonly std::string logs;   // log file for checkpointing

  //central controller chare
  mainchare [migratable] Main { 
    entry Main(CkArgMsg* msg);
    entry [reductiontarget] void energySum(double iE, double fE);  //reduction of potential energy

    // called when computes have been created
    entry [reductiontarget] void run() {
      serial {
        computeArray.doneInserting();
        CkPrintf("Computes: %d .... created\n", (5*cellArrayDimX*cellArrayDimY-3*cellArrayDimX-3*cellArrayDimY)+2);
        CkPrintf("Starting simulation .... \n\n");
        cellArray.run();
        computeArray.run();
        startBenchmarkTime = CkWallTimer();
      }
      //receive intial and final energies and compare
     when energySum(double iE1, double fE1) when energySum(double iE2, double fE2) serial {
        if(fabs(fE1 + fE2 - iE1 - iE2) > ENERGY_VAR) {
          CkPrintf("Energy value has varied significantly from %E to %E\n", iE1 + iE2, fE1 + fE2);
          CkPrintf("\nEnergy conservation test failed for maximum allowed variation of %E units.\nSIMULATION UNSUCCESSFULL\n",ENERGY_VAR);  
        } else {
          CkPrintf("\nEnergy conservation test passed for maximum allowed variation of %E units.\nSIMULATION SUCCESSFULL \n",ENERGY_VAR);
        }
        CkPrintf("Total application time %f s\n", (CkWallTimer() - startBenchmarkTime));
        CkExit();
      }
    };
  };

  //message used to send particles to computes
  message ParticleDataMsg{
    vec2 part[];
  };

  //chares to house atoms
  array [2D] Cell {
    entry Cell();  
    entry void createComputes();    //call to insert computes that I need  
    entry void receiveParticles(const std::vector<Particle> &updates);   //receive atoms that have migrated to neighboring cells to me
    entry void ResumeFromSync();    //resume from here after load balancing
    entry [reductiontarget] void startCheckpoint(); //reduction to start checkpointing
    entry void recvCheckPointDone();  //checkpointing done, resume application
    entry [reductiontarget] void reduceForces(vec2 forces[n], int n);   //receives forces from my computes on my atoms

    //function to perform iterations for Cells
    entry void run() {
      if (thisIndex.x==0 && thisIndex.y==0 ) serial {
        stepTime = CkWallTimer();
      }

      serial 
      { 
        createSection();
     //  writeCell(0);
      }

      for(stepCount = 1; stepCount <= finalStepCount; stepCount++) {
        //seurrent atom positions to my computes 
	 serial 
        { 
       CkPrintf("StepCount in Cell run %d in index (%d,%d)\n",stepCount,thisIndex.x,thisIndex.y);
     /*    if((stepCount % 50) == 0)
          {
            writeCell(stepCount);
          }*/
          sendPositions(); 
        } 

        //update properties of atoms using new force values 
        when reduceForces(vec2 forces[n], int n) serial { updateProperties(forces); }

        if ((stepCount %  MIGRATE_STEPCOUNT) == 0) {
          //send atoms that have moved beyond my cell to neighbors

          serial { migrateParticles(); } 
    
          //receive particles from my neighbors
          for(updateCount = 0; updateCount < inbrs ; updateCount++) {
            when receiveParticles(const std::vector<Particle> &updates) serial {
              for (int i = 0; i < updates.size(); ++i)
                particles.push_back(updates[i]); //add particles that have moved from neighboring cells to my cell
            }
          }
        }

        if (thisIndex.x == 0 && thisIndex.y == 0) serial {
          CkPrintf("Step %d Benchmark Time %lf ms/step\n", 
          stepCount, ((CkWallTimer() - stepTime))*1000);
          stepTime = CkWallTimer();
        }

        //periodically call load balancer
        if (stepCount >= firstLdbStep && (stepCount - firstLdbStep) % ldbPeriod == 0) {
          serial { AtSync(); }
          when ResumeFromSync() { }
        }

        //periodically checkpointing
        if (stepCount % checkptFreq == 0) {
          serial {
            //coordinate to start checkpointing
            if (thisIndex.x == 0 && thisIndex.y == 0)
              CkPrintf("[%d] CHECKPOINT at step %d\n", CkMyPe(), stepCount);
            contribute(CkCallback(CkReductionTarget(Cell,startCheckpoint),thisProxy(0,0)));
          }
          if (thisIndex.x == 0 && thisIndex.y == 0 ) {
            when startCheckpoint() serial {
              CkCallback cb(CkIndex_Cell::recvCheckPointDone(),thisProxy);
              if (checkptStrategy == 0) CkStartCheckpoint(logs.c_str(), cb);
              else CkStartMemCheckpoint(cb);
            }
          }
          when recvCheckPointDone() { }
        }

#if CMK_MEM_CHECKPOINT

        //kill one of processes to demonstrate fault tolerance
        if (stepCount == 60 && thisIndex.x == 1 && thisIndex.y == 1 ) serial {
          if (CkHasCheckpoints()) {
            CkPrintf("CkDieNow step 60\n");
            CkDieNow();
          }
        }
#endif
      }

      //everything done, reduction on kinetic energy
      serial { contribute(2*sizeof(double), energy,CkReduction::sum_double, CkCallback(CkReductionTarget(Main,energySum),mainProxy)); }
    };
  };

  //chares that do force computations for pair of cells
  array [4D] Compute {
    entry Compute();
    entry void ResumeFromSync();
    entry void calculateForces(ParticleDataMsg *msg);

    entry void run() {
      for (stepCount = 1; stepCount <= finalStepCount; stepCount++) {
        //self interaction check
        if (thisIndex.w==thisIndex.y && thisIndex.x==thisIndex.z)
          when calculateForces(ParticleDataMsg *msg) serial {// CkPrintf("Msg enter received for selfinteract for (%d,%d)&(%d,%d)\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z); 
selfInteract(msg);//CkPrintf("Msg exit received for selfinteract for (%d,%d)\n",thisIndex.w,thisIndex.x); 
	}
        else
          //receive positions from two cells
          when calculateForces(ParticleDataMsg *msg1) when calculateForces(ParticleDataMsg *msg2)
          { 
            serial {

//	 	CkPrintf("Msg enter received for interact for (%d,%d)&(%d,%d)\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z); 

              interact(msg1, msg2);// CkPrintf("Msg exit received for interact for (%d,%d)&(%d,%d)\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z); 
            }
          }
        //periodically call load balancer
        if (stepCount >= firstLdbStep && (stepCount - firstLdbStep) % ldbPeriod== 0) {
          serial { AtSync();}
          when ResumeFromSync() { }
        }
      }
      //everything done, reduction on potential energy
      serial { contribute(2*sizeof(double),energy, CkReduction::sum_double, CkCallback(CkReductionTarget(Main,energySum),mainProxy)); }
    };
  };
};
