#include "CoralBenchmark.hh"
#include "CycleTracking.hh"
#include "EnergySpectrum.hh"
#include "MC_Fast_Timer.hh"
#include "MC_Particle_Buffer.hh"
#include "MC_Processor_Info.hh"
#include "MC_SourceNow.hh"
#include "MC_Time_Info.hh"
#include "MonteCarlo.hh"
#include "NVTX_Range.hh"
#include "Parameters.hh"
#include "ParticleVault.hh"
#include "ParticleVaultContainer.hh"
#include "PopulationControl.hh"
#include "SendQueue.hh"
#include "Tallies.hh"
#include "cudaFunctions.hh"
#include "cudaUtils.hh"
#include "initMC.hh"
#include "macros.hh"
#include "qs_assert.hh"
#include "utils.hh"
#include "utilsMpi.hh"
#include <cassert>
#include <cstdio>
#include <iostream>

#include "collectMetrics.cc"

#include <pwr.h>

int run(int argc, char **argv);
void gameOver();
void cycleInit(bool loadBalance);
void cycleTracking(MonteCarlo *monteCarlo);
void cycleFinalize();

using namespace std;

MonteCarlo *mcco = NULL;

int main(int argc, char **argv)
{
  printf("HERE\n");

  PWR_Grp grp;
  PWR_Obj self;
  PWR_Cntxt cntxt;
  time_t time;
  int rc;
  double value;
  PWR_Time ts, tstart, tstop, tstart2, tstop2;
  PWR_Status status;

  Metrics pwrMetrics("main_test");

  printf("HERE2\n");

  // Get a context
  rc = PWR_CntxtInit(PWR_CNTXT_DEFAULT, PWR_ROLE_APP, "App", &cntxt);
  assert(PWR_RET_SUCCESS == rc);

  printf("HERE3\n");

  rc = PWR_CntxtGetEntryPoint(cntxt, &self);
  assert(PWR_RET_SUCCESS == rc);

  printf("HERE4\n");

  PWR_ObjType objType;
  PWR_ObjGetType(self, &objType);
  printf("I am a `%s`\n", PWR_ObjGetTypeString(objType));

  PWR_Obj parent;
  rc = PWR_ObjGetParent(self, &parent);
  assert(rc >= PWR_RET_SUCCESS);

  PWR_Grp children;
  rc = PWR_ObjGetChildren(self, &children);
  assert(rc >= PWR_RET_SUCCESS);

  int i;
  for (i = 0; i < PWR_GrpGetNumObjs(children); i++)
  {
    char name[100];
    PWR_Obj obj;
    PWR_GrpGetObjByIndx(children, i, &obj);
    PWR_ObjGetName(obj, name, 100);

    printf("child %s\n", name);
  }

  pwrMetrics.getMetrics(self);
  rc = PWR_ObjAttrGetValue(self, PWR_ATTR_FREQ, &value, &ts);
  assert(PWR_RET_SUCCESS == rc);
  printf("Frequency at time %f: %lld\n", value, ts);

  rc = PWR_ObjAttrGetValue(self, PWR_ATTR_POWER, &value, &ts);
  assert(PWR_RET_SUCCESS == rc);
  printf("Power at time %f: %lld\n", value, ts);

  // Manually set power to 15 W (Not gonna work on MacOS, at least not right
  // now). value = 15.0; printf("PWR_ObjAttrSetValue(PWR_ATTR_ENERGY)
  // value=%f\n", value); rc = PWR_ObjAttrSetValue(self, PWR_ATTR_ENERGY,
  // &value); assert(PWR_RET_SUCCESS == rc);

  PWR_AttrName name = PWR_ATTR_FREQ;

  rc = PWR_StatusCreate(cntxt, &status);
  assert(PWR_RET_SUCCESS == rc);

  rc = PWR_ObjAttrGetValues(self, 1, &name, &value, &ts, status);
  assert(PWR_RET_SUCCESS == rc);

  PWR_TimeConvert(ts, &time);

  // value = 100.10;
  // printf("PWR_ObjAttrSetValues(PWR_ATTR_ENERGY) value=%f\n", value);
  // rc = PWR_ObjAttrSetValues(self, 1, &name, &value, status);
  // assert(PWR_RET_SUCCESS == rc);
  //
  // rc = PWR_ObjAttrGetValue(self, PWR_ATTR_ENERGY, &value, &ts);
  // assert(PWR_RET_SUCCESS == rc);
  //
  // PWR_TimeConvert(ts, &time);
  // printf("PWR_ObjAttrGetValue(PWR_ATTR_ENERGY) value=%f ts=`%ld`\n", value,
  //        time);
  //
  rc = PWR_CntxtGetGrpByType(cntxt, PWR_OBJ_CORE, &grp);
  assert(PWR_RET_SUCCESS == rc);

  // value = 0.1;
  // printf("PWR_GrpAttrSetValue(PWR_ATTR_ENERGY) value=%f\n", value);
  // rc = PWR_GrpAttrSetValue(grp, PWR_ATTR_ENERGY, &value, status);
  // assert(PWR_RET_SUCCESS == rc);

  rc = PWR_ObjAttrGetValue(self, PWR_ATTR_POWER, &value, &ts);
  assert(PWR_RET_SUCCESS == rc);

  printf("PWR_ObjAttrGetValue(PWR_ATTR_ENERGY) value=%f ts=`%ld`\n", value,
         time);

  PWR_Obj node;
  rc = PWR_GrpGetObjByIndx(grp, 0, &node);
  assert(PWR_RET_SUCCESS == rc);

  /*
  PWR_Stat nodeStat;
  rc = PWR_ObjCreateStat(node, PWR_ATTR_POWER, PWR_ATTR_STAT_AVG, &nodeStat);
  assert(PWR_RET_SUCCESS == rc);

  // Start measuring node stats.
  rc = PWR_StatStart(nodeStat);

  PWR_TimePeriod statTimes;
  statTimes.start = statTimes.stop = PWR_TIME_UNINIT;
  */
  double startPower, stopPower, startFreq, stopFreq;
  rc = PWR_ObjAttrGetValue(self, PWR_ATTR_POWER, &startPower, &tstart);
  assert(PWR_RET_SUCCESS == rc);
  rc = PWR_ObjAttrGetValue(self, PWR_ATTR_FREQ, &startFreq, &tstart2);
  assert(PWR_RET_SUCCESS == rc);

  run(argc, argv);

  rc = PWR_ObjAttrGetValue(self, PWR_ATTR_POWER, &stopPower, &tstop);
  assert(PWR_RET_SUCCESS == rc);
  rc = PWR_ObjAttrGetValue(self, PWR_ATTR_FREQ, &stopFreq, &tstop2);
  assert(PWR_RET_SUCCESS == rc);

  pwrMetrics.getMetrics(self);

  printf("Power start value: %lf, time: %lld\n", startPower, tstart);
  printf("Freq start value: %lf, time: %lld\n", startFreq, tstart2);
  printf("Power stop value: %lf, time: %lld\n", stopPower, tstop);
  printf("Freq stop value: %lf, time: %lld\n", stopFreq, tstop2);

  // printf("PWR_StatGetValue(PWR_ATTR_POWER) start=%lf\n",
  //        (double)statTimes.start / 1000000000);
  //
  // printf("PWR_StatGetValue(PWR_ATTR_POWER) stop=%lf\n",
  //        (double)statTimes.stop / 1000000000);
  //
  // if (statTimes.instant != PWR_TIME_UNINIT) {
  //   printf("PWR_StatGetValue(PWR_ATTR_POWER) instant=%lf\n",
  //          (double)statTimes.instant / 1000000000);
  // }
  // PWR_StatDestroy(nodeStat);
}

int run(int argc, char **argv)
{
  mpiInit(&argc, &argv);

  Parameters params = getParameters(argc, argv);
  printParameters(params, cout);

  // mcco stores just about everything.
  mcco = initMC(params);

  int loadBalance = params.simulationParams.loadBalance;

  MC_FASTTIMER_START(MC_Fast_Timer::main); // this can be done once mcco exist.

  const int nSteps = params.simulationParams.nSteps;

  for (int ii = 0; ii < nSteps; ++ii)
  {
    cycleInit(bool(loadBalance));
    cycleTracking(mcco);
    cycleFinalize();

    mcco->fast_timer->Last_Cycle_Report(params.simulationParams.cycleTimers,
                                        mcco->processor_info->rank,
                                        mcco->processor_info->num_processors,
                                        mcco->processor_info->comm_mc_world);
  }

  MC_FASTTIMER_STOP(MC_Fast_Timer::main);

  gameOver();

  coralBenchmarkCorrectness(mcco, params);

#ifdef HAVE_UVM
  mcco->~MonteCarlo();
  gpuFree(mcco);
#else
  delete mcco;
#endif

  mpiFinalize();

  return 0;
}

void gameOver()
{
  mcco->fast_timer->Cumulative_Report(
      mcco->processor_info->rank, mcco->processor_info->num_processors,
      mcco->processor_info->comm_mc_world,
      mcco->_tallies->_balanceCumulative._numSegments);
  mcco->_tallies->_spectrum.PrintSpectrum(mcco);
}

void cycleInit(bool loadBalance)
{

  MC_FASTTIMER_START(MC_Fast_Timer::cycleInit);

  mcco->clearCrossSectionCache();

  mcco->_tallies->CycleInitialize(mcco);

  mcco->_particleVaultContainer->swapProcessingProcessedVaults();

  mcco->_particleVaultContainer->collapseProcessed();
  mcco->_particleVaultContainer->collapseProcessing();

  mcco->_tallies->_balanceTask[0]._start =
      mcco->_particleVaultContainer->sizeProcessing();

  mcco->particle_buffer->Initialize();

  MC_SourceNow(mcco);

  PopulationControl(mcco, loadBalance); // controls particle population

  RouletteLowWeightParticles(
      mcco); // Delete particles with low statistical weight

  MC_FASTTIMER_STOP(MC_Fast_Timer::cycleInit);
}

#if defined GPU_NATIVE

GLOBAL void CycleTrackingKernel(MonteCarlo *monteCarlo, int num_particles,
                                ParticleVault *processingVault,
                                ParticleVault *processedVault)
{
  int global_index = getGlobalThreadID();

  if (global_index < num_particles)
  {
    CycleTrackingGuts(monteCarlo, global_index, processingVault,
                      processedVault);
  }
}

#endif

void cycleTracking(MonteCarlo *monteCarlo)
{
  MC_FASTTIMER_START(MC_Fast_Timer::cycleTracking);

  bool done = false;

  // Determine whether or not to use GPUs if they are available (set for each
  // MPI rank)
  ExecutionPolicy execPolicy =
      getExecutionPolicy(monteCarlo->processor_info->use_gpu);

  ParticleVaultContainer &my_particle_vault =
      *(monteCarlo->_particleVaultContainer);

  // Post Inital Receives for Particle Buffer
  monteCarlo->particle_buffer->Post_Receive_Particle_Buffer(
      my_particle_vault.getVaultSize());

  // Get Test For Done Method (Blocking or non-blocking
  MC_New_Test_Done_Method::Enum new_test_done_method =
      monteCarlo->particle_buffer->new_test_done_method;

  do
  {
    int particle_count = 0; // Initialize count of num_particles processed

    while (!done)
    {
      uint64_t fill_vault = 0;

      for (uint64_t processing_vault = 0;
           processing_vault < my_particle_vault.processingSize();
           processing_vault++)
      {
        MC_FASTTIMER_START(MC_Fast_Timer::cycleTracking_Kernel);
        uint64_t processed_vault =
            my_particle_vault.getFirstEmptyProcessedVault();

        ParticleVault *processingVault =
            my_particle_vault.getTaskProcessingVault(processing_vault);
        ParticleVault *processedVault =
            my_particle_vault.getTaskProcessedVault(processed_vault);

        int numParticles = processingVault->size();

        if (numParticles != 0)
        {
          NVTX_Range trackingKernel(
              "cycleTracking_TrackingKernel"); // range ends at end of scope

          // The tracking kernel can run
          // * As a cuda kernel
          // * As an OpenMP 4.5 parallel loop on the GPU
          // * As an OpenMP 3.0 parallel loop on the CPU
          // * AS a single thread on the CPU.
          switch (execPolicy)
          {
          case gpuNative:
          {
#if defined(GPU_NATIVE)
            dim3 grid(1, 1, 1);
            dim3 block(1, 1, 1);
            int runKernel = ThreadBlockLayout(grid, block, numParticles);

            // Call Cycle Tracking Kernel
            if (runKernel)
              CycleTrackingKernel<<<grid, block>>>(
                  monteCarlo, numParticles, processingVault, processedVault);

            // Synchronize the stream so that memory is copied back before we
            // begin MPI section
            gpuPeekAtLastError();
            gpuDeviceSynchronize();
#endif
          }
          break;

          case gpuWithOpenMP:
          {
            int nthreads = 128;
            if (numParticles < 64 * 56)
              nthreads = 64;
            int nteams = (numParticles + nthreads - 1) / nthreads;
            nteams = nteams > 1 ? nteams : 1;
#ifdef HAVE_OPENMP_TARGET
#pragma omp target enter data map(to : monteCarlo[0 : 1])
#pragma omp target enter data map(to : processingVault[0 : 1])
#pragma omp target enter data map(to : processedVault[0 : 1])
#pragma omp target teams distribute parallel for num_teams(nteams) \
    thread_limit(128)
#endif
            for (int particle_index = 0; particle_index < numParticles;
                 particle_index++)
            {
              CycleTrackingGuts(monteCarlo, particle_index, processingVault,
                                processedVault);
            }
#ifdef HAVE_OPENMP_TARGET
#pragma omp target exit data map(from : monteCarlo[0 : 1])
#pragma omp target exit data map(from : processingVault[0 : 1])
#pragma omp target exit data map(from : processedVault[0 : 1])
#endif
          }
          break;

          case cpu:
#include "mc_omp_parallel_for_schedule_static.hh"
            for (int particle_index = 0; particle_index < numParticles;
                 particle_index++)
            {
              CycleTrackingGuts(monteCarlo, particle_index, processingVault,
                                processedVault);
            }
            break;
          default:
            qs_assert(false);
          } // end switch
        }

        particle_count += numParticles;

        MC_FASTTIMER_STOP(MC_Fast_Timer::cycleTracking_Kernel);

        MC_FASTTIMER_START(MC_Fast_Timer::cycleTracking_MPI);

        // Next, communicate particles that have crossed onto
        // other MPI ranks.
        NVTX_Range cleanAndComm("cycleTracking_clean_and_comm");

        SendQueue &sendQueue = *(my_particle_vault.getSendQueue());
        monteCarlo->particle_buffer->Allocate_Send_Buffer(sendQueue);

        // Move particles from send queue to the send buffers
        for (int index = 0; index < sendQueue.size(); index++)
        {
          sendQueueTuple &sendQueueT = sendQueue.getTuple(index);
          MC_Base_Particle mcb_particle;

          processingVault->getBaseParticleComm(mcb_particle,
                                               sendQueueT._particleIndex);

          int buffer =
              monteCarlo->particle_buffer->Choose_Buffer(sendQueueT._neighbor);
          monteCarlo->particle_buffer->Buffer_Particle(mcb_particle, buffer);
        }

        monteCarlo->particle_buffer->Send_Particle_Buffers(); // post MPI sends

        processingVault->clear(); // remove the invalid particles
        sendQueue.clear();

        // Move particles in "extra" vaults into the regular vaults.
        my_particle_vault.cleanExtraVaults();

        // receive any particles that have arrived from other ranks
        monteCarlo->particle_buffer->Receive_Particle_Buffers(fill_vault);

        MC_FASTTIMER_STOP(MC_Fast_Timer::cycleTracking_MPI);

      } // for loop on vaults

      MC_FASTTIMER_START(MC_Fast_Timer::cycleTracking_MPI);

      NVTX_Range collapseRange("cycleTracking_Collapse_ProcessingandProcessed");
      my_particle_vault.collapseProcessing();
      my_particle_vault.collapseProcessed();
      collapseRange.endRange();

      // Test for done - blocking on all MPI ranks
      NVTX_Range doneRange("cycleTracking_Test_Done_New");
      done = monteCarlo->particle_buffer->Test_Done_New(new_test_done_method);
      doneRange.endRange();

      MC_FASTTIMER_STOP(MC_Fast_Timer::cycleTracking_MPI);

    } // while not done: Test_Done_New()

    // Everything should be done normally.
    done = monteCarlo->particle_buffer->Test_Done_New(
        MC_New_Test_Done_Method::Blocking);

  } while (!done);

  // Make sure to cancel all pending receive requests
  monteCarlo->particle_buffer->Cancel_Receive_Buffer_Requests();
  // Make sure Buffers Memory is Free
  monteCarlo->particle_buffer->Free_Buffers();

  MC_FASTTIMER_STOP(MC_Fast_Timer::cycleTracking);
}

void cycleFinalize()
{
  MC_FASTTIMER_START(MC_Fast_Timer::cycleFinalize);

  mcco->_tallies->_balanceTask[0]._end =
      mcco->_particleVaultContainer->sizeProcessed();

  // Update the cumulative tally data.
  mcco->_tallies->CycleFinalize(mcco);

  mcco->time_info->cycle++;

  mcco->particle_buffer->Free_Memory();

  MC_FASTTIMER_STOP(MC_Fast_Timer::cycleFinalize);
}
