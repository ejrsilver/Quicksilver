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
#include <atomic>
#include <cassert>
#include <cstdio>
#include <iostream>
#include <sys/types.h>
#include <thread>
#include <unistd.h>

// #include "collectMetrics.cc"

#include <pwr.h>

int run(int argc, char **argv);
void gameOver();
void cycleInit(bool loadBalance);
void cycleTracking(MonteCarlo *monteCarlo);
void cycleFinalize();

void getGovAttrName(PWR_AttrGov gov, char *govstr, int len)
{
  if (gov == PWR_GOV_LINUX_ONDEMAND)
  {
    sprintf(govstr, "ondemand");
  }
  else if (gov == PWR_GOV_LINUX_PERFORMANCE)
  {
    sprintf(govstr, "performance");
  }
  else if (gov == PWR_GOV_LINUX_CONSERVATIVE)
  {
    sprintf(govstr, "conservative");
  }
  else if (gov == PWR_GOV_LINUX_POWERSAVE)
  {
    sprintf(govstr, "powersave");
  }
  else if (gov == PWR_GOV_LINUX_USERSPACE)
  {
    sprintf(govstr, "userspace");
  }
  else
  {
    sprintf(govstr, "schedutil");
  }
}

using namespace std;

MonteCarlo *mcco = NULL;

// variable set by main program to signal that program is still running
atomic<bool> running(true);

PWR_Obj node;
PWR_Cntxt cntxt;
int rc;
double value;
PWR_Time ts, tstart, tstop, tstart2, tstop2;
PWR_Status status;

PWR_Grp sockets;
char name[100];
uint64_t energy;
PWR_Obj socket;
PWR_ObjType socketType;
char socket_name[100];
PWR_Grp cores1;
PWR_AttrGov gov;
PWR_Obj core;
PWR_ObjType coreType;

// To create hint for sockets on CAC
uint64_t region_id_serial;
uint64_t region_id_parallel;

int main(int argc, char **argv)
{

  // Get a context
  rc = PWR_CntxtInit(PWR_CNTXT_DEFAULT, PWR_ROLE_APP, "App", &cntxt);
  assert(PWR_RET_SUCCESS == rc);
  rc = PWR_CntxtGetEntryPoint(cntxt, &node);
  assert(PWR_RET_SUCCESS == rc);

  PWR_ObjType objType;
  PWR_ObjGetType(node, &objType);
  assert(objType == PWR_OBJ_NODE);

  rc = PWR_ObjGetChildren(node, &sockets);
  assert(rc >= PWR_RET_SUCCESS);

  // Get first socket.
  PWR_GrpGetObjByIndx(sockets, 0, &socket);

  // Assert that we're reading from a socket, so that we know it has energy.
  // Leaving in frequency stuff even though right now it'll all be zeros.
  // Gonna try to add it in, but more doc reading necessary.
  // PWR_ObjGetType(socket1, &socketType);
  // assert(socketType == PWR_OBJ_SOCKET);

  PWR_ObjGetName(socket, socket_name, 100);

  // rc = PWR_ObjGetChildren(socket1, &cores1);
  // assert(rc >= PWR_RET_SUCCESS);

  // uint64_t max_freq, min_freq, init_freq, target_freq, current_freq = 0;

  // PWR_GrpGetObjByIndx(cores1, 0, &core);
  // PWR_ObjGetName(core, name, 100);

  // ############ START This stuff is just to ensure pwrAPI works as expected ############

  // PWR_ObjAttrGetValue(core, PWR_ATTR_FREQ, &init_freq, &ts);
  // assert(PWR_RET_SUCCESS == rc);
  // printf("Initial Frequency %lu\n", init_freq);

  // // If the core isn't already in userspace mode, set it.
  // gov = PWR_GOV_LINUX_USERSPACE;
  // PWR_ObjAttrSetValue(core, PWR_ATTR_GOV, &gov);
  // assert(PWR_RET_SUCCESS == rc);

  // target_freq = 2000000;
  // printf("Setting target frequency to %lu\n", target_freq);
  // PWR_ObjAttrSetValue(core, PWR_ATTR_FREQ, &target_freq);
  // assert(PWR_RET_SUCCESS == rc);

  // PWR_ObjAttrGetValue(core, PWR_ATTR_FREQ, &current_freq, &ts);
  // assert(PWR_RET_SUCCESS == rc);
  // printf("Current Frequency %lu\n", current_freq);

  // target_freq = 1200000;
  // printf("Setting target frequency to %lu\n", target_freq);
  // PWR_ObjAttrSetValue(core, PWR_ATTR_FREQ, &target_freq);
  // assert(PWR_RET_SUCCESS == rc);

  // PWR_ObjAttrGetValue(core, PWR_ATTR_FREQ, &current_freq, &ts);
  // assert(PWR_RET_SUCCESS == rc);
  // printf("Current Frequency %lu\n", current_freq);

  // gov = PWR_GOV_LINUX_POWERSAVE;
  // PWR_ObjAttrSetValue(core, PWR_ATTR_GOV, &gov);
  // assert(PWR_RET_SUCCESS == rc);

  // ############ ENDOF This stuff is just to ensure pwrAPI works as expected ############

  PWR_AppHintCreate(socket, socket_name, &region_id_parallel, PWR_REGION_PARALLEL);
  PWR_AppHintCreate(socket, socket_name, &region_id_serial, PWR_REGION_SERIAL);

  run(argc, argv);
  PWR_AppHintDestroy(&region_id_parallel);
  // set monitor thread to stop.
  running = false;
  return 0;
}

int run(int argc, char **argv)
{
  mpiInit(&argc, &argv);

  Parameters params = getParameters(argc, argv);
  printParameters(params, cout);

  // mcco stores just about everything.
  PWR_AppHintStart(&region_id_serial);
  mcco = initMC(params);
  PWR_AppHintStop(&region_id_serial);

  int loadBalance = params.simulationParams.loadBalance;

  MC_FASTTIMER_START(MC_Fast_Timer::main); // this can be done once mcco exist.

  const int nSteps = params.simulationParams.nSteps;

  for (int ii = 0; ii < nSteps; ++ii)
  {
    PWR_AppHintStart(&region_id_serial);
    cycleInit(bool(loadBalance));
    PWR_AppHintStop(&region_id_serial);

    PWR_AppHintStart(&region_id_parallel);
    cycleTracking(mcco);
    PWR_AppHintStop(&region_id_parallel);

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
