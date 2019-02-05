#include "Engines/AbcRejectionEngine.hpp"
#include "Engines/AbcSmcEngineSynthetic.hpp"
#include "Engines/AbcSmcEngineExperimental.hpp"
#include "Engines/SimulationEngineForSyntheticData.hpp"
#include "ParameterSets/SyntheticParameterSet.hpp"
#include "Other/PeriodicBoundaryConditionsConfiguration.hpp"
#include "Parallelization/Parallelization.hpp"

int main(int argc, char **argv)
{
  LaunchParallelSession(argc, argv);
  Thread thread(argc, argv);

  SyntheticParameterSet parameter_set;
  PeriodicBoundaryConditionsConfiguration pbc_config(parameter_set.Get_L(), parameter_set.Get_L());
  SimulationEngineForSyntheticData simulation_engine(parameter_set, pbc_config);
  simulation_engine.GenerateSyntheticDataWithArbitraryStationaryVelocity();

//  AbcRejectionEngine abc_rejection_engine;
//    abc_rejection_engine.PrepareSyntheticData();
//  abc_rejection_engine.RunAbcTowardsSyntheticData();

//  AbcSmcEngineSynthetic abc_smc_engine;
//  abc_smc_engine.PrepareSyntheticData();
//  abc_smc_engine.RunAbcTowardsSyntheticData();

//  AbcSmcEngineExperimental abc_smc_engine;
//  abc_smc_engine.RunAbcTowardsExperimentalData();

  FinalizeParallelSession();

  return 0;
}