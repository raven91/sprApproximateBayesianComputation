cmake_minimum_required(VERSION 3.8)
project(sprApproximateBayesianComputation)

set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES main.cpp Engines/AbcRejectionEngine.cpp Engines/AbcRejectionEngine.hpp Definitions.hpp ParameterSets/ExperimentalParameterSet.cpp ParameterSets/ExperimentalParameterSet.hpp Engines/SimulationEngineForSyntheticData.cpp Engines/SimulationEngineForSyntheticData.hpp ParameterSets/SyntheticParameterSet.cpp ParameterSets/SyntheticParameterSet.hpp Other/PeriodicBoundaryConditionsConfiguration.cpp Other/PeriodicBoundaryConditionsConfiguration.hpp Other/Vector2D.cpp Other/Vector2D.hpp Other/Matrix2D.cpp Other/Matrix2D.hpp DynamicalSystems/SprSystemWithUnitVelocity.cpp DynamicalSystems/SprSystemWithUnitVelocity.hpp DynamicalSystems/TwoParticleInteractionForce.cpp DynamicalSystems/TwoParticleInteractionForce.hpp Observers/SyntheticObserver.cpp Observers/SyntheticObserver.hpp Steppers/StochasticEulerStepperWithUnitVelocity.hpp Observers/AbcIterationObserver.cpp Observers/AbcIterationObserver.hpp Engines/AbcSmcEngineSynthetic.cpp Engines/AbcSmcEngineSynthetic.hpp Parallelization/Thread.cpp Parallelization/Thread.hpp Parallelization/Parallelization.cpp Parallelization/Parallelization.hpp DynamicalSystems/SprSystemWithArbitraryStationaryVelocity.cpp DynamicalSystems/SprSystemWithArbitraryStationaryVelocity.hpp Steppers/StochasticEulerStepperWithArbitraryStationaryVelocity.hpp Engines/AbcSmcEngineExperimental.cpp Engines/AbcSmcEngineExperimental.hpp Engines/SimulationEngineForExperimentalData.cpp Engines/SimulationEngineForExperimentalData.hpp Other/FreeBoundaryConditionsConfiguration.cpp Other/FreeBoundaryConditionsConfiguration.hpp Other/BoundaryConditionsConfiguration.cpp Other/BoundaryConditionsConfiguration.hpp ParameterSets/AbstractParameterSet.cpp ParameterSets/AbstractParameterSet.hpp DynamicalSystems/SprSystemWithFreeBoundaries.cpp DynamicalSystems/SprSystemWithFreeBoundaries.hpp)
add_executable(sprApproximateBayesianComputation ${SOURCE_FILES})

include_directories(/usr/local/include/ /Users/nikita/Documents/mpich-install/include/)
link_directories(/usr/local/lib/ /Users/nikita/Documents/mpich-install/lib/)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()
CHECK_CXX_COMPILER_FLAG("-O3" COMPILER_SUPPORTS_O3)
if(COMPILER_SUPPORTS_O3)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif()
CHECK_CXX_COMPILER_FLAG("-flto" COMPILER_SUPPORTS_FLTO)
if(COMPILER_SUPPORTS_FLTO)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
endif()