# Geant4 supports CMake 2.6.4 onwards
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

SET(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall")

# This line triggers the enabling and configuring of C and C++ compilers, sets PROJECT_SOURCE_DIR (presumably to the source directory passed to the cmake command), and sets PROJECT_BINARY_DIR (presumably to the directory where the CMake command is called)
project(WChSandBox)

# optionally disable visualisations: default to ON (build with vis enabled)
# use the "-D WITH_GEANT4_VIS = OFF" arguement when calling CMake to set this option
option(WITH_GEANT4_VIS "build with vis drivers" ON)

# Options are handled via case structures, by passing appropriate arguments to the "find_package" command. 
# You can use these options to configure datatset names and paths, CXX flags, CLHEP libraries etc.
# Find a full list here: http://geant4.web.cern.ch/geant4/UserDocumentation/UsersGuides/InstallationGuide/html/ch03s02.html 
if(WITH_GEANT4_VIS)
  find_package(Geant4 REQUIRED vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
  find_package(ROOT REQUIRED Cling TreePlayer)	# may need to add additional libraries if required; e.g. TreePlayer library for TTreeReader.

# find_package locates Geant4; or more specifically the Geant4Config.cmake file, and passes variables
# for specifying Geant4 options.
# It will look in standard locations, or its location can be specified by when calling CMake:
# "-DGeant4_DIR=/opt/Geant4/lib/Geant4-10.2.0"
# Once found it aliases the location of the Geant4Config.cmake file (packaged with Geant) to
# "Geant4_USE_FILE", which is called in the next step.

include(${Geant4_USE_FILE})
# The Geant4Config.cmake file sets a bunch of CMake variables used in the remainder of the 
# CMakeLists.txt file. These specify Geant4 directories, files, datasets and module info. 
# For example, it will set Geant4_INCLUDE_DIRS and Geant4_CXX_FLAGS. 
# It received the options given to find_package. 
# For finer control you may not call this, and set all the necessary Geant4_SOMETHING variables yourself.

include(${ROOT_USE_FILE})

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

## from example: their project was called 'event'. .. not sure needed
#ROOT_GENERATE_DICTIONARY(G__Event Event.h LINKDEF EventLinkDef.h)
##---Create a shared library with geneated dictionary
#add_library(wchsandbox SHARED Event.cxx G__Event.cxx)

# We also add the project headers
include_directories(${PROJECT_SOURCE_DIR}/include)

# Next we add the project source files. We can add them manually, which is supposedly better somehow
# or add the directory as a whole by 'globbing'
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#file({PROJECT_SOURCE_DIR}/src/DefineMRD.icc)
#file({PROJECT_SOURCE_DIR}/src/WCLitePrimaryGeneratorAction.cc)
#file({PROJECT_SOURCE_DIR}/src/MRDDetectorConstruction.cc)
#file({PROJECT_SOURCE_DIR}/src/WCLitePrimaryGeneratorMessenger.cc)
#file({PROJECT_SOURCE_DIR}/src/MRDHit.cc)
#file({PROJECT_SOURCE_DIR}/src/WCLiteRunAction.cc)
#file({PROJECT_SOURCE_DIR}/src/MRDSD.cc)
#file({PROJECT_SOURCE_DIR}/src/WCLiteStackingAction.cc)
#file({PROJECT_SOURCE_DIR}/src/SBsimMRDDB.cc)
#file({PROJECT_SOURCE_DIR}/src/WCLiteSteppingAction.cc)
#file({PROJECT_SOURCE_DIR}/src/SBsimMRDHit.cc)
#file({PROJECT_SOURCE_DIR}/src/WCLiteSteppingVerbose.cc)
#file({PROJECT_SOURCE_DIR}/src/SBsimMRDResponse.cc)
#file({PROJECT_SOURCE_DIR}/src/WCLiteTrackInformation.cc)
#file({PROJECT_SOURCE_DIR}/src/WCLiteConstructMaterials.cc)
#file({PROJECT_SOURCE_DIR}/src/WCLiteTrackingAction.cc)
#file({PROJECT_SOURCE_DIR}/src/WCLiteDetectorConstruction.cc)
#file({PROJECT_SOURCE_DIR}/src/WCLiteTrajectory.cc)
#file({PROJECT_SOURCE_DIR}/src/WCLiteEventAction.cc)
#file({PROJECT_SOURCE_DIR}/src/WCLiteTrajectoryPoint.cc)
#file({PROJECT_SOURCE_DIR}/src/mrdPMTSD.cc)
#file({PROJECT_SOURCE_DIR}/src/WCLitePhysicsList.cc)
#file({PROJECT_SOURCE_DIR}/src/mrdPMThit.cc)
#file({PROJECT_SOURCE_DIR}/src/WCLitePhysicsMessenger.cc)
#file({PROJECT_SOURCE_DIR}/src/MRDDetectorConstruction.hh)
#file({PROJECT_SOURCE_DIR}/src/WCLitePrimaryGeneratorMessenger.hh)
#file({PROJECT_SOURCE_DIR}/src/MRDHit.hh)
#file({PROJECT_SOURCE_DIR}/src/WCLiteRunAction.hh)
#file({PROJECT_SOURCE_DIR}/src/MRDSD.hh)
#file({PROJECT_SOURCE_DIR}/src/WCLiteStackingAction.hh)
#file({PROJECT_SOURCE_DIR}/src/SBsimMRDDB.hh)
#file({PROJECT_SOURCE_DIR}/src/WCLiteSteppingAction.hh)
#file({PROJECT_SOURCE_DIR}/src/SBsimMRDHit.hh)
#file({PROJECT_SOURCE_DIR}/src/WCLiteSteppingVerbose.hh)
#file({PROJECT_SOURCE_DIR}/src/SBsimMRDResponse.hh)
#file({PROJECT_SOURCE_DIR}/src/WCLiteTrackInformation.hh)
#file({PROJECT_SOURCE_DIR}/src/SBsimMRDSD.hh)
#file({PROJECT_SOURCE_DIR}/src/WCLiteTrackingAction.hh)
#file({PROJECT_SOURCE_DIR}/src/WCLiteDetectorConstruction.hh)
#file({PROJECT_SOURCE_DIR}/src/WCLiteTrajectory.hh)
#file({PROJECT_SOURCE_DIR}/src/WCLiteEventAction.hh)
#file({PROJECT_SOURCE_DIR}/src/WCLiteTrajectoryPoint.hh)
#file({PROJECT_SOURCE_DIR}/src/WCLitePhysicsList.hh)
#file({PROJECT_SOURCE_DIR}/src/WCLiteTrakectory.hh)
#file({PROJECT_SOURCE_DIR}/src/WCLitePhysicsMessenger.hh)
#file({PROJECT_SOURCE_DIR}/src/mrdPMTSD.hh)
#file({PROJECT_SOURCE_DIR}/src/WCLitePrimaryGeneratorAction.hh)
#file({PROJECT_SOURCE_DIR}/src/mrdPMThit.hh)

# Next we define the executable name and main application file, and specify its dependencies
# include all those previously specified sources and headers
add_executable(wchsandbox WChSandBox.cc ${sources} ${headers})

# and we link the application with geant4
target_link_libraries(wchsandbox ${Geant4_LIBRARIES})

message(STATUS "ROOT_LIBRARIES=${ROOT_LIBRARIES}")
target_link_libraries(wchsandbox ${ROOT_LIBRARIES})

# lastly we copy any scripts to the build directory - this helps make sure our build is completely 
# self-contained and independant from the source directory
set(WChSandBox_Scripts
  vis.mac
  MainRun.mac
  MainRunHK.mac
)

foreach(_script ${WChSandBox_Scripts})
  configure_file(
   ${PROJECT_SOURCE_DIR}/${_script}
   ${PROJECT_BINARY_DIR}/${_script}
   COPYONLY
  )
endforeach()

message(STATUS "Geant4_CXXSTD=${Geant4_CXXSTD}")

# if necessary we can also install the application
# install(TARGETS wchsandbox DESTINATION bin)

