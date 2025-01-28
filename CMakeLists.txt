cmake_minimum_required(VERSION 3.12)
project(PontederaROM VERSION 1.0.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY /root/bin_perso)

add_executable(PontederaROM PontederaROM.cpp)

target_include_directories(PontederaROM PUBLIC
    $ENV{FOAM_SRC}/OpenFOAM/lnInclude
    $ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude
    $ENV{FOAM_SRC}/TurbulenceModels/turbulenceModels/lnInclude
    $ENV{FOAM_SRC}/TurbulenceModels/incompressible/lnInclude
    $ENV{FOAM_SRC}/TurbulenceModels/compressible/lnInclude
    $ENV{FOAM_SRC}/TurbulenceModels/schemes/lnInclude
    $ENV{FOAM_SRC}/transportModels/incompressible/viscosityModels/viscosityModel
    $ENV{FOAM_SRC}/finiteVolume/lnInclude
    $ENV{FOAM_SRC}/meshTools/lnInclude
    $ENV{FOAM_SRC}/transportModels
    $ENV{FOAM_SRC}/transportModels/incompressible/lnInclude
    $ENV{FOAM_SRC}/dynamicMesh/lnInclude
    $ENV{FOAM_SRC}/dynamicFvMesh/lnInclude
    $ENV{FOAM_SRC}/sampling/lnInclude
    $ENV{FOAM_SRC}/fileFormats/lnInclude
    $ENV{FOAM_SRC}/thermophysicalModels/basic/lnInclude
    $ENV{FOAM_SRC}/thermophysicalModels/radiation/lnInclude
    $ENV{FOAM_SRC}/transportModels/compressible/lnInclude
    $ENV{FOAM_SRC}/fvOptions/lnInclude
    $ENV{FOAM_SRC}/functionObjects/forces/lnInclude


    $ENV{LIB_ITHACA_SRC}/ITHACA_THIRD_PARTY/splinter/include
    $ENV{LIB_ITHACA_SRC}/thirdparty/spectra/include
    $ENV{LIB_ITHACA_SRC}/thirdparty/Eigen
    $ENV{LIB_ITHACA_SRC}/thirdparty/redsvd
    $ENV{LIB_ITHACA_SRC}/ITHACA_CORE/ITHACAutilities
    $ENV{LIB_ITHACA_SRC}/ITHACA_CORE/ITHACAstream
    $ENV{LIB_ITHACA_SRC}/ITHACA_CORE/ITHACAPOD
    $ENV{LIB_ITHACA_SRC}/ITHACA_CORE/ITHACADMD
    $ENV{LIB_ITHACA_SRC}/ITHACA_CORE/ITHACAparallel
    $ENV{LIB_ITHACA_SRC}/ITHACA_CORE/EigenFunctions
    $ENV{LIB_ITHACA_SRC}/ITHACA_CORE/Foam2Eigen
    $ENV{LIB_ITHACA_SRC}/ITHACA_CORE/lnInclude
    $ENV{LIB_ITHACA_SRC}/ITHACA_FOMPROBLEMS/lnInclude
    $ENV{LIB_ITHACA_SRC}/ITHACA_FOMPROBLEMS/steadyNS
    $ENV{LIB_ITHACA_SRC}/ITHACA_FOMPROBLEMS/reductionProblem
    $ENV{LIB_ITHACA_SRC}/ITHACA_ROMPROBLEMS/lnInclude
    $ENV{LIB_ITHACA_SRC}/ITHACA_HR
)

target_link_directories(PontederaROM PUBLIC
    $ENV{FOAM_LIBBIN}
    $ENV{FOAM_LIBBIN}/dummy
    $ENV{FOAM_LIBBIN}/sys-openmpi
    $ENV{FOAM_USER_LIBBIN}
)

target_link_libraries(PontederaROM PUBLIC
    ITHACA_FOMPROBLEMS
    ITHACA_CORE
    ITHACA_THIRD_PARTY

    OpenFOAM
    dl
    Pstream
    finiteVolume
    meshTools
    dynamicFvMesh
    dynamicMesh
    turbulenceModels
    turbulenceModelSchemes
    incompressibleTransportModels
    incompressibleTurbulenceModels
    fvOptions
    sampling
    fluidThermophysicalModels
    radiationModels
    specie
    solverFunctionObjects
    fieldFunctionObjects
    initialisationFunctionObjects
    compressibleTurbulenceModels
    compressibleTransportModels
    combustionModels
    forces
    fileFormats
)

add_compile_options(
    -w
    -pthread
    -Wall
    -Wextra
    -Wold-style-cast
    -Wnon-virtual-dtor
    -Wno-unused-parameter
    -Wno-invalid-offsetof
    -Wno-attributes
    -Wno-unknown-pragmas
    -O3
    -ftemplate-depth-100
    -fPIC
)

add_compile_definitions(
    OPENFOAM=2106
    WM_LABEL_SIZE=32
    WM_DP
    NoRepository
)