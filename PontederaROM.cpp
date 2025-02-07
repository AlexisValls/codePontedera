#include "unsteadyNS.H"
#include "ITHACAPOD.H"
#include "ReducedUnsteadyNS.H"
#include "ITHACAstream.H"
#include <chrono>
#include<math.h>
#include<iomanip>


int main(int argc, char* argv[])
{
    // Args
    autoPtr<argList> args = autoPtr<argList> (new argList(argc, argv));
    if (!args->checkRootCase())
    {
        Foam::FatalError.exit();
    }

    // Time
    Foam::Time time(Foam::Time::controlDictName, (fileName) ".", (fileName) "."); // two last arguments stand respectively for root path and casename -> change casename later on to parametrize

    // Mesh
    Foam::fvMesh mesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            time.timeName(),
            time,
            IOobject::MUST_READ
        )
    );

    // Parameters
    ITHACAparameters* param(ITHACAparameters::getInstance(mesh,time));

    // Velocity
    PtrList<volVectorField> U = 

    // Declaration of objects
    

    // Get wanted number of modes
//    int nModes_U = param->ITHACAdict->lookupOrDefault<int>("Nmodes_U", 16);
//    int nModes_p = param->ITHACAdict->lookupOrDefault<int>("Nmodes_p", 16);


    return 0;
}