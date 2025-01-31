/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2017 by the ITHACA-FV authors
-------------------------------------------------------------------------------
License
    This file is part of ITHACA-FV
    ITHACA-FV is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    ITHACA-FV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.
    You should have received a copy of the GNU Lesser General Public License
    along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.
Description
    Example of an unsteady NS Reduction Problem
\*---------------------------------------------------------------------------*/

#include "unsteadyNS.H"
#include "ITHACAPOD.H"
#include "ReducedUnsteadyNS.H"
#include "ITHACAstream.H"
#include <chrono>
#include<math.h>
#include<iomanip>

/*---------------------------------------------------------------------------*\
                               Starting the MAIN
\*---------------------------------------------------------------------------*/

int main(int argc, char* argv[])
{
    ////////////////////////// Initialisation of instances //////////////////////
    // Construct the main unsteadyNS object
    unsteadyNS mainCase(argc, argv);

    // Read parameters from ITHACAdict file
    ITHACAparameters* para = ITHACAparameters::getInstance(mainCase._mesh(),
                             mainCase._runTime());

    /////////////////////////// Preprocessing ///////////////////////////////////



    /////////////////////////// POD computation /////////////////////////////////
    // Retrieve mean field for velocity
    autoPtr<volVectorField> meanField_U;
    ITHACAPOD::getMeanMemoryEfficient(mainCase._U(), "ITHACAoutput/Offline/", meanField_U, false);

    // Retrieve mean field for pressure
    autoPtr<volScalarField> meanField_p;
    ITHACAPOD::getMeanMemoryEfficient(mainCase._p(), "ITHACAoutput/Offline/", meanField_p, false);

/*  
    // Perform POD on velocity field
    PtrList<volVectorField> modes_U;
    ITHACAPOD::getModesMemoryEfficient
    (
        mainCase._U(),
        "ITHACAoutput/Offline/",
        modes_U,
        "U",
        false,
        false,
        false,
        10,
        false
        //meanField_U
    );


    // Create homogeneous basis functions for velocity
    example.computeLift(example.Ufield, example.liftfield, example.Uomfield);
    // Perform a POD decomposition for velocity and pressure
    ITHACAPOD::getModes(example.Uomfield, example.Umodes, example._U().name(),
                        example.podex, 0, 0,
                        NmodesUout);
    ITHACAPOD::getModes(example.Pfield, example.Pmodes, example._p().name(),
                        example.podex, 0, 0,
                        NmodesPout);
    ITHACAPOD::getModes(example.supfield, example.supmodes, example._U().name(),
                        example.podex,
                        example.supex, 1, NmodesSUPout);
    example.projectSUP("./Matrices", NmodesUproj, NmodesPproj, NmodesSUPproj);
    reducedUnsteadyNS reduced(example);
    // Set values of the reduced stuff
    reduced.nu = 0.005;
    reduced.tstart = 50;
    reduced.finalTime = 70;
    reduced.dt = 0.005;
    reduced.storeEvery = 0.005;
    reduced.exportEvery = 0.1;
    // Set the online velocity
    Eigen::MatrixXd vel_now(1, 1);
    vel_now(0, 0) = 1;
    reduced.solveOnline_sup(vel_now, 1);
    // Reconstruct the solution and export it
    reduced.reconstruct(true, "./ITHACAoutput/ReconstructionSUP/");
    */
    exit(0);
}
