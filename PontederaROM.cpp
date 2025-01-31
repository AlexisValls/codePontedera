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

template<class Type, template<class> class PatchField, class GeoMesh>
void ITHACAPOD::getModesMemoryEfficient(
    GeometricField<Type, PatchField, GeoMesh>& templateField,
    word snapshotsPath,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& modes,
    word fieldName,
    bool podex,
    bool supex,
    bool sup,
    label nmodes,
    bool correctBC,
    autoPtr<GeometricField<Type, PatchField, GeoMesh>> meanField)
{
    // Get parameters instance for POD settings
    ITHACAparameters* para(ITHACAparameters::getInstance());
    word PODkey = "POD_" + fieldName;
    word PODnorm = para->ITHACAdict->lookupOrDefault<word>(PODkey, "L2");
    
    // Verify valid norm selection
    M_Assert(PODnorm == "L2" || PODnorm == "Frobenius", 
             "The PODnorm can be only L2 or Frobenius");
    
    Info << "Performing memory efficient POD for " << fieldName 
         << " using the " << PODnorm << " norm" << endl;

    if ((podex == 0 && sup == 0) || (supex == 0 && sup == 1))
    {
        // Count number of snapshots in directory (excluding 0/ and constant/)
        fileName rootPath(".");
        Foam::Time runTime2(Foam::Time::controlDictName, rootPath, snapshotsPath);
        label nSnaps = runTime2.times().size() - 2;

        std::cout << "Found " << nSnaps << " time directories" << endl;

        // Verify we have at least one snapshot
        if (nSnaps < 1)
        {
            FatalError
                << "Error: No time directories found in " << snapshotsPath
                << exit(FatalError);
        }

        // Set number of modes based on eigensolver type
        if (para->eigensolver == "spectra")
        {
            if (nmodes == 0)
            {
                nmodes = nSnaps - 2; 
            }
            M_Assert(nmodes <= nSnaps - 2,
                    "The number of requested modes cannot be bigger than the number of snapshots - 2");
        }
        else
        {
            if (nmodes == 0)
            {
                nmodes = nSnaps;
            }
            M_Assert(nmodes <= nSnaps,
                    "The number of requested modes cannot be bigger than the number of snapshots");
        }

        // Initialize correlation matrix and boundary data structures
        Eigen::MatrixXd _corMatrix(nSnaps, nSnaps);
        _corMatrix.setZero();
        List<Eigen::MatrixXd> SnapMatrixBC;
        label NBC = templateField.boundaryField().size();
        SnapMatrixBC.resize(NBC);

        // Initialize matrices for boundary conditions
        for (label i = 0; i < NBC; i++)
        {
            SnapMatrixBC[i].resize(templateField.boundaryField()[i].size(), nSnaps);
        }

        // Build correlation matrix by processing snapshots sequentially
        for (label i = 0; i < nSnaps; i++)
        {
            // Read snapshot i
            GeometricField<Type, PatchField, GeoMesh> snapI = 
                ITHACAstream::readFieldByIndex(templateField, snapshotsPath, i);

            // Substract mean field if provided
            if (meanField)
            {
                snapI -= *meanField;
            }

            // Store boundary field data for snapshot i
            List<Eigen::VectorXd> snapIBC = Foam2Eigen::field2EigenBC(snapI);
            for (label k = 0; k < NBC; k++)
            {
                SnapMatrixBC[k].col(i) = snapIBC[k];
            }

            // Compute correlations with all subsequent snapshots
            for (label j = i; j < nSnaps; j++)
            {
                GeometricField<Type, PatchField, GeoMesh> snapJ = 
                    ITHACAstream::readFieldByIndex(templateField, snapshotsPath, j);

                // Subtract mean field if provided
                if (meanField)
                {
                    snapJ -= *meanField;
                }

                // Calculate correlation using specified norm
                if (PODnorm == "L2")
                {
                    _corMatrix(i,j) = computeInnerProduct(snapI, snapJ);
                }
                else // Frobenius norm
                {
                    _corMatrix(i,j) = computeFrobeniusInnerProduct(snapI, snapJ);
                }
                
                // Matrix is symmetric - copy value to lower triangle
                if (i != j)
                {
                    _corMatrix(j,i) = _corMatrix(i,j);
                }
            }
            
            Info << "Processed snapshot " << i + 1 << " of " << nSnaps << endl;
        }

        // Sum up correlation matrix across processors if running in parallel
        if (Pstream::parRun())
        {
            reduce(_corMatrix, sumOp<Eigen::MatrixXd>());
        }

        // Solve eigenvalue problem using selected solver
        Eigen::VectorXd eigenValues;
        Eigen::MatrixXd eigenVectors;

        Info << "####### Performing the POD using EigenDecomposition " <<
             fieldName << " #######" << endl;
        
        if (para->eigensolver == "spectra")
        {
            // Use Spectra solver for large eigenvalue problems
            std::cout << "Using Spectra EigenSolver " << std::endl;
            Spectra::DenseSymMatProd<double> op(_corMatrix);
            Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE, 
                                 Spectra::DenseSymMatProd<double>> 
                solver(&op, nmodes, nSnaps);
            
            solver.init();
            solver.compute(1000, 1e-10, Spectra::LARGEST_ALGE);
            
            M_Assert(solver.info() == Spectra::SUCCESSFUL,
                    "Eigenvalue decomposition failed");
                    
            eigenVectors = solver.eigenvectors().real();
            eigenValues = solver.eigenvalues().real();
        }
        else if (para->eigensolver == "eigen")
        {
            // Use Eigen solver for smaller problems
            std::cout << "Using Eigen EigenSolver " << std::endl;
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(_corMatrix);
            M_Assert(solver.info() == Eigen::Success,
                    "Eigenvalue decomposition failed");
                    
            eigenVectors = solver.eigenvectors().real().rowwise().reverse().leftCols(nmodes);
            eigenValues = solver.eigenvalues().real().array().reverse();
        }

        // Handle negative eigenvalues if they occur
        if (eigenValues.array().minCoeff() < 0)
        {
            eigenValues = eigenValues.array() + 2 * abs(
                                 eigenValues.array().minCoeff());
        }

        Info << "####### End of the POD for " << fieldName << " #######" << endl;

        // Construct POD modes
        modes.resize(nmodes);
        
        // Read first snapshot to get boundary conditions
        GeometricField<Type, PatchField, GeoMesh> firstSnap = 
            ITHACAstream::readFieldByIndex(templateField, snapshotsPath, 0);
                
        for (label i = 0; i < nmodes; i++)
        {                
            // Initialize mode with proper dimensions and boundary conditions
            GeometricField<Type, PatchField, GeoMesh> modeI
            (
                IOobject
                (
                    templateField.name(),
                    templateField.time().timeName(),
                    templateField.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                templateField.mesh(),
                dimensioned<Type>("zero", templateField.dimensions(), Zero),
                firstSnap.boundaryField().types()
            );

            // Construct mode as linear combination of snapshots
            for (label j = 0; j < nSnaps; j++)
            {
                GeometricField<Type, PatchField, GeoMesh> snapJ = 
                    ITHACAstream::readFieldByIndex(templateField, snapshotsPath, j);
                modeI += snapJ * eigenVectors(j,i);
            }

            // Calculate normalization factor based on selected norm
            scalar normFactor;
            if (PODnorm == "L2")
            {
                normFactor = computeInnerProduct(modeI, modeI);
            }
            else // Frobenius norm
            {
                normFactor = computeFrobeniusInnerProduct(modeI, modeI);
            }

            if (Pstream::parRun())
            {
                reduce(normFactor, sumOp<scalar>());
            }
            normFactor = Foam::sqrt(normFactor);

            // Normalize the mode
            modeI *= dimensionedScalar("normFactor", dimless, 1.0/normFactor);

            // Apply boundary conditions
            for (label k = 0; k < NBC; k++)
            {
                Eigen::VectorXd bcValues = SnapMatrixBC[k] * eigenVectors.col(i);
                bcValues = bcValues / normFactor;
                ITHACAutilities::assignBC(modeI, k, bcValues);
            }

            if (correctBC)
            {
                modeI.correctBoundaryConditions();
            }

            modes.set(i, modeI.clone());
            Info << "Constructed mode " << i + 1 << " of " << nmodes << endl;
        }

        // Save modes to appropriate directory
        if (sup)
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/supremizer/", fieldName);
        }
        else
        {
            ITHACAstream::exportFields(modes, "./ITHACAoutput/POD/", fieldName);
        }

        // Calculate and save eigenvalue data
        eigenValues = eigenValues / eigenValues.sum();
        Eigen::VectorXd cumEigenValues = eigenValues;
        for (label j = 1; j < cumEigenValues.size(); ++j)
        {
            cumEigenValues(j) += cumEigenValues(j-1);
        }

        // Export eigenvalues
        Eigen::saveMarketVector(eigenValues,
            "./ITHACAoutput/POD/Eigenvalues_" + fieldName, para->precision, para->outytpe);
        Eigen::saveMarketVector(cumEigenValues,
            "./ITHACAoutput/POD/CumEigenvalues_" + fieldName, para->precision, para->outytpe);
    }
    else
    {
        // Read existing modes instead of computing new ones
        Info << "Reading existing modes" << endl;
        if (sup)
        {
            ITHACAstream::read_fields(modes, fieldName + "sup", "./ITHACAoutput/supremizer/");
        }
        else
        {
            ITHACAstream::read_fields(modes, fieldName, "./ITHACAoutput/POD/");
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void getMeanMemoryEfficient(
    GeometricField<Type, PatchField, GeoMesh>& templateField,
    word snapshotsPath,
    autoPtr<GeometricField<Type, PatchField, GeoMesh>>& meanField,
    bool meanex)
    {
        // Count number of snapshots in directory (excluding 0/ and constant/)
        fileName rootPath(".");
        Foam::Time runTime2(Foam::Time::controlDictName, rootPath, snapshotsPath);
        label nSnaps = runTime2.times().size() - 2;
        std::cout << "Found " << nSnaps << " time directories" << endl;
                        
        // Compute mean field
        if(!meanex)
        {
            Info << "Computing the mean of snapshots" << endl;
            // Initialize mean field to zero
            *meanField = templateField*0.;

            for(int i=0; i<nSnaps; i++)
            {
                // Read snapshot i
                GeometricField<Type, PatchField, GeoMesh> snapI = ITHACAstream::readFieldByIndex(templateField, snapshotsPath, i);

                // Sum the snapshots
                *meanField += snapI;
            }
            *meanField *= 1./nSnaps;
            ITHACAstream::exportSolution(*meanField, "ITHACAoutput", "mean");
        }
        else
        {
            Info << "Reading the mean of snapshots" << endl;
            ITHACAstream::readFieldByIndex(templateField, "ITHACAoutput/mean/", 0);
        }
    }


int main(int argc, char* argv[])
{
    // Construct the main unsteadyNS object
    unsteadyNS mainCase(argc, argv);

    // Read parameters from ITHACAdict file
    ITHACAparameters* para = ITHACAparameters::getInstance(mainCase._mesh(),
                             mainCase._runTime());

    // Retrieve mean field for velocity
    autoPtr<volVectorField> meanField_U;
    //ITHACAPOD::getMeanMemoryEfficient(mainCase._U(), "ITHACAoutput/Offline/", meanField_U, false);

    // Retrieve mean field for pressure
    autoPtr<volScalarField> meanField_p;
    //ITHACAPOD::getMeanMemoryEfficient(mainCase._p(), "ITHACAoutput/Offline/", meanField_p, false);
  
    // Perform POD on velocity field
    PtrList<volVectorField> modes_U;
    ITHACAPOD::getModesMemoryEfficient
    (
        mainCase._U(),
        "/ITHACAoutput/Offline/",
        modes_U,
        "U",
        false,
        false,
        false,
        10,
        false,
        //meanField_U
    );


/*
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
