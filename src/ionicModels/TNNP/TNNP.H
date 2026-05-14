/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::TNNP

Description
    Wrapper class for calling from the cellML ionic model for the ten Tusscher,
    Noble, Noble, Panfilov, 2004, model.

SourceFiles
    TNNP.C

Author
    Philip Cardiff, UCD.
    Simão Nieto de Castro, UCD.

\*---------------------------------------------------------------------------*/

#ifndef TNNP_H
#define TNNP_H

#include "ionicModel.H"
#include "OFstream.H"
#include "volFields.H"
#include "TNNP_2004Names.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class TNNP Declaration
\*---------------------------------------------------------------------------*/

class TNNP
:
    public ionicModel
{
    // Private Data

        struct ProtectionSummaryEntry
        {
            word phase;
            word variable;
            word reason;
            label count;
            scalar firstTimeMs;
            scalar firstTimeS;
            label firstIntegrationPoint;
            scalar firstOldValue;
            scalar firstCorrectedValue;
            scalar minOldValue;
            scalar maxOldValue;
            scalar maxAbsDelta;

            ProtectionSummaryEntry();
        };

        //- States
        mutable PtrList<scalarField> STATES_;

        //- States - old time
        mutable PtrList<scalarField> STATES_OLD_;

        //- Constants
        mutable scalarField CONSTANTS_;

        //- Algebraic
        mutable PtrList<scalarField> ALGEBRAIC_;

        //- Rates
        mutable PtrList<scalarField> RATES_;

        //- Enable numerical protection for the generated TNNP algebra
        Switch numericalProtection_;

        //- Clamp gate variables to [0, 1] when numericalProtection_ is active
        Switch clampGates_;

        //- Clamp Vm in mV before evaluating the generated rates
        Switch clampVmForRates_;

        //- Write a log entry whenever a protection correction is applied
        Switch logNumericalProtection_;

        //- Floor for positive concentrations used inside log/pow divisions
        scalar concentrationFloor_;

        //- Lower Vm bound used by the rate protection [mV]
        scalar vmMinForRates_;

        //- Upper Vm bound used by the rate protection [mV]
        scalar vmMaxForRates_;

        //- Small Vm offset used to avoid the i_CaL V=0 removable singularity [mV]
        scalar vmZeroEps_;

        //- Maximum number of unique protection event types written to disk
        label maxProtectionLogEntries_;

        //- Integration point currently being processed by the ODE solver
        mutable label activeIntegrationPoint_;

        //- Number of applied numerical corrections
        mutable label protectionCorrections_;

        //- Number of unique protection event types written to disk
        mutable label protectionLogEntries_;

        //- Optional protection summary stream
        mutable autoPtr<OFstream> protectionLogPtr_;

        //- Protection summary entries, grouped by phase/variable/reason
        mutable List<ProtectionSummaryEntry> protectionSummary_;


    // Private Member Functions

        //- No copy construct
        TNNP(const TNNP&) = delete;

        //- No copy assignment
        void operator=(const TNNP&) = delete;

        //- Open the optional numerical protection log
        void openProtectionLog(const dictionary& dict) const;

        //- Log one numerical protection correction in the compact summary
        void logProtection
        (
            const scalar t,
            const label integrationPtI,
            const char* phase,
            const label stateI,
            const scalar oldValue,
            const scalar correctedValue,
            const char* reason
        ) const;

        //- Write compact numerical protection summary
        void writeProtectionSummary() const;

        //- Correct a state component and record the correction
        void correctStateComponent
        (
            scalarField& state,
            const label stateI,
            const scalar correctedValue,
            const scalar t,
            const label integrationPtI,
            const char* phase,
            const char* reason
        ) const;

        //- Apply conservative guards before evaluating the generated TNNP algebra
        void protectState
        (
            scalarField& state,
            const scalar t,
            const label integrationPtI,
            const char* phase
        ) const;


        public:

        //- Runtime type information
        TypeName("TNNP");


        //override the only available tissue in this model
        virtual List<word> supportedTissueTypes() const override;
        virtual List<word> supportedDimensions() const override
        {
            return {}; //No dimensions flag
        }



        // Constructors

            //- Construct given a dictionary, the number of integration points and a
            //  initial time step for the ODE solver
            TNNP
            (
                const dictionary& dict,
                const label nIntegrationPoints,
                const scalar initialDeltaT,
                const Switch solveVmWithinODESolver = false
            );


        //- Destructor
        virtual ~TNNP();


        // Member Functions

            virtual void writeHeader(OFstream& output) const override;
            //- Write fields to the given file given the current time
            virtual void write(const scalar t, OFstream& output) const override;




            virtual void calculateCurrent
            (
                const scalar stepStartTime,
                const scalar deltaT,
                const scalarField& Vm,
                scalarField& Im,
                Field<Field<scalar>>& states
            ) override;


            virtual void solveODE
            (
                const scalar stepStartTime,
                const scalar deltaT,
                const scalarField& Vm,
                scalarField& Im,
                Field<Field<scalar>>& states
            ) override;


            //- Number of ODE's to solve
            virtual label nEqns() const override
            {
                return NUM_STATES;
            }

            //- Derivatives of each ODE with respect to the primary unknown in that
            //  equation
            virtual void derivatives
            (
                const scalar t,
                const scalarField& y,
                scalarField& dydt
            ) const override;


            virtual void updateStatesOld(const Field<Field<scalar>>& states) const;
            virtual void resetStatesToStatesOld(Field<Field<scalar>>& states) const;

                    // Override to provide variable names for export
            virtual Foam::wordList exportedFieldNames() const override;
            virtual wordList debugPrintedNames() const override;


            // Override to export the variables
            virtual void exportStates
            (
                const Field<Field<scalar>>&,
                PtrList<volScalarField>&
            ) override;

            virtual void exportStatesIntegrationPoints
            (
                const Field<Field<scalar>>&,
                PtrList<volScalarField>&,
                const CompactListList<scalar>&
            ) override;

            // Override to print variables in solver
            void debugPrintFields
            (
                label cellI,
                scalar t1,
                scalar t2,
                scalar step
            ) const override;

            virtual void sweepCurrent
            (
                const word& currentName,
                scalar Vmin,
                scalar Vmax,
                label nPts,
                const fileName& outputFile
            ) const override;


};

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
