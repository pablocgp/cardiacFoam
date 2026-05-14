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

\*---------------------------------------------------------------------------*/

#include <math.h>
#include <cmath>
#include "TNNP.H"
#include "TNNP_2004.H"
#include "addToRunTimeSelectionTable.H"
#include "ionicModel.H"
#include"ionicModelIO.H"
#include "stimulusIO.H"
#include "volFields.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(TNNP, 0);
    addToRunTimeSelectionTable
    (
        ionicModel, TNNP, dictionary
    );
}

Foam::TNNP::ProtectionSummaryEntry::ProtectionSummaryEntry()
:
    phase(""),
    variable(""),
    reason(""),
    count(0),
    firstTimeMs(0.0),
    firstTimeS(0.0),
    firstIntegrationPoint(-1),
    firstOldValue(0.0),
    firstCorrectedValue(0.0),
    minOldValue(GREAT),
    maxOldValue(-GREAT),
    maxAbsDelta(0.0)
{}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TNNP::TNNP
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    ionicModel(dict, num, initialDeltaT, solveVmWithinODESolver),
    STATES_(num),
    STATES_OLD_(num),
    CONSTANTS_(NUM_CONSTANTS, 0.0),
    ALGEBRAIC_(num),
    RATES_(num),
    numericalProtection_
    (
        dict.lookupOrDefault<Switch>("TNNPNumericalProtection", false)
    ),
    clampGates_(dict.lookupOrDefault<Switch>("TNNPClampGates", true)),
    clampVmForRates_(dict.lookupOrDefault<Switch>("TNNPClampVmForRates", true)),
    logNumericalProtection_
    (
        dict.lookupOrDefault<Switch>("TNNPLogNumericalProtection", true)
    ),
    concentrationFloor_
    (
        dict.lookupOrDefault<scalar>("TNNPConcentrationFloor", 1e-12)
    ),
    vmMinForRates_(dict.lookupOrDefault<scalar>("TNNPVMinForRates", -200.0)),
    vmMaxForRates_(dict.lookupOrDefault<scalar>("TNNPVMaxForRates", 100.0)),
    vmZeroEps_(dict.lookupOrDefault<scalar>("TNNPVoltageZeroEps", 1e-8)),
    maxProtectionLogEntries_
    (
        dict.lookupOrDefault<label>("TNNPMaxProtectionLogEntries", 1000000)
    ),
    activeIntegrationPoint_(-1),
    protectionCorrections_(0),
    protectionLogEntries_(0),
    protectionLogPtr_(),
    protectionSummary_()
{
    if (numericalProtection_)
    {
        concentrationFloor_ = max(concentrationFloor_, SMALL);
        vmZeroEps_ = max(vmZeroEps_, SMALL);
        maxProtectionLogEntries_ = max(maxProtectionLogEntries_, label(0));

        Info<< "TNNP numerical protection enabled:" << nl
            << "    concentrationFloor = " << concentrationFloor_ << nl
            << "    clampGates = " << clampGates_ << nl
            << "    clampVmForRates = " << clampVmForRates_ << nl
            << "    Vm range = [" << vmMinForRates_
            << ", " << vmMaxForRates_ << "] mV" << nl
            << "    Vm zero epsilon = " << vmZeroEps_ << " mV" << nl
            << "    logging = " << logNumericalProtection_ << nl;

        if (logNumericalProtection_)
        {
            openProtectionLog(dict);
        }
    }

    Info<< nl << "Initialize TNNP constants:" << nl;
    forAll(STATES_, i)
    {
        STATES_.set(i,      new scalarField(NUM_STATES,     0.0));
        STATES_OLD_.set(i,  new scalarField(NUM_STATES,     0.0));
        ALGEBRAIC_.set(i,   new scalarField(NUM_ALGEBRAIC,  0.0));
        RATES_.set(i,       new scalarField(NUM_STATES,     0.0));

        // 🔑 First, set tissue using base logic + overrides
        ionicModel::setTissueFromDict(); 
        TNNPinitConsts
        (
            CONSTANTS_.data(),
            RATES_[i].data(),
            STATES_[i].data(),
            tissue(),dict
        );

        if (!utilitiesMode())
        {
            stimulusIO::loadStimulusProtocol
            (
                dict, CONSTANTS_, stim_start, stim_period_S1,stim_duration,
                stim_amplitude, nstim1, stim_period_S2, nstim2
            );
        }
    }
    Info<< CONSTANTS_ << nl;

    label i0 = rand() % STATES_.size();
    Info<< "initial states:" << nl;
    Info<< STATES_[i0] << nl;
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TNNP::~TNNP()
{
    if (numericalProtection_ && logNumericalProtection_ && protectionLogPtr_)
    {
        writeProtectionSummary();
    }

    if (numericalProtection_ && protectionCorrections_ > 0)
    {
        Info<< "TNNP numerical protection applied "
            << protectionCorrections_ << " correction(s); "
            << protectionLogEntries_ << " unique protection event type"
            << (protectionLogEntries_ == 1 ? "" : "s")
            << " written" << nl;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::TNNP::openProtectionLog(const dictionary& dict) const
{
    fileName logName
    (
        dict.lookupOrDefault<fileName>
        (
            "TNNPProtectionLog",
            "TNNP_numericalProtection_summary.dat"
        )
    );

    protectionLogPtr_.reset(new OFstream(logName));
}


void Foam::TNNP::logProtection
(
    const scalar t,
    const label integrationPtI,
    const char* phase,
    const label stateI,
    const scalar oldValue,
    const scalar correctedValue,
    const char* reason
) const
{
    ++protectionCorrections_;

    if
    (
        !logNumericalProtection_
     || !protectionLogPtr_
    )
    {
        return;
    }

    const word phaseName(phase);
    const word variableName(TNNP_STATES_NAMES[stateI]);
    const word reasonName(reason);
    const scalar absDelta = mag(oldValue - correctedValue);

    label entryI = -1;
    forAll(protectionSummary_, i)
    {
        const ProtectionSummaryEntry& e = protectionSummary_[i];
        if
        (
            e.phase == phaseName
         && e.variable == variableName
         && e.reason == reasonName
        )
        {
            entryI = i;
            break;
        }
    }

    if (entryI < 0)
    {
        if (protectionSummary_.size() >= maxProtectionLogEntries_)
        {
            return;
        }

        entryI = protectionSummary_.size();
        protectionSummary_.setSize(entryI + 1);

        ProtectionSummaryEntry& e = protectionSummary_[entryI];
        e.phase = phaseName;
        e.variable = variableName;
        e.reason = reasonName;
        e.firstTimeMs = t;
        e.firstTimeS = 1e-3*t;
        e.firstIntegrationPoint = integrationPtI;
        e.firstOldValue = oldValue;
        e.firstCorrectedValue = correctedValue;
        e.minOldValue = oldValue;
        e.maxOldValue = oldValue;
        e.maxAbsDelta = absDelta;

        protectionLogEntries_ = protectionSummary_.size();
    }

    ProtectionSummaryEntry& e = protectionSummary_[entryI];
    ++e.count;
    e.minOldValue = min(e.minOldValue, oldValue);
    e.maxOldValue = max(e.maxOldValue, oldValue);
    e.maxAbsDelta = max(e.maxAbsDelta, absDelta);
}


void Foam::TNNP::writeProtectionSummary() const
{
    if (!protectionLogPtr_)
    {
        return;
    }

    OFstream& os = protectionLogPtr_();

    os  << "# total_events " << protectionCorrections_ << nl
        << "# unique_event_types " << protectionSummary_.size() << nl
        << "# event_type = phase + variable + reason" << nl
        << nl;

    os  << "# SummaryByType" << nl
        << "# phase variable reason count first_time_ms first_time_s"
        << " first_integrationPoint first_oldValue first_correctedValue"
        << " min_oldValue max_oldValue max_abs_delta" << nl;

    forAll(protectionSummary_, i)
    {
        const ProtectionSummaryEntry& e = protectionSummary_[i];
        os  << e.phase << " "
            << e.variable << " "
            << e.reason << " "
            << e.count << " "
            << e.firstTimeMs << " "
            << e.firstTimeS << " "
            << e.firstIntegrationPoint << " "
            << e.firstOldValue << " "
            << e.firstCorrectedValue << " "
            << e.minOldValue << " "
            << e.maxOldValue << " "
            << e.maxAbsDelta << nl;
    }

    os  << nl
        << "# FirstOccurrenceByType" << nl
        << "# time_ms time_s integrationPoint phase variable oldValue"
        << " correctedValue reason" << nl;

    forAll(protectionSummary_, i)
    {
        const ProtectionSummaryEntry& e = protectionSummary_[i];
        os  << e.firstTimeMs << " "
            << e.firstTimeS << " "
            << e.firstIntegrationPoint << " "
            << e.phase << " "
            << e.variable << " "
            << e.firstOldValue << " "
            << e.firstCorrectedValue << " "
            << e.reason << nl;
    }
}


void Foam::TNNP::correctStateComponent
(
    scalarField& state,
    const label stateI,
    const scalar correctedValue,
    const scalar t,
    const label integrationPtI,
    const char* phase,
    const char* reason
) const
{
    const scalar oldValue = state[stateI];

    if (oldValue != correctedValue || !std::isfinite(oldValue))
    {
        state[stateI] = correctedValue;
        logProtection
        (
            t,
            integrationPtI,
            phase,
            stateI,
            oldValue,
            correctedValue,
            reason
        );
    }
}


void Foam::TNNP::protectState
(
    scalarField& state,
    const scalar t,
    const label integrationPtI,
    const char* phase
) const
{
    if (!numericalProtection_)
    {
        return;
    }

    if (!std::isfinite(state[V]))
    {
        correctStateComponent
        (
            state,
            V,
            -86.2,
            t,
            integrationPtI,
            phase,
            "nonfinite_V"
        );
    }

    if (clampVmForRates_)
    {
        if (state[V] < vmMinForRates_)
        {
            correctStateComponent
            (
                state,
                V,
                vmMinForRates_,
                t,
                integrationPtI,
                phase,
                "Vm_below_minimum"
            );
        }
        else if (state[V] > vmMaxForRates_)
        {
            correctStateComponent
            (
                state,
                V,
                vmMaxForRates_,
                t,
                integrationPtI,
                phase,
                "Vm_above_maximum"
            );
        }
    }

    if (mag(state[V]) < vmZeroEps_)
    {
        const scalar correctedV = (state[V] < 0.0 ? -vmZeroEps_ : vmZeroEps_);
        correctStateComponent
        (
            state,
            V,
            correctedV,
            t,
            integrationPtI,
            phase,
            "avoid_i_CaL_V_zero_singularity"
        );
    }

    const label concentrationIDs[] = {K_i, Na_i, Ca_i, Ca_SR};
    const label nConcentrationIDs = sizeof(concentrationIDs)/sizeof(label);

    for (label i = 0; i < nConcentrationIDs; ++i)
    {
        const label stateI = concentrationIDs[i];

        if (!std::isfinite(state[stateI]))
        {
            correctStateComponent
            (
                state,
                stateI,
                concentrationFloor_,
                t,
                integrationPtI,
                phase,
                "nonfinite_concentration"
            );
        }
        else if (state[stateI] < concentrationFloor_)
        {
            correctStateComponent
            (
                state,
                stateI,
                concentrationFloor_,
                t,
                integrationPtI,
                phase,
                "concentration_below_floor"
            );
        }
    }

    if (clampGates_)
    {
        const label gateIDs[] = {Xr1, Xr2, Xs, m, h, j, d, f, fCa, s, r, g};
        const label nGateIDs = sizeof(gateIDs)/sizeof(label);

        for (label i = 0; i < nGateIDs; ++i)
        {
            const label stateI = gateIDs[i];

            if (!std::isfinite(state[stateI]))
            {
                correctStateComponent
                (
                    state,
                    stateI,
                    0.0,
                    t,
                    integrationPtI,
                    phase,
                    "nonfinite_gate"
                );
            }
            else if (state[stateI] < 0.0)
            {
                correctStateComponent
                (
                    state,
                    stateI,
                    0.0,
                    t,
                    integrationPtI,
                    phase,
                    "gate_below_zero"
                );
            }
            else if (state[stateI] > 1.0)
            {
                correctStateComponent
                (
                    state,
                    stateI,
                    1.0,
                    t,
                    integrationPtI,
                    phase,
                    "gate_above_one"
                );
            }
        }
    }
}

Foam::List<Foam::word> Foam::TNNP::supportedTissueTypes() const
{
    // All three tissue variants are supported in the generated code
    return {"endocardialCells", "mCells", "epicardialCells"};
}



void Foam::TNNP::calculateCurrent
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im,
    Field<Field<scalar>>& states
)
{
    const scalar tStart = stepStartTime * 1000;
    if (Im.size() != Vm.size())
    {
        FatalErrorInFunction
            << "Im.size() != Vm.size()" << abort(FatalError);
    }
    forAll(STATES_, integrationPtI)
    {
        scalarField& STATESI    = STATES_[integrationPtI];
        scalarField& ALGEBRAICI = ALGEBRAIC_[integrationPtI];
        scalarField& RATESI     = RATES_[integrationPtI];

        // Update voltage for this integration point
        STATESI[0] = Vm[integrationPtI] * 1000;
        protectState(STATESI, tStart, integrationPtI, "calculateCurrent");

        ::TNNPcomputeVariables
        (
            tStart,
            CONSTANTS_.data(),
            RATESI.data(),
            STATESI.data(),
            ALGEBRAICI.data(),
            tissue(),
            solveVmWithinODESolver()
        );
        // Jion  is the total ionic current density used by the PDE
        Im[integrationPtI] = ALGEBRAICI[Iion_cm];

        // Copy internal states to the external buffer used by coupled solvers
        // for nonlinear residual monitoring and candidate-state resets.
        if (integrationPtI < states.size() && states[integrationPtI].size() >= NUM_STATES)
        {
            states[integrationPtI] = STATESI;
        }
    }
}


//  Solve ODE with mixed singleCell implementation and 1D-3D condition
void Foam::TNNP::solveODE
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im,
    Field<Field<scalar>>& states
)
{
    const scalar tStart = stepStartTime * 1000;
    const scalar tEnd   = (stepStartTime + deltaT) * 1000;
    const label monitorCell = 0;

    forAll(STATES_, integrationPtI)
    {
        scalarField& STATESI    = STATES_[integrationPtI];
        scalarField& ALGEBRAICI = ALGEBRAIC_[integrationPtI];
        scalarField& RATESI     = RATES_[integrationPtI];

        scalar& step = ionicModel::step()[integrationPtI];

        // If Vm is solved by the PDE, feed that Vm (in mV) into the cell model
        if (!solveVmWithinODESolver())
        {
            STATESI[0] = Vm[integrationPtI]*1000.0;
        }
        activeIntegrationPoint_ = integrationPtI;
        protectState(STATESI, tStart, integrationPtI, "solveODE_start");

        // Clamp time step (ms)
        step = min(step, deltaT * 1000.0);
        // Advance ODE system for all states
        odeSolver().solve(tStart, tEnd, STATESI, step);
        protectState(STATESI, tEnd, integrationPtI, "solveODE_end");

        // Update algebraics and rates at tEnd (includes Iion and I_stim)
        ::TNNPcomputeVariables
        (
            tEnd,
            CONSTANTS_.data(),
            RATESI.data(),
            STATESI.data(),
            ALGEBRAICI.data(),
            tissue(),
            solveVmWithinODESolver()
        );
        ::TNNPcomputeRates
        (
            tEnd,
            CONSTANTS_.data(),
            RATESI.data(),
            STATESI.data(),
            ALGEBRAICI.data(),
            tissue(),
            solveVmWithinODESolver()
        );
        if (integrationPtI == monitorCell)
        {debugPrintFields(integrationPtI, tStart, tEnd, step);}

        // Total ionic current density used by PDE
        Im[integrationPtI] = ALGEBRAICI[Iion_cm] ;

        // Copy internal states to the external buffer used by coupled solvers
        // for nonlinear residual monitoring and candidate-state resets.
        if (integrationPtI < states.size() && states[integrationPtI].size() >= NUM_STATES)
        {
            states[integrationPtI] = STATESI;
        }
    }

    activeIntegrationPoint_ = -1;
}

void Foam::TNNP::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    // Must match NUM_ALGEBRAIC from the generated TNNP code
    scalarField ALGEBRAIC_TMP(NUM_ALGEBRAIC, 0.0);

    if (numericalProtection_)
    {
        scalarField yProtected(y);
        protectState(yProtected, t, activeIntegrationPoint_, "derivatives");

        ::TNNPcomputeRates
        (
            t,
            CONSTANTS_.data(),
            dydt.data(),                              // RATES (output)
            yProtected.data(),                        // STATES (protected input)
            ALGEBRAIC_TMP.data(),                     // ALGEBRAIC (scratch)
            tissue(),
            solveVmWithinODESolver()
        );
    }
    else
    {
        ::TNNPcomputeRates
        (
            t,
            CONSTANTS_.data(),
            dydt.data(),                              // RATES (output)
            const_cast<scalarField&>(y).data(),       // STATES (input)
            ALGEBRAIC_TMP.data(),                     // ALGEBRAIC (scratch)
            tissue(),
            solveVmWithinODESolver()
        );
    }
}

void Foam::TNNP::updateStatesOld(const Field<Field<scalar>>&) const
{
    saveStateSnapshot(STATES_, STATES_OLD_);
}

void Foam::TNNP::resetStatesToStatesOld(Field<Field<scalar>>&) const
{
    restoreStateSnapshot(STATES_, STATES_OLD_);
}

// ------------------------------------------------------------------------- //
//  Writing logic in singleCell and 3D simulations

//Writing functions for singleCell implementation
Foam::wordList Foam::TNNP::exportedFieldNames() const
    {
        return ionicModelIO::exportedFieldNames
        (
            variableExport_,
            TNNP_STATES_NAMES, NUM_STATES,
            TNNP_ALGEBRAIC_NAMES, NUM_ALGEBRAIC
        );
    }

    Foam::wordList Foam::TNNP::debugPrintedNames() const
    {
        return ionicModelIO::exportedFieldNames
        (
            debugVarNames_,
            TNNP_STATES_NAMES, NUM_STATES,
            TNNP_ALGEBRAIC_NAMES, NUM_ALGEBRAIC
        );
    }

void Foam::TNNP::exportStates
(
    const Field<Field<scalar>>&,
    PtrList<volScalarField>& outFields
)
{
    ionicModelIO::exportStateFields
    (
        STATES_,ALGEBRAIC_,
        exportedFieldNames(),
        TNNP_STATES_NAMES,NUM_STATES,
        TNNP_ALGEBRAIC_NAMES,NUM_ALGEBRAIC,
        outFields
    );
}

void Foam::TNNP::exportStatesIntegrationPoints
(
    const Field<Field<scalar>>&,
    PtrList<volScalarField>& outFields,
    const CompactListList<scalar>& cellIionQuadW
)
{
    ionicModelIO::exportStateFieldsIntegrationPoints
    (
        STATES_,ALGEBRAIC_,
        exportedFieldNames(),
        TNNP_STATES_NAMES,NUM_STATES,
        TNNP_ALGEBRAIC_NAMES,NUM_ALGEBRAIC,
        outFields,
        cellIionQuadW
    );
}

void Foam::TNNP::debugPrintFields
(
    label cellI,
    scalar t1,
    scalar t2,
    scalar step
) const
{
    ionicModelIO::debugPrintFields
    (
        STATES_, ALGEBRAIC_,
        debugPrintedNames(),
        TNNP_STATES_NAMES, NUM_STATES,
        TNNP_ALGEBRAIC_NAMES, NUM_ALGEBRAIC,
        cellI,t1,t2,step
    );
}




void Foam::TNNP::writeHeader(OFstream& os) const
{
    const wordList names = exportedFieldNames();

    if (!names.empty())
    {
        ionicModelIO::writeSelectedHeader(os, names);
    }
    else
    {
        ionicModelIO::writeHeader
        (
            os,
            TNNP_STATES_NAMES, NUM_STATES,
            TNNP_ALGEBRAIC_NAMES, NUM_ALGEBRAIC
        );
    }
}

static Foam::scalar TNNP_vm(const Foam::scalarField& S)
{
    return S[0];
}

void Foam::TNNP::write(const scalar t, OFstream& os) const
{
    const wordList names = exportedFieldNames();

    if (!names.empty())
    {
        ionicModelIO::writeSelected
        (
            t, os,
            STATES_, ALGEBRAIC_,
            names,
            TNNP_STATES_NAMES, NUM_STATES,
            TNNP_ALGEBRAIC_NAMES, NUM_ALGEBRAIC
        );
    }
    else
    {
        ionicModelIO::write
        (
            t,os,
            STATES_,ALGEBRAIC_,
            RATES_,TNNP_vm
        );
    }
}


void Foam::TNNP::sweepCurrent
(
    const word& currentName,
    scalar Vmin,
    scalar Vmax,
    label nPts,
    const fileName& outputFile
) const
{
    // Retrieve dependency variables
    const auto& depMap = TNNPDependencyMap();

    if (!depMap.found(currentName))
    {
        FatalErrorInFunction
            << "Unknown current: " << currentName << nl
            << "Available currents: " << depMap.toc() << nl
            << exit(FatalError);
    }

    const wordList& deps = depMap[currentName];
    OFstream os(outputFile);
    // Write sweep header: V,<deps...>
    ionicModelIO::writeSweepHeader(os, deps);

    // Working arrays from integration point 0
    scalarField STATESI = STATES_[0];
    scalarField RATESI(NUM_STATES, 0.0);
    scalarField ALGI(NUM_ALGEBRAIC, 0.0);

    // Voltage sweep
    for (label i = 0; i < nPts; ++i)
    {
        scalar V = Vmin + (Vmax - Vmin) * scalar(i) / (nPts - 1);

        // Reset all states to baseline
        STATESI = STATES_[0];

        // Overwrite membrane voltage (dimensionless in BO2008)
        STATESI[0] = V;

        ::TNNPcomputeVariables
        (
            0.0,                       // VOI
            CONSTANTS_.data(),
            RATESI.data(),
            STATESI.data(),
            ALGI.data(),
            tissue(),
            solveVmWithinODESolver()
        );

        ionicModelIO::writeOneSweepRow
        (
            os, V, deps,STATESI,ALGI,
            TNNP_STATES_NAMES, NUM_STATES,
            TNNP_ALGEBRAIC_NAMES, NUM_ALGEBRAIC
        );
    }
    Info<< "Sweep for " << currentName
        << " written to " << outputFile << nl;
}
