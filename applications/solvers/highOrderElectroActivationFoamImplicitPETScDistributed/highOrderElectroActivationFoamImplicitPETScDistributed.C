/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

Solver
    highOrderElectroActivationFoamImplicitPETSc

Description
    Implicit monodomain solver for cardiac electrophysiology with Jacobian-Free
    Newton-Krylov (JFNK) nonlinear coupling.

    PDE solved (monodomain):

        chi * Cm * dVm/dt  =  div( D . grad Vm )  -  Iion(Vm, s)  +  Iext
                  ds/dt    =  f_ionic(Vm, s)

    where Vm is the transmembrane potential, D the conductivity tensor, chi
    the surface-to-volume ratio, Cm the membrane capacitance, Iion the
    (nonlinear) ionic current, s the vector of ionic states, and Iext an
    external stimulus current.

    Numerical features:
      * Time integration: theta-method (Backward Euler or Crank-Nicolson) or
        BDF1-4, with ESDIRK3 for the steps a multistep method cannot take
        (startup, restart at a stimulus change, a step of non-nominal
        length). The ionic states are integrated with Vm driven across the
        step (linear ramp, or the BDF/ESDIRK history interpolant).
      * Mass matrix: lumped (diagonal) or LRE-consistent high-order.
      * Diffusion operator: standard FV orthogonal stencil or LRE high-order
        face-quadrature reconstruction (Taylor expansion up to 3rd order).
      * Ionic source: optionally evaluated at LRE cell quadrature points to
        match the spatial accuracy of the diffusion operator.
      * Nonlinear solve: JFNK (matrix-free Newton + restarted GMRES with
        modified Gram-Schmidt) or Picard / diagonal-linearised variants.
      * Linear solve (non-JFNK branches): Eigen SparseLU or BiCGSTAB.
      * Robustness: if the nonlinear iteration fails to converge, Vm, Iion
        and ionic states are rolled back to the values at the beginning of
        the time step (configurable via nonlinearAcceptUnconverged).

    Diagnostics:
      * activationTime: first crossing of activationThreshold (default 0 V),
        recorded via linear interpolation across the time step.
      * P8 point and P1-P8 diagonal activation samples (Niederer 2012
        benchmark) written to postProcessing.
      * Per-time-step nonlinear residuals written to nonlinearResiduals.dat.

\*---------------------------------------------------------------------------*/

#include <petscksp.h>
#include "fvCFD.H"
// Modified for cardiacFoam: LRE (serial) -> highOrderInterp, a thin
// adapter over solids4foam's parallel movingLeastSquares. The API
// mapping, and the two places where semantics changed (quadrature
// weights, stencil storage), are documented in
// src/highOrderAdapter/highOrderInterp.H.
#include "highOrderInterp.H"
#include "globalIndex.H"
#include "syncTools.H"
// Added for cardiacFoam: needed by the coupled-patch (processor face)
// stabilisation exchange in the stiffness assemblies.
#include "PstreamBuffers.H"
#include "processorFvPatch.H"
#include "Field.H"
#include "volFields.H"
#include <cmath>
#include <chrono>
#include <algorithm>
#include <cctype>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>
#include <functional>
#include <sys/resource.h>
#include <string>
#include <vector>
#ifdef __GLIBC__
#include <malloc.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{

// Shared code of the two monodomain solvers (src/monodomainCommon). Included
// here, inside the anonymous namespace, so it has internal linkage.
#include "monodomainCommon.H"




    // ----------------------------------------------------------------- //
    // Embedded ten Tusscher--Noble--Noble--Panfilov 2004 ionic model.
    //
    // The original electro solver used cardiacFoam's runtime-selected
    // ionicModel/TNNP classes from src/.  This solver keeps the same TNNP
    // equations locally so the numerical PDE/source coupling can be kept
    // aligned with the manufactured JFNK solver without linking against
    // cardiacFoam-pc/src.  Units follow the generated CellML code:
    // Vm is stored internally in mV and time in ms.  Since 1 mV/ms = 1 V/s,
    // Iion_cm can be inserted directly into the monodomain source used by
    // the PDE, whose Vm field is in V.
    // ----------------------------------------------------------------- //

enum CONSTANTS_INDEX{
    R,              // 0  : gas constant
    T,              // 1  : temperature
    F,              // 2  : Faraday constant
    Cm,             // 3  : membrane capacitance (pF)
    V_c,            // 4  : cell volume (um^3)

    // Stimulus (S1)
    stim_start,         // 5
    stim_period_S1,     // 6
    stim_duration,      // 7
    stim_amplitude,     // 8

    // Reversal potentials & ion concentrations
    P_kna,          // 9
    K_o,            // 10
    Na_o,           // 11
    Ca_o,           // 12

    // Conductances
    g_K1,           // 13
    g_Kr,           // 14
    g_Ks,           // 15
    g_Na,           // 16
    g_bNa,          // 17
    g_CaL,          // 18
    g_bCa,          // 19
    g_to,           // 20

    // Na/K Pump
    P_NaK,          // 21
    K_mk,           // 22
    K_mNa,          // 23

    // NCX exchanger
    K_NaCa,         // 24
    K_sat,          // 25
    alpha_NCX,      // 26
    gamma_NCX,      // 27
    Km_Ca,          // 28
    Km_Nai,         // 29

    // Calcium pump + potassium pump
    g_pCa,          // 30
    K_pCa,          // 31
    g_pK,           // 32

    // Calcium dynamics (SR release / uptake)
    tau_g,          // 33
    a_rel,          // 34
    b_rel,          // 35
    c_rel,          // 36
    K_up,           // 37
    V_leak,         // 38
    Vmax_up,        // 39

    // Buffers
    Buf_c,          // 40
    K_buf_c,        // 41
    Buf_sr,         // 42
    K_buf_sr,       // 43
    V_sr,           // 44
    tau_fCa,        // 45

    // Added multi-stimulus parameters (S1/S2 protocol)
    nstim1,         // 46
    stim_period_S2, // 47
    nstim2,         // 48

    NUM_CONSTANTS
};


enum STATES_INDEX {
    V,      // membrane voltage
    K_i,    // intracellular potassium
    Na_i,   // intracellular sodium
    Ca_i,   // intracellular calcium
    Xr1,
    Xr2,
    Xs,
    m,
    h,
    j,
    d,
    f,
    fCa,
    s,
    r,
    Ca_SR,
    g,
    NUM_STATES
};


enum ALGEBRAIC_INDEX {
    Istim,
    xr1_inf,
    xr2_inf,
    xs_inf,
    m_inf,
    h_inf,
    j_inf,
    d_inf,
    f_inf,
    alpha_fCa,
    s_inf,
    r_inf,
    g_inf,
    E_Na,
    alpha_xr1,
    alpha_xr2,
    alpha_xs,
    alpha_m,
    alpha_h,
    alpha_j,
    alpha_d,
    tau_f,
    beta_fCa,
    tau_s,
    tau_r,
    d_g,
    E_K,
    beta_xr1,
    beta_xr2,
    beta_xs,
    beta_m,
    beta_h,
    beta_j,
    gama_fCa,
    beta_d,
    E_Ks,
    tau_xr1,
    tau_xr2,
    tau_xs,
    tau_m,
    tau_h,
    tau_j,
    gamma_d,
    fCa_inf,
    E_Ca,
    tau_d,
    d_fCa,
    alpha_K1,
    beta_K1,
    xK1_inf,
    i_K1,
    i_Kr,
    i_Ks,
    i_Na,
    i_b_Na,
    i_CaL,
    i_b_Ca,
    i_to,
    i_NaK,
    i_NaCa,
    i_p_Ca,
    i_p_K,
    i_rel,
    i_up,
    i_leak,
    ddt_Ca_i_total,
    ddt_Ca_sr_total,
    f_JCa_i_free,
    f_JCa_sr_free,
    Iion_cm,

    NUM_ALGEBRAIC
};

static const char* TNNP_STATES_NAMES[17] = {
    "V",        // 0: membrane voltage (millivolt)
    "K_i",      // 1: intracellular potassium (millimolar)
    "Na_i",     // 2: intracellular sodium (millimolar)
    "Ca_i",     // 3: intracellular calcium (millimolar)
    "Xr1",      // 4: Xr1 gate (dimensionless)
    "Xr2",      // 5: Xr2 gate (dimensionless)
    "Xs",       // 6: Xs gate (dimensionless)
    "m",        // 7: m gate (dimensionless)
    "h",        // 8: h gate (dimensionless)
    "j",        // 9: j gate (dimensionless)
    "d",        // 10: d gate (dimensionless)
    "f",        // 11: f gate (dimensionless)
    "fCa",      // 12: fCa gate (dimensionless)
    "s",        // 13: s gate (dimensionless)
    "r",        // 14: r gate (dimensionless)
    "Ca_SR",    // 15: sarcoplasmic reticulum calcium (millimolar)
    "g"         // 16: g gate (dimensionless)
};
static const char* TNNP_ALGEBRAIC_NAMES[70] = {
    "Istim",                    // 0
    "xr1_inf",                  // 1
    "xr2_inf",                  // 2
    "xs_inf",                   // 3
    "m_inf",                    // 4
    "h_inf",                    // 5
    "j_inf",                    // 6
    "d_inf",                    // 7
    "f_inf",                    // 8
    "alpha_fCa",               // 9
    "s_inf",                    // 10
    "r_inf",                    // 11
    "g_inf",                    // 12
    "E_Na",                     // 13
    "alpha_xr1",               // 14
    "alpha_xr2",               // 15
    "alpha_xs",                // 16
    "alpha_m",                 // 17
    "alpha_h",                 // 18
    "alpha_j",                 // 19
    "alpha_d",                 // 20
    "tau_f",                    // 21
    "beta_fCa",                // 22
    "tau_s",                    // 23
    "tau_r",                    // 24
    "d_g",                      // 25
    "E_K",                      // 26
    "beta_xr1",                // 27
    "beta_xr2",                // 28
    "beta_xs",                 // 29
    "beta_m",                  // 30
    "beta_h",                  // 31
    "beta_j",                  // 32
    "gama_fCa",                // 33
    "beta_d",                  // 34
    "E_Ks",                     // 35
    "tau_xr1",                 // 36
    "tau_xr2",                 // 37
    "tau_xs",                  // 38
    "tau_m",                   // 39
    "tau_h",                   // 40
    "tau_j",                   // 41
    "gamma_d",                 // 42
    "fCa_inf",                 // 43
    "E_Ca",                     // 44
    "tau_d",                   // 45
    "d_fCa",                   // 46
    "alpha_K1",                // 47
    "beta_K1",                 // 48
    "xK1_inf",                 // 49
    "i_K1",                     // 50
    "i_Kr",                     // 51
    "i_Ks",                     // 52
    "i_Na",                     // 53
    "i_b_Na",                   // 54
    "i_CaL",                    // 55
    "i_b_Ca",                   // 56
    "i_to",                     // 57
    "i_NaK",                    // 58
    "i_NaCa",                   // 59
    "i_p_Ca",                   // 60
    "i_p_K",                    // 61
    "i_rel",                    // 62
    "i_up",                     // 63
    "i_leak",                   // 64
    "ddt_Ca_i_total",           // 65
    "ddt_Ca_sr_total",          // 66
    "f_JCa_i_free",             // 67
    "f_JCa_sr_free",            // 68
    "Iion_cm"                   // 69
};

    // S1/S2 pacing waveform of the embedded TNNP kernel: returns the stimulus
    // current at time VOI as a train of nStim1 pulses of period stimPeriodS1
    // followed by nStim2 of period stimPeriodS2, each of stimDuration and
    // amplitude stimAmplitude. This is the CELL-level stimulus of the ionic
    // model, distinct from applyStimulus() below, which is the TISSUE-level one
    // applied to a region of the mesh.
    scalar embeddedComputeStimulus
    (
        const scalar VOI,
        const scalar stimStart,
        const scalar stimPeriodS1,
        const scalar stimDuration,
        const scalar stimAmplitude,
        const label nStim1,
        const scalar stimPeriodS2,
        const label nStim2
    )
    {
        scalar Istim = 0.0;

        if (stimPeriodS1 > 0.0 && nStim1 > 0)
        {
            const scalar tp = VOI - stimStart;
            if (tp >= 0.0 && tp <= stimPeriodS1*nStim1)
            {
                const scalar phase = tp - std::floor(tp/stimPeriodS1)*stimPeriodS1;
                if (phase >= 0.0 && phase <= stimDuration)
                {
                    Istim = -stimAmplitude;
                }
            }
        }

        if (Istim == 0.0 && stimPeriodS2 > 0.0 && nStim2 > 0)
        {
            const scalar tS1End = stimStart + stimPeriodS1*nStim1;
            const scalar tp = VOI - tS1End;
            if (tp >= 0.0 && tp <= stimPeriodS2*nStim2)
            {
                const scalar phase = tp - std::floor(tp/stimPeriodS2)*stimPeriodS2;
                if (phase >= 0.0 && phase <= stimDuration)
                {
                    Istim = -stimAmplitude;
                }
            }
        }

        return Istim;
    }

    inline Foam::scalar computeIstim(Foam::scalar t, const double* C)
    {
        return embeddedComputeStimulus
        (
            t,
            C[stim_start],
            C[stim_period_S1],
            C[stim_duration],
            C[stim_amplitude],
            Foam::label(C[nstim1]),
            C[stim_period_S2],
            Foam::label(C[nstim2])
        );
    }

void
TNNPinitConsts(double* CONSTANTS, double* RATES, double *STATES, int tissueFlag, const Foam::dictionary& stimulus)

{
STATES[V] = -86.2;
CONSTANTS[R] = 8.314;
CONSTANTS[T] = 310;
CONSTANTS[F] = 96.485;
CONSTANTS[Cm] = 185;
CONSTANTS[V_c] = 16404;


CONSTANTS[P_kna] = 0.03;
CONSTANTS[K_o] = 5.4;
CONSTANTS[Na_o] = 140;
STATES[K_i] = 138.3;
STATES[Na_i] = 11.6;
CONSTANTS[Ca_o] = 2;
STATES[Ca_i] = 0.0002;
CONSTANTS[g_K1] = 5.405;
CONSTANTS[g_Kr] = 0.096;
STATES[Xr1] = 0;
STATES[Xr2] = 1;
//Conditional tissue conductance Gks
CONSTANTS[g_Ks] = (tissueFlag == 2)
          ? 0.062
          : 0.245;

STATES[Xs] = 0;
CONSTANTS[g_Na] = 14.838;
STATES[m] = 0;
STATES[h] = 0.75;
STATES[j] = 0.75;
CONSTANTS[g_bNa] = 0.00029;
CONSTANTS[g_CaL] = 0.175;
STATES[d] = 0;
STATES[f] = 1;
STATES[fCa] = 1;
CONSTANTS[g_bCa] = 0.000592;
//Conditional tissue conductance Gto
CONSTANTS[g_to] = (tissueFlag == 3)
          ? 0.073
          : 0.294;

STATES[s] = 1;
STATES[r] = 0;
CONSTANTS[P_NaK] = 1.362;
CONSTANTS[K_mk] = 1;
CONSTANTS[K_mNa] = 40;
CONSTANTS[K_NaCa] = 1000;
CONSTANTS[K_sat] = 0.1;
CONSTANTS[alpha_NCX] = 2.5;
CONSTANTS[gamma_NCX] = 0.35;
CONSTANTS[Km_Ca] = 1.38;
CONSTANTS[Km_Nai] = 87.5;
CONSTANTS[g_pCa] = 0.825;
CONSTANTS[K_pCa] = 0.0005;
CONSTANTS[g_pK] = 0.0146;
STATES[Ca_SR] = 0.2;
STATES[g] = 1;
CONSTANTS[tau_g] = 2;
CONSTANTS[a_rel] = 0.016464;
CONSTANTS[b_rel] = 0.25;
CONSTANTS[c_rel] = 0.008232;
CONSTANTS[K_up] = 0.00025;
CONSTANTS[V_leak] = 8e-5;
CONSTANTS[Vmax_up] = 0.000425;
CONSTANTS[Buf_c] = 0.15;
CONSTANTS[K_buf_c] = 0.001;
CONSTANTS[Buf_sr] = 10;
CONSTANTS[K_buf_sr] = 0.3;
CONSTANTS[V_sr] = 1094;
CONSTANTS[tau_fCa] = 2.00000;

}

void
TNNPcomputeRates(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC,int tissueFlag, bool solveVmWithinODESolver)
{

//INa
ALGEBRAIC[m_inf] = 1.00000/std::pow(1.00000+std::exp((- 56.8600 - STATES[V])/9.03000), 2.00000);
ALGEBRAIC[alpha_m] = 1.00000/(1.00000+std::exp((- 60.0000 - STATES[V])/5.00000));
ALGEBRAIC[beta_m] = 0.100000/(1.00000+std::exp((STATES[V]+35.0000)/5.00000))+0.100000/(1.00000+std::exp((STATES[V] - 50.0000)/200.000));
ALGEBRAIC[tau_m] =  1.00000*ALGEBRAIC[alpha_m]*ALGEBRAIC[beta_m];
RATES[m] = (ALGEBRAIC[m_inf] - STATES[m])/ALGEBRAIC[tau_m];
ALGEBRAIC[h_inf] = 1.00000/std::pow(1.00000+std::exp((STATES[V]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[alpha_h] = (STATES[V]<- 40.0000 ?  0.0570000*std::exp(- (STATES[V]+80.0000)/6.80000) : 0.00000);
ALGEBRAIC[beta_h] = (STATES[V]<- 40.0000 ?  2.70000*std::exp( 0.0790000*STATES[V])+ 310000.*std::exp( 0.348500*STATES[V]) : 0.770000/( 0.130000*(1.00000+std::exp((STATES[V]+10.6600)/- 11.1000))));
ALGEBRAIC[tau_h] = 1.00000/(ALGEBRAIC[alpha_h]+ALGEBRAIC[beta_h]);
RATES[h] = (ALGEBRAIC[h_inf] - STATES[h])/ALGEBRAIC[tau_h];
ALGEBRAIC[j_inf] = 1.00000/std::pow(1.00000+std::exp((STATES[V]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[alpha_j] = (STATES[V]<- 40.0000 ? (( ( - 25428.0*std::exp( 0.244400*STATES[V]) -  6.94800e-06*std::exp( - 0.0439100*STATES[V]))*(STATES[V]+37.7800))/1.00000)/(1.00000+std::exp( 0.311000*(STATES[V]+79.2300))) : 0.00000);
ALGEBRAIC[beta_j] = (STATES[V]<- 40.0000 ? ( 0.0242400*std::exp( - 0.0105200*STATES[V]))/(1.00000+std::exp( - 0.137800*(STATES[V]+40.1400))) : ( 0.600000*std::exp( 0.0570000*STATES[V]))/(1.00000+std::exp( - 0.100000*(STATES[V]+32.0000))));
ALGEBRAIC[tau_j] = 1.00000/(ALGEBRAIC[alpha_j]+ALGEBRAIC[beta_j]);
RATES[j] = (ALGEBRAIC[j_inf] - STATES[j])/ALGEBRAIC[tau_j];

//ICaL
ALGEBRAIC[d_inf] = 1.00000/(1.00000+std::exp((- 5.00000 - STATES[V])/7.50000));
ALGEBRAIC[alpha_d] = 1.40000/(1.00000+std::exp((- 35.0000 - STATES[V])/13.0000))+0.250000;
ALGEBRAIC[gama_fCa] = 1.40000/(1.00000+std::exp((STATES[V]+5.00000)/5.00000));
ALGEBRAIC[gamma_d] = 1.00000/(1.00000+std::exp((50.0000 - STATES[V])/20.0000));
ALGEBRAIC[tau_d] =  1.00000*ALGEBRAIC[alpha_d]*ALGEBRAIC[gama_fCa]+ALGEBRAIC[gamma_d];
RATES[d] = (ALGEBRAIC[d_inf] - STATES[d])/ALGEBRAIC[tau_d];
ALGEBRAIC[f_inf] = 1.00000/(1.00000+std::exp((STATES[V]+20.0000)/7.00000));
ALGEBRAIC[tau_f] =  1125.00*std::exp(- std::pow(STATES[V]+27.0000, 2.00000)/240.000)+80.0000+165.000/(1.00000+std::exp((25.0000 - STATES[V])/10.0000));
RATES[f] = (ALGEBRAIC[f_inf] - STATES[f])/ALGEBRAIC[tau_f];

ALGEBRAIC[alpha_fCa] = 1.00000/(1.00000+std::pow(STATES[Ca_i]/0.000325000, 8.00000));
ALGEBRAIC[beta_fCa] = 0.100000/(1.00000+std::exp((STATES[Ca_i] - 0.000500000)/0.000100000));
ALGEBRAIC[beta_d] = 0.200000/(1.00000+std::exp((STATES[Ca_i] - 0.000750000)/0.000800000));
ALGEBRAIC[fCa_inf] = (ALGEBRAIC[alpha_fCa]+ALGEBRAIC[beta_fCa]+ALGEBRAIC[beta_d]+0.230000)/1.46000;
ALGEBRAIC[d_fCa] = (ALGEBRAIC[fCa_inf] - STATES[fCa])/CONSTANTS[tau_fCa];
RATES[fCa] = (ALGEBRAIC[fCa_inf]>STATES[fCa]&&STATES[V]>- 60.0000 ? 0.00000 : ALGEBRAIC[d_fCa]);

//Ito
//Conditional tissue flag for the Ito current
ALGEBRAIC[s_inf] = (tissueFlag == 1)
    ? 1.00000 / (1.00000 + std::exp((STATES[V] + 20.0000) / 5.00000))
    : 1.10000 / (1.00000 + std::exp((STATES[V] + 28.0000) / 6.00000));
ALGEBRAIC[tau_s] = (tissueFlag == 1)
     ? 85.0000*std::exp(- std::pow(STATES[V]+45.0000, 2.00000)/320.000)+5.00000/(1.00000+std::exp((STATES[V] - 20.0000)/5.00000))+3.00000
     : 1000.0000*std::exp(- std::pow(STATES[V]+67.0000, 2.00000)/1000.000)+8.00000;
RATES[s] = (ALGEBRAIC[s_inf] - STATES[s])/ALGEBRAIC[tau_s];
ALGEBRAIC[r_inf] = 1.00000/(1.00000+std::exp((20.0000 - STATES[V])/6.00000));
ALGEBRAIC[tau_r] =  9.50000*std::exp(- std::pow(STATES[V]+40.0000, 2.00000)/1800.00)+0.800000;
RATES[r] = (ALGEBRAIC[r_inf] - STATES[r])/ALGEBRAIC[tau_r];

//IKr
ALGEBRAIC[xr1_inf] = 1.00000/(1.00000+std::exp((- 26.0000 - STATES[V])/7.00000));
ALGEBRAIC[alpha_xr1] = 450.000/(1.00000+std::exp((- 45.0000 - STATES[V])/10.0000));
ALGEBRAIC[beta_xr1] = 6.00000/(1.00000+std::exp((STATES[V]+30.0000)/11.5000));
ALGEBRAIC[tau_xr1] =  1.00000*ALGEBRAIC[alpha_xr1]*ALGEBRAIC[beta_xr1];
RATES[Xr1] = (ALGEBRAIC[xr1_inf] - STATES[Xr1])/ALGEBRAIC[tau_xr1];
ALGEBRAIC[xr2_inf] = 1.00000/(1.00000+std::exp((STATES[V]+88.0000)/24.0000));
ALGEBRAIC[alpha_xr2] = 3.00000/(1.00000+std::exp((- 60.0000 - STATES[V])/20.0000));
ALGEBRAIC[beta_xr2] = 1.12000/(1.00000+std::exp((STATES[V] - 60.0000)/20.0000));
ALGEBRAIC[tau_xr2] =  1.00000*ALGEBRAIC[alpha_xr2]*ALGEBRAIC[beta_xr2];
RATES[Xr2] = (ALGEBRAIC[xr2_inf] - STATES[Xr2])/ALGEBRAIC[tau_xr2];

//IKs
ALGEBRAIC[xs_inf] = 1.00000/(1.00000+std::exp((- 5.00000 - STATES[V])/14.0000));
ALGEBRAIC[alpha_xs] = 1100.00/ std::pow((1.00000+std::exp((- 10.0000 - STATES[V])/6.00000)), 1.0 / 2);
ALGEBRAIC[beta_xs] = 1.00000/(1.00000+std::exp((STATES[V] - 60.0000)/20.0000));
ALGEBRAIC[tau_xs] =  1.00000*ALGEBRAIC[alpha_xs]*ALGEBRAIC[beta_xs];
RATES[Xs] = (ALGEBRAIC[xs_inf] - STATES[Xs])/ALGEBRAIC[tau_xs];

//SR calcium RyR gate g
ALGEBRAIC[g_inf] = (STATES[Ca_i]<0.000350000 ? 1.00000/(1.00000+std::pow(STATES[Ca_i]/0.000350000, 6.00000)) : 1.00000/(1.00000+std::pow(STATES[Ca_i]/0.000350000, 16.0000)));
ALGEBRAIC[d_g] = (ALGEBRAIC[g_inf] - STATES[g])/CONSTANTS[tau_g];
RATES[g] = (ALGEBRAIC[g_inf]>STATES[g]&&STATES[V]>- 60.0000 ? 0.00000 : ALGEBRAIC[d_g]);


ALGEBRAIC[i_NaK] = (( (( CONSTANTS[P_NaK]*CONSTANTS[K_o])/(CONSTANTS[K_o]+CONSTANTS[K_mk]))*STATES[Na_i])/(STATES[Na_i]+CONSTANTS[K_mNa]))/(1.00000+ 0.124500*std::exp(( - 0.100000*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))+ 0.0353000*std::exp(( - STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T])));
ALGEBRAIC[E_Na] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*std::log(CONSTANTS[Na_o]/STATES[Na_i]);
ALGEBRAIC[i_Na] =  CONSTANTS[g_Na]*std::pow(STATES[m], 3.00000)*STATES[h]*STATES[j]*(STATES[V] - ALGEBRAIC[E_Na]);
ALGEBRAIC[i_b_Na] =  CONSTANTS[g_bNa]*(STATES[V] - ALGEBRAIC[E_Na]);
ALGEBRAIC[i_NaCa] = ( CONSTANTS[K_NaCa]*( std::exp(( CONSTANTS[gamma_NCX]*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))*std::pow(STATES[Na_i], 3.00000)*CONSTANTS[Ca_o] -  std::exp(( (CONSTANTS[gamma_NCX] - 1.00000)*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))*std::pow(CONSTANTS[Na_o], 3.00000)*STATES[Ca_i]*CONSTANTS[alpha_NCX]))/( (std::pow(CONSTANTS[Km_Nai], 3.00000)+std::pow(CONSTANTS[Na_o], 3.00000))*(CONSTANTS[Km_Ca]+CONSTANTS[Ca_o])*(1.00000+ CONSTANTS[K_sat]*std::exp(( (CONSTANTS[gamma_NCX] - 1.00000)*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))));
RATES[Na_i] = ( - (ALGEBRAIC[i_Na]+ALGEBRAIC[i_b_Na]+ 3.00000*ALGEBRAIC[i_NaK]+ 3.00000*ALGEBRAIC[i_NaCa])*CONSTANTS[Cm])/( CONSTANTS[V_c]*CONSTANTS[F]);

ALGEBRAIC[E_K] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*std::log(CONSTANTS[K_o]/STATES[K_i]);
ALGEBRAIC[alpha_K1] = 0.100000/(1.00000+std::exp( 0.0600000*((STATES[V] - ALGEBRAIC[E_K]) - 200.000)));
ALGEBRAIC[beta_K1] = ( 3.00000*std::exp( 0.000200000*((STATES[V] - ALGEBRAIC[E_K])+100.000))+ 1.00000*std::exp( 0.100000*((STATES[V] - ALGEBRAIC[E_K]) - 10.0000)))/(1.00000+std::exp( - 0.500000*(STATES[V] - ALGEBRAIC[E_K])));
ALGEBRAIC[xK1_inf] = ALGEBRAIC[alpha_K1]/(ALGEBRAIC[alpha_K1]+ALGEBRAIC[beta_K1]);
ALGEBRAIC[i_K1] =  CONSTANTS[g_K1]*ALGEBRAIC[xK1_inf]* std::pow((CONSTANTS[K_o]/5.40000), 1.0 / 2)*(STATES[V] - ALGEBRAIC[E_K]);

ALGEBRAIC[i_to] =  CONSTANTS[g_to]*STATES[r]*STATES[s]*(STATES[V] - ALGEBRAIC[E_K]);

ALGEBRAIC[i_Kr] =  CONSTANTS[g_Kr]* std::pow((CONSTANTS[K_o]/5.40000), 1.0 / 2)*STATES[Xr1]*STATES[Xr2]*(STATES[V] - ALGEBRAIC[E_K]);
ALGEBRAIC[E_Ks] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*std::log((CONSTANTS[K_o]+ CONSTANTS[P_kna]*CONSTANTS[Na_o])/(STATES[K_i]+ CONSTANTS[P_kna]*STATES[Na_i]));
ALGEBRAIC[i_Ks] =  CONSTANTS[g_Ks]*std::pow(STATES[Xs], 2.00000)*(STATES[V] - ALGEBRAIC[E_Ks]);
ALGEBRAIC[i_CaL] = ( (( CONSTANTS[g_CaL]*STATES[d]*STATES[f]*STATES[fCa]*4.00000*STATES[V]*std::pow(CONSTANTS[F], 2.00000))/( CONSTANTS[R]*CONSTANTS[T]))*( STATES[Ca_i]*std::exp(( 2.00000*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T])) -  0.341000*CONSTANTS[Ca_o]))/(std::exp(( 2.00000*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T])) - 1.00000);
ALGEBRAIC[E_Ca] =  (( 0.500000*CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*std::log(CONSTANTS[Ca_o]/STATES[Ca_i]);
ALGEBRAIC[i_b_Ca] =  CONSTANTS[g_bCa]*(STATES[V] - ALGEBRAIC[E_Ca]);
ALGEBRAIC[i_p_K] = ( CONSTANTS[g_pK]*(STATES[V] - ALGEBRAIC[E_K]))/(1.00000+std::exp((25.0000 - STATES[V])/5.98000));
ALGEBRAIC[i_p_Ca] = ( CONSTANTS[g_pCa]*STATES[Ca_i])/(STATES[Ca_i]+CONSTANTS[K_pCa]);

ALGEBRAIC[Iion_cm] = ALGEBRAIC[i_K1]+ALGEBRAIC[i_to]+ALGEBRAIC[i_Kr]+
                ALGEBRAIC[i_Ks]+ALGEBRAIC[i_CaL]+ALGEBRAIC[i_NaK]+
                ALGEBRAIC[i_Na]+ALGEBRAIC[i_b_Na]+ALGEBRAIC[i_NaCa]+
                ALGEBRAIC[i_b_Ca]+ALGEBRAIC[i_p_K]+ALGEBRAIC[i_p_Ca];

ALGEBRAIC[Istim] = computeIstim(VOI, CONSTANTS);

RATES[V] = 0.0;
if (solveVmWithinODESolver)
    {
        RATES[V] = -ALGEBRAIC[Iion_cm] - ALGEBRAIC[Istim];
    }

RATES[K_i] = ( - ((ALGEBRAIC[i_K1]+ALGEBRAIC[i_to]+ALGEBRAIC[i_Kr]+ALGEBRAIC[i_Ks]+ALGEBRAIC[i_p_K]+ALGEBRAIC[Istim]) -  2.00000*ALGEBRAIC[i_NaK])*CONSTANTS[Cm])/( CONSTANTS[V_c]*CONSTANTS[F]);
ALGEBRAIC[i_rel] =  (( CONSTANTS[a_rel]*std::pow(STATES[Ca_SR], 2.00000))/(std::pow(CONSTANTS[b_rel], 2.00000)+std::pow(STATES[Ca_SR], 2.00000))+CONSTANTS[c_rel])*STATES[d]*STATES[g];
ALGEBRAIC[i_up] = CONSTANTS[Vmax_up]/(1.00000+std::pow(CONSTANTS[K_up], 2.00000)/std::pow(STATES[Ca_i], 2.00000));
ALGEBRAIC[i_leak] =  CONSTANTS[V_leak]*(STATES[Ca_SR] - STATES[Ca_i]);
ALGEBRAIC[ddt_Ca_i_total] = (( (- ((ALGEBRAIC[i_CaL]+ALGEBRAIC[i_b_Ca]+ALGEBRAIC[i_p_Ca]) -  2.00000*ALGEBRAIC[i_NaCa])/( 2.00000*CONSTANTS[V_c]*CONSTANTS[F]))*CONSTANTS[Cm]+ALGEBRAIC[i_leak]) - ALGEBRAIC[i_up])+ALGEBRAIC[i_rel];
ALGEBRAIC[f_JCa_i_free] = 1.00000/(1.00000+( CONSTANTS[Buf_c]*CONSTANTS[K_buf_c])/std::pow(STATES[Ca_i]+CONSTANTS[K_buf_c], 2.00000));
RATES[Ca_i] =  ALGEBRAIC[ddt_Ca_i_total]*ALGEBRAIC[f_JCa_i_free];
ALGEBRAIC[ddt_Ca_sr_total] =  (CONSTANTS[V_c]/CONSTANTS[V_sr])*(ALGEBRAIC[i_up] - (ALGEBRAIC[i_rel]+ALGEBRAIC[i_leak]));
ALGEBRAIC[f_JCa_sr_free] = 1.00000/(1.00000+( CONSTANTS[Buf_sr]*CONSTANTS[K_buf_sr])/std::pow(STATES[Ca_SR]+CONSTANTS[K_buf_sr], 2.00000));
RATES[Ca_SR] =  ALGEBRAIC[ddt_Ca_sr_total]*ALGEBRAIC[f_JCa_sr_free];
}





void
TNNPcomputeVariables(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, int tissueFlag, bool solveVmWithinODESolver)
{

ALGEBRAIC[f_inf] = 1.00000/(1.00000+std::exp((STATES[V]+20.0000)/7.00000));
ALGEBRAIC[tau_f] =  1125.00*std::exp(- std::pow(STATES[V]+27.0000, 2.00000)/240.000)+80.0000+165.000/(1.00000+std::exp((25.0000 - STATES[V])/10.0000));

//Conditional tissue flag for the Ito current
ALGEBRAIC[s_inf] = (tissueFlag == 1)
    ? 1.00000 / (1.00000 + std::exp((STATES[V] + 20.0000) / 5.00000))
    : 1.10000 / (1.00000 + std::exp((STATES[V] + 28.0000) / 6.00000));

ALGEBRAIC[tau_s] = (tissueFlag == 1)
     ? 85.0000*std::exp(- std::pow(STATES[V]+45.0000, 2.00000)/320.000)+5.00000/(1.00000+std::exp((STATES[V] - 20.0000)/5.00000))+3.00000
     : 1000.0000*std::exp(- std::pow(STATES[V]+67.0000, 2.00000)/1000.000)+8.00000;

ALGEBRAIC[r_inf] = 1.00000/(1.00000+std::exp((20.0000 - STATES[V])/6.00000));
ALGEBRAIC[tau_r] =  9.50000*std::exp(- std::pow(STATES[V]+40.0000, 2.00000)/1800.00)+0.800000;
ALGEBRAIC[g_inf] = (STATES[Ca_i]<0.000350000 ? 1.00000/(1.00000+std::pow(STATES[Ca_i]/0.000350000, 6.00000)) : 1.00000/(1.00000+std::pow(STATES[Ca_i]/0.000350000, 16.0000)));
ALGEBRAIC[d_g] = (ALGEBRAIC[g_inf] - STATES[g])/CONSTANTS[tau_g];
ALGEBRAIC[xr1_inf] = 1.00000/(1.00000+std::exp((- 26.0000 - STATES[V])/7.00000));
ALGEBRAIC[alpha_xr1] = 450.000/(1.00000+std::exp((- 45.0000 - STATES[V])/10.0000));
ALGEBRAIC[beta_xr1] = 6.00000/(1.00000+std::exp((STATES[V]+30.0000)/11.5000));
ALGEBRAIC[tau_xr1] =  1.00000*ALGEBRAIC[alpha_xr1]*ALGEBRAIC[beta_xr1];
ALGEBRAIC[xr2_inf] = 1.00000/(1.00000+std::exp((STATES[V]+88.0000)/24.0000));
ALGEBRAIC[alpha_xr2] = 3.00000/(1.00000+std::exp((- 60.0000 - STATES[V])/20.0000));
ALGEBRAIC[beta_xr2] = 1.12000/(1.00000+std::exp((STATES[V] - 60.0000)/20.0000));
ALGEBRAIC[tau_xr2] =  1.00000*ALGEBRAIC[alpha_xr2]*ALGEBRAIC[beta_xr2];
ALGEBRAIC[xs_inf] = 1.00000/(1.00000+std::exp((- 5.00000 - STATES[V])/14.0000));
ALGEBRAIC[alpha_xs] = 1100.00/ std::pow((1.00000+std::exp((- 10.0000 - STATES[V])/6.00000)), 1.0 / 2);
ALGEBRAIC[beta_xs] = 1.00000/(1.00000+std::exp((STATES[V] - 60.0000)/20.0000));
ALGEBRAIC[tau_xs] =  1.00000*ALGEBRAIC[alpha_xs]*ALGEBRAIC[beta_xs];
ALGEBRAIC[m_inf] = 1.00000/std::pow(1.00000+std::exp((- 56.8600 - STATES[V])/9.03000), 2.00000);
ALGEBRAIC[alpha_m] = 1.00000/(1.00000+std::exp((- 60.0000 - STATES[V])/5.00000));
ALGEBRAIC[beta_m] = 0.100000/(1.00000+std::exp((STATES[V]+35.0000)/5.00000))+0.100000/(1.00000+std::exp((STATES[V] - 50.0000)/200.000));
ALGEBRAIC[tau_m] =  1.00000*ALGEBRAIC[alpha_m]*ALGEBRAIC[beta_m];
ALGEBRAIC[h_inf] = 1.00000/std::pow(1.00000+std::exp((STATES[V]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[alpha_h] = (STATES[V]<- 40.0000 ?  0.0570000*std::exp(- (STATES[V]+80.0000)/6.80000) : 0.00000);
ALGEBRAIC[beta_h] = (STATES[V]<- 40.0000 ?  2.70000*std::exp( 0.0790000*STATES[V])+ 310000.*std::exp( 0.348500*STATES[V]) : 0.770000/( 0.130000*(1.00000+std::exp((STATES[V]+10.6600)/- 11.1000))));
ALGEBRAIC[tau_h] = 1.00000/(ALGEBRAIC[alpha_h]+ALGEBRAIC[beta_h]);
ALGEBRAIC[j_inf] = 1.00000/std::pow(1.00000+std::exp((STATES[V]+71.5500)/7.43000), 2.00000);
ALGEBRAIC[alpha_j] = (STATES[V]<- 40.0000 ? (( ( - 25428.0*std::exp( 0.244400*STATES[V]) -  6.94800e-06*std::exp( - 0.0439100*STATES[V]))*(STATES[V]+37.7800))/1.00000)/(1.00000+std::exp( 0.311000*(STATES[V]+79.2300))) : 0.00000);
ALGEBRAIC[beta_j] = (STATES[V]<- 40.0000 ? ( 0.0242400*std::exp( - 0.0105200*STATES[V]))/(1.00000+std::exp( - 0.137800*(STATES[V]+40.1400))) : ( 0.600000*std::exp( 0.0570000*STATES[V]))/(1.00000+std::exp( - 0.100000*(STATES[V]+32.0000))));
ALGEBRAIC[tau_j] = 1.00000/(ALGEBRAIC[alpha_j]+ALGEBRAIC[beta_j]);
ALGEBRAIC[d_inf] = 1.00000/(1.00000+std::exp((- 5.00000 - STATES[V])/7.50000));
ALGEBRAIC[alpha_d] = 1.40000/(1.00000+std::exp((- 35.0000 - STATES[V])/13.0000))+0.250000;
ALGEBRAIC[gama_fCa] = 1.40000/(1.00000+std::exp((STATES[V]+5.00000)/5.00000));
ALGEBRAIC[gamma_d] = 1.00000/(1.00000+std::exp((50.0000 - STATES[V])/20.0000));
ALGEBRAIC[tau_d] =  1.00000*ALGEBRAIC[alpha_d]*ALGEBRAIC[gama_fCa]+ALGEBRAIC[gamma_d];
ALGEBRAIC[alpha_fCa] = 1.00000/(1.00000+std::pow(STATES[Ca_i]/0.000325000, 8.00000));
ALGEBRAIC[beta_fCa] = 0.100000/(1.00000+std::exp((STATES[Ca_i] - 0.000500000)/0.000100000));
ALGEBRAIC[beta_d] = 0.200000/(1.00000+std::exp((STATES[Ca_i] - 0.000750000)/0.000800000));
ALGEBRAIC[fCa_inf] = (ALGEBRAIC[alpha_fCa]+ALGEBRAIC[beta_fCa]+ALGEBRAIC[beta_d]+0.230000)/1.46000;
ALGEBRAIC[d_fCa] = (ALGEBRAIC[fCa_inf] - STATES[fCa])/CONSTANTS[tau_fCa];
ALGEBRAIC[i_NaK] = (( (( CONSTANTS[P_NaK]*CONSTANTS[K_o])/(CONSTANTS[K_o]+CONSTANTS[K_mk]))*STATES[Na_i])/(STATES[Na_i]+CONSTANTS[K_mNa]))/(1.00000+ 0.124500*std::exp(( - 0.100000*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))+ 0.0353000*std::exp(( - STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T])));
ALGEBRAIC[E_Na] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*std::log(CONSTANTS[Na_o]/STATES[Na_i]);
ALGEBRAIC[i_Na] =  CONSTANTS[g_Na]*std::pow(STATES[m], 3.00000)*STATES[h]*STATES[j]*(STATES[V] - ALGEBRAIC[E_Na]);
ALGEBRAIC[i_b_Na] =  CONSTANTS[g_bNa]*(STATES[V] - ALGEBRAIC[E_Na]);
ALGEBRAIC[i_NaCa] = ( CONSTANTS[K_NaCa]*( std::exp(( CONSTANTS[gamma_NCX]*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))*std::pow(STATES[Na_i], 3.00000)*CONSTANTS[Ca_o] -  std::exp(( (CONSTANTS[gamma_NCX] - 1.00000)*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))*std::pow(CONSTANTS[Na_o], 3.00000)*STATES[Ca_i]*CONSTANTS[alpha_NCX]))/( (std::pow(CONSTANTS[Km_Nai], 3.00000)+std::pow(CONSTANTS[Na_o], 3.00000))*(CONSTANTS[Km_Ca]+CONSTANTS[Ca_o])*(1.00000+ CONSTANTS[K_sat]*std::exp(( (CONSTANTS[gamma_NCX] - 1.00000)*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T]))));
ALGEBRAIC[E_K] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*std::log(CONSTANTS[K_o]/STATES[K_i]);
ALGEBRAIC[alpha_K1] = 0.100000/(1.00000+std::exp( 0.0600000*((STATES[V] - ALGEBRAIC[E_K]) - 200.000)));
ALGEBRAIC[beta_K1] = ( 3.00000*std::exp( 0.000200000*((STATES[V] - ALGEBRAIC[E_K])+100.000))+ 1.00000*std::exp( 0.100000*((STATES[V] - ALGEBRAIC[E_K]) - 10.0000)))/(1.00000+std::exp( - 0.500000*(STATES[V] - ALGEBRAIC[E_K])));
ALGEBRAIC[xK1_inf] = ALGEBRAIC[alpha_K1]/(ALGEBRAIC[alpha_K1]+ALGEBRAIC[beta_K1]);
ALGEBRAIC[i_K1] =  CONSTANTS[g_K1]*ALGEBRAIC[xK1_inf]* std::pow((CONSTANTS[K_o]/5.40000), 1.0 / 2)*(STATES[V] - ALGEBRAIC[E_K]);
ALGEBRAIC[i_to] =  CONSTANTS[g_to]*STATES[r]*STATES[s]*(STATES[V] - ALGEBRAIC[E_K]);
ALGEBRAIC[i_Kr] =  CONSTANTS[g_Kr]* std::pow((CONSTANTS[K_o]/5.40000), 1.0 / 2)*STATES[Xr1]*STATES[Xr2]*(STATES[V] - ALGEBRAIC[E_K]);
ALGEBRAIC[E_Ks] =  (( CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*std::log((CONSTANTS[K_o]+ CONSTANTS[P_kna]*CONSTANTS[Na_o])/(STATES[K_i]+ CONSTANTS[P_kna]*STATES[Na_i]));
ALGEBRAIC[i_Ks] =  CONSTANTS[g_Ks]*std::pow(STATES[Xs], 2.00000)*(STATES[V] - ALGEBRAIC[E_Ks]);
ALGEBRAIC[i_CaL] = ( (( CONSTANTS[g_CaL]*STATES[d]*STATES[f]*STATES[fCa]*4.00000*STATES[V]*std::pow(CONSTANTS[F], 2.00000))/( CONSTANTS[R]*CONSTANTS[T]))*( STATES[Ca_i]*std::exp(( 2.00000*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T])) -  0.341000*CONSTANTS[Ca_o]))/(std::exp(( 2.00000*STATES[V]*CONSTANTS[F])/( CONSTANTS[R]*CONSTANTS[T])) - 1.00000);
ALGEBRAIC[E_Ca] =  (( 0.500000*CONSTANTS[R]*CONSTANTS[T])/CONSTANTS[F])*std::log(CONSTANTS[Ca_o]/STATES[Ca_i]);
ALGEBRAIC[i_b_Ca] =  CONSTANTS[g_bCa]*(STATES[V] - ALGEBRAIC[E_Ca]);
ALGEBRAIC[i_p_K] = ( CONSTANTS[g_pK]*(STATES[V] - ALGEBRAIC[E_K]))/(1.00000+std::exp((25.0000 - STATES[V])/5.98000));
ALGEBRAIC[i_p_Ca] = ( CONSTANTS[g_pCa]*STATES[Ca_i])/(STATES[Ca_i]+CONSTANTS[K_pCa]);

ALGEBRAIC[i_rel] =  (( CONSTANTS[a_rel]*std::pow(STATES[Ca_SR], 2.00000))/(std::pow(CONSTANTS[b_rel], 2.00000)+std::pow(STATES[Ca_SR], 2.00000))+CONSTANTS[c_rel])*STATES[d]*STATES[g];
ALGEBRAIC[i_up] = CONSTANTS[Vmax_up]/(1.00000+std::pow(CONSTANTS[K_up], 2.00000)/std::pow(STATES[Ca_i], 2.00000));
ALGEBRAIC[i_leak] =  CONSTANTS[V_leak]*(STATES[Ca_SR] - STATES[Ca_i]);
ALGEBRAIC[ddt_Ca_i_total] = (( (- ((ALGEBRAIC[i_CaL]+ALGEBRAIC[i_b_Ca]+ALGEBRAIC[i_p_Ca]) -  2.00000*ALGEBRAIC[i_NaCa])/( 2.00000*CONSTANTS[V_c]*CONSTANTS[F]))*CONSTANTS[Cm]+ALGEBRAIC[i_leak]) - ALGEBRAIC[i_up])+ALGEBRAIC[i_rel];
ALGEBRAIC[f_JCa_i_free] = 1.00000/(1.00000+( CONSTANTS[Buf_c]*CONSTANTS[K_buf_c])/std::pow(STATES[Ca_i]+CONSTANTS[K_buf_c], 2.00000));
ALGEBRAIC[ddt_Ca_sr_total] =  (CONSTANTS[V_c]/CONSTANTS[V_sr])*(ALGEBRAIC[i_up] - (ALGEBRAIC[i_rel]+ALGEBRAIC[i_leak]));
ALGEBRAIC[f_JCa_sr_free] = 1.00000/(1.00000+( CONSTANTS[Buf_sr]*CONSTANTS[K_buf_sr])/std::pow(STATES[Ca_SR]+CONSTANTS[K_buf_sr], 2.00000));

ALGEBRAIC[Iion_cm] = ALGEBRAIC[i_K1]+ALGEBRAIC[i_to]+ALGEBRAIC[i_Kr]+
                ALGEBRAIC[i_Ks]+ALGEBRAIC[i_CaL]+ALGEBRAIC[i_NaK]+
                ALGEBRAIC[i_Na]+ALGEBRAIC[i_b_Na]+ALGEBRAIC[i_NaCa]+
                ALGEBRAIC[i_b_Ca]+ALGEBRAIC[i_p_K]+ALGEBRAIC[i_p_Ca];

ALGEBRAIC[Istim] = computeIstim(VOI, CONSTANTS);
}



    // ten Tusscher-Noble-Noble-Panfilov (2004) ionic model, embedded in this
    // translation unit rather than taken from the ionicModels library.
    //
    // It owns one state vector per INTEGRATION POINT, not per cell: the
    // high-order path evaluates Iion at the cell quadrature points, so the ODE
    // is integrated there. The scratch arrays cover two disjoint banks - Gauss
    // points and cell centres - because the front-aware hybrid needs both alive
    // at once; see solveODE() and the note on cellSlotOffset_ above.
    class EmbeddedTNNPModel
    {
        dictionary dict_;
        label nPoints_;

        // Added for cardiacFoam: number of scratch slots reserved for
        // cell-centred integration, and the index at which that second bank
        // starts. The per-point scratch below is a single array covering two
        // disjoint point sets: Gauss-point slots occupy [0, nPoints_) and
        // cell-centre slots occupy [cellSlotOffset_, cellSlotOffset_+nCells_).
        // See solveODE() for why the two banks must not overlap.
        label nCellSlots_;
        label cellSlotOffset_;

        scalarField constants_;
        Field<Field<scalar>> algebraic_;
        Field<Field<scalar>> rates_;
        Field<Field<scalar>> oldStates_;
        wordList exportedNames_;
        wordList debugNames_;
        label tissue_;
        word odeSolverName_;
        scalar odeAbsTol_;
        scalar odeRelTol_;
        label odeMaxSteps_;
        scalar odeInitialStepMs_;
        Switch numericalProtection_;
        Switch clampGates_;
        Switch clampVmForRates_;
        Switch logNumericalProtection_;
        scalar concentrationFloor_;
        scalar vmMinForRates_;
        scalar vmMaxForRates_;
        scalar vmZeroEps_;
        label maxProtectionLogEntries_;
        word protectionLogName_;
        label protectionCorrections_;

        // Added for cardiacFoam: WHAT the protection corrects, by state and by
        // reason, with the first occurrence of each pair. The total alone
        // (capped at maxProtectionLogEntries_) said "1000000" and nothing else.
        // With TNNPProtectionDiagnostics on and the protection OFF, the same
        // tests run and are recorded but no state is touched: it shows what
        // leaves the physical range in an unprotected run, before it ends in a
        // SIGFPE. The first occurrence of each pair is printed at once, because
        // a SIGFPE never reaches the destructor that writes the summary.
        enum ProtectionReason
        {
            prNonFinite, prVmRange, prVmZero, prConcentration, prGate,
            nProtectionReasons
        };
        Switch protectionDiagnostics_;
        List<label> protectionCount_;
        List<scalar> protectionFirstTimeMs_;
        List<label> protectionFirstPoint_;
        List<scalar> protectionFirstValue_;
        List<scalar> protectionMinValue_;
        List<scalar> protectionMaxValue_;
        List<word> protectionFirstPhase_;

        static const char* protectionReasonName(const label r)
        {
            static const char* names[] =
                {"nonFinite", "VmRange", "VmZero", "concentrationFloor", "gateRange"};
            return names[r];
        }

        // Added for cardiacFoam: OpenMP over the ionic ODEs. Every point is an
        // independent initial-value problem sharing no state, so the loop is
        // embarrassingly parallel - this is the same treatment the MMS solver
        // already applies in updateIntegrationPointStatesODE.
        //
        // Off by default: it interacts with the MPI rank count, and oversubscribing
        // a node is slower than either alone.
        Switch parallelODE_;
        label parallelODEThreshold_;

        static label tissueFlag(const word& tissueName)
        {
            if (tissueName == "epicardialCells") return 1;
            if (tissueName == "mCells") return 2;
            if (tissueName == "endocardialCells") return 3;
            if (tissueName == "myocyte") return 4;

            FatalErrorInFunction
                << "Unsupported TNNP tissue '" << tissueName << "'. Valid options are "
                << "epicardialCells, mCells, endocardialCells, myocyte."
                << exit(FatalError);

            return 1;
        }

        static label stateIndex(const word& name)
        {
            for (label i = 0; i < NUM_STATES; ++i)
            {
                if (name == TNNP_STATES_NAMES[i]) return i;
            }
            return -1;
        }

        static label algebraicIndex(const word& name)
        {
            for (label i = 0; i < NUM_ALGEBRAIC; ++i)
            {
                if (name == TNNP_ALGEBRAIC_NAMES[i]) return i;
            }
            return -1;
        }

        void loadStimulusProtocol()
        {
            constants_[stim_start] = dict_.lookupOrDefault<scalar>("stim_start", 0.0);
            constants_[stim_period_S1] = dict_.lookupOrDefault<scalar>("stim_period_S1", 0.0);
            constants_[stim_duration] = dict_.lookupOrDefault<scalar>("stim_duration", 0.0);
            constants_[stim_amplitude] = dict_.lookupOrDefault<scalar>("stim_amplitude", 0.0);
            constants_[nstim1] = dict_.lookupOrDefault<scalar>("nstim1", 0.0);
            constants_[stim_period_S2] = dict_.lookupOrDefault<scalar>("stim_period_S2", 0.0);
            constants_[nstim2] = dict_.lookupOrDefault<scalar>("nstim2", 0.0);
        }

        // Modified for cardiacFoam: the only shared mutable state on the ODE
        // path. The bound test and the increment must be one indivisible step,
        // so this is "critical" rather than "atomic" - an atomic increment would
        // still race with the test. It only executes when a correction actually
        // fires, which is rare, so it does not serialise the loop.
        void recordCorrection
        (
            const label stateI,
            const label reason,
            const scalar tMs,
            const label pointI,
            const char* phase,
            const scalar value,
            const scalar VmMilliVolt
        )
        {
            #ifdef _OPENMP
            #pragma omp critical(cardiacFoamProtectionCount)
            #endif
            {
                if (protectionCorrections_ < labelMax)
                {
                    ++protectionCorrections_;
                }
                const label k = stateI*nProtectionReasons + reason;
                if (protectionCount_[k] == 0)
                {
                    protectionFirstTimeMs_[k] = tMs;
                    protectionFirstPoint_[k] = pointI;
                    protectionFirstValue_[k] = value;
                    protectionFirstPhase_[k] = phase;
                    Pout<< "TNNP protection "
                        << (numericalProtection_ ? "[corrected]" : "[diagnostic only]")
                        << ": first " << protectionReasonName(reason)
                        << " of " << TNNP_STATES_NAMES[stateI]
                        << " at t = " << tMs << " ms, point " << pointI
                        << ", phase " << phase << ", value " << value
                        << " (Vm = " << VmMilliVolt << " mV)" << endl;
                }
                if (protectionCount_[k] < labelMax)
                {
                    ++protectionCount_[k];
                }
                if (std::isfinite(value))
                {
                    protectionMinValue_[k] = min(protectionMinValue_[k], value);
                    protectionMaxValue_[k] = max(protectionMaxValue_[k], value);
                }
            }
        }

        void protectState
        (
            scalarField& state,
            const scalar tMs,
            const label pointI,
            const char* phase
        )
        {
            if (!numericalProtection_ && !protectionDiagnostics_) return;

            // With the protection off (diagnostics only) the test is recorded
            // and the state is left alone, so the checks below see the raw
            // values.
            auto correct = [&](const label i, const scalar value, const label reason)
            {
                if (state[i] != value)
                {
                    recordCorrection(i, reason, tMs, pointI, phase, state[i], state[V]);
                    if (numericalProtection_)
                    {
                        state[i] = value;
                    }
                }
            };

            if (!std::isfinite(state[V]))
            {
                correct(V, -86.2, prNonFinite);
            }
            if (clampVmForRates_)
            {
                correct(V, min(max(state[V], vmMinForRates_), vmMaxForRates_), prVmRange);
            }
            if (mag(state[V]) < vmZeroEps_)
            {
                correct(V, state[V] < 0.0 ? -vmZeroEps_ : vmZeroEps_, prVmZero);
            }

            const label concentrationIDs[] = {K_i, Na_i, Ca_i, Ca_SR};
            for (const label id : concentrationIDs)
            {
                if (!std::isfinite(state[id]))
                {
                    correct(id, concentrationFloor_, prNonFinite);
                }
                else if (state[id] < concentrationFloor_)
                {
                    correct(id, concentrationFloor_, prConcentration);
                }
            }

            if (clampGates_)
            {
                const label gateIDs[] = {Xr1, Xr2, Xs, m, h, j, d, f, fCa, s, r, g};
                for (const label id : gateIDs)
                {
                    if (!std::isfinite(state[id]))
                    {
                        correct(id, 0.0, prNonFinite);
                    }
                    else
                    {
                        correct
                        (
                            id, min(max(state[id], scalar(0.0)), scalar(1.0)), prGate
                        );
                    }
                }
            }
        }

        // Added for cardiacFoam: the membrane potential [mV] the driver
        // prescribes at time tMs of an integration over [t0Ms, t0Ms + dtMs].
        static scalar drivenVmMilliVolt
        (
            const VmDriverNodes& nodes,
            const scalar tMs,
            const scalar t0Ms,
            const scalar dtMs
        )
        {
            const scalar s = (dtMs > SMALL) ? (tMs - t0Ms)/dtMs : 1.0;
            return 1000.0*interpolateVmDriver(nodes, s);
        }

        // Modified for cardiacFoam: with a driver, the rates are evaluated with
        // V replaced by the driven Vm(tMs) instead of whatever the state vector
        // carries. V is not integrated here (its rate is zero off the
        // single-cell path), so the state entry only ever held the value set at
        // the start of the step; the driver is what makes Vm vary across it.
        // Without a driver nothing changes, bit for bit.
        void computeRates
        (
            const scalar tMs,
            const scalarField& y,
            scalarField& dydt,
            scalarField& algebraic,
            const label pointI,
            const char* phase,
            const VmDriverNodes* nodes = nullptr,
            const scalar t0Ms = 0.0,
            const scalar dtMs = 0.0
        )
        {
            if (numericalProtection_)
            {
                scalarField yProtected(y);
                if (nodes)
                {
                    yProtected[V] = drivenVmMilliVolt(*nodes, tMs, t0Ms, dtMs);
                }
                protectState(yProtected, tMs, pointI, phase);
                TNNPcomputeRates
                (
                    tMs,
                    constants_.data(),
                    dydt.data(),
                    yProtected.data(),
                    algebraic.data(),
                    tissue_,
                    false
                );
            }
            else
            {
                scalarField yCopy(y);
                if (nodes)
                {
                    yCopy[V] = drivenVmMilliVolt(*nodes, tMs, t0Ms, dtMs);
                }
                if (protectionDiagnostics_)
                {
                    // Records only: the protection is off, nothing is changed.
                    protectState(yCopy, tMs, pointI, phase);
                }
                TNNPcomputeRates
                (
                    tMs,
                    constants_.data(),
                    dydt.data(),
                    yCopy.data(),
                    algebraic.data(),
                    tissue_,
                    false
                );
            }
        }

        void computeVariables
        (
            const scalar tMs,
            scalarField& y,
            scalarField& rates,
            scalarField& algebraic
        )
        {
            TNNPcomputeVariables
            (
                tMs,
                constants_.data(),
                rates.data(),
                y.data(),
                algebraic.data(),
                tissue_,
                false
            );
        }

        void eulerStep
        (
            const scalar tMs,
            const scalar hMs,
            scalarField& y,
            scalarField& rates,
            scalarField& algebraic,
            const label pointI,
            const VmDriverNodes* nodes,
            const scalar t0Ms,
            const scalar dtMs
        )
        {
            computeRates(tMs, y, rates, algebraic, pointI, "Euler", nodes, t0Ms, dtMs);
            forAll(y, i)
            {
                y[i] += hMs*rates[i];
            }
        }

        void rk4Step
        (
            const scalar tMs,
            const scalar hMs,
            scalarField& y,
            scalarField& rates,
            scalarField& algebraic,
            const label pointI,
            const VmDriverNodes* nodes,
            const scalar t0Ms,
            const scalar dtMs
        )
        {
            scalarField k1(NUM_STATES, 0.0), k2(NUM_STATES, 0.0);
            scalarField k3(NUM_STATES, 0.0), k4(NUM_STATES, 0.0);
            scalarField yt(NUM_STATES, 0.0);

            computeRates(tMs, y, k1, algebraic, pointI, "RK4_k1", nodes, t0Ms, dtMs);
            forAll(y, i) yt[i] = y[i] + 0.5*hMs*k1[i];
            computeRates(tMs + 0.5*hMs, yt, k2, algebraic, pointI, "RK4_k2", nodes, t0Ms, dtMs);
            forAll(y, i) yt[i] = y[i] + 0.5*hMs*k2[i];
            computeRates(tMs + 0.5*hMs, yt, k3, algebraic, pointI, "RK4_k3", nodes, t0Ms, dtMs);
            forAll(y, i) yt[i] = y[i] + hMs*k3[i];
            computeRates(tMs + hMs, yt, k4, algebraic, pointI, "RK4_k4", nodes, t0Ms, dtMs);

            forAll(y, i)
            {
                y[i] += (hMs/6.0)*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
            }
            rates = k4;
        }

        void rkf45Trial
        (
            const scalar tMs,
            const scalar hMs,
            const scalarField& y,
            scalarField& yFifth,
            scalarField& rates,
            scalarField& algebraic,
            const label pointI,
            scalar& err,
            const VmDriverNodes* nodes,
            const scalar t0Ms,
            const scalar dtMs
        )
        {
            scalarField k1(NUM_STATES, 0.0), k2(NUM_STATES, 0.0), k3(NUM_STATES, 0.0);
            scalarField k4(NUM_STATES, 0.0), k5(NUM_STATES, 0.0), k6(NUM_STATES, 0.0);
            scalarField yt(NUM_STATES, 0.0), yFourth(NUM_STATES, 0.0);

            // Added for cardiacFoam: an intermediate stage, or the result, that
            // leaves the domain where the rates are defined is a REJECTED trial
            // (err = GREAT, so the caller shrinks h), not a state to evaluate.
            // The rates take the log of K_i, Na_i and Ca_i (the Nernst
            // potentials): a non-positive one is FE_INVALID, which OpenFOAM
            // traps as a SIGFPE BEFORE the error estimate could reject the
            // trial. Found on a GP run, 2-D hexa 0.5 mm, without the numerical
            // protection: stage k5 at one Gauss point gave Na_i = -401.6 mM
            // with Vm = -102.7 mV at t = 12.31 ms, and the run died there at
            // dt = 1e-4 and at dt = 5e-5 alike. No state is modified: this is
            // the step-size control doing its job, which the protection used to
            // hide by clamping the trial state instead.
            auto inDomain = [](const scalarField& s)
            {
                forAll(s, i)
                {
                    if (!std::isfinite(s[i])) return false;
                }
                return s[K_i] > 0.0 && s[Na_i] > 0.0 && s[Ca_i] > 0.0;
            };
            auto reject = [&]()
            {
                yFifth = y;
                err = GREAT;
            };

            computeRates(tMs, y, k1, algebraic, pointI, "RKF45_k1", nodes, t0Ms, dtMs);
            forAll(y, i) yt[i] = y[i] + hMs*(1.0/5.0)*k1[i];
            if (!inDomain(yt)) return reject();
            computeRates(tMs + hMs/5.0, yt, k2, algebraic, pointI, "RKF45_k2", nodes, t0Ms, dtMs);
            forAll(y, i) yt[i] = y[i] + hMs*((3.0/40.0)*k1[i] + (9.0/40.0)*k2[i]);
            if (!inDomain(yt)) return reject();
            computeRates(tMs + 3.0*hMs/10.0, yt, k3, algebraic, pointI, "RKF45_k3", nodes, t0Ms, dtMs);
            forAll(y, i) yt[i] = y[i] + hMs*((3.0/10.0)*k1[i] - (9.0/10.0)*k2[i] + (6.0/5.0)*k3[i]);
            if (!inDomain(yt)) return reject();
            computeRates(tMs + 3.0*hMs/5.0, yt, k4, algebraic, pointI, "RKF45_k4", nodes, t0Ms, dtMs);
            forAll(y, i) yt[i] = y[i] + hMs*((-11.0/54.0)*k1[i] + (5.0/2.0)*k2[i] - (70.0/27.0)*k3[i] + (35.0/27.0)*k4[i]);
            if (!inDomain(yt)) return reject();
            computeRates(tMs + hMs, yt, k5, algebraic, pointI, "RKF45_k5", nodes, t0Ms, dtMs);
            forAll(y, i) yt[i] = y[i] + hMs*((1631.0/55296.0)*k1[i] + (175.0/512.0)*k2[i] + (575.0/13824.0)*k3[i] + (44275.0/110592.0)*k4[i] + (253.0/4096.0)*k5[i]);
            if (!inDomain(yt)) return reject();
            computeRates(tMs + 7.0*hMs/8.0, yt, k6, algebraic, pointI, "RKF45_k6", nodes, t0Ms, dtMs);

            err = 0.0;
            forAll(y, i)
            {
                yFifth[i] = y[i] + hMs*((37.0/378.0)*k1[i] + (250.0/621.0)*k3[i] + (125.0/594.0)*k4[i] + (512.0/1771.0)*k6[i]);
                yFourth[i] = y[i] + hMs*((2825.0/27648.0)*k1[i] + (18575.0/48384.0)*k3[i] + (13525.0/55296.0)*k4[i] + (277.0/14336.0)*k5[i] + 0.25*k6[i]);
                const scalar scale = odeAbsTol_ + odeRelTol_*max(mag(yFifth[i]), mag(yFourth[i]));
                err = max(err, mag(yFifth[i] - yFourth[i])/max(scale, SMALL));
            }
            if (!inDomain(yFifth)) return reject();
            rates = k6;
        }

        // Modified for cardiacFoam: nodes, when given, is the Vm driver of this
        // point (see VmDriverNodes); its last node is VmVolt, the candidate, so
        // state[V] - set to VmVolt below and never integrated - still carries
        // the end-of-step value that the final computeVariables needs. Without
        // nodes the rates see VmVolt for the whole step, as before.
        void advancePoint
        (
            const scalar tStartMs,
            const scalar dtMs,
            const scalar VmVolt,
            const VmDriverNodes* nodes,
            scalarField& state,
            scalarField& algebraic,
            scalarField& rates,
            const label pointI
        )
        {
            if (state.size() != NUM_STATES)
            {
                state.setSize(NUM_STATES, 0.0);
                TNNPinitConsts(constants_.data(), rates.data(), state.data(), tissue_, dict_);
                loadStimulusProtocol();
            }

            state[V] = 1000.0*VmVolt;
            protectState(state, tStartMs, pointI, "solveODE_start");

            if (dtMs <= SMALL)
            {
                computeVariables(tStartMs, state, rates, algebraic);
                return;
            }

            if (odeSolverName_ == "Euler")
            {
                eulerStep(tStartMs, dtMs, state, rates, algebraic, pointI, nodes, tStartMs, dtMs);
                protectState(state, tStartMs + dtMs, pointI, "solveODE_end");
                computeVariables(tStartMs + dtMs, state, rates, algebraic);
                return;
            }

            if (odeSolverName_ == "RK4")
            {
                rk4Step(tStartMs, dtMs, state, rates, algebraic, pointI, nodes, tStartMs, dtMs);
                protectState(state, tStartMs + dtMs, pointI, "solveODE_end");
                computeVariables(tStartMs + dtMs, state, rates, algebraic);
                return;
            }

            // Modified for cardiacFoam: odeSolverName_ is validated in the
            // constructor (RKF45, RK4, Euler, one spelling each), so reaching
            // this point means RKF45. The aliases forwardEuler, RKT45 and rkf45
            // are gone.

            scalar tau = 0.0;
            // Modified for cardiacFoam, as in the MMS solver: every integration
            // starts from the same first sub-step, stateODEInitialStep.
            //
            // It used to start from the LAST accepted step of the previous
            // integration of the same scratch slot (stepMs_, seeded once with
            // dt). That made the answer depend on history the state itself does
            // not carry: on which slot a point maps to (so on the integration
            // mode and on the partition), on whether the step was rolled back
            // (the history was not restored), and on what was integrated
            // before - so two evaluations of the same Vm candidate differed at
            // the RKF45 tolerance, which is noise inside JFNK's finite-
            // difference Jacobian. Now the result is a function of its inputs.
            scalar h =
                odeInitialStepMs_ > SMALL
              ? min(dtMs, odeInitialStepMs_)
              : dtMs;
            const scalar hMin = max(dtMs*1.0e-12, SMALL);
            scalarField yTrial(NUM_STATES, 0.0);

            for (label subStep = 0; subStep < max(odeMaxSteps_, label(1)); ++subStep)
            {
                if (tau >= dtMs - SMALL) break;
                h = min(h, dtMs - tau);

                scalar err = GREAT;
                rkf45Trial
                (
                    tStartMs + tau,
                    h,
                    state,
                    yTrial,
                    rates,
                    algebraic,
                    pointI,
                    err,
                    nodes,
                    tStartMs,
                    dtMs
                );

                if (err <= 1.0 || h <= hMin)
                {
                    state = yTrial;
                    tau += h;
                }

                const scalar factor = min
                (
                    scalar(4.0),
                    max(scalar(0.1), scalar(0.84)*std::pow(1.0/(err + SMALL), 0.25))
                );
                h = max(hMin, min(dtMs - tau, h*factor));
            }

            if (tau < dtMs - 10.0*SMALL)
            {
                WarningInFunction
                    << "Embedded TNNP RKF45 reached maxSteps before completing dt. "
                    << "tau = " << tau << " ms, dt = " << dtMs << " ms" << nl;
            }

            protectState(state, tStartMs + dtMs, pointI, "solveODE_end");
            computeVariables(tStartMs + dtMs, state, rates, algebraic);
        }

    public:

        // Added for cardiacFoam: the "no slot map" default for solveODE().
        // A function-local static is used rather than labelUList::null() so
        // that .empty() is always well defined.
        static const labelList& emptyScratchIDs()
        {
            static const labelList empty;
            return empty;
        }

        // Added for cardiacFoam: threading of the per-point ODE loop, from the
        // shared OpenMP controls (openMPControls.H).
        void setOpenMP(const bool useOpenMP, const label threshold)
        {
            parallelODE_ = useOpenMP;
            parallelODEThreshold_ = threshold;
        }

        // Added for cardiacFoam: first scratch slot of the cell-centre bank.
        // Callers integrating cell-centred points alongside Gauss points build
        // their slot map as cellSlotOffset() + cellI.
        label cellSlotOffset() const
        {
            return cellSlotOffset_;
        }

        label nCellSlots() const
        {
            return nCellSlots_;
        }

        // Modified for cardiacFoam: takes nCellSlots (= mesh.nCells()) and sizes
        // the per-point scratch to nIntegrationPoints + nCellSlots, so that the
        // front-aware hybrid's Gauss-point pass and cell-centre pass get
        // disjoint scratch banks instead of both starting at slot 0.
        // oldStates_ is deliberately NOT extended: it is only ever assigned
        // whole-array (oldStates_ = states), so it resizes itself to the state
        // array in use and is not indexed by scratch slot.
        EmbeddedTNNPModel
        (
            const dictionary& dict,
            const label nIntegrationPoints,
            const label nCellSlots,
            const scalar initialDeltaT
        )
        :
            dict_(dict),
            nPoints_(nIntegrationPoints),
            nCellSlots_(nCellSlots),
            cellSlotOffset_(nIntegrationPoints),
            constants_(NUM_CONSTANTS, 0.0),
            algebraic_
            (
                nIntegrationPoints + nCellSlots,
                Field<scalar>(NUM_ALGEBRAIC, 0.0)
            ),
            rates_
            (
                nIntegrationPoints + nCellSlots,
                Field<scalar>(NUM_STATES, 0.0)
            ),
            oldStates_(nIntegrationPoints, Field<scalar>(NUM_STATES, 0.0)),
            exportedNames_(),
            debugNames_(),
            tissue_(tissueFlag(dict.lookupOrDefault<word>("tissue", "epicardialCells"))),
            odeSolverName_(dict.lookupOrDefault<word>("solver", dict.lookupOrDefault<word>("stateODESolver", "RKF45"))),
            odeAbsTol_(max(dict.lookupOrDefault<scalar>("absTol", dict.lookupOrDefault<scalar>("stateODEAbsTol", 1.0e-9)), scalar(SMALL))),
            odeRelTol_(max(dict.lookupOrDefault<scalar>("relTol", dict.lookupOrDefault<scalar>("stateODERelTol", 1.0e-6)), scalar(SMALL))),
            odeMaxSteps_(dict.lookupOrDefault<label>("maxSteps", dict.lookupOrDefault<label>("stateODEMaxSteps", 10000))),
            // Added for cardiacFoam, as in the MMS solver: the first sub-step
            // of EVERY RKF45 integration, in seconds like the MMS key (1e-5 s
            // by default), stored in ms because this integrator works in ms.
            odeInitialStepMs_
            (
                1000.0*dict.lookupOrDefault<scalar>("stateODEInitialStep", 1.0e-5)
            ),
            numericalProtection_(dict.lookupOrDefault<Switch>("TNNPNumericalProtection", false)),
            clampGates_(dict.lookupOrDefault<Switch>("TNNPClampGates", true)),
            clampVmForRates_(dict.lookupOrDefault<Switch>("TNNPClampVmForRates", true)),
            logNumericalProtection_(dict.lookupOrDefault<Switch>("TNNPLogNumericalProtection", true)),
            concentrationFloor_(max(dict.lookupOrDefault<scalar>("TNNPConcentrationFloor", 1.0e-12), scalar(SMALL))),
            vmMinForRates_(dict.lookupOrDefault<scalar>("TNNPVMinForRates", -200.0)),
            vmMaxForRates_(dict.lookupOrDefault<scalar>("TNNPVMaxForRates", 100.0)),
            vmZeroEps_(max(dict.lookupOrDefault<scalar>("TNNPVoltageZeroEps", 1.0e-8), scalar(SMALL))),
            maxProtectionLogEntries_(dict.lookupOrDefault<label>("TNNPMaxProtectionLogEntries", 1000000)),
            protectionLogName_(dict.lookupOrDefault<word>("TNNPProtectionLog", "TNNP_numericalProtection_summary.dat")),
            protectionCorrections_(0),
            protectionDiagnostics_
            (
                dict.lookupOrDefault<Switch>("TNNPProtectionDiagnostics", false)
            ),
            protectionCount_(NUM_STATES*nProtectionReasons, label(0)),
            protectionFirstTimeMs_(NUM_STATES*nProtectionReasons, scalar(0)),
            protectionFirstPoint_(NUM_STATES*nProtectionReasons, label(-1)),
            protectionFirstValue_(NUM_STATES*nProtectionReasons, scalar(0)),
            protectionMinValue_(NUM_STATES*nProtectionReasons, GREAT),
            protectionMaxValue_(NUM_STATES*nProtectionReasons, -GREAT),
            protectionFirstPhase_(NUM_STATES*nProtectionReasons),
            // Added for cardiacFoam. The key name matches the one the solver
            // report documents and the one the tutorial driver already writes
            // into this dictionary - which until now nothing read.
            // Set by setOpenMP() from the shared OpenMP controls
            // (stateODEUseOpenMP, stateODEOpenMPThreshold), as in the MMS
            // solver; they used to be TNNP keys of their own.
            parallelODE_(false),
            parallelODEThreshold_(256)
        {
            // Added for cardiacFoam: validated here, at start-up, instead of
            // at the first ODE integration.
            requireOneOf("stateODESolver", odeSolverName_, {"RKF45", "RK4", "Euler"});

            if (dict_.found("exportedVariables"))
            {
                dict_.lookup("exportedVariables") >> exportedNames_;
            }
            if (dict_.found("debugPrintVariables"))
            {
                dict_.lookup("debugPrintVariables") >> debugNames_;
            }

            Info<< "Using embedded TNNP ionic model" << nl
                << "    tissue flag = " << tissue_ << nl
                << "    ODE solver = " << odeSolverName_ << nl
                << "    absTol = " << odeAbsTol_ << nl
                << "    relTol = " << odeRelTol_ << nl
                << "    maxSteps = " << odeMaxSteps_ << nl;

            if (numericalProtection_)
            {
                Info<< "    numerical protection enabled" << nl;
            }
        }

        ~EmbeddedTNNPModel()
        {
            if ((numericalProtection_ && logNumericalProtection_) || protectionDiagnostics_)
            {
                // Modified for cardiacFoam: protectionLogName_ is a bare
                // file name, and mpirun starts every rank in the case
                // directory, so all ranks used to open and truncate the SAME
                // file. Each non-master rank now writes inside its own
                // processorN. Every rank still constructs the stream - see the
                // deadlock note on postProcessingDir().
                const fileName protectionLogPath
                (
                    Pstream::parRun() && !Pstream::master()
                  ? fileName("processor" + Foam::name(Pstream::myProcNo()))
                        /protectionLogName_
                  : fileName(protectionLogName_)
                );

                OFstream os(protectionLogPath);
                os  << "# Embedded TNNP numerical protection summary" << nl
                    << "corrections " << protectionCorrections_ << nl
                    << "# mode "
                    << (numericalProtection_ ? "corrected" : "diagnostic only") << nl
                    << "# state reason count first_time_ms first_point"
                    << " first_phase first_value min_value max_value" << nl;
                forAll(protectionCount_, k)
                {
                    if (protectionCount_[k] == 0) continue;
                    os  << TNNP_STATES_NAMES[k/nProtectionReasons] << ' '
                        << protectionReasonName(k % nProtectionReasons) << ' '
                        << protectionCount_[k] << ' '
                        << protectionFirstTimeMs_[k] << ' '
                        << protectionFirstPoint_[k] << ' '
                        << protectionFirstPhase_[k] << ' '
                        << protectionFirstValue_[k] << ' '
                        << protectionMinValue_[k] << ' '
                        << protectionMaxValue_[k] << nl;
                }
            }
        }

        label nEqns() const
        {
            return NUM_STATES;
        }

        wordList exportedFieldNames() const
        {
            return exportedNames_;
        }

        // Initialise the ionic states on the point set the caller owns.
        //
        // The caller sizes the state array to the point set the ODE is
        // integrated on, and that is NOT always nPoints_: with
        // stateIntegrationMode = cellCentredReconstruct the ODE is driven by
        // the cell-centred Vm, so the array holds mesh.nCells() entries while
        // nPoints_ counts the Iion Gauss points. Resizing to nPoints_ here
        // used to make the size guard in solveODE()/calculateCurrent()
        // permanently unsatisfiable in that mode, so the guard re-fired on
        // every call and silently reset the whole cell model to its resting
        // state at every time step (no state history could ever accumulate).
        //
        // The per-point scratch (rates_, algebraic_) is sized in the
        // constructor to nPoints_, the largest point count the model is ever
        // driven with, and is left alone here.
        void initialiseStates(Field<Field<scalar>>& states)
        {
            if (states.empty())
            {
                states.setSize(nPoints_);
            }

            // The scratch is indexed by the point set being integrated, which
            // may be larger than states.size(); keep it covering all nPoints_
            // entries regardless of how many states the caller holds.
            forAll(rates_, pointI)
            {
                rates_[pointI].setSize(NUM_STATES, 0.0);
                algebraic_[pointI].setSize(NUM_ALGEBRAIC, 0.0);
            }

            forAll(states, pointI)
            {
                states[pointI].setSize(NUM_STATES, 0.0);
                TNNPinitConsts
                (
                    constants_.data(),
                    rates_[pointI].data(),
                    states[pointI].data(),
                    tissue_,
                    dict_
                );
                loadStimulusProtocol();
                computeVariables(0.0, states[pointI], rates_[pointI], algebraic_[pointI]);
            }
            oldStates_ = states;
        }

        void updateStatesOld(const Field<Field<scalar>>& states)
        {
            oldStates_ = states;
        }

        void resetStatesToStatesOld(Field<Field<scalar>>& states)
        {
            states = oldStates_;
        }

        void calculateCurrent
        (
            const scalar stepStartTime,
            const scalar,
            const scalarField& Vm,
            scalarField& Im,
            Field<Field<scalar>>& states
        )
        {
            const scalar tStartMs = 1000.0*stepStartTime;
            if (Im.size() != Vm.size())
            {
                FatalErrorInFunction
                    << "Im.size() != Vm.size()" << abort(FatalError);
            }

            // A mismatch here is a programming error, not something to repair
            // on the fly: re-initialising would discard the accumulated state
            // history. The state array must already be sized to the point set
            // it is integrated on (mesh.nCells() for cellCentredReconstruct,
            // the Iion Gauss points otherwise).
            if (states.size() != Vm.size())
            {
                FatalErrorInFunction
                    << "states.size() = " << states.size()
                    << " != Vm.size() = " << Vm.size() << nl
                    << "The ionic state array must be sized to the point set "
                    << "it is integrated on." << abort(FatalError);
            }

            forAll(Vm, pointI)
            {
                scalarField& y = states[pointI];
                y[V] = 1000.0*Vm[pointI];
                protectState(y, tStartMs, pointI, "calculateCurrent");
                computeVariables(tStartMs, y, rates_[pointI], algebraic_[pointI]);
                Im[pointI] = algebraic_[pointI][Iion_cm];
            }
        }

        void solveODE
        (
            const scalar stepStartTime,
            const scalar deltaT,
            const scalarField& Vm,
            scalarField& Im,
            Field<Field<scalar>>& states,
            const labelUList& scratchIDs = emptyScratchIDs(),
            const VmDriverStages* driver = nullptr
        )
        {
            const scalar tStartMs = 1000.0*stepStartTime;
            const scalar dtMs = 1000.0*deltaT;

            // Added for cardiacFoam: the Vm driver, when given, holds one field
            // per node except the last, which is Vm itself (the candidate). A
            // missing or mis-sized node is a silently wrong ODE driver rather
            // than a crash, so it is checked here, before the parallel loop.
            if (driver)
            {
                if (driver->nNodes < 2 || driver->nNodes > driver->s.size())
                {
                    FatalErrorInFunction
                        << "Vm driver with " << driver->nNodes << " nodes."
                        << abort(FatalError);
                }

                for (label j = 0; j < driver->nNodes - 1; ++j)
                {
                    if
                    (
                        !driver->field[j]
                     || driver->field[j]->size() != Vm.size()
                    )
                    {
                        FatalErrorInFunction
                            << "Vm driver node " << j << " of "
                            << driver->nNodes << " is not set or is not sized"
                            << " to the " << Vm.size() << " points integrated."
                            << abort(FatalError);
                    }
                }
            }
            if (Im.size() != Vm.size())
            {
                FatalErrorInFunction
                    << "Im.size() != Vm.size()" << abort(FatalError);
            }

            // Added for cardiacFoam: the caller may pass an explicit scratch
            // slot for each point. This is what keeps the front-aware hybrid
            // correct. That path calls solveODE twice per iteration, once with
            // a compacted array of nF front Gauss points and once with a
            // compacted array of nS smooth cells, both indexed from 0. Without
            // a slot map both passes write algebraic_/rates_[0..n),
            // so the second silently overwrites the first's scratch even
            // though the two arrays describe different physical points.
            if (!scratchIDs.empty())
            {
                if (scratchIDs.size() != Vm.size())
                {
                    FatalErrorInFunction
                        << "scratchIDs.size() = " << scratchIDs.size()
                        << " != Vm.size() = " << Vm.size()
                        << abort(FatalError);
                }

                forAll(scratchIDs, i)
                {
                    if (scratchIDs[i] < 0 || scratchIDs[i] >= rates_.size())
                    {
                        FatalErrorInFunction
                            << "scratchIDs[" << i << "] = " << scratchIDs[i]
                            << " out of range [0, " << rates_.size() << ")"
                            << abort(FatalError);
                    }
                }
            }

            // A mismatch here is a programming error, not something to repair
            // on the fly: re-initialising would discard the accumulated state
            // history. The state array must already be sized to the point set
            // it is integrated on (mesh.nCells() for cellCentredReconstruct,
            // the Iion Gauss points otherwise).
            if (states.size() != Vm.size())
            {
                FatalErrorInFunction
                    << "states.size() = " << states.size()
                    << " != Vm.size() = " << Vm.size() << nl
                    << "The ionic state array must be sized to the point set "
                    << "it is integrated on." << abort(FatalError);
            }

            // Modified for cardiacFoam: OpenMP across the points.
            //
            // Safe because every iteration touches only data of its own point:
            // states[pointI] and Im[pointI] are distinct entries, and the
            // scratch is addressed by SLOT, with the slot map injective (the
            // identity when absent, and the front/smooth maps of the hybrid are
            // disjoint by construction - see the two scratch banks above). The
            // only shared mutable state on this path is the protection counter,
            // which recordCorrection() now serialises.
            //
            // Written as a canonical for-loop rather than forAll because OpenMP
            // cannot parallelise the macro's form. Iteration order does not
            // affect the result: the points do not communicate.
            //
            // Caveat worth stating: a FatalError raised inside the loop escapes
            // a parallel region, which OpenMP does not define. It aborts the
            // process either way, but the message may interleave. The MMS
            // solver's OpenMP loops carry the same exposure.
            const label nPointsLocal = Vm.size();

            #ifdef _OPENMP
            #pragma omp parallel for schedule(static) \
                if (parallelODE_ && nPointsLocal >= parallelODEThreshold_)
            #endif
            for (label pointI = 0; pointI < nPointsLocal; ++pointI)
            {
                // The scratch is addressed by slot, not by position in the
                // passed array. With no slot map the two coincide, so every
                // other call site is unaffected.
                const label slotI =
                    scratchIDs.empty() ? pointI : scratchIDs[pointI];

                // Per-point nodes, local to the iteration, so the OpenMP loop
                // shares nothing.
                VmDriverNodes nodes;
                if (driver)
                {
                    nodes.nNodes = driver->nNodes;
                    for (label j = 0; j < driver->nNodes - 1; ++j)
                    {
                        nodes.s[j] = driver->s[j];
                        nodes.v[j] = (*driver->field[j])[pointI];
                    }
                    nodes.s[driver->nNodes - 1] = driver->s[driver->nNodes - 1];
                    nodes.v[driver->nNodes - 1] = Vm[pointI];
                }

                advancePoint
                (
                    tStartMs,
                    dtMs,
                    Vm[pointI],
                    driver ? &nodes : nullptr,
                    states[pointI],
                    algebraic_[slotI],
                    rates_[slotI],
                    slotI
                );
                Im[pointI] = algebraic_[slotI][Iion_cm];
            }
        }

        void exportStates
        (
            const Field<Field<scalar>>& states,
            PtrList<volScalarField>& outFields
        )
        {
            forAll(outFields, outI)
            {
                const word& name = exportedNames_[outI];
                const label sI = stateIndex(name);
                const label aI = algebraicIndex(name);
                scalarField& fld = outFields[outI].primitiveFieldRef();

                forAll(fld, cellI)
                {
                    if (sI >= 0 && cellI < states.size())
                    {
                        fld[cellI] = states[cellI][sI];
                    }
                    else if (aI >= 0 && cellI < algebraic_.size())
                    {
                        fld[cellI] = algebraic_[cellI][aI];
                    }
                    else
                    {
                        fld[cellI] = 0.0;
                    }
                }
                outFields[outI].correctBoundaryConditions();
            }
        }

        void exportStatesIntegrationPoints
        (
            const Field<Field<scalar>>& states,
            PtrList<volScalarField>& outFields,
            const CompactListList<scalar>& cellIionQuadW
        )
        {
            forAll(outFields, outI)
            {
                const word& name = exportedNames_[outI];
                const label sI = stateIndex(name);
                const label aI = algebraicIndex(name);
                scalarField& fld = outFields[outI].primitiveFieldRef();

                label integrationPointI = 0;
                forAll(cellIionQuadW, cellI)
                {
                    scalar wSum = 0.0;
                    scalar value = 0.0;
                    forAll(cellIionQuadW[cellI], qI)
                    {
                        const scalar w = cellIionQuadW[cellI][qI];
                        wSum += w;
                        if (sI >= 0 && integrationPointI < states.size())
                        {
                            value += w*states[integrationPointI][sI];
                        }
                        else if (aI >= 0 && integrationPointI < algebraic_.size())
                        {
                            value += w*algebraic_[integrationPointI][aI];
                        }
                        ++integrationPointI;
                    }
                    fld[cellI] = value/max(wSum, SMALL);
                }
                outFields[outI].correctBoundaryConditions();
            }
        }
    };

    // Representative cell size: the geometric mean edge of a cell if the domain
    // were divided uniformly, i.e. (active volume / total cells)^(1/dim). Used
    // only to size the stable time step, so a global average is the right
    // notion. The cell count is reduced, so the value does not depend on the
    // decomposition.
    scalar characteristicDx(const fvMesh& mesh)
    {
        const boundBox bb(mesh.points());
        const vector span = bb.max() - bb.min();
        const label dim = max(mesh.nGeometricD(), label(1));

        scalar activeMeasure = 1.0;
        if (dim >= 1)
        {
            activeMeasure *= max(span.x(), SMALL);
        }
        if (dim >= 2)
        {
            activeMeasure *= max(span.y(), SMALL);
        }
        if (dim >= 3)
        {
            activeMeasure *= max(span.z(), SMALL);
        }

        const scalar nCellsGlobal = scalar(returnReduce(mesh.nCells(), sumOp<label>()));

        return std::pow(activeMeasure/max(nCellsGlobal, scalar(1.0)), 1.0/scalar(dim));
    }

    // Explicit-diffusion stability limit dt = CFL*dx^2/(dim*D_eff), with
    // D_eff = max(D)/(chi*Cm). Reported (and optionally used) as a reference
    // even though the solver is implicit, because it is the scale at which the
    // nonlinear coupling stops being stiff. Reduced with min so all ranks agree.
    scalar computeStableDeltaT
    (
        const volTensorField& conductivity,
        const scalar chiVal,
        const scalar CmVal,
        const scalar CFL,
        const scalar dx,
        const label dim
    )
    {
        const scalar Dxx = max(conductivity.component(tensor::XX)).value();
        const scalar Dyy = max(conductivity.component(tensor::YY)).value();
        const scalar Dzz = max(conductivity.component(tensor::ZZ)).value();
        const scalar Dmax = max(Dxx, max(Dyy, Dzz));
        const scalar DefMax = Dmax/(chiVal*CmVal);

        scalar dtStable = GREAT;
        if (dim > 0 && DefMax > SMALL)
        {
            dtStable = CFL*sqr(dx)/(scalar(dim)*DefMax);
        }

        reduce(dtStable, minOp<scalar>());
        return dtStable;
    }

    // Tissue-level external stimulus over the step [t0, t0 + dt]: in the cells
    // of every stimulus region, stimulusIntensity times the FRACTION of the step
    // that falls inside the region's window [tStart, tStart + duration], and
    // zero everywhere else. Regions and their start times come from
    // constant/stimulusProtocol. Returns the fraction per region, which is
    // what the BDF restart compares from one step to the next.
    //
    // Modified for cardiacFoam: the current used to be switched on for the
    // WHOLE step whenever t0 lay in the CLOSED window. The step starting exactly
    // at tStart + duration was therefore still stimulated, so a 2 ms stimulus
    // at dt = 0.1 ms injected the charge of 2.1 ms: an O(dt) error in the
    // stimulus itself, whatever the time scheme, and a first-order floor under
    // every temporal convergence study. Averaging over the step makes the
    // injected charge exact for any dt and any alignment: with the window on
    // step boundaries it is exactly the windowed current, and a window that
    // ends inside a step contributes that step's covered fraction. The value
    // is constant within the step, so every scheme - the theta family, BDF and
    // the ESDIRK stages - sees the same forcing.
    //
    // Where two windows cover the same cell the larger fraction wins, which is
    // the old "set, not add" semantics carried over.
    List<scalar> applyStimulus
    (
        const scalar t0,
        const scalar dt,
        volScalarField& externalStimulusCurrent,
        const List<labelList>& stimulusCellIDsList,
        const List<scalar>& stimulusStartTimes,
        const scalar stimulusIntensity,
        const scalar stimulusDuration
    )
    {
        scalarField& extI = externalStimulusCurrent.primitiveFieldRef();
        extI = 0.0;

        List<scalar> fraction(stimulusCellIDsList.size(), 0.0);

        forAll(stimulusCellIDsList, bI)
        {
            const scalar tStart = stimulusStartTimes[bI];
            const scalar tEnd = tStart + stimulusDuration;
            const scalar overlap = min(t0 + dt, tEnd) - max(t0, tStart);

            // Modified for cardiacFoam: t0 is runTime.value(), which accumulates
            // round-off, so a step meant to start (or end) exactly on a window
            // edge can overlap it by a sliver. That sliver changed the fraction
            // by more than the 1e-12 that triggers a BDF history discard, so
            // the window's end discarded the history TWICE: measured in E02 at
            // dt = 6.25e-6 only (t = 0.002 and 0.00200625), which gave BDF4
            // three extra ESDIRK3 steps on the finest dt of a convergence
            // ladder. Overlaps within 1e-9 of the step - the tolerance
            // clippedStep uses - are snapped to 0 or to the full step. Any
            // other fraction is unchanged, bit for bit.
            if (dt <= 0.0 || overlap <= 1.0e-9*dt)
            {
                continue;
            }

            fraction[bI] =
                overlap >= (1.0 - 1.0e-9)*dt ? scalar(1) : overlap/dt;

            const scalar current = fraction[bI]*stimulusIntensity;
            const labelList& stimulusCellIDs = stimulusCellIDsList[bI];
            forAll(stimulusCellIDs, cI)
            {
                const label cellI = stimulusCellIDs[cI];
                if (mag(current) > mag(extI[cellI]))
                {
                    extI[cellI] = current;
                }
            }
        }

        externalStimulusCurrent.correctBoundaryConditions();

        return fraction;
    }


    // Worst relative L2 change over all state variables between two nonlinear
    // iterations, also filling the per-state residuals for reporting. The
    // states converge at different rates, so the stopping test uses the maximum
    // rather than an average.
    scalar maxStateRelativeL2Difference
    (
        const Field<Field<scalar>>& current,
        const Field<Field<scalar>>& previous,
        scalarField& stateResiduals
    )
    {
        stateResiduals = 0.0;

        // Modified for cardiacFoam: this function is COLLECTIVE - it reduces
        // twice per state below - so every rank has to reach every reduction,
        // and reach it the same number of times.
        //
        // What was here before returned 0.0 early when a rank's state array came
        // up empty, which is exactly what happens to a rank whose partition has
        // no cells. That rank skipped both reductions while the others blocked
        // in them: not a wrong number, a HANG. Reproduced on a 2-cell mesh with
        // `simple (4 1 1)` and with `scotch` at np=8 - both killed at the
        // timeout, while the same case at np=1 finishes in one second (test
        // E19).
        //
        // So the empty case is now handled by letting the accumulation loops
        // below contribute nothing: `forAll(current, ...)` over an empty field
        // iterates zero times, num and den stay 0, and the rank still takes part
        // in every reduce.
        //
        // The trip count needs the same treatment. `current[0].size()` is
        // rank-local and cannot even be read on an empty rank, so it is reduced:
        // otherwise two ranks could run this loop a different number of times
        // and desynchronise the reductions for the same reason.
        label nStates = stateResiduals.size();
        if (!current.empty())
        {
            nStates = min(current[0].size(), stateResiduals.size());
        }
        reduce(nStates, minOp<label>());

        scalar maxResidual = 0.0;

        for (label stateI = 0; stateI < nStates; ++stateI)
        {
            scalar num = 0.0;
            scalar den = 0.0;

            forAll(current, pointI)
            {
                if (stateI >= current[pointI].size() || stateI >= previous[pointI].size())
                {
                    continue;
                }

                const scalar c = current[pointI][stateI];
                const scalar diff = c - previous[pointI][stateI];
                num += diff*diff;
                den += c*c;
            }

            reduce(num, sumOp<scalar>());
            reduce(den, sumOp<scalar>());

            stateResiduals[stateI] = std::sqrt(num)/(std::sqrt(den) + SMALL);
            maxResidual = max(maxResidual, stateResiduals[stateI]);
        }

        return maxResidual;
    }

    // Added for cardiacFoam: directory for this run's postProcessing output.
    //
    // In a parallel run runTime.path() is <case>/processorN, so every one of
    // these files used to be written inside a processor directory where the
    // tutorial driver does not look for them.
    //
    // EVERY rank must still construct the stream. Guarding with
    // "if (!Pstream::master()) return;" before opening DEADLOCKS: OpenFOAM's
    // file handler can be collective, so the master blocks inside the
    // collective open while the other ranks have already moved on. That failure
    // was diagnosed in the MMS solver as a run which completed every time step,
    // printed its timings and then hung. Only the destination differs here.
    fileName postProcessingDir(const Time& runTime)
    {
        const fileName base =
        (
            Pstream::master()
          ? runTime.rootPath()/runTime.globalCaseName()
          : runTime.path()
        );

        return base/"postProcessing"/"highOrderElectroActivationFoamImplicitPETScDistributed";
    }











    // The cell-centred ionic states at the Iion quadrature points, all states
    // in one walk of the stencil, through the shared reconstruction core
    // (monodomainReconstruction.H):
    //
    //   statesCells [cellI][stateI]  ->  statesIP [ipI][stateI]
    //
    // Only the cells with writeCellMask true are written (all when it is
    // null). limitToNeighbourRange applies the Barth-Jespersen-style limiter
    // (stateReconstructLimiter), which the stiff TNNP states need across an
    // under-resolved front.
    void reconstructStatesAtIionIntegrationPoints
    (
        const fvMesh& mesh,
        const Field<Field<scalar>>& statesCells,
        const label nStates,
        const highOrderInterp& LREInterp_states,
        const highOrderInterp& LREInterp_Iion,
        Field<Field<scalar>>& statesIP,
        const boolList* writeCellMask,
        const bool limitToNeighbourRange,
        const Switch useOpenMP,
        const label openMPThreshold
    )
    {
        const CompactListList<point>& qp = LREInterp_Iion.cellQuadPoints();

        // One contiguous column per state: the core reads each field by cell.
        List<scalarField> columns(nStates, scalarField(statesCells.size()));
        List<const scalarField*> fields(nStates);
        for (label stateI = 0; stateI < nStates; ++stateI)
        {
            scalarField& col = columns[stateI];
            forAll(col, cellI)
            {
                col[cellI] = statesCells[cellI][stateI];
            }
            fields[stateI] = &col;
        }

        auto write = [&statesIP](const label ip, const label k, const scalar v)
        {
            statesIP[ip][k] = v;
        };

        reconstructCellFieldsAtIntegrationPoints
        (
            mesh,
            fields,
            LREInterp_states,
            qp,
            write,
            writeCellMask,
            limitToNeighbourRange,
            useOpenMP,
            openMPThreshold
        );
    }


    // This is the discrete cell-average operator and it is the natural
    // companion to reconstructVmAtIionIntegrationPoints(): together they
    // form a Gauss-quadrature evaluation of  (1/|Omega|) int_Omega Iion(Vm).
    void averageIionIntegrationPointsToCells
    (
        const scalarField& IionIntegrationPoints,
        const highOrderInterp& LREInterp_Iion,
        volScalarField& Iion
    )
    {
        const fvMesh& mesh = Iion.mesh();
        const CompactListList<scalar>& cellIionQuadW =
            LREInterp_Iion.cellQuadWeightPhysical();

        scalarField& IionCells = Iion.primitiveFieldRef();
        label integrationPointI = 0;

        forAll(mesh.cells(), cellI)
        {
            scalar iBar = 0.0;
            scalar wSum = 0.0;

            forAll(cellIionQuadW[cellI], qI)
            {
                const scalar w = cellIionQuadW[cellI][qI];
                iBar += w*IionIntegrationPoints[integrationPointI];
                wSum += w;
                ++integrationPointI;
            }

            // Normalise by sum of weights so the result is an unbiased
            // average regardless of whether quadrature weights sum to 1.
            IionCells[cellI] = iBar/max(wSum, SMALL);
        }

        Iion.correctBoundaryConditions();
    }

    // Quadrature-average the per-Gauss-point state vectors to cell centres,
    // component by component. Companion to averageIionIntegrationPointsToCells
    // for the Field<Field<scalar>> state carrier; used by the front-aware
    // hybrid to build cell-centred states from the Gauss-point carrier.
    void averageStatesIntegrationPointsToCells
    (
        const Field<Field<scalar>>& statesIP,
        const label nStates,
        const highOrderInterp& LREInterp_Iion,
        Field<Field<scalar>>& statesCells,
        const boolList* cellMask = nullptr
    )
    {
        const CompactListList<scalar>& cellIionQuadW =
            LREInterp_Iion.cellQuadWeightPhysical();

        label integrationPointI = 0;
        forAll(cellIionQuadW, cellI)
        {
            const bool updateThisCell =
                (cellMask == nullptr) || (*cellMask)[cellI];

            if (!updateThisCell)
            {
                integrationPointI += cellIionQuadW[cellI].size();
                continue;
            }

            Field<scalar>& sc = statesCells[cellI];
            sc = 0.0;
            scalar wSum = 0.0;

            forAll(cellIionQuadW[cellI], qI)
            {
                const scalar w = cellIionQuadW[cellI][qI];
                const Field<scalar>& sq = statesIP[integrationPointI];
                for (label k = 0; k < nStates; ++k)
                {
                    sc[k] += w*sq[k];
                }
                wSum += w;
                ++integrationPointI;
            }

            const scalar invW = 1.0/max(wSum, SMALL);
            for (label k = 0; k < nStates; ++k)
            {
                sc[k] *= invW;
            }
        }
    }

    // Flag "front" cells for the front-aware hybrid: those whose sub-cell Vm
    // range (over their Iion Gauss points) exceeds a threshold, or which are
    // stimulated. This classification is computed ONCE per time step (from the
    // committed Vm) and frozen for the whole nonlinear solve -- recomputing it
    // per nonlinear iteration makes the front/smooth switch flip between
    // iterations and turns the Picard/Newton fixed point into a period-2 limit
    // cycle that never converges.
    void flagFrontCells
    (
        const fvMesh& mesh,
        const scalarField& VmIntegrationPoints,
        const highOrderInterp& LREInterp_Iion,
        const volScalarField& externalStimulusCurrent,
        const scalar vmRangeThreshold,
        boolList& frontCell
    )
    {
        const CompactListList<point>& cellQP = LREInterp_Iion.cellQuadPoints();
        label ip = 0;
        forAll(mesh.cells(), cellI)
        {
            scalar vMin = GREAT, vMax = -GREAT;
            forAll(cellQP[cellI], qI)
            {
                const scalar v = VmIntegrationPoints[ip];
                vMin = min(vMin, v);
                vMax = max(vMax, v);
                ++ip;
            }
            frontCell[cellI] =
                (vMax - vMin > vmRangeThreshold)
             || (mag(externalStimulusCurrent[cellI]) > SMALL);
        }
    }

    // Grow the set of front cells by nLayers cell layers (face neighbours),
    // so the moving front stays inside the Gauss-resolved band for the whole
    // time step and the front/smooth transition is padded.
    // Added for cardiacFoam: opt-in front-cell count diagnostic.
    bool debugFrontCellCount = false;

    // Grow the front flag by nLayers layers of face neighbours, so that the
    // expensive per-Gauss-point ODE integration is applied to a band around the
    // propagating front rather than exactly on it.
    void dilateFrontCells
    (
        const fvMesh& mesh,
        const label nLayers,
        boolList& frontCell
    )
    {
        if (nLayers <= 0)
        {
            return;
        }

        const labelListList& cellCells = mesh.cellCells();

        // Added for cardiacFoam: marker field used to carry the front flag
        // across processor boundaries. mesh.cellCells() stops at a processor
        // face, so without this the front band simply failed to grow across a
        // partition. Unlike the stabilisation gap this is an O(1) error, not a
        // shrinking surface effect: a cell classified front integrates a
        // different set of ODEs from one classified smooth.
        //
        // Constructing the field with zeroGradient does NOT put zeroGradient on
        // the processor patches: fvPatchField::New overrides the requested type
        // with the patch's own constraint type, so processor patches receive a
        // real processorFvPatchScalarField and correctBoundaryConditions()
        // performs the halo exchange. One exchange per dilation layer.
        autoPtr<volScalarField> markerPtr;
        if (Pstream::parRun())
        {
            markerPtr.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        "frontCellMarker",
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar("zero", dimless, 0.0),
                    zeroGradientFvPatchScalarField::typeName
                )
            );
        }

        for (label layer = 0; layer < nLayers; ++layer)
        {
            boolList grown(frontCell);
            forAll(frontCell, cellI)
            {
                if (!frontCell[cellI])
                {
                    continue;
                }
                const labelList& nbr = cellCells[cellI];
                forAll(nbr, j)
                {
                    grown[nbr[j]] = true;
                }
            }

            if (Pstream::parRun())
            {
                volScalarField& marker = markerPtr();
                scalarField& m = marker.primitiveFieldRef();
                forAll(frontCell, cellI)
                {
                    m[cellI] = frontCell[cellI] ? 1.0 : 0.0;
                }
                marker.correctBoundaryConditions();

                forAll(mesh.boundary(), patchI)
                {
                    const fvPatch& p = mesh.boundary()[patchI];
                    if (!p.coupled())
                    {
                        continue;
                    }

                    const scalarField nbrMarker
                    (
                        marker.boundaryField()[patchI].patchNeighbourField()
                    );
                    const labelUList& faceCells = p.faceCells();

                    forAll(faceCells, i)
                    {
                        if (nbrMarker[i] > 0.5)
                        {
                            grown[faceCells[i]] = true;
                        }
                    }
                }
            }

            frontCell = grown;
        }

        // Added for cardiacFoam: the sharpest available check that the front
        // classification is partition-independent. It is an integer, so it must
        // match EXACTLY between rank counts; before the exchange above it
        // differed by roughly the number of cells adjacent to a cut.
        if (debugFrontCellCount)
        {
            label nFront = 0;
            forAll(frontCell, cellI)
            {
                if (frontCell[cellI]) ++nFront;
            }
            reduce(nFront, sumOp<label>());
            Info<< "    [front] cells = " << nFront << endl;
        }
    }

    // Detect the first time Vm crosses the activation threshold in each
    // cell and record it via linear interpolation across the time step:
    //
    //     w = (Vthr - Vm^n) / (Vm^{n+1} - Vm^n)   in [0, 1]
    //     activationTime = t_n + w * dt
    //
    // Once a cell is marked activated its flag in calculateActivationTime
    // is cleared so subsequent threshold re-crossings (e.g. recovery and
    // re-excitation) are ignored.
    void updateActivationTimes
    (
        const scalarField& VmOld,
        const volScalarField& Vm,
        volScalarField& activationTime,
        boolList& calculateActivationTime,
        const scalar t0,
        const scalar dt,
        const scalar activationThreshold
    )
    {
        const scalarField& VmI = Vm.primitiveField();
        scalarField& activationTimeI = activationTime.primitiveFieldRef();

        forAll(activationTimeI, cellI)
        {
            if (!calculateActivationTime[cellI])
            {
                continue;
            }

            if (VmOld[cellI] < activationThreshold && VmI[cellI] >= activationThreshold)
            {
                const scalar denom = VmI[cellI] - VmOld[cellI];
                scalar w = 0.0;
                if (mag(denom) > VSMALL)
                {
                    w = (activationThreshold - VmOld[cellI])/denom;
                }

                w = min(max(w, scalar(0.0)), scalar(1.0));
                activationTimeI[cellI] = t0 + w*dt;
                calculateActivationTime[cellI] = false;
            }
        }

        activationTime.correctBoundaryConditions();
    }

    // Index of the cell whose centre is nearest a sample point, or -1 on the
    // ranks that do not own it. Used to locate benchmark probe points.
    label nearestCellToPoint
    (
        const fvMesh& mesh,
        const point& samplePoint
    )
    {
        const vectorField& C = mesh.C();
        label bestCellI = -1;
        scalar bestD2 = GREAT;

        forAll(C, cellI)
        {
            const scalar d2 = magSqr(C[cellI] - samplePoint);
            if (d2 < bestD2)
            {
                bestD2 = d2;
                bestCellI = cellI;
            }
        }

        // Modified for cardiacFoam: the search above is over LOCAL cells only,
        // so every rank used to believe it owned the sample point. Downstream
        // that meant each rank triggered the stop condition on ITS OWN nearest
        // cell, clipped effectiveEndTime differently, and left the time loop at
        // a different step - some ranks dropping out of the collective solve
        // while the others blocked in it. A hang, preceded by corrupted physics,
        // because effectiveEndTime also sets the final step's dt.
        //
        // Exactly one rank now returns a cell; the others return -1. Ties are
        // broken by lowest rank, so the choice is deterministic and independent
        // of the decomposition. reduce() is preferred over Pstream::broadcast,
        // which broadcasts from masterNo() rather than from an arbitrary root.
        if (Pstream::parRun())
        {
            scalar globalBestD2 = bestD2;
            reduce(globalBestD2, minOp<scalar>());

            label winner =
                (bestCellI >= 0 && bestD2 <= globalBestD2)
              ? Pstream::myProcNo()
              : Pstream::nProcs();
            reduce(winner, minOp<label>());

            if (winner >= Pstream::nProcs())
            {
                FatalErrorInFunction
                    << "Could not find a cell for point " << samplePoint
                    << " on any rank" << exit(FatalError);
            }

            return (winner == Pstream::myProcNo()) ? bestCellI : -1;
        }

        if (bestCellI < 0)
        {
            FatalErrorInFunction
                << "Could not find a cell for point " << samplePoint
                << exit(FatalError);
        }

        return bestCellI;
    }


    // Sample a field at an arbitrary point by inverse-distance weighting over
    // its k nearest cell centres, with an exact-hit shortcut. This is how the
    // Niederer benchmark reads activation times at points that are not cell
    // centres.
    scalar sampleIDW
    (
        const volScalarField& field,
        const point& samplePoint,
        const label requestedNeighbours,
        const scalar requestedPower
    )
    {
        const vectorField& C = field.mesh().C();
        const scalarField& values = field.primitiveField();

        // Modified for cardiacFoam: the k nearest cells must be the k nearest
        // GLOBALLY. The search below is over local cells, so on N ranks this
        // used to return "the k nearest on my rank" - N different answers, none
        // of them right unless the sample point happened to sit well inside one
        // partition.
        //
        // Each rank builds its own top-k, the candidates are gathered, and every
        // rank merges the same set with the same ordering, so all ranks agree by
        // construction. Ties are broken by global cell ID, which reproduces the
        // serial result exactly: the serial loop tests d2 < bestD2 strictly, so
        // among equal distances it keeps the LOWEST cell index, and in serial
        // the global ID is the cell index.
        //
        // Note the early return on an exact hit has moved AFTER the merge. Left
        // inside the search loop it would make one rank return immediately while
        // the others entered the gather, which deadlocks.
        const label nGlobalCells = returnReduce(C.size(), sumOp<label>());
        if (nGlobalCells == 0)
        {
            return 0.0;
        }

        const label nNbrs =
            min(max(requestedNeighbours, label(1)), nGlobalCells);
        const scalar power = max(requestedPower, scalar(SMALL));

        List<scalar> bestD2(nNbrs, GREAT);
        List<label> bestId(nNbrs, -1);
        List<scalar> bestVal(nNbrs, 0.0);

        forAll(C, cellI)
        {
            const scalar d2 = magSqr(C[cellI] - samplePoint);

            for (label i = 0; i < nNbrs; ++i)
            {
                if (d2 < bestD2[i])
                {
                    for (label j = nNbrs - 1; j > i; --j)
                    {
                        bestD2[j] = bestD2[j - 1];
                        bestId[j] = bestId[j - 1];
                        bestVal[j] = bestVal[j - 1];
                    }

                    bestD2[i] = d2;
                    bestId[i] = gCol(cellI);
                    bestVal[i] = values[cellI];
                    break;
                }
            }
        }

        // Added for cardiacFoam: merge every rank's candidates identically.
        if (Pstream::parRun())
        {
            List<List<scalar>> allD2(Pstream::nProcs());
            List<List<label>> allId(Pstream::nProcs());
            List<List<scalar>> allVal(Pstream::nProcs());
            allD2[Pstream::myProcNo()] = bestD2;
            allId[Pstream::myProcNo()] = bestId;
            allVal[Pstream::myProcNo()] = bestVal;
            Pstream::allGatherList(allD2);
            Pstream::allGatherList(allId);
            Pstream::allGatherList(allVal);

            DynamicList<label> order;
            DynamicList<scalar> d2Flat;
            DynamicList<label> idFlat;
            DynamicList<scalar> valFlat;
            forAll(allD2, procI)
            {
                forAll(allD2[procI], i)
                {
                    if (allId[procI][i] < 0)
                    {
                        continue;
                    }
                    order.append(d2Flat.size());
                    d2Flat.append(allD2[procI][i]);
                    idFlat.append(allId[procI][i]);
                    valFlat.append(allVal[procI][i]);
                }
            }

            // Sort by (distance, global cell ID). The ID makes the order total,
            // so every rank produces the same list.
            labelList idx(order);
            std::stable_sort
            (
                idx.begin(),
                idx.end(),
                [&](const label a, const label b)
                {
                    if (d2Flat[a] != d2Flat[b]) return d2Flat[a] < d2Flat[b];
                    return idFlat[a] < idFlat[b];
                }
            );

            bestD2 = List<scalar>(nNbrs, GREAT);
            bestId = List<label>(nNbrs, -1);
            bestVal = List<scalar>(nNbrs, 0.0);
            const label nKeep = min(nNbrs, label(idx.size()));
            for (label i = 0; i < nKeep; ++i)
            {
                bestD2[i] = d2Flat[idx[i]];
                bestId[i] = idFlat[idx[i]];
                bestVal[i] = valFlat[idx[i]];
            }
        }

        // Exact hit on a cell centre: return that cell's value. Deferred to
        // here so that every rank has already taken part in the gather.
        if (bestId[0] >= 0 && bestD2[0] < VSMALL)
        {
            return bestVal[0];
        }

        scalar wSum = 0.0;
        scalar vSum = 0.0;
        forAll(bestId, i)
        {
            if (bestId[i] < 0)
            {
                continue;
            }

            const scalar d = std::sqrt(max(bestD2[i], scalar(VSMALL)));
            const scalar w = 1.0/std::pow(d, power);
            wSum += w;
            vSum += w*bestVal[i];
        }

        return vSum/max(wSum, scalar(SMALL));
    }

    // Write the benchmark activation-time outputs: the value at the P8 corner
    // point and a profile sampled along the diagonal of the slab, both by
    // inverse-distance weighting, into the solver's postProcessing directory.
    void writeActivationSamples
    (
        const Time& runTime,
        const volScalarField& activationTime,
        const point& diagonalStart,
        const point& diagonalEnd,
        const point& p8Point,
        const label nDiagonalSamplesRequested,
        const label sampleNeighbours,
        const scalar samplePower,
        const bool verbose = true
    )
    {
        // Modified for cardiacFoam: see postProcessingDir().
        const fileName sampleDir(postProcessingDir(runTime));
        mkDir(sampleDir);

        const label nDiagonalSamples = max(nDiagonalSamplesRequested, label(2));

        const scalar p8Time = sampleIDW
        (
            activationTime,
            p8Point,
            sampleNeighbours,
            samplePower
        );

        OFstream p8File(sampleDir/"P8_activationTime.dat");
        p8File<< "# label x y z activationTime_s activationTime_ms" << nl;
        p8File<< "P8 "
              << p8Point.x() << ' ' << p8Point.y() << ' ' << p8Point.z() << ' '
              << p8Time << ' ' << 1000.0*p8Time << nl;

        OFstream diagFile(sampleDir/"diagonalActivationTime.dat");
        diagFile<< "# index x y z distance_m distance_mm activationTime_s activationTime_ms" << nl;

        const vector diag = diagonalEnd - diagonalStart;
        for (label i = 0; i < nDiagonalSamples; ++i)
        {
            const scalar s = scalar(i)/scalar(nDiagonalSamples - 1);
            const point p = diagonalStart + s*diag;
            const scalar distance = mag(p - diagonalStart);
            const scalar activation = sampleIDW
            (
                activationTime,
                p,
                sampleNeighbours,
                samplePower
            );

            diagFile<< i << ' '
                    << p.x() << ' ' << p.y() << ' ' << p.z() << ' '
                    << distance << ' ' << 1000.0*distance << ' '
                    << activation << ' ' << 1000.0*activation << nl;
        }

        if (verbose)
        {
            Info<< "P8 activation time = " << p8Time << " s ("
                << 1000.0*p8Time << " ms)" << nl
                << "Activation samples written to " << sampleDir << nl;
        }
    }
}

int main(int argc, char* argv[])
{
#ifdef __GLIBC__
    mallopt(M_ARENA_MAX, 2);
    mallopt(M_MMAP_THRESHOLD, 64*1024);
    mallopt(M_TRIM_THRESHOLD, 64*1024);
#endif

    #include "setRootCaseLists.H"
    PetscSession petscSession(argc, argv);
    #include "createTime.H"
    #include "createMesh.H"
    // Global cell addressing and row offset, before anything assembles a
    // matrix (shared code, see monodomainBase.H).
    initialiseGlobalCells(mesh);

    if (Pstream::parRun())
    {
        Info<< "Distributed: " << gCellsPtr->totalSize() << " cells over "
            << Pstream::nProcs() << " ranks" << endl;
    }

    #include "createFields.H"

    // Added for cardiacFoam: reject Vm boundary conditions the assembly cannot
    // represent, at startup, instead of producing a silently wrong answer.
    //
    // The patch dispatchers key off the Vm patch field's TYPE NAME and only
    // branch on the four types below. Anything else used to fall out of the
    // bottom of the `if` - and it did so in four different ways, because the
    // orthogonal dispatcher has no `else` at all while the high-order one drops
    // into the Dirichlet-style flux loop. Measured before this check:
    // `fixedGradient` with a 50 V/m applied flux and `calculated` both produced
    // a trace IDENTICAL to `zeroGradient`, exit 0, no warning. The requested
    // boundary condition never reached the operator.
    //
    // Checked on every rank: `processor` patches only exist in parallel, and a
    // rank-local list would disagree between ranks.
    //
    // Modified for cardiacFoam: fixedValue and fixedVoltage are no longer in the
    // list. The operators assemble their diagonal term -a/V, so the boundary IS
    // clamped - but to zero, because the +a*Vm_b/V lift on the right-hand side
    // is never formed: Vm.boundaryField() is only ever read for its type, so the
    // patch VALUE cannot reach the system. Measured on the low-order path with
    // walls at -40 mV: max(Vm) -0.0198 V at step 5 against -0.0667 insulated,
    // 47 mV of error and growing, exit 0, no warning. On the high-order path the
    // run died inside solids4foam instead, an accidental abort rather than a
    // check. A comment in the orthogonal assembly claimed the boundary value
    // "enters through the RHS elsewhere"; there was no elsewhere. Refused here
    // until the lift is implemented (the MMS solver has one, built from its
    // exact solution, which is the natural model for a port).
    {
        const wordHashSet supportedVmPatchTypes
        ({
            word("empty"),
            zeroGradientFvPatchScalarField::typeName,
            word("processor")
        });

        wordHashSet offending;
        forAll(mesh.boundary(), patchI)
        {
            const fvPatch& patch = mesh.boundary()[patchI];
            const word bcType =
                patch.lookupPatchField<volScalarField, scalar>("Vm").type();

            if (!supportedVmPatchTypes.found(bcType))
            {
                offending.insert(bcType);
            }
        }

        if (!offending.empty())
        {
            FatalErrorInFunction
                << "Unsupported boundary condition(s) on Vm: "
                << offending.sortedToc() << nl
                << "The assembly only implements " << supportedVmPatchTypes.sortedToc()
                << nl
                << "Any other type is silently ignored by the operator, which"
                << " produces a plausible but wrong answer rather than an error."
                << nl
                << "fixedValue/fixedVoltage are refused because the Dirichlet"
                << " value never reaches the right-hand side: the boundary would"
                << " be clamped to 0 V, not to the value requested."
                << exit(FatalError);
        }
    }

    // Added for cardiacFoam: report the chunk size and how many chunks it
    // actually produces on THIS rank. Without the second number a chunking test
    // cannot tell whether it exercised the multi-chunk path at all - the count
    // is local, so a mesh that needs two chunks in serial can need one per rank
    // in parallel.
    Info<< "Assembly chunking     : faceChunk = " << assemblyFaceChunk << endl;
    Pout<< "Assembly chunking     : nInternalFaces = " << mesh.nInternalFaces()
        << " -> " << (mesh.nInternalFaces() + assemblyFaceChunk - 1)
                     /assemblyFaceChunk
        << " face chunk(s)" << endl;

    const auto tStartTotal = std::chrono::steady_clock::now();
    const auto tStartSetup = tStartTotal;

    const label nStates = ionicModel->nEqns();
    EmbeddedTNNPModel* tnnpModel = &ionicModel();

    // Front-aware hybrid: Gauss-point ODEs in "front" cells, cell-centre ODE +
    // reconstruction elsewhere. It carries the states at the Gauss points (like
    // gaussPointODE), so it is excluded from statesAtCells below.
    const bool useFrontHybrid =
        frontAwareHybrid
     && useHighOrder_Iion
     && mesh.nGeometricD() > 1;

    // stateIntegrationMode = cellCentredReconstruct keeps the persistent states
    // at cell centres (1 ODE/cell) and reconstructs them at the Iion Gauss
    // points on demand; gaussPointODE (and the hybrid) keep them at the Gauss
    // points.
    const bool statesAtCells =
        reconstructStatesFromCellCentres
     && !useFrontHybrid
     && useHighOrder_Iion
     && mesh.nGeometricD() > 1;

    Field<Field<scalar>> states
    (
        statesAtCells ? mesh.nCells() : totalIionIntegrationPoints,
        Field<scalar>(nStates, 0.0)
    );
    ionicModel->initialiseStates(states);

    // Transient IP-sized buffer + scratch field used when reconstructing the
    // cell-centred states at the Iion Gauss points (cellCentredReconstruct).
    Field<Field<scalar>> statesIP
    (
        statesAtCells ? totalIionIntegrationPoints : 0,
        Field<scalar>(nStates, 0.0)
    );

    // Front-aware hybrid working set (per-cell front flag and a reusable
    // cell-centred state buffer for the smooth-cell ODE + reconstruction).
    boolList frontCell(useFrontHybrid ? mesh.nCells() : 0, false);

    // Added for cardiacFoam: the per-cell count of mode switches (front <->
    // smooth) for the end-of-run diagnostic, printed only with the hybrid. A
    // cell the front simply crosses switches twice; more means it went back
    // and forth.
    boolList frontCellLast(useFrontHybrid ? mesh.nCells() : 0, false);
    labelList frontSwitches(useFrontHybrid ? mesh.nCells() : 0, 0);
    bool frontClassified = false;
    Field<Field<scalar>> statesCellBuf
    (
        useFrontHybrid ? mesh.nCells() : 0,
        Field<scalar>(nStates, 0.0)
    );

    // Added for cardiacFoam: the hybrid's CELL-CENTRED states, kept between
    // steps as pure cellCentredReconstruct keeps its own. A smooth cell
    // integrates its ODE from here. It used to start, every step, from the
    // average of its Gauss-point states - but those are a reconstruction of the
    // cell states, and averaging a reconstruction is not the identity (on a
    // square cell it adds (h^2/24) Lap(s)): a small diffusion of the states
    // applied once per STEP, so stronger the smaller dt, which showed as a
    // floor of ~2.2e-3 in the hybrid's dt convergence (E02) that no time scheme
    // got below. Only the front cells keep their states at the Gauss points.
    // For a front cell the entry is the average of its Gauss-point states,
    // which is what it starts from the one time it leaves the band.
    Field<Field<scalar>> hybridCellStates;

    scalarField VmIntegrationPoints(totalIionIntegrationPoints, 0.0);
    scalarField IionIntegrationPoints(totalIionIntegrationPoints, 0.0);

    const label dim = max(mesh.nGeometricD(), label(1));
    const scalar dx = characteristicDx(mesh);
    const scalar CFL = timeIntegrationProperties.lookupOrDefault<scalar>("CFL", 0.1);
    const scalar dtExplicitReference = computeStableDeltaT
    (
        conductivity,
        chi.value(),
        Cm.value(),
        CFL,
        dx,
        dim
    );

    scalar dt = runTime.deltaTValue();
    const TimeSchemeCoeffs timeScheme = timeSchemeCoeffs(implicitScheme);
    const scalar theta = timeScheme.theta;
    const label bdfOrder = timeScheme.bdfOrder;
    const bool historyStateODEDriver = (stateODEDriver == "history");
    const EsdirkTableau esdirkTab = esdirk3Tableau();

    // Invariant rather than a reachable branch today (ported from the MMS
    // solver): ESDIRK3 covers every BDF offered. Kept because adding BDF5 would
    // break it silently - a startup of insufficient order does not fail, it
    // caps the observed order without any error.
    if (bdfOrder > 1 && esdirkTab.order + 1 < bdfOrder)
    {
        FatalErrorInFunction
            << "bdfStartup = " << bdfStartup << " supplies starting values of"
            << " order " << (esdirkTab.order + 1) << ", below the O(dt^"
            << bdfOrder << ") that " << implicitScheme << " requires."
            << exit(FatalError);
    }
    const word massMatrixMode = normalizedMassMatrixType(massMatrixType);
    const word& memoryOptimizationMode = memoryOptimization;

    const bool memoryOptimizationEffective = memoryOptimizationIsEffective
    (
        memoryOptimizationMode,
        useHighOrder_Vm,
        dim,
        mesh,
        memoryOptimizationCellThreshold
    );

    const bool useAdaptiveTripletReserve =
        memoryOptimizationEffective
     && memoryOptimizationAdaptiveTripletReserve;

    const label stiffnessTripletsPerFaceReserve =
            useAdaptiveTripletReserve
        ? (
                stabilisationAlpha > SMALL
              ? memoryOptimizationStabilisedTripletsPerFace
              : memoryOptimizationFluxTripletsPerFace
          )
        : label(800);

    const bool trimHeapAfterLargeSetup =
        memoryOptimizationEffective && memoryOptimizationTrimHeap;

    const bool compactMassAssembly =
        memoryOptimizationEffective && memoryOptimizationCompactMassAssembly;

    Info<< "Running high-order electro activation implicit PETSc solver" << nl
        << "Dimension = " << dim << nl
        << "dx = " << dx << nl
        << "requested dt = " << dt << nl
        << "explicit stable dt reference = " << dtExplicitReference << nl
        << "implicitScheme = " << implicitScheme << nl
        // Recorded even where inert (theta schemes, BDF1): the observed order
        // of a BDF run means nothing without its startup and its ODE driver.
        << "bdfStartup = " << bdfStartup << nl
        << "stateODEDriver = " << stateODEDriver << nl
        << "stateODEGaussPointHistory = " << stateODEGaussPointHistory << nl
        << "massMatrix = " << massMatrixMode << nl
        << "linearSolverBackend = " << linearSolverBackend << nl
        << "PETSc linear KSP = " << petscLinearKspType
        << ", PC = " << petscLinearPcType << nl
        << "useHighOrder_Vm = " << useHighOrder_Vm << nl
        << "useHighOrder_Iion = " << useHighOrder_Iion << nl
        << "stateIntegrationMode = " << stateIntegrationMode << nl
        << "frontAwareHybrid = " << frontAwareHybrid
        << (frontAwareHybrid ? " (hybridFrontVmRange = " : "")
        << (frontAwareHybrid ? Foam::name(hybridFrontVmRange) : word(""))
        << (frontAwareHybrid ? " V, hybridFrontDilation = " : "")
        << (frontAwareHybrid ? Foam::name(hybridFrontDilation) : word(""))
        << (frontAwareHybrid ? ")" : "") << nl
        << "stateReconstructLimiter = " << stateReconstructLimiter << nl
        << "stabilisationAlpha = " << stabilisationAlpha << nl
        << "memoryOptimization = " << memoryOptimization
        << " (effective: "
        << (memoryOptimizationEffective ? "true" : "false") << ")" << nl
        << "stiffnessTripletsPerFaceReserve = "
        << stiffnessTripletsPerFaceReserve << nl
        << "compactMassAssembly = "
        << (compactMassAssembly ? "true" : "false") << nl
        << "Iion integration points = " << totalIionIntegrationPoints << nl
        << "nonlinearMethod = " << nonlinearMethod << nl
        << "JFNK linear backend = " << jfnkLinearSolverBackend << nl
        << "JFNK PETSc KSP = " << jfnkPetscKspType
        << ", PC = " << jfnkPetscPcType << nl
        << "state ODE solver = " << timeIntegrationProperties.lookupOrDefault<word>("solver", "default") << nl
        << "nonlinearRelaxation = " << nonlinearRelaxation << nl
        << "nonlinearStatesTolerance = " << nonlinearStatesTolerance << nl
        << "nonlinearRequireStatesConvergence = "
        << nonlinearRequireStatesConvergence << nl
        << "activationThreshold = " << activationThreshold << " V" << nl
        << "stopAfterPointActivation = " << stopAfterPointActivation << nl
        << "stopActivationPoint = " << stopActivationPoint << nl
        << "stopDelayAfterActivation = " << stopDelayAfterActivation << " s" << nl
        << "implicitNonlinearIterations = " << implicitNonlinearIterations << nl
        << "implicitLinearSolver = " << implicitLinearSolver << nl << endl;

    if (statesAtCells)
    {
        Info<< "LREInterp_states order = " << LREInterp_statesPtr().order()
            << " (states evolved at cell centres, reconstructed at Gauss points)"
            << endl;
    }
    if (useFrontHybrid)
    {
        Info<< "LREInterp_states order = " << LREInterp_statesPtr().order()
            << " (front-aware hybrid: Gauss-point ODEs on the front,"
            << " cell-centre + reconstruction elsewhere)" << endl;

        hybridCellStates.setSize(mesh.nCells(), Field<scalar>(nStates, 0.0));
        averageStatesIntegrationPointsToCells
        (
            states, nStates, LREInterp_IionPtr(), hybridCellStates
        );
    }

    const label stopActivationCellI =
        stopAfterPointActivation
      ? nearestCellToPoint(mesh, stopActivationPoint)
      : -1;

    scalar effectiveEndTime = runTime.endTime().value();
    bool stopPointTriggered = false;
    scalar stopPointActivationTime = -GREAT;

    // Modified for cardiacFoam: nearestCellToPoint returns the cell on exactly
    // ONE rank and -1 on every other, and this banner used to index mesh.C()
    // with it on all of them. Info is a null stream off the master, but its
    // ARGUMENTS are still evaluated, so the losing ranks read mesh.C()[-1]:
    // undefined behaviour, which in an Opt build occasionally decoded as a
    // signalling double and raised SIGFPE in mag() - measured at dx = 0.125 mm
    // on 4 ranks, and it blocked that refinement level.
    //
    // Only the owner reads the mesh now, and its values reach the log whatever
    // rank it is: the others contribute -GREAT to a max-reduction. That also
    // fixes the second defect of the old banner, which printed the master's
    // own -1 whenever the owner was another rank.
    if (stopAfterPointActivation)
    {
        label stopCellGlobal = -1;
        vector stopCellCentre(-GREAT, -GREAT, -GREAT);
        scalar stopCellDistance = -GREAT;

        if (stopActivationCellI >= 0)
        {
            stopCellGlobal = gCellsPtr->toGlobal(stopActivationCellI);
            stopCellCentre = mesh.C()[stopActivationCellI];
            stopCellDistance =
                mag(mesh.C()[stopActivationCellI] - stopActivationPoint);
        }

        reduce(stopCellGlobal, maxOp<label>());
        reduce(stopCellCentre, maxOp<vector>());
        reduce(stopCellDistance, maxOp<scalar>());

        Info<< "Stop point nearest cell = " << stopCellGlobal
            << " (global), centre = " << stopCellCentre
            << ", distance = " << stopCellDistance
            << " m" << nl << endl;
    }

    const wordList exportNames = ionicModel->exportedFieldNames();
    if (!exportNames.empty())
    {
        Info<< "Exporting ionic fields: " << exportNames << nl;
    }

    SpMat M;
    SpMat K;
    SpMat L;

    if (profileTimings)
    {
        logMemoryCheckpoint("before mass assembly");
    }

    if (!useHighOrder_Vm || dim <= 1 || massMatrixMode == "lumped")
    {
        if (massMatrixMode == "consistent" && (!useHighOrder_Vm || dim <= 1))
        {
            Info<< "Consistent mass requested without high-order Vm; "
                << "using the diagonal/lumped mass matrix." << nl;
        }
        assembleDiagonalMassMatrix(mesh, 1.0, M);
    }
    else
    {
        assembleConsistentMassMatrixHO
        (
            mesh,
            LREInterp_VmPtr(),
            1.0,
            compactMassAssembly,
            assemblyCellChunk,
            M
        );
    }

    if (profileTimings)
    {
        logMemoryCheckpoint("after mass assembly");
    }

#ifdef __GLIBC__
    if (trimHeapAfterLargeSetup)
    {
        malloc_trim(0);
    }
#endif

    if (profileTimings)
    {
        logMemoryCheckpoint("before stiffness assembly");
    }

    if (useHighOrder_Vm && dim > 1)
    {
        assembleHighOrderStiffnessMatrix
        (
            mesh,
            conductivity,
            LREInterp_VmPtr(),
            stabilisationAlpha,
            stiffnessTripletsPerFaceReserve,
            assemblyFaceChunk,
            K
        );
    }
    else
    {
        assembleStandardOrthogonalStiffnessMatrix
        (
            mesh,
            conductivity,
            LREInterp_VmPtr(),
            stabilisationAlpha,
            K
        );
    }

    if (profileTimings)
    {
        logMemoryCheckpoint("after stiffness assembly");
    }

    // Added for cardiacFoam, ported from the MMS solver: with CF_DUMP_STENCILS
    // set, write the MLS stencil MEMBERSHIP - one line per cell, "globalCell :
    // global columns" - so that a serial and a parallel build can be diffed
    // cell by cell. This is one level below CF_DUMP_K: it says whether a
    // differing operator row comes from different stencil members or from
    // different coefficients over the same members, which is the question an
    // entry-by-entry operator comparison (E15) cannot answer.
    //
    // Unlike the MMS, this solver builds up to three interpolators from three
    // LRECoeffs blocks, each with its own Nn: Vm (operator and consistent mass
    // matrix), Iion and states (quadrature points). Any of them can lose a
    // member at a processor cut, so every one that exists is dumped. Vm keeps
    // the MMS file names, so the MMS readers work unchanged; the other two carry
    // their tag. The Vm interpolator always exists (see createFields.H); the
    // other two only on the high-order Iion path.
    if (getenv("CF_DUMP_STENCILS"))
    {
        const word rank(Foam::name(Pstream::myProcNo()));

        // Cell centres alongside, at full precision: the selection is a
        // nearest-N by distance, so a centroid that differs in the last bits
        // between a serial and a decomposed mesh is enough to reorder a tie.
        OFstream osc("cellCentre_dump_rank" + rank + ".dat");
        osc.precision(17);
        forAll(mesh.C(), cellI)
        {
            osc << gCellsPtr->toGlobal(cellI) << ' ' << mesh.C()[cellI].x()
                << ' ' << mesh.C()[cellI].y() << ' ' << mesh.C()[cellI].z()
                << nl;
        }

        auto dumpStencils = [&](const word& tag, const highOrderInterp& interp)
        {
            const globalIndex& gc = interp.globalCells();
            const CompactListList<label>& cellSt = interp.globalCellStencils();

            OFstream os("stencil_" + tag + "dump_rank" + rank + ".dat");
            forAll(mesh.C(), cellI)
            {
                os  << gc.toGlobal(cellI) << " :";
                const UList<label> st = cellSt[cellI];
                forAll(st, i)
                {
                    os  << ' ' << st[i];
                }
                os  << nl;
            }

            // Face stencils, keyed by LOCAL face index: the reader maps them
            // with processorN/constant/polyMesh/faceProcAddressing.
            const CompactListList<label>& faceSt = interp.globalFaceStencils();

            OFstream osf("faceStencil_" + tag + "dump_rank" + rank + ".dat");
            forAll(faceSt, faceI)
            {
                osf << faceI << " :";
                const UList<label> st = faceSt[faceI];
                forAll(st, i)
                {
                    osf << ' ' << st[i];
                }
                osf << nl;
            }
        };

        if (LREInterp_VmPtr.valid())
        {
            dumpStencils(word::null, LREInterp_VmPtr());
        }
        if (LREInterp_IionPtr.valid())
        {
            dumpStencils("Iion_", LREInterp_IionPtr());
        }
        if (LREInterp_statesPtr.valid())
        {
            dumpStencils("states_", LREInterp_statesPtr());
        }
    }

    // Added for cardiacFoam, ported from the MMS solver: with CF_DUMP_K set in
    // the environment, write K as (global row, global column, value) triples,
    // one file per rank.
    //
    // This is the check that a row-sum or nnz fingerprint cannot make: mapping
    // the parallel indices back through processorN/constant/polyMesh/
    // cellProcAddressing gives an entry-by-entry comparison against the serial
    // assembly. E15 does that already for the operator handed to KSP through
    // PETSc's option database; this dumps K alone, before M and dt are mixed
    // in, which is where a stiffness-only change (the non-orthogonal
    // correction) has to be measured. Off by default, no cost.
    if (getenv("CF_DUMP_K"))
    {
        OFstream os
        (
            "K_dump_rank" + Foam::name(Pstream::myProcNo()) + ".dat"
        );
        os.precision(17);
        for (label r = 0; r < K.outerSize(); ++r)
        {
            for (SpMat::InnerIterator it(K, r); it; ++it)
            {
                os  << (globalRowStart() + r) << ' ' << it.col() << ' '
                    << it.value() << nl;
            }
        }
    }

#ifdef __GLIBC__
    if (trimHeapAfterLargeSetup)
    {
        malloc_trim(0);
    }
#endif

    L = K;
    L *= 1.0/(chi.value()*Cm.value());

    autoPtr<OFstream> nonlinearResidualFilePtr;
    if (writeNonlinearResiduals)
    {
        // Modified for cardiacFoam: see postProcessingDir().
        const fileName residualDir(postProcessingDir(runTime));
        mkDir(residualDir);
        nonlinearResidualFilePtr.reset
        (
            new OFstream(residualDir/"nonlinearResiduals.dat")
        );
        nonlinearResidualFilePtr()
            << "# time_s step nonlinearMethod stateODESolver iter linearIterations linearError "
            << "coupledResidual Vm_relL2 Iion_relL2 maxState_relL2 "
            << "lineSearchIters converged";
        for (label stateI = 0; stateI < nStates; ++stateI)
        {
            nonlinearResidualFilePtr() << " state" << stateI << "_relL2";
        }
        nonlinearResidualFilePtr() << nl;
    }

    // In-loop timers of the shared reconstructions (monodomainBase.H).
    FineGrainedTimings fineTimings;

    PetscShellKspSolver jfnkPetscShellSolver;

    // The JFNK preconditioner and the controls of the JFNK linear solve,
    // shared with the MMS solver (monodomainKrylov.H).
    JfnkPreconditioner jfnkPc
    (
        jfnkPreconditioner,
        usesPetscBackend(jfnkLinearSolverBackend),
        jfnkPreconditionerKspType,
        jfnkPreconditionerPcType,
        jfnkPreconditionerTolerance,
        jfnkPreconditionerMaxIterations,
        jfnkMaxKrylovIterations,
        jfnkPreconditionerOptionsPrefix,
        petscUseOptions,
        jfnkPreconditionerFillFactor,
        jfnkPreconditionerDropTolerance,
        jfnkPreconditionerUpdateFrequency
    );

    JfnkKrylovControls jfnkKrylov;
    jfnkKrylov.petsc = usesPetscBackend(jfnkLinearSolverBackend);
    jfnkKrylov.kspType = jfnkPetscKspType;
    jfnkKrylov.pcType = jfnkPetscPcType;
    jfnkKrylov.restart = jfnkPetscRestart;
    jfnkKrylov.maxKrylovIterations = jfnkMaxKrylovIterations;
    jfnkKrylov.maxRestarts = jfnkMaxRestarts;
    jfnkKrylov.tolerance = jfnkLinearTolerance;
    jfnkKrylov.optionsPrefix = jfnkPetscOptionsPrefix;
    jfnkKrylov.useOptions = petscUseOptions;

    const auto tEndSetup = std::chrono::steady_clock::now();
    const auto tStartLoop = tEndSetup;

    label nSteps = 0;

    // ----------------------------------------------------------------- //
    // State for the linear-extrapolation Newton initial guess.
    //
    //   Vm_init = Vm_n + (Vm_n - Vm_{n-1}) * (dt / dt_prev)
    //
    // Requires Vm from the previous time step (VmTwoStepsAgo) and the
    // previous time-step length (dtPrev). On the first step we have no
    // history so we fall back to x_0 = Vm_n (legacy behaviour).
    // ----------------------------------------------------------------- //
    scalarField VmTwoStepsAgo(mesh.nCells(), 0.0);
    scalar dtPrev = 0.0;
    bool hasPrevStep = false;

    // Live activation output: write the diagonal/P8 samples once up front (so
    // the file exists from t=0 with every point still pending) and then again
    // whenever a new point crosses the activation threshold, so the run can be
    // monitored as it progresses. writeActivationSamples() overwrites the files.
    label nActivatedPrev = 0;
    forAll(calculateActivationTime, cellI)
    {
        if (!calculateActivationTime[cellI])
        {
            ++nActivatedPrev;
        }
    }
    reduce(nActivatedPrev, sumOp<label>());

    // Activation-failure tracking (see abortOnActivationFailure). The status is
    // written to postProcessing/.../activationStatus.dat at the end.
    scalar tLastActivationChange = runTime.value();
    word activationStatus("RUNNING");
    bool activationAborted = false;

    writeActivationSamples
    (
        runTime,
        activationTime,
        diagonalStart,
        diagonalEnd,
        p8Point,
        nDiagonalSamples,
        sampleNeighbours,
        samplePower,
        false
    );

    // ----------------------------------------------------------------------- //
    // Phase B: persistent linear-solver state, reused across time steps and
    // nonlinear iterations. AImplicit/BImplicit only depend on the (constant)
    // M, L and theta and on dt, which is fixed except on the final clipped
    // step; they are therefore assembled once and rebuilt only when dt changes.
    // The factorisations (PETSc KSP or Eigen) live for the whole run and are
    // refreshed according to the MatrixUpdate policy computed each solve.
    // ----------------------------------------------------------------------- //
    SpMat AImplicit;
    SpMat BImplicit;

    // Added for cardiacFoam: PETSc-backed operators for the two products that
    // used to be plain Eigen expressions. They are rebuilt only when the
    // matrices are (on the first step and whenever dt changes), so the extra
    // PETSc assembly is not on the per-iteration path.
    DistributedMatVec applyAImplicitOp;
    DistributedMatVec applyBImplicitOp;

    auto applyAImplicit = [&](const EigVec& x) { return applyAImplicitOp(x); };
    auto applyBImplicit = [&](const EigVec& x) { return applyBImplicitOp(x); };
    scalar assembledDt = -1.0;
    bool linearSolverNeedsRebuild = true;
    PetscKspMatrixSolver persistentPetscSolver;
    PersistentEigenSolvers persistentEigenSolvers;

    // Added for cardiacFoam, ported from the MMS solver: the BDF history, i.e.
    // the k-1 levels OLDER than V^n that a BDF-k step needs. VmBdfOlder[j] holds
    // V^{n-1-j}, so index 0 is the most recent of them; V^n itself is VmOld.
    // Empty for the theta family and BDF1. Sized and zeroed up front so that a
    // level read before it is filled is a wrong answer rather than an
    // out-of-bounds read - which bdfValidLevels makes impossible anyway.
    List<scalarField> VmBdfOlder(max(bdfOrder - 1, label(0)));
    forAll(VmBdfOlder, j)
    {
        VmBdfOlder[j].setSize(mesh.nCells(), 0.0);
    }

    // How many of those levels hold a value on the CURRENT smooth piece of the
    // solution at the nominal spacing. A BDF step needs all k-1; until then
    // the step is taken by ESDIRK3. Reset to zero by everything that breaks
    // the history: the start of the run, a change in the stimulus (a
    // discontinuity in the source, across which the history polynomial would
    // interpolate a kink), a step of non-nominal length, and a rolled-back step.
    label bdfValidLevels = 0;

    // The stimulus fractions of the previous step, to detect that change.
    List<scalar> stimulusFractionPrev;

    // The step BDF assumes. A step shorter than this (the clipped last step,
    // the stop point) is taken by ESDIRK3.
    const scalar nominalDt = runTime.deltaTValue();

    // The ESDIRK stage operator M/(gamma dt) - L, one for every implicit stage
    // because the method is SINGLY diagonally implicit, with its own cached
    // solvers so that the BDF operator's factorisation is not thrown away. Only
    // built for BDF2-4.
    SpMat AEsdirk;
    DistributedMatVec applyEsdirkOp;
    PetscKspMatrixSolver esdirkPetscSolver;
    PersistentEigenSolvers esdirkEigenSolvers;
    bool esdirkSolverNeedsRebuild = true;

    // L as a distributed product, for the stage right-hand sides F = L Y + s.
    DistributedMatVec applyLOp;
    if (bdfOrder > 1)
    {
        applyLOp.reset(L);
    }

    // Modified for cardiacFoam: stop within a small fraction of the NOMINAL
    // step of the end time, not within SMALL.
    //
    // The accumulated time is a floating-point sum of dt and can land just
    // BELOW the end time by more than SMALL (1e-15): 160 steps of 3.2e-3 sum to
    // 0.51199999999999868, 1.3e-15 short of 0.512. With the old test that
    // remainder took one more pass through the loop, clipped to dt = 1.3e-15:
    // an operator M/dt fifteen orders larger than the real one, a
    // refactorisation, a degenerate solve and an extra step in the count. The
    // MMS solver hit the same accumulation and took a whole extra step.
    //
    // The MMS fix (stop within half a step) does not carry over, because dt is
    // not constant here: the last step is clipped to land on the end time, and
    // the stop-point trigger moves effectiveEndTime to an arbitrary instant. A
    // half-step test would drop a legitimate partial last step. The tolerance
    // is therefore relative to the nominal dt, and small: the accumulation
    // error is ~ nSteps*eps*t, about 1e-11 s for 1e5 steps over 1 s, i.e.
    // ~1e-6 of dt = 1e-5 s, so 1e-3*dt leaves three orders of margin while the
    // most it can drop is a remainder of 1e-3*dt, which no physics resolves.
    // Taken once, before the loop, from the nominal dt: setDeltaT() below
    // changes deltaTValue() on the clipped last step.
    const scalar endTimeTolerance = 1.0e-3*runTime.deltaTValue();

    // Nonlinear iterations that stopped at the residual floor instead of the
    // tolerance (nonlinearIterationExit, monodomainBase.H).
    label nResidualFloorExits = 0;

    // The floor of the denominator of the Iion increment test
    // (relativeL2Difference, monodomainBase.H): with nonlinearIionIncrementFloor
    // on, the largest global L2 norm of Iion over the accepted steps so far;
    // off, zero, which is the plain relative test. The net ionic current
    // crosses zero on the plateau, where ||Iion|| of the instant is no scale at
    // all. A running maximum is uniform across ranks and independent of dt to
    // the accuracy of the solution, but it loosens the test everywhere past
    // the peak of the upstroke, not only at the crossing.
    scalar IionIncrementFloor =
        nonlinearIionIncrementFloor
      ? Foam::sqrt(gSumSqr(Iion.primitiveField()))
      : 0.0;

    // Rolled-back steps: in total, and in a row. A rolled-back step advances
    // time with the solution frozen, so a run of them silently stops the
    // dynamics - and the run still ends with exit 0.
    label nRolledBackSteps = 0;
    label nConsecutiveRollbacks = 0;

    // Modified for cardiacFoam: the time of step n is loopStartTime + n*dt,
    // set after every nominal step, not the running sum of dt. The sum drifts:
    // at dt = 1.25e-5 it reached the end time short by 2e-9*dt, which clipped
    // the last step to a non-nominal length and handed it to ESDIRK3; and a
    // window edge of the stimulus fell a sliver off a step boundary (see
    // applyStimulus).
    const scalar loopStartTime = runTime.value();

    while (runTime.value() < effectiveEndTime - endTimeTolerance)
    {
        const scalar t0 = runTime.value();
        const scalar remaining = effectiveEndTime - t0;

        // A remainder within the clippedStep tolerance of the nominal step is
        // the nominal step.
        dt = remaining < (1.0 - 1.0e-9)*nominalDt ? remaining : nominalDt;
        runTime.setDeltaT(dt);

        // ----------------------------------------------------------------- //
        // Build the theta-scheme implicit operators.
        //
        //   The monodomain PDE in semi-discrete form is
        //
        //       dVm/dt = L * Vm + s(Vm, states)
        //
        //   where L = K / (chi * Cm) is the discrete diffusion operator and
        //   s = -Iion/(chi*Cm) + Iext/(chi*Cm) lumps the reactive sources.
        //   Applying the theta rule (theta=1: BE, theta=0.5: CN) and dividing
        //   by dt yields
        //
        //       AImplicit * Vm^{n+1} = BImplicit * Vm^n + source-blend
        //
        //   with  AImplicit = M/dt - theta*L
        //         BImplicit = M/dt + (1-theta)*L
        //
        //   Under BDF-k (theta = 1) A carries a_0/dt instead of 1/dt, and Vm^n
        //   is replaced by the history vector -sum_{j>=1} a_j V^{n+1-j}; see
        //   TimeSchemeCoeffs. The ESDIRK stage operator M/(gamma dt) - L is
        //   assembled alongside for BDF2-4.
        //
        //   M is either the lumped (diagonal) or consistent high-order mass
        //   matrix; L already carries the 1/(chi*Cm) scaling.
        //
        //   Phase B: AImplicit/BImplicit are constant while dt is unchanged, so
        //   they are only reassembled on the first step and whenever dt changes
        //   (the final clipped step). A dt change also forces the persistent
        //   factorisation to be rebuilt on the next solve.
        // ----------------------------------------------------------------- //
        if (mag(dt - assembledDt) > SMALL)
        {
            // Modified for cardiacFoam: the only thing BDF changes in A is the
            // scalar in front of M. a_0 is 1 for the theta family and BDF1, so
            // backwardEuler and BDF1 build the same matrix bit for bit.
            AImplicit = M;
            AImplicit *= (timeScheme.a[0]/dt);
            AImplicit -= theta*L;

            if (bdfOrder > 1)
            {
                AEsdirk = M;
                AEsdirk *= (1.0/(esdirkTab.gamma*dt));
                AEsdirk -= L;
                applyEsdirkOp.reset(AEsdirk);
                esdirkSolverNeedsRebuild = true;
            }

            BImplicit = M;
            BImplicit *= (1.0/dt);
            if (theta < 1.0 - SMALL)
            {
                // Crank-Nicolson (or any theta < 1) needs the explicit half of
                // the diffusion operator on the right-hand side. For pure
                // Backward Euler (theta=1) this term vanishes and is skipped.
                BImplicit += (1.0 - theta)*L;
            }

            // Added for cardiacFoam: mirror the freshly assembled matrices
            // into PETSc. Must follow every reassembly, which is exactly here.
            reportOperatorFingerprint("AImplicit", AImplicit);
            reportOperatorFingerprint("BImplicit", BImplicit);

            applyAImplicitOp.reset(AImplicit);
            applyBImplicitOp.reset(BImplicit);

            assembledDt = dt;
            linearSolverNeedsRebuild = true;
        }

        volScalarField VmOld
        (
            IOobject
            (
                "VmOld",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Vm
        );

        const scalarField VmOldValues(Vm.primitiveField());
        Field<Field<scalar>> statesOld(states);

        // The hybrid's cell-centred states at t^n, and the ones the iterate
        // being solved ends with. The latter is written only by evaluations of
        // statesGuess (statesGuessPtr), never by the auxiliary ones - the
        // diagonalIion perturbation, the JFNK probes - which share the buffer.
        const Field<Field<scalar>> hybridCellStatesOld(hybridCellStates);
        Field<Field<scalar>> hybridCellStatesGuess(hybridCellStates);
        const Field<Field<scalar>>* statesGuessPtr = nullptr;

        volScalarField IionOld
        (
            IOobject
            (
                "IionOld",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Iion
        );

        volScalarField IionGuess
        (
            IOobject
            (
                "IionGuess",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Iion
        );

        const List<scalar> stimulusFraction = applyStimulus
        (
            t0,
            dt,
            externalStimulusCurrent,
            stimulusCellIDsList,
            stimulusStartTimes,
            stimulusIntensity.value(),
            stimulusDuration.value()
        );

        // Added for cardiacFoam: does BDF take this step, or ESDIRK3?
        //
        // A multistep method interpolates through its history, so any
        // discontinuity in the source between the levels it uses is
        // interpolated as if it were smooth: the k-1 steps after it carry an
        // O(dt) error that caps the global order, with no warning. The stimulus
        // is exactly such a discontinuity. Its per-step value changes only
        // where a window opens or closes, and there the history is DISCARDED:
        // V^n, on the boundary, is a valid node of the new smooth piece, the
        // older levels are not. The fractions depend on t0 and dt alone, so the
        // decision is uniform across ranks.
        //
        // A step of non-nominal length (the clipped last step, the stop point)
        // cannot be a BDF step either: the coefficients assume equal spacing.
        const bool clippedStep = mag(dt - nominalDt) > 1.0e-9*nominalDt;

        if (bdfOrder > 1)
        {
            bool stimulusChanged = false;
            if (stimulusFractionPrev.size() == stimulusFraction.size())
            {
                forAll(stimulusFraction, bI)
                {
                    if
                    (
                        mag(stimulusFraction[bI] - stimulusFractionPrev[bI])
                      > 1.0e-12
                    )
                    {
                        stimulusChanged = true;
                    }
                }
            }
            stimulusFractionPrev = stimulusFraction;

            if (stimulusChanged && bdfValidLevels > 0)
            {
                Info<< "BDF history discarded at t = " << t0
                    << ": the stimulus changed. The next " << (bdfOrder - 1)
                    << " step(s) are taken by ESDIRK3." << endl;
                bdfValidLevels = 0;
            }
        }

        const bool esdirkThisStep =
            bdfOrder > 1 && (bdfValidLevels < bdfOrder - 1 || clippedStep);

        // Added for cardiacFoam: V^n at the Iion integration points, kept for the
        // whole step. It is the first node of the linear Vm driver on the paths
        // that integrate the ODEs at those points (gaussPointODE and the front
        // cells of the hybrid); VmIntegrationPoints itself is overwritten with
        // every candidate.
        scalarField VmOldIntegrationPoints;

        if (useHighOrder_Iion)
        {
            reconstructVmAtIionIntegrationPoints
            (
                VmOld,
                useHighOrder_Vm,
                LREInterp_VmPtr(),
                LREInterp_IionPtr(),
                VmIntegrationPoints,
                stateODEUseOpenMP,
                stateODEOpenMPThreshold,
                profileTimings,
                fineTimings
            );

            VmOldIntegrationPoints = VmIntegrationPoints;

            if (statesAtCells)
            {
                // statesOld lives at cell centres -> reconstruct to the Gauss
                // points before evaluating Iion there.
                reconstructStatesAtIionIntegrationPoints
                (
                    mesh,
                    statesOld,
                    nStates,
                    LREInterp_statesPtr(),
                    LREInterp_IionPtr(),
                    statesIP,
                    nullptr,
                    stateReconstructLimiter,
                    stateODEUseOpenMP,
                    stateODEOpenMPThreshold
                );
            }

            ionicModel->calculateCurrent
            (
                t0,
                dt,
                VmIntegrationPoints,
                IionIntegrationPoints,
                statesAtCells ? statesIP : statesOld
            );

            averageIionIntegrationPointsToCells
            (
                IionIntegrationPoints,
                LREInterp_IionPtr(),
                IionOld
            );
        }
        else
        {
            ionicModel->calculateCurrent
            (
                t0,
                dt,
                VmOld.internalField(),
                IionOld,
                statesOld
            );
            IionOld.correctBoundaryConditions();
        }

        // Front-aware hybrid: classify the front cells ONCE per time step (from
        // the committed Vm, reconstructed above into VmIntegrationPoints) and
        // freeze the flag for the whole nonlinear solve. Dilating a couple of
        // cell layers covers the front's motion within the step and keeps the
        // fixed-point map continuous (no front/smooth flip between iterations).
        if (useFrontHybrid)
        {
            flagFrontCells
            (
                mesh,
                VmIntegrationPoints,
                LREInterp_IionPtr(),
                externalStimulusCurrent,
                hybridFrontVmRange,
                frontCell
            );
            dilateFrontCells(mesh, hybridFrontDilation, frontCell);

            if (frontClassified)
            {
                forAll(frontCell, cellI)
                {
                    if (frontCell[cellI] != frontCellLast[cellI])
                    {
                        ++frontSwitches[cellI];
                    }
                }
            }
            frontCellLast = frontCell;
            frontClassified = true;
        }

        scalarField sourceOld(mesh.nCells(), 0.0);
        scalarField sourceGuess(mesh.nCells(), 0.0);
        forAll(sourceOld, cellI)
        {
            sourceOld[cellI] =
               -IionOld[cellI]
              + externalStimulusCurrent[cellI]/(chi.value()*Cm.value());
        }

        // Modified for cardiacFoam, ported from the MMS solver: the history
        // vector handed to B. For the theta family it is V^n; for BDF-k it is
        // -sum_{j>=1} a_j V^{n+1-j}, formed HERE as one field so that a BDF step
        // still costs a single mass-matrix product. For BDF1, a_1 = -1 and
        // there are no older levels, so the scaling is by exactly 1.0 and the
        // result is bit-identical to backwardEuler. An ESDIRK step does not use
        // it (it forms B V^n itself).
        EigVec Vn = fieldToEigVec(VmOld);

        if (bdfOrder > 0 && !esdirkThisStep)
        {
            Vn *= -timeScheme.a[1];

            forAll(VmBdfOlder, j)
            {
                const scalar aj = timeScheme.a[j + 2];
                const scalarField& Vj = VmBdfOlder[j];

                forAll(Vj, cellI)
                {
                    Vn[cellI] -= aj*Vj[cellI];
                }
            }
        }

        volScalarField VmGuess
        (
            IOobject
            (
                "VmGuess",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            VmOld
        );

        // ------------------------------------------------------------- //
        // Newton initial guess: linear extrapolation Vm_n + (Vm_n -
        // Vm_{n-1}) * (dt / dt_prev), clamped to a physical bracket.
        // First time step has no history -> identity (Vm_n).
        // ------------------------------------------------------------- //
        if (jfnkInitGuessOrder >= 1 && hasPrevStep && dtPrev > SMALL)
        {
            scalarField& Vg = VmGuess.primitiveFieldRef();
            const scalarField& Vn_ = VmOld.primitiveField();
            const scalar ratio = dt / dtPrev;
            forAll(Vg, cellI)
            {
                scalar v = Vn_[cellI] + ratio*(Vn_[cellI] - VmTwoStepsAgo[cellI]);
                if (v < jfnkInitGuessVmMin) v = jfnkInitGuessVmMin;
                if (v > jfnkInitGuessVmMax) v = jfnkInitGuessVmMax;
                Vg[cellI] = v;
            }
            VmGuess.correctBoundaryConditions();
        }

        // Snapshot Vm_n into VmTwoStepsAgo *before* Newton mutates the
        // field, so the next iteration of the time loop sees the correct
        // "previous step" Vm even after rollback.
        VmTwoStepsAgo = VmOld.primitiveField();
        dtPrev = dt;
        hasPrevStep = true;

        if (tnnpModel)
        {
            tnnpModel->updateStatesOld(statesOld);
        }

        volScalarField sourceDerivative
        (
            IOobject
            (
                "sourceDerivative",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("sourceDerivative", dimless/dimTime, 0.0),
            Vm.boundaryField().types()
        );

        // Added for cardiacFoam, ported from the MMS solver: the interval the
        // state ODEs are integrated over, measured from t0. For every one-step
        // scheme and for BDF it is the whole step. ESDIRK is the exception:
        // stage i lands at t0 + c_i dt, so its states are integrated over
        // c_i dt. Only the ESDIRK block assigns it, and it restores dt before
        // leaving - a stale value would integrate the states over the wrong
        // interval, which a convergence study reads as a bad order.
        scalar stageDt = dt;

        // Added for cardiacFoam, ported from the MMS solver: the Vm driver of
        // the state ODEs on the CELL-indexed paths (cellCentredReconstruct and
        // the low-order Iion path). Node 0 is V^n; the last node is the
        // candidate, which solveODE takes from its Vm argument. The two-node
        // linear blend set here is what backwardEuler, crankNicolson and BDF1
        // use; BDF2-4 and the ESDIRK stages raise nNodes further down unless
        // stateODEDriver is linear.
        //
        // The Gauss-point paths (gaussPointODE and the front cells of the
        // hybrid) keep the linear blend always, as the MMS does: the
        // higher-order interpolant there would mean storing every history level
        // and every stage value at every quadrature point.
        VmDriverStages cellDriver;
        cellDriver.nNodes = 2;
        cellDriver.s[0] = 0.0;
        cellDriver.s[1] = 1.0;
        cellDriver.field[0] = &VmOldValues;

        // A BDF step drives its cell-indexed state ODEs with the polynomial
        // through its OWN history levels (ported from the MMS solver): nodes
        // at s = -(k-1) ... -1, 0, 1 in units of dt from t0, the candidate
        // last. The ODEs are integrated over s in [0, 1], inside the node hull,
        // so this interpolates rather than extrapolates. The levels are already
        // stored, so it costs only the interpolation. An ESDIRK step sets its
        // own driver below.
        if (bdfOrder > 1 && !esdirkThisStep && historyStateODEDriver)
        {
            cellDriver.nNodes = bdfOrder + 1;

            // Oldest first, so the candidate lands on the last slot.
            for (label j = 0; j < bdfOrder - 1; ++j)
            {
                // VmBdfOlder[m] holds V^{n-1-m}, so walk it backwards.
                const label m = bdfOrder - 2 - j;
                cellDriver.s[j] = -scalar(m + 1);
                cellDriver.field[j] = &VmBdfOlder[m];
            }

            cellDriver.s[bdfOrder - 1] = 0.0;
            cellDriver.field[bdfOrder - 1] = &VmOldValues;

            cellDriver.s[bdfOrder] = 1.0;
            cellDriver.field[bdfOrder] = nullptr;
        }

        // The linear blend at the Gauss points, for the paths above.
        VmDriverStages gaussDriver;
        gaussDriver.nNodes = 2;
        gaussDriver.s[0] = 0.0;
        gaussDriver.s[1] = 1.0;
        gaussDriver.field[0] = &VmOldIntegrationPoints;

        // Added for cardiacFoam: with stateODEGaussPointHistory the Gauss-point
        // paths take the same history polynomial as cellDriver. The BDF levels
        // at the Gauss points are reconstructed from the stored cell levels by
        // the operator that gave VmOldIntegrationPoints, so each is what that
        // step had there and nothing is kept between steps. The ESDIRK block
        // does the same with its stage values.
        const bool gaussHistory =
            stateODEGaussPointHistory && historyStateODEDriver
         && useHighOrder_Iion && (useFrontHybrid || !statesAtCells);
        autoPtr<volScalarField> VmLevelScratchPtr;
        List<scalarField> VmBdfOlderIntegrationPoints;
        if (gaussHistory)
        {
            VmLevelScratchPtr.reset(new volScalarField("VmLevelScratch", VmOld));
        }
        if (gaussHistory && bdfOrder > 1 && !esdirkThisStep)
        {
            VmBdfOlderIntegrationPoints.setSize(bdfOrder - 1);
            forAll(VmBdfOlderIntegrationPoints, m)
            {
                VmBdfOlderIntegrationPoints[m].setSize
                (
                    VmOldIntegrationPoints.size()
                );
                reconstructVmLevelAtIionIntegrationPoints
                (
                    VmBdfOlder[m],
                    VmLevelScratchPtr(),
                    useHighOrder_Vm,
                    LREInterp_VmPtr(),
                    LREInterp_IionPtr(),
                    VmBdfOlderIntegrationPoints[m],
                    stateODEUseOpenMP,
                    stateODEOpenMPThreshold,
                    profileTimings,
                    fineTimings
                );
            }

            gaussDriver.nNodes = bdfOrder + 1;
            for (label j = 0; j < bdfOrder - 1; ++j)
            {
                const label m = bdfOrder - 2 - j;
                gaussDriver.s[j] = -scalar(m + 1);
                gaussDriver.field[j] = &VmBdfOlderIntegrationPoints[m];
            }
            gaussDriver.s[bdfOrder - 1] = 0.0;
            gaussDriver.field[bdfOrder - 1] = &VmOldIntegrationPoints;
            gaussDriver.s[bdfOrder] = 1.0;
            gaussDriver.field[bdfOrder] = nullptr;
        }

        // Evaluate the coupled (Vm, states, Iion) triple for a given
        // membrane-potential candidate.
        //
        // When resetStatesToTimeStart is true the ionic state variables are
        // restored to their values at the beginning of the current time step
        // BEFORE integrating the ODE. This is required during nonlinear
        // iterations: each Newton/Picard step must integrate from t_n to
        // t_{n+1} starting from the same state, otherwise the iterates
        // drift along the ODE trajectory instead of converging.
        auto evaluateNonlinearFields =
        [&](
            const volScalarField& VmCandidate,
            Field<Field<scalar>>& statesCandidate,
            volScalarField& IionCandidate,
            const bool resetStatesToTimeStart
        )
        {
            if (resetStatesToTimeStart)
            {
                // Generic reset: works for every ionic model since `states`
                // is the externally-visible state vector field.
                statesCandidate = statesOld;

                if (tnnpModel)
                {
                    // TNNP additionally keeps its own internal state snapshot
                    // (STATES_OLD_) used for stiff-gating protection; refresh
                    // it here so the next ODE step starts from a consistent
                    // pair of (external, internal) states.
                    tnnpModel->resetStatesToStatesOld(statesCandidate);
                }
            }

            if (useHighOrder_Iion)
            {
                reconstructVmAtIionIntegrationPoints
                (
                    VmCandidate,
                    useHighOrder_Vm,
                    LREInterp_VmPtr(),
                    LREInterp_IionPtr(),
                    VmIntegrationPoints,
                    stateODEUseOpenMP,
                    stateODEOpenMPThreshold,
                    profileTimings,
                    fineTimings
                );

                // Protect TNNP from unphysical probes coming from the
                // matrix-free JFNK Jacobian. VmIntegrationPoints is a
                // working buffer for the Iion evaluation so we clamp it in
                // place; this does not affect VmCandidate or any persistent
                // field.
                if (jfnkClampODEInput)
                {
                    clampPhysical
                    (
                        VmIntegrationPoints,
                        jfnkInitGuessVmMin,
                        jfnkInitGuessVmMax
                    );
                }

                if (useFrontHybrid)
                {
                    // ------------------------------------------------------- //
                    // Front-aware hybrid. statesCandidate is Gauss-sized.
                    //   front cells  : one ODE per Gauss point (local Vm),
                    //   smooth cells : one ODE at the cell centre + LRE
                    //                  reconstruction of the states to Gauss.
                    // VmIntegrationPoints is already (optionally) clamped above.
                    // ------------------------------------------------------- //
                    const CompactListList<point>& cellQP =
                        LREInterp_IionPtr().cellQuadPoints();

                    // (a) frontCell is FROZEN for the whole time step (computed
                    //     once before the nonlinear loop from the committed Vm,
                    //     then dilated). It must NOT be recomputed here: doing so
                    //     per iteration flips the front/smooth switch and makes
                    //     the fixed point a non-convergent period-2 limit cycle.

                    // (b) statesCellBuf <- the cell-centred states at t^n:
                    //     a smooth cell starts its ODE from its OWN state, as
                    //     in cellCentredReconstruct (see hybridCellStates).
                    //     Continuing an integration without the reset keeps
                    //     the old average.
                    if (resetStatesToTimeStart)
                    {
                        statesCellBuf = hybridCellStatesOld;
                    }
                    else
                    {
                        averageStatesIntegrationPointsToCells
                        (
                            statesCandidate, nStates, LREInterp_IionPtr(),
                            statesCellBuf
                        );
                    }

                    // (c) Front cells: gather their Gauss points, advance one
                    //     ODE per point with the local Vm, scatter back.
                    DynamicList<label> frontIPs;
                    {
                        label ip = 0;
                        forAll(mesh.cells(), cellI)
                        {
                            const label nq = cellQP[cellI].size();
                            if (frontCell[cellI])
                            {
                                for (label k = 0; k < nq; ++k)
                                {
                                    frontIPs.append(ip + k);
                                }
                            }
                            ip += nq;
                        }
                    }
                    if (frontIPs.size() > 0)
                    {
                        const label nF = frontIPs.size();
                        scalarField VmF(nF), ImF(nF, 0.0);
                        Field<Field<scalar>> sF(nF);
                        forAll(frontIPs, i)
                        {
                            VmF[i] = VmIntegrationPoints[frontIPs[i]];
                            sF[i] = statesCandidate[frontIPs[i]];
                        }
                        // Modified for cardiacFoam: frontIPs are the true Gauss
                        // point indices of these front points, so passing them
                        // as the scratch map puts each point's algebraic_,
                        // and rates_ in the slot that belongs to it.
                        // Previously this pass used slots 0..nF-1, which the
                        // smooth-cell pass below then overwrote.
                        // Linear Vm driver at these Gauss points, compacted
                        // like VmF (see gaussDriver).
                        // Every stored node of gaussDriver (one with the
                        // linear blend, the whole history polynomial under
                        // stateODEGaussPointHistory), compacted like VmF; the
                        // last node, the candidate, comes from VmF itself.
                        VmDriverStages frontDriver(gaussDriver);
                        List<scalarField> frontNodes(gaussDriver.nNodes);
                        for (label j = 0; j < gaussDriver.nNodes; ++j)
                        {
                            if (gaussDriver.field[j])
                            {
                                const scalarField& full = *gaussDriver.field[j];
                                frontNodes[j].setSize(nF);
                                forAll(frontIPs, i)
                                {
                                    frontNodes[j][i] = full[frontIPs[i]];
                                }
                                frontDriver.field[j] = &frontNodes[j];
                            }
                        }

                        ionicModel->solveODE
                        (
                            t0, stageDt, VmF, ImF, sF, frontIPs, &frontDriver
                        );
                        forAll(frontIPs, i)
                        {
                            statesCandidate[frontIPs[i]] = sF[i];
                        }
                    }

                    // (d) Smooth cells: one ODE at the cell-centre Vm, from the
                    //     cell-average old state; result overwrites statesCellBuf.
                    boolList smoothCell(frontCell.size());
                    forAll(frontCell, cellI)
                    {
                        smoothCell[cellI] = !frontCell[cellI];
                    }
                    DynamicList<label> smoothCells;
                    forAll(smoothCell, cellI)
                    {
                        if (smoothCell[cellI]) smoothCells.append(cellI);
                    }
                    if (smoothCells.size() > 0)
                    {
                        const label nS = smoothCells.size();
                        const scalarField& VmC = VmCandidate.internalField();
                        scalarField VmS(nS), ImS(nS, 0.0);
                        Field<Field<scalar>> sS(nS);
                        forAll(smoothCells, i)
                        {
                            scalar v = VmC[smoothCells[i]];
                            if (jfnkClampODEInput)
                            {
                                if (v < jfnkInitGuessVmMin) v = jfnkInitGuessVmMin;
                                else if (v > jfnkInitGuessVmMax) v = jfnkInitGuessVmMax;
                            }
                            VmS[i] = v;
                            sS[i] = statesCellBuf[smoothCells[i]];
                        }

                        // Added for cardiacFoam: these are cell-centred points,
                        // so they take slots from the second scratch bank. The
                        // map is injective and disjoint from the front pass's
                        // Gauss slots, so the two passes no longer alias.
                        labelList smoothSlots(nS);
                        forAll(smoothCells, i)
                        {
                            smoothSlots[i] =
                                ionicModel->cellSlotOffset() + smoothCells[i];
                        }

                        // Linear Vm driver at these cell centres, compacted
                        // like VmS.
                        scalarField VmOldS(nS);
                        forAll(smoothCells, i)
                        {
                            VmOldS[i] = VmOldValues[smoothCells[i]];
                        }
                        VmDriverStages smoothDriver;
                        smoothDriver.field[0] = &VmOldS;

                        ionicModel->solveODE
                        (
                            t0, stageDt, VmS, ImS, sS, smoothSlots,
                            &smoothDriver
                        );
                        forAll(smoothCells, i)
                        {
                            statesCellBuf[smoothCells[i]] = sS[i];
                        }
                    }

                    // (e) Front cells' new cell state = average of their advanced
                    //     Gauss states (needed as a smooth reconstruction stencil).
                    averageStatesIntegrationPointsToCells
                    (
                        statesCandidate, nStates, LREInterp_IionPtr(),
                        statesCellBuf, &frontCell
                    );

                    // (f) Reconstruct the new cell states to the Gauss points in
                    //     smooth cells only (front Gauss states already advanced).
                    reconstructStatesAtIionIntegrationPoints
                    (
                        mesh,
                        statesCellBuf,
                        nStates,
                        LREInterp_statesPtr(),
                        LREInterp_IionPtr(),
                        statesCandidate,
                        &smoothCell,
                        stateReconstructLimiter,
                        stateODEUseOpenMP,
                        stateODEOpenMPThreshold
                    );

                    // (g) Iion at all Gauss points from the assembled states.
                    ionicModel->calculateCurrent
                    (
                        t0,
                        dt,
                        VmIntegrationPoints,
                        IionIntegrationPoints,
                        statesCandidate
                    );

                    // The cell-centred states this evaluation ends with: the
                    // integrated ones in smooth cells, the average of the
                    // advanced Gauss states in front cells. Kept only for the
                    // iterate itself.
                    if (&statesCandidate == statesGuessPtr)
                    {
                        hybridCellStatesGuess = statesCellBuf;
                    }
                }
                else if (statesAtCells)
                {
                    // cellCentredReconstruct: one ODE per cell (driven by the
                    // cell-centred Vm), then reconstruct the states at the
                    // Gauss points and evaluate Iion there (compute-only).
                    if (jfnkClampODEInput)
                    {
                        scalarField VmClamped(VmCandidate.internalField());
                        clampPhysical
                        (
                            VmClamped,
                            jfnkInitGuessVmMin,
                            jfnkInitGuessVmMax
                        );
                        ionicModel->solveODE
                        (
                            t0, stageDt, VmClamped, IionCandidate,
                            statesCandidate,
                            EmbeddedTNNPModel::emptyScratchIDs(), &cellDriver
                        );
                    }
                    else
                    {
                        ionicModel->solveODE
                        (
                            t0,
                            stageDt,
                            VmCandidate.internalField(),
                            IionCandidate,
                            statesCandidate,
                            EmbeddedTNNPModel::emptyScratchIDs(),
                            &cellDriver
                        );
                    }

                    reconstructStatesAtIionIntegrationPoints
                    (
                        mesh,
                        statesCandidate,
                        nStates,
                        LREInterp_statesPtr(),
                        LREInterp_IionPtr(),
                        statesIP,
                        nullptr,
                        stateReconstructLimiter,
                        stateODEUseOpenMP,
                        stateODEOpenMPThreshold
                    );

                    ionicModel->calculateCurrent
                    (
                        t0,
                        dt,
                        VmIntegrationPoints,
                        IionIntegrationPoints,
                        statesIP
                    );
                }
                else
                {
                    // gaussPointODE: integrate the ODE independently at every
                    // Gauss point.
                    ionicModel->solveODE
                    (
                        t0,
                        stageDt,
                        VmIntegrationPoints,
                        IionIntegrationPoints,
                        statesCandidate,
                        EmbeddedTNNPModel::emptyScratchIDs(),
                        &gaussDriver
                    );
                }

                averageIionIntegrationPointsToCells
                (
                    IionIntegrationPoints,
                    LREInterp_IionPtr(),
                    IionCandidate
                );
            }
            else
            {
                if (jfnkClampODEInput)
                {
                    // VmCandidate is const, so we copy into a clamped
                    // scratch buffer before invoking the ODE. The copy is
                    // O(nCells) and negligible compared to the integrator
                    // sub-steps inside solveODE.
                    scalarField VmClamped(VmCandidate.internalField());
                    clampPhysical
                    (
                        VmClamped,
                        jfnkInitGuessVmMin,
                        jfnkInitGuessVmMax
                    );
                    ionicModel->solveODE
                    (
                        t0,
                        stageDt,
                        VmClamped,
                        IionCandidate,
                        statesCandidate,
                        EmbeddedTNNPModel::emptyScratchIDs(),
                        &cellDriver
                    );
                }
                else
                {
                    ionicModel->solveODE
                    (
                        t0,
                        stageDt,
                        VmCandidate.internalField(),
                        IionCandidate,
                        statesCandidate,
                        EmbeddedTNNPModel::emptyScratchIDs(),
                        &cellDriver
                    );
                }
                IionCandidate.correctBoundaryConditions();
            }
        };

        // Build the right-hand side of the implicit system for a given Iion
        // candidate.
        //
        //   rhs = BImplicit * Vm^n  +  theta*src^{n+1}  +  (1-theta)*src^n
        //
        // No explicit dt factor is needed on the source terms: after writing
        // (Vm^{n+1} - Vm^n)/dt = ... and dividing both sides by dt, the
        // sources naturally appear with weights theta and (1-theta) only.
        // The reaction source is  src = -Iion/(chi*Cm) + Iext/(chi*Cm).
        auto rhsFromIion = [&](const volScalarField& IionCandidate)
        {
            // Modified for cardiacFoam: see applyAImplicit.
            EigVec rhs = applyBImplicit(Vn);
            forAll(sourceGuess, cellI)
            {
                const scalar sourceNp1 =
                   -IionCandidate[cellI]
                  + externalStimulusCurrent[cellI]/(chi.value()*Cm.value());

                rhs[cellI] +=
                    theta*sourceNp1
                  + (1.0 - theta)*sourceOld[cellI];
            }
            return rhs;
        };

        auto computeSourceDerivative =
        [&](
            const volScalarField& VmLinearisationPoint,
            const volScalarField& IionLinearisationPoint
        )
        {
            volScalarField VmPlus
            (
                IOobject
                (
                    "VmPlus",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                VmLinearisationPoint
            );

            VmPlus.primitiveFieldRef() += diagonalIionEpsilon;
            VmPlus.correctBoundaryConditions();

            Field<Field<scalar>> statesPlus(statesOld);
            volScalarField IionPlus
            (
                IOobject
                (
                    "IionPlus",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                IionLinearisationPoint
            );

            evaluateNonlinearFields(VmPlus, statesPlus, IionPlus, true);

            scalarField& deriv = sourceDerivative.primitiveFieldRef();
            forAll(deriv, cellI)
            {
                deriv[cellI] =
                  -(IionPlus[cellI] - IionLinearisationPoint[cellI])
                   /max(diagonalIionEpsilon, SMALL);
            }
            sourceDerivative.correctBoundaryConditions();
        };

        // Modified for cardiacFoam: nonlinearMethod is validated at start-up in
        // createFields.H (requireOneOf), one spelling per method, so exactly
        // one of these is true. The aliases picard, jfnk, diagonal and
        // localDiagonal are gone.
        //
        // The validation exists because the failure it prevents is not a bad
        // run but a bad STUDY: a typo in a sweep script used to collapse the
        // whole method matrix onto Picard while the directories were still
        // named after JFNK and diagonalIion and every run exited 0. Measured,
        // `nonlinearMethod Picrad` produced a trace bit-identical to Picard.
        const bool useDiagonalIion = nonlinearMethod == "diagonalIion";
        const bool useJFNK = nonlinearMethod == "JFNK";

        Field<Field<scalar>> statesGuess(statesOld);
        statesGuessPtr = &statesGuess;
        evaluateNonlinearFields(VmGuess, statesGuess, IionGuess, true);

        scalarField stateResiduals(nStates, 0.0);

        auto writeAndPrintResiduals =
        [&](
            const label corr,
            const label linearIterations,
            const scalar linearError,
            const scalar coupledResidual,
            const scalar VmResidual,
            const scalar IionResidual,
            const scalar maxStateResidual,
            const scalarField& stateResidualValues,
            const bool converged,
            const label lineSearchIters = -1  // -1 = N/A (Picard / early exit)
        )
        {
            if (logNonlinearIterations)
            {
                Info<< "Time = " << t0 << "," << nl
                    << "       Nonlinear method:          " << nonlinearMethod << nl
                    << "       Implicit electro solver:   " << implicitLinearSolver
                    << "; iterations = " << linearIterations
                    << ", estimated error = " << linearError << nl
                    << "       NonLinSolver:              iter = " << corr
                    << ", Vm residual = " << VmResidual
                    << ", Iion residual = " << IionResidual
                    << ", max state residual = " << maxStateResidual
                    << ", coupled residual = " << coupledResidual;
                if (lineSearchIters >= 0)
                {
                    Info<< ", lineSearchIters = " << lineSearchIters;
                }
                Info<< ", converged = " << converged << nl
                    << "       State residuals:";

                forAll(stateResidualValues, stateI)
                {
                    Info<< " s" << stateI << '=' << stateResidualValues[stateI];
                }
                Info<< nl;
            }

            if (writeNonlinearResiduals && nonlinearResidualFilePtr.valid())
            {
                nonlinearResidualFilePtr()
                    << t0 << ' ' << (nSteps + 1) << ' ' << nonlinearMethod << ' '
                    << timeIntegrationProperties.lookupOrDefault<word>("solver", "default") << ' '
                    << corr << ' ' << linearIterations << ' ' << linearError << ' '
                    << coupledResidual << ' '
                    << VmResidual << ' ' << IionResidual << ' '
                    << maxStateResidual << ' ' << lineSearchIters << ' ' << converged;
                forAll(stateResidualValues, stateI)
                {
                    nonlinearResidualFilePtr() << ' ' << stateResidualValues[stateI];
                }
                nonlinearResidualFilePtr() << nl;
            }
        };

        // Track whether the nonlinear solve converged in either branch
        // (JFNK or Picard/diagonal). Used after both loops to decide whether
        // to commit the iterate or roll back to the time-step-start state.
        bool nonlinearConverged = false;
        label nonlinearIters = 0;

        // Added for cardiacFoam, ported from the MMS solver: the ESDIRK3 step,
        // for the steps BDF-k cannot take (see esdirkThisStep above).
        //
        // Stage 1 is explicit, Y_1 = V^n. For i >= 2 the stage solves
        //
        //   A_g Y_i = (1/g) [ B V^n + sum_{j<i} a_ij F_j ] + s(Y_i)
        //
        // with A_g = M/(g dt) - L, B = M/dt (theta = 1 under BDF, so BImplicit
        // is exactly that), F_j = L Y_j + s(Y_j) and s the reaction source
        // -Iion + Iext/(chi Cm). That falls out of dividing the stage equation
        // M (Y_i - V^n)/dt = sum_{j<=i} a_ij F_j by g. Stiffly accurate: the
        // last stage IS V^{n+1}, with its states, so no combination step.
        //
        // Every stage is solved by a Picard fixed point, whatever
        // nonlinearMethod says, as on the MMS side: JFNK and diagonalIion would
        // each need their per-stage Jacobian machinery threaded through, and
        // these are at most k-1 steps per restart. Relaxation and stopping
        // criteria are the Picard branch's own.
        if (esdirkThisStep)
        {
            const label nS = esdirkTab.nStages;
            const scalar gam = esdirkTab.gamma;
            const scalar chiCm = chi.value()*Cm.value();

            // The reaction source s(Y) at the cells, for the stage equations.
            auto stageSource = [&](const volScalarField& IionCandidate)
            {
                EigVec src(mesh.nCells());
                forAll(IionCandidate, cellI)
                {
                    src[cellI] =
                       -IionCandidate[cellI]
                      + externalStimulusCurrent[cellI]/chiCm;
                }
                return src;
            };

            List<EigVec> stageF(nS);
            const EigVec BVn = applyBImplicit(fieldToEigVec(VmOld));

            // The stage VALUES, kept so the state ODEs of later stages can be
            // driven by a polynomial through them. Sized up front: a zero-length
            // field behind a live driver pointer would be read out of bounds.
            List<scalarField> stageVm(nS);
            forAll(stageVm, j)
            {
                stageVm[j].setSize(mesh.nCells(), 0.0);
            }
            stageVm[0] = VmOldValues;

            // The stage values at the Gauss points, for the Gauss-point paths
            // under stateODEGaussPointHistory; each is reconstructed once, when
            // its stage is accepted.
            List<scalarField> stageVmIntegrationPoints(gaussHistory ? nS : 0);
            if (gaussHistory)
            {
                stageVmIntegrationPoints[0] = VmOldIntegrationPoints;
            }

            Info<< "ESDIRK3 step at t = " << t0 << " (dt = " << dt << "): "
                << (clippedStep ? "non-nominal step length" : "building the BDF history")
                << ", stages solved by Picard" << endl;

            // Stage 1, explicit: Y_1 = V^n with the states at t^n, which is what
            // stageDt = 0 delivers (the ODE integrator returns over a zero
            // interval).
            stageDt = 0.0;
            VmGuess.primitiveFieldRef() = VmOldValues;
            VmGuess.correctBoundaryConditions();
            evaluateNonlinearFields(VmGuess, statesGuess, IionGuess, true);

            stageF[0] =
                applyLOp(fieldToEigVec(VmGuess)) + stageSource(IionGuess);

            nonlinearConverged = true;

            for (label i = 1; i < nS; ++i)
            {
                stageDt = esdirkTab.c[i]*dt;

                // Drive the cell-indexed state ODEs with the polynomial through
                // the stage values already known, in time normalised by THIS
                // stage's interval: s_j = c_j/c_i, so s_0 = 0 and s_i = 1. The
                // last node is the stage being solved. The abscissae need not be
                // monotone; [0, 1] lies inside their hull because c_0 = 0 and
                // c_i are both nodes.
                if (historyStateODEDriver)
                {
                    cellDriver.nNodes = i + 1;
                    for (label j = 0; j < i; ++j)
                    {
                        cellDriver.s[j] = esdirkTab.c[j]/esdirkTab.c[i];
                        cellDriver.field[j] = &stageVm[j];
                    }
                    cellDriver.s[i] = 1.0;
                    cellDriver.field[i] = nullptr;

                    if (gaussHistory)
                    {
                        gaussDriver.nNodes = i + 1;
                        for (label j = 0; j < i; ++j)
                        {
                            gaussDriver.s[j] = esdirkTab.c[j]/esdirkTab.c[i];
                            gaussDriver.field[j] = &stageVmIntegrationPoints[j];
                        }
                        gaussDriver.s[i] = 1.0;
                        gaussDriver.field[i] = nullptr;
                    }
                }
                else
                {
                    // The ramp from stage 0 (= V^n) to the stage being solved.
                    cellDriver.nNodes = 2;
                    cellDriver.s[0] = 0.0;
                    cellDriver.s[1] = 1.0;
                    cellDriver.field[0] = &stageVm[0];
                    cellDriver.field[1] = nullptr;
                }

                EigVec rhsBase = BVn;
                for (label j = 0; j < i; ++j)
                {
                    rhsBase += esdirkTab.a[i][j]*stageF[j];
                }
                rhsBase /= gam;

                // Fields consistent with the current stage guess and interval.
                evaluateNonlinearFields(VmGuess, statesGuess, IionGuess, true);

                bool stageConverged = false;

                scalar coupledResidualPrevIter = GREAT;

                for
                (
                    label corr = 0;
                    corr < max(implicitNonlinearIterations, label(1));
                    ++corr
                )
                {
                    const scalarField VmPrevious(VmGuess.primitiveField());
                    const scalarField IionPrevious(IionGuess.primitiveField());
                    Field<Field<scalar>> statesPrevious(statesGuess);

                    const EigVec rhs = rhsBase + stageSource(IionGuess);

                    label linearIterations = 0;
                    scalar linearError = GREAT;

                    MatrixUpdate updatePolicy = MatrixUpdate::Reuse;
                    if (esdirkSolverNeedsRebuild)
                    {
                        updatePolicy = MatrixUpdate::Rebuild;
                        esdirkSolverNeedsRebuild = false;
                    }

                    const EigVec Ysol = solveSparseSystem
                    (
                        AEsdirk,
                        rhs,
                        linearSolverBackend,
                        implicitLinearSolver,
                        petscLinearKspType,
                        petscLinearPcType,
                        petscLinearRestart,
                        petscLinearOptionsPrefix,
                        petscUseOptions,
                        implicitTolerance,
                        implicitMaxIterations,
                        linearIterations,
                        linearError,
                        esdirkPetscSolver,
                        esdirkEigenSolvers,
                        updatePolicy,
                        petscReusePreconditioner
                    );

                    const EigVec Yprev = fieldToEigVec(VmGuess);
                    eigVecToField
                    (
                        Yprev + nonlinearRelaxation*(Ysol - Yprev), VmGuess
                    );
                    VmGuess.correctBoundaryConditions();
                    evaluateNonlinearFields
                    (
                        VmGuess, statesGuess, IionGuess, true
                    );

                    const EigVec rhsNew = rhsBase + stageSource(IionGuess);
                    const scalar coupledResidual = relativeL2Norm
                    (
                        applyEsdirkOp(fieldToEigVec(VmGuess)) - rhsNew, rhsNew
                    );
                    const scalar VmResidual =
                        relativeL2Difference(VmGuess.primitiveField(), VmPrevious);
                    const scalar IionResidual =
                        relativeL2Difference
                        (
                            IionGuess.primitiveField(), IionPrevious, IionIncrementFloor
                        );
                    const scalar maxStateResidual = maxStateRelativeL2Difference
                    (
                        statesGuess,
                        statesPrevious,
                        stateResiduals
                    );

                    const bool converged = nonlinearIterationStops
                    (
                        nonlinearIterationExit
                        (
                            corr,
                            implicitMinNonlinearIterations,
                            coupledResidual,
                            coupledResidualPrevIter,
                            nonlinearTolerance,
                            VmResidual,
                            nonlinearVmTolerance,
                            IionResidual,
                            nonlinearIionTolerance,
                            nonlinearRequireStatesConvergence,
                            maxStateResidual,
                            nonlinearStatesTolerance
                        ),
                        coupledResidual,
                        nonlinearTolerance,
                        nResidualFloorExits
                    );
                    coupledResidualPrevIter = coupledResidual;

                    writeAndPrintResiduals
                    (
                        corr + 1,
                        linearIterations,
                        linearError,
                        coupledResidual,
                        VmResidual,
                        IionResidual,
                        maxStateResidual,
                        stateResiduals,
                        converged
                    );

                    ++nonlinearIters;
                    if (converged)
                    {
                        stageConverged = true;
                        break;
                    }
                }

                if (!stageConverged)
                {
                    nonlinearConverged = false;
                }

                // The fields are already evaluated at the accepted stage value
                // (the last thing the loop did), and F_i feeds every later stage.
                stageF[i] =
                    applyLOp(fieldToEigVec(VmGuess)) + stageSource(IionGuess);

                stageVm[i] = VmGuess.primitiveField();

                if (gaussHistory && i + 1 < nS)
                {
                    stageVmIntegrationPoints[i].setSize
                    (
                        VmOldIntegrationPoints.size()
                    );
                    reconstructVmLevelAtIionIntegrationPoints
                    (
                        stageVm[i],
                        VmLevelScratchPtr(),
                        useHighOrder_Vm,
                        LREInterp_VmPtr(),
                        LREInterp_IionPtr(),
                        stageVmIntegrationPoints[i],
                        stateODEUseOpenMP,
                        stateODEOpenMPThreshold,
                        profileTimings,
                        fineTimings
                    );
                }
            }

            // Restore the one-step defaults. stageVm is local to this block, so
            // a driver left pointing into it would dangle - and a dangling read
            // is a wrong answer, not a crash.
            stageDt = dt;
            cellDriver.nNodes = 2;
            cellDriver.s[0] = 0.0;
            cellDriver.s[1] = 1.0;
            cellDriver.field[0] = &VmOldValues;
            cellDriver.field[1] = nullptr;
            gaussDriver.nNodes = 2;
            gaussDriver.s[0] = 0.0;
            gaussDriver.s[1] = 1.0;
            gaussDriver.field[0] = &VmOldIntegrationPoints;
            gaussDriver.field[1] = nullptr;
        }
        else if (useJFNK)
        {
            EigVec x = fieldToEigVec(VmGuess);

            jfnkPc.beginStep();

            // IMPORTANT: residualFor and standardResidual must use an explicit
            // `-> EigVec` return type. Without it, the deduced return type is
            // Eigen's expression-template (CwiseBinaryOp<..., Product<...>>),
            // which holds *references* to the temporaries produced by
            // fieldToEigVec(...) and rhsFromIion(...). Those temporaries die
            // when the lambda returns, leaving the expression with dangling
            // pointers; later materialisation reads garbage sizes and
            // triggers std::bad_alloc (140 TB allocations have been observed).
            // Returning by value as EigVec forces the expression to be
            // evaluated *inside* the lambda while the temporaries are alive.
            auto residualFor =
            [&](const EigVec& xCandidate, const bool updateGuess) -> EigVec
            {
                // Standard JfNK residual: F(Vm) = A*Vm - rhs(Iion(Vm))
                // This avoids an inner linear solve inside every matVec evaluation.
                auto standardResidual =
                [&](volScalarField& VmCandidate,
                    Field<Field<scalar>>& statesCandidate,
                    volScalarField& IionCandidate) -> EigVec
                {
                    evaluateNonlinearFields(VmCandidate, statesCandidate, IionCandidate, true);
                    // Modified for cardiacFoam: AImplicit has LOCAL rows and GLOBAL
                    // columns, so Eigen cannot form this product - it would
                    // need the whole global vector. PETSc's MatMult does the
                    // halo exchange internally.
                    return applyAImplicit(fieldToEigVec(VmCandidate))
                         - rhsFromIion(IionCandidate);
                };

                if (updateGuess)
                {
                    eigVecToField(xCandidate, VmGuess);
                    VmGuess.correctBoundaryConditions();
                    return standardResidual(VmGuess, statesGuess, IionGuess);
                }

                volScalarField VmTmp
                (
                    IOobject
                    (
                        "VmTmp",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    VmGuess
                );
                volScalarField IionTmp
                (
                    IOobject
                    (
                        "IionTmp",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    IionGuess
                );
                Field<Field<scalar>> statesTmp(statesOld);

                eigVecToField(xCandidate, VmTmp);
                VmTmp.correctBoundaryConditions();
                return standardResidual(VmTmp, statesTmp, IionTmp);
            };

            // The increments measured at the END of the previous corr, which are
            // what the check at the top of the next one is entitled to look at.
            //
            // That check used to recompute them against a copy taken two lines
            // further up, with only residualFor(x, true) in between - and that
            // call re-evaluates the same x the bottom of the previous iteration
            // already evaluated, with the states reset to statesOld both times.
            // So the "increments" were the difference between two evaluations of
            // one point: zero, or RKF45 step-history noise. Four of the five
            // tolerance terms could never fail, the test reduced to
            // `coupledResidual <= nonlinearTolerance`, and JFNK returned the
            // jfnkInitGuessOrder extrapolation without taking a single Newton
            // step whenever it landed under that tolerance. Same defect as the
            // MMS solver, where it was found first: on its M08 in 3-D, 24 of 50
            // timesteps exited that way.
            //
            // GREAT on entry means the first pass can never exit there, so at
            // least one Newton step is always taken. stateResiduals needs no
            // copy: only the bottom-of-loop call writes it, so on an early exit
            // it still holds the per-state values that go with
            // maxStateIncrPrevIter.
            scalar VmIncrPrevIter = GREAT;
            scalar IionIncrPrevIter = GREAT;
            scalar maxStateIncrPrevIter = GREAT;

            // Note for cardiacFoam: the MMS solver refuses its Eigen JFNK path
            // in parallel because its GMRES does not reduce its inner products.
            // This one does - solveGMRES goes through gDot/gNorm and caps the
            // Krylov dimension uniformly across ranks - and its ILUT
            // preconditioner is factorised on the local diagonal block, so the
            // Eigen JFNK path here IS parallel and deliberately has no such
            // guard. Only the assembled Eigen path (solveSparseSystem) does.

            scalar coupledResidualPrevIter = GREAT;

            for
            (
                label corr = 0;
                corr < max(implicitNonlinearIterations, label(1));
                ++corr
            )
            {
                const scalarField VmPrevious(VmGuess.primitiveField());
                const scalarField IionPrevious(IionGuess.primitiveField());
                Field<Field<scalar>> statesPrevious(statesGuess);

                const EigVec R = residualFor(x, true);
                const scalar coupledResidual = relativeL2Norm(R, x);

                // The increments that decide whether this iteration can be
                // skipped are the ones the PREVIOUS one measured; recomputing
                // them here would measure zero, see the note above the loop.
                const scalar VmResidualInitial = VmIncrPrevIter;
                const scalar IionResidualInitial = IionIncrPrevIter;
                const scalar maxStateResidualInitial = maxStateIncrPrevIter;

                // Early-exit check: skip the GMRES solve when the current
                // iterate meets the coupled tolerance AND the increments of the
                // previous iteration met theirs. It cannot fire at corr == 0.
                const bool residualConverged =
                    corr + 1 >= implicitMinNonlinearIterations
                 && coupledResidual <= nonlinearTolerance
                 && VmResidualInitial <= nonlinearVmTolerance
                 && IionResidualInitial <= nonlinearIionTolerance
                 && (
                        !nonlinearRequireStatesConvergence
                     || maxStateResidualInitial <= nonlinearStatesTolerance
                    );

                if (residualConverged)
                {
                    writeAndPrintResiduals
                    (
                        corr + 1,
                        0,
                        0.0,
                        coupledResidual,
                        VmResidualInitial,
                        IionResidualInitial,
                        maxStateResidualInitial,
                        stateResiduals,
                        true
                    );
                    nonlinearConverged = true;
                    nonlinearIters = corr + 1;
                    break;
                }

                // Finite-difference approximation to the Jacobian-vector
                // product:    J*v ~= (F(x + eps*v) - F(x)) / eps
                // The forward-difference perturbation follows the standard
                // JFNK formula  eps = sqrt(eps_mach) * (1 + ||x||) / ||v|| ,
                // which balances truncation and round-off error and is
                // invariant under scaling of v.
                //
                // The explicit `-> EigVec` return type is critical: see the
                // long comment on residualFor() above. Without it, the
                // returned CwiseBinaryOp captures references to R (safe,
                // it lives outside) AND to the temporary returned by
                // residualFor(..., false), causing dangling-pointer
                // bad_alloc inside Eigen when the result is finally evaluated.
                const scalar epsScale = jfnkEpsilon*(1.0 + gNorm(x));
                auto matVec = [&](const EigVec& v) -> EigVec
                {
                    const scalar eps = epsScale/max(gNorm(v), SMALL);
                    return (residualFor(x + eps*v, false) - R)/eps;
                };

                // Solve  J * delta = -R  with restarted GMRES.
                // The matVec closure above re-evaluates the nonlinear
                // residual every Krylov step, so each iteration triggers a
                // full ionic-ODE integration; this dominates the cost.
                label gmresIterations = 0;
                scalar gmresError = GREAT;
                EigVec delta(R.size());

                // Shared with the MMS solver (monodomainKrylov.H): the
                // preconditioner is rebuilt at the first iteration of the
                // step and every jfnkPreconditionerUpdateFrequency after it.
                if (jfnkPc.needsUpdate(corr))
                {
                    if (jfnkPc.usesSourceDerivative())
                    {
                        computeSourceDerivative(VmGuess, IionGuess);
                    }
                    jfnkPc.update
                    (
                        AImplicit,
                        jfnkPc.usesSourceDerivative()
                      ? &sourceDerivative.primitiveField()
                      : nullptr,
                        theta,
                        corr
                    );
                }

                delta = solveJfnkLinearSystem
                (
                    matVec,
                    R,
                    jfnkPc,
                    jfnkPetscShellSolver,
                    jfnkKrylov,
                    gmresIterations,
                    gmresError
                );

                // ----------------------------------------------------- //
                // Newton update with optional Armijo backtracking.
                //
                // The unconditional damped step  x <- x + omega * delta
                // is fragile when the linearisation of F at the current
                // iterate fails to capture the true nonlinearity (typical
                // case: cardiac AP upstroke with large dt). We optionally
                // wrap the step in a backtracking line search that halves
                // alpha until the sufficient-decrease Armijo condition
                //
                //     ||F(x + alpha*delta)||^2
                //         <= (1 - 2*c*alpha) * ||F(x)||^2
                //
                // is satisfied. This guarantees monotone descent in the
                // residual norm and globalises Newton.
                // ----------------------------------------------------- //
                label lineSearchIters = 0;
                EigVec Rnew(R.size());
                scalar coupledResidualNew = GREAT;

                if (jfnkLineSearch)
                {
                    // Snapshot of the current iterate; restored on each
                    // backtracking trial.
                    const EigVec xBackup = x;
                    const scalar Rold2 = gSquaredNorm(R);
                    scalar alpha = nonlinearRelaxation;

                    while (true)
                    {
                        // x = xBackup + alpha * delta
                        x = xBackup + alpha*delta;

                        // F(x_new); commits VmGuess/statesGuess/IionGuess
                        Rnew = residualFor(x, true);
                        const scalar Rnew2 = gSquaredNorm(Rnew);

                        // Armijo sufficient-decrease test
                        if (Rnew2 <= (1.0 - 2.0*jfnkArmijoC*alpha)*Rold2)
                        {
                            break;  // accept trial
                        }

                        ++lineSearchIters;
                        if (lineSearchIters >= jfnkLineSearchMaxIter)
                        {
                            break;  // budget exhausted
                        }

                        alpha *= 0.5;
                        if (alpha < jfnkLineSearchAlphaMin)
                        {
                            // alpha floored; do one final trial with the
                            // floor and accept whatever comes out.
                            alpha = jfnkLineSearchAlphaMin;
                            x = xBackup + alpha*delta;
                            Rnew = residualFor(x, true);
                            break;
                        }
                    }
                    coupledResidualNew = relativeL2Norm(Rnew, x);
                }
                else
                {
                    // Legacy path: unconditional damped step.
                    x += nonlinearRelaxation*delta;
                    Rnew = residualFor(x, true);
                    coupledResidualNew = relativeL2Norm(Rnew, x);
                }

                const scalar VmResidual =
                    relativeL2Difference(VmGuess.primitiveField(), VmPrevious);
                const scalar IionResidual =
                    relativeL2Difference
                    (
                        IionGuess.primitiveField(), IionPrevious, IionIncrementFloor
                    );
                const scalar maxStateResidual = maxStateRelativeL2Difference
                (
                    statesGuess,
                    statesPrevious,
                    stateResiduals
                );

                // Shared stopping test (monodomainBase.H); this used to be an
                // inline copy of it, without the residual-floor exit.
                const bool converged = nonlinearIterationStops
                (
                    nonlinearIterationExit
                    (
                        corr,
                        implicitMinNonlinearIterations,
                        coupledResidualNew,
                        coupledResidualPrevIter,
                        nonlinearTolerance,
                        VmResidual,
                        nonlinearVmTolerance,
                        IionResidual,
                        nonlinearIionTolerance,
                        nonlinearRequireStatesConvergence,
                        maxStateResidual,
                        nonlinearStatesTolerance
                    ),
                    coupledResidualNew,
                    nonlinearTolerance,
                    nResidualFloorExits
                );
                coupledResidualPrevIter = coupledResidualNew;

                writeAndPrintResiduals
                (
                    corr + 1,
                    gmresIterations,
                    gmresError,
                    coupledResidualNew,
                    VmResidual,
                    IionResidual,
                    maxStateResidual,
                    stateResiduals,
                    converged,
                    lineSearchIters
                );

                // What the top-of-loop check of the next iteration may look at.
                VmIncrPrevIter = VmResidual;
                IionIncrPrevIter = IionResidual;
                maxStateIncrPrevIter = maxStateResidual;

                nonlinearIters = corr + 1;
                if (converged)
                {
                    nonlinearConverged = true;
                    break;
                }
            }
        }
        else
        {
            scalar coupledResidualPrevIter = GREAT;

            for
            (
                label corr = 0;
                corr < max(implicitNonlinearIterations, label(1));
                ++corr
            )
            {
                const scalarField VmPrevious(VmGuess.primitiveField());
                const scalarField IionPrevious(IionGuess.primitiveField());
                Field<Field<scalar>> statesPrevious(statesGuess);

                evaluateNonlinearFields(VmGuess, statesGuess, IionGuess, true);
                EigVec rhs = rhsFromIion(IionGuess);
                SpMat ACurrent = AImplicit;

                if (useDiagonalIion)
                {
                    computeSourceDerivative(VmGuess, IionGuess);
                    const EigVec sourceDerivativeVec = fieldToEigVec(sourceDerivative);
                    const EigVec VmLinearisationPoint = fieldToEigVec(VmGuess);

                    addDiagonalToMatrix
                    (
                        sourceDerivative.primitiveField(),
                       -theta,
                        ACurrent
                    );

                    rhs -= theta
                       * sourceDerivativeVec.cwiseProduct(VmLinearisationPoint);
                }

                label linearIterations = 0;
                scalar linearError = GREAT;

                // Phase B: decide how the persistent solver treats the matrix.
                //   - first solve, or dt changed  -> Rebuild (full setup)
                //   - diagonalIion (diagonal moves every iteration) -> ValuesOnly
                //   - Picard with a constant operator -> Reuse (solve only)
                MatrixUpdate updatePolicy;
                if (linearSolverNeedsRebuild)
                {
                    updatePolicy = MatrixUpdate::Rebuild;
                    linearSolverNeedsRebuild = false;
                }
                else if (useDiagonalIion)
                {
                    updatePolicy = MatrixUpdate::ValuesOnly;
                }
                else
                {
                    updatePolicy = MatrixUpdate::Reuse;
                }

                const EigVec Vsol = solveSparseSystem
                (
                    ACurrent,
                    rhs,
                    linearSolverBackend,
                    implicitLinearSolver,
                    petscLinearKspType,
                    petscLinearPcType,
                    petscLinearRestart,
                    petscLinearOptionsPrefix,
                    petscUseOptions,
                    implicitTolerance,
                    implicitMaxIterations,
                    linearIterations,
                    linearError,
                    persistentPetscSolver,
                    persistentEigenSolvers,
                    updatePolicy,
                    petscReusePreconditioner
                );

                const EigVec Vprev = fieldToEigVec(VmGuess);
                const EigVec Vrelaxed =
                    Vprev + nonlinearRelaxation*(Vsol - Vprev);

                eigVecToField(Vrelaxed, VmGuess);
                VmGuess.correctBoundaryConditions();
                evaluateNonlinearFields(VmGuess, statesGuess, IionGuess, true);

                const EigVec nonlinearResidualVec =
                    applyAImplicit(fieldToEigVec(VmGuess)) - rhsFromIion(IionGuess);
                const scalar coupledResidual =
                    relativeL2Norm(nonlinearResidualVec, rhsFromIion(IionGuess));
                const scalar VmResidual =
                    relativeL2Difference(VmGuess.primitiveField(), VmPrevious);
                const scalar IionResidual =
                    relativeL2Difference
                    (
                        IionGuess.primitiveField(), IionPrevious, IionIncrementFloor
                    );
                const scalar maxStateResidual = maxStateRelativeL2Difference
                (
                    statesGuess,
                    statesPrevious,
                    stateResiduals
                );

                const bool converged = nonlinearIterationStops
                (
                    nonlinearIterationExit
                    (
                        corr,
                        implicitMinNonlinearIterations,
                        coupledResidual,
                        coupledResidualPrevIter,
                        nonlinearTolerance,
                        VmResidual,
                        nonlinearVmTolerance,
                        IionResidual,
                        nonlinearIionTolerance,
                        nonlinearRequireStatesConvergence,
                        maxStateResidual,
                        nonlinearStatesTolerance
                    ),
                    coupledResidual,
                    nonlinearTolerance,
                    nResidualFloorExits
                );
                coupledResidualPrevIter = coupledResidual;

                writeAndPrintResiduals
                (
                    corr + 1,
                    linearIterations,
                    linearError,
                    coupledResidual,
                    VmResidual,
                    IionResidual,
                    maxStateResidual,
                    stateResiduals,
                    converged
                );

                nonlinearIters = corr + 1;
                if (converged)
                {
                    nonlinearConverged = true;
                    break;
                }
            }
        }

        // ----------------------------------------------------------------- //
        // Commit or roll back the nonlinear iterate.
        //
        // If the Newton/Picard loop exhausted its iteration budget without
        // reaching the convergence tolerances, we have two unsafe options
        // available:
        //
        //   (a) accept the last iterate, which may be far from the implicit
        //       solution and pollute every subsequent time step;
        //   (b) discard the time step's progress and reset to t_n values.
        //
        // We choose (b) by default: print a clear warning and restore Vm,
        // Iion and the ionic states from their time-step-start snapshots.
        // The user can opt back into (a) with `nonlinearAcceptUnconverged`
        // in spatialIntegrationProperties.
        // ----------------------------------------------------------------- //
        if (!nonlinearConverged && !nonlinearAcceptUnconverged)
        {
            WarningInFunction
                << "Nonlinear solver did not converge at t = " << t0
                << " after " << nonlinearIters << " iterations."
                << " Rolling back Vm, Iion and ionic states to the values"
                << " at the beginning of the time step." << nl
                << "Set 'nonlinearAcceptUnconverged true;' in"
                << " spatialIntegrationProperties to accept unconverged"
                << " iterates instead." << endl;

            VmGuess.primitiveFieldRef() = VmOldValues;
            VmGuess.correctBoundaryConditions();
            IionGuess.primitiveFieldRef() = IionOld.primitiveField();
            IionGuess.correctBoundaryConditions();
            statesGuess = statesOld;

            if (tnnpModel)
            {
                tnnpModel->resetStatesToStatesOld(statesGuess);
            }

            // The next step starts from the same state, and away from a
            // stimulus edge the cell model is autonomous: it is the same
            // nonlinear problem and fails the same way. Stop instead of
            // carrying a frozen solution to the end time.
            ++nRolledBackSteps;
            if (++nConsecutiveRollbacks > nonlinearMaxConsecutiveRollbacks)
            {
                FatalErrorInFunction
                    << nConsecutiveRollbacks << " consecutive time steps"
                    << " rolled back (nonlinearMaxConsecutiveRollbacks = "
                    << nonlinearMaxConsecutiveRollbacks << "), the last at t = "
                    << t0 << ". The solution has been frozen since the first"
                    << " of them; tighten dt or loosen the nonlinear"
                    << " tolerances." << exit(FatalError);
            }
        }
        else if (!nonlinearConverged && nonlinearAcceptUnconverged)
        {
            WarningInFunction
                << "Nonlinear solver did not converge at t = " << t0
                << " after " << nonlinearIters << " iterations."
                << " Accepting the unconverged iterate"
                << " (nonlinearAcceptUnconverged=true)." << endl;
        }

        if (nonlinearConverged || nonlinearAcceptUnconverged)
        {
            nConsecutiveRollbacks = 0;
        }

        // Commit the (possibly rolled-back) iterate into the persistent
        // fields used by subsequent time steps and post-processing.
        Vm.primitiveFieldRef() = VmGuess.primitiveField();
        Vm.correctBoundaryConditions();
        Iion.primitiveFieldRef() = IionGuess.primitiveField();
        Iion.correctBoundaryConditions();
        states = statesGuess;

        if (nonlinearIionIncrementFloor)
        {
            IionIncrementFloor = max
            (
                IionIncrementFloor,
                Foam::sqrt(gSumSqr(Iion.primitiveField()))
            );
        }

        // The hybrid's cell-centred states move on with the step, unless it
        // was rolled back.
        if (useFrontHybrid && (nonlinearConverged || nonlinearAcceptUnconverged))
        {
            hybridCellStates = hybridCellStatesGuess;
        }

        // Added for cardiacFoam, ported from the MMS solver: rotate the BDF
        // history once the step is final. V^n becomes V^{n-1} and the oldest
        // level falls off the end.
        //
        // Two things RESTART it instead (a departure from the MMS, which
        // rotates regardless): a rolled-back step, which advanced time but not
        // the solution, so the history would hold two identical levels one dt
        // apart; and a step of non-nominal length, after which the levels are
        // no longer equally spaced. A restart costs k-1 ESDIRK3 steps; a
        // history that is silently inconsistent costs the order.
        if (bdfOrder > 1)
        {
            const bool rolledBack =
                !nonlinearConverged && !nonlinearAcceptUnconverged;

            if (rolledBack || clippedStep)
            {
                bdfValidLevels = 0;
            }
            else
            {
                for (label j = VmBdfOlder.size() - 1; j >= 1; --j)
                {
                    VmBdfOlder[j] = VmBdfOlder[j - 1];
                }
                VmBdfOlder[0] = VmOldValues;

                bdfValidLevels = min(bdfValidLevels + 1, bdfOrder - 1);
            }
        }

        if (outFields.size())
        {
            if (useHighOrder_Iion && !statesAtCells)
            {
                ionicModel->exportStatesIntegrationPoints
                (
                    states,
                    outFields,
                    LREInterp_IionPtr().cellQuadWeightPhysical()
                );
            }
            else
            {
                // cellCentredReconstruct (and the low-order path) keep the
                // authoritative states at cell centres.
                ionicModel->exportStates(states, outFields);
            }
        }

        if (useHighOrder_Vm && dim > 1)
        {
            computeHighOrderLaplacian
            (
                Vm,
                conductivity,
                LREInterp_VmPtr(),
                fluxVm_HO,
                lapVm
            );
        }
        else
        {
            lapVm = fvc::laplacian(conductivity, Vm);
        }
        rhsVm = lapVm/(chi*Cm) - Iion + externalStimulusCurrent/(chi*Cm);

        updateActivationTimes
        (
            VmOldValues,
            Vm,
            activationTime,
            calculateActivationTime,
            t0,
            dt,
            activationThreshold
        );

        // Refresh the live activation samples whenever new points crossed the
        // threshold this step (files are overwritten; quiet to avoid log spam).
        label nActivated = 0;
        forAll(calculateActivationTime, cellI)
        {
            if (!calculateActivationTime[cellI])
            {
                ++nActivated;
            }
        }
        reduce(nActivated, sumOp<label>());

        if (nActivated > nActivatedPrev)
        {
            nActivatedPrev = nActivated;
            tLastActivationChange = runTime.value();
            writeActivationSamples
            (
                runTime,
                activationTime,
                diagonalStart,
                diagonalEnd,
                p8Point,
                nDiagonalSamples,
                sampleNeighbours,
                samplePower,
                false
            );
        }

        // Early-abort on activation failure (recorded in activationStatus.dat).
        if (abortOnActivationFailure && !stopPointTriggered)
        {
            const scalar tNow = runTime.value();
            if (nActivated == 0)
            {
                if (tNow > activationIgnitionTimeout)
                {
                    activationStatus = "FAILED_NO_IGNITION";
                    activationAborted = true;
                }
            }
            else if
            (
                stopAfterPointActivation
             && (tNow - tLastActivationChange) > activationStallTime
            )
            {
                activationStatus = "FAILED_CONDUCTION_BLOCK";
                activationAborted = true;
            }

            if (activationAborted)
            {
                WarningInFunction
                    << "Activation failure detected (" << activationStatus
                    << ") at t = " << tNow << " s, nActivated = " << nActivated
                    << "; aborting the run early. See activationStatus.dat."
                    << nl << endl;
                break;
            }
        }

        // Modified for cardiacFoam: only ONE rank owns the stop cell now, so
        // the test can no longer be gated on stopActivationCellI >= 0 - the
        // other ranks would never trigger. The owner's value is broadcast with
        // a max-reduction (every other rank contributes -GREAT), which makes
        // stopPointTriggered, stopPointActivationTime and effectiveEndTime
        // uniform BY CONSTRUCTION, and with them the loop's exit step,
        // FAILED_CONDUCTION_BLOCK and the final activationStatus.
        if (stopAfterPointActivation && !stopPointTriggered)
        {
            scalar actTime =
                (stopActivationCellI >= 0)
              ? activationTime[stopActivationCellI]
              : -GREAT;
            reduce(actTime, maxOp<scalar>());

            if (actTime > SMALL)
            {
                stopPointTriggered = true;
                stopPointActivationTime = actTime;
                effectiveEndTime = min
                (
                    effectiveEndTime,
                    stopPointActivationTime + stopDelayAfterActivation
                );

                Info<< "Stop activation point crossed threshold at t = "
                    << stopPointActivationTime << " s; solver will stop at t = "
                    << effectiveEndTime << " s" << nl << endl;
            }
        }

        activationVelocity = fvc::grad
        (
            1.0/(activationTime + dimensionedScalar("SMALL", dimTime, SMALL))
        );

        ++nSteps;
        ++runTime;

        // Only the last step can be clipped (or the one that reaches the stop
        // point), so up to it every step is nominal and n*dt is exact.
        runTime.setTime
        (
            clippedStep ? t0 + dt : loopStartTime + nSteps*nominalDt,
            runTime.timeIndex()
        );

        if (nSteps <= 5 || nSteps % stepLogInterval == 0)
        {
            Info<< "Step " << nSteps
                << " time = " << runTime.value()
                << " min(Vm) = " << gMin(Vm.primitiveField())
                << " max(Vm) = " << gMax(Vm.primitiveField()) << nl;
        }

        runTime.write();
    }

    activationVelocity = fvc::grad
    (
        1.0/(activationTime + dimensionedScalar("SMALL", dimTime, SMALL))
    );
    activationTime.write();
    activationVelocity.write();
    Vm.write();
    Iion.write();
    externalStimulusCurrent.write();
    lapVm.write();
    rhsVm.write();
    forAll(outFields, i)
    {
        outFields[i].write();
    }

    if (stopAfterPointActivation && !stopPointTriggered)
    {
        Info<< "Stop activation point did not cross the threshold before t = "
            << runTime.value() << " s" << nl << endl;
    }

    writeActivationSamples
    (
        runTime,
        activationTime,
        diagonalStart,
        diagonalEnd,
        p8Point,
        nDiagonalSamples,
        sampleNeighbours,
        samplePower
    );

    // Finalise and record the activation status so a driver can tell a healthy
    // run from a non-igniting / conduction-blocked one without parsing the log.
    if (!activationAborted)
    {
        if (stopPointTriggered)
        {
            activationStatus = "OK";
        }
        else if (nActivatedPrev == 0)
        {
            activationStatus = "FAILED_NO_IGNITION";
        }
        else
        {
            activationStatus = "INCOMPLETE";
        }
    }
    {
        // Modified for cardiacFoam: see postProcessingDir().
        const fileName statusDir(postProcessingDir(runTime));
        mkDir(statusDir);
        OFstream statusFile(statusDir/"activationStatus.dat");
        statusFile
            << "# status finalTime_s nActivated stopPointActivated" << nl
            << activationStatus << ' ' << runTime.value() << ' '
            << nActivatedPrev << ' ' << (stopPointTriggered ? 1 : 0) << nl;
        Info<< "Activation status: " << activationStatus
            << " (t = " << runTime.value() << " s, nActivated = "
            << nActivatedPrev << ")" << nl << endl;
    }

    const auto tEndTotal = std::chrono::steady_clock::now();
    const scalar setupWallTime =
        std::chrono::duration<scalar>(tEndSetup - tStartSetup).count();
    const scalar loopWallTime =
        std::chrono::duration<scalar>(tEndTotal - tStartLoop).count();
    const scalar totalWallTime =
        std::chrono::duration<scalar>(tEndTotal - tStartTotal).count();
    const scalar peakRSSMB = currentPeakRSSMB();

    const fileName performanceDir
    (
        postProcessingDir(runTime)   // Modified for cardiacFoam
    );
    mkDir(performanceDir);
    OFstream perfFile(performanceDir/"solverPerformance.dat");
    perfFile
        << "# solver finalTime_s nSteps setupWall_s loopWall_s totalWall_s peakRSS_MB "
        << "memoryOptimization memoryOptimizationEffective "
        << "linearSolverBackend petscLinearKspType petscLinearPcType "
        << "jfnkLinearSolverBackend jfnkPetscKspType jfnkPetscPcType"
        << nl;
    perfFile
        << "highOrderElectroActivationFoamImplicitPETSc "
        << runTime.value() << ' '
        << nSteps << ' '
        << setupWallTime << ' '
        << loopWallTime << ' '
        << totalWallTime << ' '
        << peakRSSMB << ' '
        << memoryOptimization << ' '
        << (memoryOptimizationEffective ? "true" : "false") << ' '
        << linearSolverBackend << ' '
        << petscLinearKspType << ' '
        << petscLinearPcType << ' '
        << jfnkLinearSolverBackend << ' '
        << jfnkPetscKspType << ' '
        << jfnkPetscPcType
        << nl;

    Info<< "Solver performance written to " << performanceDir
        << "/solverPerformance.dat" << nl
        << "Wall time: setup = " << setupWallTime
        << " s, loop = " << loopWallTime
        << " s, total = " << totalWallTime << " s" << nl
        << "Peak RSS = " << peakRSSMB << " MB" << nl << endl;

    runTime.printExecutionTime(Info);

    if (useFrontHybrid)
    {
        label nSwitch = 0;
        label nMany = 0;
        label maxSwitch = 0;
        forAll(frontSwitches, cellI)
        {
            nSwitch += frontSwitches[cellI];
            if (frontSwitches[cellI] > 2) ++nMany;
            maxSwitch = max(maxSwitch, frontSwitches[cellI]);
        }
        reduce(nSwitch, sumOp<label>());
        reduce(nMany, sumOp<label>());
        reduce(maxSwitch, maxOp<label>());
        Info<< "Hybrid front band: " << nSwitch << " mode switches; " << nMany
            << " cells switched more than twice (max " << maxSwitch << ")" << nl;
    }

    Info<< "Nonlinear iterations stopped at the residual floor: "
        << nResidualFloorExits << nl;

    Info<< "Rolled-back time steps: " << nRolledBackSteps << nl;

    Info<< "End" << nl << endl;

    return 0;
}
