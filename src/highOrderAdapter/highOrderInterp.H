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

Class
    Foam::highOrderInterp

Description
    Added for cardiacFoam: a thin compatibility adapter that exposes the small
    part of the old, serial LRE API that the cardiac solvers actually use, on
    top of solids4foam's parallel movingLeastSquares.

    The old library (Castrillo et al., Comput Struct 268:106829, 2022) lived in
    solids4foam as numerics/LRE and had no processor-boundary support: its
    stencil construction returns before the (commented out) MPI exchange, and
    processor patches raise NotImplemented. movingLeastSquares replaces it with
    a genuinely parallel implementation in which stencils carry global cell IDs
    and every evaluation entry point performs its own halo exchange, so a caller
    needs no changes to run in parallel.

    This adapter exists so the ~7000-line cardiac solvers can be moved across in
    one place rather than at every call site, and so the dictionary translation
    lives in a single documented location. It is deliberately temporary: once
    the cardiac solvers are ported and verified, the call sites can be pointed
    straight at movingLeastSquares and this header deleted.

    IMPORTANT - quadrature weights changed meaning.
    LRE normalised its quadrature weights (face weights summed to 1 over a face,
    cell weights to 1 over a cell), so callers multiplied by |Sf| or by the cell
    volume themselves. fvMeshQuadrature returns PHYSICAL weights: face weights
    sum to |Sf| and cell weights sum to the cell volume. To stop that difference
    from silently rescaling a result, the weight accessors here are deliberately
    NOT called faceQuadWeight()/cellQuadWeight(): a caller ported from LRE gets
    a compile error and has to visit the site.

    Not provided, because the cardiac solvers never used them: QRInterpCoeffs,
    assembleStiffnessMatrixCoeffs, assembleLumpedMassMatrix,
    surfaceIntegrateVector, and every vector-field overload.

Authors
    Pablo Castrillo, UCD.

\*---------------------------------------------------------------------------*/

#ifndef highOrderInterp_H
#define highOrderInterp_H

#include "fvMesh.H"
#include "volFields.H"
#include "CompactListList.H"
#include "movingLeastSquares.H"
#include "movingLeastSquaresStencil.H"
#include "fvMeshQuadrature.H"
#include "symmTensor3rdOrder.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class highOrderInterp Declaration
\*---------------------------------------------------------------------------*/

class highOrderInterp
{
    // Private data

        //- The wrapped parallel moving-least-squares engine
        movingLeastSquares mls_;


    // Private Member Functions

        //- Translate an LRE-style coefficients dictionary into the keys
        //  movingLeastSquares expects.
        //
        //      N                    -> polynomialOrder
        //      Nn                   -> faceStencilExtraCells
        //                              cellStencilExtraCells
        //      weightFunction       -> weightFunctionCoeffs/type
        //      k                    -> weightFunctionCoeffs/k
        //      calcConditionNumber  -> calcConditionNumber (unchanged)
        //      haloDepthScale       -> haloDepthScale (new, optional)
        //
        //  Both libraries define Nn as the number of stencil points ADDED on
        //  top of the minimum needed for the polynomial order, so it carries
        //  over unchanged.
        //
        //  Dropped, because movingLeastSquares has no equivalent and does not
        //  need one: maxStencilSize and nLayers only ever sized a DynamicList
        //  in LRE's dead legacy path; useQRDecomposition and useGlobalStencils
        //  selected between LRE code paths of which only one worked.
        static dictionary translateDict(const dictionary& lreDict);

        //- Disallow default bitwise copy construct and assignment
        highOrderInterp(const highOrderInterp&);
        void operator=(const highOrderInterp&);


public:

    // Public typedefs

        //- Fully symmetric third-order tensor, ten independent components in
        //  the order xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz. LRE used
        //  a FixedList<scalar, 10> with exactly this layout, so callers that
        //  only pass it to cubicForm() need no change.
        typedef symmTensor3rdOrder symmTensor3Order;


    // Constructors

        //- Construct from mesh, per-patch stencil inclusion flags and an
        //  LRE-style coefficients dictionary
        highOrderInterp
        (
            const fvMesh& mesh,
            const boolList& includePatchInStencils,
            const dictionary& dict
        );


    // Member Functions

        //- Return the wrapped engine, for callers that want the native API
        const movingLeastSquares& mls() const
        {
            return mls_;
        }

        //- Polynomial order of the local reconstruction
        label order() const
        {
            return mls_.polynomialOrder();
        }

        //- Contract a symmetric third-order tensor with a vector three times,
        //  T_ijk d_i d_j d_k. Identical to LRE::cubicForm, including the twoD
        //  flag that suppresses the out-of-plane component.
        inline static scalar cubicForm
        (
            const symmTensor3Order& T,
            const vector& d,
            const bool twoD
        )
        {
            return Foam::cubicForm
            (
                T,
                twoD ? vector(d.x(), d.y(), 0.0) : d
            );
        }


    // Differential operators at cell centres

        //- Gradient of a scalar field at the cell centres
        tmp<volVectorField> grad(const volScalarField& s) const
        {
            return mls_.grad(s);
        }

        //- Second derivative (Hessian) of a scalar field at the cell centres
        tmp<volSymmTensorField> hessian(const volScalarField& s) const
        {
            return mls_.secondGrad(s);
        }

        //- Third derivative of a scalar field at the cell centres
        autoPtr<List<symmTensor3Order>> thirdDeriv(const volScalarField& s) const
        {
            return mls_.thirdGrad(s);
        }

        //- Gradient of a scalar field at every face quadrature point, indexed
        //  [faceI][quadPointI]. Replaces LRE::gradScalarFaceQuad; the return
        //  type is a CompactListList rather than a List<List>, which indexes
        //  identically but stores contiguously.
        autoPtr<CompactListList<vector>> gradScalarFaceQuad
        (
            const volScalarField& s
        ) const;


    // Quadrature

        //- Face quadrature points, indexed [faceI][quadPointI]
        const CompactListList<point>& faceQuadPoints() const
        {
            return mls_.quadrature().faceQuadPoints();
        }

        //- Cell quadrature points, indexed [cellI][quadPointI]
        const CompactListList<point>& cellQuadPoints() const
        {
            return mls_.quadrature().cellQuadPoints();
        }

        //- Face quadrature weights. PHYSICAL: they sum to |Sf| over a face, so
        //  do NOT multiply by the face area as was needed with LRE.
        const CompactListList<scalar>& faceQuadWeightPhysical() const
        {
            return mls_.quadrature().faceQuadWeights();
        }

        //- Cell quadrature weights. PHYSICAL: they sum to the cell volume.
        //  Callers that renormalise by the row sum are unaffected.
        const CompactListList<scalar>& cellQuadWeightPhysical() const
        {
            return mls_.quadrature().cellQuadWeights();
        }


    // Stencils and coefficients, for callers assembling their own operators

        //- Per-cell stencil of GLOBAL cell indices
        const CompactListList<label>& globalCellStencils() const
        {
            return mls_.stencilData().cellsStencil();
        }

        //- Per-face stencil of GLOBAL cell indices
        const CompactListList<label>& globalFaceStencils() const
        {
            return mls_.stencilData().facesStencil();
        }

        //- Global cell addressing, for translating local to global indices
        const globalIndex& globalCells() const
        {
            return mls_.stencilData().globalCells();
        }

        //- Cell-centre gradient coefficients. Row cellI follows
        //  globalCellStencils()[cellI] and carries one extra trailing entry,
        //  at index globalCellStencils()[cellI].size(), for the cell itself.
        const CompactListList<vector>& QRGradCoeffs() const
        {
            return mls_.cellGradCoeffsData();
        }

        //- Cell-centre second derivative coefficients, same row layout.
        //  Empty when polynomialOrder < 2, in which case there is no second
        //  derivative to speak of.
        //
        //  LRE allocated these arrays unconditionally and left them unused at
        //  low order, so callers could bind the reference and then guard only
        //  the *use* with order() >= 2. movingLeastSquares instead never builds
        //  them below the required order and aborts inside calcCellCoeffs() if
        //  they are requested. Returning an empty list here preserves the LRE
        //  contract, so a caller that binds unconditionally still works.
        const CompactListList<symmTensor>& cellHessianCoeffs() const
        {
            static const CompactListList<symmTensor> empty;
            return order() >= 2 ? mls_.cellSecondGradCoeffsData() : empty;
        }

        //- Cell-centre third derivative coefficients, same row layout.
        //  Empty when polynomialOrder < 3; see cellHessianCoeffs().
        const CompactListList<symmTensor3Order>& cellThirdDerivCoeffs() const
        {
            static const CompactListList<symmTensor3Order> empty;
            return order() >= 3 ? mls_.cellThirdGradCoeffsData() : empty;
        }

        //- Face gradient coefficients, indexed [faceI][quadPointI][stencilI].
        //  Columns follow globalFaceStencils()[faceI]; a Dirichlet boundary
        //  face carries one extra trailing ghost column at index
        //  globalFaceStencils()[faceI].size().
        const List<CompactListList<vector>>& QRGradFaceGPCoeffs() const
        {
            return mls_.faceGradCoeffsData();
        }
};


// * * * * * * * * * * * * * * * * * Implementation  * * * * * * * * * * * * //
// Header-only on purpose: this is a thin, temporary adapter, so keeping it out
// of the build system means a consumer only needs an -I flag.

inline dictionary highOrderInterp::translateDict(const dictionary& lreDict)
{
    dictionary out;

    const label N = readLabel(lreDict.lookup("N"));
    const label Nn = readLabel(lreDict.lookup("Nn"));

    out.add("polynomialOrder", N);
    out.add("faceStencilExtraCells", Nn);
    out.add("cellStencilExtraCells", Nn);

    dictionary wfc;
    wfc.add
    (
        "type",
        lreDict.lookupOrDefault<word>("weightFunction", word("one"))
    );
    if (lreDict.found("k"))
    {
        wfc.add("k", readScalar(lreDict.lookup("k")));
    }
    out.add("weightFunctionCoeffs", wfc);

    out.add
    (
        "calcConditionNumber",
        lreDict.lookupOrDefault<Switch>("calcConditionNumber", false)
    );

    // Not an LRE key: movingLeastSquares defaults it to polynomialOrder*2.5.
    // Exposed here because it controls how far the stencil search reaches into
    // a neighbouring rank, and may need raising for high p on many ranks.
    if (lreDict.found("haloDepthScale"))
    {
        out.add("haloDepthScale", readScalar(lreDict.lookup("haloDepthScale")));
    }

    // LRE could treat a zeroGradient patch by mirroring the whole face stencil
    // across the boundary, opted into with this key. movingLeastSquares has the
    // mirroring machinery (it uses it for symmetry patches) but no equivalent
    // opt-in, and the cardiac solvers do not need one: they impose a
    // homogeneous Neumann condition by setting the boundary flux to zero, which
    // is the correct finite-volume treatment and is exact for the manufactured
    // solution. Warn rather than fail, so existing dictionaries keep working.
    if (lreDict.found("mirrorPatches"))
    {
        WarningInFunction
            << "Ignoring 'mirrorPatches': movingLeastSquares has no mirrored"
            << " zeroGradient treatment." << nl
            << "    Boundary faces on those patches get a one-sided stencil."
            << " For a homogeneous Neumann condition imposed as zero boundary"
            << " flux this makes no difference to the result." << endl;
    }

    return out;
}


inline highOrderInterp::highOrderInterp
(
    const fvMesh& mesh,
    const boolList& includePatchInStencils,
    const dictionary& dict
)
:
    mls_(mesh, includePatchInStencils, translateDict(dict))
{}


inline autoPtr<CompactListList<vector>>
highOrderInterp::gradScalarFaceQuad(const volScalarField& s) const
{
    // Size the result from the quadrature: one row per face, one entry per
    // quadrature point on that face. Empty patches contribute zero-length rows.
    const CompactListList<point>& faceQuadPts = faceQuadPoints();

    labelList rowSizes(faceQuadPts.size(), 0);
    forAll(faceQuadPts, faceI)
    {
        rowSizes[faceI] = faceQuadPts[faceI].size();
    }

    autoPtr<CompactListList<vector>> resultPtr
    (
        new CompactListList<vector>(rowSizes)
    );
    CompactListList<vector>& result = resultPtr();

    forAll(result, faceI)
    {
        forAll(result[faceI], qpI)
        {
            result[faceI][qpI] = vector::zero;
        }
    }

    // fGrad performs its own halo exchange, so this is parallel-correct with
    // nothing required of the caller. Note the two-argument form is used
    // deliberately: the autoPtr-returning overload in movingLeastSquares calls
    // a member that does not exist and does not compile if instantiated.
    mls_.fGrad(s, result);

    return resultPtr;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
