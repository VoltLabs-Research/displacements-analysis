#include <volt/displacements_engine.h>
#include <volt/analysis/cutoff_neighbor_finder.h>
#include <volt/core/progress.h>
#include <volt/displacements_log.h>

#include <stdexcept>
#include <vector>
#include <cmath>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

namespace Volt{

using namespace Particles;

DisplacementsEngine::DisplacementsEngine(
    ParticleProperty* positions,
    const SimulationCell& cell,
    ParticleProperty* refPositions,
    const SimulationCell& refCell,
    ParticleProperty* identifiers,
    ParticleProperty* refIdentifiers,
    bool useMinimumImageConvention,
    AffineMapping affineMapping
)
    : _positions(positions),
      _refPositions(refPositions),
      _identifiers(identifiers),
      _refIdentifiers(refIdentifiers),
      _simCell(cell),
      _simCellRef(refCell),
      _useMinimumImageConvention(useMinimumImageConvention),
      _affineMapping(affineMapping){}

void DisplacementsEngine::perform(DisplacementsLog& log){
    if(!_positions || !_refPositions){
        throw std::runtime_error("DisplacementsEngine: null input properties.");
    }

    const std::size_t n = _positions->size();
    if(n == 0){
        _displacementProperty.reset();
        _displacementMagnitudeProperty.reset();
        return;
    }

    _displacementProperty = std::make_shared<ParticleProperty>(n, ParticleProperty::DisplacementProperty, 3, true);
    _displacementMagnitudeProperty = std::make_shared<ParticleProperty>(n, ParticleProperty::DisplacementMagnitudeProperty, 1, true);

    const bool slipVectors = _computeSlipVector && _slipCutoff > 0.0;
    ProgressSteps steps({
        {"Matching atoms to the reference", 1.0},
        {"Computing displacements", 1.0},
        {"Computing slip vectors", 2.0, slipVectors},
    });
    const ReferenceMatch match = ReferenceMatch::build(
        _identifiers, _refIdentifiers, n, _refPositions->size(), ReferenceMatch::MissingAtoms::Rejected
    );
    log.dump(match);
    const std::vector<std::size_t>& currentToRef = match.currentToReference;
    const std::vector<std::size_t>& refToCurrent = match.referenceToCurrent;
    steps.next();

    const Point3* pos = _positions->constDataPoint3();
    const Point3* refPos = _refPositions->constDataPoint3();

    Vector3* outU = _displacementProperty->dataVector3();
    double* outMag = _displacementMagnitudeProperty->dataDouble();

    if(!pos || !refPos || !outU || !outMag){
        throw std::runtime_error("DisplacementsEngine: null data pointers.");
    }

    const auto pbcFlags = _simCellRef.pbcFlags();
    const auto refCellMatrix = _simCellRef.matrix();

    auto doRange = [&](std::size_t begin, std::size_t end){
        if(_affineMapping != AffineMapping::NoMapping){
            const AffineTransformation reducedToAbsolute =
                (_affineMapping == AffineMapping::ToReferenceCell) ? _simCellRef.matrix() : _simCell.matrix();
         
            const AffineTransformation currRecip = _simCell.inverseMatrix();
            const AffineTransformation refRecip = _simCellRef.inverseMatrix();

            for(std::size_t i = begin; i < end; i++){
                const std::size_t j = currentToRef[i];

                Point3 reducedCurr = currRecip * pos[i];
                Point3 reducedRef = refRecip * refPos[j];

                Vector3 delta = reducedCurr - reducedRef;

                if(_useMinimumImageConvention){
                    for(std::size_t k = 0; k < 3; k++){
                        if(pbcFlags[k]){
                            delta[k] -= std::floor(delta[k] + 0.5);
                        }
                    }
                }

                Vector3 u = reducedToAbsolute * delta;
                outU[i] = u;
                outMag[i] = u.length();
            }
        }else{
            for(std::size_t i = begin; i < end; i++){
                const std::size_t j = currentToRef[i];
                Vector3 u = pos[i] - refPos[j];

                if(_useMinimumImageConvention){
                    for(std::size_t k = 0; k < 3; k++){
                        if(pbcFlags[k]){
                            const Vector3 ck = refCellMatrix.column(k);

                            while((u + ck).squaredLength() < u.squaredLength()){
                                u += ck;
                            }

                            while((u - ck).squaredLength() < u.squaredLength()){
                                u -= ck;
                            }
                        }
                    }
                }

                outU[i] = u;
                outMag[i] = u.length();
            }
        }
    };

    tbb::parallel_for(
        tbb::blocked_range<std::size_t>(0, n, 1024),
        [&](const tbb::blocked_range<std::size_t>& r){
            doRange(r.begin(), r.end());
        }
    );
    log.dump(*this);
    steps.next();

    if(slipVectors){
        computeSlipVectors(currentToRef, refToCurrent);
        log.dump(*_slipVectorProperty);
    }
}

void DisplacementsEngine::computeSlipVectors(
    const std::vector<std::size_t>& currentToRefIndexMap,
    const std::vector<std::size_t>& refToCurrentIndexMap
){
    const std::size_t n = _positions->size();
    _slipVectorProperty = std::make_shared<ParticleProperty>(n, ParticleProperty::DisplacementProperty, 3, true);
    _slipVectorMagnitudeProperty = std::make_shared<ParticleProperty>(n, ParticleProperty::DisplacementMagnitudeProperty, 1, true);

    CutoffNeighborFinder neighborFinder;
    if(!neighborFinder.prepare(_slipCutoff, _refPositions, _simCellRef)){
        throw std::runtime_error("DisplacementsEngine: failed to prepare the slip-vector neighbor list.");
    }

    const Vector3* displacements = _displacementProperty->dataVector3();
    Vector3* outSlip = _slipVectorProperty->dataVector3();
    double* outSlipMag = _slipVectorMagnitudeProperty->dataDouble();
    const double thresholdSquared = _slipThreshold * _slipThreshold;

    ProgressReporter progress("Computing slip vectors", n);
    tbb::parallel_for(
        tbb::blocked_range<std::size_t>(0, n, 256),
        [&](const tbb::blocked_range<std::size_t>& r){
            for(std::size_t i = r.begin(); i < r.end(); i++){
                const std::size_t refIndex = currentToRefIndexMap[i];
                Vector3 accumulated(0.0, 0.0, 0.0);
                int slippedNeighbors = 0;

                for(CutoffNeighborFinder::Query neighborQuery(neighborFinder, refIndex);
                    !neighborQuery.atEnd(); neighborQuery.next()){
                    const std::size_t neighborCurrentIndex = refToCurrentIndexMap[neighborQuery.current()];
                    if(neighborCurrentIndex == ReferenceMatch::kNoCounterpart) continue;

                    const Vector3 relativeSlip = displacements[neighborCurrentIndex] - displacements[i];
                    if(relativeSlip.squaredLength() <= thresholdSquared) continue;

                    accumulated += relativeSlip;
                    slippedNeighbors++;
                }

                if(slippedNeighbors > 0){
                    const Vector3 slip = accumulated * (-1.0 / slippedNeighbors);
                    outSlip[i] = slip;
                    outSlipMag[i] = slip.length();
                }
            }
            progress.advance(r.size());
        }
    );
}

}
