#include <volt/displacements_engine.h>
#include <volt/analysis/cutoff_neighbor_finder.h>

#include <stdexcept>
#include <limits>
#include <unordered_map>
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
    AffineMappingType affineMapping
)
    : _positions(positions),
      _refPositions(refPositions),
      _identifiers(identifiers),
      _refIdentifiers(refIdentifiers),
      _simCell(cell),
      _simCellRef(refCell),
      _useMinimumImageConvention(useMinimumImageConvention),
      _affineMapping(affineMapping){}

void DisplacementsEngine::buildParticleMapping(
    std::vector<std::size_t>& currentToRefIndexMap,
    std::vector<std::size_t>& refToCurrentIndexMap,
    bool requireCompleteCurrentToRefMapping,
    bool requireCompleteRefToCurrentMapping
) const{
    const std::size_t nCurr = _positions->size();
    const std::size_t nRef = _refPositions->size();

    currentToRefIndexMap.resize(nCurr);
    refToCurrentIndexMap.resize(nRef);

    const bool haveIds = (_identifiers && _refIdentifiers);
    
    if(haveIds){
        if(_identifiers->size() != nCurr){
            throw std::runtime_error("DisplacementsEngine: identifiers size != positions size.");
        }

        if(_refIdentifiers->size() != nRef){
            throw std::runtime_error("DisplacementsEngine: refIdentifiers size != refPositions size.");
        }

        // Map ID -> index (reference)
        std::unordered_map<int, std::size_t> refMap;
        refMap.reserve(nRef * 2);

        for(std::size_t i = 0; i < nRef; i++){
            const int id = _refIdentifiers->getInt(i);
            auto [it, inserted] = refMap.emplace(id, i);
            if(!inserted){
                throw std::runtime_error("DisplacementsEngine: duplicate particle identifier in reference configuration.");
            }
        }

        // Check duplicates in current + build currentMap
        std::unordered_map<int, std::size_t> currMap;
        currMap.reserve(nCurr * 2);

        for(std::size_t i = 0; i < nCurr; i++){
            const int id = _identifiers->getInt(i);
            auto [it, inserted] = currMap.emplace(id, i);
            if(!inserted) {
                throw std::runtime_error("DisplacementsEngine: duplicate particle identifier in current configuration.");
            }
        }

        // Build current -> ref
        for(std::size_t i = 0; i < nCurr; i++){
            const int id = _identifiers->getInt(i);
            auto it = refMap.find(id);
            if(it != refMap.end()){
                currentToRefIndexMap[i] = it->second;
            }else if(requireCompleteCurrentToRefMapping){
                throw std::runtime_error("DisplacementsEngine: particle ID exists in current but not in reference.");
            }else{
                currentToRefIndexMap[i] = std::numeric_limits<std::size_t>::max();
            }
        }

        // Build ref -> current
        for(std::size_t i = 0; i < nRef; i++){
            const int id = _refIdentifiers->getInt(i);
            auto it = currMap.find(id);
            if(it != currMap.end()){
                refToCurrentIndexMap[i] = it->second;
            }else if(requireCompleteRefToCurrentMapping){
                throw std::runtime_error("DisplacementsEngine: particle ID exists in reference but not in current.");
            }else{
                // OJO: aca estabas escribiendo en el vector equivocado
                refToCurrentIndexMap[i] = std::numeric_limits<std::size_t>::max();
            }
        }
    }else{
        // NO IDs: require same number of particles and assume same ordering
        if(nCurr != nRef){
            throw std::runtime_error("DisplacementsEngine: positions and refPositions size mismatch and no identifiers present.");
        }

        for(std::size_t i = 0; i < nCurr; i++){
            currentToRefIndexMap[i] = i;
            refToCurrentIndexMap[i] = i;
        }
    }
}

void DisplacementsEngine::perform(){
    if(!_positions || !_refPositions){
        throw std::runtime_error("DisplacementsEngine: null input properties.");
    }

    const std::size_t n = _positions->size();
    if(n == 0){
        _displacementProperty.reset();
        _displacementMagnitudeProperty.reset();
        return;
    }

    // Allocate output properties
    _displacementProperty = std::make_shared<ParticleProperty>(n, ParticleProperty::DisplacementProperty, 3, true);
    _displacementMagnitudeProperty = std::make_shared<ParticleProperty>(n, ParticleProperty::DisplacementMagnitudeProperty, 1, true);

    // Build mapping
    std::vector<std::size_t> currentToRef;
    std::vector<std::size_t> refToCurrent;
    buildParticleMapping(currentToRef, refToCurrent, /*requireCompleteCurrentToRefMapping=*/true, /*requireCompleteRefToCurrentMapping=*/false);

    // Accessors
    const Point3* pos = _positions->constDataPoint3();
    const Point3* refPos = _refPositions->constDataPoint3();

    Vector3* outU = _displacementProperty->dataVector3();
    double* outMag = _displacementMagnitudeProperty->dataDouble();

    if(!pos || !refPos || !outU || !outMag){
        throw std::runtime_error("DisplacementsEngine: null data pointers.");
    }

    // PBC and cell matrices
    const auto pbcFlags = _simCellRef.pbcFlags();
    const auto refCellMatrix = _simCellRef.matrix();

    auto doRange = [&](std::size_t begin, std::size_t end){
        if(_affineMapping != AffineMappingType::NoMapping){
            // Reduced coordinate delta + MIC + Mapping to absolute
            const AffineTransformation reducedToAbsolute =
                (_affineMapping == AffineMappingType::ToReferenceCell) ? _simCellRef.matrix() : _simCell.matrix();
         
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
            // Direct displacement in absolute coords + MIC by cell columns
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

    if(_computeSlipVector && _slipCutoff > 0.0){
        computeSlipVectors(currentToRef, refToCurrent);
    }
}

void DisplacementsEngine::computeSlipVectors(
    const std::vector<std::size_t>& currentToRefIndexMap,
    const std::vector<std::size_t>& refToCurrentIndexMap
){
    const std::size_t n = _positions->size();
    _slipVectorProperty = std::make_shared<ParticleProperty>(n, ParticleProperty::DisplacementProperty, 3, true);
    _slipVectorMagnitudeProperty = std::make_shared<ParticleProperty>(n, ParticleProperty::DisplacementMagnitudeProperty, 1, true);

    // Neighbors come from the reference configuration: the slip vector compares
    // each bond against its pre-deformation counterpart, so the bond list must
    // predate the deformation as well.
    CutoffNeighborFinder neighborFinder;
    if(!neighborFinder.prepare(_slipCutoff, _refPositions, _simCellRef)){
        throw std::runtime_error("DisplacementsEngine: failed to prepare the slip-vector neighbor list.");
    }

    const Vector3* displacements = _displacementProperty->dataVector3();
    Vector3* outSlip = _slipVectorProperty->dataVector3();
    double* outSlipMag = _slipVectorMagnitudeProperty->dataDouble();
    const double thresholdSquared = _slipThreshold * _slipThreshold;
    constexpr std::size_t unmapped = std::numeric_limits<std::size_t>::max();

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
                    if(neighborCurrentIndex == unmapped) continue;

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
        }
    );
}

}
