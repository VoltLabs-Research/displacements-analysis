#pragma once

#include <volt/core/volt.h>
#include <volt/core/simulation_cell.h>
#include <volt/core/particle_property.h>

#include <memory>

namespace Volt{

class DisplacementsEngine{
public:
    enum class AffineMappingType{
        NoMapping = 0,
        ToReferenceCell,
        ToCurrentCell
    };

    DisplacementsEngine(
        Particles::ParticleProperty* positions,
        const SimulationCell& cell,
        Particles::ParticleProperty* refPositions,
        const SimulationCell& refCell,
        Particles::ParticleProperty* identifiers,
        Particles::ParticleProperty* refIdentifiers,
        bool useMinimumImageConvention = true,
        AffineMappingType affineMapping = AffineMappingType::NoMapping
    );

    void perform();

    // Zimmerman slip vector: for each atom, the negative average relative
    // displacement of reference-frame neighbors whose displacement difference
    // exceeds the threshold. Disabled unless configured before perform().
    void setSlipVectorOptions(bool compute, double cutoff, double threshold){
        _computeSlipVector = compute;
        _slipCutoff = cutoff;
        _slipThreshold = threshold;
    }

    std::shared_ptr<Particles::ParticleProperty> displacements() const{
        return _displacementProperty;
    }

    std::shared_ptr<Particles::ParticleProperty> displacementMagnitudes() const{
        return _displacementMagnitudeProperty;
    }

    std::shared_ptr<Particles::ParticleProperty> slipVectors() const{
        return _slipVectorProperty;
    }

    std::shared_ptr<Particles::ParticleProperty> slipVectorMagnitudes() const{
        return _slipVectorMagnitudeProperty;
    }

private:
    void buildParticleMapping(
        std::vector<std::size_t>& currentToRefIndexMap,
        std::vector<std::size_t>& refToCurrentIndexMap,
        bool requireCompleteCurrentToRefMapping,
        bool requireCompleteRefToCurrentMapping
    ) const;

    void computeSlipVectors(
        const std::vector<std::size_t>& currentToRefIndexMap,
        const std::vector<std::size_t>& refToCurrentIndexMap
    );

    Particles::ParticleProperty* _positions = nullptr;
    Particles::ParticleProperty* _refPositions = nullptr;
    Particles::ParticleProperty* _identifiers = nullptr;
    Particles::ParticleProperty* _refIdentifiers = nullptr;

    SimulationCell _simCell;
    SimulationCell _simCellRef;

    bool _useMinimumImageConvention = true;
    AffineMappingType _affineMapping = AffineMappingType::NoMapping;

    bool _computeSlipVector = false;
    double _slipCutoff = 0.0;
    double _slipThreshold = 0.5;

    std::shared_ptr<Particles::ParticleProperty> _displacementProperty;
    std::shared_ptr<Particles::ParticleProperty> _displacementMagnitudeProperty;
    std::shared_ptr<Particles::ParticleProperty> _slipVectorProperty;
    std::shared_ptr<Particles::ParticleProperty> _slipVectorMagnitudeProperty;
};

}
