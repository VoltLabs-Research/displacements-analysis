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

    std::shared_ptr<Particles::ParticleProperty> displacements() const{
        return _displacementProperty;
    }

    std::shared_ptr<Particles::ParticleProperty> displacementMagnitudes() const{
        return _displacementMagnitudeProperty;
    }

private:
    void buildParticleMapping(
        std::vector<std::size_t>& currentToRefIndexMap,
        std::vector<std::size_t>& refToCurrentIndexMap,
        bool requireCompleteCurrentToRefMapping,
        bool requireCompleteRefToCurrentMapping
    ) const;

    Particles::ParticleProperty* _positions = nullptr;
    Particles::ParticleProperty* _refPositions = nullptr;
    Particles::ParticleProperty* _identifiers = nullptr;
    Particles::ParticleProperty* _refIdentifiers = nullptr;

    SimulationCell _simCell;
    SimulationCell _simCellRef;

    bool _useMinimumImageConvention = true;
    AffineMappingType _affineMapping = AffineMappingType::NoMapping;

    std::shared_ptr<Particles::ParticleProperty> _displacementProperty;
    std::shared_ptr<Particles::ParticleProperty> _displacementMagnitudeProperty;
};

}
