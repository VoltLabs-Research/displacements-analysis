#pragma once

#include <volt/core/volt.h>
#include <volt/core/simulation_cell.h>
#include <volt/core/particle_property.h>
#include <volt/analysis/reference_configuration.h>

#include <cstddef>
#include <limits>
#include <memory>
#include <vector>

namespace Volt{

class DisplacementsLog;

class DisplacementsEngine{
public:
    // The displacement and slip vector magnitudes over the whole frame, for the summary table and the log.
    struct Summary{
        std::size_t atoms = 0;
        bool hasDisplacements = false;
        double totalMagnitude = 0.0;
        double maxMagnitude = 0.0;
        double minMagnitude = std::numeric_limits<double>::max();
        bool hasSlipVectors = false;
        double totalSlipMagnitude = 0.0;
        double maxSlipMagnitude = 0.0;
        std::size_t slippedAtoms = 0;
    };

    DisplacementsEngine(
        Particles::ParticleProperty* positions,
        const SimulationCell& cell,
        Particles::ParticleProperty* refPositions,
        const SimulationCell& refCell,
        Particles::ParticleProperty* identifiers,
        Particles::ParticleProperty* refIdentifiers,
        bool useMinimumImageConvention = true,
        AffineMapping affineMapping = AffineMapping::NoMapping
    );

    void perform(DisplacementsLog& log);

    void setSlipVectorOptions(bool compute, double cutoff, double threshold){
        _computeSlipVector = compute;
        _slipCutoff = cutoff;
        _slipThreshold = threshold;
    }

    bool useMinimumImageConvention() const{
        return _useMinimumImageConvention;
    }

    AffineMapping affineMapping() const{
        return _affineMapping;
    }

    const SimulationCell& referenceCell() const{
        return _simCellRef;
    }

    bool computesSlipVector() const{
        return _computeSlipVector;
    }

    double slipCutoff() const{
        return _slipCutoff;
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
    AffineMapping _affineMapping = AffineMapping::NoMapping;

    bool _computeSlipVector = false;
    double _slipCutoff = 0.0;
    double _slipThreshold = 0.5;

    std::shared_ptr<Particles::ParticleProperty> _displacementProperty;
    std::shared_ptr<Particles::ParticleProperty> _displacementMagnitudeProperty;
    std::shared_ptr<Particles::ParticleProperty> _slipVectorProperty;
    std::shared_ptr<Particles::ParticleProperty> _slipVectorMagnitudeProperty;
};

}
