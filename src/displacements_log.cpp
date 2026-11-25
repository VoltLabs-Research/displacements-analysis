#include <volt/displacements_log.h>

#include <spdlog/spdlog.h>

#include <string>

namespace Volt{

using namespace Particles;

void DisplacementsLog::dump(const DisplacementsEngine& engine){
    std::string minimumImage;
    if(engine.useMinimumImageConvention()){
        const auto pbcFlags = engine.referenceCell().pbcFlags();
        for(std::size_t k = 0; k < 3; k++){
            if(!pbcFlags[k]) continue;
            if(!minimumImage.empty()) minimumImage += ' ';
            minimumImage += "xyz"[k];
        }
        if(minimumImage.empty()) minimumImage = "no periodic axis";
    }else{
        minimumImage = "off";
    }
    const auto affineMapping = engine.affineMapping();
    const char* mapping = affineMapping == AffineMapping::ToReferenceCell ? "to reference cell"
        : affineMapping == AffineMapping::ToCurrentCell ? "to current cell" : "none";
    spdlog::info(
        "Displacements computed for {} atoms in {} (minimum image convention: {}; affine mapping: {})",
        engine.displacements()->size(), phaseElapsed(), minimumImage, mapping
    );
    if(engine.computesSlipVector() && !(engine.slipCutoff() > 0.0)){
        spdlog::warn("Slip vectors skipped: slip_cutoff must be positive (got {})", engine.slipCutoff());
    }
}

void DisplacementsLog::dump(const ParticleProperty&){
    spdlog::info("Slip vectors computed in {}", phaseElapsed());
}

void DisplacementsLog::dump(const DisplacementsEngine::Summary& summary){
    const double atoms = static_cast<double>(summary.atoms);
    if(summary.hasDisplacements){
        spdlog::info("Displacement magnitude: min {:.4f}, mean {:.4f}, max {:.4f}",
            summary.minMagnitude, summary.totalMagnitude / atoms, summary.maxMagnitude);
    }
    if(summary.hasSlipVectors){
        spdlog::info("Slip vector magnitude: mean {:.4f}, max {:.4f}; {} atoms with nonzero slip ({:.1f}%)",
            summary.totalSlipMagnitude / atoms, summary.maxSlipMagnitude, summary.slippedAtoms,
            percentOf(summary.slippedAtoms, summary.atoms));
    }
}

}
