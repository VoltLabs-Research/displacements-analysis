#include <volt/displacements_service.h>
#include <volt/displacements_engine.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <spdlog/spdlog.h>

namespace Volt{

using namespace Volt::Particles;

DisplacementsService::DisplacementsService()
    : _hasReference(false),
      _useMinimumImageConvention(true),
      _affineMapping(DisplacementsEngine::AffineMappingType::NoMapping){}

void DisplacementsService::setReferenceFrame(const LammpsParser::Frame &ref){
    _referenceFrame = ref;
    _hasReference = true;
}

void DisplacementsService::setOptions(bool useMinimumImageConvention, DisplacementsEngine::AffineMappingType affineMapping){
    _useMinimumImageConvention = useMinimumImageConvention;
    _affineMapping = affineMapping;
}

json DisplacementsService::compute(const LammpsParser::Frame& currentFrame, const std::string &outputFilename){
    auto startTime = std::chrono::high_resolution_clock::now();

    if(currentFrame.natoms <= 0){
        return AnalysisResult::failure("Invalid number of atoms");
    }

    const LammpsParser::Frame &refFrame = _hasReference ? _referenceFrame : currentFrame;

    if(currentFrame.natoms != refFrame.natoms){
        return AnalysisResult::failure("Atom count mismatch between current and reference frames");
    }

    auto positions = FrameAdapter::createPositionPropertyShared(currentFrame);
    if(!positions){
        return AnalysisResult::failure("Failed to create position property");
    }

    auto refPositions = FrameAdapter::createPositionPropertyShared(refFrame);
    if(!refPositions){
        return AnalysisResult::failure("Failed to create reference position property");
    }

    auto identifiers = FrameAdapter::createIdentifierProperty(currentFrame);
    auto refIdentifiers = FrameAdapter::createIdentifierProperty(refFrame);
    if(!identifiers || !refIdentifiers){
        return AnalysisResult::failure("Failed to create identifier properties");
    }

    spdlog::info("Starting displacement analysis (MIC = {}, affineMapping = {})...",
        _useMinimumImageConvention ? "true" : "false",
        static_cast<int>(_affineMapping)
    );

    DisplacementsEngine engine(
        positions.get(),
        currentFrame.simulationCell,
        refPositions.get(),
        refFrame.simulationCell,
        identifiers.get(),
        refIdentifiers.get(),
        _useMinimumImageConvention,
        _affineMapping
    );

    engine.perform();

    json result = AnalysisResult::success();
    AnalysisResult::addTiming(result, startTime);

    if(!outputFilename.empty()){
        // TODO: Implement msgpack export in standalone package
        spdlog::warn("File output not yet implemented in standalone package");
    }
    result["displacements"] = json::array();

    return result;
}

}
