#include <volt/displacements_service.h>
#include <volt/displacements_engine.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/utilities/json_utils.h>
#include <spdlog/spdlog.h>
#include <limits>

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

    if(currentFrame.natoms <= 0)
        return AnalysisResult::failure("Invalid number of atoms");

    const LammpsParser::Frame &refFrame = _hasReference ? _referenceFrame : currentFrame;

    if(currentFrame.natoms != refFrame.natoms)
        return AnalysisResult::failure("Atom count mismatch between current and reference frames");

    auto positions = FrameAdapter::createPositionPropertyShared(currentFrame);
    if(!positions)
        return AnalysisResult::failure("Failed to create position property");

    auto refPositions = FrameAdapter::createPositionPropertyShared(refFrame);
    if(!refPositions)
        return AnalysisResult::failure("Failed to create reference position property");

    auto identifiers = FrameAdapter::createIdentifierProperty(currentFrame);
    auto refIdentifiers = FrameAdapter::createIdentifierProperty(refFrame);
    if(!identifiers || !refIdentifiers)
        return AnalysisResult::failure("Failed to create identifier properties");

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

    auto U = engine.displacements();
    auto Umag = engine.displacementMagnitudes();

    const size_t n = static_cast<size_t>(currentFrame.natoms);
    double totalMag = 0.0;
    double maxMag = 0.0;
    double minMag = std::numeric_limits<double>::max();

    if(U && Umag){
        for(size_t i = 0; i < n; i++){
            double m = Umag->getDouble(i);
            totalMag += m;
            if(m > maxMag) maxMag = m;
            if(m < minMag) minMag = m;
        }
    }

    json result;
    result["main_listing"] = {
        { "average_displacement_magnitude", n > 0 ? totalMag / n : 0.0 },
        { "max_displacement_magnitude", maxMag },
        { "min_displacement_magnitude", (minMag == std::numeric_limits<double>::max()) ? 0.0 : minMag }
    };

    json perAtom = json::array();
    for(size_t i = 0; i < n; i++){
        json a;
        a["id"] = currentFrame.ids[i];

        if(U){
            Vector3 u = U->dataVector3()[i];
            a["displacement"] = {u.x(), u.y(), u.z()};
        } else {
            a["displacement"] = {0.0, 0.0, 0.0};
        }

        a["magnitude"] = Umag ? Umag->getDouble(i) : 0.0;
        perAtom.push_back(a);
    }
    result["per-atom-properties"] = perAtom;

    if(!outputFilename.empty()){
        const std::string outputPath = outputFilename + "_displacements.msgpack";
        if(JsonUtils::writeJsonMsgpackToFile(result, outputPath, false)){
            spdlog::info("Displacements msgpack written to {}", outputPath);
        }else{
            spdlog::warn("Could not write displacements msgpack: {}", outputPath);
        }
    }

    return result;
}

}

