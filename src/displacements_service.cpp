#include <volt/displacements_service.h>
#include <volt/displacements_engine.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/utilities/msgpack_atom_writer.h>
#include <spdlog/spdlog.h>

#include <fstream>
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

    if(!outputFilename.empty()){
        // _displacements.msgpack: summary only (no per-atom array)
        const std::string outputPath = outputFilename + "_displacements.msgpack";
        if(JsonUtils::writeJsonMsgpackToFile(result, outputPath, false)){
            spdlog::info("Displacements msgpack written to {}", outputPath);
        }else{
            spdlog::warn("Could not write displacements msgpack: {}", outputPath);
        }

        // _atoms.msgpack: streaming, no DOM
        auto fieldWriter = [&](MsgpackWriter& w, std::size_t i, int& count){
            count = 2;
            w.write_key("displacement"); w.write_array_header(3);
            const Vector3 u = U ? U->dataVector3()[i] : Vector3(0.0, 0.0, 0.0);
            w.write_double(u.x()); w.write_double(u.y()); w.write_double(u.z());
            w.write_key("magnitude"); w.write_double(Umag ? Umag->getDouble(i) : 0.0);
        };

        const std::string atomsPath = outputFilename + "_atoms.msgpack";
        streamAtomsToFile(atomsPath, currentFrame,
            [](std::size_t){ return std::string("All"); },
            fieldWriter
        );
        spdlog::info("Exported atoms data to: {}", atomsPath);
    }

    return result;
}

}
