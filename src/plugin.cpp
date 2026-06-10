#include <volt/plugin/plugin_main.h>
#include <volt/plugin/output_serializer.h>
#include <volt/displacements_engine.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>

#include <limits>

namespace Volt {

using namespace Particles;

static DisplacementsEngine::AffineMappingType parseAffineMapping(const std::string& s) {
    if (s == "toReferenceCell") return DisplacementsEngine::AffineMappingType::ToReferenceCell;
    if (s == "toCurrentCell") return DisplacementsEngine::AffineMappingType::ToCurrentCell;
    return DisplacementsEngine::AffineMappingType::NoMapping;
}

} // namespace Volt

static const Volt::Plugin::PluginDescriptor descriptor{
    .name = "volt-displacements",
    .description = "Displacements Analysis",
    .options = {
        {"--mic", "bool", "Use minimum image convention", "true"},
        {"--affine_mapping", "string", "noMapping|toReferenceCell|toCurrentCell", "noMapping"},
    },
    .needsReferenceFrame = true
};

VOLT_PLUGIN_MAIN(descriptor,
    [](const auto& opts, const Volt::LammpsParser::Frame& frame,
       const Volt::LammpsParser::Frame* refFramePtr,
       const std::string& outputBase) -> Volt::Plugin::json {
        using namespace Volt;
        using namespace Volt::Particles;

        if (frame.natoms <= 0)
            return AnalysisResult::failure("Invalid number of atoms");

        const LammpsParser::Frame& refFrame = refFramePtr ? *refFramePtr : frame;

        if (frame.natoms != refFrame.natoms)
            return AnalysisResult::failure("Atom count mismatch between current and reference frames");

        auto positions = FrameAdapter::createPositionPropertyShared(frame);
        if (!positions) return AnalysisResult::failure("Failed to create position property");

        auto refPositions = FrameAdapter::createPositionPropertyShared(refFrame);
        if (!refPositions) return AnalysisResult::failure("Failed to create reference position property");

        auto identifiers = FrameAdapter::createIdentifierProperty(frame);
        auto refIdentifiers = FrameAdapter::createIdentifierProperty(refFrame);
        if (!identifiers || !refIdentifiers)
            return AnalysisResult::failure("Failed to create identifier properties");

        bool mic = CLI::getBool(opts, "--mic", true);
        auto affineMapping = parseAffineMapping(CLI::getString(opts, "--affine_mapping", "noMapping"));

        DisplacementsEngine engine(
            positions.get(), frame.simulationCell,
            refPositions.get(), refFrame.simulationCell,
            identifiers.get(), refIdentifiers.get(),
            mic, affineMapping
        );
        engine.perform();

        auto U = engine.displacements();
        auto Umag = engine.displacementMagnitudes();

        const std::size_t n = static_cast<std::size_t>(frame.natoms);
        double totalMag = 0.0, maxMag = 0.0;
        double minMag = std::numeric_limits<double>::max();

        if (U && Umag) {
            for (std::size_t i = 0; i < n; ++i) {
                double m = Umag->getDouble(i);
                totalMag += m;
                if (m > maxMag) maxMag = m;
                if (m < minMag) minMag = m;
            }
        }

        nlohmann::json result;
        result["main_listing"] = {
            {"average_displacement_magnitude", n > 0 ? totalMag / static_cast<double>(n) : 0.0},
            {"max_displacement_magnitude", maxMag},
            {"min_displacement_magnitude", (minMag == std::numeric_limits<double>::max()) ? 0.0 : minMag}
        };

        if (!outputBase.empty()) {
            Plugin::serializePluginOutput(outputBase, frame, result, {
                .summaryFileSuffix = "_displacements",
                .bucketResolver = [](std::size_t) { return std::string("All"); },
                .perAtomColumnWriter = [&U, &Umag](ColumnarAtomWriter& w, std::size_t i) {
                    const Vector3 u = U ? U->dataVector3()[i] : Vector3(0.0, 0.0, 0.0);
                    w.field("displacement", std::vector<double>{u.x(), u.y(), u.z()});
                    w.field("magnitude", Umag ? Umag->getDouble(i) : 0.0);
                }
            });
        }

        return result;
    })
