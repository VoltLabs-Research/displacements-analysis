#include <volt/plugin/plugin_main.h>
#include <volt/plugin/output_serializer.h>
#include <volt/displacements_engine.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/structures/crystal_structure_types.h>

#include <limits>

namespace Volt {

using namespace Particles;

static DisplacementsEngine::AffineMappingType parseAffineMapping(const std::string& s) {
    if (s == "toReferenceCell") return DisplacementsEngine::AffineMappingType::ToReferenceCell;
    if (s == "toCurrentCell") return DisplacementsEngine::AffineMappingType::ToCurrentCell;
    if (s != "noMapping" && s != "none")
        spdlog::warn("Unknown affine_mapping '{}' (expected noMapping, toReferenceCell or toCurrentCell); no affine mapping applied", s);
    return DisplacementsEngine::AffineMappingType::NoMapping;
}

}

static const Volt::Plugin::PluginDescriptor descriptor{
    .name = "volt-displacements",
    .description = "Displacements Analysis",
    .options = {
        {"--mic", "bool", "Use minimum image convention", "true"},
        {"--affine_mapping", "string", "noMapping|toReferenceCell|toCurrentCell", "noMapping"},
        {"--compute_slip_vector", "bool", "Compute the Zimmerman slip vector per atom", "true"},
        {"--slip_cutoff", "double", "Reference-frame neighbor cutoff distance for the slip vector", "3.5"},
        {"--slip_threshold", "double", "Relative displacement above which a neighbor counts as slipped", "0.5"},
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
        if (!refFramePtr)
            spdlog::warn("No reference frame given (--reference); the current frame is its own reference, so every displacement is zero");

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
        engine.setSlipVectorOptions(
            CLI::getBool(opts, "--compute_slip_vector", true),
            CLI::getDouble(opts, "--slip_cutoff", 3.5),
            CLI::getDouble(opts, "--slip_threshold", 0.5)
        );
        engine.perform();

        auto U = engine.displacements();
        auto Umag = engine.displacementMagnitudes();
        auto S = engine.slipVectors();
        auto Smag = engine.slipVectorMagnitudes();

        const std::size_t n = static_cast<std::size_t>(frame.natoms);
        double totalMag = 0.0, maxMag = 0.0;
        double minMag = std::numeric_limits<double>::max();
        double totalSlipMag = 0.0, maxSlipMag = 0.0;

        if (U && Umag) {
            for (std::size_t i = 0; i < n; ++i) {
                double m = Umag->getDouble(i);
                totalMag += m;
                if (m > maxMag) maxMag = m;
                if (m < minMag) minMag = m;
            }
            spdlog::info("Displacement magnitude: min {:.4f}, mean {:.4f}, max {:.4f}",
                minMag, totalMag / static_cast<double>(n), maxMag);
        }

        if (Smag) {
            std::size_t slippedAtoms = 0;
            for (std::size_t i = 0; i < n; ++i) {
                double m = Smag->getDouble(i);
                totalSlipMag += m;
                if (m > maxSlipMag) maxSlipMag = m;
                if (m > 0.0) ++slippedAtoms;
            }
            spdlog::info("Slip vector magnitude: mean {:.4f}, max {:.4f}; {} atoms with nonzero slip ({:.1f}%)",
                totalSlipMag / static_cast<double>(n), maxSlipMag, slippedAtoms, percentOf(slippedAtoms, n));
        }

        nlohmann::json result;
        result["main_listing"] = {
            {"average_displacement_magnitude", n > 0 ? totalMag / static_cast<double>(n) : 0.0},
            {"max_displacement_magnitude", maxMag},
            {"min_displacement_magnitude", (minMag == std::numeric_limits<double>::max()) ? 0.0 : minMag}
        };
        if (Smag) {
            result["main_listing"]["average_slip_vector_magnitude"] = n > 0 ? totalSlipMag / static_cast<double>(n) : 0.0;
            result["main_listing"]["max_slip_vector_magnitude"] = maxSlipMag;
        }

        if (!outputBase.empty()) {
            const auto* structureCol = frame.findAtomProperty("structure_id");
            if (!structureCol) structureCol = frame.findAtomProperty("structure_type");
            const auto structureTypeAt = [structureCol](std::size_t i) -> int {
                return (structureCol && i < structureCol->size()) ? structureCol->ints[i] : -1;
            };

            Plugin::serializePluginOutput(outputBase, frame, result, {
                .summaryFileSuffix = "_displacements",
                .bucketResolver = [&](std::size_t i) {
                    const int s = structureTypeAt(i);
                    return s < 0 ? std::string("All") : std::string(structureTypeName(s));
                },
                .perAtomColumnWriter = [&U, &Umag, &S, &Smag](ColumnarAtomWriter& w, std::size_t i) {
                    const Vector3 u = U ? U->dataVector3()[i] : Vector3(0.0, 0.0, 0.0);
                    w.field("displacement", std::vector<double>{u.x(), u.y(), u.z()});
                    w.field("magnitude", Umag ? Umag->getDouble(i) : 0.0);
                    if (S && Smag) {
                        const Vector3 s = S->dataVector3()[i];
                        w.field("slip_vector", std::vector<double>{s.x(), s.y(), s.z()});
                        w.field("slip_vector_magnitude", Smag->getDouble(i));
                    }
                },
                .resolveStructureId = structureCol ? StructureIdResolver(structureTypeAt) : StructureIdResolver{}
            });
        }

        return result;
    })
