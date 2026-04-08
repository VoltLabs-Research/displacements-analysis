#include <volt/cli/common.h>
#include <volt/displacements_service.h>
#include <oneapi/tbb/global_control.h>

using namespace Volt;
using namespace Volt::CLI;

static void showUsage(const std::string& name){
    printUsageHeader(name, "Volt - Displacements Analysis");
    std::cerr
        << "  --reference <file>            Reference LAMMPS dump file.\n"
        << "                                If omitted, current frame is used (≈ zero displacement).\n"
        << "  --mic                         Use minimum image convention. [default: true]\n"
        << "  --affineMapping <mode>        Affine mapping mode: noMapping|toReferenceCell|toCurrentCell [default: noMapping]\n"
        << "  --threads <int>               Max worker threads (TBB/OMP). [default: auto]\n";
    printHelpOption();
}

static DisplacementsEngine::AffineMappingType parseAffineMapping(const std::string& s){
    if(s == "noMapping") return DisplacementsEngine::AffineMappingType::NoMapping;
    if(s == "toReferenceCell") return DisplacementsEngine::AffineMappingType::ToReferenceCell;
    if(s == "toCurrentCell") return DisplacementsEngine::AffineMappingType::ToCurrentCell;

    spdlog::warn("Unknown affineMapping '{}', defaulting to 'none'.", s);
    return DisplacementsEngine::AffineMappingType::NoMapping;
}

int main(int argc, char* argv[]){
    if(argc < 2){
        showUsage(argv[0]);
        return 1;
    }

    std::string filename, outputBase;
    auto opts = parseArgs(argc, argv, filename, outputBase);

    if(hasOption(opts, "--help") || filename.empty()){
        showUsage(argv[0]);
        return filename.empty() ? 1 : 0;
    }

    const int requestedThreads = std::max(1, getInt(opts, "--threads", std::thread::hardware_concurrency() > 0
        ? static_cast<int>(std::thread::hardware_concurrency())
        : 1));
    oneapi::tbb::global_control parallelControl(
        oneapi::tbb::global_control::max_allowed_parallelism,
        static_cast<std::size_t>(requestedThreads)
    );
    initLogging("volt-displacements");
    spdlog::info("Using {} threads (OneTBB)", requestedThreads);

    LammpsParser::Frame frame;
    if(!parseFrame(filename, frame)) return 1;

    // Parse reference frame if provided
    std::string refFile = getString(opts, "--reference");
    LammpsParser::Frame refFrame;
    bool hasReference = false;

    if(!refFile.empty()){
        spdlog::info("Parsing reference file: {}", refFile);
        LammpsParser refParser;
        if(!refParser.parseFile(refFile, refFrame)){
            spdlog::error("Failed to parse reference file: {}", refFile);
            return 1;
        }
        if(refFrame.natoms != frame.natoms){
            spdlog::error("Atom count mismatch: current={} reference={}", frame.natoms, refFrame.natoms);
            return 1;
        }
        hasReference = true;
        spdlog::info("Reference loaded: {} atoms", refFrame.natoms);
    }

    outputBase = deriveOutputBase(filename, outputBase);
    spdlog::info("Output base: {}", outputBase);

    // Options
    bool mic = getBool(opts, "--mic", true);
    auto affineMapping = parseAffineMapping(getString(opts, "--affineMapping", "none"));

    DisplacementsService analyzer;
    analyzer.setOptions(mic, affineMapping);

    if(hasReference){
        analyzer.setReferenceFrame(refFrame);
    }

    spdlog::info("Starting displacements analysis...");
    json result = analyzer.compute(frame, outputBase);

    if(result.value("is_failed", false)){
        spdlog::error("Analysis failed: {}", result.value("error", "Unknown error"));
        return 1;
    }

    spdlog::info("Displacements analysis completed.");
    return 0;
}
