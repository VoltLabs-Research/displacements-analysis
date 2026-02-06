#pragma once

#include <volt/core/volt.h>
#include <volt/core/lammps_parser.h>
#include <volt/core/particle_property.h>
#include <volt/displacements_engine.h>
#include <nlohmann/json.hpp>
#include <memory>
#include <string>

namespace Volt{
using json = nlohmann::json;
 
class DisplacementsService{
public:
    DisplacementsService();

    void setReferenceFrame(const LammpsParser::Frame &ref);
    void setOptions(bool useMinimumImageConvention, DisplacementsEngine::AffineMappingType affineMapping);

    json compute(const LammpsParser::Frame &currentFrame, const std::string &outputFilename);

private:
    bool _hasReference;
    LammpsParser::Frame _referenceFrame;

    bool _useMinimumImageConvention;
    DisplacementsEngine::AffineMappingType _affineMapping;
};

}
