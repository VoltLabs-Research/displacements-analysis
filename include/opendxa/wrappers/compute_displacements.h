#pragma once

#include <opendxa/core/opendxa.h>
#include <opendxa/core/lammps_parser.h>
#include <opendxa/core/particle_property.h>
#include <opendxa/analysis/compute_displacements.h>
#include <nlohmann/json.hpp>
#include <memory>
#include <string>

namespace OpenDXA{
using json = nlohmann::json;
 
class DisplacementsWrapper{
public:
    DisplacementsWrapper();

    void setReferenceFrame(const LammpsParser::Frame &ref);
    void setOptions(bool useMinimumImageConvention, ComputeDisplacements::AffineMappingType affineMapping);

    json compute(const LammpsParser::Frame &currentFrame, const std::string &outputFilename);

private:
    std::shared_ptr<Particles::ParticleProperty> createPositionProperty(const LammpsParser::Frame& frame);
    std::shared_ptr<Particles::ParticleProperty> createIdentifierProperty(const LammpsParser::Frame& frame);

    bool _hasReference;
    LammpsParser::Frame _referenceFrame;

    bool _useMinimumImageConvention;
    ComputeDisplacements::AffineMappingType _affineMapping;
};

}