#pragma once

#include <volt/core/particle_property.h>
#include <volt/core/run_log.h>
#include <volt/displacements_engine.h>

namespace Volt{

// The execution log of one displacements analysis. Each phase hands it what it built, and the log writes what that
// holds and how long the phase took. Scripts parse these lines, so their wording stays as it is.
class DisplacementsLog : public RunLog{
public:
    using RunLog::dump;

    // The displacements, once computed; it also warns when the slip vectors that were asked for are skipped.
    void dump(const DisplacementsEngine& engine);
    void dump(const Particles::ParticleProperty& slipVectors);
    void dump(const DisplacementsEngine::Summary& summary);
};

}
