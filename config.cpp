/***************************************************************************
 * SeisComP spectral moment-magnitude plugin (mwspec)                      *
 *                                                                         *
 * Shared configuration reader for the amplitude and magnitude processors. *
 *                                                                         *
 * GNU Affero General Public License Usage - see LICENSE.                   *
 ***************************************************************************/


#define SEISCOMP_COMPONENT MwSpec

#include <seiscomp/logging/log.h>
#include <seiscomp/config/config.h>

#include "mwspec.h"


namespace Seiscomp {
namespace Magnitudes {
namespace MwSpec {


namespace {

// Reads a list of model layer strings from the global configuration. Each
// list item is "depth vp vs qp qap qs qas density".
bool readModelLayers(const Config::Config *cfg, const std::string &key,
                     std::vector<std::string> &layers) {
	if ( !cfg ) {
		return false;
	}
	try {
		layers = cfg->getStrings(key);
		return !layers.empty();
	}
	catch ( ... ) {}
	return false;
}

}


bool readMwSpecConfig(const Processing::Settings &settings,
                      const std::string &prefix, MwSpecConfig &out) {
	const Config::Config *cfg = settings.localConfiguration;

	// The spectral model layers and the phase are physically shared by the
	// amplitude and magnitude processors, so they are read from one canonical
	// key (under "magnitudes.Mw(spec)") regardless of the caller's prefix. This
	// guarantees both sides use a consistent model and phase.
	const std::string canonical = std::string("magnitudes.") + MWSPEC_TYPE;

	// --- phase (canonical, with per-prefix fallback) ----------------------
	{
		std::string phase;
		if ( (settings.getValue(phase, canonical + ".phase") ||
		      settings.getValue(phase, prefix + ".phase")) && !phase.empty() ) {
			char p = static_cast<char>(::toupper(phase[0]));
			if ( p == 'P' || p == 'S' ) {
				out.phase = p;
			}
			else {
				SEISCOMP_ERROR("%s.phase: must be P or S, got '%s'",
				               prefix.c_str(), phase.c_str());
				return false;
			}
		}
	}

	// --- layered spectral model ------------------------------------------
	{
		std::vector<std::string> layers;
		if ( readModelLayers(cfg, canonical + ".model", layers) ||
		     readModelLayers(cfg, prefix + ".model", layers) ) {
			if ( !out.model.setLayers(layers) ) {
				SEISCOMP_ERROR("%s.model: malformed layer definition", prefix.c_str());
				return false;
			}
		}
		else {
			// No layered model configured: build a single uniform layer from the
			// scalar parameters so the attenuation (Q0, Qalpha) and the source
			// velocity/density are DIRECTLY configurable. Defaults reproduce the
			// previous Seisan single-layer fallback (vp=6, vs=3.5, density=3,
			// Q0=0 i.e. no attenuation correction). For reliable Mw, set Q0/Qalpha
			// to a regionally-calibrated value (the analogue of Seisan CODAQ/QLG).
			double vp = 6.0, vs = 3.5, density = 3.0, q0 = 0.0, qalpha = 0.0;
			settings.getValue(vp,      canonical + ".vp");
			settings.getValue(vs,      canonical + ".vs");
			settings.getValue(density, canonical + ".density");
			settings.getValue(q0,      canonical + ".Q0");
			settings.getValue(qalpha,  canonical + ".Qalpha");
			std::string layer =
				"0.0 " + std::to_string(vp) + " " + std::to_string(vs) + " " +
				std::to_string(q0) + " " + std::to_string(qalpha) + " " +
				std::to_string(q0) + " " + std::to_string(qalpha) + " " +
				std::to_string(density);
			out.model.setLayers({layer});
		}

		double kp = 0.01, ks = 0.02;
		bool haveKp = settings.getValue(kp, prefix + ".kappaP");
		bool haveKs = settings.getValue(ks, prefix + ".kappaS");
		if ( haveKp || haveKs ) {
			out.model.setKappa(kp, ks);
		}

		double q1p = 1.0, q1s = 1.0;
		bool haveQ1p = settings.getValue(q1p, prefix + ".qBelow1HzP");
		bool haveQ1s = settings.getValue(q1s, prefix + ".qBelow1HzS");
		if ( haveQ1p || haveQ1s ) {
			out.model.setQBelow1Hz(q1p, q1s);
		}

		double qCorner = 0.0;
		if ( settings.getValue(qCorner, prefix + ".qCorner") ) {
			out.model.setQCorner(qCorner);
		}
	}

	// --- moment computation ----------------------------------------------
	settings.getValue(out.radiation,        prefix + ".radiation");
	settings.getValue(out.freeSurface,      prefix + ".freeSurface");
	settings.getValue(out.geoDepth1,        prefix + ".geoDepth1");
	settings.getValue(out.geoDepth2,        prefix + ".geoDepth2");
	settings.getValue(out.herkijDistanceKm, prefix + ".herkijDistance");
	settings.getValue(out.minDistanceDeg,   prefix + ".minimumDistance");
	settings.getValue(out.maxDistanceDeg,   prefix + ".maximumDistance");

	// --- spectral measurement --------------------------------------------
	settings.getValue(out.signalPreTime,  prefix + ".signalPreTime");
	settings.getValue(out.signalDuration, prefix + ".signalDuration");
	settings.getValue(out.noiseGap,       prefix + ".noiseGap");
	settings.getValue(out.fixedFmin,      prefix + ".fmin");
	settings.getValue(out.fixedFmax,      prefix + ".fmax");
	settings.getValue(out.nfreq,          prefix + ".nfreq");
	settings.getValue(out.minSNR,         prefix + ".minSNR");
	settings.getValue(out.maxResidual,    prefix + ".maxResidual");
	settings.getValue(out.calibration,    prefix + ".calibration");
	settings.getValue(out.applyTaper,     prefix + ".taper");

	// --- grid search tunables --------------------------------------------
	settings.getValue(out.fit.ngridF,  prefix + ".gridF");
	settings.getValue(out.fit.ngridOm, prefix + ".gridOm");
	settings.getValue(out.fit.nloop,   prefix + ".gridLoops");
	settings.getValue(out.fit.dkappa,  prefix + ".dkappa");
	settings.getValue(out.fit.norm,    prefix + ".norm");

	if ( out.geoDepth2 <= out.geoDepth1 ) {
		SEISCOMP_ERROR("%s.geoDepth2 (%f) must be > geoDepth1 (%f)",
		               prefix.c_str(), out.geoDepth2, out.geoDepth1);
		return false;
	}
	if ( out.nfreq < 8 ) {
		SEISCOMP_ERROR("%s.nfreq (%d) too small", prefix.c_str(), out.nfreq);
		return false;
	}

	return true;
}


}
}
}
