/***************************************************************************
 * SeisComP spectral moment-magnitude plugin (mwspec)                      *
 *                                                                         *
 * Magnitude processor: Omega0 -> seismic moment -> Mw, following Seisan    *
 * automag.for. Mw(spec) is itself a moment magnitude, so estimateMw is     *
 * the identity.                                                            *
 *                                                                         *
 * Copyright (C) 2026 Mustafa Comoglu (Geoscience Australia)               *
 * GNU Affero General Public License Usage - see LICENSE.                   *
 ***************************************************************************/


#define SEISCOMP_COMPONENT MwSpec

#include <seiscomp/logging/log.h>
#include <seiscomp/datamodel/stationmagnitude.h>
#include <seiscomp/math/geo.h>

#include "mwspec.h"
#include "version.h"


namespace Seiscomp {
namespace Magnitudes {
namespace MwSpec {


REGISTER_MAGNITUDEPROCESSOR(MagnitudeProcessor_MwSpec, MWSPEC_TYPE);


MagnitudeProcessor_MwSpec::MagnitudeProcessor_MwSpec()
: Processing::MagnitudeProcessor(MWSPEC_TYPE) {}


bool MagnitudeProcessor_MwSpec::setup(const Processing::Settings &settings) {
	if ( !Processing::MagnitudeProcessor::setup(settings) ) {
		return false;
	}

	_cfg = MwSpecConfig();
	if ( !readMwSpecConfig(settings, std::string("magnitudes.") + type(), _cfg) ) {
		return false;
	}

	// Distance gate (enforced by the base before computeMagnitude). Lets users
	// drop teleseismic station magnitudes, e.g. maximumDistance = 12 deg.
	_minimumDistanceDeg = _cfg.minDistanceDeg;
	_maximumDistanceDeg = _cfg.maxDistanceDeg;

	SEISCOMP_DEBUG("%s: phase=%c radiation=%.3f freeSurface=%.2f model layers=%zu "
	               "dist=%.1f-%.1f deg",
	               type().c_str(), _cfg.phase, _cfg.radiation, _cfg.freeSurface,
	               _cfg.model.layerCount(), _cfg.minDistanceDeg, _cfg.maxDistanceDeg);
	return true;
}


MagnitudeProcessor_MwSpec::Status
MagnitudeProcessor_MwSpec::computeMagnitude(
	double amplitude, const std::string &unit,
	double period, double /*snr*/,
	double delta, double depth,
	const DataModel::Origin *,
	const DataModel::SensorLocation *,
	const DataModel::Amplitude *,
	const Locale *,
	double &value) {

	_lastValid = false;

	if ( amplitude <= 0.0 ) {
		return AmplitudeOutOfRange;
	}

	// Omega0 is carried in nm*s (same convention as Mwp).
	double omega0 = amplitude;
	if ( !convertAmplitude(omega0, unit, MWSPEC_AMP_UNIT) ) {
		return InvalidAmplitudeUnit;
	}

	if ( depth < 0.0 ) {
		depth = 0.0;
	}

	const double epicentralKm = Math::Geo::deg2km(delta);
	const SourceParams sp = _cfg.model.paramsAt(depth, _cfg.phase);

	// With the empirical attenuation table the geometric spreading was already
	// removed on the amplitude side (Omega0 is a source level), so use R=1 here
	// and let the one-time `calibration` offset set the absolute level.
	double geoDistKm;
	if ( _cfg.useAttenTable ) {
		geoDistKm = 1.0;
	}
	else {
		geoDistKm = geoSpreadingDistance(_cfg.phase, epicentralKm, depth,
		                                 _cfg.geoDepth1, _cfg.geoDepth2,
		                                 _cfg.herkijDistanceKm);
		if ( geoDistKm <= 0.0 ||
		     geoDistKm >= 0.5 * std::numeric_limits<double>::max() ) {
			return DistanceOutOfRange;
		}
	}

	const double cornerFreq = (period > 0.0) ? (1.0 / period) : 0.0;

	const MomentResult mr = momentFromOmega(omega0, cornerFreq,
	                                        sp.velocity(_cfg.phase), sp.density,
	                                        geoDistKm, _cfg.radiation,
	                                        _cfg.freeSurface);
	if ( mr.m0 <= 0.0 ) {
		return Error;
	}

	value = mr.mw;

	_lastValid = true;
	_last = mr;
	_lastCornerFreq = cornerFreq;
	_lastGeoDistKm = geoDistKm;

	SEISCOMP_DEBUG("%s: Om0=%g nm*s fc=%.3f Hz R=%.1f km v=%.2f rho=%.2f "
	               "-> log M0=%.2f Mw=%.2f",
	               type().c_str(), omega0, cornerFreq, geoDistKm,
	               sp.velocity(_cfg.phase), sp.density, mr.logM0, value);

	return OK;
}


void MagnitudeProcessor_MwSpec::finalizeMagnitude(DataModel::StationMagnitude *magnitude) const {
	if ( !magnitude ) {
		return;
	}

	try {
		magnitude->creationInfo().setVersion(MWSPEC_VERSION);
	}
	catch ( ... ) {
		DataModel::CreationInfo ci;
		ci.setVersion(MWSPEC_VERSION);
		magnitude->setCreationInfo(ci);
	}

	if ( !_lastValid ) {
		return;
	}

	setComment(magnitude, "M0", _last.m0, "%.3e");
	if ( _lastCornerFreq > 0.0 ) {
		setComment(magnitude, "fc", _lastCornerFreq);
		setComment(magnitude, "sourceRadius", _last.sourceRadius, "%.0f");
		setComment(magnitude, "stressDrop", _last.stressDrop * 0.1, "%.3g");  // bar -> MPa
	}
	if ( !_cfg.useAttenTable ) {
		setComment(magnitude, "geoDistance", _lastGeoDistKm, "%.1f");
	}
}


}
}
}
