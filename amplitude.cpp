/***************************************************************************
 * SeisComP spectral moment-magnitude plugin (mwspec)                      *
 *                                                                         *
 * Amplitude processor: deconvolves the phase window to ground displacement *
 * (handled by the base via setDataUnit(Meter)), builds the displacement    *
 * amplitude spectrum, corrects it for Q and kappa attenuation, selects the *
 * usable S/N band and fits the Brune omega-square model. The result is the *
 * spectral flat level Omega0 (carried in nm*s) and the corner frequency    *
 * (carried as the amplitude "period"). Ported from Seisan SPEC/AUTOMAG.    *
 *                                                                         *
 * GNU Affero General Public License Usage - see LICENSE.                   *
 ***************************************************************************/


#define SEISCOMP_COMPONENT MwSpec

#include <seiscomp/logging/log.h>
#include <seiscomp/core/datetime.h>
#include <seiscomp/datamodel/origin.h>
#include <seiscomp/datamodel/sensorlocation.h>
#include <seiscomp/math/fft.h>
#include <seiscomp/math/geo.h>
#include <seiscomp/math/mean.h>
#include <seiscomp/math/windows/cosine.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "mwspec.h"


namespace Seiscomp {
namespace Magnitudes {
namespace MwSpec {


// The registered processor is the component combiner (see combiner.cpp); the
// worker below is instantiated as its per-component child.


namespace {

const double PI = 3.141592654;

int nextPow2(int n) {
	int p = 1;
	while ( p < n ) {
		p <<= 1;
	}
	return p;
}

/**
 * Displacement amplitude spectrum of a real segment (already in metres).
 * Fills freq[k] [Hz] and ampNmS[k] [nm·s] for the usable FFT bins
 * (k = 1 .. fftn/2-1), using U(f) = |FFT|/fsamp (the same normalisation as
 * Seisan spec_value), converted to nanometre·second.
 */
bool displacementSpectrum(const double *seg, int npts, double fsamp,
                          bool taper, std::vector<double> &freq,
                          std::vector<double> &ampNmS) {
	if ( npts < 8 || fsamp <= 0.0 ) {
		return false;
	}

	std::vector<double> work(seg, seg + npts);

	// Demean.
	double mean = Math::Statistics::mean(npts, work.data());
	for ( int i = 0; i < npts; ++i ) {
		work[i] -= mean;
	}

	// Quarter-sine taper, 5% of the samples at each end (Seisan SPEC
	// `applytaper`). Using Seisan's exact taper shape makes the spectral level
	// reproduce Seisan bit-for-bit (Math::fft already zero-pads to the next
	// power of two, as Seisan's prepare() does); see the validation in
	// docs/SEISAN_SPECTRUM_COMPARISON.md. A Hann/cosine taper instead reads
	// ~x0.82 (~0.05 Mw) low in the Omega0 band.
	if ( taper ) {
		const int ntap = std::max(2, static_cast<int>(0.05 * npts));
		const double arg = (PI / 2.0) / (ntap - 1);
		for ( int i = 0; i < ntap; ++i ) {
			const double w = std::sin(i * arg);
			work[i] *= w;
			work[npts - 1 - i] *= w;
		}
	}

	Math::ComplexArray spec;
	Math::fft(spec, npts, work.data());
	if ( spec.empty() ) {
		return false;
	}

	const int fftn = nextPow2(npts);
	const double df = fsamp / fftn;
	const int nbins = fftn / 2;
	const double dt = 1.0 / fsamp;

	freq.clear();
	ampNmS.clear();
	freq.reserve(nbins);
	ampNmS.reserve(nbins);

	// Skip bin 0 (DC / packed Nyquist). Go up to fftn/2-1 (valid in both the
	// FFTW and the in-tree real-FFT packings).
	for ( int k = 1; k < nbins; ++k ) {
		if ( k >= static_cast<int>(spec.size()) ) {
			break;
		}
		double amp = std::abs(spec[k]) * dt * 1.0e9;  // metres·s -> nm·s
		freq.push_back(k * df);
		ampNmS.push_back(amp);
	}

	return !freq.empty();
}

/**
 * Linear interpolation of log10(amplitude) versus log10(frequency), clamped
 * to the endpoints. @p freq is ascending. Returns the interpolated log10
 * amplitude at @p f, or a very small level if inputs are degenerate.
 */
double logLogInterp(const std::vector<double> &freq,
                    const std::vector<double> &amp, double f) {
	const size_t n = std::min(freq.size(), amp.size());
	if ( n == 0 || f <= 0.0 ) {
		return -30.0;
	}
	if ( f <= freq[0] ) {
		return std::log10(std::max(amp[0], 1.0e-30));
	}
	if ( f >= freq[n - 1] ) {
		return std::log10(std::max(amp[n - 1], 1.0e-30));
	}

	// Binary search for the bracketing bins.
	size_t lo = 0, hi = n - 1;
	while ( hi - lo > 1 ) {
		size_t mid = (lo + hi) / 2;
		if ( freq[mid] <= f ) {
			lo = mid;
		}
		else {
			hi = mid;
		}
	}

	const double lf0 = std::log10(freq[lo]);
	const double lf1 = std::log10(freq[hi]);
	const double la0 = std::log10(std::max(amp[lo], 1.0e-30));
	const double la1 = std::log10(std::max(amp[hi], 1.0e-30));
	const double t = (lf1 == lf0) ? 0.0 : (std::log10(f) - lf0) / (lf1 - lf0);
	return la0 + t * (la1 - la0);
}

}  // namespace


// ---------------------------------------------------------------------------
AmplitudeProcessor_MwSpec::AmplitudeProcessor_MwSpec()
: Processing::AmplitudeProcessor(MWSPEC_TYPE) {
	setUsedComponent(Vertical);
	setDataUnit(Meter);          // deconvolve to ground displacement
	setUnit(MWSPEC_AMP_UNIT);    // Omega0 carried in nm*s
	_enableResponses = true;

	// We perform our own spectral S/N test; disable the scalar SNR gate.
	setMinSNR(0);

	// Full depth/distance range (local to teleseismic). The geometric
	// spreading handles the distance dependence in the magnitude processor.
	setMinDepth(-5);
	setMaxDepth(800);
	setMinDist(0);
	setMaxDist(180);

	applyConfig();
}


void AmplitudeProcessor_MwSpec::applyConfig() {
	// Component depends on the phase: P on the vertical, S on a horizontal.
	setUsedComponent(_cfg.phase == 'S' ? FirstHorizontal : Vertical);

	// Trigger-relative windows. The noise window has the same length as the
	// signal window and sits ahead of it, separated by noiseGap.
	const double sigStart = -_cfg.signalPreTime;
	const double sigEnd   = _cfg.signalDuration;
	const double noiEnd   = sigStart - _cfg.noiseGap;
	const double noiStart = noiEnd - (_cfg.signalDuration + _cfg.signalPreTime);

	setSignalStart(sigStart);
	setSignalEnd(sigEnd);
	setNoiseStart(noiStart);
	setNoiseEnd(noiEnd);
}


bool AmplitudeProcessor_MwSpec::setup(const Processing::Settings &settings) {
	if ( !Processing::AmplitudeProcessor::setup(settings) ) {
		return false;
	}

	_cfg = MwSpecConfig();
	if ( !readMwSpecConfig(settings, std::string("amplitudes.") + type(), _cfg) ) {
		return false;
	}

	// Spectral Mw REQUIRES full instrument-response removal to ground
	// displacement. Re-assert these after the base setup so a config key or a
	// reset cannot silently leave us measuring counts/velocity.
	_enableResponses = true;
	setDataUnit(Meter);

	applyConfig();

	SEISCOMP_DEBUG("%s: phase=%c preTime=%.1f duration=%.1f nfreq=%d minSNR=%.1f",
	               type().c_str(), _cfg.phase, _cfg.signalPreTime,
	               _cfg.signalDuration, _cfg.nfreq, _cfg.minSNR);
	return true;
}


int AmplitudeProcessor_MwSpec::capabilities() const {
	return Processing::AmplitudeProcessor::capabilities();
}


Processing::AmplitudeProcessor::IDList
AmplitudeProcessor_MwSpec::capabilityParameters(Capability cap) const {
	return Processing::AmplitudeProcessor::capabilityParameters(cap);
}


bool AmplitudeProcessor_MwSpec::setParameter(Capability cap,
                                             const std::string &value) {
	return Processing::AmplitudeProcessor::setParameter(cap, value);
}


void AmplitudeProcessor_MwSpec::prepareData(DoubleArray &data) {
	// The base class deconvolves to the configured data unit (displacement).
	// Validate the metadata required for that here so failures are explicit.
	const Processing::Stream &sc = _streamConfig[targetComponent()];

	if ( sc.gain == 0.0 ) {
		setStatus(MissingGain, 1);
		return;
	}

	SignalUnit unit;
	if ( !unit.fromString(sc.gainUnit.c_str()) ) {
		setStatus(IncompatibleUnit, 2);
		return;
	}

	if ( _enableResponses ) {
		Processing::Sensor *sensor = sc.sensor();
		if ( !sensor || !sensor->response() ) {
			setStatus(MissingResponse, 1);
			return;
		}
	}

	Processing::AmplitudeProcessor::prepareData(data);
}


bool AmplitudeProcessor_MwSpec::computeAmplitude(
	const DoubleArray &data,
	size_t i1, size_t i2,
	size_t si1, size_t si2,
	double /*offset*/,
	AmplitudeIndex *dt,
	AmplitudeValue *amplitude,
	double *period, double *snr) {

	const double fsamp = _stream.fsamp;
	if ( fsamp <= 0.0 ) {
		setStatus(Error, 1);
		return false;
	}

	const int dataSize = static_cast<int>(data.size());
	const int sigStart = static_cast<int>(i1);
	const int sigEnd   = static_cast<int>(i2);
	const int nsig     = sigEnd - sigStart;
	if ( nsig < 16 || sigStart < 0 || sigEnd > dataSize ) {
		setStatus(Error, 2);
		return false;
	}

	// Noise: same length, ending noiseGap before the signal window.
	const int gap = static_cast<int>(_cfg.noiseGap * fsamp);
	int noiEnd   = sigStart - gap;
	int noiStart = noiEnd - nsig;
	if ( noiStart < 0 ) {
		noiStart = 0;
	}
	const int nnoise = noiEnd - noiStart;
	const bool haveNoise = (nnoise >= 16 && noiEnd <= dataSize && noiEnd > noiStart);

	// --- spectra ----------------------------------------------------------
	std::vector<double> sFreq, sAmp, nFreq, nAmp;
	if ( !displacementSpectrum(data.typedData() + sigStart, nsig, fsamp,
	                           _cfg.applyTaper, sFreq, sAmp) ) {
		setStatus(Error, 3);
		return false;
	}
	if ( haveNoise ) {
		displacementSpectrum(data.typedData() + noiStart, nnoise, fsamp,
		                     _cfg.applyTaper, nFreq, nAmp);
	}

	// --- source/attenuation parameters at the hypocentre ------------------
	double depth = 0.0;
	double travelTime = 0.0;
	if ( _environment.hypocenter ) {
		try { depth = _environment.hypocenter->depth().value(); }
		catch ( ... ) {}
		if ( _trigger ) {
			try {
				travelTime = (*_trigger - _environment.hypocenter->time().value()).length();
			}
			catch ( ... ) {}
		}
	}
	if ( travelTime < 0.0 ) {
		travelTime = 0.0;
	}

	const SourceParams sp = _cfg.model.paramsAt(depth, _cfg.phase);

	// Hypocentral distance [km], needed only for the empirical attenuation
	// table (which folds geometric spreading + anelastic into one term).
	double rhyp = 0.0;
	if ( _cfg.useAttenTable ) {
		if ( _environment.hypocenter && _environment.receiver ) {
			try {
				const double elat = _environment.hypocenter->latitude().value();
				const double elon = _environment.hypocenter->longitude().value();
				const double slat = _environment.receiver->latitude();
				const double slon = _environment.receiver->longitude();
				double dDeg, az, baz;
				Math::Geo::delazi(elat, elon, slat, slon, &dDeg, &az, &baz);
				const double epiKm = Math::Geo::deg2km(dDeg);
				rhyp = std::sqrt(epiKm * epiKm + depth * depth);
			}
			catch ( ... ) {}
		}
		if ( rhyp <= 0.0 ) {
			setStatus(Error, 7);   // need geometry for the attenuation table
			return false;
		}
	}

	// --- log-spaced evaluation frequencies --------------------------------
	const double windowLen = (nsig / fsamp);
	double flow = (_cfg.fixedFmin > 0.0) ? _cfg.fixedFmin
	                                      : std::max(0.05, 2.0 / windowLen);
	double fhigh = (_cfg.fixedFmax > 0.0) ? _cfg.fixedFmax : fsamp / 2.5;
	if ( fhigh > fsamp / 2.5 ) {
		fhigh = fsamp / 2.5;
	}
	if ( fhigh <= flow ) {
		setStatus(Error, 4);
		return false;
	}

	const int nf = _cfg.nfreq;
	std::vector<double> farray(nf), logSig(nf), logNoise(nf);
	const double llow = std::log10(flow);
	const double lhigh = std::log10(fhigh);

	for ( int i = 0; i < nf; ++i ) {
		const double f = std::pow(10.0, llow + (lhigh - llow) * i / (nf - 1));
		farray[i] = f;

		// Attenuation correction (Seisan get_om_f0): divide the spectrum by
		// exp(-pi f tt / Q(f)) * exp(-pi kappa f), i.e. boost it back up.
		double corrLog10 = 0.0;  // log10 of the multiplicative correction
		if ( _cfg.useAttenTable ) {
			// Empirical table: one additive log10 term (geometric spreading +
			// anelastic). The magnitude processor then skips geometric
			// spreading (R=1). Near-site kappa stays separate below.
			corrLog10 = _cfg.attenTable.correction(f, rhyp);
			if ( sp.kappa != 0.0 ) {
				corrLog10 += (PI * sp.kappa * f) / std::log(10.0);
			}
		}
		else {
			if ( sp.q0 > 0.0 ) {
				const double q = evalQ(sp.q0, f, sp.qalpha, sp.qcorner);
				if ( q > 0.0 ) {
					corrLog10 += (PI * f * travelTime / q) / std::log(10.0);
				}
			}
			if ( sp.kappa != 0.0 ) {
				corrLog10 += (PI * sp.kappa * f) / std::log(10.0);
			}
		}

		logSig[i] = logLogInterp(sFreq, sAmp, f) + corrLog10 + _cfg.calibration;

		if ( haveNoise && !nFreq.empty() ) {
			logNoise[i] = logLogInterp(nFreq, nAmp, f) + corrLog10 + _cfg.calibration;
		}
		else {
			logNoise[i] = -30.0;  // no noise estimate -> effectively infinite S/N
		}
	}

	// --- frequency band ---------------------------------------------------
	double fmin, fmax;
	double snrLog10 = 3.0;
	if ( _cfg.fixedFmin > 0.0 && _cfg.fixedFmax > 0.0 ) {
		fmin = _cfg.fixedFmin;
		fmax = _cfg.fixedFmax;
	}
	else {
		SNBand band = selectSNBand(farray, logSig, logNoise);
		if ( !band.ok ) {
			setStatus(LowSNR, 0);
			return false;
		}
		fmin = band.fmin;
		fmax = band.fmax;
		snrLog10 = band.maxSNRLog10;
	}

	// --- Brune fit --------------------------------------------------------
	BruneFit fit = bruneGridSearch(farray, logSig, fmin, fmax, _cfg.fit);
	if ( !fit.ok ) {
		setStatus(Error, 5);
		return false;
	}
	if ( fit.cornerFreq <= 0.0 || fit.residual > _cfg.maxResidual ) {
		setStatus(Error, 6);
		return false;
	}

	const double snrLinear = std::pow(10.0, snrLog10);
	if ( _cfg.minSNR > 0.0 && snrLinear < _cfg.minSNR ) {
		setStatus(LowSNR, snrLinear);
		return false;
	}

	// --- outputs ----------------------------------------------------------
	// SeisComP's deconvolveFFT removes only the normalised response shape; the
	// overall sensitivity (gain) is NOT removed and must be divided out here,
	// exactly as the ML/MN/A5_2 processors do. Without this the spectral level
	// (and hence M0) is too large by the gain (~1e8..1e10).
	const double gain = std::fabs(_streamConfig[targetComponent()].gain);
	if ( gain == 0.0 ) {
		setStatus(MissingGain, 0);
		return false;
	}

	// Omega0 (linear flat level) in nm*s, in true ground-displacement units.
	amplitude->value = std::pow(10.0, fit.omega0Log10) / gain;

	// Carry the corner frequency as a "period": the base divides the value by
	// fsamp to obtain seconds, so fsamp/fc yields the corner period 1/fc.
	*period = fsamp / fit.cornerFreq;

	*snr = snrLinear;

	dt->index = 0.5 * (sigStart + sigEnd);
	dt->begin = sigStart - dt->index;
	dt->end   = sigEnd - dt->index;

	SEISCOMP_DEBUG("%s.%s.%s %c: Om0=%g nm*s fc=%.3f Hz band=%.2f-%.2f Hz "
	               "res=%.3f tt=%.1fs gain=%g",
	               _environment.networkCode.c_str(),
	               _environment.stationCode.c_str(),
	               _environment.locationCode.c_str(),
	               _cfg.phase, amplitude->value, fit.cornerFreq,
	               fmin, fmax, fit.residual, travelTime, gain);

	(void)si1; (void)si2;
	return true;
}


}
}
}
