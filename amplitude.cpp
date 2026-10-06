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
 * Copyright (C) 2026 Mustafa Comoglu (Geoscience Australia)               *
 * GNU Affero General Public License Usage - see LICENSE.                   *
 ***************************************************************************/


#define SEISCOMP_COMPONENT MwSpec

#include <seiscomp/logging/log.h>
#include <seiscomp/core/datetime.h>
#include <seiscomp/datamodel/amplitude.h>
#include <seiscomp/datamodel/arrival.h>
#include <seiscomp/datamodel/origin.h>
#include <seiscomp/datamodel/pick.h>
#include <seiscomp/datamodel/sensorlocation.h>
#include <seiscomp/math/fft.h>
#include <seiscomp/math/geo.h>
#include <seiscomp/math/mean.h>
#include <seiscomp/math/windows/cosine.h>
#include <seiscomp/seismology/ttt.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>

#include "mwspec.h"
#include "version.h"


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


void writeArray(std::ostream &os, const char *name, const std::vector<double> &v) {
	os << "\"" << name << "\":[";
	for ( size_t i = 0; i < v.size(); ++i ) {
		os << (i ? "," : "") << (std::isfinite(v[i]) ? v[i] : -30.0);
	}
	os << "]";
}

}  // namespace


/**
 * Everything the spectrum viewer needs about one measurement. Written as
 * JSON to $MWSPEC_DUMP_DIR when the measurement ends, also when it is
 * rejected (status says why), and kept as the processor's spectral
 * diagnostics. Log arrays include the Q/kappa correction and `calibration`
 * but not the gain: subtract log10(gain) for nm*s.
 */
struct SpectrumDump {
	std::string path;
	std::string stream, phase, onset, status = "error";
	std::string signalBegin, signalEnd, noiseBegin, noiseEnd;
	double gain = 0, calibration = 0, travelTime = 0;
	std::vector<double> freq, logSig, logNoise, logCorr;
	bool   fitted = false;
	double om0Log10 = 0, fc = 0, fmin = 0, fmax = 0, residual = 0, snrLog10 = 0;
	double deltaKappa = 0;
	bool   manualBand = false;
	OPT(Core::TimeWindow) signalWindow, noiseWindow;
#ifdef MWSPEC_SPECTRAL_DIAGNOSTICS
	Processing::SpectralDiagnostics *diag = nullptr;

	void toDiagnostics(Processing::SpectralDiagnostics &d) const {
		using Processing::SpectralCurve;
		using Processing::SpectralParameter;

		d.clear();
		d.status = status;

		const std::string cha = stream.substr(stream.rfind('.') + 1);
		const char comp = cha.empty() ? 0 : cha.back();

		if ( signalWindow ) {
			d.windows.push_back({ SpectralCurve::Signal, comp, *signalWindow });
		}
		if ( noiseWindow ) {
			d.windows.push_back({ SpectralCurve::Noise, comp, *noiseWindow });
		}
		const double g = gain != 0 ? std::log10(std::fabs(gain)) : 0.0;

		auto curve = [&](SpectralCurve::Role role, const std::string &label,
		                 const std::string &unit) {
			SpectralCurve c;
			c.role = role;
			c.label = cha + " " + label;
			c.component = comp;
			c.unit = unit;
			return c;
		};

		if ( !freq.empty() ) {
			SpectralCurve sig = curve(SpectralCurve::Signal, "corrected", "nm*s");
			SpectralCurve raw = curve(SpectralCurve::Other, "raw", "nm*s");
			SpectralCurve noise = curve(SpectralCurve::Noise, "noise", "nm*s");
			SpectralCurve corr = curve(SpectralCurve::Correction, "Q/kappa correction", "");
			for ( size_t i = 0; i < freq.size(); ++i ) {
				sig.freq.push_back(freq[i]);
				sig.value.push_back(std::pow(10.0, logSig[i] - g));
				raw.freq.push_back(freq[i]);
				raw.value.push_back(std::pow(10.0, logSig[i] - logCorr[i] - calibration - g));
				if ( logNoise[i] > -29 ) {
					noise.freq.push_back(freq[i]);
					noise.value.push_back(std::pow(10.0, logNoise[i] - g));
				}
				corr.freq.push_back(freq[i]);
				corr.value.push_back(std::pow(10.0, logCorr[i]));
			}
			d.curves.push_back(std::move(sig));
			d.curves.push_back(std::move(raw));
			d.curves.push_back(std::move(noise));
			d.curves.push_back(std::move(corr));
		}

		if ( fitted ) {
			SpectralCurve model = curve(SpectralCurve::Model, "Brune model", "nm*s");
			for ( double f : freq ) {
				model.freq.push_back(f);
				model.value.push_back(std::pow(10.0, om0Log10 - g
				                               - std::log10(1.0 + (f / fc) * (f / fc))
				                               - PI * deltaKappa * f));
			}
			d.curves.push_back(std::move(model));

			d.bands.emplace_back(fmin, fmax);

			auto param = [&](const std::string &id, double value,
			                 const std::string &unit, bool marker = false) {
				SpectralParameter p;
				p.id = id;
				p.value = value;
				p.unit = unit;
				p.frequencyMarker = marker;
				d.parameters.push_back(p);
			};

			if ( manualBand ) {
				param("manual band", 1, "");
			}
			param("Omega0", std::pow(10.0, om0Log10 - g), "nm*s");
			param("fc", fc, "Hz", true);
			param("residual", residual, "");
			param("SNR", std::pow(10.0, snrLog10), "");
			if ( deltaKappa != 0 ) {
				param("deltaKappa", deltaKappa, "s");
			}
			param("travelTime", travelTime, "s");
		}
	}
#endif

	~SpectrumDump() {
#ifdef MWSPEC_SPECTRAL_DIAGNOSTICS
		if ( diag ) {
			toDiagnostics(*diag);
		}
#endif
		if ( path.empty() ) {
			return;
		}
		std::ofstream os(path);
		os.precision(7);
		os << "{\"stream\":\"" << stream << "\",\"phase\":\"" << phase
		   << "\",\"onset\":\"" << onset << "\",\"status\":\"" << status
		   << "\",\"signalWindow\":[\"" << signalBegin << "\",\"" << signalEnd
		   << "\"],\"noiseWindow\":[\"" << noiseBegin << "\",\"" << noiseEnd
		   << "\"],\"gain\":" << gain << ",\"calibration\":" << calibration
		   << ",\"travelTime\":" << travelTime << ",";
		writeArray(os, "freq", freq); os << ",";
		writeArray(os, "logSig", logSig); os << ",";
		writeArray(os, "logNoise", logNoise); os << ",";
		writeArray(os, "logCorr", logCorr);
		if ( fitted ) {
			os << ",\"fit\":{\"om0Log10\":" << om0Log10 << ",\"fc\":" << fc
			   << ",\"fmin\":" << fmin << ",\"fmax\":" << fmax
			   << ",\"residual\":" << residual << ",\"snrLog10\":" << snrLog10 << "}";
		}
		os << "}\n";
	}
};


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

	if ( const char *dir = std::getenv("MWSPEC_DUMP_DIR") ) {
		_dumpDir = dir;
	}

	applyConfig();
}


void AmplitudeProcessor_MwSpec::applyConfig() {
	// Component depends on the phase: P on the vertical, S on a horizontal.
	setUsedComponent(_cfg.phase == 'S' ? FirstHorizontal : Vertical);
	applyWindows();
}


void AmplitudeProcessor_MwSpec::applyWindows() {
	// Trigger-relative windows. The noise window has the same length as the
	// signal window and ends noiseGap before the phase onset at the trigger;
	// for S it therefore stays ahead of P while the signal moves to S.
	const double sigStart = -_cfg.signalPreTime;
	const double sigEnd   = std::max(_cfg.signalDuration, _signalEnd - _signalShift);
	const double noiEnd   = sigStart - _cfg.noiseGap;
	const double noiStart = noiEnd - (sigEnd - sigStart);

	setSignalStart(_signalShift + sigStart);
	setSignalEnd(_signalShift + sigEnd);
	setNoiseStart(noiStart);
	setNoiseEnd(noiEnd);
}


double AmplitudeProcessor_MwSpec::sOnsetShift(const DataModel::Origin *hypocenter,
                                              const DataModel::SensorLocation *receiver,
                                              const DataModel::Pick *pick) {
	_onsetSource.clear();
	if ( _cfg.phase != 'S' ) {
		return 0.0;
	}
	if ( _cfg.sOnset == MwSpecConfig::SOnsetTrigger || !_trigger || !hypocenter ) {
		_onsetSource = "trigger";
		return 0.0;
	}

	// 1. The earliest S-type pick of this station associated with the origin.
	if ( _cfg.sOnset == MwSpecConfig::SOnsetAuto && pick ) {
		const auto &wid = pick->waveformID();
		OPT(Core::Time) best;
		for ( size_t i = 0; i < hypocenter->arrivalCount(); ++i ) {
			const DataModel::Arrival *arr = hypocenter->arrival(i);
			std::string phase;
			try { phase = arr->phase().code(); }
			catch ( ... ) { continue; }
			if ( phase.empty() || ::toupper(phase[0]) != 'S' ) {
				continue;
			}
			const DataModel::Pick *sp = DataModel::Pick::Find(arr->pickID());
			if ( !sp || sp->waveformID().networkCode() != wid.networkCode() ||
			     sp->waveformID().stationCode() != wid.stationCode() ) {
				continue;
			}
			try {
				const Core::Time t = sp->time().value();
				if ( !best || t < *best ) {
					best = t;
				}
			}
			catch ( ... ) {}
		}
		if ( best ) {
			_onsetSource = "pick";
			return (*best - *_trigger).length();
		}
	}

	// 2. The first S-type phase of the travel-time table, configured as for
	//    the SeisComP amplitude time-window expressions.
	if ( receiver ) {
		try {
			TravelTimeTableInterfacePtr ttt = TravelTimeTableInterfaceFactory::Create(
				config().ttInterface.empty() ? "LOCSAT" : config().ttInterface.c_str());
			if ( ttt && ttt->setModel(config().ttModel.empty() ? "iasp91" : config().ttModel) ) {
				double elev = 0.0;
				try { elev = receiver->elevation(); }
				catch ( ... ) {}
				TravelTimeList *tts = ttt->compute(
					hypocenter->latitude().value(), hypocenter->longitude().value(),
					_srcDepthKm, receiver->latitude(), receiver->longitude(), elev);
				if ( tts ) {
					tts->sortByTime();
					for ( const auto &tt : *tts ) {
						if ( !tt.phase.empty() && tt.phase[0] == 'S' ) {
							const Core::Time onset = *_originTime + Core::TimeSpan(tt.time);
							delete tts;
							_onsetSource = "ttt";
							return (onset - *_trigger).length();
						}
					}
					delete tts;
				}
			}
		}
		catch ( std::exception &e ) {
			SEISCOMP_WARNING("%s: S travel time failed: %s", type().c_str(), e.what());
		}
	}

	_onsetSource = "trigger";
	return 0.0;
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


#ifdef MWSPEC_SPECTRAL_DIAGNOSTICS
bool AmplitudeProcessor_MwSpec::setSpectralBand(double fmin, double fmax) {
	if ( fmin <= 0.0 ) {
		_bandFmin = _bandFmax = 0.0;
		return true;
	}

	if ( fmax <= fmin ) {
		return false;
	}

	_bandFmin = fmin;
	_bandFmax = fmax;
	return true;
}
#endif


int AmplitudeProcessor_MwSpec::capabilities() const {
#ifdef MWSPEC_SPECTRAL_DIAGNOSTICS
	return Processing::AmplitudeProcessor::capabilities() | Spectrum;
#else
	return Processing::AmplitudeProcessor::capabilities();
#endif
}


Processing::AmplitudeProcessor::IDList
AmplitudeProcessor_MwSpec::capabilityParameters(Capability cap) const {
	return Processing::AmplitudeProcessor::capabilityParameters(cap);
}


bool AmplitudeProcessor_MwSpec::setParameter(Capability cap,
                                             const std::string &value) {
	return Processing::AmplitudeProcessor::setParameter(cap, value);
}


void AmplitudeProcessor_MwSpec::setEnvironment(
		const DataModel::Origin *hypocenter,
		const DataModel::SensorLocation *receiver,
		const DataModel::Pick *pick) {
	Processing::AmplitudeProcessor::setEnvironment(hypocenter, receiver, pick);

	_srcDepthKm = 0.0;
	_originTime = Core::None;
	_rhypKm = 0.0;
	_signalShift = 0.0;
	_signalEnd = 0.0;
	_epiKm = 0.0;
	_onsetSource.clear();

	if ( !hypocenter ) {
		applyWindows();
		return;
	}

	try { _srcDepthKm = hypocenter->depth().value(); }
	catch ( ... ) {}
	try { _originTime = hypocenter->time().value(); }
	catch ( ... ) {}

	if ( receiver ) {
		try {
			double dDeg, az, baz;
			Math::Geo::delazi(hypocenter->latitude().value(),
			                  hypocenter->longitude().value(),
			                  receiver->latitude(), receiver->longitude(),
			                  &dDeg, &az, &baz);
			_epiKm = Math::Geo::deg2km(dDeg);
			_rhypKm = std::sqrt(_epiKm * _epiKm + _srcDepthKm * _srcDepthKm);
		}
		catch ( ... ) {}
	}

	// Move the S signal window to the S onset (needs the trigger, which
	// scamp and scolv set before the environment).
	_signalShift = _originTime ? sOnsetShift(hypocenter, receiver, pick) : 0.0;

	// S at regional distances: make the window reach past the Lg arrival.
	if ( _cfg.phase == 'S' && _cfg.lgVelocity > 0.0 && _epiKm > 0.0 &&
	     _originTime && _trigger &&
	     Math::Geo::km2deg(_epiKm) <= _cfg.lgMaxDistanceDeg ) {
		const Core::Time lgEnd = *_originTime +
			Core::TimeSpan(_epiKm / _cfg.lgVelocity + _cfg.lgMargin);
		_signalEnd = (lgEnd - *_trigger).length();
	}
	applyWindows();
}


void AmplitudeProcessor_MwSpec::writeDiagnostics(DataModel::Amplitude *amplitude,
                                                 const std::string &suffix) const {
	if ( !_fit.valid ) {
		return;
	}

	setComment(amplitude, "Om0" + suffix, _fit.omega0);
	setComment(amplitude, "fc" + suffix, _fit.cornerFreq);
	setComment(amplitude, "fmin" + suffix, _fit.fmin);
	setComment(amplitude, "fmax" + suffix, _fit.fmax);
	setComment(amplitude, "fitResidual" + suffix, _fit.residual, "%.3f");
	if ( _fit.deltaKappa != 0.0 ) {
		setComment(amplitude, "deltaKappa" + suffix, _fit.deltaKappa);
	}
	if ( !_cfg.useAttenTable ) {
		setComment(amplitude, "travelTime" + suffix, _fit.travelTime, "%.2f");
	}
	if ( !_onsetSource.empty() ) {
		setComment(amplitude, "sOnset" + suffix, _onsetSource);
	}
}


void AmplitudeProcessor_MwSpec::finalizeAmplitude(DataModel::Amplitude *amplitude) const {
	if ( !amplitude ) {
		return;
	}

	amplitude->setMethodID(std::string("Brune/") + _cfg.phase +
	                       (_cfg.useAttenTable ? "/table" : "/Q"));

	try {
		amplitude->creationInfo().setVersion(MWSPEC_VERSION);
	}
	catch ( ... ) {
		DataModel::CreationInfo ci;
		ci.setVersion(MWSPEC_VERSION);
		amplitude->setCreationInfo(ci);
	}

	writeDiagnostics(amplitude, "");
}


void AmplitudeProcessor_MwSpec::prepareData(DoubleArray &data) {
	// The base class deconvolves to the configured data unit (displacement).
	// Validate the metadata required for that here so failures are explicit.
	const Processing::Stream &sc = _streamConfig[targetComponent()];

#ifdef MWSPEC_SPECTRAL_DIAGNOSTICS
	_diag.clear();
#endif

	if ( sc.gain == 0.0 ) {
		setStatus(MissingGain, 1);
#ifdef MWSPEC_SPECTRAL_DIAGNOSTICS
		_diag.status = "missing gain";
#endif
		return;
	}

	SignalUnit unit;
	if ( !unit.fromString(sc.gainUnit.c_str()) ) {
		setStatus(IncompatibleUnit, 2);
#ifdef MWSPEC_SPECTRAL_DIAGNOSTICS
		_diag.status = "incompatible gain unit " + sc.gainUnit;
#endif
		return;
	}

	if ( _enableResponses ) {
		Processing::Sensor *sensor = sc.sensor();
		if ( !sensor || !sensor->response() ) {
			setStatus(MissingResponse, 1);
#ifdef MWSPEC_SPECTRAL_DIAGNOSTICS
			_diag.status = "missing response";
#endif
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

	_fit = FitDiagnostics();

	const double fsamp = _stream.fsamp;
	if ( fsamp <= 0.0 ) {
		setStatus(Error, 1);
		return false;
	}

	const int dataSize = static_cast<int>(data.size());
	const int sigStart = static_cast<int>(i1);
	const int sigEnd   = static_cast<int>(i2);
	const int nsig     = sigEnd - sigStart;

	const Core::Time dataStart = dataTimeWindow().startTime();
	auto timeAt = [&](int idx) {
		return (dataStart + Core::TimeSpan(idx / fsamp)).iso();
	};

	SpectrumDump dump;
	dump.stream = _environment.networkCode + "." + _environment.stationCode + "." +
	              _environment.locationCode + "." +
	              _streamConfig[targetComponent()].code();
	if ( !_dumpDir.empty() ) {
		dump.path = _dumpDir + "/" + dump.stream + ".json";
	}
	dump.phase = std::string(1, _cfg.phase);
	dump.onset = _onsetSource;
	dump.gain = _streamConfig[targetComponent()].gain;
	dump.calibration = _cfg.calibration;
#ifdef MWSPEC_SPECTRAL_DIAGNOSTICS
	dump.diag = &_diag;
#endif

	if ( nsig < 16 || sigStart < 0 || sigEnd > dataSize ) {
		setStatus(Error, 2);
		dump.status = "signal window incomplete";
		return false;
	}

	// Noise: same length as the signal, ending noiseGap before the phase
	// onset. Without a signal shift that is just ahead of the signal window;
	// an S window moved to the S onset keeps its noise ahead of P.
	int noiEnd;
	if ( _signalShift == 0.0 ) {
		noiEnd = sigStart - static_cast<int>(_cfg.noiseGap * fsamp);
	}
	else {
		const double dt0 = (*_trigger - dataStart).length();
		noiEnd = static_cast<int>((dt0 + config().noiseEnd) * fsamp + 0.5);
	}
	int noiStart = noiEnd - nsig;
	if ( noiStart < 0 ) {
		noiStart = 0;
	}
	const int nnoise = noiEnd - noiStart;
	const bool haveNoise = (nnoise >= 16 && noiEnd <= dataSize && noiEnd > noiStart);

	dump.signalBegin = timeAt(sigStart);
	dump.signalEnd = timeAt(sigEnd);
	dump.signalWindow = Core::TimeWindow(dataStart + Core::TimeSpan(sigStart / fsamp),
	                                     dataStart + Core::TimeSpan(sigEnd / fsamp));
	if ( haveNoise ) {
		dump.noiseBegin = timeAt(noiStart);
		dump.noiseEnd = timeAt(noiEnd);
		dump.noiseWindow = Core::TimeWindow(dataStart + Core::TimeSpan(noiStart / fsamp),
		                                    dataStart + Core::TimeSpan(noiEnd / fsamp));
	}

	// --- spectra ----------------------------------------------------------
	std::vector<double> sFreq, sAmp, nFreq, nAmp;
	if ( !displacementSpectrum(data.typedData() + sigStart, nsig, fsamp,
	                           _cfg.applyTaper, sFreq, sAmp) ) {
		setStatus(Error, 3);
		dump.status = "spectrum failed";
		return false;
	}
	if ( haveNoise ) {
		displacementSpectrum(data.typedData() + noiStart, nnoise, fsamp,
		                     _cfg.applyTaper, nFreq, nAmp);
	}

	// --- source/attenuation parameters at the hypocentre ------------------
	// Values were copied in setEnvironment(); the Origin itself may be gone.
	// The travel time of the analysed phase (P: the trigger; S: its onset).
	double travelTime = 0.0;
	if ( _originTime && _trigger ) {
		travelTime = (*_trigger - *_originTime).length() + _signalShift;
	}
	if ( travelTime < 0.0 ) {
		travelTime = 0.0;
	}
	dump.travelTime = travelTime;

	const SourceParams sp = _cfg.model.paramsAt(_srcDepthKm, _cfg.phase);

	// Hypocentral distance [km], needed only for the empirical attenuation
	// table (which folds geometric spreading + anelastic into one term).
	const double rhyp = _rhypKm;
	if ( _cfg.useAttenTable && rhyp <= 0.0 ) {
		setStatus(Error, 7);   // need geometry for the attenuation table
		dump.status = "no geometry for the attenuation table";
		return false;
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
		dump.status = "empty frequency range";
		return false;
	}

	const int nf = _cfg.nfreq;
	std::vector<double> farray(nf), logSig(nf), logNoise(nf), logCorr(nf);
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

		logCorr[i] = corrLog10;
		logSig[i] = logLogInterp(sFreq, sAmp, f) + corrLog10 + _cfg.calibration;

		if ( haveNoise && !nFreq.empty() ) {
			logNoise[i] = logLogInterp(nFreq, nAmp, f) + corrLog10 + _cfg.calibration;
		}
		else {
			logNoise[i] = -30.0;  // no noise estimate -> effectively infinite S/N
		}
	}

	dump.freq = farray;
	dump.logSig = logSig;
	dump.logNoise = logNoise;
	dump.logCorr = logCorr;

	// --- frequency band ---------------------------------------------------
	double fmin, fmax;
	double snrLog10 = 3.0;
	if ( _bandFmin > 0.0 ) {
		// Band set by the analyst in the review window
		fmin = std::max(_bandFmin, farray.front());
		fmax = std::min(_bandFmax, farray.back());
		if ( fmax <= fmin ) {
			setStatus(Error, 8);
			dump.status = "manual band outside the spectrum";
			return false;
		}
		dump.manualBand = true;
		if ( haveNoise ) {
			snrLog10 = -30.0;
			for ( int i = 0; i < nf; ++i ) {
				if ( farray[i] >= fmin && farray[i] <= fmax ) {
					snrLog10 = std::max(snrLog10, logSig[i] - logNoise[i]);
				}
			}
		}
	}
	else if ( _cfg.fixedFmin > 0.0 && _cfg.fixedFmax > 0.0 ) {
		fmin = _cfg.fixedFmin;
		fmax = _cfg.fixedFmax;
	}
	else {
		SNBand band = selectSNBand(farray, logSig, logNoise);
		if ( !band.ok ) {
			setStatus(LowSNR, 0);
			dump.status = "no usable S/N band";
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
		dump.status = "Brune fit failed";
		return false;
	}

	dump.fitted = true;
	dump.om0Log10 = fit.omega0Log10;
	dump.fc = fit.cornerFreq;
	dump.fmin = fmin;
	dump.fmax = fmax;
	dump.residual = fit.residual;
	dump.snrLog10 = snrLog10;
	dump.deltaKappa = fit.deltaKappa;

	if ( fit.cornerFreq <= 0.0 || fit.residual > _cfg.maxResidual ) {
		setStatus(Error, 6);
		dump.status = "fit residual above maxResidual";
		return false;
	}

	const double snrLinear = std::pow(10.0, snrLog10);
	if ( _cfg.minSNR > 0.0 && snrLinear < _cfg.minSNR ) {
		setStatus(LowSNR, snrLinear);
		dump.status = "SNR below minSNR";
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
		dump.status = "missing gain";
		return false;
	}
	dump.status = "ok";

	// Omega0 (linear flat level) in nm*s, in true ground-displacement units.
	amplitude->value = std::pow(10.0, fit.omega0Log10) / gain;

	// Carry the corner frequency as a "period": the base divides the value by
	// fsamp to obtain seconds, so fsamp/fc yields the corner period 1/fc.
	*period = fsamp / fit.cornerFreq;

	*snr = snrLinear;

	dt->index = 0.5 * (sigStart + sigEnd);
	dt->begin = sigStart - dt->index;
	dt->end   = sigEnd - dt->index;

	_fit.valid      = true;
	_fit.omega0     = amplitude->value;
	_fit.cornerFreq = fit.cornerFreq;
	_fit.fmin       = fmin;
	_fit.fmax       = fmax;
	_fit.residual   = fit.residual;
	_fit.deltaKappa = fit.deltaKappa;
	_fit.travelTime = travelTime;

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
