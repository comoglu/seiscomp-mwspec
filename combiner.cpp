/***************************************************************************
 * SeisComP spectral moment-magnitude plugin (mwspec)                      *
 *                                                                         *
 * Component combiner: the registered "Mw(spec)" amplitude processor.      *
 * For a P phase it runs one worker on the vertical. For an S phase it runs *
 * two workers (N and E horizontals) and combines their Omega0 — by default *
 * the vector sum sqrt(N^2 + E^2), i.e. the total horizontal S-wave motion, *
 * the standard horizontal combination for spectral moment magnitudes.      *
 * Modelled on the gempa MLc/MLh two-horizontal proxy.                      *
 *                                                                         *
 * Copyright (C) 2026 Mustafa Comoglu (Geoscience Australia)               *
 * GNU Affero General Public License Usage - see LICENSE.                   *
 ***************************************************************************/


#define SEISCOMP_COMPONENT MwSpec

#include <seiscomp/logging/log.h>
#include <seiscomp/client/application.h>
#include <seiscomp/core/optional.h>
#include <seiscomp/datamodel/amplitude.h>

#include <cctype>
#include <cmath>
#include <exception>
#include <functional>

#include "mwspec.h"


namespace Seiscomp {
namespace Magnitudes {
namespace MwSpec {


REGISTER_AMPLITUDEPROCESSOR(AmplitudeProcessor_MwSpecCombiner, MWSPEC_TYPE);


namespace {

// Average of two AmplitudeTime windows (envelope of both).
Processing::AmplitudeProcessor::AmplitudeTime
avgTime(const Processing::AmplitudeProcessor::AmplitudeTime &a,
        const Processing::AmplitudeProcessor::AmplitudeTime &b) {
	Processing::AmplitudeProcessor::AmplitudeTime t;
	t.reference = a.reference;
	t.begin = std::min(a.begin, b.begin);
	t.end   = std::max(a.end, b.end);
	return t;
}

/**
 * Phase from the module configuration, for use before setup(). scamp and
 * scolv choose which streams to subscribe from usedComponent() *before*
 * calling setup(), so an S configuration must already report Horizontal
 * here. A phase given only in a station binding cannot be seen this early.
 */
char configuredPhase() {
	if ( !SCCoreApp ) {
		return 'P';
	}
	for ( const char *prefix : {"magnitudes.", "amplitudes."} ) {
		try {
			const std::string p = SCCoreApp->configGetString(std::string(prefix) + MWSPEC_TYPE + ".phase");
			if ( !p.empty() ) {
				return ::toupper(p[0]) == 'S' ? 'S' : 'P';
			}
		}
		catch ( ... ) {}
	}
	return 'P';
}

}  // namespace


// ---------------------------------------------------------------------------
AmplitudeProcessor_MwSpecCombiner::AmplitudeProcessor_MwSpecCombiner()
: Processing::AmplitudeProcessor(MWSPEC_TYPE) {
	setUnit(MWSPEC_AMP_UNIT);
	setUsedComponent(configuredPhase() == 'S' ? Horizontal : Vertical);

	_c0.setPublishFunction(std::bind(&AmplitudeProcessor_MwSpecCombiner::newAmplitude,
	                                 this, std::placeholders::_1, std::placeholders::_2));
	_c1.setPublishFunction(std::bind(&AmplitudeProcessor_MwSpecCombiner::newAmplitude,
	                                 this, std::placeholders::_1, std::placeholders::_2));

	AmplitudeProcessor_MwSpecCombiner::reset();
}


// ---------------------------------------------------------------------------
int AmplitudeProcessor_MwSpecCombiner::capabilities() const {
	return _c0.capabilities() | Combiner;
}


Processing::AmplitudeProcessor::IDList
AmplitudeProcessor_MwSpecCombiner::capabilityParameters(Capability cap) const {
	if ( cap == Combiner ) {
		IDList l;
		l.push_back("Vector sum");
		l.push_back("Average");
		l.push_back("Geometric mean");
		l.push_back("Max");
		l.push_back("Min");
		return l;
	}
	return _c0.capabilityParameters(cap);
}


bool AmplitudeProcessor_MwSpecCombiner::setParameter(Capability cap,
                                                     const std::string &value) {
	if ( cap == Combiner ) {
		if      ( value == "Vector sum" )     _combiner = CombineVectorSum;
		else if ( value == "Average" )        _combiner = CombineAverage;
		else if ( value == "Geometric mean" ) _combiner = CombineGeometricMean;
		else if ( value == "Max" )            _combiner = CombineMax;
		else if ( value == "Min" )            _combiner = CombineMin;
		else return false;
		return true;
	}
	_c1.setParameter(cap, value);
	return _c0.setParameter(cap, value);
}


std::string AmplitudeProcessor_MwSpecCombiner::parameter(Capability cap) const {
	if ( cap == Combiner ) {
		switch ( _combiner ) {
			case CombineVectorSum:     return "Vector sum";
			case CombineAverage:       return "Average";
			case CombineGeometricMean: return "Geometric mean";
			case CombineMax:           return "Max";
			case CombineMin:           return "Min";
		}
	}
	return _c0.parameter(cap);
}


// ---------------------------------------------------------------------------
void AmplitudeProcessor_MwSpecCombiner::reset() {
	Processing::AmplitudeProcessor::reset();
	_results[0] = _results[1] = Core::None;
	_c0.reset();
	_c1.reset();
}


// ---------------------------------------------------------------------------
bool AmplitudeProcessor_MwSpecCombiner::setup(const Processing::Settings &settings) {
	// Propagate any aliased type to the children.
	_c0._type = _type;
	_c1._type = _type;

	// Read the phase so we know whether to run one (P) or two (S) workers.
	MwSpecConfig cfg;
	if ( !readMwSpecConfig(settings, std::string("amplitudes.") + type(), cfg) ) {
		return false;
	}
	_phase = cfg.phase;

	// The streams were chosen from usedComponent() before setup(); a phase
	// that differs (e.g. set per binding) would never receive its data.
	if ( (_phase == 'S') != (usedComponent() == Horizontal) ) {
		SEISCOMP_ERROR("%s: phase %c conflicts with the module-level phase; set "
		               "magnitudes.%s.phase in the module configuration, not "
		               "per binding", type().c_str(), _phase, type().c_str());
		return false;
	}

	// Optional combiner override (S only).
	std::string c;
	if ( settings.getValue(c, std::string("amplitudes.") + type() + ".combiner") && !c.empty() ) {
		if      ( c == "vector_sum" )     _combiner = CombineVectorSum;
		else if ( c == "average" )        _combiner = CombineAverage;
		else if ( c == "geometric_mean" ) _combiner = CombineGeometricMean;
		else if ( c == "max" )            _combiner = CombineMax;
		else if ( c == "min" )            _combiner = CombineMin;
		else {
			SEISCOMP_ERROR("%s.combiner: unknown '%s' (use vector_sum|average|"
			               "geometric_mean|max|min)", type().c_str(), c.c_str());
			return false;
		}
	}

	if ( _phase == 'S' ) {
		_nActive = 2;
		setUsedComponent(Horizontal);
		_c0.setUsedComponent(FirstHorizontal);
		_c1.setUsedComponent(SecondHorizontal);
		_c0.streamConfig(FirstHorizontalComponent)  = streamConfig(FirstHorizontalComponent);
		_c1.streamConfig(SecondHorizontalComponent) = streamConfig(SecondHorizontalComponent);
	}
	else {
		_nActive = 1;
		setUsedComponent(Vertical);
		_c0.setUsedComponent(Vertical);
		_c0.streamConfig(VerticalComponent) = streamConfig(VerticalComponent);
	}

	if ( !Processing::AmplitudeProcessor::setup(settings) ) {
		return false;
	}

	// The worker's own setup() reads the shared config and (via applyConfig)
	// resets the used component to phase default; re-assert the S split after.
	if ( !_c0.setup(settings) ) {
		return false;
	}
	if ( _nActive == 2 ) {
		if ( !_c1.setup(settings) ) {
			return false;
		}
		_c0.setUsedComponent(FirstHorizontal);
		_c1.setUsedComponent(SecondHorizontal);
	}
	else {
		_c0.setUsedComponent(Vertical);
	}

	return true;
}


// ---------------------------------------------------------------------------
void AmplitudeProcessor_MwSpecCombiner::setTrigger(const Core::Time &trigger) {
	Processing::AmplitudeProcessor::setTrigger(trigger);
	// scamp and scolv set the trigger before setup(), i.e. before _nActive
	// is known: always forward it to both workers.
	_c0.setTrigger(trigger);
	_c1.setTrigger(trigger);
}


void AmplitudeProcessor_MwSpecCombiner::setEnvironment(
		const DataModel::Origin *hypocenter,
		const DataModel::SensorLocation *receiver,
		const DataModel::Pick *pick) {
	_c0.setEnvironment(hypocenter, receiver, pick);
	if ( _nActive == 2 ) {
		_c1.setEnvironment(hypocenter, receiver, pick);
	}
}


void AmplitudeProcessor_MwSpecCombiner::computeTimeWindow() {
	_c0.computeTimeWindow();

	if ( _nActive == 1 ) {
		if ( _c0.isFinished() ) {
			setStatus(_c0.status(), _c0.statusValue());
			setTimeWindow(Core::TimeWindow());
			return;
		}
		setTimeWindow(_c0.timeWindow());
		return;
	}

	_c1.computeTimeWindow();
	if ( _c0.isFinished() ) {
		setStatus(_c0.status(), _c0.statusValue());
		setTimeWindow(Core::TimeWindow());
		return;
	}
	if ( _c1.isFinished() ) {
		setStatus(_c1.status(), _c1.statusValue());
		setTimeWindow(Core::TimeWindow());
		return;
	}
	setTimeWindow(_c0.timeWindow() | _c1.timeWindow());
}


void AmplitudeProcessor_MwSpecCombiner::close() const {}


// ---------------------------------------------------------------------------
bool AmplitudeProcessor_MwSpecCombiner::feed(const Record *record) {
	// No exception may leave feed(): StreamApplication::readRecords() answers
	// an exception from storeRecord() with `delete rec`, although the workers
	// already hold that record as their last record. The next record then
	// releases a freed Record (SIGSEGV in WaveformProcessor::store, see
	// CRASH_DIAGNOSIS_2026-07-22.md). Fail this processor instead.
	try {
		return feedWorkers(record);
	}
	catch ( std::exception &e ) {
		SEISCOMP_ERROR("%s %s: exception while processing data: %s",
		               type().c_str(), record->streamID().c_str(), e.what());
	}
	catch ( ... ) {
		SEISCOMP_ERROR("%s %s: unknown exception while processing data",
		               type().c_str(), record->streamID().c_str());
	}

	setStatus(Error, 10);
	return false;
}


bool AmplitudeProcessor_MwSpecCombiner::feedWorkers(const Record *record) {
	if ( status() > WaveformProcessor::Finished ) {
		return false;
	}

	if ( _nActive == 1 ) {
		if ( _c0.isFinished() ) {
			return false;
		}
		if ( record->channelCode() == _streamConfig[VerticalComponent].code() ) {
			_c0.feed(record);
			if ( _c0.status() == InProgress ) {
				setStatus(WaveformProcessor::InProgress, _c0.statusValue());
			}
			else if ( _c0.isFinished() && !isFinished() ) {
				// A non-Finished terminal status is propagated here; a Finished
				// status is propagated by newAmplitude() after it emits.
				if ( _c0.status() != Finished ) {
					setStatus(_c0.status(), _c0.statusValue());
				}
			}
		}
		return true;
	}

	// S: route each horizontal to its worker.
	if ( _c0.isFinished() && _c1.isFinished() ) {
		return false;
	}
	if ( record->channelCode() == _streamConfig[FirstHorizontalComponent].code() ) {
		if ( !_c0.isFinished() ) {
			_c0.feed(record);
			if ( _c0.status() == InProgress ) {
				setStatus(WaveformProcessor::InProgress, _c0.statusValue());
			}
			else if ( _c0.isFinished() && _c1.isFinished() && !isFinished() ) {
				setStatus(_c0.status() == Finished ? _c1.status() : _c0.status(),
				          _c0.status() == Finished ? _c1.statusValue() : _c0.statusValue());
			}
		}
	}
	else if ( record->channelCode() == _streamConfig[SecondHorizontalComponent].code() ) {
		if ( !_c1.isFinished() ) {
			_c1.feed(record);
			if ( _c1.status() == InProgress ) {
				setStatus(WaveformProcessor::InProgress, _c1.statusValue());
			}
			else if ( _c1.isFinished() && _c0.isFinished() && !isFinished() ) {
				setStatus(_c1.status() == Finished ? _c0.status() : _c1.status(),
				          _c1.status() == Finished ? _c0.statusValue() : _c1.statusValue());
			}
		}
	}

	return true;
}


// ---------------------------------------------------------------------------
void AmplitudeProcessor_MwSpecCombiner::newAmplitude(
		const AmplitudeProcessor *proc,
		const AmplitudeProcessor::Result &res) {

	if ( isFinished() ) {
		return;
	}

	// Single-component (P): forward directly.
	if ( _nActive == 1 ) {
		setStatus(Finished, 100.0);
		emitAmplitude(res);
		return;
	}

	const int idx = (proc == &_c0) ? 0 : 1;
	_results[idx] = ComponentResult();
	_results[idx]->value  = res.amplitude;
	_results[idx]->time   = res.time;
	_results[idx]->snr    = res.snr;
	_results[idx]->period = res.period;

	if ( !_results[0] || !_results[1] ) {
		return;   // wait for the other horizontal
	}

	setStatus(Finished, 100.0);

	const double vN = _results[0]->value.value;
	const double vE = _results[1]->value.value;

	Result out;
	out.record    = res.record;
	out.component = Horizontal;
	out.snr       = 0.5 * (_results[0]->snr + _results[1]->snr);
	out.time      = avgTime(_results[0]->time, _results[1]->time);
	out.period    = (_results[0]->period > 0 && _results[1]->period > 0)
	                ? 0.5 * (_results[0]->period + _results[1]->period)
	                : -1.0;

	switch ( _combiner ) {
		case CombineVectorSum:
			out.amplitude.value = std::sqrt(vN * vN + vE * vE);
			break;
		case CombineGeometricMean:
			out.amplitude.value = std::sqrt(vN * vE);
			break;
		case CombineAverage:
			out.amplitude.value = 0.5 * (vN + vE);
			break;
		case CombineMax:
			out.amplitude = (vN >= vE) ? _results[0]->value : _results[1]->value;
			out.period    = (vN >= vE) ? _results[0]->period : _results[1]->period;
			out.snr       = (vN >= vE) ? _results[0]->snr : _results[1]->snr;
			break;
		case CombineMin:
			out.amplitude = (vN <= vE) ? _results[0]->value : _results[1]->value;
			out.period    = (vN <= vE) ? _results[0]->period : _results[1]->period;
			out.snr       = (vN <= vE) ? _results[0]->snr : _results[1]->snr;
			break;
	}

	emitAmplitude(out);
}


// ---------------------------------------------------------------------------
const Processing::AmplitudeProcessor *
AmplitudeProcessor_MwSpecCombiner::componentProcessor(Component comp) const {
	switch ( comp ) {
		case VerticalComponent:
			return _nActive == 1 ? &_c0 : nullptr;
		case FirstHorizontalComponent:
			return _nActive == 2 ? &_c0 : nullptr;
		case SecondHorizontalComponent:
			return _nActive == 2 ? &_c1 : nullptr;
		default:
			break;
	}
	return nullptr;
}


const DoubleArray *
AmplitudeProcessor_MwSpecCombiner::processedData(Component comp) const {
	switch ( comp ) {
		case VerticalComponent:
			return _nActive == 1 ? _c0.processedData(comp) : nullptr;
		case FirstHorizontalComponent:
			return _nActive == 2 ? _c0.processedData(comp) : nullptr;
		case SecondHorizontalComponent:
			return _nActive == 2 ? _c1.processedData(comp) : nullptr;
		default:
			break;
	}
	return nullptr;
}


// ---------------------------------------------------------------------------
void AmplitudeProcessor_MwSpecCombiner::reprocess(OPT(double) searchBegin,
                                                  OPT(double) searchEnd) {
	setStatus(WaitingForData, 0);
	_results[0] = _results[1] = Core::None;

	_c0.reprocess(searchBegin, searchEnd);
	if ( _nActive == 2 ) {
		_c1.reprocess(searchBegin, searchEnd);
	}

	if ( !isFinished() ) {
		if ( _c0.status() > Finished ) {
			setStatus(_c0.status(), _c0.statusValue());
		}
		else if ( _nActive == 2 && _c1.status() > Finished ) {
			setStatus(_c1.status(), _c1.statusValue());
		}
	}
}


// ---------------------------------------------------------------------------
void AmplitudeProcessor_MwSpecCombiner::finalizeAmplitude(DataModel::Amplitude *amplitude) const {
	if ( !amplitude ) {
		return;
	}

	_c0.finalizeAmplitude(amplitude);
	if ( _nActive == 1 ) {
		return;
	}

	// S: the worker wrote unsuffixed values of N only; replace them by
	// per-component values plus the corner frequency that the magnitude uses
	// (the inverse of the combined period).
	for ( const char *id : {"Om0", "fc", "fmin", "fmax", "fitResidual",
	                        "deltaKappa", "travelTime", "sOnset"} ) {
		amplitude->removeComment(DataModel::CommentIndex(id));
	}
	_c0.writeDiagnostics(amplitude, ".N");
	_c1.writeDiagnostics(amplitude, ".E");

	try {
		const double period = amplitude->period().value();
		if ( period > 0.0 ) {
			setComment(amplitude, "fc", 1.0 / period);
		}
	}
	catch ( ... ) {}
}


// The combiner itself never computes a raw amplitude; the workers do.
bool AmplitudeProcessor_MwSpecCombiner::computeAmplitude(
		const DoubleArray &, size_t, size_t, size_t, size_t,
		double, AmplitudeIndex *, AmplitudeValue *, double *, double *) {
	return false;
}


#ifdef MWSPEC_SPECTRAL_DIAGNOSTICS
bool AmplitudeProcessor_MwSpecCombiner::setSpectralBand(double fmin, double fmax) {
	bool ok = _c0.setSpectralBand(fmin, fmax);
	if ( _nActive == 2 ) {
		ok = _c1.setSpectralBand(fmin, fmax) && ok;
	}
	return ok;
}


const Processing::SpectralDiagnostics *
AmplitudeProcessor_MwSpecCombiner::spectralDiagnostics() const {
	if ( _nActive == 1 ) {
		return _c0.spectralDiagnostics();
	}

	// S: both horizontals in one view, parameters suffixed with the
	// component like the amplitude comments (fc.N, fc.E, ...).
	_diag.clear();
	const AmplitudeProcessor_MwSpec *workers[2] = { &_c0, &_c1 };
	const char *suffix[2] = { ".N", ".E" };
	std::string status;

	for ( int i = 0; i < 2; ++i ) {
		const Processing::SpectralDiagnostics *d = workers[i]->spectralDiagnostics();
		if ( !d || d->empty() ) {
			continue;
		}

		_diag.curves.insert(_diag.curves.end(), d->curves.begin(), d->curves.end());
		_diag.bands.insert(_diag.bands.end(), d->bands.begin(), d->bands.end());
		_diag.windows.insert(_diag.windows.end(), d->windows.begin(), d->windows.end());
		for ( auto p : d->parameters ) {
			p.id += suffix[i];
			_diag.parameters.push_back(p);
		}

		if ( !status.empty() ) {
			status += ", ";
		}
		status += std::string(suffix[i] + 1) + ": " + (d->status.empty() ? "?" : d->status);
	}

	_diag.status = status;
	return &_diag;
}
#endif


}
}
}
