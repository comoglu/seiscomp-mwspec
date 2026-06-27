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
 * GNU Affero General Public License Usage - see LICENSE.                   *
 ***************************************************************************/


#define SEISCOMP_COMPONENT MwSpec

#include <seiscomp/logging/log.h>
#include <seiscomp/core/optional.h>

#include <cmath>
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

}  // namespace


// ---------------------------------------------------------------------------
AmplitudeProcessor_MwSpecCombiner::AmplitudeProcessor_MwSpecCombiner()
: Processing::AmplitudeProcessor(MWSPEC_TYPE) {
	setUnit(MWSPEC_AMP_UNIT);
	setUsedComponent(Vertical);   // overridden in setup() once the phase is known

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
	_c0.setTrigger(trigger);
	if ( _nActive == 2 ) {
		_c1.setTrigger(trigger);
	}
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


// The combiner itself never computes a raw amplitude; the workers do.
bool AmplitudeProcessor_MwSpecCombiner::computeAmplitude(
		const DoubleArray &, size_t, size_t, size_t, size_t,
		double, AmplitudeIndex *, AmplitudeValue *, double *, double *) {
	return false;
}


}
}
}
