/***************************************************************************
 * SeisComP spectral moment-magnitude plugin (mwspec)                      *
 *                                                                         *
 * Per-station moment magnitude Mw(spec) from a Brune omega-square fit of   *
 * the displacement spectrum, porting Seisan SPEC/AUTOMAG.                  *
 *                                                                         *
 * GNU Affero General Public License Usage - see LICENSE.                   *
 ***************************************************************************/


#ifndef SEISCOMP_MAGNITUDES_MWSPEC_PLUGIN_H
#define SEISCOMP_MAGNITUDES_MWSPEC_PLUGIN_H


#include <seiscomp/core/plugin.h>
#include <seiscomp/core/version.h>
#include <seiscomp/processing/amplitudeprocessor.h>
#include <seiscomp/processing/magnitudeprocessor.h>

#include <string>
#include <vector>

#include "brune.h"


namespace Seiscomp {
namespace Magnitudes {
namespace MwSpec {


//! Registered amplitude/magnitude type string.
#define MWSPEC_TYPE "Mw(spec)"

//! Unit of the carried spectral flat level (Omega0), as for Mwp.
#define MWSPEC_AMP_UNIT "nm*s"


/**
 * Generic frequency x distance attenuation table — an empirical alternative to
 * the parametric Q model, for agencies with a regionally-calibrated attenuation
 * (e.g. a GIT result). A(f, rhyp) is the ADDITIVE log10 path correction
 * (geometric spreading + anelastic + any regional term) that brings the
 * observed displacement spectrum to the source:  log10(S) = log10(obs) + A.
 * Near-site kappa stays a separate per-station config. Loaded from a CSV named
 * by magnitudes.Mw(spec).attenuationTable; see the plugin's make_atten_table.py
 * for the format and a worked example.
 */
struct AttenTable {
	std::vector<double> logFreqs;          //!< log10(grid frequencies), ascending
	std::vector<double> dists;             //!< grid hypocentral distances [km], ascending
	std::vector<std::vector<double>> A;    //!< A[i][j] = correction at (freq i, dist j)

	bool empty() const {
		return logFreqs.empty() || dists.empty() || A.empty();
	}

	//! Bilinear interpolation (linear in log10 f, linear in distance), clamped
	//! to the grid edges. @p f in Hz, @p rhyp in km. Returns additive log10.
	double correction(double f, double rhyp) const;
};


/**
 * Configuration shared by the amplitude and magnitude processors. Both read
 * the same keys (under their respective "amplitudes."/"magnitudes." prefix)
 * so that a single binding describes one consistent spectral model.
 */
struct MwSpecConfig {
	char     phase            = 'P';     //!< 'P' or 'S'
	SpecModel model;                      //!< layered velocity/Q/density model

	// Moment computation
	double   radiation        = 0.6;     //!< average radiation pattern
	double   freeSurface      = 2.0;     //!< free-surface amplification factor
	double   geoDepth1        = 50.0;    //!< [km] body->surface spreading transition top
	double   geoDepth2        = 100.0;   //!< [km] body->surface spreading transition bottom
	double   herkijDistanceKm = 100.0;   //!< [km] Herrmann-Kijko transition distance
	double   minDistanceDeg   = 0.0;     //!< accept station mags from this distance [deg]
	double   maxDistanceDeg   = 180.0;   //!< ...up to this distance (set ~12 to drop teleseismic)

	// Spectral measurement (amplitude side)
	double   signalPreTime    = 1.0;     //!< [s] window lead before the phase onset
	double   signalDuration   = 25.0;    //!< [s] signal window length after the onset
	double   noiseGap         = 1.0;     //!< [s] gap between noise and signal windows
	double   fixedFmin        = 0.0;     //!< [Hz] fixed band low edge (0 = automatic)
	double   fixedFmax        = 0.0;     //!< [Hz] fixed band high edge (0 = automatic)
	int      nfreq            = 100;     //!< number of log-spaced spectral samples
	double   minSNR           = 2.0;     //!< minimum amplitude SNR to accept
	double   maxResidual      = 1.0;     //!< maximum Brune-fit residual to accept
	double   calibration      = 0.0;     //!< log10 additive level calibration (Seisan match)
	bool     applyTaper       = true;    //!< cosine taper before the FFT
	BruneFitOptions fit;                  //!< grid-search tunables

	// Optional empirical attenuation table (alternative to the parametric Q +
	// geometric spreading). When a path is configured the amplitude processor
	// corrects the spectrum with the table (geometric + anelastic) and the
	// magnitude processor skips its own geometric spreading (R=1); the absolute
	// level is then set once via `calibration`. See AttenTable.
	std::string attenTablePath;           //!< magnitudes.Mw(spec).attenuationTable
	bool        useAttenTable = false;    //!< true when attenTablePath is set+loaded
	AttenTable  attenTable;               //!< parsed grid (amplitude side uses it)
};


/**
 * Reads MwSpecConfig from a station binding. @p prefix is "amplitudes.Mw(spec)"
 * or "magnitudes.Mw(spec)". Missing keys keep their defaults; the function only
 * fails on a malformed spectral model definition.
 */
bool readMwSpecConfig(const Processing::Settings &settings,
                      const std::string &prefix, MwSpecConfig &out);


// ---------------------------------------------------------------------------
//  Amplitude processor: builds the displacement spectrum, fits the Brune
//  model and emits Omega0 (flat level, nm*s) with the corner frequency as the
//  carried "period".
// ---------------------------------------------------------------------------
class SC_SYSTEM_CLIENT_API AmplitudeProcessor_MwSpec : public Processing::AmplitudeProcessor {
	public:
		AmplitudeProcessor_MwSpec();

	public:
		bool setup(const Processing::Settings &settings) override;
		int capabilities() const override;
		Processing::AmplitudeProcessor::IDList
		    capabilityParameters(Capability cap) const override;
		bool setParameter(Capability cap, const std::string &value) override;

	protected:
		void prepareData(DoubleArray &data) override;

		bool computeAmplitude(const DoubleArray &data,
		                      size_t i1, size_t i2,
		                      size_t si1, size_t si2,
		                      double offset,
		                      AmplitudeIndex *dt,
		                      AmplitudeValue *amplitude,
		                      double *period, double *snr) override;

	private:
		void applyConfig();

	private:
		MwSpecConfig _cfg;

	friend class AmplitudeProcessor_MwSpecCombiner;
};


// ---------------------------------------------------------------------------
//  Component combiner (registered as "Mw(spec)"). For P it runs a single
//  worker on the vertical; for S it runs two workers on the horizontals and
//  combines their Omega0 (default: vector sum sqrt(N^2 + E^2) = total S-wave
//  horizontal motion). This is the entry point scamp/scolv instantiates.
// ---------------------------------------------------------------------------
class SC_SYSTEM_CLIENT_API AmplitudeProcessor_MwSpecCombiner : public Processing::AmplitudeProcessor {
	public:
		AmplitudeProcessor_MwSpecCombiner();

	public:
		int capabilities() const override;
		Processing::AmplitudeProcessor::IDList
		    capabilityParameters(Capability cap) const override;
		bool setParameter(Capability cap, const std::string &value) override;
		std::string parameter(Capability cap) const override;

		void reset() override;
		bool setup(const Processing::Settings &settings) override;

		void setTrigger(const Core::Time &trigger) override;
		void setEnvironment(const DataModel::Origin *hypocenter,
		                    const DataModel::SensorLocation *receiver,
		                    const DataModel::Pick *pick) override;
		void computeTimeWindow() override;
		void close() const override;
		bool feed(const Record *record) override;

		const AmplitudeProcessor *componentProcessor(Component comp) const override;
		const DoubleArray *processedData(Component comp) const override;

		void reprocess(OPT(double) searchBegin, OPT(double) searchEnd) override;

	protected:
		bool computeAmplitude(const DoubleArray &data,
		                      size_t i1, size_t i2,
		                      size_t si1, size_t si2,
		                      double offset,
		                      AmplitudeIndex *dt,
		                      AmplitudeValue *amplitude,
		                      double *period, double *snr) override;

	private:
		void newAmplitude(const AmplitudeProcessor *proc,
		                  const AmplitudeProcessor::Result &res);

		struct ComponentResult {
			Processing::AmplitudeProcessor::AmplitudeValue value;
			Processing::AmplitudeProcessor::AmplitudeTime  time;
			double snr;
			double period;
		};

		enum CombineMode {
			CombineVectorSum,      //!< sqrt(N^2 + E^2): total horizontal S motion
			CombineAverage,
			CombineGeometricMean,
			CombineMax,
			CombineMin
		};

		mutable AmplitudeProcessor_MwSpec _c0;   //!< vertical (P) or N (S)
		mutable AmplitudeProcessor_MwSpec _c1;   //!< E (S only)
		char     _phase    = 'P';
		int      _nActive  = 1;                  //!< 1 for P, 2 for S
		CombineMode _combiner = CombineVectorSum;
		OPT(ComponentResult) _results[2];
};


// ---------------------------------------------------------------------------
//  Magnitude processor: turns Omega0 into Mw using the geometric spreading,
//  source velocity and density at the hypocentre.
// ---------------------------------------------------------------------------
class SC_SYSTEM_CLIENT_API MagnitudeProcessor_MwSpec : public Processing::MagnitudeProcessor {
	public:
		MagnitudeProcessor_MwSpec();

	public:
		void setDefaults() override {}
		bool setup(const Processing::Settings &settings) override;

		Status computeMagnitude(double amplitude, const std::string &unit,
		                        double period, double snr,
		                        double delta, double depth,
		                        const DataModel::Origin *hypocenter,
		                        const DataModel::SensorLocation *receiver,
		                        const DataModel::Amplitude *,
		                        const Locale *,
		                        double &value) override;

	private:
		MwSpecConfig _cfg;
};


}
}
}


#endif
