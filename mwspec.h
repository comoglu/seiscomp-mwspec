/***************************************************************************
 * SeisComP spectral moment-magnitude plugin (mwspec)                      *
 *                                                                         *
 * Per-station moment magnitude Mw(spec) from a Brune omega-square fit of   *
 * the displacement spectrum, porting Seisan SPEC/AUTOMAG.                  *
 *                                                                         *
 * Copyright (C) 2026 Mustafa Comoglu (Geoscience Australia)               *
 * GNU Affero General Public License Usage - see LICENSE.                   *
 ***************************************************************************/


#ifndef SEISCOMP_MAGNITUDES_MWSPEC_PLUGIN_H
#define SEISCOMP_MAGNITUDES_MWSPEC_PLUGIN_H


#include <seiscomp/core/plugin.h>
#include <seiscomp/core/version.h>
#include <seiscomp/datamodel/comment.h>
#include <seiscomp/processing/amplitudeprocessor.h>
#include <seiscomp/processing/magnitudeprocessor.h>

// Spectral diagnostics for the amplitude review window (SeisComP with the
// SpectralDiagnosticsProvider interface); without it the plugin builds as before.
#if __has_include(<seiscomp/processing/spectraldiagnostics.h>)
#include <seiscomp/processing/spectraldiagnostics.h>
#define MWSPEC_SPECTRAL_DIAGNOSTICS 1
#endif

#include <cstdio>
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

	//! Where the S signal window starts (amplitudes are triggered on P):
	//! the associated S pick, else the first theoretical S-type arrival
	//! (Auto); always theoretical (TravelTime); or the P trigger (Trigger,
	//! the pre-0.7 behaviour).
	enum SOnset { SOnsetAuto, SOnsetTravelTime, SOnsetTrigger };
	SOnset   sOnset           = SOnsetAuto;

	//! S only: the window ends no earlier than the Lg arrival, epicentral
	//! distance / lgVelocity + lgMargin after the origin time, so that it
	//! covers Lg at regional distances where the S onset is Sn. 0 = off.
	//! Applies up to lgMaxDistance (Lg is a regional, continental phase).
	double   lgVelocity       = 3.0;     //!< [km/s] slowest Lg group velocity
	double   lgMargin         = 10.0;    //!< [s] after the Lg arrival
	double   lgMaxDistanceDeg = 20.0;    //!< [deg]
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


/**
 * Sets comment @p id of an Amplitude or StationMagnitude to the number
 * @p value (printed with @p fmt), replacing an existing comment with the same
 * id. The fit diagnostics travel to the database/QuakeML this way, so scolv
 * and external tools can show them; see the README for the comment ids.
 */
template <typename T>
void setComment(T *obj, const std::string &id, const std::string &text) {
	DataModel::Comment *c = obj->comment(DataModel::CommentIndex(id));
	if ( c ) {
		c->setText(text);
		return;
	}
	DataModel::CommentPtr nc = new DataModel::Comment;
	nc->setId(id);
	nc->setText(text);
	obj->add(nc.get());
}

template <typename T>
void setComment(T *obj, const std::string &id, double value,
                const char *fmt = "%.4g") {
	char buf[64];
	std::snprintf(buf, sizeof(buf), fmt, value);
	setComment(obj, id, std::string(buf));
}


/**
 * Per-component spectral fit result kept by the amplitude worker after a
 * successful computeAmplitude(), for finalizeAmplitude().
 */
struct FitDiagnostics {
	bool   valid      = false;
	double omega0     = 0.0;   //!< flat level [nm*s], gain-corrected
	double cornerFreq = 0.0;   //!< [Hz]
	double fmin       = 0.0;   //!< fitted band low edge [Hz]
	double fmax       = 0.0;   //!< fitted band high edge [Hz]
	double residual   = 0.0;   //!< Brune-fit misfit
	double deltaKappa = 0.0;   //!< fitted delta-kappa
	double travelTime = 0.0;   //!< [s] used for the Q correction
};


// ---------------------------------------------------------------------------
//  Amplitude processor: builds the displacement spectrum, fits the Brune
//  model and emits Omega0 (flat level, nm*s) with the corner frequency as the
//  carried "period".
// ---------------------------------------------------------------------------
class SC_SYSTEM_CLIENT_API AmplitudeProcessor_MwSpec : public Processing::AmplitudeProcessor
#ifdef MWSPEC_SPECTRAL_DIAGNOSTICS
                                                     , public Processing::SpectralDiagnosticsProvider
#endif
{
	public:
		AmplitudeProcessor_MwSpec();

	public:
		bool setup(const Processing::Settings &settings) override;
		int capabilities() const override;
		Processing::AmplitudeProcessor::IDList
		    capabilityParameters(Capability cap) const override;
		bool setParameter(Capability cap, const std::string &value) override;

		void setEnvironment(const DataModel::Origin *hypocenter,
		                    const DataModel::SensorLocation *receiver,
		                    const DataModel::Pick *pick) override;

		void finalizeAmplitude(DataModel::Amplitude *amplitude) const override;

		//! Diagnostics of the last successful fit (valid == false otherwise).
		const FitDiagnostics &fitDiagnostics() const { return _fit; }

#ifdef MWSPEC_SPECTRAL_DIAGNOSTICS
		//! Spectra, Brune model, band and fit values of the last measurement
		//! (also when it was rejected; status says why).
		const Processing::SpectralDiagnostics *spectralDiagnostics() const override {
			return &_diag;
		}

		//! A band set by the user replaces the automatic S/N band selection
		bool canSetSpectralBand() const override { return true; }
		bool setSpectralBand(double fmin, double fmax) override;
#endif

		//! Writes the methodID, version and fit-diagnostic comments.
		//! @p suffix is appended to the comment ids (e.g. ".N" for S workers).
		void writeDiagnostics(DataModel::Amplitude *amplitude,
		                      const std::string &suffix) const;

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
		//! Sets the trigger-relative signal/noise windows (signal shifted by
		//! _signalShift). Unlike applyConfig() it leaves the component alone,
		//! which the combiner assigns per worker.
		void applyWindows();
		//! S only: seconds from the trigger (P) to the S onset, or 0.
		double sOnsetShift(const DataModel::Origin *hypocenter,
		                   const DataModel::SensorLocation *receiver,
		                   const DataModel::Pick *pick);

	private:
		MwSpecConfig   _cfg;
		FitDiagnostics _fit;
#ifdef MWSPEC_SPECTRAL_DIAGNOSTICS
		Processing::SpectralDiagnostics _diag;
#endif

		double      _signalShift = 0.0;   //!< [s] signal window offset from trigger
		double      _bandFmin = 0.0;      //!< [Hz] user fit band, 0 = automatic
		double      _bandFmax = 0.0;      //!< [Hz] user fit band
		std::string _onsetSource;         //!< "pick", "ttt" or "trigger" (S only)
		std::string _dumpDir;             //!< MWSPEC_DUMP_DIR: write spectra as JSON

		// Hypocentre values copied in setEnvironment(). The base class keeps
		// only a raw Origin pointer, and scamp does not keep a messaging
		// origin alive until our (long) window completes, so it must not be
		// dereferenced in computeAmplitude().
		double          _srcDepthKm = 0.0;
		OPT(Core::Time) _originTime;
		double          _rhypKm     = 0.0;   //!< 0 = geometry unavailable
		double          _epiKm      = 0.0;   //!< epicentral distance, 0 = unknown
		double          _signalEnd  = 0.0;   //!< [s] S: window end from the trigger, 0 = default

	friend class AmplitudeProcessor_MwSpecCombiner;
};


// ---------------------------------------------------------------------------
//  Component combiner (registered as "Mw(spec)"). For P it runs a single
//  worker on the vertical; for S it runs two workers on the horizontals and
//  combines their Omega0 (default: vector sum sqrt(N^2 + E^2) = total S-wave
//  horizontal motion). This is the entry point scamp/scolv instantiates.
// ---------------------------------------------------------------------------
class SC_SYSTEM_CLIENT_API AmplitudeProcessor_MwSpecCombiner : public Processing::AmplitudeProcessor
#ifdef MWSPEC_SPECTRAL_DIAGNOSTICS
                                                             , public Processing::SpectralDiagnosticsProvider
#endif
{
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

		void finalizeAmplitude(DataModel::Amplitude *amplitude) const override;

#ifdef MWSPEC_SPECTRAL_DIAGNOSTICS
		//! The diagnostics of the active worker(s), merged for S.
		const Processing::SpectralDiagnostics *spectralDiagnostics() const override;
		bool canSetSpectralBand() const override { return true; }
		//! Sets the band of all active workers
		bool setSpectralBand(double fmin, double fmax) override;
#endif

	protected:
		bool computeAmplitude(const DoubleArray &data,
		                      size_t i1, size_t i2,
		                      size_t si1, size_t si2,
		                      double offset,
		                      AmplitudeIndex *dt,
		                      AmplitudeValue *amplitude,
		                      double *period, double *snr) override;

	private:
		bool feedWorkers(const Record *record);

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
#ifdef MWSPEC_SPECTRAL_DIAGNOSTICS
		mutable Processing::SpectralDiagnostics _diag;
#endif
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

		void finalizeMagnitude(DataModel::StationMagnitude *magnitude) const override;

	private:
		MwSpecConfig _cfg;

		// Source parameters of the last successful computeMagnitude(). scmag
		// and scolv call finalizeMagnitude() right after it for the same
		// amplitude, which is how they reach the StationMagnitude comments.
		bool         _lastValid = false;
		MomentResult _last;
		double       _lastCornerFreq = 0.0;
		double       _lastGeoDistKm  = 0.0;
};


}
}
}


#endif
