/***************************************************************************
 * SeisComP spectral moment-magnitude plugin (mwspec)                      *
 *                                                                         *
 * Faithful port of the Seisan SPEC/AUTOMAG Brune-model spectral fitting   *
 * and moment magnitude computation (Lars Ottemoeller).                    *
 *                                                                         *
 * References in the Seisan source tree:                                   *
 *   PRO/automag.for       - moment/Mw block (om, f0 -> M0 -> Mw)          *
 *   LIB/auto_proc.for      - get_om_f0, grid_om_fo, eval_spec_fit1, eval_q*
 *   LIB/spec_dist.for       - geometric spreading                         *
 *   LIB/libsei.for          - get_att_vel (layered model interpolation)   *
 *                                                                         *
 * GNU Affero General Public License Usage - see LICENSE.                   *
 ***************************************************************************/


#ifndef SEISCOMP_MAGNITUDES_MWSPEC_BRUNE_H
#define SEISCOMP_MAGNITUDES_MWSPEC_BRUNE_H


#include <string>
#include <vector>


namespace Seiscomp {
namespace Magnitudes {
namespace MwSpec {


/**
 * One layer of the Seisan spectral attenuation/velocity model.
 * Mirrors a "SPEC MODEL h,vp,vs,qp,qap,qs,qas,d" line in SEISAN.DEF.
 * All velocities in km/s, density in g/cm^3, depth in km.
 */
struct SpecLayer {
	double depth   = 0.0;   //!< layer reference depth [km]
	double vp      = 6.0;   //!< P velocity [km/s]
	double vs      = 3.5;   //!< S velocity [km/s]
	double qp      = 0.0;   //!< Q0 for P (0 => no Q correction)
	double qap     = 0.0;   //!< Q alpha for P
	double qs      = 0.0;   //!< Q0 for S
	double qas     = 0.0;   //!< Q alpha for S
	double density = 2.6;   //!< density [g/cm^3]
};


/**
 * Phase- and depth-dependent material/attenuation parameters,
 * the result of interpolating the layered model at the source depth.
 */
struct SourceParams {
	double vp       = 6.0;
	double vs       = 3.5;
	double density  = 2.6;
	double q0       = 0.0;   //!< Q0 for the requested phase
	double qalpha   = 0.0;   //!< Q alpha for the requested phase
	double qcorner  = 0.0;   //!< Q corner frequency (0 => power law Q0*f^alpha)
	double kappa    = 0.0;   //!< near-surface attenuation for the phase
	double qBelow1Hz = 1.0;  //!< constant Q below 1 Hz (informational)

	//! Returns vp for phase 'P', otherwise vs.
	double velocity(char phase) const { return phase == 'P' ? vp : vs; }
};


/**
 * The Seisan spectral model: a stack of depth layers plus the global
 * kappa and "Q below 1 Hz" constants. Loaded from configuration; falls
 * back to Seisan's single-layer defaults (vp=6.0, vs=3.5, density=3.0).
 */
class SpecModel {
	public:
		SpecModel();

		//! Resets to Seisan's hard-coded single-layer default model.
		void setDefaults();

		/**
		 * Parses "SPEC MODEL" style layer definitions. Each entry is a
		 * whitespace-separated "depth vp vs qp qap qs qas density".
		 * Returns false (and keeps previous state) on a malformed entry.
		 */
		bool setLayers(const std::vector<std::string> &layerDefs);

		void setKappa(double kappaP, double kappaS);
		void setQBelow1Hz(double qpBelow1Hz, double qsBelow1Hz);
		void setQCorner(double qCorner) { _qCorner = qCorner; }

		//! Number of configured layers (always >= 1 after construction).
		size_t layerCount() const { return _layers.size(); }

		/**
		 * Interpolates the layered model at the given source depth [km]
		 * for the given phase ('P' or 'S'), exactly as Seisan get_att_vel:
		 * linear interpolation in depth, clamped to the layer range.
		 */
		SourceParams paramsAt(double depthKm, char phase) const;

	private:
		std::vector<SpecLayer> _layers;
		double _kappaP      = 0.01;
		double _kappaS      = 0.02;
		double _qpBelow1Hz  = 1.0;
		double _qsBelow1Hz  = 1.0;
		double _qCorner     = 0.0;
};


/**
 * Seisan eval_q: frequency-dependent quality factor.
 *   qcorner == 0 : Q(f) = q0 * f^qalpha
 *   qcorner != 0 : Q(f) = q0 * (1 + (f/qcorner)^qalpha)
 */
double evalQ(double q0, double f, double qalpha, double qcorner);


/**
 * Seisan spec_dist: equivalent geometric-spreading distance R [km] for the
 * moment computation (i.e. 1/spreading-factor). For P (and near/deep S) this
 * is the hypocentral distance; for regional S beyond herkijDistance it follows
 * the Herrmann-Kijko surface-wave spreading, interpolated across the
 * geoDepth1..geoDepth2 band; beyond 30 degrees a teleseismic approximation.
 *
 * @param epicentralKm  epicentral distance [km]
 * @param depthKm       hypocentral depth [km]
 * Returns the equivalent distance R [km] (>0), or a very large value if the
 * spreading factor is ~0 (teleseismic), which the caller must guard.
 */
double geoSpreadingDistance(char phase, double epicentralKm, double depthKm,
                            double geoDepth1, double geoDepth2,
                            double herkijDistanceKm);


/** Result of a Brune ω² spectral fit. */
struct BruneFit {
	bool   ok          = false;
	double omega0Log10 = 0.0;  //!< log10 of the displacement flat level (Seisan "om")
	double cornerFreq  = 0.0;  //!< corner frequency f0 [Hz]
	double deltaKappa  = 0.0;  //!< fitted delta-kappa (0 unless dkappa search enabled)
	double residual    = 0.0;  //!< L-norm misfit (lower is better)
	double fmin        = 0.0;  //!< low edge of the frequency band actually fitted [Hz]
	double fmax        = 0.0;  //!< high edge of the frequency band actually fitted [Hz]
};


/** Result of the signal-to-noise frequency-band selection. */
struct SNBand {
	bool   ok      = false;
	size_t minIdx  = 0;    //!< index into farray of the low-frequency edge
	size_t maxIdx  = 0;    //!< index into farray of the high-frequency edge
	double fmin    = 0.0;  //!< low-frequency edge [Hz]
	double fmax    = 0.0;  //!< high-frequency edge [Hz]
	double maxSNRLog10 = 0.0; //!< peak log10(signal/noise) over the band
};


/**
 * Selects the usable frequency band from the (attenuation-corrected) log10
 * signal and noise spectral levels, porting the S/N logic of Seisan
 * get_om_f0: requires the peak log10 S/N > 0.4, locates the spectral minimum
 * as the high-frequency edge, walks up from the low frequencies to the first
 * sustained good-S/N point, and enforces the band-width sanity checks (>= 0.3
 * log units wide, trimmed to <= 3 log units). Returns ok=false if no usable
 * band exists.
 */
SNBand selectSNBand(const std::vector<double> &farray,
                    const std::vector<double> &logSignal,
                    const std::vector<double> &logNoise);


/** Tunables for the grid search (Seisan defaults: 100, 100, 5, norm 1). */
struct BruneFitOptions {
	int    ngridF = 100;  //!< grid nodes along corner frequency
	int    ngridOm = 100; //!< grid nodes along flat level
	int    nloop  = 5;    //!< zoom iterations (step shrinks ~1/loop^2)
	double dkappa = 0.0;  //!< if >0, also search delta-kappa in {-dkappa,0,+dkappa}
	double norm   = 1.0;  //!< residual norm (Seisan uses 1 = L1)
};


/**
 * L-norm misfit between the Brune model and the observed log10 spectrum,
 * Seisan eval_spec_fit1 (here returned directly as the misfit, i.e. lower is
 * better; Seisan returns 30 - misfit so the grid search maximizes it).
 *
 *   model(f) = om - log10(1 + (f/f0)^2) - pi*k*f
 */
double bruneMisfit(const std::vector<double> &farray,
                   const std::vector<double> &logSpec,
                   double fmin, double fmax,
                   double f0, double om, double k, double norm);


/**
 * Grid-search Brune fit, Seisan grid_om_fo. The inputs are the discrete
 * log-spaced frequency array and the matching attenuation-corrected log10
 * displacement spectrum, plus the [fmin,fmax] band determined from S/N.
 */
BruneFit bruneGridSearch(const std::vector<double> &farray,
                         const std::vector<double> &logSpec,
                         double fmin, double fmax,
                         const BruneFitOptions &opt);


/** Seismic-moment results derived from a Brune fit. */
struct MomentResult {
	double m0         = 0.0;  //!< seismic moment [N·m]
	double logM0      = 0.0;  //!< log10(M0)
	double mw         = 0.0;  //!< moment magnitude
	double sourceRadius = 0.0; //!< Brune source radius [m]
	double stressDrop = 0.0;  //!< stress drop [Pa]
};


/**
 * Seismic moment and Mw from the spectral flat level, reproducing Seisan
 * automag.for bit-for-bit (Seisan units: density g/cm^3, velocity km/s,
 * distance km, flat level nm·s):
 *
 *   M0 = 4π · ρ · Ω0 · (1/(freeSurface·radiation)) · R · v³ · 1e15/1e9   [N·m]
 *   Mw = (2/3)·log10(M0) − 6.06
 *
 * Seisan hard-codes 1/(freeSurface·radiation) = 0.833 = 1/(2·0.6); here both
 * factors are configurable. The 1e15 = (1e3 density)(1e9 v³)(1e3 R) unit
 * conversion and the 1e-9 converts Ω0 from nm·s to m·s.
 *
 * @param omega0      displacement flat level Ω0 [nm·s] (linear, NOT log10)
 * @param geoDistKm   equivalent geometric-spreading distance R [km]
 * @param density     source density [g/cm^3]
 * @param velocity    source velocity for the phase [km/s]
 */
MomentResult momentFromOmega(double omega0, double cornerFreq,
                             double velocity, double density,
                             double geoDistKm,
                             double radiation, double freeSurface);


}
}
}


#endif
