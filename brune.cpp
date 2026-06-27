/***************************************************************************
 * SeisComP spectral moment-magnitude plugin (mwspec)                      *
 *                                                                         *
 * Faithful port of the Seisan SPEC/AUTOMAG Brune-model spectral fitting   *
 * and moment-magnitude computation (Lars Ottemoeller).                    *
 *                                                                         *
 * GNU Affero General Public License Usage - see LICENSE.                   *
 ***************************************************************************/


#include "brune.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <sstream>


namespace Seiscomp {
namespace Magnitudes {
namespace MwSpec {


namespace {

const double PI = 3.141592654;  // matches Seisan's literal

// Linear interpolation matching Seisan lin_interp: if the abscissae coincide,
// fall back to the lower sample (avoids division by zero at a layer boundary).
double linInterp(double x, double x1, double x2, double y1, double y2) {
	if ( x2 == x1 ) {
		return y1;
	}
	return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
}

}


// ---------------------------------------------------------------------------
SpecModel::SpecModel() {
	setDefaults();
}


void SpecModel::setDefaults() {
	// Seisan's single-layer fallback (get_att_vel: vp=6.0, vs=3.5, density=3.0,
	// no Q correction). kappa defaults follow SEISAN.DEF (0.01 P, 0.02 S).
	_layers.assign(1, SpecLayer());
	_layers[0].depth   = 0.0;
	_layers[0].vp      = 6.0;
	_layers[0].vs      = 3.5;
	_layers[0].qp      = 0.0;
	_layers[0].qap     = 0.0;
	_layers[0].qs      = 0.0;
	_layers[0].qas     = 0.0;
	_layers[0].density = 3.0;
	_kappaP     = 0.01;
	_kappaS     = 0.02;
	_qpBelow1Hz = 1.0;
	_qsBelow1Hz = 1.0;
	_qCorner    = 0.0;
}


bool SpecModel::setLayers(const std::vector<std::string> &layerDefs) {
	if ( layerDefs.empty() ) {
		return false;
	}

	std::vector<SpecLayer> layers;
	layers.reserve(layerDefs.size());

	for ( const std::string &def : layerDefs ) {
		std::istringstream iss(def);
		SpecLayer l;
		if ( !(iss >> l.depth >> l.vp >> l.vs >> l.qp >> l.qap
		           >> l.qs >> l.qas >> l.density) ) {
			return false;  // malformed line -> keep previous model
		}
		layers.push_back(l);
	}

	// Keep layers ordered by depth, as the interpolation assumes.
	std::sort(layers.begin(), layers.end(),
	          [](const SpecLayer &a, const SpecLayer &b) {
		return a.depth < b.depth;
	});

	_layers.swap(layers);
	return true;
}


void SpecModel::setKappa(double kappaP, double kappaS) {
	_kappaP = kappaP;
	_kappaS = kappaS;
}


void SpecModel::setQBelow1Hz(double qpBelow1Hz, double qsBelow1Hz) {
	_qpBelow1Hz = qpBelow1Hz;
	_qsBelow1Hz = qsBelow1Hz;
}


SourceParams SpecModel::paramsAt(double depthKm, char phase) const {
	SourceParams p;
	const bool isP = (phase == 'P');

	p.kappa      = isP ? _kappaP : _kappaS;
	p.qBelow1Hz  = isP ? _qpBelow1Hz : _qsBelow1Hz;
	p.qcorner    = _qCorner;

	const size_t n = _layers.size();

	// Single layer: return it directly (Seisan get_att_vel nspec_model==1).
	if ( n == 1 ) {
		const SpecLayer &l = _layers[0];
		p.vp = l.vp; p.vs = l.vs; p.density = l.density;
		p.q0     = isP ? l.qp  : l.qs;
		p.qalpha = isP ? l.qap : l.qas;
		return p;
	}

	// Find the first layer whose reference depth >= source depth (Seisan loop).
	size_t k = n - 1;
	for ( size_t i = 0; i < n; ++i ) {
		if ( depthKm <= _layers[i].depth ) {
			k = i;
			break;
		}
	}
	size_t j = (k > 0) ? k - 1 : 0;

	// Below the deepest layer: clamp to it. Above the shallowest: clamp to it.
	if ( depthKm > _layers[n - 1].depth ) {
		k = n - 1; j = k;
	}
	if ( depthKm < _layers[0].depth ) {
		k = 0; j = 0;
	}

	const SpecLayer &lj = _layers[j];
	const SpecLayer &lk = _layers[k];

	p.density = linInterp(depthKm, lj.depth, lk.depth, lj.density, lk.density);
	p.vs      = linInterp(depthKm, lj.depth, lk.depth, lj.vs, lk.vs);
	p.vp      = linInterp(depthKm, lj.depth, lk.depth, lj.vp, lk.vp);

	if ( isP ) {
		p.q0     = linInterp(depthKm, lj.depth, lk.depth, lj.qp, lk.qp);
		p.qalpha = linInterp(depthKm, lj.depth, lk.depth, lj.qap, lk.qap);
	}
	else {
		p.q0     = linInterp(depthKm, lj.depth, lk.depth, lj.qs, lk.qs);
		p.qalpha = linInterp(depthKm, lj.depth, lk.depth, lj.qas, lk.qas);
	}

	return p;
}


// ---------------------------------------------------------------------------
double evalQ(double q0, double f, double qalpha, double qcorner) {
	// Seisan eval_q.
	if ( qcorner == 0.0 ) {
		return q0 * std::pow(f, qalpha);
	}
	return q0 * (1.0 + std::pow(f / qcorner, qalpha));
}


// ---------------------------------------------------------------------------
double geoSpreadingDistance(char phase, double epicentralKm, double depthKm,
                            double geoDepth1, double geoDepth2,
                            double herkijDistanceKm) {
	// Port of Seisan spec_dist, returning R = 1/factor (the equivalent
	// spreading distance the moment formula multiplies by).
	const double degKm = 111.2;

	// Teleseismic (> 30 deg): factor = 0.0048 / (30 + delta_deg).
	if ( epicentralKm >= 30.0 * degKm ) {
		double factor = 0.0048 / (30.0 + epicentralKm / degKm);
		if ( factor <= 0.0 ) {
			return std::numeric_limits<double>::max();
		}
		return 1.0 / factor;
	}

	double sdistance = std::sqrt(depthKm * depthKm + epicentralKm * epicentralKm);
	if ( sdistance < 1.0e-5 ) {
		sdistance = 1.0e-5;
	}

	char p = (phase == 'P' || phase == 'S') ? phase : 'S';

	// P wave: always hypocentral spreading.
	if ( p == 'P' ) {
		return sdistance;
	}

	// S wave near the source or deep: hypocentral spreading.
	if ( epicentralKm < herkijDistanceKm || depthKm >= geoDepth2 ) {
		return sdistance;
	}

	// S wave, regional distance, shallow: Herrmann-Kijko surface spreading.
	if ( epicentralKm >= herkijDistanceKm && depthKm < geoDepth1 ) {
		return std::sqrt(herkijDistanceKm * epicentralKm);
	}

	// S wave, regional distance, transitional depth: interpolate between
	// Herrmann-Kijko and body-wave spreading.
	if ( epicentralKm >= herkijDistanceKm &&
	     depthKm >= geoDepth1 && depthKm < geoDepth2 ) {
		double frac = (depthKm - geoDepth1) / (geoDepth2 - geoDepth1);
		double factor = (1.0 - frac) * std::sqrt(epicentralKm * herkijDistanceKm)
		                + frac * sdistance;
		return factor;  // already a distance (Seisan inverts factor=1/factor)
	}

	return sdistance;
}


// ---------------------------------------------------------------------------
SNBand selectSNBand(const std::vector<double> &farray,
                    const std::vector<double> &logSignal,
                    const std::vector<double> &logNoise) {
	// Port of the signal-to-noise band logic embedded in Seisan get_om_f0.
	// All arrays are assumed ascending in frequency and the same length.
	SNBand band;
	const int n = static_cast<int>(std::min({farray.size(), logSignal.size(),
	                                          logNoise.size()}));
	if ( n < 3 ) {
		return band;
	}

	// (1) Peak log10 S/N; reject if the signal never reaches ~2.5x the noise.
	double mx = -std::numeric_limits<double>::max();
	for ( int i = 0; i < n; ++i ) {
		const double d = logSignal[i] - logNoise[i];
		if ( d > mx ) {
			mx = d;
		}
	}
	band.maxSNRLog10 = mx;
	if ( mx <= 0.4 ) {
		return band;
	}

	// (2) Spectral minimum -> high-frequency edge.
	int maxIdx = -1;
	double lo = std::numeric_limits<double>::max();
	for ( int i = 0; i < n; ++i ) {
		if ( logSignal[i] < lo ) {
			lo = logSignal[i];
			maxIdx = i;
		}
	}
	if ( maxIdx < 0 ) {
		return band;
	}
	if ( maxIdx == 0 ) {
		maxIdx = 1;
	}

	// Refine: do not trust frequencies where the S/N has decayed for good.
	{
		int k = 0;
		for ( int i = maxIdx / 2; i <= maxIdx; ++i ) {
			const double d = logSignal[i] - logNoise[i];
			if ( d < 0.3 ) {
				++k;
			}
			else {
				--k;
				if ( k < 0 ) {
					k = 0;
				}
			}
			if ( k > 2 ) {
				maxIdx = i;
				break;
			}
		}
	}

	// (3) Low-frequency edge: first sustained good-S/N point.
	int ind = -1;
	int cnt = 0;
	for ( int i = 0; i < n; ++i ) {
		const double d = logSignal[i] - logNoise[i];
		const bool good = (d - mx / 2.0 > 0.0) || (d > 1.0);
		if ( good && ind < 0 ) {
			ind = i;
			++cnt;
		}
		else if ( ind >= 0 && good ) {
			++cnt;
		}
		if ( ind >= 0 && (i - ind > 5) && cnt < 2 ) {
			ind = -1;
			cnt = 0;
		}
	}
	if ( ind >= 0 && (logSignal[ind] - logNoise[ind] < 0.3) && ind < maxIdx ) {
		ind = -1;
	}
	if ( ind < 0 || cnt < 5 ) {
		return band;
	}
	int minIdx = ind;

	if ( minIdx > maxIdx ) {
		maxIdx = n - 1;
	}

	// (6) Reject bands narrower than a factor of 2 (0.3 log units).
	if ( std::log10(farray[maxIdx]) - std::log10(farray[minIdx]) < 0.3 ) {
		return band;
	}

	// (7) Trim to at most 3 log units of bandwidth.
	if ( std::log10(farray[maxIdx]) - std::log10(farray[minIdx]) > 3.0 ) {
		for ( int i = minIdx + 1; i <= maxIdx; ++i ) {
			if ( std::log10(farray[i]) - std::log10(farray[minIdx]) > 3.0 ) {
				maxIdx = i;
				break;
			}
		}
	}

	// (8) Mean S/N over the band must be appreciable.
	double res = 0.0;
	const int m = maxIdx - minIdx + 1;
	for ( int i = minIdx; i <= maxIdx; ++i ) {
		res += (logSignal[i] - logNoise[i]) / m;
	}
	if ( res <= 0.2 ) {
		return band;
	}

	band.ok     = true;
	band.minIdx = static_cast<size_t>(minIdx);
	band.maxIdx = static_cast<size_t>(maxIdx);
	band.fmin   = farray[minIdx];
	band.fmax   = farray[maxIdx];
	return band;
}


// ---------------------------------------------------------------------------
double bruneMisfit(const std::vector<double> &farray,
                   const std::vector<double> &logSpec,
                   double fmin, double fmax,
                   double f0, double om, double k, double norm) {
	// Port of Seisan eval_spec_fit1 (returns the misfit itself; Seisan returns
	// 30 - misfit so its grid search maximizes).
	double res = 0.0;
	int cnt = 0;
	const size_t n = std::min(farray.size(), logSpec.size());

	for ( size_t i = 0; i < n; ++i ) {
		const double f = farray[i];
		if ( f >= fmin && f <= fmax && f0 > 0.0 && f0 < fmax ) {
			double x = om - std::log10(1.0 + (f / f0) * (f / f0));  // source
			x = x - PI * k * f;                                     // delta-kappa
			++cnt;
			res += std::pow(std::fabs(x - logSpec[i]), norm);
		}
	}

	if ( cnt > 0 ) {
		return std::pow(res / static_cast<double>(cnt), 1.0 / norm);
	}
	return std::numeric_limits<double>::max();
}


// ---------------------------------------------------------------------------
BruneFit bruneGridSearch(const std::vector<double> &farray,
                         const std::vector<double> &logSpec,
                         double fmin, double fmax,
                         const BruneFitOptions &opt) {
	BruneFit fit;
	fit.fmin = fmin;
	fit.fmax = fmax;

	if ( farray.empty() || farray.size() != logSpec.size() ||
	     fmax <= fmin || opt.ngridF < 2 || opt.ngridOm < 2 || opt.nloop < 1 ) {
		return fit;
	}

	const int xmax = opt.ngridF;
	const int ymax = opt.ngridOm;
	const int nkappa = (opt.dkappa != 0.0) ? 3 : 1;

	double best = -std::numeric_limits<double>::max();
	double f0 = 0.0, om = 0.0, ka = 0.0;

	for ( int loop = 1; loop <= opt.nloop; ++loop ) {
		const double loop2 = static_cast<double>(loop) * loop;
		const double xstep = (fmax - fmin) / xmax / loop2;
		const double ystep = 15.0 / ymax / loop2;

		double xstart, ystart;
		if ( loop == 1 ) {
			xstart = fmin;
			ystart = -3.0;
		}
		else {
			xstart = f0 - (xmax * xstep / 2.0);
			if ( xstart <= fmin ) {
				xstart = fmin;
			}
			ystart = om - (ymax * ystep / 2.0);
		}

		for ( int l = 1; l <= nkappa; ++l ) {
			const double k = static_cast<double>(l - 2) * opt.dkappa;
			for ( int i = 1; i <= xmax; ++i ) {
				const double x = xstart + (i - 1) * xstep;  // corner frequency
				for ( int jj = 1; jj <= ymax; ++jj ) {
					const double y = ystart + (jj - 1) * ystep;  // flat level
					// Seisan maximizes (30 - misfit); equivalently minimize misfit.
					const double misfit = bruneMisfit(farray, logSpec,
					                                  fmin, fmax, x, y, k, opt.norm);
					const double score = -misfit;
					if ( score > best ) {
						best = score;
						f0 = x;
						om = y;
						ka = k;
					}
				}
			}
		}
	}

	fit.ok          = (best > -std::numeric_limits<double>::max());
	fit.omega0Log10 = om;
	fit.cornerFreq  = f0;
	fit.deltaKappa  = ka;
	fit.residual    = -best;  // the misfit at the best node
	return fit;
}


// ---------------------------------------------------------------------------
MomentResult momentFromOmega(double omega0, double cornerFreq,
                             double velocity, double density,
                             double geoDistKm,
                             double radiation, double freeSurface) {
	MomentResult r;

	if ( omega0 <= 0.0 || velocity <= 0.0 || density <= 0.0 ||
	     geoDistKm <= 0.0 || radiation <= 0.0 || freeSurface <= 0.0 ) {
		return r;
	}

	// Seisan automag.for, generalised so radiation/freeSurface are explicit
	// (Seisan hard-codes 0.833 = 1/(2.0*0.6)).
	const double kk = 1.0 / (freeSurface * radiation);

	double m0 = 4.0 * PI * density * omega0 * kk * geoDistKm
	            * std::pow(1000.0, 5)   // unit conversion (density, v^3, R)
	            / 1.0e9;                 // omega0 nm·s -> m·s
	m0 *= velocity * velocity * velocity;

	r.m0    = m0;
	r.logM0 = std::log10(m0);
	r.mw    = 2.0 / 3.0 * r.logM0 - 6.06;

	// Brune source dimension and stress drop (Seisan automag.for), with the
	// radius reported in metres.
	if ( cornerFreq > 0.0 ) {
		const double radiusKm = 0.37 * velocity / cornerFreq;
		r.sourceRadius = radiusKm * 1000.0;
		if ( radiusKm > 0.0 ) {
			r.stressDrop = (0.44 * m0) / (1.0e14 * radiusKm * radiusKm * radiusKm);
		}
	}

	return r;
}


}
}
}
