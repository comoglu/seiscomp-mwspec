// Standalone unit tests for the framework-independent Brune/Seisan core.
// Build:  g++ -std=c++17 -I. test_brune.cpp brune.cpp -o /tmp/test_brune && /tmp/test_brune
//
// These tests check the science port for internal consistency and against
// hand-computable Seisan reference values. They do NOT touch SeisComP.

#include "../brune.h"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace Seiscomp::Magnitudes::MwSpec;

static int g_fail = 0;

static void check(bool cond, const char *name, double got = 0, double exp = 0) {
	if ( cond ) {
		printf("  ok   %s\n", name);
	}
	else {
		printf("  FAIL %s (got %.6g, expected %.6g)\n", name, got, exp);
		++g_fail;
	}
}

static bool close(double a, double b, double tol) { return std::fabs(a - b) <= tol; }


int main() {
	printf("== evalQ ==\n");
	// Q(f) = q0 * f^qalpha for qcorner==0
	check(close(evalQ(500, 1.0, 0.7, 0.0), 500.0, 1e-9), "Q at 1Hz = q0");
	check(close(evalQ(500, 4.0, 0.7, 0.0), 500 * std::pow(4.0, 0.7), 1e-6), "Q power law");
	// Q(f) = q0*(1+(f/fc)^qalpha) for qcorner!=0
	check(close(evalQ(500, 2.0, 1.0, 2.0), 500 * (1 + 1.0), 1e-6), "Q corner form");

	printf("== geoSpreadingDistance ==\n");
	// P: always hypocentral distance sqrt(epi^2 + depth^2)
	check(close(geoSpreadingDistance('P', 100, 10, 50, 100, 100),
	            std::sqrt(100.0 * 100 + 10 * 10), 1e-6), "P hypocentral");
	// S, shallow, regional: Herrmann-Kijko sqrt(herkij*epi)
	check(close(geoSpreadingDistance('S', 200, 5, 50, 100, 100),
	            std::sqrt(100.0 * 200), 1e-6), "S Herrmann-Kijko");
	// S, near (< herkij): hypocentral
	check(close(geoSpreadingDistance('S', 50, 5, 50, 100, 100),
	            std::sqrt(50.0 * 50 + 25), 1e-6), "S near = hypocentral");
	// Teleseismic (> 30 deg): R = 1/factor, factor=0.0048/(30+deg)
	{
		double epi = 40 * 111.2, deg = 40.0;
		double expR = 1.0 / (0.0048 / (30 + deg));
		check(close(geoSpreadingDistance('P', epi, 10, 50, 100, 100), expR, 1e-3),
		      "teleseismic spreading");
	}

	printf("== SpecModel ==\n");
	{
		SpecModel m;  // defaults: single layer vp=6, vs=3.5, density=3.0
		SourceParams p = m.paramsAt(10.0, 'P');
		check(close(p.vp, 6.0, 1e-9) && close(p.density, 3.0, 1e-9), "default single layer");
		check(close(p.kappa, 0.01, 1e-9), "default kappa P");
		SourceParams ps = m.paramsAt(10.0, 'S');
		check(close(ps.kappa, 0.02, 1e-9), "default kappa S");

		// Two-layer interpolation (SEISAN.DEF style lines).
		std::vector<std::string> layers = {
			"3.0 5.8 3.2 500 0.7 400 0.7 2.6",
			"15.0 6.8 3.9 500 0.7 400 0.7 2.9"
		};
		SpecModel m2;
		bool okLayers = m2.setLayers(layers);
		check(okLayers, "setLayers parse");
		SourceParams q = m2.paramsAt(9.0, 'P');  // halfway between 3 and 15
		check(close(q.vp, 0.5 * (5.8 + 6.8), 1e-6), "vp interpolation midpoint", q.vp);
		check(close(q.density, 0.5 * (2.6 + 2.9), 1e-6), "density interpolation", q.density);
	}

	printf("== bruneGridSearch (recover synthetic source) ==\n");
	{
		const double om0 = 2.3;   // log10 flat level (nm*s)
		const double fc0 = 1.7;   // corner freq Hz
		std::vector<double> farray, logSpec;
		const int n = 100;
		double flo = 0.1, fhi = 20.0;
		for ( int i = 0; i < n; ++i ) {
			double f = std::pow(10.0, std::log10(flo) +
			                    (std::log10(fhi) - std::log10(flo)) * i / (n - 1));
			farray.push_back(f);
			logSpec.push_back(om0 - std::log10(1.0 + (f / fc0) * (f / fc0)));
		}
		BruneFitOptions opt;  // Seisan defaults 100/100/5
		BruneFit fit = bruneGridSearch(farray, logSpec, flo, fhi, opt);
		check(fit.ok, "fit ok");
		check(close(fit.omega0Log10, om0, 0.05), "recover om0", fit.omega0Log10, om0);
		check(close(fit.cornerFreq, fc0, 0.15), "recover fc", fit.cornerFreq, fc0);
		check(fit.residual < 0.02, "low residual on synthetic", fit.residual, 0);
	}

	printf("== selectSNBand ==\n");
	{
		const int n = 100;
		std::vector<double> farray, logSig, logNoise;
		double flo = 0.1, fhi = 20.0;
		for ( int i = 0; i < n; ++i ) {
			double f = std::pow(10.0, std::log10(flo) +
			                    (std::log10(fhi) - std::log10(flo)) * i / (n - 1));
			farray.push_back(f);
			// Brune signal well above a flat noise floor in the mid band.
			logSig.push_back(2.0 - std::log10(1.0 + (f / 2.0) * (f / 2.0)));
			logNoise.push_back(-2.0);
		}
		SNBand band = selectSNBand(farray, logSig, logNoise);
		check(band.ok, "band found");
		check(band.fmin < band.fmax && band.fmin >= flo, "band ordered/in range");
	}

	printf("== momentFromOmega (Seisan formula, exact) ==\n");
	{
		// Reproduce Seisan automag arithmetic by hand:
		// M0 = 4*pi*rho*10^om * (1/(fs*rad)) * R * 1000^5/1e9 * v^3
		double om = 2.0, rho = 2.6, R = 100.0, v = 6.0, rad = 0.6, fs = 2.0;
		double omega0 = std::pow(10.0, om);  // linear nm*s
		double kk = 1.0 / (fs * rad);
		double expM0 = 4.0 * 3.141592654 * rho * omega0 * kk * R
		               * std::pow(1000.0, 5) / 1e9 * v * v * v;
		double expMw = 2.0 / 3.0 * std::log10(expM0) - 6.06;

		MomentResult mr = momentFromOmega(omega0, 1.7, v, rho, R, rad, fs);
		check(close(mr.m0, expM0, expM0 * 1e-9), "M0 matches Seisan formula", mr.m0, expM0);
		check(close(mr.mw, expMw, 1e-6), "Mw matches", mr.mw, expMw);
		// Brune radius r = 0.37*v/fc (km) -> metres
		check(close(mr.sourceRadius, 0.37 * v / 1.7 * 1000.0, 1e-3), "source radius");
		printf("    [info] om=2.0 nm*s, R=100km, v=6 -> M0=%.3e Nm, Mw=%.2f\n",
		       mr.m0, mr.mw);
	}

	printf("\n%s (%d failure%s)\n", g_fail ? "TESTS FAILED" : "ALL TESTS PASSED",
	       g_fail, g_fail == 1 ? "" : "s");
	return g_fail ? 1 : 0;
}
