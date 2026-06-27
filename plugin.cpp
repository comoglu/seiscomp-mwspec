/***************************************************************************
 * SeisComP spectral moment-magnitude plugin (mwspec)                      *
 *                                                                         *
 * Plugin entry point. Registers the Mw(spec) amplitude and magnitude       *
 * processors (a port of Seisan SPEC/AUTOMAG spectral Mw).                  *
 *                                                                         *
 * GNU Affero General Public License Usage - see LICENSE.                   *
 ***************************************************************************/


#include <seiscomp/core/plugin.h>

#include "version.h"


ADD_SC_PLUGIN(
	"Spectral moment magnitude Mw(spec): per-station Brune omega-square fit of "
	"the displacement spectrum (port of Seisan SPEC/AUTOMAG).",
	"gempa/Seisan port",
	MWSPEC_VERSION_MAJOR, MWSPEC_VERSION_MINOR, MWSPEC_VERSION_PATCH
)
