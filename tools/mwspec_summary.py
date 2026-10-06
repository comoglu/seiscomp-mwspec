#!/usr/bin/env python3
############################################################################
# Mw(spec) origin summary for the scolv origin panel                       #
#                                                                          #
# scolv runs the scripts listed in display.origin.addons for each origin,  #
# writes the Origin (SeisComP binary archive) to stdin and shows the first #
# line printed on stdout. This prints the network Mw(spec) and the median  #
# source parameters of the contributing station magnitudes, read from the  #
# plugin's StationMagnitude comments (M0, fc, stressDrop).                 #
#                                                                          #
# Copyright (C) 2026 Mustafa Comoglu (Geoscience Australia)                #
# GNU Affero General Public License Usage - see LICENSE.                   #
############################################################################

import os
import statistics
import sys
import tempfile

import seiscomp.client
import seiscomp.datamodel as dm
import seiscomp.io

TYPE = "Mw(spec)"


def read_origin(data):
    # BinaryArchive reads from a file only
    with tempfile.NamedTemporaryFile(suffix=".bin") as f:
        f.write(data)
        f.flush()
        ar = seiscomp.io.BinaryArchive()
        if not ar.open(f.name):
            return None
        obj = ar.readObject()
        ar.close()
    return dm.Origin.Cast(obj)


def comments(obj):
    out = {}
    for i in range(obj.commentCount()):
        c = obj.comment(i)
        try:
            out[c.id()] = float(c.text())
        except ValueError:
            pass
    return out


def network_magnitude(origin):
    for i in range(origin.magnitudeCount()):
        mag = origin.magnitude(i)
        if mag.type() == TYPE:
            return mag
    return None


def summarize(origin, stamags, mag):
    """stamags: list of (publicID, {comment id: value})."""
    used = None
    if mag is not None and mag.stationMagnitudeContributionCount() > 0:
        used = set()
        for i in range(mag.stationMagnitudeContributionCount()):
            c = mag.stationMagnitudeContribution(i)
            try:
                if c.weight() <= 0:
                    continue
            except ValueError:
                pass
            used.add(c.stationMagnitudeID())
    values = [v for pid, v in stamags if used is None or pid in used]

    def median(key):
        xs = [v[key] for v in values if key in v]
        return (statistics.median(xs), xs) if xs else (None, [])

    parts = []
    if mag is not None:
        n = mag.stationCount() if _has(mag.stationCount) else len(values)
        parts.append(f"{TYPE} {mag.magnitude().value():.2f} ({n} sta)")
    fc, _ = median("fc")
    if fc is not None:
        parts.append(f"fc {fc:.2g} Hz")
    sd, sds = median("stressDrop")
    if sd is not None:
        if len(sds) >= 4:
            q = statistics.quantiles(sds, n=4)
            parts.append(f"Δσ {sd:.2g} MPa [{q[0]:.2g}–{q[2]:.2g}]")
        else:
            parts.append(f"Δσ {sd:.2g} MPa")
    m0, _ = median("M0")
    if m0 is not None:
        parts.append(f"M0 {m0:.2e} N·m")
    return " · ".join(parts) if parts else "-"


def _has(getter):
    try:
        getter()
        return True
    except ValueError:
        return False


def from_memory(origin):
    stamags = []
    for i in range(origin.stationMagnitudeCount()):
        sm = origin.stationMagnitude(i)
        if sm.type() == TYPE:
            stamags.append((sm.publicID(), comments(sm)))
    return stamags


class DatabaseReader(seiscomp.client.Application):
    """Loads the Mw(spec) station magnitudes (with comments) of an origin."""

    def __init__(self, origin_id):
        super().__init__(1, [sys.argv[0]])
        self.setMessagingEnabled(False)
        self.setDatabaseEnabled(True, False)
        self.setLoggingToStdErr(False)
        self._origin_id = origin_id
        self.origin = None
        self.stamags = []

    def run(self):
        q = self.query()
        self.origin = dm.Origin.Cast(q.getObject(dm.Origin.TypeInfo(), self._origin_id))
        if self.origin is None:
            return True
        q.loadMagnitudes(self.origin)
        for i in range(self.origin.magnitudeCount()):
            q.loadStationMagnitudeContributions(self.origin.magnitude(i))
        q.loadStationMagnitudes(self.origin)
        for i in range(self.origin.stationMagnitudeCount()):
            sm = self.origin.stationMagnitude(i)
            if sm.type() == TYPE:
                q.loadComments(sm)
                self.stamags.append((sm.publicID(), comments(sm)))
        return True


def main():
    sys.stdout.reconfigure(encoding="utf-8")  # scolv reads the line as UTF-8
    origin = read_origin(sys.stdin.buffer.read())
    if origin is None:
        print("-")
        return 0

    stamags = from_memory(origin)
    mag = network_magnitude(origin)

    # scolv usually has the station magnitudes but not their comments; fall
    # back to the database (committed origins) in that case.
    if not any(v for _, v in stamags):
        app = DatabaseReader(origin.publicID())
        app()
        if app.origin is not None and app.stamags:
            stamags = app.stamags
            mag = network_magnitude(app.origin) or mag

    print(summarize(origin, stamags, mag))
    return 0


if __name__ == "__main__":
    sys.exit(main())
