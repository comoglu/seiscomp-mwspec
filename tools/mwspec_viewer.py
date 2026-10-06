#!/usr/bin/env python3
############################################################################
# Mw(spec) spectrum viewer — scolv commandMenuAction popup                 #
#                                                                          #
# Launched by scolv as  mwspec_viewer.py <OriginID> [<EventID>]. Recomputes #
# Mw(spec) for the origin with the installed plugin (scamp + scmag offline, #
# Mw(spec) only) and shows, per station, the spectra the plugin actually   #
# fitted: Q/kappa-corrected and raw displacement spectrum, pre-P noise,    #
# the fitted Brune model, the fitted band and fc. The spectra come from    #
# the plugin itself (MWSPEC_DUMP_DIR), so nothing is re-implemented here.  #
#                                                                          #
#   scolv.cfg:                                                             #
#     olv.commandMenuAction.mwspec.enable  = true                          #
#     olv.commandMenuAction.mwspec.command = @DATADIR@/client/mwspec_viewer.py
#     olv.commandMenuAction.mwspec.text    = "Mw(spec) spectra"            #
#                                                                          #
# Waveforms are read with the launching scolv's recordstream (its -I, else #
# its config) unless -I is given here.                                     #
# Headless check: mwspec_viewer.py <OriginID> --png out.png [--ep x.xml]   #
#                                                                          #
# Copyright (C) 2026 Mustafa Comoglu (Geoscience Australia)                #
# GNU Affero General Public License Usage - see LICENSE.                   #
############################################################################

import argparse
import glob
import json
import math
import os
import shutil
import statistics
import subprocess
import sys
import tempfile
import traceback
import xml.etree.ElementTree as ET

import seiscomp.config
import seiscomp.system
from seiscomp.client import Application

TYPE = "Mw(spec)"
COMP_COLORS = {"Z": "#1f77b4", "N": "#d62728", "E": "#2ca02c"}
SEISCOMP = os.path.join(os.environ.get("SEISCOMP_ROOT", "/opt/seiscomp"), "bin", "seiscomp")


# ---------------------------------------------------------------------------
#  Computation (no GUI)
# ---------------------------------------------------------------------------

class StationResult:
    def __init__(self, sid):
        self.sid = sid               # NET.STA
        self.distance_km = None
        self.dumps = []              # plugin JSON dumps, one per component
        self.mw = None
        self.used = False            # contributes to the network magnitude
        self.amp = {}                # amplitude comments
        self.mag = {}                # station magnitude comments
        self.snr = None

    @property
    def status(self):
        if self.mw is not None:
            return "used" if self.used else "not used"
        reasons = sorted({d["status"] for d in self.dumps if d["status"] != "ok"})
        return ", ".join(reasons) if reasons else "no result"

    @property
    def onset(self):
        return ",".join(sorted({d["onset"] for d in self.dumps if d.get("onset")})) or "-"


def _parent_args():
    """argv of the launching process (scolv starts the command directly)."""
    try:
        with open(f"/proc/{os.getppid()}/cmdline", "rb") as f:
            return [a.decode("utf-8", "replace") for a in f.read().split(b"\0") if a]
    except OSError:
        return []


def scolv_recordstream():
    """The recordstream of the scolv instance that launched us: its -I /
    --record-url argument, else the config of that scolv (or its alias)."""
    name = "scolv"
    args = _parent_args()
    if args and os.path.basename(args[0]).startswith("scolv"):
        name = os.path.basename(args[0])
        for i, a in enumerate(args[1:], 1):
            if a in ("-I", "--record-url") and i + 1 < len(args):
                return args[i + 1]
            if a.startswith("--record-url="):
                return a.split("=", 1)[1]
            if a.startswith("-I") and len(a) > 2:
                return a[2:]
    cfg = seiscomp.config.Config()
    seiscomp.system.Environment.Instance().initConfig(cfg, name)
    try:
        return cfg.getString("recordstream")
    except Exception:
        return ""


def run(cmd, env=None, log=print):
    log("$ " + " ".join(cmd[2:] if cmd[:2] == [SEISCOMP, "exec"] else cmd))
    p = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8",
                       errors="replace", env=env, timeout=1800)
    if p.returncode != 0:
        raise RuntimeError(f"{cmd[2] if len(cmd) > 2 else cmd[0]} failed "
                           f"(exit {p.returncode}):\n{p.stderr[-2000:]}")
    return p.stdout


def _ns(root):
    return root.tag.split("}")[0] + "}" if root.tag.startswith("{") else ""


def _comments(el, ns):
    out = {}
    for c in el.findall(ns + "comment"):
        cid, text = c.findtext(ns + "id"), c.findtext(ns + "text")
        try:
            out[cid] = float(text)
        except (TypeError, ValueError):
            out[cid] = text
    return out


def keep_only_origin(path, origin_id):
    """scamp --ep processes every origin in the file: keep only ours."""
    tree = ET.parse(path)
    root = tree.getroot()
    ns = _ns(root)
    if ns:
        ET.register_namespace("", ns[1:-1])
    for ep in root.iter(ns + "EventParameters"):
        for child in list(ep):
            if child.tag == ns + "event" or \
               (child.tag == ns + "origin" and child.get("publicID") != origin_id):
                ep.remove(child)
    tree.write(path, xml_declaration=True, encoding="UTF-8")


def compute(origin_id, database, recordstream, workdir, ep=None, log=print):
    """Runs the plugin for one origin; returns (summary dict, [StationResult])."""
    dump_dir = os.path.join(workdir, "spectra")
    cfg_dir = os.path.join(workdir, "config")
    os.makedirs(dump_dir)
    os.makedirs(cfg_dir)

    # Restrict scamp/scmag to Mw(spec) without touching the system config:
    # a private user-config dir (keeping the user's global.cfg, if any).
    user_global = os.path.join(os.path.expanduser("~/.seiscomp"), "global.cfg")
    if os.path.exists(user_global):
        shutil.copy(user_global, cfg_dir)
    with open(os.path.join(cfg_dir, "scamp.cfg"), "w") as f:
        f.write(f"amplitudes = {TYPE}\n")
    with open(os.path.join(cfg_dir, "scmag.cfg"), "w") as f:
        f.write(f"magnitudes = {TYPE}\n")
    env = dict(os.environ, SEISCOMP_LOCAL_CONFIG=cfg_dir, MWSPEC_DUMP_DIR=dump_dir)

    origin_xml = os.path.join(workdir, "origin.xml")
    if ep:
        shutil.copy(ep, origin_xml)
    else:
        run([SEISCOMP, "exec", "scxmldump", "-d", database, "-O", origin_id,
             "-P", "-o", origin_xml], log=log)

    log("Measuring Mw(spec) amplitudes (fetching waveforms)…")
    amp_xml = os.path.join(workdir, "amplitudes.xml")
    keep_only_origin(origin_xml, origin_id)
    out = run([SEISCOMP, "exec", "scamp", "-d", database, "--ep", origin_xml,
               "-I", recordstream, "--force"], env=env, log=log)
    with open(amp_xml, "w") as f:
        f.write(out)

    log("Computing magnitudes…")
    out = run([SEISCOMP, "exec", "scmag", "-d", database, "--ep", amp_xml],
              env=env, log=log)
    mag_xml = os.path.join(workdir, "magnitudes.xml")
    with open(mag_xml, "w") as f:
        f.write(out)

    return collect(origin_id, origin_xml, mag_xml, dump_dir)


def collect(origin_id, origin_xml, mag_xml, dump_dir):
    stations = {}

    def station(sid):
        if sid not in stations:
            stations[sid] = StationResult(sid)
        return stations[sid]

    for path in sorted(glob.glob(os.path.join(dump_dir, "*.json"))):
        with open(path) as f:
            d = json.load(f)
        d["component"] = d["stream"].split(".")[-1][-1:]
        net, sta = d["stream"].split(".")[:2]
        station(f"{net}.{sta}").dumps.append(d)

    # Distances from the origin's arrivals (via their picks' stations)
    root = ET.parse(origin_xml).getroot()
    ns = _ns(root)
    pick_sta = {}
    for p in root.iter(ns + "pick"):
        w = p.find(ns + "waveformID")
        pick_sta[p.get("publicID")] = f"{w.get('networkCode')}.{w.get('stationCode')}"
    summary = {"originID": origin_id}
    for o in root.iter(ns + "origin"):
        if origin_id and o.get("publicID") != origin_id:
            continue
        summary["time"] = o.findtext(f"{ns}time/{ns}value")
        summary["latitude"] = o.findtext(f"{ns}latitude/{ns}value")
        summary["longitude"] = o.findtext(f"{ns}longitude/{ns}value")
        summary["depth"] = o.findtext(f"{ns}depth/{ns}value")
        for a in o.findall(ns + "arrival"):
            sid = pick_sta.get(a.findtext(ns + "pickID"))
            dist = a.findtext(ns + "distance")
            if sid and dist and sid in stations and stations[sid].distance_km is None:
                stations[sid].distance_km = float(dist) * 111.195
        break

    # Results of scmag (amplitudes, station and network magnitudes)
    root = ET.parse(mag_xml).getroot()
    ns = _ns(root)
    amp_sta = {}
    for a in root.iter(ns + "amplitude"):
        if a.findtext(ns + "type") != TYPE:
            continue
        w = a.find(ns + "waveformID")
        sid = f"{w.get('networkCode')}.{w.get('stationCode')}"
        amp_sta[a.get("publicID")] = sid
        st = station(sid)
        st.amp = _comments(a, ns)
        snr = a.findtext(ns + "snr")
        st.snr = float(snr) if snr else None
    smag_sta = {}
    for o in root.iter(ns + "origin"):
        if origin_id and o.get("publicID") != origin_id:
            continue
        for sm in o.findall(ns + "stationMagnitude"):
            if sm.findtext(ns + "type") != TYPE:
                continue
            sid = amp_sta.get(sm.findtext(ns + "amplitudeID"))
            if not sid:
                w = sm.find(ns + "waveformID")
                sid = f"{w.get('networkCode')}.{w.get('stationCode')}"
            st = station(sid)
            st.mw = float(sm.findtext(f"{ns}magnitude/{ns}value"))
            st.mag = _comments(sm, ns)
            smag_sta[sm.get("publicID")] = sid
        for m in o.findall(ns + "magnitude"):
            if m.findtext(ns + "type") != TYPE:
                continue
            summary["mw"] = float(m.findtext(f"{ns}magnitude/{ns}value"))
            for c in m.findall(ns + "stationMagnitudeContribution"):
                w = c.findtext(ns + "weight")
                sid = smag_sta.get(c.findtext(ns + "stationMagnitudeID"))
                if sid and (w is None or float(w) > 0):
                    stations[sid].used = True
        break

    used = [s for s in stations.values() if s.used] or \
           [s for s in stations.values() if s.mw is not None]
    summary["n"] = len(used)
    for key in ("fc", "stressDrop", "M0"):
        xs = [s.mag[key] for s in used if isinstance(s.mag.get(key), float)]
        summary[key] = statistics.median(xs) if xs else None
    phases = {d["phase"] for s in stations.values() for d in s.dumps}
    summary["phase"] = ",".join(sorted(phases)) or "?"

    result = sorted(stations.values(),
                    key=lambda s: (s.distance_km is None, s.distance_km or 0, s.sid))
    return summary, result


def summary_text(s):
    parts = [f"{TYPE} {s['mw']:.2f} ({s['n']} sta, {s['phase']})" if s.get("mw") is not None
             else f"{TYPE}: no network magnitude ({s['phase']})"]
    if s.get("fc"):
        parts.append(f"fc {s['fc']:.2g} Hz")
    if s.get("stressDrop"):
        parts.append(f"Δσ {s['stressDrop']:.2g} MPa")
    if s.get("M0"):
        parts.append(f"M0 {s['M0']:.2e} N·m")
    return "   ·   ".join(parts)


# ---------------------------------------------------------------------------
#  Plotting (shared by GUI and --png)
# ---------------------------------------------------------------------------

def plot_station(ax, st, legend=True):
    ax.clear()
    ax.set_xscale("log")
    ax.set_yscale("log")
    for d in sorted(st.dumps, key=lambda d: d["component"]):
        comp = d["component"]
        col = COMP_COLORS.get(comp, "#555555")
        g = math.log10(abs(d["gain"])) if d["gain"] else 0.0
        f = d["freq"]
        sig = [10 ** (v - g) for v in d["logSig"]]
        raw = [10 ** (s - c - d["calibration"] - g)
               for s, c in zip(d["logSig"], d["logCorr"])]
        noise = [10 ** (v - g) for v in d["logNoise"] if v > -29]
        ax.plot(f, sig, color=col, lw=1.6, label=f"{comp} corrected")
        ax.plot(f, raw, color=col, lw=0.8, ls=":", label=f"{comp} raw")
        if len(noise) == len(f):
            ax.plot(f, noise, color=col, lw=0.9, alpha=0.35, label=f"{comp} noise (pre-P)")
        fit = d.get("fit")
        if fit:
            om = fit["om0Log10"] - g
            model = [10 ** (om - math.log10(1 + (x / fit["fc"]) ** 2)) for x in f]
            ax.plot(f, model, color="black", lw=1.1, ls="--",
                    label=f"{comp} Brune fc={fit['fc']:.2f} Hz")
            ax.axvspan(fit["fmin"], fit["fmax"], color=col, alpha=0.06)
            ax.axvline(fit["fc"], color=col, lw=0.8, ls="-.")

    dist = f"{st.distance_km:.0f} km" if st.distance_km is not None else "? km"
    title = f"{st.sid}   {dist}   "
    if st.mw is not None:
        title += f"Mw {st.mw:.2f}"
        if isinstance(st.mag.get("stressDrop"), float):
            title += f"   Δσ {st.mag['stressDrop']:.2g} MPa"
        if not st.used:
            title += "   (not used)"
    else:
        title += st.status
    ax.set_title(title, fontsize=9, color="black" if st.mw is not None else "#b00020")
    ax.set_xlabel("Frequency (Hz)", fontsize=8)
    ax.set_ylabel("Displacement spectrum (nm·s)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.grid(True, which="both", alpha=0.2)
    if legend and st.dumps:
        ax.legend(fontsize=6, loc="lower left")


def save_overview(path, summary, stations, ncols=4):
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    n = max(1, len(stations))
    nrows = math.ceil(n / ncols)
    fig = Figure(figsize=(4.2 * ncols, 3.3 * nrows + 0.6))
    FigureCanvasAgg(fig)
    fig.suptitle(summary_text(summary), fontsize=11)
    for i, st in enumerate(stations):
        plot_station(fig.add_subplot(nrows, ncols, i + 1), st, legend=(i == 0))
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=90)


# ---------------------------------------------------------------------------
#  GUI
# ---------------------------------------------------------------------------

def make_gui(app):
    from PyQt5.QtCore import Qt, QThread, pyqtSignal
    from PyQt5.QtGui import QFont
    from PyQt5.QtWidgets import (QFileDialog, QHBoxLayout, QHeaderView, QLabel,
                                 QMainWindow, QPushButton, QSplitter,
                                 QTableWidget, QTableWidgetItem, QVBoxLayout,
                                 QWidget, QAbstractItemView)
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
    from matplotlib.figure import Figure

    class Worker(QThread):
        progress = pyqtSignal(str)
        done = pyqtSignal(object, object)
        failed = pyqtSignal(str)

        def run(self):
            try:
                summary, stations = compute(app.origin_id, app.databaseURI(),
                                            app.rs_url, app.workdir, app.args.ep,
                                            log=self.progress.emit)
                self.done.emit(summary, stations)
            except Exception as e:
                traceback.print_exc()
                self.failed.emit(str(e))

    class NumItem(QTableWidgetItem):
        def __init__(self, value, fmt):
            super().__init__("" if value is None else fmt.format(value))
            self._v = value if value is not None else float("inf")
            self.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)

        def __lt__(self, other):
            return self._v < getattr(other, "_v", float("inf"))

    COLS = [("Station", None), ("Dist km", "{:.0f}"), ("Mw", "{:.2f}"),
            ("fc Hz", "{:.3g}"), ("Δσ MPa", "{:.2g}"), ("M0 N·m", "{:.2e}"),
            ("Band Hz", None), ("Residual", "{:.2f}"), ("SNR", "{:.0f}"),
            ("S onset", None), ("Status", None)]

    class Window(QMainWindow):
        def __init__(self):
            super().__init__()
            self.setWindowTitle(f"Mw(spec) spectra — {app.origin_id}")
            self.resize(1300, 760)
            self.stations = []
            self.summary = None
            central = QWidget()
            self.setCentralWidget(central)
            root = QVBoxLayout(central)
            hdr = QHBoxLayout()
            self.lbl = QLabel("Loading…")
            self.lbl.setFont(QFont("sans-serif", 10, QFont.Bold))
            hdr.addWidget(self.lbl, 1)
            self.btn_png = QPushButton("Save PNG…")
            self.btn_png.setEnabled(False)
            self.btn_png.clicked.connect(self.save_png)
            hdr.addWidget(self.btn_png)
            root.addLayout(hdr)

            split = QSplitter(Qt.Horizontal)
            self.table = QTableWidget(0, len(COLS))
            self.table.setHorizontalHeaderLabels([c for c, _ in COLS])
            self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
            self.table.setSelectionMode(QAbstractItemView.SingleSelection)
            self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
            self.table.verticalHeader().setVisible(False)
            self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            self.table.itemSelectionChanged.connect(self.show_selected)
            split.addWidget(self.table)
            self.fig = Figure(figsize=(6, 5))
            self.canvas = FigureCanvasQTAgg(self.fig)
            self.ax = self.fig.add_subplot(111)
            split.addWidget(self.canvas)
            split.setSizes([650, 650])
            root.addWidget(split, 1)
            self.status = QLabel("")
            self.status.setFont(QFont("monospace", 8))
            root.addWidget(self.status)

        def start(self):
            self.worker = Worker()
            self.worker.progress.connect(self.status.setText)
            self.worker.done.connect(self.loaded)
            self.worker.failed.connect(self.failed)
            self.worker.start()

        def failed(self, msg):
            self.lbl.setText("Mw(spec) computation failed")
            self.status.setText(msg.splitlines()[0] if msg else "")
            self.ax.clear()
            self.ax.text(0.02, 0.98, msg[-1500:], va="top", fontsize=7, family="monospace",
                         transform=self.ax.transAxes)
            self.canvas.draw_idle()

        def loaded(self, summary, stations):
            self.summary, self.stations = summary, stations
            self.lbl.setText(summary_text(summary))
            self.status.setText(f"{len(stations)} stations — spectra by the installed "
                                f"Mw(spec) plugin; click a row (headers sort)")
            self.btn_png.setEnabled(True)
            self.table.setSortingEnabled(False)
            self.table.setRowCount(len(stations))
            for r, st in enumerate(stations):
                fmin = st.amp.get("fmin", st.amp.get("fmin.N"))
                fmax = st.amp.get("fmax", st.amp.get("fmax.N"))
                res = st.amp.get("fitResidual", st.amp.get("fitResidual.N"))
                vals = [st.sid, st.distance_km, st.mw, st.mag.get("fc"),
                        st.mag.get("stressDrop"), st.mag.get("M0"),
                        f"{fmin:.2g}–{fmax:.2g}" if isinstance(fmin, float) and isinstance(fmax, float) else "",
                        res if isinstance(res, float) else None, st.snr,
                        st.onset, st.status]
                for c, ((_, fmt), v) in enumerate(zip(COLS, vals)):
                    item = NumItem(v, fmt) if fmt else QTableWidgetItem(str(v))
                    item.setData(Qt.UserRole, r)
                    if st.mw is None:
                        item.setForeground(Qt.darkRed)
                    elif not st.used:
                        item.setForeground(Qt.darkGray)
                    self.table.setItem(r, c, item)
            self.table.setSortingEnabled(True)
            if stations:
                self.table.selectRow(0)

        def show_selected(self):
            items = self.table.selectedItems()
            if not items:
                return
            st = self.stations[items[0].data(Qt.UserRole)]
            plot_station(self.ax, st)
            self.fig.tight_layout()
            self.canvas.draw_idle()

        def save_png(self):
            path, _ = QFileDialog.getSaveFileName(
                self, "Save all station spectra", f"mwspec_{app.origin_id.replace('/', '_')}.png",
                "PNG (*.png)")
            if path:
                save_overview(path, self.summary, self.stations)
                self.status.setText(f"Saved {path}")

    return Window()


# ---------------------------------------------------------------------------
#  Application
# ---------------------------------------------------------------------------

class ViewerApp(Application):
    def __init__(self, args):
        self.args = args
        self.origin_id = args.origin
        self.rs_url = args.recordstream
        self.workdir = None
        self.window = None
        self._qapp = None
        if not args.png:
            from PyQt5.QtWidgets import QApplication
            self._qapp = QApplication.instance() or QApplication(sys.argv[:1])
        sc_args = [sys.argv[0], "--logging.file=false", "--logging.level=1"]
        if args.database:
            sc_args += ["-d", args.database]
        Application.__init__(self, len(sc_args), sc_args)
        self.setMessagingEnabled(False)
        self.setDatabaseEnabled(True, True)
        self.setDaemonEnabled(False)

    def run(self):
        self.rs_url = self.rs_url or scolv_recordstream()
        if not self.rs_url:
            print("No recordstream: set it in scolv or pass -I", file=sys.stderr)
            return False
        self.workdir = tempfile.mkdtemp(prefix="mwspec_viewer_")
        try:
            if self.args.png:
                summary, stations = compute(self.origin_id, self.databaseURI(),
                                            self.rs_url, self.workdir, self.args.ep,
                                            log=lambda m: print(m, file=sys.stderr))
                save_overview(self.args.png, summary, stations)
                print(summary_text(summary))
                return True
            from PyQt5.QtCore import QTimer
            self.window = make_gui(self)
            self.window.show()
            QTimer.singleShot(100, self.window.start)
            return self._qapp.exec_() == 0
        finally:
            if self.args.keep:
                print(f"work files kept in {self.workdir}", file=sys.stderr)
            else:
                shutil.rmtree(self.workdir, ignore_errors=True)


def main():
    p = argparse.ArgumentParser(description="Mw(spec) spectrum viewer (scolv popup)")
    p.add_argument("origin", help="origin publicID (scolv passes it)")
    p.add_argument("event", nargs="?", help="event publicID (ignored)")
    p.add_argument("-I", "--recordstream", help="RecordStream URL (default: scolv's)")
    p.add_argument("-d", "--database", help="database URI (default: configured)")
    p.add_argument("--ep", help="read origin+picks from this XML instead of the database")
    p.add_argument("--png", help="no GUI: write all station spectra to this PNG")
    p.add_argument("--keep", action="store_true", help="keep the work directory")
    args = p.parse_args()
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
    return 0 if ViewerApp(args)() else 1


if __name__ == "__main__":
    sys.exit(main())
