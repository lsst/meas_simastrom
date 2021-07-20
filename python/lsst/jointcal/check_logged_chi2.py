# This file is part of jointcal.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
Extract chi2 and degrees of freedom values logged by one or more jointcal runs,
print warnings about oddities, and make plots.
"""

import argparse
import dataclasses
import itertools
import os.path
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import seaborn as sns  # noqa: E402
sns.set_style("ticks", {"legend.frameon": True})
sns.set_context("talk")


@dataclasses.dataclass
class Chi2Data:
    """Store the chi2 values read in from a jointcal log file.
    """
    kind: list()
    raw: np.ndarray
    ndof: np.ndarray
    reduced: np.ndarray
    init_count: int = dataclasses.field(init=False)

    def __post_init__(self):
        # ensure the array types are correct
        self.raw = np.array(self.raw, dtype=np.float64)
        self.ndof = np.array(self.ndof, dtype=np.int)
        self.reduced = np.array(self.reduced, dtype=np.float64)
        self.init_count = self._find_init()

    def _find_init(self):
        """Return the index of the first "fit step", after initialization.

        NOTE
        ----
        There are never more than ~25 items in the list, so search optimization
        is not worth the trouble.
        """
        # Logs pre-DM-25779
        if "Fit prepared" in self.kind:
            return self.kind.index("Fit prepared") + 1
        # Logs post-DM-25779
        elif "Fit iteration 0" in self.kind:
            return self.kind.index("Fit iteration 0")
        else:
            raise RuntimeError(f"Cannot find end of initialization sequence in {self.kind}")


class LogParser:
    """Parse a jointcal logfile to extract chi2 values and plot them.

    Call the instance with the path to a file to check it for anamolous chi2
    and output plots to your current directory.

    Parameters
    ----------
    plot : `bool`
        Make plots for each file (saved to the current working directory)?
    verbose : `bool`
        Print extra updates during processing?
    longlog : `bool`
        Use a regex designed for the long
    """
    def __init__(self, plot=True, verbose=True, longlog=True):
        # This regular expression extracts the chi2 values, and the "kind" of
        # chi2 (e.g. "Initial", "Fit iteration").
        # Chi2 values in the log look like this, for example for the "normal" format:
        # jointcal INFO: Initial chi2/ndof : 2.50373e+16/532674=4.7003e+10
        if longlog is False:
            chi2_re = r"jointcal INFO: (?P<kind>.+) chi2/ndof : (?P<chi2>.+)/(?P<ndof>.+)=(?P<reduced_chi2>.+)"  # noqa: E501
        else:
            # and for the "longlog" format:
            # INFO  2021-07-20T11:08:55.715-0500 jointcal ()(jointcal.py:1391)- Initial chi2/ndof : 4.44678e+16/508342=8.74762e+10  # noqa: E501
            chi2_re = r"INFO (.*) jointcal \(\)(\(jointcal\.py\:[0-9]*\))- (?P<kind>.+) chi2/ndof : (?P<chi2>.+)/(?P<ndof>.+)=(?P<reduced_chi2>.+)"  # noqa: E501
        self.matcher = re.compile(chi2_re)
        self.plot = plot
        self.verbose = verbose

        # Reuse the Figure to speed up plotting and save memory.
        self.fig = plt.figure(figsize=(15, 8))

        # How to find the beginning and end of the relevant parts of the log
        # to scan for chi2 values.
        self.section_start = {"astrometry": "=== Starting astrometric fitting...",
                              "photometry": "=== Starting photometric fitting..."}
        # Scan for astrometry until the "====== Now processing photometry..." line,
        # and for photometry until the end of the file.
        self.section_end = {"astrometry": "====== Now processing photometry...",
                            "photometry": None}

    def __call__(self, logfile):
        """Parse logfile to extract chi2 values and generate and save plots.

        The plot output is written to the current directory, with the name
        derived from the basename of ``logfile``.

        Parameters
        ----------
        logfile : `str`
            The filename of the jointcal log to process.

        Returns
        -------
        astrometry, photometry : `Chi2Data`
            The chi2 values that were found in the log.
        """
        title = os.path.basename(logfile)
        if self.verbose:
            print("Processing:", title)

        with open(logfile) as opened_log:
            # Astrometry is always run first, so we can scan for that until the
            # end of that section, and then continue scanning for photometry.
            astrometry = self._extract_chi2(opened_log, "astrometry")
            increased = self._find_chi2_increase(astrometry, title, "astrometry")
            photometry = self._extract_chi2(opened_log, "photometry")
            increased |= self._find_chi2_increase(photometry, title, "photometry")

        if astrometry is None and photometry is None and self.verbose:
            print(f"WARNING: No chi2 values found in {logfile}.")

        if increased or self.plot:
            self._plot(astrometry, photometry, title)
            plotfile = f"{os.path.splitext(title)[0]}.png"
            plt.savefig(plotfile, bbox_inches="tight")
            print("Saved plot:", plotfile)

        return astrometry, photometry

    def _find_chi2_increase(self, chi2Data, title, label, threshold=1):
        """Return True and print a message if the raw chi2 increases
        markedly.
        """
        if chi2Data is None:
            return False
        diff = np.diff(chi2Data.raw)
        ratio = diff/chi2Data.raw[:-1]
        if np.any(ratio > threshold):
            increased = np.where(ratio > threshold)[0]
            print(f"{title} has increasing {label} chi2:")
            for x in zip(chi2Data.raw[increased], chi2Data.raw[increased + 1],
                         ratio[increased], diff[increased]):
                print(f"{x[0]:.6} -> {x[1]:.6} (ratio: {x[2]:.6}, diff: {x[3]:.6})")
            return True
        return False

    def _extract_chi2(self, opened_log, section, verbose=True):
        """Return the values extracted from the chi2 statements in the logfile.
        """
        if verbose:
            print(f"Scanning {section}")
        start = self.section_start[section]
        end = self.section_end[section]
        kind = []
        chi2 = []
        ndof = []
        reduced = []
        # Skip over lines until we get to the section start line.
        for line in opened_log:
            if start in line:
                break

        for line in opened_log:
            # Stop parsing at the section end line, or end of file if end is None.
            if end is not None and end in line:
                break
            if "chi2" in line:
                match = self.matcher.search(line)
                if match is not None:
                    kind.append(match.group("kind"))
                    chi2.append(match.group("chi2"))
                    ndof.append(match.group("ndof"))
                    reduced.append(match.group("reduced_chi2"))
                    print(f"Found: {kind[-1]}: {chi2[-1]}/{ndof[-1]} = {reduced[-1]}")

        # No chi2 values were found (e.g., photometry wasn't run).
        if len(kind) == 0:
            return None

        return Chi2Data(kind, np.array(chi2, dtype=np.float64),
                        np.array(ndof, dtype=int), np.array(reduced, dtype=np.float64))

    def _plot(self, astrometry, photometry, title):
        """Generate plots of chi2 values.

        Parameters
        ----------
        astrometry : `Chi2Data` or None
            The as-read astrometry data, or None if there is none to plot.
        photometry : `Chi2Data` or None
            The as-read photometry data, or None if there is none to plot.
        title : `str`
            Title for the whole plot.
        """
        palette = itertools.cycle(sns.color_palette())

        self.fig.clf()
        ax0, ax1 = self.fig.subplots(ncols=2, gridspec_kw={"wspace": 0.05})

        self.fig.suptitle(title)
        # Use a log scale if any of the chi2 values are very large.
        if max(getattr(astrometry, "raw", [0])) > 100 or max(getattr(photometry, "raw", [0])) > 100:
            ax0.set_yscale("log")
        ax1.yaxis.set_label_position("right")
        ax1.yaxis.tick_right()

        if astrometry is not None:
            patch1, patch2 = self._plot_axes(ax0, ax1, astrometry, palette, label="astrometry")

        if photometry is not None:
            patch3, patch4 = self._plot_axes(ax0, ax1, photometry, palette, label="photometry")

        # Let matplotlib figure out the best legend location: if there is data
        # in the "upper right", we definitely want to see it.
        handles, labels = ax0.get_legend_handles_labels()
        ax1.legend(handles, labels)

    def _plot_axes(self, ax0, ax1, chi2Data, palette, label=""):
        """Make the chi2 and degrees of freedom subplots."""
        xrange = np.arange(0, len(chi2Data.raw), dtype=float)

        # mark chi2=1
        ax0.axhline(1, color='grey', ls='--')
        # mark the separation between initialization and iteration
        ax0.axvline(chi2Data.init_count-0.5, color='grey', lw=0.9)
        color = next(palette)
        patch1 = ax0.plot(xrange[:chi2Data.init_count], chi2Data.raw[:chi2Data.init_count], '*', ms=10,
                          label=f"{label} pre-init", color=color)
        patch2 = ax0.plot(xrange[chi2Data.init_count:], chi2Data.raw[chi2Data.init_count:], 'o', ms=10,
                          label=f"{label} post-init", color=color)
        patch1 = ax0.plot(xrange[:chi2Data.init_count], chi2Data.reduced[:chi2Data.init_count], '*',
                          markerfacecolor="none", ms=10, color=color)
        patch2 = ax0.plot(xrange[chi2Data.init_count:], chi2Data.reduced[chi2Data.init_count:], 'o',
                          markerfacecolor="none", ms=10, label=f"{label} reduced", color=color)

        ax0.set_xlabel("Iteration #", fontsize=20)
        ax0.set_ylabel(r"$\chi ^2$", fontsize=20)

        # mark the separation between initialization and iteration
        ax1.axvline(chi2Data.init_count-0.5, color='grey', lw=0.9)
        ax1.plot(xrange[:chi2Data.init_count], chi2Data.ndof[:chi2Data.init_count], '*', ms=10,
                 label="pre-init", color=color)
        ax1.plot(xrange[chi2Data.init_count:], chi2Data.ndof[chi2Data.init_count:], 'o', ms=10,
                 label="post-init", color=color)

        ax1.set_xlabel("Iteration #", fontsize=20)
        ax1.set_ylabel("# degrees of freedom", fontsize=20)

        return patch1[0], patch2[0]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("files", metavar="files", nargs="+",
                        help="Log file(s) to extract chi2 values from.")
    parser.add_argument("--plot", action="store_true",
                        help="Generate a plot PNG for each log file, otherwise just for questionable ones.")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Print extra information during processing.")
    return parser.parse_args()


def main():
    args = parse_args()
    log_parser = LogParser(plot=args.plot, verbose=args.verbose)
    for file in args.files:
        log_parser(file)
