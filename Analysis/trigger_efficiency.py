import logging
import os
import warnings

import awkward as ak
import coffea.hist
import coffea.processor
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import numpy

from coffea.nanoevents import NanoAODSchema

warnings.filterwarnings("ignore", category=FutureWarning)
hep.style.use("CMS")


def build_objs(events, opt):
    """Return a record array (pt, eta, charge) per event for the given object type."""
    if opt == "mu":
        gp = events.GenPart
        mu = gp[
            (abs(gp.pdgId) == 13)
            & gp.hasFlags(["isLastCopy"])
            & gp.hasFlags(["isTauDecayProduct"])
        ]
        return ak.zip({
            "pt": mu.pt,
            "eta": mu.eta,
            "charge": -numpy.sign(mu.pdgId),
        })

    elif opt == "tauh":
        tau = events.GenVisTau
        return ak.zip({
            "pt": tau.pt,
            "eta": tau.eta,
            "charge": tau.charge,
        })

    raise ValueError("%s not a valid option." % opt)


def find_best_pair(obj1, obj2, same_objs):
    """Per event: opposite-charge pair passing pt/eta cuts with the largest pt sum.
    Returns (event_mask, best_pairs) where best_pairs has one record (a, b) per selected event."""
    if same_objs:
        pairs = ak.combinations(obj1, 2, fields=["a", "b"])
    else:
        pairs = ak.cartesian({"a": obj1, "b": obj2})

    good = (
        (pairs.a.pt > 20) & (abs(pairs.a.eta) < 2.5)
        & (pairs.b.pt > 20) & (abs(pairs.b.eta) < 2.5)
        & (pairs.a.charge * pairs.b.charge < 0)
    )
    pairs = pairs[good]

    sel = ak.num(pairs, axis=1) > 0
    pairs = pairs[sel]

    best_idx = ak.argmax(pairs.a.pt + pairs.b.pt, axis=1, keepdims=True)
    best = ak.flatten(pairs[best_idx], axis=1)

    return ak.to_numpy(sel), best


class MyProcessor(coffea.processor.ProcessorABC):

    def __init__(self, genObj1opt, genObj2opt, triggers=None):
        self.genObj1opt = genObj1opt
        self.genObj2opt = genObj2opt
        self.triggers = list(triggers) if triggers else []

        dataset_axis = coffea.hist.Cat("dataset", "Dataset")

        def h(name, nbins, lo, hi):
            return coffea.hist.Hist("Counts", dataset_axis, coffea.hist.Bin(name, name, nbins, lo, hi))

        self._accumulator = coffea.processor.dict_accumulator({
            "MET_pt": h("MET_pt", 100, 0, 1000),
            "GenMET_pt": h("GenMET_pt", 100, 0, 1000),
        })

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.accumulator.identity()
        dataset = events.metadata["dataset"]

        # Trigger selection
        if self.triggers:
            mask = numpy.zeros(len(events), dtype=bool)
            for trig in self.triggers:
                if trig in events.HLT.fields:
                    mask |= ak.to_numpy(events.HLT[trig]).astype(bool)
                else:
                    logging.warning("Trigger %s not found in file, skipping.", trig)
            events = events[mask]

        if len(events) == 0:
            return output

        obj1 = build_objs(events, self.genObj1opt)
        obj2 = build_objs(events, self.genObj2opt)

        sel, _ = find_best_pair(obj1, obj2, self.genObj1opt == self.genObj2opt)

        if not numpy.any(sel):
            return output

        events = events[sel]

        output["MET_pt"].fill(dataset=dataset, MET_pt=ak.to_numpy(events.MET.pt))
        output["GenMET_pt"].fill(dataset=dataset, GenMET_pt=ak.to_numpy(events.GenMET.pt))

        return output

    def postprocess(self, accumulator):
        return accumulator


def get_hist(output, histName, dataset):
    h = output[histName].integrate("dataset", dataset)
    vals = h.values()[()]
    edges = h.axes()[0].edges()
    return vals, edges


def plot_efficiency(den, num, edges, outfile, xrange, xtitle, legendtitle, lumiText):
    frac = num.sum() / den.sum() if den.sum() > 0 else 0.0

    fig, (ax, rax) = plt.subplots(
        2, 1, figsize=(10, 10), sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
    )

    hep.histplot(den, bins=edges, ax=ax, histtype="step", color="tab:blue", linewidth=2, label="Before HLT")
    hep.histplot(num, bins=edges, ax=ax, histtype="step", color="tab:red", linewidth=2,
                 label="After HLT (frac = %0.2f)" % frac)

    ax.set_ylabel("Counts")
    ax.set_xlim(*xrange)
    ax.set_ylim(0, max(den.max(), 1) * 1.5)
    ax.legend(title=legendtitle, loc="upper right")
    hep.cms.text("Simulation Preliminary", ax=ax)
    ax.set_title(lumiText, loc="right", fontsize=20)

    centers = 0.5 * (edges[1:] + edges[:-1])
    with numpy.errstate(divide="ignore", invalid="ignore"):
        eff = numpy.where(den > 0, num / den, numpy.nan)
        err = numpy.sqrt(eff * (1 - eff) / numpy.where(den > 0, den, numpy.nan))
    rax.errorbar(centers, eff, yerr=err, xerr=0.5 * numpy.diff(edges), fmt="o", color="black", markersize=3)
    rax.set_ylim(0, 1.1)
    rax.set_ylabel("Efficiency")
    rax.set_xlabel(xtitle)

    fig.savefig(outfile)
    fig.savefig(outfile.replace(".pdf", ".png"), dpi=150)
    plt.close(fig)


def main():

    datasets = {
        "stau500_lsp1_ctau10mm": [
            "/data/dust/user/vashisht/CMSSW_17_0_0_pre3/src/Run3_LLStaus/Analysis/stau500_ctau10mm/0000/ExoNanoAOD_1.root"
        ],
    }

    triggers = [
        "HLT_PFMET120_PFMHT120_IDTight",
        "HLT_PFMET130_PFMHT130_IDTight",
        "HLT_PFMET140_PFMHT140_IDTight",
        "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
        "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight",
        "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight",
        "HLT_IsoMu24",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130",
        "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140",
        "HLT_MET105_IsoTrk50",
        "HLT_MET120_IsoTrk50",
        "HLT_DoubleMediumChargedIsoDisplacedPFTauHPS36_Trk1_eta2p1_v1",
    ]

    genObj1opt, genObj2opt = "mu", "tauh"
    extraTag = "incl-HLT-MET-IsoTrk-PFTau"

    treename = "Events"
    executor = coffea.processor.iterative_executor
    executor_args = {"schema": NanoAODSchema}

    output_den = coffea.processor.run_uproot_job(
        datasets, treename,
        MyProcessor(genObj1opt, genObj2opt),
        executor, executor_args,
    )

    output_num = coffea.processor.run_uproot_job(
        datasets, treename,
        MyProcessor(genObj1opt, genObj2opt, triggers=triggers),
        executor, executor_args,
    )

    d_hist = {
        "MET_pt":    {"xrange": (0, 500), "xtitle": r"$p_{T}^{miss}$ [GeV]"},
        "GenMET_pt": {"xrange": (0, 500), "xtitle": r"$p_{T,gen}^{miss}$ [GeV]"},
    }

    tag = "%s-%s%s" % (genObj1opt, genObj2opt, "_%s" % extraTag if extraTag else "")

    for dataset in datasets:

        outdir = "output/plots/%s/%s" % (dataset, tag)
        os.makedirs(outdir, exist_ok=True)

        for histName, cfg in d_hist.items():

            print("*" * 10, dataset, histName, "*" * 10)

            den, edges = get_hist(output_den, histName, dataset)
            num, _ = get_hist(output_num, histName, dataset)

            print("den integral: %f, num integral: %f" % (den.sum(), num.sum()))

            plot_efficiency(
                den, num, edges,
                outfile="%s/%s.pdf" % (outdir, histName),
                xrange=cfg["xrange"],
                xtitle=cfg["xtitle"],
                legendtitle="%s\n(%s)" % (dataset, tag),
                lumiText="2018 (13 TeV)",
            )

    return 0


if __name__ == "__main__":
    main()