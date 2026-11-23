import numpy as np
import uproot
import matplotlib.pyplot as plt
from scipy.stats import beta

# --------------------------------------------------
# ROOT files
# --------------------------------------------------

D_DEN = "data_IsoMuonly/hists/"      # only isomu24 trigger
D_NUM = "data_IsoMuonandpfmet/hists/"          # with PFMET120 trigger

DEN = D_DEN + "Cut_012_atleast_two_loose_jets_final_min_MET_MHT.root"
NUM = D_NUM + "Cut_013_atleast_two_loose_jets_final_min_MET_MHT.root"

# groups = {
#     "Signal (all stau samples)": [
#         f"SMS-TStauStau_MStau-500_ctau-{ctau}mm_mLSP-1"
#         # for mass in (400, 500, 600)
#         for ctau in (1, 10, 100)
#     ],

#     "W+jets": [
#         f"WtoLNu-4Jets_Bin-{i}J"
#         for i in ( 3, 4)
#     ],
# }


# --------------------------------------------------
# Datasets (data only)
# --------------------------------------------------

groups = {
    "Muon": ["Muon"],
}

# Plots to produce: (x-axis max in GeV, output filename)
PLOTS = [
    (1000, "data_minMHT_MET_turnon_0-1000.png"),
    (2000, "data_minMHT_MET_turnon_0-2000.png"),
]

# Set to True to print per-bin numbers for debugging
DEBUG = True


def get(path, datasets):
    f = uproot.open(path)

    vals = None
    var = None
    edges = None

    for ds in datasets:
        hist_path = f"{ds}/nominal/hist"
        h = f[hist_path]

        h_vals = h.values(flow=False)
        h_var = h.variances(flow=False)
        h_edges = h.axis().edges()

        # Data: unweighted counts. If no sumw2 is stored, use Poisson (var = N).
        if h_var is None:
            h_var = h_vals.copy()

        if vals is None:
            vals = h_vals.copy()
            var = h_var.copy()
            edges = h_edges.copy()
        else:
            if not np.allclose(edges, h_edges):
                raise ValueError("Inconsistent histogram binning!")

            vals += h_vals
            var += h_var

    f.close()
    return vals, var, edges


# --------------------------------------------------
# Calculate efficiency for each group (done once)
# --------------------------------------------------

results = {}

for label, datasets in groups.items():
    # n: denominator
    # k: numerator
    n, n_var, edges = get(DEN, datasets)
    k, k_var, edges_num = get(NUM, datasets)

    if not np.allclose(edges, edges_num):
        raise ValueError("Numerator and denominator binning differ!")

    eff = np.divide(
        k, n,
        out=np.zeros_like(k, dtype=float),
        where=n > 0
    )
    eff = np.clip(eff, 0, 1)

    # Effective number of events (equals n for unweighted data)
    n_eff = np.divide(
        n**2, n_var,
        out=np.zeros_like(n, dtype=float),
        where=n_var > 0
    )

    # Enforce 0 <= k_eff <= n_eff (avoids floating-point overshoot)
    k_eff = np.clip(eff * n_eff, 0, n_eff)

    lo = np.zeros_like(eff)
    hi = np.ones_like(eff)

    # Bins with essentially no statistics are treated as "empty"
    valid = n_eff >= 1

    kv, nv = k_eff[valid], n_eff[valid]
    lo[valid] = beta.ppf(0.16, kv, nv - kv + 1)
    hi[valid] = beta.ppf(0.84, kv + 1, nv - kv)

    # Edge cases (use tolerance, not ==)
    lo[valid & (k_eff < 1e-9)] = 0
    hi[valid & (k_eff > n_eff - 1e-9)] = 1

    if DEBUG:
        print(f"\n=== {label} ===")
        for lo_e, hi_e, nn, kk, ne in zip(edges[:-1], edges[1:], n, k, n_eff):
            e = kk / nn if nn > 0 else 0
            flag = "  <-- k>n!" if kk > nn * (1 + 1e-9) else ""
            empty = "  [EMPTY]" if ne < 1 else ""
            print(f"{lo_e:6.0f}-{hi_e:6.0f}  n={nn:10.3g}  k={kk:10.3g}  "
                  f"n_eff={ne:8.2f}  eff={e:.3f}{flag}{empty}")

    results[label] = {
        "x": 0.5 * (edges[1:] + edges[:-1]),
        "xerr": 0.5 * np.diff(edges),
        "eff": eff,
        "lo": lo,
        "hi": hi,
        "valid": valid,
    }


# --------------------------------------------------
# Plotting
# --------------------------------------------------

def make_plot(xmax, outname):
    plt.figure(figsize=(8, 6))

    empty_labeled = False

    for i, (label, r) in enumerate(results.items()):
        color = f"C{i}"
        v = r["valid"]
        e = ~v

        # Normal bins with efficiency + Clopper-Pearson-style errors
        plt.errorbar(
            r["x"][v],
            r["eff"][v],
            xerr=r["xerr"][v],
            yerr=[
                np.clip(r["eff"][v] - r["lo"][v], 0, None),
                np.clip(r["hi"][v] - r["eff"][v], 0, None),
            ],
            fmt="o",
            color=color,
            capsize=3,
            label=label,
        )

        # Empty bins: black X at eff = 0, labeled once
        if e.any():
            plt.errorbar(
                r["x"][e],
                np.zeros(e.sum()),
                xerr=r["xerr"][e],
                fmt="x",
                color="black",
                markersize=8,
                markeredgewidth=2,
                linestyle="none",
                label=None if empty_labeled else "Empty bin (no events)",
            )
            empty_labeled = True

    plt.axhline(1, linestyle="--", color="gray")
    # plt.axhline(0.85, linestyle="--", color="red", label="85% efficiency")
    plt.axhline(0.90, linestyle="--", color="blue", label="90% efficiency")
    plt.axhline(0.95, linestyle="--", color="green", label="95% efficiency")

    plt.xlabel("Minimum of MET and MHT [GeV]")
    plt.ylabel("Trigger efficiency")
    plt.title("PFMET120_PFMHT120_IDTight Trigger Turn-on")

    plt.xlim(0, xmax)
    plt.ylim(0, 1.1)
    plt.grid(alpha=0.3)
    plt.legend(loc="lower right")

    plt.tight_layout()
    plt.savefig(outname, dpi=300)


for xmax, outname in PLOTS:
    make_plot(xmax, outname)

plt.show()