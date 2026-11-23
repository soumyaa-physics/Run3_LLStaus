# Analysis (STau 2018 UL setup) - Manual

## Steps for RUN-3 analysis:

Using the stau2024_signal.json config file i run the processor: stau_processor_signal.py using:
```sh
python3 -m pepper.runproc --debug stau_processor_signal.py configs/proc_2024/stau2024_signal.json
```
you can get cutflows and histograms

## Trigger Efficiency Study

The MET trigger efficiency is measured relative to the `IsoMu24` reference trigger:

$$
\epsilon_{\mathrm{MET}} = \frac{N(\mathrm{IsoMu24} \cap \mathrm{MET})}{N(\mathrm{IsoMu24})}
$$

The trigger is defined in `configs/proc_2024/stau2024_trigger_study.json`. The numerator and denominator differ only by the MET trigger cut in `stau_processor_trigger_study.py` which encodes the trigger whose efficiency you want to calculate:

```python
selector.add_cut("MET_trigger_pass", self.met_trigger_pass)
```

| Sample | MET trigger cut |
|---|---|
| Denominator (`IsoMu24` only) | **Commented out** |
| Numerator (`IsoMu24` and `MET`) | **Uncommented** |

### Running the Processor

Run the processor once for each configuration above.

Debug mode:

```bash
python3 -m pepper.runproc stau_processor_trigger_study.py configs/proc_2024/stau2024_trigger_study.json --debug
```

Condor mode:

```bash
python3 -m pepper.runproc stau_processor_trigger_study.py configs/proc_2024/stau2024_trigger_study.json --condor 200 --output withMETcut_MHT
```

### Trigger Turn-On Plots

Use `plot_turnon.py` to produce the turn-on plots. Choose one final offline selection:

- `atleast_two_loose_jets_final`
- `two_tight_jets`

**Use the same offline cut for the numerator and denominator** so the comparison is consistent.