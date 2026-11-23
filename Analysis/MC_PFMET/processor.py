# to run: python3 -m pepper.runproc stau_processor_trigger_study.py configs/proc_2024/stau2024_trigger_study.json --debug
#   --condor 200 --output prod_v1
import pepper
import logging
import awkward as ak
from functools import partial
import awkward as ak
import numpy as np
import pepper
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stau_processor_base import BaseProcessor  

logger = logging.getLogger(__name__)


class Processor(BaseProcessor):
    """Signal-region processor. Shared functions come from BaseProcessor."""
    def zero_handler(func):
        def _function(self, data, *args, **kwargs):
            if len(data) > 0: return func(self, data, *args, **kwargs)
            else: return ak.Array([])
        return _function

    def __init__(self, config, eventdir):
        # Base __init__ sets histogram_format and calls the pepper init
        super().__init__(config, eventdir)

        if "pileup_reweighting" not in config:
            logger.error("No pileup reweigthing specified")
        if "DY_ZptLO_weights" not in config:
            logger.error("No DY k-factor specified")
        if "MET_trigger_sfs" not in config:
            logger.error("No MET trigger scale factors specified")

        if "DY_jet_reweight" not in config or len(config["DY_jet_reweight"]) == 0:
            logger.error("No DY jet-binned sample stitching is specified")
        if "W_jet_reweight" not in config or len(config["W_jet_reweight"]) == 0:
            logger.error("No WjetsToLNu jet-binned sample stitching is specified")

        if "jet_fake_rate" not in config and len(config["jet_fake_rate"]) == 0:
            self.predict_jet_fakes = False
            logger.warning("No jet fake rate specified")
        else:
            self.predict_jet_fakes = config["predict_yield"]

        if "jet_veto_map" not in config:
            logger.error("No jet veto map is specified")

        if ("jet_uncertainty" not in config and config["compute_systematics"]):
            logger.warning("No jet uncertainty specified")

        if ("jet_resolution" not in config or "jet_ressf" not in config):
            logger.warning("No jet resolution or no jet resolution scale "
                           "factor specified. This is necessary for "
                           "smearing, even if not computing systematics")
        if "jet_correction_mc" not in config and (
                ("jet_resolution" in config and "jet_ressf" in config) or
                ("reapply_jec" in config and config["reapply_jec"])):
            raise pepper.config.ConfigError(
                "Need jet_correction_mc for propagating jet "
                "smearing/variation to MET or because reapply_jec is true")
        if ("jet_correction_data" not in config and "reapply_jec" in config
                and config["reapply_jec"]):
            raise pepper.config.ConfigError(
                "Need jet_correction_data because reapply_jec is true")

        if config["propagate_eff_factors"] and "signal_eff_factors" in config:
            self.propagate_eff_factors = True
        else:
            self.propagate_eff_factors = False

    def process_selection(self, selector, dsname, is_mc, filler):
        era = self.get_era(selector.data, is_mc)

        if self.config["stau_properties_study"]:
            assert dsname.startswith("SMS-TStauStau")
            self.signal_process_study(selector, dsname, is_mc)
            return

        # Triggers, MET filters, PV 
        selector.add_cut("MET", self.MET_cut)
        self.add_trigger(selector, dsname, is_mc)
        self.add_quality_cuts(selector, is_mc)
        selector.add_cut("MET_trigger_pass", self.met_trigger_pass)
        # one tight tag muon
        selector.set_column("Muon_tag", self.muons_tag)
        selector.add_cut("single_muon_tag", self.single_muon_tag_cut)
        selector.add_cut("Lead_pass_trigger", self.lead_pass_trigger)
        # selector.add_cut("muon_sfs", partial(self.get_muon_sfs, is_mc=is_mc)) # this does not work

        # lepton selection
        logger.debug("Running lepton veto")
        selector.set_column("muon_veto", self.muons_veto)
        selector.add_cut("loose_muon_veto", self.loose_muon_veto_cut)
        selector.set_column("electron_veto", self.electron_veto)
        selector.add_cut("loose_electron_veto", self.loose_electron_veto_cut)
        
        # selector.set_column("mt_muon", partial(self.mt, name="Muon_tag"))
        # selector.add_cut("mt_muon", self.mt_muon_cut)

        # if is_mc:
        #     selector.set_column("sum_ll_gen", self.sum_ll_gen)
        selector.add_cut("dy_gen_sfs", partial(self.get_dy_gen_sfs, is_mc=is_mc, dsname=dsname))
        selector.set_column("Jet_select", self.jet_selection)
        if self.config["b-veto"]:
            selector.add_cut("b-veto", partial(self.b_tagged_jet_cut, name="Jet_select"))

        # ---------------------- jet selections - to match baseline region ----------------------
        selector.set_column("Jet_select", self.jet_basis)
        selector.set_multiple_columns(self.missing_energy)

        selector.set_column("Jet_select", self.getloose_jets)
        selector.set_column("PfCands", self.pfcand_valid)
        selector.set_column("Jet_lead_pfcand", partial(self.get_matched_pfCands, match_object="Jet_select", dR=0.4))
        selector.set_column("Jet_select", self.set_jet_dxy)
        # selector.add_cut("mt_muon2", self.mt_muon_cut)
        selector.add_cut("atleast_two_loose_jets", self.has_two_jets)

        selector.set_column("sum_jj", self.sum_jj)
        # selector.set_multiple_columns(self.missing_energy) # sets HT_valid, HT_miss_valid, HT and HT_miss (all defined in base processor)
        selector.set_multiple_columns(self.mt_jets)
        selector.set_column("dphi_jet1_jet2", self.dphi_jet1_jet2)
        selector.set_column("dr_jet1_jet2", self.dr_jet1_jet2)

        selector.add_cut("dphi_min_cut", self.dphi_min_cut)
        selector.set_column("mt2_j1_j2_MET", self.get_mt2)
        selector.set_column("binning_schema", self.binning_schema)

        selector.set_multiple_columns(self.set_njets_pass)
        if self.config["compute_systematics"] and is_mc and self.propagate_eff_factors:
            selector.add_cut("signal_eff_sfs", partial(self.signal_eff_unc, dsname=dsname))
        if self.config["predict_yield"] and not is_mc:
            selector.set_multiple_columns(partial(self.predict_yield, weight=selector.systematics["weight"]))

        selector.add_cut("atleast_two_loose_jets_final", self.has_two_jets)
        selector.set_column("Jet_select", self.gettight_jets)
        selector.add_cut("two_tight_jets", self.has_two_jets)

        # selector.add_cut("ht_cut", self.ht_cut)

    @zero_handler
    def met_trigger_pass(self, data):
        iso_mu24 = data["HLT"]["IsoMu24"]
        met_pass = None
        for path in self.config["met_trigger_paths"]:
            trigger_name = path.removeprefix("HLT_")
            bit = data["HLT"][trigger_name]

            logger.info(
                f"Trigger {trigger_name}: {ak.sum(bit)} events pass"
            )

            met_pass = bit if met_pass is None else (met_pass | bit)

        # Require BOTH IsoMu24 and the MET trigger OR
        passed = iso_mu24 & met_pass

        logger.info(
            f"IsoMu24: {ak.sum(iso_mu24)} events pass"
        )
        logger.info(
            f"Combined IsoMu24 AND MET trigger OR: {ak.sum(passed)} events pass"
        )

        return passed

    @zero_handler
    def muons_tag(self, data):
        muons = data["Muon"]
        is_good = (
            (muons.pt > self.config["muon_pt_min"])
            & (muons.eta < self.config["muon_eta_max"])
            & (muons.eta > self.config["muon_eta_min"])
            & (muons[self.config["muon_ID"]] == 1)
            & (muons.pfIsoId >= self.config["muon_pfIsoId"])
            & (abs(muons.dxy) <= self.config["muon_absdxy"])
            & (abs(muons.dz) <= self.config["muon_absdz"])
        )
        return muons[is_good]

    @zero_handler
    def single_muon_tag_cut(self, data):
        return ak.num(data["Muon_tag"]) == 1

    @zero_handler
    def lead_pass_trigger(self, data):
        lead_muons = data["Muon_tag"][:, :1]
        trig_muons = data["TrigObj"]
        trig_muons = trig_muons[
            (abs(trig_muons.id) == 13)
            & (trig_muons.filterBits >= 8)
        ]
        matches, dRlist = lead_muons.nearest(trig_muons, return_metric=True, threshold=0.2)
        has_matched = ~ak.is_none(matches, axis=1)
        return has_matched[:, 0]  # since lead_muons has form: [ [muon1], [muon2], ....]
    
    @zero_handler
    def get_muon_sfs(self, data, is_mc):
        weight = np.ones(len(data))
        if is_mc:
            muons = data["Muon_tag"]
            id_iso_sfs, systematics = \
                self.muon_sfs(muons, sfs_name_config="muon_sf")
            weight *= ak.to_numpy(id_iso_sfs)
            muon_leading = muons[:, :1]  # trigger_sfs are applied only to leading muon
            mu_trigger_sfs, systematics_trig = \
                self.muon_sfs(muon_leading, sfs_name_config="muon_sf_trigger")
            weight *= ak.to_numpy(mu_trigger_sfs)
            systematics.update(systematics_trig)
            return weight, systematics
        else:
            return weight
    
    @zero_handler
    def muons_veto(self, data):
        muons = data["Muon"]

        is_good = (
            (muons.pt > self.config["muon_veto_pt_min"])
            & (muons.eta < self.config["muon_veto_eta_max"])
            & (muons.eta > self.config["muon_veto_eta_min"])
            & (muons[self.config["muon_veto_ID"]] == 1)
            & (muons.pfIsoId >= self.config["muon_veto_pfIsoId"])
        )

        is_tag = (
            (muons.pt > self.config["muon_pt_min"])
            & (muons.eta < self.config["muon_eta_max"])
            & (muons.eta > self.config["muon_eta_min"])
            & (muons[self.config["muon_ID"]] == 1)
            & (muons.pfIsoId >= self.config["muon_pfIsoId"])
            & (abs(muons.dxy) <= self.config["muon_absdxy"])
            & (abs(muons.dz) <= self.config["muon_absdz"])
        )

        return muons[(is_good & (~is_tag))]

    @zero_handler
    def loose_muon_veto_cut(self, data):
        return ak.num(data["muon_veto"]) == 0
    
    @zero_handler
    def electron_veto(self, data):
        ele = data["Electron"]
        is_good = (
            (ele.pt > self.config["elec_veto_pt"])
            & (ele.eta < self.config["elec_veto_eta_max"])
            & (ele.eta > self.config["elec_veto_eta_min"])
            & (ele[self.config["elec_veto"]] == 1)
            & (ele[self.config["elec_ID"]] == 1)
        )
        return ele[is_good]


    @zero_handler
    def loose_electron_veto_cut(self, data):
        return ak.num(data["electron_veto"]) == 0

    @zero_handler
    def mt_muon_cut(self, data):
        return ((data["mt_muon"] > self.config["mt_muon_min"]) & (data["mt_muon"] < self.config["mt_muon_max"]))

    @zero_handler
    def b_tagged_jet_cut(self, data, name="Jet"):  # base uses jets.btagged
        jets = data[name]
        # Jet_btagDeepFlavB satisfies the Medium (>0.2783) WP:
        # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18
        year = self.config["year"]
        tagger = "deepjet"
        wp = "medium"
        wptuple = pepper.scale_factors.BTAG_WP_CUTS[tagger][year]
        if not hasattr(wptuple, wp):
            raise pepper.config.ConfigError(
                "Invalid working point \"{}\" for {} in year {}".format(
                    wp, tagger, year))
        b_tagged_idx = jets["btagDeepFlavB"] > getattr(wptuple, wp)
        return ak.num(jets[b_tagged_idx], axis=1) == 0

# ---------------------- jet selections - to match baseline region ----------------------

    @zero_handler
    def jet_selection(self, data):
        jets = data["JetCHS"]

        # Apply jet kinematic selection
        jets = jets[
            (self.config["jet_eta_min"] < jets.eta)
            & (jets.eta < self.config["jet_eta_max"])
            & (self.config["jet_pt_min"] < jets.pt)
            # & (jets.jetId >= self.config["jet_jetId"])
        ]

        muons = data["Muon_tag"]

        # Calculate delta R between every jet and every muon
        mu_b = muons[:, None, :]

        deta = jets.eta[:, :, None] - mu_b.eta
        dphi = (jets.phi[:, :, None] - mu_b.phi + np.pi) % (2 * np.pi) - np.pi

        dR = np.sqrt(deta**2 + dphi**2)

        # Keep jets with NO muon within the veto cone
        has_muon_match = ak.any(
            dR < self.config["tag_muon_veto_dR"],
            axis=-1
        )

        isoJets = jets[~has_muon_match]

        return isoJets

    @zero_handler
    def jet_basis(self, data):
        jets = data["Jet_select"]
        n_before = ak.num(jets, axis=1)
        logger.info(
            f"Selecting CHS jets: "
            f"exactly 2 = {ak.sum(n_before == 2)}, "
            f">=2 = {ak.sum(n_before >= 2)}"
        )
        return jets

    @zero_handler
    def set_jet_dxy(self, data):
        jets = data["Jet_select"]
        n_before = ak.num(jets, axis=1)
        logger.info(
            f"Before dxy cut: "
            f"exactly 2 = {ak.sum(n_before == 2)}, "
            f">=2 = {ak.sum(n_before >= 2)}"
        )
        jets = jets[jets.dxy >= self.config["jet_dxy_min"]]
        n_after = ak.num(jets, axis=1)
        logger.info(
            f"After dxy cut: "
            f"exactly 2 = {ak.sum(n_after == 2)}, "
            f">=2 = {ak.sum(n_after >= 2)}"
        )
        return jets

    @zero_handler
    def ht_cut(self, data):
        return data["HT_valid"] > 150

    @zero_handler
    def jet_veto(self, data, jet_name):
        # mask events with at least one jet in veto map (might be too tight)
        jets = data[jet_name]
        config = self.config["jet_veto_map"]
        mask_per_jet = config(eta=jets.eta, phi=jets.phi)
        mask = ak.any(mask_per_jet, axis=-1)
        return ~mask

    @zero_handler
    def has_jets(self, data):
        return ak.num(data["Jet_select"]) > 0

    @zero_handler
    def has_two_jets(self, data):
        jets = data["Jet_select"]
        n = ak.num(jets)
        logger.info(f"n jets: events={len(data)}, >=2 jets: {ak.sum(n >= 2)}")
        return n >= 2
    
    @zero_handler
    def gettight_jets(self, data):
        jets = data["Jet_select"]
        jets = jets[(jets.disTauTag_score1 >= self.config["tight_thr"])]
        n_after = ak.num(jets, axis=1)
        logger.info(
            f"after tight cut: "
            f"exactly 2 = {ak.sum(n_after == 2)}, "
            f">=2 = {ak.sum(n_after >= 2)}"
        )
        return jets


    @zero_handler
    def signal_eff_unc(self, data, dsname):
        weights = np.ones(len(data))
        if not dsname.startswith("SMS-TStauStau"):
            return weights
        else:
            factor = self.config["signal_eff_factors"]
            systematics = {}
            RT1 = data["tight_bin1"]
            RT2 = data["tight_bin2"]
            weights[RT1] = factor["central"]
            weights[RT2] = factor["central"]**2
            if self.config["compute_systematics"]:
                weights_up = np.ones(len(data))
                weights_up[RT1] = factor["up"]
                weights_up[RT2] = factor["up"]**2
                weights_down = np.ones(len(data))
                weights_down[RT1] = factor["down"]
                weights_down[RT2] = factor["down"]**2
                systematics["signal_eff_unc_up"] = weights_up / weights
                systematics["signal_eff_unc_down"] = weights_down / weights
            return weights, systematics
