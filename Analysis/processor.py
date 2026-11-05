# This file illustrates how to implement a processor, realizing the selection
# steps and outputting histograms and a cutflow with efficiencies.
# Here we create a very simplified version of the ttbar-to-dilep processor.
# One can run this processor using
# 'python3 -m pepper.runproc --debug example_processor.py example_config.json'
# Above command probably will need a little bit of time before all cuts are
# applied once. This is because a chunk of events are processed simultaneously.
# You change adjust the number of events in a chunk and thereby the memory
from nis import match
from unittest import result
import pepper
import awkward as ak
from functools import partial
import numpy as np
import mt2
import numba as nb
import coffea
import logging
import os
from copy import copy

from coffea.nanoevents import NanoAODSchema

logger = logging.getLogger(__name__)

class Processor(pepper.ProcessorBasicPhysics):
    # We use the ConfigTTbarLL instead of its base Config, to use some of its
    # predefined extras
    config_class = pepper.ConfigSTau
    
    def zero_handler(func):
        def _function(self, data, *args, **kwargs):
            if len(data) > 0: return func(self, data, *args, **kwargs)
            else: return ak.Array([])
        return _function
    
    def __init__(self, config, eventdir):
        # Initialize the class, maybe overwrite some config variables and
        # load additional files if needed
        # Can set and modify configuration here as well
        config["histogram_format"] = "root"
        # Need to call parent init to make histograms and such ready
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

        # Triggers
        pos_triggers, neg_triggers = pepper.misc.get_trigger_paths_for(
            dsname,
            is_mc,
            self.config["dataset_trigger_map"],
            self.config["dataset_trigger_order"])
        selector.add_cut("Trigger", partial(
            self.passing_trigger, pos_triggers, neg_triggers))
            
        print("Positive triggers:", pos_triggers)
        print("Negative triggers:", neg_triggers)
        # if is_mc and ( dsname.startswith("DYJetsToLL_M-50") or \
        #             dsname.startswith("DY1JetsToLL_M-50") or \
        #             dsname.startswith("DY2JetsToLL_M-50") or \
        #             dsname.startswith("DY3JetsToLL_M-50") or \
        #             dsname.startswith("DY4JetsToLL_M-50") ):
        #     selector.add_cut("DY jet reweighting",
        #         partial(self.do_dy_jet_reweighting))
        # if is_mc and ( dsname.startswith("WJetsToLNu") or \
        #             dsname.startswith("W1JetsToLNu") or \
        #             dsname.startswith("W2JetsToLNu") or \
        #             dsname.startswith("W3JetsToLNu") or \
        #             dsname.startswith("W4JetsToLNu") ):
        #     selector.systematics["weight"] = \
        #         ak.full_like(selector.systematics["weight"], 1.0)
        #     selector.add_cut("W jet reweighting",
        #         partial(self.do_w_jet_reweighting))

        # if is_mc and "pileup_reweighting" in self.config:
        #     selector.add_cut("Pileup reweighting", partial(
        #         self.do_pileup_reweighting, dsname))
        
        # if is_mc and self.config["year"] in ("2016", "2017", "ul2016pre",
        #                                 "ul2016post", "ul2017"):
        #     selector.add_cut("L1Prefiring", self.add_l1_prefiring_weights)

        selector.add_cut("MET filters", partial(self.met_filters, is_mc))
        
        # PV cut
        selector.add_cut("PV", lambda data: data["PV"].npvsGood > 0)
        
        # lepton selection
        logger.debug(f"Running lepton veto")
        selector.set_column("Muon", self.select_muons)
        selector.set_column("Electron", self.select_electrons)
        selector.add_cut("muon_veto", self.muon_veto)
        selector.add_cut("elec_veto", self.elec_veto)

        selector.set_column("Lepton", lambda data: ak.Array([[]] * len(data)))
        selector.set_column("Jet_select", self.jet_basis)

        if self.config["dxy_cut_study"]:
            self.dxy_cut_study(selector, dsname, is_mc)
            return

        selector.set_column("Jet_select", self.getloose_jets)
        selector.set_column("PfCands", self.pfcand_valid)
        selector.set_column("Jet_lead_pfcand", partial(self.get_matched_pfCands, match_object="Jet_select", dR=0.4))
        selector.set_column("Jet_select", self.set_jet_dxy)
        selector.add_cut("two_loose_jets", self.has_two_jets)

        selector.set_column("sum_jj", self.sum_jj)
        selector.set_multiple_columns(self.missing_energy)
        selector.set_multiple_columns(self.mt_jets)
        selector.set_column("dphi_jet1_jet2", self.dphi_jet1_jet2)
        selector.set_column("dr_jet1_jet2", self.dr_jet1_jet2)

        selector.add_cut("dphi_min_cut", self.dphi_min_cut)
        selector.set_column("mt2_j1_j2_MET", self.get_mt2)
        selector.set_column("binning_schema", self.binning_schema)
        
        # Tagger part for calculating scale factors
        # Scale factors should be calculated -
        # before cuts on the number of the jets
        selector.set_multiple_columns(self.set_njets_pass)
        if self.config["compute_systematics"] and is_mc and self.propagate_eff_factors:
            selector.add_cut("signal_eff_sfs", partial(self.signal_eff_unc, dsname=dsname))
        # selector.set_multiple_columns(self.set_njets_pass_finebin)
        if self.config["predict_yield"] and not is_mc:
            selector.set_multiple_columns(partial(self.predict_yield, weight=selector.systematics["weight"]))
        
        ## Dummy basemark to have cut where regions were not categories yet
        # selector.set_cat("control_region", {"RT0", "RT1", "RT2", "INC"})
        # selector.set_multiple_columns(partial(self.categories_bins))
        selector.add_cut("two_loose_jets_final", lambda data: ak.Array(np.ones(len(data))))
        selector.add_cut("ht_cut", self.ht_cut)


    def dxy_cut_study(self, selector, dsname, is_mc):
        
        selector.set_column("Jet_select", self.gettight_jets)
        selector.set_column("Jet_tau", partial(self.jet_tau, is_mc=is_mc))
        
        selector.add_cut("cut_separator1", lambda data: np.ones(len(data)))
        
        selector.add_cut("one_tight", lambda data: ak.num(data["Jet_select"]) >= 1)
        # Jet_select was not redefined yet
        selector.add_cut("cut_separato2", lambda data: np.ones(len(data)))

        selector.set_column("PfCands", self.pfcand_valid)
        selector.set_column("Jet_lead_pfcand", partial(self.get_matched_pfCands, match_object="Jet_select", dR=0.4))
        selector.set_column("Jet_select", self.set_jet_dxy)
        selector.set_column("Jet_tau", partial(self.jet_tau, is_mc=is_mc))

        selector.add_cut("after_define_dxy", lambda data: np.ones(len(data)))
        selector.add_cut("two_tight", lambda data: ak.num(data["Jet_select"]) >= 2)

    ### Signal process study ---->

    def signal_process_study(self, selector, dsname, is_mc):
        # Define lifetime of GenVisTau:
        selector.set_column("GenVisTau", self.gen_vis_tau)
        selector.set_column("Jet_select", self.jet_selection)
        # match GenVisTau to Jet_select
        selector.set_column("Jet_select", partial(self.jet_tauvis_match, is_mc=is_mc))
        selector.set_column("PfCands", self.pfcand_valid)
        selector.set_column("Jet_lead_pfcand", partial(self.get_matched_pfCands, match_object="Jet_select", dR=0.4))
        selector.set_column("Jet_select", self.set_jet_dxy)
        selector.add_cut("checkpoint", lambda data: np.ones(len(data)))
        selector.set_column("jet_score", self.jet_lead_score)

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

    def gen_vis_tau(self, data):
        tau = data.GenVisTau
        # define trevel distance and transverse travel distance of GenVisTau
        genpart = tau.parent
        assert ak.all(abs(genpart.pdgId) == 15)
        # calculate travel distance and transverse travel distance wrt. to the mother vertex
        tau["travel"] = np.sqrt(genpart.vx**2 + genpart.vy**2 + genpart.vz**2)
        tau["tr_travel"] = np.sqrt(genpart.vx**2 + genpart.vy**2)
        return tau

    @zero_handler
    def jet_tauvis_match(self, data, is_mc):
        jets = data["Jet_select"]
        if not is_mc: return jets
        tau_vis = data.GenVisTau[ ((data.GenVisTau.pt > 30) & (abs(data.GenVisTau.eta) < 2.4) &
                                  (data.GenVisTau.parent.hasFlags(["fromHardProcess"])))
                                ]
        matches_h, dRlist = jets.nearest(tau_vis, return_metric=True, threshold=0.3)
        # add travel distance
        jets["travel"] = matches_h.travel
        jets["tr_travel"] = matches_h.tr_travel
        tau_jet = jets[~ak.is_none(matches_h, axis=1)]

        return tau_jet
    
    ### Signal process study end <---

    def categories_bins(self, data):
        if len(data) == 0:
            return {
                "RT0": ak.Array([]),
                "RT1":ak.Array([]),
                "RT2":ak.Array([]),
                "INC":ak.Array([])}
        jets = data["Jet_select"]
        jets = jets[(jets.disTauTag_score1 >= self.config["tight_thr"])]
        n_tight = ak.num(jets)
        RT0 = (n_tight == 0)
        RT1 = (n_tight == 1)
        RT2 = (n_tight == 2)
        return {"RT0":RT0, "RT1":RT1, "RT2":RT2, "INC":np.ones(len(data))}

    @zero_handler
    def sum_ll_gen(self, data):
        part = data["GenPart"]
        part = part[ part.hasFlags("isLastCopy")
                & (part.hasFlags("fromHardProcess")
                & ((abs(part["pdgId"]) == 11) 
                | (abs(part["pdgId"]) == 13)
                | (abs(part["pdgId"]) == 12)
                | (abs(part["pdgId"]) == 14)
                | (abs(part["pdgId"]) == 15)
                | (abs(part["pdgId"]) == 16)))
                | (part.hasFlags("isDirectHardProcessTauDecayProduct"))
        ]
        sum_p4 = part.sum(axis=1) # sum over all particles in event
        return sum_p4

    @zero_handler
    def get_dy_gen_sfs(self, data, is_mc, dsname):
        weight = np.ones(len(data))
        if is_mc and dsname.startswith("DY"):
            z_boson = data["sum_ll_gen"]
            dy_gen_sfs = self.config["DY_lo_sfs"](mass=z_boson.mass, pt=z_boson.pt)
            weight *= ak.to_numpy(dy_gen_sfs)
            return weight
        else:
            return weight

    @zero_handler
    def do_dy_jet_reweighting(self, data):
        njet = data["LHE"]["Njets"]
        weights = self.config["DY_jet_reweight"][njet]
        return weights
    
    @zero_handler
    def do_w_jet_reweighting(self, data):
        njet = data["LHE"]["Njets"]
        weights = self.config["W_jet_reweight"][njet]
        return weights

    @zero_handler
    def skim_jets(self, data):
        jets = data["Jet_select"]
        file_name = os.path.basename(data.metadata["filename"])
        dataset_name = data.metadata["dataset"]
        prefix = f"evn_{data.metadata['entrystart']}_{data.metadata['entrystop']}_"
        path = self.config["skim_path"]
        path = f"{path}/{dataset_name}"
        print(path)
        print(f"{path}/"+prefix+file_name+".parquet")
        if not os.path.exists(path):
            os.makedirs(path)
        ak.to_parquet(jets, f"{path}/"+prefix+file_name+".parquet")
        return np.ones(len(data))
    
    @zero_handler
    def MET_trigger_sfs(self, data, met_name="MET"):
        met_pt = data[met_name].pt
        scale_factors = self.config["MET_trigger_sfs"]
        met_thr = self.config["MET_trigger_plateau"]
        sfs = scale_factors(pt=met_pt)
        systematics = {}
        if self.config["compute_systematics"]:
            sfs_up = scale_factors(variation="up", pt=met_pt)
            sfs_down = scale_factors(variation="down", pt=met_pt)
            # ratio of the up/down to nominal and for 120-250 GeV - 2x stat. unc. of the MET_trigger_sfs is used, for > 250 - 1x stat. unc.
            systematics["MET_trigger_sfs_up"] = ak.where(met_pt > met_thr, sfs_up / sfs, 2 * sfs_up / sfs)
            systematics["MET_trigger_sfs_down"] = ak.where(met_pt > met_thr, sfs_down / sfs, 2 * sfs_down / sfs)
        return sfs, systematics
    
    @zero_handler
    def HEM_veto(self, data, is_mc):
        weight = np.ones(len(data), dtype=np.float32)
        jets = data["Jet"]
        elctron = data["Electron"]
        electron_in15or16_hem = ( (elctron.pt > 20) & (elctron.eta > -3.0) & (elctron.eta < -1.3) & (elctron.phi > -1.57) & (elctron.phi < -0.87) )
        jet_in15or16_hem = ( (jets.pt > 20) & (jets.eta > -3.2) & (jets.eta < -1.3) & (jets.phi > -1.77) & (jets.phi < -0.67) )
        in_hem = (ak.any(electron_in15or16_hem, axis=-1) | ak.any(jet_in15or16_hem, axis=-1))
        if is_mc:
            weight[in_hem] = (1-0.66)
        else:
            issue_period = (data.run >= 319077)
            weight[in_hem & issue_period] = 0.0
        return weight   
    
    @zero_handler
    def jet_lead(self, data):
        # select only first two pt leading jets
        return data["Jet_select"][:,:2]

    @zero_handler
    def jet_passed(self, data):
        jets = data["Jet_select"]
        jets = jets[jets.disTauTag_score1 >= self.config["tight_thr"]]
        return jets
    
    @zero_handler
    def jet_tau(self, data, is_mc):
        jets = data["Jet_select"]
        if not is_mc: return jets
        tau_vis = data.GenVisTau[ ((data.GenVisTau.pt > 30) & (abs(data.GenVisTau.eta) < 2.4) &
                                  (data.GenVisTau.parent.hasFlags(["fromHardProcess"])))
                                ]
        matches_h, dRlist = jets.nearest(tau_vis, return_metric=True, threshold=0.3)
        tau_jet = jets[~ak.is_none(matches_h, axis=1)]
        return tau_jet
    
    @zero_handler
    def jet_taupass(self, data, is_mc):
        jets = data["Jet_pass"]
        if not is_mc: return jets
        tau_vis = data.GenVisTau[ ((data.GenVisTau.pt > 30) & (abs(data.GenVisTau.eta) < 2.4) &
                                  (data.GenVisTau.parent.hasFlags(["fromHardProcess"])))
                                ]
        matches_h, dRlist = jets.nearest(tau_vis, return_metric=True, threshold=0.3)
        tau_jet = jets[~ak.is_none(matches_h, axis=1)]
        return tau_jet
    
    @zero_handler
    def jet_taupass2(self, data, is_mc):
        jets = data["Jet_pass"]
        if not is_mc: return jets
        tau_vis = data.GenVisTau[ ((data.GenVisTau.pt > 30) & (abs(data.GenVisTau.eta) < 2.4) &
                                  (data.GenVisTau.parent.hasFlags(["fromHardProcess"])))
                                ]
        matches_h, dRlist = jets.nearest(tau_vis, return_metric=True, threshold=0.3)
        two_jets_pass = (ak.num(jets) == 2)
        tau_jet = jets[ (~ak.is_none(matches_h, axis=1)) & two_jets_pass]
        return tau_jet
    
    @zero_handler
    def jet_lead_score(self, data):
        # select only first two biggest score jets
        jets = data["Jet_select"]
        sort_idx = ak.argsort(jets.disTauTag_score1, axis=-1, ascending=False)
        jets = jets[sort_idx]
        return jets[:,:2]
    
    @zero_handler    
    def get_mt2(self, data):
        jet1 = data["Jet_select"][:,0]
        jet2 = data["Jet_select"][:,1]
        met = data["PFMET"]
        jet1_px = jet1.pt * np.cos(jet1.phi)
        jet1_py = jet1.pt * np.sin(jet1.phi)

        jet2_px = jet2.pt * np.cos(jet2.phi)
        jet2_py = jet2.pt * np.sin(jet2.phi)

        met_px = met.pt * np.cos(met.phi)
        met_py = met.pt * np.sin(met.phi)
        return mt2.mt2(
            jet1.mass, jet1_px, jet1_py,
            jet2.mass, jet2_px, jet2_py,
            met_px, met_py,
            0, 0
        )
        
    @zero_handler
    def missing_energy(self, data):

        # Selected jets
        jets = data["Jet_select"]

        HT_valid = ak.sum(jets.pt, axis=-1)

        px = ak.sum(jets.pt * np.cos(jets.phi), axis=-1)
        py = ak.sum(jets.pt * np.sin(jets.phi), axis=-1)

        HT_miss_valid = np.sqrt(px**2 + py**2)

        # All jets
        jets = data["Jet"]

        HT = ak.sum(jets.pt, axis=-1)

        px = ak.sum(jets.pt * np.cos(jets.phi), axis=-1)
        py = ak.sum(jets.pt * np.sin(jets.phi), axis=-1)

        HT_miss = np.sqrt(px**2 + py**2)

        return {
            "HT_valid": HT_valid,
            "HT_miss_valid": HT_miss_valid,
            "HT": HT,
            "HT_miss": HT_miss
        }
        
    def delta_phi(self, phi1_ak, phi2_ak):
        phi1 = np.array(phi1_ak)
        phi2 = np.array(phi2_ak)
        assert phi1.shape == phi2.shape
        d = phi1 - phi2
        indx_pos = d>np.pi
        d[indx_pos] -= np.pi*2
        indx_neg = d<=-np.pi
        d[indx_neg] += np.pi*2
        return d
    
    @zero_handler
    def mt_jets(self, data):

        jet1 = data["Jet_select"][:, 0]
        jet2 = data["Jet_select"][:, 1]

        met = data["PFMET"]

        # Transverse mass of leading jet and MET
        delta_phi1 = self.delta_phi(jet1.phi, met.phi)

        mt_j1 = np.sqrt(
            2 * jet1.pt * met.pt * (1.0 - np.cos(delta_phi1))
        )

        # Transverse mass of subleading jet and MET
        delta_phi2 = self.delta_phi(jet2.phi, met.phi)

        mt_j2 = np.sqrt(
            2 * jet2.pt * met.pt * (1.0 - np.cos(delta_phi2))
        )

        return {
            "mt_jet1": mt_j1,
            "mt_jet2": mt_j2,
            "mt_sum": mt_j1 + mt_j2
        }
        
    @zero_handler
    def dphi_jet1_jet2(self, data):
        return self.delta_phi(data["Jet_select"][:,0].phi,
                              data["Jet_select"][:,1].phi)
    @zero_handler
    def dr_jet1_jet2(self, data):

        jets = data["Jet_select"]

        jet1 = jets[:, 0]
        jet2 = jets[:, 1]

        delta_eta = jet1.eta - jet2.eta
        delta_phi = self.delta_phi(jet1.phi, jet2.phi)

        return np.sqrt(delta_eta**2 + delta_phi**2)
    
    @zero_handler
    def MET_cut(self, data):
        return data["MET"].pt > self.config["MET_cut"]
    
    @zero_handler
    def select_muons(self, data):
        muons = data["Muon"]
        is_good = (
              (muons.pt > self.config["muon_veto_pt_min"])
            & (muons.eta < self.config["muon_veto_eta_max"])
            & (muons.eta > self.config["muon_veto_eta_min"])
            & (muons[self.config["muon_veto_ID"]] == 1)
            & (muons.pfIsoId >= self.config["muon_veto_pfIsoId"])
            )
        return muons[is_good]
    
    @zero_handler
    def select_electrons(self, data):
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
    def muon_veto(self, data):
        return ak.num(data["Muon"]) == 0
    
    @zero_handler
    def elec_veto(self, data):
        return ak.num(data["Electron"]) == 0

    @zero_handler
    def jet_basis(self,data):
        jets = data["JetCHS"]
        return jets
    
    @zero_handler
    def jet_selection(self, data, from_collection="JetCHS"):
        jets = data[from_collection]
        jets = jets[(
            (self.config["jet_eta_min"] < jets.eta)
            & (jets.eta < self.config["jet_eta_max"])
            & (self.config["jet_pt_min"] < jets.pt)
            # & (jets.jetId >= self.config["jet_jetId"] )
            )]
        return jets
    
    @zero_handler
    def has_jets(self, data):
        return ak.num(data["Jet_select"]) > 0

    @zero_handler
    def has_two_jets(self, data):
        jets = data["Jet_select"]
        n = ak.num(jets)
        logger.info(f"n jetsv: events={len(data)}, >=2 jets: {ak.sum(n>=2)}")
        return n >= 2
    
    @zero_handler
    def has_small_dxy(self, data):
        jets = data["Jet_select"]
        dxy1 = (jets.dxy[:,0] < 0.4)
        dxy2 = (jets.dxy[:,1] < 0.4)
        return (dxy1 | dxy2)

    @zero_handler
    def getloose_jets(self, data):
        jets = data["Jet_select"]
        jets = jets[(jets.disTauTag_score1 >= self.config["loose_thr"])]
        return jets

    @zero_handler
    def gettight_jets(self, data):
        jets = data["Jet_select"]
        jets = jets[(jets.disTauTag_score1 >= self.config["tight_thr"])]
        return jets

    @zero_handler
    def pfcand_valid(self, data):
        pfCands = data["PFCand"]
        is_good = (
            (pfCands.pt > self.config["pfcand_pt"])
            & (pfCands.eta < self.config["pfcand_eta_max"])
            & (pfCands.eta > self.config["pfcand_eta_min"])
            # & (pfCands[self.config["track"]])
        )
        pfCands_selected = pfCands[is_good]
        sort_idx = ak.argsort(pfCands_selected.pt, axis=-1, ascending=False)
        return pfCands_selected[sort_idx]

    @zero_handler
    def match_jet_to_pfcand(self, data, jet_name=None, pf_name=None, dR=0.4):
        jets = data[jet_name]
        pfcands = data[pf_name]

        # shape: events x jets x pfcands
        pf_b = pfcands[:, None, :]                       # broadcast over jets
        deta = jets.eta[:, :, None] - pf_b.eta
        dphi = (jets.phi[:, :, None] - pf_b.phi + np.pi) % (2 * np.pi) - np.pi
        dr = np.sqrt(deta**2 + dphi**2)

        # pfcands duplicated per jet, then masked by dR
        pf_per_jet = ak.broadcast_arrays(jets.eta[:, :, None], pf_b)[1]
        return pf_per_jet[dr < dR]
    
    @zero_handler
    def get_matched_pfCands(self, data, match_object, dR=0.4):
        pfCands = self.match_jet_to_pfcand(data, jet_name=match_object, pf_name="PfCands", dR=dR)
        pfCands_lead = ak.firsts(pfCands, axis=-1)
        # pfCands_lead["dxysig"] = pfCands_lead.dxy / pfCands_lead.dxyError
        # pfCands_lead["ip3d"] = np.sqrt(pfCands_lead.dxy**2 + pfCands_lead.dz**2)
        # pfCands_lead["dxy_weight"] = ak.mean(pfCands.dxy, weight=pfCands.pt, axis=-1)
        # pfCands_lead["dxysig_weight"] = ak.mean(pfCands.dxy / pfCands.dxyError, weight=pfCands.pt, axis=-1)
        return pfCands_lead
    
    @zero_handler
    def set_jet_dxy(self, data):
        jets = data["Jet_select"]
        # JetCHS already contains these quantities
        jets["dxy"] = np.abs(jets.dxy)
        jets["dz"] = np.abs(jets.dz)
        jets["dxy_err"] = np.abs(jets.dxyerror)
        jets["dz_err"] = np.abs(jets.dzerror)
        jets["dxysig"] = jets.dxy / jets.dxy_err

        jets = jets[
            jets.dxy >= self.config["jet_dxy_min"]
        ]
        # # Mask jets with dxy nan (no selected pfcands matching)
        # bad_jets = ak.is_none(data["Jet_lead_pfcand"].dxy, axis=-1)
        # jets = ak.mask(jets, ~bad_jets) # mask bad jets to keep coorect shape
        # lead_pf = data["Jet_lead_pfcand"]
        # jets["dz"] = np.abs(lead_pf.dz)
        # jets["dxy"] = np.abs(lead_pf.dxy)
        # jets["dxy_weight"] = np.abs(lead_pf.dxy_weight)
        # jets["dxysig"] = np.abs(lead_pf.dxysig)
        # jets["dxysig_weight"] = np.abs(lead_pf.dxysig_weight)
        # jets["ip3d"] = lead_pf.ip3d
        # #### Extra variables for the jet dxy cut study #### begin
        # jets["dz_err"] = np.abs(lead_pf.dzError)
        # jets["dxy_err"] = np.abs(lead_pf.dxyError)
        # jets["vxy"] = np.sqrt(lead_pf.vx**2 + lead_pf.vy**2)
        # jets["vz"] = np.abs(lead_pf.vz)
        # jets["fromPV"] = lead_pf.fromPV
        # jets["lostInnerHits"] = lead_pf.lostInnerHits
        # #### end
        # jets = jets[~bad_jets] # remove bad jets
        # jets = jets[jets.dxy >= self.config["jet_dxy_min"]]
        return jets
    
    @zero_handler
    def sum_jj(self, data):

        jets = data['Jet_select']

        j1 = jets[:, 0]
        j2 = jets[:, 1]

        mt1 = np.sqrt(j1.mass**2 + j1.pt**2)
        mt2 = np.sqrt(j2.mass**2 + j2.pt**2)

        mjj2 = (
            j1.mass**2 + j2.mass**2
            + 2 * (
                mt1 * mt2 * np.cosh(j1.eta - j2.eta)
                - j1.pt * j2.pt * np.cos(j1.phi - j2.phi)
            )
        )

        return np.sqrt(np.maximum(mjj2, 0))

    @zero_handler
    def dphi_min_cut(self, data):
        return abs(data["dphi_jet1_jet2"]) > self.config["dphi_j1j2"]

    
    @zero_handler
    def b_tagged_jet_cut(self, data, name="Jet"):
        # To leading jets are excluded! 
        jets = data[name]
        # Jet_btagDeepFlavB satisfies the Medium (>0.2783) WP:
        # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18
        # b_tagged_idx = (jet_not_signal.btagDeepFlavB > 0.2783)
        # return ak.num(jet_not_signal[b_tagged_idx], axis=-1) == 0
        # print(jets.btagged)
        # print(ak.num(jets[jets.btagged], axis=1))
        return ak.num(jets[jets.btagged], axis=1) == 0
    
    @zero_handler
    def set_njets_pass(self, data):
        jets_score = data["Jet_select"].disTauTag_score1
        n_pass = []
        for score in self.config["score_pass"]:
            jets_pass = jets_score[(jets_score>=score)]
            passed = ak.num(jets_pass, axis=1)
            n_pass.append( passed )
        n_pass = ak.from_regular(
            np.stack(n_pass, axis=1), axis=-1)
        tight_wp = self.config["score_pass"].index(self.config["tight_thr"])
        tight_pass = n_pass[:,tight_wp]
        tight_bin0 = (n_pass[:,tight_wp] == 0)
        tight_bin1 = (n_pass[:,tight_wp] == 1)
        tight_bin2 = (n_pass[:,tight_wp] == 2)
        return {
            "n_pass" : n_pass,
            "n_pass_score_bin" : ak.local_index(n_pass, axis=1),
            "tight_pass" : tight_pass,
            "tight_bin0" : tight_bin0,
            "tight_bin1" : tight_bin1,
            "tight_bin2" : tight_bin2
        }
    
    @zero_handler
    def set_njets_pass_finebin(self, data):
        jets_score = data["Jet_select"].disTauTag_score1
        n_pass = []
        for score in self.config["score_pass_finebin"]:
            jets_pass = jets_score[(jets_score>=score)]
            passed = ak.num(jets_pass, axis=1)
            n_pass.append( passed )
        n_pass = ak.from_regular(
            np.stack(n_pass, axis=1), axis=-1)
        # print(n_pass)
        # print(ak.local_index(n_pass, axis=1))
        return {
            "n_pass_finebin" : n_pass,
            "n_pass_score_bin_finebin" : ak.local_index(n_pass, axis=1),
        }
    
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

    @zero_handler
    def predict_yield(self, data, weight=None):

        jets = data["Jet_select"]
        
        weights_bin0to1 = []
        weights_bin0to2 = []
        weights_bin1to2 = []
        
        for score in self.config["score_pass"]:
            
            # fake =  self.config["jet_fake_rate"](jet_pt=jets.pt)
            # fake =  self.config["jet_fake_rate"](jet_dxy=jets.dxy)
            fake =  self.config["jet_fake_rate"](jet_pt=jets.pt, jet_dxy=jets.dxy)
            
            # from bin 0 to bin 1 and 2
            events_0tag = (ak.num(jets[(jets.disTauTag_score1 >= score)]) == 0)
            masked_jets = ak.broadcast_arrays(jets.disTauTag_score1, True)[1]
            masked_jets = events_0tag * masked_jets
            fake0 = ak.mask(fake, masked_jets)
            f_1, f_2 = fake0[:,0], fake0[:,1]
            from0to1 = ( f_1*(1-f_2) + f_2*(1-f_1) ) / ((1-f_2)*(1-f_1))
            from0to2 = ( f_1*f_2 ) / ((1-f_2)*(1-f_1))
            from0to1 = ak.fill_none(from0to1, 0.0)
            from0to2 = ak.fill_none(from0to2, 0.0)
            weights_bin0to1.append(from0to1)
            weights_bin0to2.append(from0to2)
            
            # from bin 1 to bin 2
            events_1tag = (ak.num(jets[(jets.disTauTag_score1 >= score)]) == 1)
            masked_jets = ak.broadcast_arrays(jets.disTauTag_score1, True)[1]
            masked_jets = events_1tag * masked_jets
            fake1 = ak.mask(fake, masked_jets)
            f_1, f_2 = fake1[:,0], fake1[:,1]
            from1to2 = ( f_1*f_2 ) / (f_1*(1-f_2) + f_2*(1-f_1))
            from1to2 = ak.fill_none(from1to2, 0.0)
            weights_bin1to2.append(from1to2)
            
        yield_bin0to1 = ak.from_regular(np.stack(weights_bin0to1, axis=1), axis=-1)
        yield_bin0to2 = ak.from_regular(np.stack(weights_bin0to2, axis=1), axis=-1)
        yield_bin1to2 = ak.from_regular(np.stack(weights_bin1to2, axis=1), axis=-1)
            
        # now we need to each predicted yield assign cooresponding score bin
        score_bin = ak.local_index(yield_bin0to1, axis=1) + 1 # +1 because we skip first bin
        
        # One of the WP (called tight) is used for the final analysis
        # the weight for this will be saved saparetly to have one weight per event
        tight_wp = self.config["score_pass"].index(self.config["tight_thr"])

        return_cols = {
            "yield_bin0to1" : weight*yield_bin0to1,
            "yield_bin1to2" : weight*yield_bin1to2,
            "yield_bin0to2" : weight*yield_bin0to2,
            "tight_yield_bin0to1" : weight*yield_bin0to1[:,tight_wp],
            "tight_yield_bin1to2" : weight*yield_bin1to2[:,tight_wp],
            "tight_yield_bin0to2" : weight*yield_bin0to2[:,tight_wp],
            "score_bin"     : score_bin
        }

        if self.config["compute_systematics"]:
            for sys in ["stat_up","stat_down", "sys_up", "sys_down"]:
                
                weights_bin0to1 = []
                weights_bin0to2 = []
                weights_bin1to2 = []
                
                for score in self.config["score_pass"]:
            
                    if sys == "stat_up":
                        fake = self.config["jet_fake_rate"](jet_pt=jets.pt, jet_dxy=jets.dxy, variation="up")
                    elif sys == "stat_down":
                        fake = self.config["jet_fake_rate"](jet_pt=jets.pt, jet_dxy=jets.dxy, variation="down")
                    elif sys == "sys_up": # +10% systematic uncertainty for the fake rate
                        fake = self.config["jet_fake_rate"](jet_pt=jets.pt, jet_dxy=jets.dxy) * 1.1
                    elif sys == "sys_down": # -10% systematic uncertainty for the fake rate
                        fake = self.config["jet_fake_rate"](jet_pt=jets.pt, jet_dxy=jets.dxy) * 0.9
                    else:
                        raise ValueError("Unknown systematic uncertainty")
                    
                    # from bin 0 to bin 1 and 2
                    events_0tag = (ak.num(jets[(jets.disTauTag_score1 >= score)]) == 0)
                    masked_jets = ak.broadcast_arrays(jets.disTauTag_score1, True)[1]
                    masked_jets = events_0tag * masked_jets
                    fake0 = ak.mask(fake, masked_jets)
                    f_1, f_2 = fake0[:,0], fake0[:,1]
                    from0to1 = ( f_1*(1-f_2) + f_2*(1-f_1) ) / ((1-f_2)*(1-f_1))
                    from0to2 = ( f_1*f_2 ) / ((1-f_2)*(1-f_1))
                    from0to1 = ak.fill_none(from0to1, 0.0)
                    from0to2 = ak.fill_none(from0to2, 0.0)
                    weights_bin0to1.append(from0to1)
                    weights_bin0to2.append(from0to2)
                    
                    # from bin 1 to bin 2
                    events_1tag = (ak.num(jets[(jets.disTauTag_score1 >= score)]) == 1)
                    masked_jets = ak.broadcast_arrays(jets.disTauTag_score1, True)[1]
                    masked_jets = events_1tag * masked_jets
                    fake1 = ak.mask(fake, masked_jets)
                    f_1, f_2 = fake1[:,0], fake1[:,1]
                    from1to2 = ( f_1*f_2 ) / (f_1*(1-f_2) + f_2*(1-f_1))
                    from1to2 = ak.fill_none(from1to2, 0.0)
                    weights_bin1to2.append(from1to2)
                    
                yield_bin0to1 = ak.from_regular(np.stack(weights_bin0to1, axis=1), axis=-1)
                yield_bin0to2 = ak.from_regular(np.stack(weights_bin0to2, axis=1), axis=-1)
                yield_bin1to2 = ak.from_regular(np.stack(weights_bin1to2, axis=1), axis=-1)

                tight_yield_bin0to1 = yield_bin0to1[:,tight_wp]
                tight_yield_bin1to2 = yield_bin1to2[:,tight_wp]
                tight_yield_bin0to2 = yield_bin0to2[:,tight_wp]

                return_cols["tight_yield_bin0to1_"+sys] = weight*tight_yield_bin0to1
                return_cols["tight_yield_bin1to2_"+sys] = weight*tight_yield_bin1to2
                return_cols["tight_yield_bin0to2_"+sys] = weight*tight_yield_bin0to2
        
        # print(return_cols.keys())
        return return_cols
        
    # Adding gen-jet flavor to every jet:
    def gen_lep(self, data):
        if len(data) == 0:
            return {"gen_mu":ak.Array([]),
                    "gen_ele":ak.Array([])}
        gen_mu = data.GenPart[
            (abs(data.GenPart.pdgId) == 13)
            & data.GenPart.hasFlags(["isHardProcess"])
            & data.GenPart.hasFlags(["isFirstCopy"])
        ]
        gen_ele = data.GenPart[
            (abs(data.GenPart.pdgId) == 11)
            & data.GenPart.hasFlags(["isHardProcess"])
            & data.GenPart.hasFlags(["isFirstCopy"])
        ]
        return {
            "gen_mu"  : gen_mu,
            "gen_ele" : gen_ele
        }
        
    @zero_handler
    def jet_updated_flavour(self, data):
        
        jets = data["Jet_select"]
        
        # add hadronic tau flavour:
        tau_vis = data.GenVisTau[ ((data.GenVisTau.pt > 30) & (abs(data.GenVisTau.eta) < 2.4) &
                                  (data.GenVisTau.parent.hasFlags(["fromHardProcess"])))
                                ]
        matches_tauhad, _ = jets.nearest(tau_vis, return_metric=True, threshold=0.4)
        matches_mu, _  = jets.nearest(data["gen_mu"], return_metric=True, threshold=0.4)
        matches_ele, _ = jets.nearest(data["gen_ele"], return_metric=True, threshold=0.4)
        
        
        updFlavour = ak.where(~ak.is_none(matches_tauhad, axis=1), 15, jets.partonFlavour)
        updFlavour = ak.where(~ak.is_none(matches_mu, axis=1), 13, updFlavour)
        updFlavour = ak.where(~ak.is_none(matches_ele, axis=1), 11, updFlavour)
        
        return np.abs(updFlavour)

    @zero_handler
    def binning_schema(self, data):
        jets = data["Jet_select"]
        # declear variables for binning
        met = data["PFMET"].pt
        jet2_pt = jets[:,1].pt
        mt2 = data["mt2_j1_j2_MET"]
        # create empty binning
        bins = np.full((len(met)), np.nan)
        # B1 = (jet2_pt < 50) & (met >= 250)
        # B2 = (jet2_pt < 50) & (met < 250) & (mt2 < 100)
        # B3 = (jet2_pt < 50) & (met < 250) & (mt2 >= 100)
        # B4 = (jet2_pt >= 50) & (jet2_pt < 100)
        # B5 = (jet2_pt >= 100)
        # bins[B1] = 1
        # bins[B2] = 2
        # bins[B3] = 3
        # bins[B4] = 4
        # bins[B5] = 5
        B1 = (jet2_pt < 50) & (met < 250) & (mt2 < 100)
        B2 = (jet2_pt < 50) & (met < 250) & (mt2 >= 100)
        B3 = (jet2_pt < 50) & (met >= 250)
        B4 = (jet2_pt >= 50) & (jet2_pt < 100) & (met < 250) & (mt2 < 100)
        B5 = (jet2_pt >= 50) & (jet2_pt < 100) & (met < 250) & (mt2 >= 100)
        B6 = (jet2_pt >= 50) & (jet2_pt < 100) & (met >= 250)
        B7 = (jet2_pt >= 100) & (mt2 < 100)
        B8 = (jet2_pt >= 100) & (mt2 >= 100)
        bins[B1] = 1
        bins[B2] = 2
        bins[B3] = 3
        bins[B4] = 4
        bins[B5] = 5
        bins[B6] = 6
        bins[B7] = 7
        bins[B8] = 8

        return bins