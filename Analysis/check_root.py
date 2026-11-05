import uproot

f = uproot.open("/pnfs/desy.de/cms/tier2/store/user/svashish/LLSTauProduction_RUN3/SUS-RunIIISummer2024_EXONANO_stau700_lsp1_ctau10m_v3/SUS-Run3Summer2024-stau700_lsp1_ctau10m/NANOAOD/260820_162504/0000/SUS-Run3Summer2024-stau700_lsp1_ctau10m_ExoNanoAOD_100.root")

print(f.keys())

tree = f["Events"]

print(f"\nNumber of events: {tree.num_entries}")
print(f"Number of branches: {len(tree.keys())}\n")

for branch in tree.keys():
    print(branch)
