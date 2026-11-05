# GROMACS Protein–Ligand MD Automation

Bash/Python scripts that take a docked protein–ligand complex through system setup, equilibration and production MD in **GROMACS 2023–2026**, then run a complete, quality-controlled analysis that ends in Excel workbooks, 600-dpi figures and a QC report.

The workflow follows the standard GROMACS protein–ligand protocol (CHARMM27 force field, TIP3P water, SwissParam ligand parameters). What this repository adds is the automation, the robustness checks and the analysis/reporting pipeline.

---

## What it does

```
REC.pdb + SwissParam ligand ──► setup_md_swissparam.sh ──► run_md_protocol.sh ──► run_md_analysis.sh
  (from docking)                  topology, box, water,      EM → NVT → NPT →        RMSD, RMSF, Rg, H-bonds,
                                  ions, restraints, index,   production MD           SASA, contacts, DSSP, PCA,
                                  EM.tpr + grompp tests      (resumes from .cpt)     interaction energy → Excel,
                                                                                     PNG figures, QC report
```

| Script | Purpose |
|---|---|
| `scripts/setup_md_swissparam.sh` | Builds the solvated, neutralised complex from `REC.pdb` and the SwissParam output (`unk.itp`, `unk.pdb`). Renames the ligand to `LIG`, merges coordinates, edits `topol.top`, builds the box, water, ions, ligand position restraints and the `Protein_LIG` index group, writes `EM.tpr`, and test-runs `grompp` for NVT, NPT and MD so input errors show up before any simulation starts. |
| `scripts/run_md_protocol.sh` | Runs EM → NVT → NPT → production MD, skips stages that already finished, resumes an interrupted production run from `MD.cpt`, and starts the analysis when MD is done. Can wait for another running `mdrun` (`--wait`). |
| `scripts/run_md_analysis.sh` | Analyses a finished run. The simulation folder is only read; everything goes to `MD_Analysis_<N>ns/`. |
| `scripts/make_results.py` | Builds the Excel files, figures, CSV copies and QC report from the GROMACS outputs. Called by `run_md_analysis.sh`. |
| `scripts/energy_pick.py` | Picks `gmx energy` terms by name, never by menu number. |
| `scripts/package_results.sh` | Collects every QC-passed analysis under a root folder (`<root>/<target>/<simulation>/MD_Analysis_<N>ns`) into one folder with results, the complete analysis and the simulation inputs, and zips it for sharing. |
| `mdp_examples/` | Example `.mdp` files (CHARMM settings, V-rescale, C-rescale, 200 ns production). Check them against your own protocol before use. |

## Key design choices

- **Groups by name, never by number.** Index groups are built with `gmx select`; energy terms are chosen from the `gmx energy` menu by name.
- **Ligand position restraints live inside `LIG.itp`.** In a multi-chain protein, restraints placed at the end of `topol.top` would attach to the wrong molecule.
- **Robust PBC treatment for multi-chain complexes.** First frame made whole → protein + ligand clustered → frames every 250 ps → `-pbc nojump` against the clustered reference → `-pbc mol -center -ur compact` for ligand/complex analyses. This avoids split chains and false "ligand left the pocket" events.
- **SwissParam topologies are handled.** SwissParam `[ atomtypes ]` often lack atomic numbers, so `gmx hbond` would see a ligand with no donors or acceptors. Analysis-only copies of the topology get atomic numbers from the atom masses; your original files are never changed.
- **Interaction energy by rerun.** Energy groups are not supported on GPU runs, so protein–ligand Coulomb and Lennard-Jones energies are recomputed with `mdrun -rerun` on the CPU.
- **Nothing is invented.** Every value comes from GROMACS output. The molecular surface (SES, `SurfMol`) cannot be computed with GROMACS and is left empty on purpose.
- **Quality control.** Before any file is trusted, the QC step checks: identical 250 ps time grid across all analyses, no duplicate frames, RMSD ≈ 0 at t = 0, no chain-split jumps in RMSD or Rg, integer non-negative H-bond counts, non-negative SASA, and that every Excel cell matches the raw `.xvg` value. Ligand RMSD > 1 nm and ligand leaving contact are reported as warnings: they are scientific results, not errors.

## Requirements

- Linux, GROMACS **2023 or newer** (uses the selection-based `gmx hbond` and built-in `gmx dssp`); developed with GROMACS 2026.3 (conda-forge, CUDA)
- Python 3.9+ with `numpy`, `pandas`, `matplotlib`, `openpyxl`
- A ligand topology from [SwissParam](https://www.swissparam.ch/) for the setup script (the analysis script also works with CGenFF, ACPYPE and other ligand topologies)

```bash
conda create -n gmx -c conda-forge gromacs numpy pandas matplotlib openpyxl
conda activate gmx
```

## Usage

Start every session like this (`set +H` stops Bash from altering selections that contain `!`; `LC_NUMERIC=C` keeps decimal points in locales that use commas):

```bash
set +H
conda activate gmx
export LC_NUMERIC=C LC_ALL=C
```

**1. Prepare the inputs** (UCSF Chimera/ChimeraX): from the docked complex, save the ligand with hydrogens as `LIG.mol2` and the DockPrep-processed receptor as `REC.pdb`. Parameterise the ligand on SwissParam and unzip the result into the simulation folder together with `REC.pdb` and your `.mdp` files (`ions.mdp`, `EM.mdp`, `NVT.mdp`, `NPT.mdp`, `MD.mdp`).

**2. Build the system**

```bash
bash scripts/setup_md_swissparam.sh "/path/to/simulation_folder"
# optional: copy force field / water from an earlier run and compare .mdp files
bash scripts/setup_md_swissparam.sh "/path/to/new_folder" "/path/to/reference_folder"
# optional settings
CONC=0.15 BOXD=1.2 BOXTYPE=dodecahedron bash scripts/setup_md_swissparam.sh "/path/to/folder"
```

**3. Run the simulation and the analysis** (keeps running after the terminal closes)

```bash
nohup setsid bash scripts/run_md_protocol.sh "/path/to/simulation_folder" --wait > /dev/null 2>&1 < /dev/null &
tail -f "/path/to/simulation_folder/protocol_log.txt"
```

**4. Analyse an existing run on its own**

```bash
bash scripts/run_md_analysis.sh "/path/to/simulation_folder" 200      # length in ns
LIGAND=UNL bash scripts/run_md_analysis.sh "/path/to/folder" 100       # name the ligand moleculetype explicitly
```

**5. Package all QC-passed results for sharing**

```bash
bash scripts/package_results.sh "/path/to/Simulation_root" 200
```

Check the real trajectory length first: `gmx check -f MD.xtc 2>&1 | tr '\r' '\n' | grep -E "Last frame|Coords"`. The analysis stops if the trajectory is shorter than requested.

## Output (`MD_Analysis_<N>ns/`)

| Folder | Contents |
|---|---|
| `Excel/` | `analysis_<sim>.xlsx` (27-column standard layout in Å, Å², kJ/mol, %), `Analysis_2_<sim>.xlsx` (+ ligand columns, interaction energy, column definitions), `Protein_Ligand_MD_Analysis.xlsx` (detail sheets in nm), `Summary_Statistics.xlsx` (mean, SD, median, min/max with time, last-window statistics, RMSF statistics, binding summary) |
| `Figures_PNG/` | 19 figures at 600 dpi: RMSD, ligand RMSD, RMSF (with chain boundaries and contact residues), Rg, H-bonds, SASA and buried SASA, minimum distance and contacts, secondary structure, interaction energy, energies, temperature/pressure/density, box, PCA, free-energy landscape, contact occupancy |
| `Data_CSV/`, `Raw_XVG/` | CSV copies of every table; raw GROMACS outputs |
| `Structures/` | Complex at 0 ns and at the final frame; input PDBs |
| `Reports/` | `QC_REPORT.txt`, `gmx_commands.txt` (exact commands), `Methods_parameters.txt` (values read from your `.mdp`, topology and log files), `pipeline_log.txt` |

## Limitations

- One ligand molecule type per system; membrane systems (lipids) need `LIGAND=` set explicitly.
- `setup_md_swissparam.sh` expects a protein-only `REC.pdb` (no cofactors, metals or waters) and SwissParam output.
- Protein–ligand interaction energy is short-range Coulomb + Lennard-Jones, **not** a binding free energy. Use MM/PBSA or alchemical methods for affinities.
- One trajectory is one sample. Run replicates before drawing conclusions.
- Membrane proteins (GPCRs, channels, transporters) must be simulated in a lipid bilayer (e.g. built with CHARMM-GUI Membrane Builder); the setup script here builds water-only systems.

## Use so far

The analysis pipeline has been run on eight 200 ns protein–ligand simulations (1.6 µs in total, GROMACS 2026.3 on one RTX 4090), all passing the automatic QC checks.

## Credits

- Core GROMACS protein–ligand protocol: Dr. Dweipayan Goswami ([GROMACS-Protein-Ligand](https://github.com/DweipayanG/GROMACS-Protein-Ligand)) and J. A. Lemkul's GROMACS tutorials (mdtutorials.com).
- Ligand parameters: SwissParam (Zoete et al., *J. Comput. Chem.* 2011; Bugnon et al., *J. Chem. Inf. Model.* 2023).
- GROMACS: Abraham et al., *SoftwareX* 2015.
- `sort_mol2_bonds.pl` and `cgenff_charmm2gmx.py`, used in some ligand preparation workflows, are third-party scripts and are **not** included here; get them from their original sources.

Automation and analysis scripts written by Md. Al Hadi Islam with AI assistance (Claude, Anthropic).

## Citation

If you use these scripts, please cite this repository (see `CITATION.cff`) together with GROMACS, the force field and the ligand parameterisation tool you used.

## License

MIT, see `LICENSE`.
