# 🧬 JE2 Proteome Quest

A gamified, zero-install explorer for the Bio253 (UZH/ETHZ, HS2026) *Staphylococcus aureus* JE2 proteome and phosphoproteome. Built in the Wednesday vibe-coding sprint.

**Play:** clone the repo and double-click `index.html`. Nothing to install. You need internet for the plotting and 3D libraries and for AlphaFold.

```bash
git clone https://github.com/Gigubot/VibeCoding-25.11.2026.git
cd VibeCoding-25.11.2026
# open index.html in Chrome, Edge or Firefox
```

## Two data modes

| Mode | What you see | Where the data come from |
|---|---|---|
| **DEMO** (default) | 2,300 proteins, 12 contrasts, **simulated values** | Generated in `app.js`. Some public protein identities in `demo_data.js` (JE2 locus tags, NCTC8325 UniProt IDs) make it realistic |
| **REAL** | **Every protein (or phosphosite) and every contrast** of the course data | Click **📂 Load course folder** and choose your local course `data/` or `data/derived/` folder |

The course data are unpublished. In REAL mode the files are decompressed and read **inside your browser**. They are never uploaded and never written into this repo. `.gitignore` blocks data files.

Files used in REAL mode, all from the course `derived/` folder unless noted:

| File | Needed for |
|---|---|
| `dea_total.tsv.gz` / `dea_phospho_STYH.tsv.gz` / `dea_phospho_STY.tsv.gz` | Volcano and all games. **Required** for the chosen dataset |
| `protein_annotation.tsv.gz` | Gene names, products, UniProt orthologs (used by AlphaFold), protein families |
| `abundances_<dataset>.tsv.gz` | Per-sample dot plots on protein cards |
| `sites_STYH.tsv.gz` | Phosphosites highlighted on the sequence and the 3D model |
| `raw/resources/JE2_proteins.fasta` | Real JE2 sequences (otherwise the NCTC8325 reference from UniProt is shown) |

## What's inside

**Live volcano**
- Contrast picker, FDR and |log2 effect| sliders, and draggable dashed cutoff lines.
- Up/down counts written in plain English.
- Toggle for `lod_imputed` estimates, CSV export of hits, SVG export through the camera icon.

**Games and XP**

| Game | Teaches |
|---|---|
| 🎯 **Sign Sniper**: which group has more? After a 5-answer streak, interaction questions unlock | `A_vs_B` means `diff = log2(A) − log2(B)` |
| 🪤 **LOD Trap**: clicking a purple × asks what `lod_imputed` means | An imputed estimate is presence/absence evidence, not a precise fold change |
| 📉 **False-discovery meter**: expected false hits ≤ FDR cutoff × hits | Multiple testing. The "Reformed p-hacker" badge rewards going back to FDR ≤ 0.1 |
| 🕵️ **Mystery Protein**: five clues, and fewer clues earn more XP | Linking a position on the volcano to biology |
| 🧬 **Family Hunt**: catch the members of a protein family that pass your cutoffs | Families are defined from product names, so it works for every protein. An empty hunt shows what GSEA can find and thresholds can't |
| 📖 **Protein Pokédex**: every protein card you open is collected; 🧊 marks a viewed structure, 👑 marks Stk1/Stp1 | Exploring the dataset |

**Protein card:**
- Per-sample values for all 27 samples (NA, never 0).
- Every contrast with its meaning in words.
- Links to UniProt, AlphaFold DB and STRING.
- Sequence viewer with regex motif search (for example `[ST]P` or `R..S`) and measured phosphosites.
- **AlphaFold 3D model** (3Dmol.js), coloured by pLDDT, with phosphosites as pink spheres.

**Engine:** 6 levels, 14 badges, and a leaderboard for each player name, saved in this browser (`localStorage`).

### The two formulas

- Simple contrast: `diff = log2(A) − log2(B)`. For `TSB_vs_SN_at_anc`, positive means higher in TSB, which is **lower in PASN**.
- Interaction `interaction_A_vs_B_at_TSB_vs_SN`: `(A_TSB − B_TSB) − (A_SN − B_SN) = R_B − R_A`, with `R_g = log2(SN_g) − log2(TSB_g)`. Positive means B's change from TSB to SN is more positive than A's. It does **not** simply mean "higher in B".

### Caveats shown in the app

- AlphaFold models are the **NCTC8325 reference ortholog**, so JE2 residue numbers can differ. The sequence view flags mismatched sites.
- The `lod_imputed` flag in the course export applies to the whole model, so check the per-sample dots.
- PASN is cell-free spent medium of *P. aeruginosa*, not co-culture. The `stp` mutant also has strongly reduced Stk1.

## Files

| File | Role |
|---|---|
| `index.html` | Layout and styles. Loads Plotly and 3Dmol.js from their CDNs |
| `app.js` | Data layer (demo generator and local loader), volcano, games, protein card |
| `demo_data.js` | Public protein identities used only to label the simulated demo |

## Sprint checklist / next ideas

- [x] Repo plus spec, working skeleton, games, local data loader
- [ ] Compare mode: two volcanoes side by side (e.g. anc vs stk PASN response)
- [ ] Phospho vs parent protein: site effect minus protein effect
- [ ] Shared leaderboard (e.g. a GitHub Gist) instead of per-browser storage
- [ ] Sound effects (Web Audio) on streaks and level-ups
- [ ] Daily Mystery Protein with a fixed random seed per date

Built with AI assistance (Claude). Course data © FGCZ/authors, unpublished; not included. Code: MIT.
