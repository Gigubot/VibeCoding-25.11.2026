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
| **DEMO** (default) | ~2,300 detected proteins from the public *S. aureus* NCTC8325 reference proteome (UniProt UP000008816), all **27 systematic contrasts**, **simulated values** | Protein identities, pathways and keywords are fetched live from UniProt; values are generated in `app.js`. Offline fallback: `demo_data.js` |
| **REAL** | **Every protein (or phosphosite) and every contrast (33)** of the course data | Click **📂 Load course folder** and choose your local course `data/` or `data/derived/` folder |

The course data are unpublished. In REAL mode the files are decompressed and read **inside your browser**. They are never uploaded and never written into this repo. `.gitignore` blocks data files.

Files used in REAL mode, all from the course `derived/` folder unless noted:

| File | Needed for |
|---|---|
| `dea_total.tsv.gz` / `dea_phospho_STYH.tsv.gz` / `dea_phospho_STY.tsv.gz` | Volcano and all games. **Required** for the chosen dataset |
| `protein_annotation.tsv.gz` | Gene names, products, UniProt orthologs (`staa8_uniprot`, used for pathways and AlphaFold) |
| `abundances_<dataset>.tsv.gz` | Per-sample dot plots on protein cards |
| `sites_STYH.tsv.gz` | Phosphosites highlighted on the sequence and the 3D model |
| `raw/resources/JE2_proteins.fasta` | Real JE2 sequences (otherwise the NCTC8325 reference from UniProt is shown) |

## What's inside (v0.2)

**Live volcano**
- Contrast picker with all contrasts grouped by type (medium effects, genotype differences, interactions), FDR and |log2 effect| sliders, and draggable dashed cutoff lines.
- **Show-on-plot switches** with counts: ▲ up hits · ▼ down hits · not selected · ✕ `lod_imputed` estimates. Click a switch to hide or show that group.
- Up/down counts in plain English, CSV export of hits, SVG export via the camera icon.

**Games, sound and memes**

| Game | How it plays | Teaches |
|---|---|---|
| 🎯 **Sign Sniper** | Arcade: 3 lives, 8-second timer per shot, combos multiply points, keyboard ← →, wrong or slow shots cost XP, high score. Interaction questions (×2) unlock at combo ×4 | `A_vs_B` means `diff = log2(A) − log2(B)` |
| 🧪 **Family Catcher** | Protein bubbles fall; catch members of a family (UniProt keyword such as *Purine biosynthesis*, *Ribosomal protein*, *Virulence*) with a beaker and dodge impostors. Colour = direction in the current contrast, gold ring = passes your cutoffs (bonus), 🛡 shield power-up, rising speed | Which proteins belong together, and what the current contrast does to them |
| 🗺️ **Pathway Explorer** | Category → pathway → protein → **AlphaFold 3D model**, with ▲/▼ hit counts per pathway. The pathway is highlighted on the volcano | Pathway-level view of the data |
| 🪤 **LOD Trap** | Clicking a purple × asks what `lod_imputed` means | An imputed estimate is presence/absence evidence, not a precise fold change |
| 📉 **False-discovery meter** | Expected false hits ≤ FDR cutoff × hits, live | Multiple testing ("Reformed p-hacker" badge) |
| 📖 **Pokédex** | Every protein you open is collected; 🧊 structure viewed; 👑 Stk1/Stp1 | Exploring the dataset |

- **Sound effects** synthesised with the Web Audio API (no audio files): XP chime rising with combo, level-up fanfare, badge bell, XP-loss "wah-wah", level-down, buzzer, countdown ticks, game over, meme whoosh. 🔊/🔇 button in the header.
- **Flying memes** on streaks and high scores (text and emoji only, e.g. 🧥 "THE MORE YOU BUY, THE MORE YOU SAVE", 📈 STONKS, 💀 SKILL ISSUE).
- **Engine:** 6 levels (you can also level down), 16 badges, and a leaderboard with high scores for each player name, saved in this browser (`localStorage`).

**Protein card:** per-sample values for all 27 samples (NA, never 0), every contrast with its meaning in words, UniProt pathways and keywords, links to UniProt/AlphaFold DB/STRING, sequence viewer with regex motif search and measured phosphosites, and the AlphaFold 3D model coloured by pLDDT with phosphosites as pink spheres.

### The two formulas

- Simple contrast: `diff = log2(A) − log2(B)`. For `TSB_vs_SN_at_anc`, positive means higher in TSB, which is **lower in PASN**.
- Interaction `interaction_A_vs_B_at_TSB_vs_SN`: `(A_TSB − B_TSB) − (A_SN − B_SN) = R_B − R_A`, with `R_g = log2(SN_g) − log2(TSB_g)`. Positive means B's change from TSB to SN is more positive than A's. It does **not** simply mean "higher in B".

### Caveats shown in the app

- AlphaFold models are the **NCTC8325 reference ortholog**, so JE2 residue numbers can differ. The sequence view flags mismatched sites.
- The `lod_imputed` flag in the course export applies to the whole model, so check the per-sample dots.
- PASN is cell-free spent medium of *P. aeruginosa*, not co-culture. The `stp` mutant also has strongly reduced Stk1.
- Pathways and keywords come from the UniProt NCTC8325 reference proteome, matched to JE2 proteins through the course ortholog column; proteins without an ortholog have no pathway or structure.

## Files

| File | Role |
|---|---|
| `index.html` | Layout and styles. Loads Plotly and 3Dmol.js from their CDNs |
| `app.js` | Data layer (UniProt annotation, demo generator, local loader), volcano, sound engine, games, protein card |
| `demo_data.js` | Offline fallback: public protein identities used only when UniProt is unreachable |

## Sprint checklist / next ideas

- [x] v0.1: repo, volcano, games, local data loader
- [x] v0.2: sound effects, arcade Sign Sniper with memes, Family Catcher, Pathway Explorer with 3D, all contrasts, show/hide switches
- [ ] Compare mode: two volcanoes side by side (e.g. anc vs stk PASN response)
- [ ] Phospho vs parent protein: site effect minus protein effect
- [ ] Shared leaderboard (e.g. a GitHub Gist) instead of per-browser storage
- [ ] Daily Mystery Protein with a fixed random seed per date

Built with AI assistance (Claude). Course data © FGCZ/authors, unpublished; not included. Code: MIT.
