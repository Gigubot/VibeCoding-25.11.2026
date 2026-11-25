/* JE2 Proteome Quest v0.2 — a gamified explorer for the Bio253 S. aureus JE2 data set.
 * Runs entirely in the browser. Course data are read from YOUR disk and never uploaded.
 * Sign rule used everywhere: contrast A_vs_B  ->  diff = log2(A) - log2(B).
 */
'use strict';

/* ============================ 1. utilities ============================ */
const $ = s => document.querySelector(s);
const esc = s => String(s ?? '').replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
const fmt = (x, d = 2) => (x == null || Number.isNaN(x)) ? 'NA' : (x > 0 ? '+' : '') + Number(x).toFixed(d);
const fmtP = x => (x == null || Number.isNaN(x)) ? 'NA' : (x < 0.001 ? x.toExponential(1) : x.toFixed(3));
const pick = (arr, r = Math.random) => arr[Math.floor(r() * arr.length)];
function mulberry32(a) { return function () { a |= 0; a = a + 0x6D2B79F5 | 0; let t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; }; }
function hashStr(s) { let h = 2166136261; for (const c of s) { h ^= c.charCodeAt(0); h = Math.imul(h, 16777619); } return h >>> 0; }
function gauss(r) { let u = 0, v = 0; while (!u) u = r(); while (!v) v = r(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); }
function erfc(x) { // Numerical Recipes approximation, |error| < 1.2e-7
  const z = Math.abs(x), t = 1 / (1 + 0.5 * z);
  const r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
  return x >= 0 ? r : 2 - r;
}
const pTwoSided = z => Math.min(1, erfc(Math.abs(z) / Math.SQRT2));
function bh(ps) { // Benjamini-Hochberg adjusted p-values
  const n = ps.length, idx = ps.map((p, i) => [p, i]).sort((a, b) => a[0] - b[0]); const q = new Array(n); let min = 1;
  for (let k = n - 1; k >= 0; k--) { const [p, i] = idx[k]; min = Math.min(min, p * n / (k + 1)); q[i] = min; } return q;
}

/* ============================ 2. design + contrast language ============================ */
const GENOTYPES = ['anc', 'stk', 'stp'];
const MEDIA = ['TSB', 'SNpulse', 'SN'];
const GROUPS = GENOTYPES.flatMap(g => MEDIA.map(m => `${g}_${m}`));
const GENO_COL = { anc: '#0072B2', stk: '#E69F00', stp: '#009E73' };
const NICE = { anc: 'ancestor (anc)', stk: 'kinase mutant (stk)', stp: 'phosphatase mutant (stp)', TSB: 'TSB', SN: 'PASN (SN, ~15 h)', SNpulse: 'PASN pulse (15 min)' };
const MEDIUM_PAIRS = [['TSB', 'SN'], ['TSB', 'SNpulse'], ['SNpulse', 'SN']];
const GENO_PAIRS = [['stk', 'anc'], ['stp', 'anc'], ['stp', 'stk']];
const DEMO_CONTRASTS = [
  ...MEDIUM_PAIRS.flatMap(([a, b]) => GENOTYPES.map(g => `${a}_vs_${b}_at_${g}`)),
  ...GENO_PAIRS.flatMap(([a, b]) => MEDIA.map(m => `${a}_vs_${b}_at_${m}`)),
  ...GENO_PAIRS.flatMap(([a, b]) => MEDIUM_PAIRS.map(([m1, m2]) => `interaction_${a}_vs_${b}_at_${m1}_vs_${m2}`))];

function parseContrast(c) {
  let m = c.match(/^interaction_([A-Za-z]+)_vs_([A-Za-z]+)_at_([A-Za-z]+)_vs_([A-Za-z]+)$/);
  if (m) return { type: 'int', A: m[1], B: m[2], M1: m[3], M2: m[4] };
  m = c.match(/^([A-Za-z]+)_vs_([A-Za-z]+)_at_([A-Za-z]+)$/);
  if (m) return { type: 'simple', A: m[1], B: m[2], ctx: m[3] };
  m = c.match(/^([A-Za-z]+)_vs_([A-Za-z]+)$/);
  if (m) return { type: 'simple', A: m[1], B: m[2], ctx: '' };
  return { type: 'other' };
}
function contrastGroup(c) {
  const p = parseContrast(c);
  if (p.type === 'int') return 'Interactions (genotype × medium)';
  if (p.type === 'simple') return MEDIA.includes(p.A) ? 'Medium effects (within a genotype)' : 'Genotype differences (within a medium)';
  return 'Other contrasts';
}
/** Labels for a positive / negative effect, plus the formula and its plain-English reading. */
function contrastInfo(c) {
  const p = parseContrast(c);
  if (p.type === 'int') {
    const R = g => `R<sub>${g}</sub>`;
    return {
      ...p, pos: `${p.B}'s ${p.M1}→${p.M2} change is larger`, neg: `${p.A}'s ${p.M1}→${p.M2} change is larger`,
      formula: `effect = (${p.A}<sub>${p.M1}</sub> − ${p.B}<sub>${p.M1}</sub>) − (${p.A}<sub>${p.M2}</sub> − ${p.B}<sub>${p.M2}</sub>) = ${R(p.B)} − ${R(p.A)}, &nbsp; with R<sub>g</sub> = log2(${p.M2}<sub>g</sub>) − log2(${p.M1}<sub>g</sub>)`,
      english: `Positive: the change from ${NICE[p.M1] || p.M1} to ${NICE[p.M2] || p.M2} is more positive in ${p.B} than in ${p.A}. Negative: the reverse. It does NOT simply mean "higher in ${p.B}".`
    };
  }
  if (p.type === 'simple') {
    const where = p.ctx ? ` (in ${p.ctx})` : '';
    const pasn = (p.A === 'TSB' && /^SN/.test(p.B)) ? ' = lower in PASN' : '';
    return {
      ...p, pos: `higher in ${p.A}${where}${pasn}`, neg: `higher in ${p.B}${where}`,
      formula: `diff = log2(${p.A}) − log2(${p.B})${p.ctx ? ` &nbsp; within ${p.ctx}` : ''}`,
      english: `Positive: more protein in ${NICE[p.A] || p.A}${pasn}. Negative: more in ${NICE[p.B] || p.B}. |diff| = 1 means a 2-fold difference.`
    };
  }
  return { ...p, pos: 'positive effect', neg: 'negative effect', formula: 'diff = log2(first group) − log2(second group)', english: 'Check the contrast definition in the course DEA report.' };
}

/* ============================ 3. biology annotation (UniProt, public) ============================ */
/* Name-based families: a fallback that works offline for every protein. */
const FAMILIES = [
  ['Ribosomal proteins', /ribosomal protein/i],
  ['Aminoacyl-tRNA ligases', /tRNA ligase|tRNA synthetase/i],
  ['Purine biosynthesis', /phosphoribosyl(amine|amino|formyl|glycinamide)|IMP dehydrogenase|IMP cyclohydrolase|adenylosuccinate|GMP synthase|pur operon/i],
  ['ATP synthase', /ATP synthase/i],
  ['Glycolysis', /glyceraldehyde-3-phosphate dehydrogenase|phosphoglycerate (kinase|mutase)|enolase|pyruvate kinase|fructose-bisphosphate|triose-?phosphate isomerase|glucose-6-phosphate isomerase|6-phosphofructokinase/i],
  ['Respiration', /quinol oxidase|cytochrome|NADH dehydrogenase|menaquinone|succinate dehydrogenase/i],
  ['Oxidative stress', /catalase|superoxide dismutase|peroxiredoxin|alkyl hydroperoxide|thioredoxin|glutathione/i],
  ['ABC transporters', /ABC transporter/i],
  ['Cell wall', /penicillin-binding|peptidoglycan|transglycosylase|autolysin|amidase|teichoic|sortase/i],
  ['Iron & heme uptake', /siderophore|ferric|ferrous|ferrichrome|heme|iron/i],
  ['Virulence & toxins', /toxin|hemolysin|leukocidin|immunoglobulin-binding|coagulase|fibronectin-binding|clumping factor|staphylokinase|lipase|nuclease/i],
];
const familyOf = product => { for (const [name, re] of FAMILIES) if (re.test(product || '')) return name; return ''; };

/* UniProt "biological process" keywords, grouped into explorer categories. */
const KW_CATEGORIES = {
  '⚙️ Metabolism': ['Amino-acid biosynthesis', 'Branched-chain amino acid biosynthesis', 'Aromatic amino acid biosynthesis', 'Arginine biosynthesis', 'Histidine biosynthesis', 'Lysine biosynthesis', 'Methionine biosynthesis', 'Diaminopimelate biosynthesis',
    'Purine biosynthesis', 'Pyrimidine biosynthesis', 'Glycolysis', 'Gluconeogenesis', 'Tricarboxylic acid cycle', 'Carbohydrate metabolism', 'Lactose metabolism', 'One-carbon metabolism',
    'Lipid metabolism', 'Lipid biosynthesis', 'Fatty acid biosynthesis', 'Fatty acid metabolism', 'Phospholipid biosynthesis', 'Menaquinone biosynthesis', 'Porphyrin biosynthesis', 'Heme biosynthesis',
    'Molybdenum cofactor biosynthesis', 'Queuosine biosynthesis', 'Carotenoid biosynthesis', 'Nitrate assimilation'],
  '⚡ Energy & respiration': ['Electron transport', 'ATP synthesis', 'Hydrogen ion transport'],
  '🧬 Information processing': ['Protein biosynthesis', 'Ribosomal protein', 'Ribosome biogenesis', 'Aminoacyl-tRNA synthetase', 'tRNA processing', 'rRNA processing', 'Transcription', 'Transcription regulation',
    'DNA replication', 'DNA repair', 'DNA recombination', 'SOS response'],
  '🧱 Cell envelope & division': ['Cell wall biogenesis/degradation', 'Peptidoglycan synthesis', 'Teichoic acid biosynthesis', 'Cell shape', 'Cell division', 'Septation', 'Peptidoglycan-anchor'],
  '🚚 Transport': ['Transport', 'Ion transport', 'Sugar transport', 'Amino-acid transport', 'Iron transport', 'Nickel transport', 'Sodium transport', 'Phosphotransferase system', 'Protein transport'],
  '📡 Signalling & stress': ['Two-component regulatory system', 'Kinase', 'Protein phosphatase', 'Phosphoprotein', 'Stress response', 'Chaperone', 'Antibiotic resistance', 'Toxin-antitoxin system'],
  '☠️ Virulence': ['Virulence', 'Toxin', 'Cytolysis', 'Hemolysis', 'Secreted'],
};
const KW_TO_CAT = new Map(Object.entries(KW_CATEGORIES).flatMap(([cat, list]) => list.map(k => [k, cat])));
const PW_CAT = '🧪 Metabolic pathways (UniProt)';
const FAM_CAT = '🏷️ Name-based families';

const PROTEOME = new Map(); // NCTC8325 reference proteome UP000008816 (public UniProt data)
const PROTEOME_URL = 'https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000008816&fields=accession,gene_primary,gene_oln,protein_name,cc_pathway,keyword&format=tsv';
async function loadProteome() {
  let txt = null;
  try { txt = localStorage.getItem('je2q-proteome-v1'); } catch (e) { /* no storage */ }
  if (!txt) {
    const r = await fetch(PROTEOME_URL); if (!r.ok) throw new Error(`UniProt ${r.status}`); txt = await r.text();
    try { localStorage.setItem('je2q-proteome-v1', txt); } catch (e) { /* too big or blocked: fine */ }
  }
  for (const l of txt.split('\n').slice(1)) {
    if (!l) continue; const [acc, gene, oln, pname, pw, kw] = l.split('\t');
    const seen = new Set(); const pathways = [];
    for (const s of (pw || '').split(/PATHWAY: /).slice(1)) {
      const parts = s.replace(/\{[^}]*\}\.?/g, '').split(';').map(x => x.trim().replace(/\.$/, '')).filter(Boolean);
      if (parts[0] && !seen.has(parts[0])) { seen.add(parts[0]); pathways.push({ top: parts[0], sub: parts.slice(1).join(' · ') }); }
    }
    PROTEOME.set(acc, { acc, gene: gene || '', oln: (oln || '').split(/[\s;]/)[0], name: (pname || '').split(' (')[0], pathways, keywords: (kw || '').split(';').filter(Boolean) });
  }
  return PROTEOME.size;
}
function annotate(e) {
  const p = PROTEOME.get(e.uniprot);
  e.pathways = p ? p.pathways : [];
  e.processes = p ? p.keywords.filter(k => KW_TO_CAT.has(k)) : [];
  if (p && !e.product) e.product = p.name;
  e.family = familyOf(e.product);
}

/* ============================ 4. data layer ============================ */
const DATA = {
  mode: 'demo', dataset: 'demo', isSite: false,
  entities: new Map(), rows: new Map(), index: new Map(), contrasts: [],
  abund: new Map(), sites: new Map(), seqs: new Map(), files: null, demoMeans: null,
};
let EXPLORER = new Map(); // category -> Map(pathway -> {ids:Set, sub:Map})
let CATCH = new Map();    // group -> [ids]
function indexRows() {
  DATA.index = new Map(); DATA.contrasts = [...DATA.rows.keys()];
  for (const [c, arr] of DATA.rows) DATA.index.set(c, new Map(arr.map(r => [r.id, r])));
}

/* ---- 4a. simulated demo on top of the real reference proteome ---- */
const KW_SHIFT = { 'Ribosomal protein': { R: -1.2, stk: 0, stp: 0 }, 'Purine biosynthesis': { R: -0.8, stk: -1.0, stp: -0.2 }, 'Electron transport': { R: 0.6, stk: 0.9, stp: 0 },
  'Stress response': { R: 1.0, stk: 0, stp: 0.4 }, 'Iron transport': { R: 1.4, stk: 0, stp: 0 }, 'Virulence': { R: -1.3, stk: 0.3, stp: -0.5 },
  'Aminoacyl-tRNA synthetase': { R: -0.4, stk: 0, stp: -0.6 }, 'Amino-acid biosynthesis': { R: 0.9, stk: 0.2, stp: 0 }, 'Glycolysis': { R: 0.7, stk: 0, stp: 0.3 } };
const FAM_SHIFT = { 'Ribosomal proteins': KW_SHIFT['Ribosomal protein'], 'Purine biosynthesis': KW_SHIFT['Purine biosynthesis'], 'Respiration': KW_SHIFT['Electron transport'],
  'Oxidative stress': KW_SHIFT['Stress response'], 'Iron & heme uptake': KW_SHIFT['Iron transport'], 'Virulence & toxins': KW_SHIFT['Virulence'], 'Aminoacyl-tRNA ligases': KW_SHIFT['Aminoacyl-tRNA synthetase'] };
const SYN_TEMPLATES = [i => `50S ribosomal protein L${i}`, i => `30S ribosomal protein S${i}`, () => 'ATP synthase subunit', () => 'ABC transporter ATP-binding protein', () => 'penicillin-binding protein',
  () => 'thioredoxin', () => 'iron ABC transporter substrate-binding protein', () => 'leukocidin family toxin', () => 'cytochrome aa3 quinol oxidase subunit', () => 'enolase',
  () => 'hypothetical protein', () => 'hypothetical protein', () => 'transcriptional regulator', () => 'MFS transporter', () => 'membrane protein', () => 'hydrolase', () => 'oxidoreductase', () => 'hypothetical protein'];

function demoList(r) {
  const list = [];
  if (PROTEOME.size) {
    for (const p of PROTEOME.values()) {
      if (r() < 0.2 && !/^Q2FZ6[45]$/.test(p.acc)) continue; // ~20 % not detected, like a real run
      const special = { Q2FZ64: 'Stk1', Q2FZ65: 'Stp1' }[p.acc];
      list.push({ id: p.acc, key: p.acc, uniprot: p.acc, locus: p.oln, name: special || p.gene || p.oln || p.acc, product: p.name });
    }
  } else { // offline fallback: public identities + synthetic fillers
    (window.DEMO_PROTEINS || []).forEach(([uniprot, locus, name, product]) => list.push({ id: locus, key: locus, name, locus, product, uniprot }));
    for (let i = 1; list.length < 2300; i++) { const id = `DEMO_${String(i).padStart(4, '0')}`; list.push({ id, key: id, name: id, locus: id, product: SYN_TEMPLATES[i % SYN_TEMPLATES.length](i % 30 + 1) + ' (simulated)', uniprot: '' }); }
  }
  return list;
}
function buildDemo() {
  const r = mulberry32(253); const list = demoList(r); const means = new Map();
  DATA.entities = new Map(); DATA.rows = new Map(); DATA.abund = new Map(); DATA.sites = new Map(); DATA.seqs = new Map();
  for (const p of list) {
    annotate(p);
    const fs = p.processes.map(k => KW_SHIFT[k]).find(Boolean) || FAM_SHIFT[p.family] || { R: 0, stk: 0, stp: 0 };
    const mu = 22 + gauss(r) * 1.5;
    const Ranc = fs.R * (0.6 + 0.8 * r()) + gauss(r) * 0.35 + (r() < 0.2 ? gauss(r) * 1.7 : 0);
    const R = { anc: Ranc, stk: Ranc + fs.stk * (0.6 + 0.8 * r()) + gauss(r) * 0.25, stp: Ranc + fs.stp * (0.6 + 0.8 * r()) + gauss(r) * 0.25 };
    const g = { anc: 0, stk: gauss(r) * 0.25 + (r() < 0.06 ? gauss(r) * 1.6 : 0), stp: gauss(r) * 0.25 + (r() < 0.06 ? gauss(r) * 1.6 : 0) };
    p.missing = [];
    if (p.name === 'Stk1') { g.stk = -6; g.stp = -2.5; p.missing = ['stk_TSB', 'stk_SNpulse', 'stk_SN']; } // transposon in stk1; less Stk1 in stp (polar effect)
    if (p.name === 'Stp1') { g.stp = -6; p.missing = ['stp_TSB', 'stp_SNpulse', 'stp_SN']; }
    const m = {};
    for (const gt of GENOTYPES) { m[`${gt}_TSB`] = mu + g[gt]; m[`${gt}_SN`] = mu + g[gt] + R[gt]; m[`${gt}_SNpulse`] = mu + g[gt] + 0.3 * R[gt] + gauss(r) * 0.15; }
    if (!p.missing.length && r() < 0.04) { const grp = pick(GROUPS, r); m[grp] -= 3.5; p.missing = [grp]; }
    p.lod = p.missing.length > 0; p.se = 0.16 + 0.28 * r();
    means.set(p.id, m); DATA.entities.set(p.id, p);
  }
  for (const c of DEMO_CONTRASTS) {
    const pc = parseContrast(c); const arr = [];
    for (const p of list) {
      const m = means.get(p.id); let diff, s = p.se;
      if (pc.type === 'int') { diff = (m[`${pc.A}_${pc.M1}`] - m[`${pc.B}_${pc.M1}`]) - (m[`${pc.A}_${pc.M2}`] - m[`${pc.B}_${pc.M2}`]); s *= 2; }
      else if (MEDIA.includes(pc.A)) { diff = m[`${pc.ctx}_${pc.A}`] - m[`${pc.ctx}_${pc.B}`]; s *= Math.SQRT2; }
      else { diff = m[`${pc.A}_${pc.ctx}`] - m[`${pc.B}_${pc.ctx}`]; s *= Math.SQRT2; }
      diff += gauss(r) * 0.05;
      arr.push({ id: p.id, diff, p: pTwoSided(diff / s), fdr: 1, lod: p.lod });
    }
    const q = bh(arr.map(x => x.p)); arr.forEach((x, i) => { x.fdr = q[i]; });
    DATA.rows.set(c, arr);
  }
  DATA.demoMeans = means; DATA.mode = 'demo'; DATA.isSite = false; DATA.dataset = 'demo';
  indexRows();
}
function demoAbundance(id) {
  const p = DATA.entities.get(id), m = DATA.demoMeans.get(id); const r = mulberry32(hashStr(id)); const out = [];
  for (const g of GENOTYPES) for (const md of MEDIA) for (let rep = 1; rep <= 3; rep++) {
    const grp = `${g}_${md}`; out.push({ genotype: g, medium: md, replicate: rep, value: p.missing.includes(grp) ? null : m[grp] + gauss(r) * p.se * 1.2 });
  }
  return out;
}

/* ---- 4b. the real course data, read locally from the user's folder ---- */
async function readText(file) {
  if (file.name.endsWith('.gz')) return new Response(file.stream().pipeThrough(new DecompressionStream('gzip'))).text();
  return file.text();
}
function parseTSV(text) {
  const lines = text.split(/\r?\n/); const head = lines[0].split('\t'); const out = [];
  for (let i = 1; i < lines.length; i++) {
    if (!lines[i]) continue; const v = lines[i].split('\t'); const o = {};
    for (let j = 0; j < head.length; j++) o[head[j]] = v[j] === undefined || v[j] === 'NA' ? '' : v[j];
    out.push(o);
  }
  return out;
}
const findFile = (files, re) => [...files].find(f => re.test(f.webkitRelativePath || f.name));
const num = s => (s === '' || s == null) ? NaN : +s;
const isTrue = s => /^(TRUE|true|1)$/.test(s || '');

async function loadCourseFolder(files, dataset) {
  DATA.files = files;
  const status = msg => { $('#loadStatus').textContent = msg; };
  const deaFile = findFile(files, new RegExp(`dea_${dataset}\\.tsv(\\.gz)?$`));
  if (!deaFile) { status(`Could not find dea_${dataset}.tsv.gz — pick the course "data" or "derived" folder.`); return false; }
  status('Reading protein annotation…');
  const ann = new Map();
  const annFile = findFile(files, /protein_annotation\.tsv(\.gz)?$/);
  if (annFile) for (const a of parseTSV(await readText(annFile))) ann.set(a.protein_key, a);
  status(`Reading ${deaFile.name} (all proteins, all contrasts)…`);
  const dea = parseTSV(await readText(deaFile));
  const ents = new Map(), rows = new Map();
  const isSite = 'site_id' in (dea[0] || {}) || 'site' in (dea[0] || {});
  for (const d of dea) {
    if (isTrue(d.is_contaminant)) continue;
    const id = d.site_id || d.site || d.protein_key;
    if (!ents.has(id)) {
      const a = ann.get(d.protein_key) || d;
      const base = a.staa8_gene || a.locus_tag || d.protein_key;
      const site = isSite && d.residue ? { residue: d.residue, position: +d.position } : null;
      const e = { id, key: d.protein_key, locus: a.locus_tag || '', product: a.product || '', uniprot: a.staa8_uniprot || '', name: site ? `${base} ${site.residue}${site.position}` : base, site };
      annotate(e); ents.set(id, e);
    }
    const diff = num(d.diff); if (Number.isNaN(diff)) continue;
    if (!rows.has(d.contrast)) rows.set(d.contrast, []);
    const fdr = num(d.FDR);
    rows.get(d.contrast).push({ id, diff, p: num(d['p.value']), fdr: Number.isNaN(fdr) ? 1 : fdr, lod: d.estimate_type === 'lod_imputed' });
  }
  DATA.entities = ents; DATA.rows = rows; DATA.isSite = isSite; DATA.dataset = dataset; DATA.mode = 'real';
  DATA.abund = new Map(); DATA.sites = new Map(); DATA.seqs = new Map();
  indexRows();
  const abFile = findFile(files, new RegExp(`abundances_${dataset}\\.tsv(\\.gz)?$`));
  if (abFile) {
    status(`Reading ${abFile.name} (per-sample values)…`);
    for (const a of parseTSV(await readText(abFile))) {
      const id = a.site_id || a.site || a.protein_key; if (!ents.has(id)) continue;
      if (!DATA.abund.has(id)) DATA.abund.set(id, []);
      const v = num(a.normalized_abundance);
      DATA.abund.get(id).push({ genotype: a.genotype, medium: a.medium, replicate: a.replicate, value: Number.isNaN(v) ? null : v });
    }
  }
  const addSite = (key, residue, position, prob) => {
    if (!key || !residue) return; const list = DATA.sites.get(key) || []; if (list.some(s => s.position === position && s.residue === residue)) return;
    list.push({ residue, position, prob, label: `${residue}${position}` }); DATA.sites.set(key, list);
  };
  if (isSite) for (const e of ents.values()) if (e.site) addSite(e.key, e.site.residue, e.site.position, NaN);
  const siteFile = findFile(files, /sites_STYH\.tsv(\.gz)?$/);
  if (siteFile) {
    status('Reading phosphosites (sites_STYH)…'); const seenSite = new Set();
    for (const s of parseTSV(await readText(siteFile))) { if (seenSite.has(s.site_id)) continue; seenSite.add(s.site_id); addSite(s.protein_key, s.residue, +s.position, num(s.best_loc_prob)); }
  }
  const fasta = findFile(files, /JE2_proteins\.fa(sta)?$/);
  if (fasta) {
    status('Reading JE2 protein sequences…');
    let key = null, locus = null, seq = [];
    const flush = () => { if (key) { const s = seq.join(''); DATA.seqs.set(key, s); if (locus) DATA.seqs.set(locus, s); } };
    for (const line of (await readText(fasta)).split(/\r?\n/)) {
      if (line.startsWith('>')) { flush(); key = line.slice(1).split(/\s/)[0].replace(/\|/g, '_'); locus = (line.match(/locus_tag=([^\]\s]+)/) || [])[1] || null; seq = []; }
      else seq.push(line.trim());
    }
    flush();
  }
  status(`Loaded ${ents.size.toLocaleString()} ${isSite ? 'phosphosites' : 'proteins'} × ${rows.size} contrasts from your disk (nothing uploaded).`);
  return true;
}

/* ============================ 5. sound effects (Web Audio, synthesised) ============================ */
const SFX = {
  ctx: null,
  ac() {
    if (!this.ctx) { const C = window.AudioContext || window.webkitAudioContext; if (!C) return null; this.ctx = new C(); }
    if (this.ctx.state === 'suspended') this.ctx.resume(); return this.ctx;
  },
  tone(freq, t0, dur, { type = 'triangle', vol = 0.16, slide = 0 } = {}) {
    const ac = this.ac(); if (!ac || GAME.muted) return;
    const o = ac.createOscillator(), g = ac.createGain(), t = ac.currentTime + t0;
    o.type = type; o.frequency.setValueAtTime(freq, t);
    if (slide) o.frequency.exponentialRampToValueAtTime(Math.max(40, freq * slide), t + dur);
    g.gain.setValueAtTime(0.0001, t); g.gain.exponentialRampToValueAtTime(vol, t + 0.012); g.gain.exponentialRampToValueAtTime(0.0001, t + dur);
    o.connect(g).connect(ac.destination); o.start(t); o.stop(t + dur + 0.03);
  },
  noise(t0, dur, { vol = 0.12, from = 400, to = 4000 } = {}) {
    const ac = this.ac(); if (!ac || GAME.muted) return;
    const len = Math.floor(ac.sampleRate * dur), buf = ac.createBuffer(1, len, ac.sampleRate), d = buf.getChannelData(0);
    for (let i = 0; i < len; i++) d[i] = Math.random() * 2 - 1;
    const src = ac.createBufferSource(), f = ac.createBiquadFilter(), g = ac.createGain(), t = ac.currentTime + t0;
    src.buffer = buf; f.type = 'bandpass'; f.Q.value = 1.2; f.frequency.setValueAtTime(from, t); f.frequency.exponentialRampToValueAtTime(to, t + dur);
    g.gain.setValueAtTime(vol, t); g.gain.exponentialRampToValueAtTime(0.0001, t + dur);
    src.connect(f).connect(g).connect(ac.destination); src.start(t); src.stop(t + dur);
  },
  play(name, n = 0) {
    const up = 1 + 0.06 * Math.min(n, 12); // pitch rises with combo
    switch (name) {
      case 'xp': this.tone(880 * up, 0, 0.09); this.tone(1320 * up, 0.07, 0.13); break;
      case 'bigxp': [784, 988, 1319, 1568].forEach((f, i) => this.tone(f, i * 0.07, 0.16)); break;
      case 'catch': this.tone(1047 * up, 0, 0.06, { type: 'square', vol: 0.07 }); this.tone(1568 * up, 0.05, 0.1, { type: 'square', vol: 0.07 }); break;
      case 'levelup': [523, 659, 784, 1047].forEach((f, i) => this.tone(f, i * 0.1, 0.22, { type: 'square', vol: 0.09 })); [1047, 1319, 1568].forEach(f => this.tone(f, 0.42, 0.7, { vol: 0.1 })); this.noise(0.42, 0.6, { vol: 0.05, from: 3000, to: 9000 }); break;
      case 'badge': this.tone(1568, 0, 0.7, { type: 'sine', vol: 0.14 }); this.tone(2093, 0, 0.7, { type: 'sine', vol: 0.07 }); this.tone(2637, 0.12, 0.6, { type: 'sine', vol: 0.06 }); break;
      case 'lose': this.tone(392, 0, 0.2, { type: 'sawtooth', vol: 0.09, slide: 0.75 }); this.tone(294, 0.18, 0.35, { type: 'sawtooth', vol: 0.09, slide: 0.6 }); break;
      case 'leveldown': [440, 349, 262, 196].forEach((f, i) => this.tone(f, i * 0.14, 0.25, { type: 'sawtooth', vol: 0.08 })); break;
      case 'wrong': this.tone(140, 0, 0.28, { type: 'square', vol: 0.1 }); this.tone(105, 0, 0.28, { type: 'square', vol: 0.08 }); break;
      case 'tick': this.tone(1400, 0, 0.035, { type: 'square', vol: 0.05 }); break;
      case 'start': [392, 523, 659, 784].forEach((f, i) => this.tone(f, i * 0.06, 0.12, { type: 'square', vol: 0.07 })); break;
      case 'gameover': [523, 392, 330, 262].forEach((f, i) => this.tone(f, i * 0.18, 0.3)); break;
      case 'whoosh': this.noise(0, 0.45, { vol: 0.08, from: 300, to: 3000 }); break;
      case 'shield': this.tone(660, 0, 0.3, { type: 'sine', vol: 0.12, slide: 2 }); break;
      case 'click': this.tone(700, 0, 0.04, { type: 'square', vol: 0.04 }); break;
    }
  },
};

/* ============================ 6. memes (text + emoji, no images of real people) ============================ */
const MEMES = {
  good: ['📈 STONKS', '🧑‍🍳 LET HIM COOK', '✨ +1000 AURA', '🗿 SIGMA SCIENTIST', '🦛 MOO DENG APPROVES', '😎 just a chill guy', '🖐️ 6 7 🖐️', '🐐 GOAT', '🧢 NO CAP',
    '💅 it\'s giving significance', '🧠 BIG BRAIN', '🚀 GPU GO BRRR', '🔥 ON FIRE', '🍿 PLOT TWIST: IT WAS THE KINASE', '💯 FDR-CERTIFIED'],
  big: ['🧥 THE MORE YOU BUY, THE MORE YOU SAVE 🧥', '🧥 LEATHER-JACKET MODE: ON 🟩', '🟩 ONE MORE GPU FOR THIS VOLCANO 🟩', '👑 PROTEOME LEGEND 👑'],
  bad: ['💀 SKILL ISSUE', '🤡', '😭 EMOTIONAL DAMAGE', '🫠 OHIO MOMENT', '📉 NOT STONKS', '🔥 this is fine 🔥', '🙈 SIGN BACKWARDS'],
};
let lastMeme = 0;
function memeFly(kind = 'good') {
  const now = performance.now(); if (now - lastMeme < 500 && kind !== 'big') return; lastMeme = now;
  const d = document.createElement('div'); d.className = `meme ${kind}`; d.textContent = pick(MEMES[kind]);
  if (kind !== 'big') { d.style.top = `${12 + Math.random() * 60}%`; d.style.setProperty('--rot', `${(Math.random() - 0.5) * 16}deg`); d.style.animationDuration = `${2.2 + Math.random()}s`; }
  $('#memes').appendChild(d); SFX.play(kind === 'bad' ? 'wrong' : 'whoosh');
  setTimeout(() => d.remove(), 3600);
}

/* ============================ 7. game engine: XP, levels, badges ============================ */
const TITLES = ['Pipette Rookie', 'Spectrum Sleuth', 'Volcano Ranger', 'FDR Guardian', 'Phospho Master', 'Proteome Legend'];
const BADGES = {
  first_card: ['🔬', 'First look', 'Open your first protein card'],
  raw_first: ['📊', 'Raw-data respecter', 'Look at per-sample values for 5 proteins'],
  sniper5: ['🎯', 'Sharp shooter', 'Combo ×5 in Sign Sniper'],
  sniper25: ['🏹', 'Sign sensei', '25 correct Sign Sniper answers'],
  sniper500: ['🧥', 'Leather-jacket sniper', '500 points in one Sign Sniper game'],
  int_master: ['🧠', 'Interaction thinker', 'Answer 5 interaction questions correctly'],
  lod3: ['🪤', 'Trap spotter', 'Recognise 3 LOD-imputed estimates'],
  reformed: ['😇', 'Reformed p-hacker', 'Loosen FDR above 0.25, then return to ≤ 0.1'],
  catcher300: ['🧪', 'Beaker boss', '300 points in one Family Catcher game'],
  combo10: ['🔥', 'Catch streak', 'Combo ×10 in Family Catcher'],
  explorer5: ['🗺️', 'Pathfinder', 'Visit 5 pathways'],
  explorer15: ['🧭', 'Cartographer', 'Visit 15 pathways'],
  struct3: ['🧊', 'Structural biologist', 'View 3 AlphaFold structures'],
  dex25: ['📚', 'Collector', '25 entries in your Pokédex'],
  legend: ['👑', 'Legendary pair', 'Collect both Stk1 and Stp1'],
  real: ['💾', 'Real scientist', 'Load the real course data'],
};
const STORE_KEY = 'je2-proteome-quest-v1';
const GAME = { players: {}, current: 'Player 1', muted: false };
function loadStore() { try { const s = JSON.parse(localStorage.getItem(STORE_KEY) || 'null'); if (s) Object.assign(GAME, s); } catch (e) { /* storage unavailable */ } }
function saveStore() { try { localStorage.setItem(STORE_KEY, JSON.stringify(GAME)); } catch (e) { /* ignore */ } }
function me() {
  const p = GAME.players[GAME.current] ||= {};
  const defaults = { xp: 0, badges: [], dex: [], structs: [], visited: [], raw: 0, sniperOk: 0, intOk: 0, lodOk: 0, wasHacker: false, sniperBest: 0, catcherBest: 0 };
  for (const [k, v] of Object.entries(defaults)) if (p[k] === undefined) p[k] = Array.isArray(v) ? [] : v;
  return p;
}
const level = xp => Math.min(TITLES.length, 1 + Math.floor(Math.sqrt(Math.max(0, xp) / 60)));
function toast(html, kind = '') { const t = document.createElement('div'); t.className = `toast ${kind}`; t.innerHTML = html; $('#toasts').appendChild(t); setTimeout(() => t.remove(), 3200); }
function addXP(n, why, { combo = 0, quiet = false } = {}) {
  const p = me(); const before = level(p.xp); p.xp = Math.max(0, p.xp + n); saveStore(); renderHUD(true);
  const after = level(p.xp);
  if (n < 0) { SFX.play('lose'); toast(`−${-n} XP · ${esc(why)}`, 'bad'); }
  else if (!quiet) { SFX.play(n >= 40 ? 'bigxp' : 'xp', combo); toast(`+${n} XP · ${esc(why)}`); }
  else SFX.play('xp', combo);
  if (after > before) { setTimeout(() => SFX.play('levelup'), 250); toast(`⬆️ LEVEL UP! Lv ${after}: <b>${TITLES[after - 1]}</b>`, 'gold'); memeFly('big'); }
  if (after < before) { setTimeout(() => SFX.play('leveldown'), 250); toast(`⬇️ Level down… Lv ${after}: ${TITLES[after - 1]}`, 'bad'); }
}
function award(id) {
  const p = me(); if (p.badges.includes(id)) return; p.badges.push(id); saveStore();
  const [icon, name] = BADGES[id]; setTimeout(() => SFX.play('badge'), 150); toast(`${icon} Badge unlocked: <b>${name}</b>`, 'gold'); addXP(50, 'badge', { quiet: true }); renderBadges();
}
function renderHUD(pulse = false) {
  const p = me(), L = level(p.xp), lo = 60 * (L - 1) ** 2, hi = 60 * L ** 2;
  $('#hudLevel').textContent = `Lv ${L} · ${TITLES[L - 1]}`;
  $('#hudXP').textContent = `${p.xp} XP`;
  $('#hudBar').style.width = `${L >= TITLES.length ? 100 : Math.round(100 * (p.xp - lo) / (hi - lo))}%`;
  $('#hudBadges').textContent = `${p.badges.length}/${Object.keys(BADGES).length} 🏅`;
  if (pulse) { const h = $('#hudXP'); h.classList.remove('pulse'); void h.offsetWidth; h.classList.add('pulse'); }
  $('#muteBtn').textContent = GAME.muted ? '🔇' : '🔊';
}
function renderBadges() {
  const p = me();
  $('#badgeGrid').innerHTML = Object.entries(BADGES).map(([id, [icon, name, how]]) =>
    `<div class="badge ${p.badges.includes(id) ? 'got' : ''}" title="${esc(how)}"><span>${icon}</span><b>${esc(name)}</b><small>${esc(how)}</small></div>`).join('');
  const board = Object.entries(GAME.players).map(([n, q]) => [n, q.xp || 0, (q.badges || []).length, q.sniperBest || 0, q.catcherBest || 0]).sort((a, b) => b[1] - a[1]);
  $('#leaderboard').innerHTML = board.map(([n, xp, b, s, c], i) => `<tr class="${n === GAME.current ? 'me' : ''}"><td>${['🥇', '🥈', '🥉'][i] || i + 1}</td><td>${esc(n)}</td><td>${xp}</td><td>${b}</td><td>${s}</td><td>${c}</td></tr>`).join('');
}

/* ============================ 8. live volcano, visibility switches, false-discovery meter ============================ */
const UI = { contrast: '', fdr: 0.1, lfc: 1, show: { up: true, down: true, ns: true, lod: true }, highlight: new Set(), plotted: false, trapDone: new Set(), px: { cat: null, pw: null, id: null }, tab: 'tSniper' };
const rowsNow = () => DATA.rows.get(UI.contrast) || [];
const isHit = r => r.fdr < UI.fdr && Math.abs(r.diff) > UI.lfc;
const yOf = r => -Math.log10(Math.max(r.fdr, 1e-300));
const catOf = r => r.lod ? 'lod' : isHit(r) ? (r.diff > 0 ? 'up' : 'down') : 'ns';

function renderToggles(counts) {
  const info = contrastInfo(UI.contrast);
  const defs = [['up', '▲', info.pos, '#d55e00'], ['down', '▼', info.neg, '#0072b2'], ['ns', '●', 'not selected', '#9aa3ad'], ['lod', '✕', 'lod_imputed estimates', '#7b5ea7']];
  $('#toggles').innerHTML = defs.map(([k, sym, label, col]) =>
    `<button class="chip ${UI.show[k] ? 'on' : 'off'}" data-k="${k}" role="switch" aria-checked="${UI.show[k]}" title="Click to ${UI.show[k] ? 'hide' : 'show'} these points">
      <span class="sw"><span></span></span><span class="sym" style="color:${col}">${sym}</span> ${esc(label)} <b>${counts[k]}</b></button>`).join('');
  $('#toggles').querySelectorAll('.chip').forEach(b => b.onclick = () => { UI.show[b.dataset.k] = !UI.show[b.dataset.k]; SFX.play('click'); renderVolcano(); });
}
function renderVolcano() {
  const groups = { up: [], down: [], ns: [], lod: [] }; const hl = [];
  for (const r of rowsNow()) { groups[catOf(r)].push(r); if (UI.highlight.has(r.id)) hl.push(r); }
  renderToggles({ up: groups.up.length, down: groups.down.length, ns: groups.ns.length, lod: groups.lod.length });
  renderCounts(groups);
  if (!window.Plotly) { $('#volcano').innerHTML = '<p class="warn">Plotly could not load (offline?). The games still work.</p>'; return; }
  const hover = r => { const e = DATA.entities.get(r.id); return `<b>${esc(e.name)}</b><br>${esc(e.product).slice(0, 60)}<br>diff ${fmt(r.diff)} · FDR ${fmtP(r.fdr)}${r.lod ? '<br>⚠ lod_imputed' : ''}`; };
  const tr = (arr, marker, visible = true) => ({ type: 'scattergl', mode: 'markers', visible, x: arr.map(r => r.diff), y: arr.map(yOf), customdata: arr.map(r => r.id), text: arr.map(hover), hovertemplate: '%{text}<extra></extra>', marker });
  const traces = [
    tr(groups.ns, { size: 5, color: '#c3c9d2' }, UI.show.ns),
    tr(groups.up, { size: 7, color: '#d55e00' }, UI.show.up),
    tr(groups.down, { size: 7, color: '#0072b2' }, UI.show.down),
    tr(groups.lod, { size: 8, symbol: 'x-thin', line: { width: 1.5, color: '#7b5ea7' }, color: '#7b5ea7' }, UI.show.lod),
    tr(hl, { size: 15, color: 'rgba(0,0,0,0)', line: { width: 2.5, color: '#f5b700' } }),
  ];
  const yFdr = -Math.log10(UI.fdr);
  const layout = {
    margin: { l: 55, r: 10, t: 10, b: 45 }, hovermode: 'closest', dragmode: 'pan', uirevision: UI.contrast, showlegend: false,
    xaxis: { title: 'log2 effect (diff)', zeroline: true, zerolinecolor: '#999' }, yaxis: { title: '−log10(FDR)' },
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: '#fff',
    shapes: [
      { type: 'line', xref: 'x', yref: 'paper', x0: -UI.lfc, x1: -UI.lfc, y0: 0, y1: 1, line: { color: '#e69f00', width: 2, dash: 'dash' } },
      { type: 'line', xref: 'x', yref: 'paper', x0: UI.lfc, x1: UI.lfc, y0: 0, y1: 1, line: { color: '#e69f00', width: 2, dash: 'dash' } },
      { type: 'line', xref: 'paper', yref: 'y', x0: 0, x1: 1, y0: yFdr, y1: yFdr, line: { color: '#009e73', width: 2, dash: 'dash' } }],
  };
  const config = { responsive: true, displaylogo: false, edits: { shapePosition: true }, modeBarButtonsToRemove: ['lasso2d', 'select2d'], toImageButtonOptions: { format: 'svg', filename: `volcano_${UI.contrast}` } };
  Plotly.react('volcano', traces, layout, config);
  if (!UI.plotted) {
    UI.plotted = true; const v = $('#volcano');
    v.on('plotly_click', ev => onPointClick(ev.points[0].customdata));
    v.on('plotly_relayout', ev => { // dragging a cutoff line updates the sliders
      let changed = false;
      for (const k of Object.keys(ev)) {
        const m = k.match(/^shapes\[(\d)\]\.([xy])0$/); if (!m) continue;
        if (m[1] !== '2' && m[2] === 'x') { UI.lfc = Math.min(6, Math.abs(+ev[k])); changed = true; }
        if (m[1] === '2' && m[2] === 'y') { UI.fdr = Math.min(0.5, Math.max(1e-6, 10 ** -(+ev[k]))); changed = true; }
      }
      if (changed) { syncSliders(); renderVolcano(); }
    });
  }
}
function syncSliders() {
  $('#fdrSlider').value = UI.fdr; $('#fdrVal').textContent = UI.fdr < 0.001 ? UI.fdr.toExponential(1) : UI.fdr.toFixed(3);
  $('#lfcSlider').value = UI.lfc; $('#lfcVal').textContent = UI.lfc.toFixed(2);
}
function renderCounts(groups) {
  const info = contrastInfo(UI.contrast); const hits = groups.up.length + groups.down.length;
  $('#counts').innerHTML = `<div class="count up">▲ ${groups.up.length}<small>${esc(info.pos)}</small></div><div class="count down">▼ ${groups.down.length}<small>${esc(info.neg)}</small></div>`;
  $('#formula').innerHTML = `<code>${info.formula}</code><p>${esc(info.english)}</p>`;
  // false-discovery meter: with BH at level q, expected false discoveries among H hits ≤ q·H
  const exp = UI.fdr * hits; const pct = Math.min(100, 100 * UI.fdr / 0.5); const p = me(); let verdict;
  if (UI.fdr > 0.25 || UI.lfc < 0.3) { verdict = '🚨 p-hacking zone: many of these hits may be noise.'; if (!p.wasHacker) { p.wasHacker = true; SFX.play('wrong'); } }
  else if (UI.fdr <= 0.1 && UI.lfc >= 1) { verdict = '😇 Responsible cutoffs (course default: FDR < 0.1, |diff| > 1).'; if (p.wasHacker) { p.wasHacker = false; award('reformed'); } }
  else verdict = '🤔 Exploratory cutoffs — say so when you report them.';
  saveStore();
  $('#meter').innerHTML = `<div class="meterbar"><div style="width:${pct}%"></div></div>
    <p><b>Expected false hits ≤ ${exp.toFixed(1)}</b> of ${hits} &nbsp;<span class="muted">(FDR cutoff × hits = ${UI.fdr.toFixed(3)} × ${hits}; LOD-imputed estimates not counted)</span></p><p>${verdict}</p>`;
}
function onPointClick(id) {
  const r = DATA.index.get(UI.contrast)?.get(id);
  if (r && r.lod && !UI.trapDone.has(id)) return lodTrap(id);
  openCard(id);
}

/* ============================ 9. games ============================ */
/* 9a. LOD Trap — clicking an imputed estimate triggers a question before the card opens */
function lodTrap(id) {
  const e = DATA.entities.get(id);
  const opts = [['A precise fold change — I can trust the number', false],
    ['Presence/absence evidence: at least one group mean was imputed from the detection limit. Check the per-sample values!', true],
    ['It means the protein is a contaminant', false]].sort(() => Math.random() - 0.5);
  showModal(`<h3>🪤 LOD Trap!</h3><p><b>${esc(e.name)}</b> is flagged <code>estimate_type = lod_imputed</code>. What does that tell you?</p>
    ${opts.map(([t, ok]) => `<button class="opt" data-ok="${ok}">${esc(t)}</button>`).join('')}`, box => {
    box.querySelectorAll('.opt').forEach(b => b.onclick = () => {
      UI.trapDone.add(id);
      if (b.dataset.ok === 'true') { const p = me(); p.lodOk++; addXP(15, 'trap spotted'); if (p.lodOk >= 3) award('lod3'); b.classList.add('ok'); memeFly('good'); }
      else { b.classList.add('bad'); addXP(-5, 'fell into the LOD trap'); memeFly('bad'); }
      setTimeout(() => { closeModal(); openCard(id); }, 1000);
    });
  });
}

/* 9b. Sign Sniper — arcade: 3 lives, 8 s per shot, combos, memes */
const SN = { running: false, T: 8000 };
function sniperIntro(extra = '') {
  clearTimeout(SN.timer); clearInterval(SN.tick); SN.running = false;
  $('#sniper').innerHTML = `<div class="arcade-intro">${extra}<div class="bigicon">🎯</div><h4>Sign Sniper</h4>
    <p>A protein, a contrast, a fold change. Shoot the right side before the timer runs out!</p>
    <p>❤️❤️❤️ 3 lives · ⏱ 8 s per shot · combos multiply points · 🧠 interaction questions (×2) unlock at combo ×4</p>
    <p>Your best: <b>${me().sniperBest}</b></p><button id="snStart" class="big">▶ PLAY</button><p class="muted">Keyboard: ← left target · → right target</p></div>`;
  $('#snStart').onclick = sniperStart;
}
function sniperStart() { Object.assign(SN, { running: true, score: 0, lives: 3, combo: 0, correct: 0 }); SFX.play('start'); sniperNext(); }
function sniperQuestion() {
  const ints = DATA.contrasts.filter(c => parseContrast(c).type === 'int'); const simples = DATA.contrasts.filter(c => parseContrast(c).type === 'simple');
  const isInt = ints.length > 0 && SN.combo >= 4 && Math.random() < 0.4; const pool = isInt ? ints : simples;
  for (let t = 0; t < 40 && pool.length; t++) {
    const c = pick(pool); const cand = (DATA.rows.get(c) || []).filter(x => !x.lod && x.fdr < 0.1 && Math.abs(x.diff) > 1);
    if (cand.length) return { c, r: pick(cand), info: contrastInfo(c), isInt };
  }
  return null;
}
function sniperNext() {
  if (!SN.running) return; const q = sniperQuestion();
  if (!q) { sniperIntro('<p class="bad">No clear hits to quiz on in this dataset.</p>'); return; }
  SN.q = q; SN.t0 = performance.now(); SN.answered = false; const e = DATA.entities.get(q.r.id);
  const hearts = '❤️'.repeat(SN.lives) + '🤍'.repeat(Math.max(0, 3 - SN.lives));
  $('#sniper').innerHTML = `<div class="arcade-hud"><span>${hearts}</span><span>⭐ ${SN.score}</span><span class="combo ${SN.combo >= 2 ? 'hot' : ''}">combo ×${Math.max(1, Math.min(SN.combo, 5))}</span></div>
    <div class="timer"><div style="animation-duration:${SN.T}ms"></div></div>
    <p class="q">${q.isInt ? '<span class="tag hard">🧠 INTERACTION ×2</span> ' : ''}<b>${esc(e.name)}</b> <span class="muted">${esc(e.product).slice(0, 70)}</span></p>
    <p class="contrastline"><code>${esc(q.c)}</code> &nbsp; diff = <b class="bigdiff">${fmt(q.r.diff)}</b></p>
    <div class="targets"><button class="target" data-side="pos">◀ ${esc(q.info.pos)}</button><button class="target" data-side="neg">${esc(q.info.neg)} ▶</button></div>
    <div id="snFeedback"></div>`;
  $('#sniper').querySelectorAll('.target').forEach(b => b.onclick = () => sniperAnswer(b.dataset.side));
  SN.timer = setTimeout(() => sniperAnswer(null), SN.T);
  clearInterval(SN.tick); SN.tick = setInterval(() => { if (SN.T - (performance.now() - SN.t0) < 3200 && !SN.answered) SFX.play('tick'); }, 1000);
}
function sniperAnswer(side) {
  if (!SN.running || SN.answered) return; SN.answered = true; clearTimeout(SN.timer); clearInterval(SN.tick);
  const q = SN.q; const correctSide = q.r.diff > 0 ? 'pos' : 'neg'; const ok = side === correctSide; const p = me();
  $('#sniper').querySelectorAll('.target').forEach(b => { b.disabled = true; b.classList.add(b.dataset.side === correctSide ? 'ok' : 'bad'); });
  $('#sniper .timer div').style.animationPlayState = 'paused';
  let fb;
  if (ok) {
    SN.combo++; SN.correct++; p.sniperOk++; if (q.isInt) p.intOk++;
    const left = Math.max(0, SN.T - (performance.now() - SN.t0));
    const pts = Math.round((10 + left / 400) * Math.min(SN.combo, 5) * (q.isInt ? 2 : 1)); SN.score += pts;
    addXP(Math.min(30, 4 + 2 * SN.combo) * (q.isInt ? 2 : 1), 'Sign Sniper', { combo: SN.combo, quiet: true });
    if (SN.combo % 5 === 0) memeFly('big'); else if (SN.combo >= 2 && Math.random() < 0.7) memeFly('good');
    if (SN.combo >= 5) award('sniper5'); if (p.sniperOk >= 25) award('sniper25'); if (p.intOk >= 5) award('int_master');
    fb = `<div class="fb ok">✅ +${pts} pts${left > SN.T / 2 ? ' ⚡ quick-draw bonus' : ''}</div>`;
  } else {
    SN.lives--; SN.combo = 0; addXP(-5, side === null ? 'too slow' : 'wrong side'); memeFly('bad');
    fb = `<div class="fb bad">${side === null ? '⏰ Too slow!' : '❌ Missed!'} <code>${q.info.formula}</code> → ${fmt(q.r.diff)} is ${q.r.diff > 0 ? 'positive' : 'negative'}: <b>${esc(q.r.diff > 0 ? q.info.pos : q.info.neg)}</b></div>`;
  }
  saveStore(); $('#snFeedback').innerHTML = fb;
  setTimeout(() => (SN.lives <= 0 ? sniperOver() : sniperNext()), ok ? 850 : 2300);
}
function sniperOver() {
  SN.running = false; SFX.play('gameover'); const p = me(); let extra = '';
  const bonus = Math.round(SN.score / 10); if (bonus) addXP(bonus, 'game bonus');
  if (SN.score > p.sniperBest) { p.sniperBest = SN.score; saveStore(); extra = '<p class="record">🏆 NEW HIGH SCORE!</p>'; setTimeout(() => memeFly('big'), 600); }
  if (SN.score >= 500) award('sniper500');
  renderBadges();
  sniperIntro(`<div class="gameover"><h3>GAME OVER</h3><p>⭐ ${SN.score} points · ${SN.correct} correct</p>${extra}</div>`);
}

/* 9c. Family Catcher — falling-bubble arcade on a canvas */
const CG = { running: false };
function catcherIntro(extra = '') {
  const groups = [...CATCH.entries()].filter(([, ids]) => ids.length <= 260).sort((a, b) => a[0].localeCompare(b[0]));
  const sel = CG.group && CATCH.has(CG.group) ? CG.group : (CATCH.has('Ribosomal protein') ? 'Ribosomal protein' : groups[0]?.[0]);
  $('#catcher').innerHTML = `<div class="arcade-intro">${extra}<div class="bigicon">🧪</div><h4>Family Catcher</h4>
    <p>Proteins rain down. Catch members of your family with the beaker, dodge the impostors!</p>
    <div class="row" style="justify-content:center"><select id="cgGroup">${groups.map(([g, ids]) => `<option value="${esc(g)}" ${g === sel ? 'selected' : ''}>${esc(g)} (${ids.length})</option>`).join('')}</select></div>
    <p class="muted" id="cgExamples"></p>
    <p>❤️❤️❤️ 3 lives · ⭐ gold ring = significant in the current contrast (bonus) · 🛡 shield blocks one mistake · speed keeps rising</p>
    <p>Your best: <b>${me().catcherBest}</b></p><button id="cgStart" class="big">▶ PLAY</button><p class="muted">Move with the mouse, finger or ← →</p></div>`;
  const showEx = () => { const ids = CATCH.get($('#cgGroup').value) || []; $('#cgExamples').innerHTML = 'Members look like: ' + ids.slice(0, 7).map(id => `<b>${esc(DATA.entities.get(id).name)}</b>`).join(', ') + '…'; };
  if (groups.length) { $('#cgGroup').onchange = showEx; showEx(); $('#cgStart').onclick = () => catcherStart($('#cgGroup').value); }
  else $('#cgStart').disabled = true;
}
function catcherStart(group) {
  const members = new Set(CATCH.get(group) || []); if (!members.size) return;
  const idx = DATA.index.get(UI.contrast);
  const decoys = [...DATA.entities.keys()].filter(id => !members.has(id));
  $('#catcher').innerHTML = `<canvas id="cgCanvas" tabindex="0"></canvas><div class="row"><button id="cgStop" class="ghost">■ Stop</button><span class="muted">Family: <b>${esc(group)}</b> · contrast <code>${esc(UI.contrast)}</code></span></div>`;
  const cv = $('#cgCanvas'), dpr = window.devicePixelRatio || 1; const W = cv.clientWidth || 460, H = 400;
  cv.width = W * dpr; cv.height = H * dpr; cv.style.height = `${H}px`; const ctx = cv.getContext('2d'); ctx.scale(dpr, dpr);
  Object.assign(CG, { running: true, group, members, decoys, idx, ctx, W, H, px: W / 2, pw: 96, bubbles: [], floaters: [], score: 0, lives: 3, combo: 0, best: 0, shield: false, caught: new Map(), t0: performance.now(), last: performance.now(), nextSpawn: 0 });
  const move = x => { CG.px = Math.max(CG.pw / 2, Math.min(W - CG.pw / 2, x)); };
  cv.onmousemove = ev => move(ev.offsetX);
  cv.ontouchmove = ev => { const r = cv.getBoundingClientRect(); move(ev.touches[0].clientX - r.left); ev.preventDefault(); };
  $('#cgStop').onclick = () => catcherOver(true);
  cv.focus(); SFX.play('start'); requestAnimationFrame(catcherLoop);
}
function catcherSpawn() {
  const isShield = Math.random() < 0.04; const member = !isShield && Math.random() < 0.45;
  const id = isShield ? null : member ? pick([...CG.members]) : pick(CG.decoys); const e = id ? DATA.entities.get(id) : null;
  const r = id ? CG.idx?.get(id) : null;
  CG.bubbles.push({ id, member, shield: isShield, x: 36 + Math.random() * (CG.W - 72), y: -34, r: 31,
    label: e ? e.name.slice(0, 11) : '🛡', sub: e ? (e.product || '').slice(0, 18) : 'SHIELD', diff: r?.diff, sig: r ? isHit(r) : false, wob: Math.random() * 6 });
}
function catcherLoop(ts) {
  if (!CG.running) return;
  const dt = Math.min(0.05, (ts - CG.last) / 1000); CG.last = ts; const el = (ts - CG.t0) / 1000;
  const speed = Math.min(270, 70 + el * 3.2), interval = Math.max(360, 1050 - el * 14);
  if (ts > CG.nextSpawn) { catcherSpawn(); CG.nextSpawn = ts + interval; }
  const { ctx, W, H } = CG; const by = H - 34;
  for (const b of CG.bubbles) { b.y += speed * dt; b.x += Math.sin(ts / 400 + b.wob) * 0.4; }
  for (const b of [...CG.bubbles]) {
    const inX = Math.abs(b.x - CG.px) < CG.pw / 2 + b.r * 0.5;
    if (b.y + b.r >= by - 6 && b.y - b.r <= by + 10 && inX) { CG.bubbles.splice(CG.bubbles.indexOf(b), 1); catcherCatch(b); continue; }
    if (b.y - b.r > H) { CG.bubbles.splice(CG.bubbles.indexOf(b), 1); if (b.member) { if (CG.combo >= 3) floater(b.x, H - 60, `missed ${b.label}`, '#ff8a80'); CG.combo = 0; } }
  }
  // draw
  const g = ctx.createLinearGradient(0, 0, 0, H); g.addColorStop(0, '#0f2233'); g.addColorStop(1, '#173b4f'); ctx.fillStyle = g; ctx.fillRect(0, 0, W, H);
  for (const b of CG.bubbles) {
    ctx.beginPath(); ctx.arc(b.x, b.y, b.r, 0, Math.PI * 2);
    ctx.fillStyle = b.shield ? '#4fc3f7' : b.diff == null ? '#7f8c99' : b.diff > 0 ? 'rgba(213,94,0,.88)' : 'rgba(0,114,178,.88)'; ctx.fill();
    if (b.sig) { ctx.lineWidth = 4; ctx.strokeStyle = '#f5b700'; ctx.stroke(); }
    ctx.fillStyle = '#fff'; ctx.textAlign = 'center'; ctx.font = 'bold 12px system-ui'; ctx.fillText(b.label, b.x, b.y + 1);
    ctx.font = '8px system-ui'; ctx.fillStyle = 'rgba(255,255,255,.85)'; ctx.fillText(b.sub, b.x, b.y + 13);
  }
  // beaker
  ctx.fillStyle = CG.shield ? 'rgba(79,195,247,.35)' : 'rgba(255,255,255,.12)'; ctx.strokeStyle = '#e8f1f8'; ctx.lineWidth = 3;
  ctx.beginPath(); ctx.moveTo(CG.px - CG.pw / 2, by - 10); ctx.lineTo(CG.px - CG.pw / 2 + 8, by + 22); ctx.lineTo(CG.px + CG.pw / 2 - 8, by + 22); ctx.lineTo(CG.px + CG.pw / 2, by - 10); ctx.fill(); ctx.stroke();
  ctx.fillStyle = '#29d391'; ctx.fillRect(CG.px - CG.pw / 2 + 10, by + 8, CG.pw - 20, 12);
  // floaters + HUD
  ctx.textAlign = 'center';
  for (const f of [...CG.floaters]) { f.t += dt; f.y -= 40 * dt; ctx.globalAlpha = Math.max(0, 1 - f.t / 1.2); ctx.fillStyle = f.c; ctx.font = 'bold 14px system-ui'; ctx.fillText(f.s, f.x, f.y); ctx.globalAlpha = 1; if (f.t > 1.2) CG.floaters.splice(CG.floaters.indexOf(f), 1); }
  ctx.textAlign = 'left'; ctx.font = 'bold 15px system-ui'; ctx.fillStyle = '#fff';
  ctx.fillText(`${'❤️'.repeat(CG.lives)}${'🤍'.repeat(Math.max(0, 3 - CG.lives))}  ⭐ ${CG.score}  combo ×${CG.combo}${CG.shield ? '  🛡' : ''}`, 10, 22);
  ctx.textAlign = 'right'; ctx.font = '12px system-ui'; ctx.fillStyle = '#9fd8ff'; ctx.fillText(`catch: ${CG.group}`, W - 10, 22);
  if (CG.lives <= 0) return catcherOver();
  requestAnimationFrame(catcherLoop);
}
function floater(x, y, s, c) { CG.floaters.push({ x, y, s, c, t: 0 }); }
function catcherCatch(b) {
  if (b.shield) { CG.shield = true; SFX.play('shield'); floater(b.x, b.y, '🛡 SHIELD!', '#4fc3f7'); return; }
  if (b.member) {
    CG.combo++; CG.best = Math.max(CG.best, CG.combo); const mult = 1 + Math.floor(CG.combo / 5); const pts = (10 + (b.sig ? 10 : 0)) * mult; CG.score += pts;
    CG.caught.set(b.id, (CG.caught.get(b.id) || 0) + 1); SFX.play('catch', CG.combo);
    floater(b.x, b.y, `+${pts}${b.sig ? ' ⭐' : ''}`, b.sig ? '#f5b700' : '#29d391');
    if (CG.combo % 5 === 0) { memeFly(CG.combo % 15 === 0 ? 'big' : 'good'); addXP(5 * mult, `catch streak ×${CG.combo}`, { combo: CG.combo, quiet: true }); }
    if (CG.combo >= 10) award('combo10');
  } else if (CG.shield) { CG.shield = false; SFX.play('shield'); floater(b.x, b.y, `🛡 blocked ${b.label}`, '#4fc3f7'); }
  else {
    CG.lives--; CG.combo = 0; addXP(-3, `${b.label} is not ${CG.group}`); floater(b.x, b.y, `✗ ${b.label}: impostor!`, '#ff6b6b');
    if (CG.lives > 0) memeFly('bad'); else setTimeout(() => catcherOver(), 0);
  }
}
function catcherOver(stopped = false) {
  if (!CG.running) return; CG.running = false; SFX.play('gameover'); const p = me(); let extra = '';
  const xp = Math.round(CG.score / 5); if (xp) addXP(xp, 'Family Catcher');
  if (CG.score > p.catcherBest) { p.catcherBest = CG.score; saveStore(); if (CG.score > 0) { extra = '<p class="record">🏆 NEW HIGH SCORE!</p>'; setTimeout(() => memeFly('big'), 500); } }
  if (CG.score >= 300) award('catcher300'); renderBadges();
  const caught = [...CG.caught.keys()].map(id => DATA.entities.get(id)).filter(Boolean);
  catcherIntro(`<div class="gameover"><h3>${stopped ? 'STOPPED' : 'GAME OVER'}</h3><p>⭐ ${CG.score} points · best combo ×${CG.best} · ${caught.length} different ${esc(CG.group)} proteins caught</p>${extra}
    ${caught.length ? `<p class="muted">Learn your catch (click to open):</p><div class="chips">${caught.map(e => `<button class="mini" data-id="${esc(e.id)}">${esc(e.name)}</button>`).join('')}</div>` : ''}</div>`);
  $('#catcher').querySelectorAll('.mini').forEach(b => b.onclick = () => openCard(b.dataset.id));
}

/* 9d. Pathway Explorer — category → pathway → protein → AlphaFold 3D */
function computeGroups() {
  EXPLORER = new Map(); CATCH = new Map();
  const add = (cat, pw, id, sub = '') => {
    if (!cat || !pw) return; if (!EXPLORER.has(cat)) EXPLORER.set(cat, new Map()); const m = EXPLORER.get(cat);
    if (!m.has(pw)) m.set(pw, { ids: new Set(), sub: new Map() }); const o = m.get(pw); o.ids.add(id); if (sub) o.sub.set(id, sub);
  };
  for (const e of DATA.entities.values()) {
    for (const p of e.pathways || []) add(PW_CAT, p.top, e.id, p.sub);
    for (const k of e.processes || []) add(KW_TO_CAT.get(k), k, e.id);
    if (e.family) add(FAM_CAT, e.family, e.id);
  }
  for (const [cat, m] of EXPLORER) { for (const [pw, o] of m) if (o.ids.size < 2) m.delete(pw); if (!m.size) EXPLORER.delete(cat); }
  for (const [cat, m] of EXPLORER) for (const [pw, o] of m) if (o.ids.size >= 8 && cat !== PW_CAT && !CATCH.has(pw)) CATCH.set(pw, [...o.ids]);
}
function pwStats(ids) {
  const idx = DATA.index.get(UI.contrast); let up = 0, down = 0;
  for (const id of ids) { const r = idx?.get(id); if (r && !r.lod && isHit(r)) { if (r.diff > 0) up++; else down++; } }
  return { up, down };
}
function renderExplorer() {
  const P = UI.px; const box = $('#explorer');
  const crumbs = ['<a data-l="0">🗺️ All</a>']; if (P.cat) crumbs.push(`<a data-l="1">${esc(P.cat)}</a>`); if (P.pw) crumbs.push(`<a data-l="2">${esc(P.pw)}</a>`); if (P.id) crumbs.push(`<span>${esc(DATA.entities.get(P.id)?.name)}</span>`);
  let html = `<div class="crumbs">${crumbs.join(' › ')}</div><p class="muted">▲/▼ = hits in the current contrast <code>${esc(UI.contrast)}</code> with your cutoffs.</p>`;
  if (!P.cat) {
    html += `<div class="pxgrid">${[...EXPLORER.entries()].map(([cat, m]) => { const ids = new Set([...m.values()].flatMap(o => [...o.ids])); const s = pwStats(ids);
      return `<button class="pxcard" data-cat="${esc(cat)}"><b>${esc(cat)}</b><small>${m.size} pathways · ${ids.size} proteins · <span class="u">▲${s.up}</span> <span class="d">▼${s.down}</span></small></button>`; }).join('')}</div>`;
    if (!EXPLORER.size) html += '<p class="warn">No pathway annotation available (UniProt unreachable?).</p>';
  } else if (!P.pw) {
    const list = [...EXPLORER.get(P.cat).entries()].map(([pw, o]) => ({ pw, n: o.ids.size, ...pwStats(o.ids) })).sort((a, b) => (b.up + b.down) - (a.up + a.down) || b.n - a.n);
    const visited = new Set(me().visited);
    html += list.map(x => `<button class="pxrow" data-pw="${esc(x.pw)}"><span>${visited.has(x.pw) ? '✅' : '🆕'} <b>${esc(x.pw)}</b> <small class="muted">${x.n} proteins</small></span>
      <span class="minibar"><i class="u" style="width:${100 * x.up / x.n}%"></i><i class="d" style="width:${100 * x.down / x.n}%"></i></span><small><span class="u">▲${x.up}</span> <span class="d">▼${x.down}</span></small></button>`).join('');
  } else if (!P.id) {
    const o = EXPLORER.get(P.cat).get(P.pw); const idx = DATA.index.get(UI.contrast);
    const list = [...o.ids].map(id => ({ e: DATA.entities.get(id), r: idx?.get(id), sub: o.sub.get(id) })).sort((a, b) => (a.r?.fdr ?? 2) - (b.r?.fdr ?? 2));
    html += `<p><button id="pxHL" class="ghost">✨ Show all ${list.length} on the volcano</button> <span class="muted">Pick a protein to see its 3D structure.</span></p>` + list.map(({ e, r, sub }) => {
      const cls = r ? (r.lod ? 'lod' : isHit(r) ? (r.diff > 0 ? 'u' : 'd') : 'n') : 'n';
      return `<button class="pxrow" data-id="${esc(e.id)}"><span><b>${esc(e.name)}</b> <small class="muted">${esc(e.product).slice(0, 48)}${sub ? ` · ${esc(sub).slice(0, 40)}` : ''}</small></span>
        <span class="diffchip ${cls}">${r ? fmt(r.diff) : 'NA'}${r?.lod ? ' ✕' : ''}</span><small>FDR ${r ? fmtP(r.fdr) : 'NA'}</small></button>`; }).join('');
  } else {
    const e = DATA.entities.get(P.id), r = DATA.index.get(UI.contrast)?.get(P.id), info = contrastInfo(UI.contrast);
    html += `<div class="pxprotein"><h3>${esc(e.name)} ${isLegend(e) ? '👑' : ''}</h3><p>${esc(e.product)}</p>
      <p>${r ? `In <code>${esc(UI.contrast)}</code>: diff <b>${fmt(r.diff)}</b>, FDR ${fmtP(r.fdr)} → <b>${esc(r.diff > 0 ? info.pos : info.neg)}</b>${r.lod ? ' <span class="tag warn">lod_imputed</span>' : ''}` : 'Not measured in this contrast.'}</p>
      <div id="pxViewer" class="viewer"></div><p id="pxMsg" class="muted"></p>
      <div class="row"><button id="pxCard">Open full protein card</button><button id="pxBack" class="ghost">◀ Back to ${esc(P.pw)}</button></div></div>`;
  }
  box.innerHTML = html;
  box.querySelectorAll('.crumbs a').forEach(a => a.onclick = () => { const l = +a.dataset.l; if (l < 1) P.cat = null; if (l < 2) P.pw = null; P.id = null; SFX.play('click'); renderExplorer(); });
  box.querySelectorAll('[data-cat]').forEach(b => b.onclick = () => { P.cat = b.dataset.cat; SFX.play('click'); renderExplorer(); });
  box.querySelectorAll('[data-pw]').forEach(b => b.onclick = () => {
    P.pw = b.dataset.pw; const p = me(); SFX.play('click');
    if (!p.visited.includes(P.pw)) { p.visited.push(P.pw); saveStore(); addXP(10, `discovered ${P.pw}`); if (p.visited.length >= 5) award('explorer5'); if (p.visited.length >= 15) award('explorer15'); }
    UI.highlight = new Set(EXPLORER.get(P.cat).get(P.pw).ids); renderVolcano(); renderExplorer();
  });
  box.querySelectorAll('.pxrow[data-id]').forEach(b => b.onclick = () => { P.id = b.dataset.id; SFX.play('click'); collect(P.id); renderExplorer(); });
  if ($('#pxHL')) $('#pxHL').onclick = () => { UI.highlight = new Set(EXPLORER.get(P.cat).get(P.pw).ids); renderVolcano(); $('#volcano').scrollIntoView({ behavior: 'smooth' }); };
  if (P.id) {
    const e = DATA.entities.get(P.id); UI.highlight = new Set([P.id]); renderVolcano();
    $('#pxCard').onclick = () => openCard(P.id); $('#pxBack').onclick = () => { P.id = null; renderExplorer(); };
    if (e.uniprot) showStructure(e, $('#pxViewer'), $('#pxMsg')); else $('#pxMsg').textContent = 'No UniProt ortholog mapped — no AlphaFold model available.';
  }
}

/* 9e. Pokédex + protein card */
const isLegend = e => !!e && (/^(stk1|stp1|pknB)$/i.test(e.name) || /^B7H15_0635[05]$/.test(e.locus) || /^Q2FZ6[45]$/.test(e.uniprot));
function collect(id) {
  const p = me(); if (p.dex.includes(id)) return; p.dex.push(id); saveStore(); addXP(5, 'new Pokédex entry');
  if (p.dex.length === 1) award('first_card'); if (p.dex.length >= 25) award('dex25');
  const leg = p.dex.map(i => DATA.entities.get(i)).filter(isLegend).map(e => e.uniprot || e.locus || e.name.toLowerCase());
  if (new Set(leg).size >= 2) award('legend');
  renderDex();
}
function renderDex() {
  const p = me(); const ents = p.dex.map(id => DATA.entities.get(id)).filter(Boolean);
  $('#dexCount').textContent = `${ents.length} of ${DATA.entities.size.toLocaleString()} collected in this data set · 🧊 = structure viewed · 👑 = legendary`;
  $('#dexGrid').innerHTML = ents.map(e => `<button class="dex ${p.structs.includes(e.id) ? 'star' : ''} ${isLegend(e) ? 'legend' : ''}" data-id="${esc(e.id)}">
    <b>${esc(e.name)}</b><small>${esc(e.processes?.[0] || e.family || e.product).slice(0, 34)}</small></button>`).join('') || '<p class="muted">Open protein cards (volcano, browser, explorer) to fill your Pokédex.</p>';
  $('#dexGrid').querySelectorAll('.dex').forEach(b => b.onclick = () => openCard(b.dataset.id));
}
function openCard(id) {
  const e = DATA.entities.get(id); if (!e) return; collect(id);
  const rows = DATA.contrasts.map(c => [c, DATA.index.get(c).get(id)]).filter(x => x[1]);
  const up = e.uniprot; const tags = [...(e.pathways || []).map(p => p.top), ...(e.processes || [])].slice(0, 8);
  $('#card').innerHTML = `<div class="cardhead"><div><h2>${esc(e.name)} ${isLegend(e) ? '👑' : ''}</h2>
      <p>${esc(e.product)}</p><p class="muted">${esc(e.locus || e.key)}${up ? ` · UniProt ${esc(up)} (NCTC8325 reference)` : ''}${DATA.mode === 'demo' ? ' · <span class="tag warn">simulated values</span>' : ''}</p>
      <p>${tags.map(t => `<span class="tag">${esc(t)}</span>`).join(' ')}</p></div><button class="ghost" id="cardClose">✕</button></div>
    <div class="links">${up ? `<a target="_blank" rel="noopener" href="https://www.uniprot.org/uniprotkb/${encodeURIComponent(up)}/entry">UniProt</a>
      <a target="_blank" rel="noopener" href="https://alphafold.ebi.ac.uk/entry/${encodeURIComponent(up)}">AlphaFold DB</a>
      <a target="_blank" rel="noopener" href="https://string-db.org/cgi/network?identifiers=${encodeURIComponent(up)}&species=93061">STRING</a>` : '<span class="muted">No UniProt ortholog mapped.</span>'}</div>
    <div class="cardgrid"><div><h4>Per-sample values <span class="muted">(look before trusting a fold change)</span></h4><div id="cardPlot"></div></div>
      <div><h4>All ${rows.length} contrasts</h4><div class="tablewrap"><table><tr><th>contrast</th><th>diff</th><th>FDR</th><th>meaning</th></tr>
      ${rows.map(([c, r]) => { const i = contrastInfo(c); const hit = r.fdr < 0.1 && Math.abs(r.diff) > 1; return `<tr class="${hit ? 'hit' : ''}"><td><code>${esc(c)}</code></td><td>${fmt(r.diff)}${r.lod ? ' <span title="lod_imputed">⚠</span>' : ''}</td><td>${fmtP(r.fdr)}</td><td>${esc(r.diff > 0 ? i.pos : i.neg)}</td></tr>`; }).join('')}
      </table></div></div></div>
    <h4>Sequence & AlphaFold structure</h4>
    <div class="row"><button id="structBtn" ${up ? '' : 'disabled'}>🧊 Load sequence + 3D model</button>
      <input id="motif" placeholder="motif regex, e.g. R..S or [ST]P"><span id="structMsg" class="muted"></span></div>
    <div class="seqgrid"><div id="seq"></div><div id="viewer" class="viewer"></div></div>`;
  $('#cardClose').onclick = () => { $('#card').innerHTML = ''; };
  $('#structBtn').onclick = () => showStructure(e, $('#viewer'), $('#structMsg'), true);
  $('#motif').oninput = () => renderSeq(e);
  if (DATA.seqs.has(e.key) || DATA.seqs.has(e.locus)) renderSeq(e);
  plotSamples(e);
  $('#card').scrollIntoView({ behavior: 'smooth', block: 'start' });
}
function plotSamples(e) {
  const vals = DATA.mode === 'demo' ? demoAbundance(e.id) : (DATA.abund.get(e.id) || []);
  if (!vals.length) { $('#cardPlot').innerHTML = '<p class="muted">Per-sample values need abundances_*.tsv.gz in the loaded folder.</p>'; return; }
  const p = me(); p.raw++; saveStore(); if (p.raw >= 5) award('raw_first');
  if (!window.Plotly) return;
  const traces = GENOTYPES.map((g, gi) => {
    const v = vals.filter(x => x.genotype === g && x.value != null);
    return { type: 'scatter', mode: 'markers', name: g, x: v.map(x => MEDIA.indexOf(x.medium) + (gi - 1) * 0.22 + (Math.random() - 0.5) * 0.08), y: v.map(x => x.value), marker: { size: 9, color: GENO_COL[g], line: { width: 1, color: '#fff' } }, hovertemplate: `${g} %{y:.2f}<extra></extra>` };
  });
  const missing = vals.filter(x => x.value == null).length;
  Plotly.react('cardPlot', traces, { height: 260, margin: { l: 45, r: 10, t: 25, b: 35 }, title: { text: missing ? `${missing}/27 samples not detected (NA, not 0)` : 'detected in all 27 samples', font: { size: 12 } },
    xaxis: { tickvals: [0, 1, 2], ticktext: MEDIA, range: [-0.5, 2.5] }, yaxis: { title: 'log2 abundance' }, legend: { orientation: 'h', y: -0.25 }, paper_bgcolor: 'rgba(0,0,0,0)' }, { displaylogo: false, responsive: true });
}
const STRUCT_CACHE = new Map();
async function fetchAlphaFold(acc) {
  if (STRUCT_CACHE.has(acc)) return STRUCT_CACHE.get(acc);
  const res = await fetch(`https://alphafold.ebi.ac.uk/api/prediction/${encodeURIComponent(acc)}`);
  if (!res.ok) throw new Error(`AlphaFold has no model for ${acc} (${res.status})`);
  const entry = (await res.json())[0];
  if (!entry.uniprotSequence) { const f = await fetch(`https://rest.uniprot.org/uniprotkb/${encodeURIComponent(acc)}.fasta`); entry.uniprotSequence = (await f.text()).split('\n').slice(1).join(''); }
  entry.pdbText = await (await fetch(entry.pdbUrl)).text();
  STRUCT_CACHE.set(acc, entry); return entry;
}
async function showStructure(e, viewerEl, msgEl, withSeq = false) {
  msgEl.textContent = '🧊 Fetching AlphaFold model…';
  try {
    const entry = await fetchAlphaFold(e.uniprot); e.refSeq = entry.uniprotSequence;
    if (withSeq) renderSeq(e);
    if (!window.$3Dmol) { msgEl.textContent = '3Dmol.js could not load — sequence only.'; return; }
    if (!viewerEl.isConnected) return; // user navigated away
    viewerEl.innerHTML = '';
    const viewer = $3Dmol.createViewer(viewerEl, { backgroundColor: 'white' });
    viewer.addModel(entry.pdbText, 'pdb');
    const plddt = b => b > 90 ? '#0053d6' : b > 70 ? '#65cbf3' : b > 50 ? '#ffdb13' : '#ff7d45';
    viewer.setStyle({}, { cartoon: { colorfunc: a => plddt(a.b) } });
    const sites = sitesFor(e).filter(s => s.matchesRef);
    for (const s of sites) viewer.addStyle({ resi: s.position }, { sphere: { radius: 1.8, color: '#d6007a' } });
    viewer.zoomTo(); viewer.spin('y', 0.6); viewer.render();
    viewerEl.onclick = () => viewer.spin(false);
    const p = me(); if (!p.structs.includes(e.id)) { p.structs.push(e.id); saveStore(); addXP(20, 'structure unlocked'); if (p.structs.length >= 3) award('struct3'); renderDex(); }
    msgEl.textContent = `${entry.modelEntityId || 'AlphaFold'} · ${entry.uniprotSequence.length} aa · colour = pLDDT confidence (blue high → orange low)${sites.length ? ` · ${sites.length} phosphosite(s) as pink spheres` : ''} · click to stop spinning`;
  } catch (err) { msgEl.textContent = String(err.message || err); }
}
function sitesFor(e) {
  const list = DATA.sites.get(e.key) || DATA.sites.get(e.locus) || [];
  const je2 = DATA.seqs.get(e.key) || DATA.seqs.get(e.locus);
  return list.map(s => ({ ...s, matchesJE2: je2 ? je2[s.position - 1] === s.residue : null, matchesRef: e.refSeq ? e.refSeq[s.position - 1] === s.residue : false }));
}
function renderSeq(e) {
  const je2 = DATA.seqs.get(e.key) || DATA.seqs.get(e.locus); const seq = je2 || e.refSeq; if (!seq || !$('#seq')) return;
  const sites = sitesFor(e); const siteAt = new Map(sites.map(s => [s.position, s]));
  const motifMark = new Set(); let motifN = 0; const pat = $('#motif')?.value;
  if (pat) { try { const re = new RegExp(pat, 'gi'); let m; while ((m = re.exec(seq)) && motifN < 500) { for (let i = m.index; i < m.index + Math.max(1, m[0].length); i++) motifMark.add(i); motifN++; if (!m[0].length) re.lastIndex++; } } catch (err) { motifN = -1; } }
  let html = '';
  for (let i = 0; i < seq.length; i++) {
    if (i % 60 === 0) html += `${i ? '\n' : ''}<span class="pos">${String(i + 1).padStart(5)}</span> `; else if (i % 10 === 0) html += ' ';
    const s = siteAt.get(i + 1); const cls = [s ? 'site' : '', motifMark.has(i) ? 'motif' : ''].join(' ').trim();
    html += cls ? `<span class="${cls}" title="${s ? `${s.residue}${s.position}${Number.isFinite(s.prob) ? ` · loc. prob ${s.prob.toFixed(2)}` : ''}` : ''}">${seq[i]}</span>` : seq[i];
  }
  const note = je2 ? 'JE2 sequence from your FASTA.' : 'Reference (NCTC8325) sequence from UniProt/AlphaFold — JE2 numbering can differ.';
  const mism = sites.filter(s => !je2 && !s.matchesRef).map(s => s.label);
  $('#seq').innerHTML = `<p class="muted">${seq.length} aa · ${note} ${sites.length ? `· ${sites.length} measured phosphosite(s) in pink` : ''}${mism.length ? ` · ⚠ residue mismatch in reference at ${esc(mism.join(', '))}` : ''}${motifN > 0 ? ` · motif matches: ${motifN}` : motifN < 0 ? ' · invalid regex' : ''}</p><pre>${html}</pre>`;
}

/* ============================ 10. modal, controls, wiring ============================ */
function showModal(html, wire) { $('#modal').innerHTML = `<div class="modalbox">${html}</div>`; $('#modal').classList.add('open'); wire($('#modal .modalbox')); }
function closeModal() { $('#modal').classList.remove('open'); $('#modal').innerHTML = ''; }
function setContrast(c) { UI.contrast = c; $('#contrastSel').value = c; }
function refreshAll() {
  const byGroup = new Map(); for (const c of DATA.contrasts) { const g = contrastGroup(c); if (!byGroup.has(g)) byGroup.set(g, []); byGroup.get(g).push(c); }
  $('#contrastSel').innerHTML = [...byGroup.entries()].map(([g, cs]) => `<optgroup label="${esc(g)} (${cs.length})">${cs.map(c => `<option>${esc(c)}</option>`).join('')}</optgroup>`).join('');
  setContrast(DATA.contrasts.includes('TSB_vs_SN_at_anc') ? 'TSB_vs_SN_at_anc' : DATA.contrasts[0]);
  $('#nContrasts').textContent = `${DATA.contrasts.length} contrasts`;
  $('#dataBadge').className = `databadge ${DATA.mode}`;
  $('#dataBadge').textContent = DATA.mode === 'demo' ? `DEMO · ${DATA.entities.size.toLocaleString()} proteins · simulated values` : `REAL · ${DATA.dataset} · local only`;
  UI.highlight = new Set(); UI.px = { cat: null, pw: null, id: null };
  computeGroups();
  renderVolcano(); sniperIntro(); if (CG.running) catcherOver(true); else catcherIntro(); renderExplorer(); renderDex(); renderBadges(); renderHUD(); renderSearch();
  $('#card').innerHTML = '';
}
function renderSearch() {
  const q = ($('#search').value || '').toLowerCase(); const onlyHits = $('#onlyHits').checked;
  const list = rowsNow().filter(r => { const e = DATA.entities.get(r.id); return (!q || `${e.name} ${e.locus} ${e.product} ${e.uniprot}`.toLowerCase().includes(q)) && (!onlyHits || isHit(r)); });
  list.sort((a, b) => a.fdr - b.fdr);
  $('#results').innerHTML = `<p class="muted">${list.length.toLocaleString()} match(es)${list.length > 200 ? ' · showing 200 lowest FDR' : ''}</p><table><tr><th>name</th><th>product</th><th>diff</th><th>FDR</th></tr>` +
    list.slice(0, 200).map(r => { const e = DATA.entities.get(r.id); return `<tr data-id="${esc(r.id)}" class="${isHit(r) ? 'hit' : ''}"><td>${esc(e.name)}</td><td>${esc(e.product).slice(0, 55)}</td><td>${fmt(r.diff)}${r.lod ? ' ⚠' : ''}</td><td>${fmtP(r.fdr)}</td></tr>`; }).join('') + '</table>';
  $('#results').querySelectorAll('tr[data-id]').forEach(tr => tr.onclick = () => openCard(tr.dataset.id));
}
function exportHits() {
  const info = contrastInfo(UI.contrast);
  const lines = [['name', 'locus_tag', 'product', 'uniprot', 'contrast', 'diff', 'FDR', 'lod_imputed', 'direction'].join(',')];
  for (const r of rowsNow().filter(isHit)) { const e = DATA.entities.get(r.id); lines.push([e.name, e.locus, e.product, e.uniprot, UI.contrast, r.diff.toFixed(4), r.fdr.toExponential(3), r.lod, r.diff > 0 ? info.pos : info.neg].map(v => `"${String(v).replace(/"/g, '""')}"`).join(',')); }
  const a = document.createElement('a'); a.href = URL.createObjectURL(new Blob([lines.join('\n')], { type: 'text/csv' })); a.download = `hits_${UI.contrast}_FDR${UI.fdr}_LFC${UI.lfc}.csv`; a.click();
}
function switchTab(id) {
  if (UI.tab === 'tCatcher' && id !== 'tCatcher' && CG.running) catcherOver(true);
  if (UI.tab === 'tSniper' && id !== 'tSniper' && SN.running) sniperIntro('<p class="muted">Game stopped — press play to start again.</p>');
  UI.tab = id; SFX.play('click');
  document.querySelectorAll('.tabs button').forEach(x => x.classList.toggle('active', x.dataset.tab === id));
  document.querySelectorAll('.tabpane').forEach(p => { p.hidden = p.id !== id; });
  if (id === 'tExplore') renderExplorer();
  if (id === 'tCatcher' && !CG.running) catcherIntro();
}

async function init() {
  loadStore(); me();
  $('#player').value = GAME.current;
  $('#player').onchange = () => { GAME.current = $('#player').value.trim() || 'Player 1'; me(); saveStore(); renderHUD(); renderBadges(); renderDex(); sniperIntro(); catcherIntro(); };
  $('#muteBtn').onclick = () => { GAME.muted = !GAME.muted; saveStore(); renderHUD(); if (!GAME.muted) SFX.play('xp'); };
  $('#contrastSel').onchange = () => { UI.contrast = $('#contrastSel').value; UI.highlight = new Set(); renderVolcano(); renderSearch(); if (UI.tab === 'tExplore') renderExplorer(); };
  $('#fdrSlider').oninput = () => { UI.fdr = +$('#fdrSlider').value; syncSliders(); renderVolcano(); renderSearch(); };
  $('#lfcSlider').oninput = () => { UI.lfc = +$('#lfcSlider').value; syncSliders(); renderVolcano(); renderSearch(); };
  $('#resetCut').onclick = () => { UI.fdr = 0.1; UI.lfc = 1; UI.show = { up: true, down: true, ns: true, lod: true }; syncSliders(); renderVolcano(); renderSearch(); };
  $('#exportBtn').onclick = exportHits;
  $('#search').oninput = renderSearch; $('#onlyHits').onchange = renderSearch;
  $('#modal').onclick = ev => { if (ev.target.id === 'modal') closeModal(); };
  document.querySelectorAll('.tabs button').forEach(b => b.onclick = () => switchTab(b.dataset.tab));
  document.addEventListener('keydown', ev => {
    if (/INPUT|SELECT|TEXTAREA/.test(document.activeElement?.tagName)) return;
    if (UI.tab === 'tSniper' && SN.running && !SN.answered) { if (ev.key === 'ArrowLeft') { sniperAnswer('pos'); ev.preventDefault(); } if (ev.key === 'ArrowRight') { sniperAnswer('neg'); ev.preventDefault(); } }
    if (UI.tab === 'tCatcher' && CG.running) { if (ev.key === 'ArrowLeft') { CG.px = Math.max(CG.pw / 2, CG.px - 26); ev.preventDefault(); } if (ev.key === 'ArrowRight') { CG.px = Math.min(CG.W - CG.pw / 2, CG.px + 26); ev.preventDefault(); } }
  });
  document.addEventListener('pointerdown', () => SFX.ac(), { once: true }); // browsers need a click before audio
  $('#folder').onchange = async ev => {
    const files = ev.target.files; if (!files.length) return;
    if (await loadCourseFolder(files, $('#datasetSel').value)) { refreshAll(); award('real'); }
  };
  $('#datasetSel').onchange = async () => {
    if (!DATA.files) { $('#loadStatus').textContent = 'Choose your course data folder first.'; return; }
    if (await loadCourseFolder(DATA.files, $('#datasetSel').value)) refreshAll();
  };
  $('#demoBtn').onclick = () => { buildDemo(); refreshAll(); $('#loadStatus').textContent = 'Demo data: simulated values on the real S. aureus reference proteome.'; };
  syncSliders(); renderHUD();
  $('#loadStatus').textContent = 'Loading the S. aureus NCTC8325 reference proteome from UniProt (public)…';
  try { const n = await loadProteome(); $('#loadStatus').textContent = `Demo: ${n.toLocaleString()} UniProt reference proteins with simulated values. Pick your local course data or derived folder to play with the real data — files are read in your browser, never uploaded.`; }
  catch (err) { $('#loadStatus').textContent = 'UniProt unreachable (offline?) — using the small offline demo. Pathways fall back to name-based families.'; }
  buildDemo(); refreshAll();
}
window.addEventListener('DOMContentLoaded', init);
