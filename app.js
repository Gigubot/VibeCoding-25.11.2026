/* JE2 Proteome Quest — a gamified explorer for the Bio253 S. aureus JE2 data set.
 * Runs entirely in the browser. Course data are read from YOUR disk and never uploaded.
 * Sign rule used everywhere: contrast A_vs_B  ->  diff = log2(A) - log2(B).
 */
'use strict';

/* ======================= 1. small utilities ======================= */
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

/* ======================= 2. design + contrast language ======================= */
const GENOTYPES = ['anc', 'stk', 'stp'];
const MEDIA = ['TSB', 'SNpulse', 'SN'];
const GROUPS = GENOTYPES.flatMap(g => MEDIA.map(m => `${g}_${m}`));
const GENO_COL = { anc: '#0072B2', stk: '#E69F00', stp: '#009E73' };
const NICE = { anc: 'ancestor (anc)', stk: 'kinase mutant (stk)', stp: 'phosphatase mutant (stp)', TSB: 'TSB', SN: 'PASN (SN, ~15 h)', SNpulse: 'PASN pulse (15 min)' };

function parseContrast(c) {
  let m = c.match(/^interaction_([A-Za-z]+)_vs_([A-Za-z]+)_at_([A-Za-z]+)_vs_([A-Za-z]+)$/);
  if (m) return { type: 'int', A: m[1], B: m[2], M1: m[3], M2: m[4] };
  m = c.match(/^([A-Za-z]+)_vs_([A-Za-z]+)_at_([A-Za-z]+)$/);
  if (m) return { type: 'simple', A: m[1], B: m[2], ctx: m[3] };
  m = c.match(/^([A-Za-z]+)_vs_([A-Za-z]+)$/);
  if (m) return { type: 'simple', A: m[1], B: m[2], ctx: '' };
  return { type: 'other' };
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

/* Protein families — defined only from product names, so they work for every protein. */
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

/* ======================= 3. data layer ======================= */
const DATA = {
  mode: 'demo', dataset: 'total', isSite: false,
  entities: new Map(),   // id -> {id,key,name,locus,product,uniprot,family,site}
  rows: new Map(),       // contrast -> [{id,diff,p,fdr,lod}]
  index: new Map(),      // contrast -> Map(id -> row)
  contrasts: [],
  abund: new Map(),      // id -> [{genotype,medium,replicate,value}]
  sites: new Map(),      // protein_key or locus_tag -> [{residue,position,prob,label}]
  seqs: new Map(),       // protein_key or locus_tag -> JE2 sequence (from FASTA)
  files: null,           // the user's folder (FileList) for switching datasets
  demoMeans: null,
};

function indexRows() {
  DATA.index = new Map(); DATA.contrasts = [...DATA.rows.keys()];
  for (const [c, arr] of DATA.rows) DATA.index.set(c, new Map(arr.map(r => [r.id, r])));
}

/* ---- 3a. simulated demo: consistent group means for every protein, then contrasts ---- */
const DEMO_CONTRASTS = ['TSB_vs_SN_at_anc', 'TSB_vs_SN_at_stk', 'TSB_vs_SN_at_stp', 'TSB_vs_SNpulse_at_anc', 'SNpulse_vs_SN_at_anc',
  'stk_vs_anc_at_TSB', 'stp_vs_anc_at_TSB', 'stk_vs_anc_at_SN', 'stp_vs_anc_at_SN',
  'interaction_stk_vs_anc_at_TSB_vs_SN', 'interaction_stp_vs_anc_at_TSB_vs_SN', 'interaction_stp_vs_stk_at_TSB_vs_SN'];
const DEMO_TEMPLATES = [
  i => `50S ribosomal protein L${i}`, i => `30S ribosomal protein S${i}`,
  i => `${pick(['alanine', 'glycine', 'serine', 'leucine', 'valine', 'proline', 'lysine', 'arginine', 'tyrosine', 'histidine'], mulberry32(i))}--tRNA ligase`,
  () => 'ATP synthase subunit', () => 'ABC transporter ATP-binding protein', () => 'ABC transporter permease',
  () => 'penicillin-binding protein', () => 'thioredoxin', () => 'iron ABC transporter substrate-binding protein',
  () => 'leukocidin family toxin', () => 'cytochrome aa3 quinol oxidase subunit', () => 'enolase',
  () => 'hypothetical protein', () => 'hypothetical protein', () => 'hypothetical protein', () => 'transcriptional regulator',
  () => 'MFS transporter', () => 'amino acid permease', () => 'GNAT family N-acetyltransferase', () => 'DUF family protein',
  () => 'hypothetical protein', () => 'membrane protein', () => 'hydrolase', () => 'oxidoreductase', () => 'lipoprotein',
  () => 'hypothetical protein', () => 'transporter', () => 'SDR family oxidoreductase', () => 'DNA-binding protein', () => 'hypothetical protein'];

function buildDemo() {
  const r = mulberry32(253);
  const list = []; const seen = new Set();
  (window.DEMO_PROTEINS || []).forEach(([uniprot, locus, name, product]) => {
    if (seen.has(locus)) return; seen.add(locus);
    list.push({ id: locus, key: locus, name, locus, product, uniprot, site: null });
  });
  for (let i = 1; list.length < 2300; i++) {
    const id = `DEMO_${String(i).padStart(4, '0')}`;
    list.push({ id, key: id, name: id, locus: id, product: DEMO_TEMPLATES[i % DEMO_TEMPLATES.length](i % 30 + 1) + ' (simulated)', uniprot: '', site: null });
  }
  // family-level "biology" so that games have something to discover
  const FAM_SHIFT = { 'Ribosomal proteins': { R: -1.2, stk: 0, stp: 0 }, 'Purine biosynthesis': { R: -0.8, stk: -1.0, stp: -0.2 },
    'Respiration': { R: 0.6, stk: 0.9, stp: 0 }, 'Oxidative stress': { R: 1.0, stk: 0, stp: 0.4 }, 'Iron & heme uptake': { R: 1.4, stk: 0, stp: 0 },
    'Virulence & toxins': { R: -1.3, stk: 0.3, stp: -0.5 }, 'Aminoacyl-tRNA ligases': { R: -0.4, stk: 0, stp: -0.6 } };
  const means = new Map();
  DATA.entities = new Map(); DATA.rows = new Map(); DATA.abund = new Map(); DATA.sites = new Map(); DATA.seqs = new Map();
  for (const p of list) {
    p.family = familyOf(p.product);
    const fs = FAM_SHIFT[p.family] || { R: 0, stk: 0, stp: 0 };
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
    p.lod = p.missing.length > 0;     // like the course export: the flag is model-wide
    p.se = 0.16 + 0.28 * r();
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
    const grp = `${g}_${md}`;
    out.push({ genotype: g, medium: md, replicate: rep, value: p.missing.includes(grp) ? null : m[grp] + gauss(r) * p.se * 1.2 });
  }
  return out;
}

/* ---- 3b. the real course data, read locally from the user's folder ---- */
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
      const gene = a.staa8_gene || '';
      const base = gene || a.locus_tag || d.protein_key;
      const site = isSite && d.residue ? { residue: d.residue, position: +d.position } : null;
      ents.set(id, { id, key: d.protein_key, locus: a.locus_tag || '', product: a.product || '', uniprot: a.staa8_uniprot || '',
        name: site ? `${base} ${site.residue}${site.position}` : base, family: familyOf(a.product), site });
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
  // phosphosites: from a phospho DEA itself, or from sites_STYH.tsv.gz
  const addSite = (key, residue, position, prob) => {
    if (!key || !residue) return; const list = DATA.sites.get(key) || []; if (list.some(s => s.position === position && s.residue === residue)) return;
    list.push({ residue, position, prob, label: `${residue}${position}` }); DATA.sites.set(key, list);
  };
  if (isSite) for (const e of ents.values()) if (e.site) addSite(e.key, e.site.residue, e.site.position, NaN);
  const siteFile = findFile(files, /sites_STYH\.tsv(\.gz)?$/);
  if (siteFile) {
    status('Reading phosphosites (sites_STYH)…');
    const seenSite = new Set();
    for (const s of parseTSV(await readText(siteFile))) {
      if (seenSite.has(s.site_id)) continue; seenSite.add(s.site_id);
      addSite(s.protein_key, s.residue, +s.position, num(s.best_loc_prob));
    }
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

/* ======================= 4. game engine: XP, levels, badges ======================= */
const TITLES = ['Pipette Rookie', 'Spectrum Sleuth', 'Volcano Ranger', 'FDR Guardian', 'Phospho Master', 'Proteome Legend'];
const BADGES = {
  first_card: ['🔬', 'First look', 'Open your first protein card'],
  raw_first: ['📊', 'Raw-data respecter', 'Look at per-sample values for 5 proteins'],
  sniper5: ['🎯', 'Sharp shooter', '5-answer streak in Sign Sniper'],
  sniper25: ['🏹', 'Sign sensei', '25 correct Sign Sniper answers'],
  int_master: ['🧠', 'Interaction thinker', 'Answer 5 interaction questions correctly'],
  lod3: ['🪤', 'Trap spotter', 'Recognise 3 LOD-imputed estimates'],
  reformed: ['😇', 'Reformed p-hacker', 'Loosen FDR above 0.25, then return to ≤ 0.1'],
  family: ['🧬', 'Family reunion', 'Complete a Family Hunt'],
  detective: ['🕵️', 'Detective', 'Solve a Mystery Protein with ≤ 2 clues'],
  dex10: ['📖', 'Collector', '10 entries in your Pokédex'],
  dex50: ['📚', 'Archivist', '50 entries in your Pokédex'],
  struct3: ['🧊', 'Structural biologist', 'View 3 AlphaFold structures'],
  legend: ['👑', 'Legendary pair', 'Collect both Stk1 and Stp1'],
  real: ['💾', 'Real scientist', 'Load the real course data'],
};
const STORE_KEY = 'je2-proteome-quest-v1';
const GAME = { players: {}, current: 'Player 1' };
function loadStore() { try { const s = JSON.parse(localStorage.getItem(STORE_KEY) || 'null'); if (s) Object.assign(GAME, s); } catch (e) { /* storage unavailable: keep in memory */ } }
function saveStore() { try { localStorage.setItem(STORE_KEY, JSON.stringify(GAME)); } catch (e) { /* ignore */ } }
function me() {
  const p = GAME.players[GAME.current] ||= {};
  for (const [k, v] of Object.entries({ xp: 0, badges: [], dex: [], structs: [], raw: 0, sniperOk: 0, streak: 0, intOk: 0, lodOk: 0, wasHacker: false })) if (p[k] === undefined) p[k] = Array.isArray(v) ? [] : v;
  return p;
}
const level = xp => Math.min(TITLES.length, 1 + Math.floor(Math.sqrt(xp / 60)));
function toast(html, kind = '') { const t = document.createElement('div'); t.className = `toast ${kind}`; t.innerHTML = html; $('#toasts').appendChild(t); setTimeout(() => t.remove(), 3800); }
function addXP(n, why) {
  const p = me(); const before = level(p.xp); p.xp += n; saveStore(); renderHUD();
  toast(`+${n} XP · ${esc(why)}`);
  if (level(p.xp) > before) toast(`⬆️ Level ${level(p.xp)}: <b>${TITLES[level(p.xp) - 1]}</b>`, 'gold');
}
function award(id) {
  const p = me(); if (p.badges.includes(id)) return; p.badges.push(id); saveStore();
  const [icon, name] = BADGES[id]; toast(`${icon} Badge unlocked: <b>${name}</b>`, 'gold'); addXP(50, 'badge'); renderBadges();
}
function renderHUD() {
  const p = me(), L = level(p.xp), lo = 60 * (L - 1) ** 2, hi = 60 * L ** 2;
  $('#hudLevel').textContent = `Lv ${L} · ${TITLES[L - 1]}`;
  $('#hudXP').textContent = `${p.xp} XP`;
  $('#hudBar').style.width = `${L >= TITLES.length ? 100 : Math.round(100 * (p.xp - lo) / (hi - lo))}%`;
  $('#hudBadges').textContent = `${p.badges.length}/${Object.keys(BADGES).length} 🏅`;
}
function renderBadges() {
  const p = me();
  $('#badgeGrid').innerHTML = Object.entries(BADGES).map(([id, [icon, name, how]]) =>
    `<div class="badge ${p.badges.includes(id) ? 'got' : ''}" title="${esc(how)}"><span>${icon}</span><b>${esc(name)}</b><small>${esc(how)}</small></div>`).join('');
  const board = Object.entries(GAME.players).map(([n, q]) => [n, q.xp || 0, (q.badges || []).length]).sort((a, b) => b[1] - a[1]);
  $('#leaderboard').innerHTML = board.map(([n, xp, b], i) => `<tr class="${n === GAME.current ? 'me' : ''}"><td>${['🥇', '🥈', '🥉'][i] || i + 1}</td><td>${esc(n)}</td><td>${xp}</td><td>${b}</td></tr>`).join('');
}

/* ======================= 5. live volcano + false-discovery meter ======================= */
const UI = { contrast: '', fdr: 0.1, lfc: 1, showLod: true, game: 'explore', hunt: null, highlight: new Set(), plotted: false };
const rowsNow = () => DATA.rows.get(UI.contrast) || [];
const isHit = r => r.fdr < UI.fdr && Math.abs(r.diff) > UI.lfc;
const yOf = r => -Math.log10(Math.max(r.fdr, 1e-300));

function renderVolcano() {
  if (!window.Plotly) { $('#volcano').innerHTML = '<p class="warn">Plotly could not load (offline?). The games still work.</p>'; renderCounts(); return; }
  const info = contrastInfo(UI.contrast);
  const groups = { up: [], down: [], ns: [], lod: [], hl: [] };
  for (const r of rowsNow()) {
    if (r.lod && !UI.showLod) continue;
    if (UI.highlight.has(r.id)) groups.hl.push(r);
    if (r.lod) groups.lod.push(r); else groups[isHit(r) ? (r.diff > 0 ? 'up' : 'down') : 'ns'].push(r);
  }
  const hover = r => { const e = DATA.entities.get(r.id); return `<b>${esc(e.name)}</b><br>${esc(e.product).slice(0, 60)}<br>diff ${fmt(r.diff)} · FDR ${fmtP(r.fdr)}${r.lod ? '<br>⚠ lod_imputed' : ''}`; };
  const tr = (arr, name, marker) => ({ type: 'scattergl', mode: 'markers', name, x: arr.map(r => r.diff), y: arr.map(yOf), customdata: arr.map(r => r.id), text: arr.map(hover), hovertemplate: '%{text}<extra></extra>', marker });
  const traces = [
    tr(groups.ns, 'not selected', { size: 5, color: '#c3c9d2' }),
    tr(groups.up, info.pos, { size: 7, color: '#d55e00' }),
    tr(groups.down, info.neg, { size: 7, color: '#0072b2' }),
    tr(groups.lod, 'lod_imputed (×)', { size: 8, symbol: 'x-thin', line: { width: 1.5, color: '#7b5ea7' }, color: '#7b5ea7' }),
    tr(groups.hl, 'highlighted', { size: 14, color: 'rgba(0,0,0,0)', line: { width: 2.5, color: '#f5b700' } }),
  ];
  const yFdr = -Math.log10(UI.fdr);
  const layout = {
    margin: { l: 55, r: 10, t: 10, b: 45 }, hovermode: 'closest', dragmode: 'pan', uirevision: UI.contrast,
    xaxis: { title: 'log2 effect (diff)', zeroline: true, zerolinecolor: '#999' }, yaxis: { title: '−log10(FDR)' },
    legend: { orientation: 'h', y: -0.2, font: { size: 11 } }, paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: '#fff',
    shapes: [
      { type: 'line', xref: 'x', yref: 'paper', x0: -UI.lfc, x1: -UI.lfc, y0: 0, y1: 1, line: { color: '#e69f00', width: 2, dash: 'dash' } },
      { type: 'line', xref: 'x', yref: 'paper', x0: UI.lfc, x1: UI.lfc, y0: 0, y1: 1, line: { color: '#e69f00', width: 2, dash: 'dash' } },
      { type: 'line', xref: 'paper', yref: 'y', x0: 0, x1: 1, y0: yFdr, y1: yFdr, line: { color: '#009e73', width: 2, dash: 'dash' } }],
  };
  const config = { responsive: true, displaylogo: false, edits: { shapePosition: true }, modeBarButtonsToRemove: ['lasso2d', 'select2d'],
    toImageButtonOptions: { format: 'svg', filename: `volcano_${UI.contrast}` } };
  Plotly.react('volcano', traces, layout, config);
  if (!UI.plotted) {
    UI.plotted = true;
    const v = $('#volcano');
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
  renderCounts();
}
function syncSliders() {
  $('#fdrSlider').value = UI.fdr; $('#fdrVal').textContent = UI.fdr < 0.001 ? UI.fdr.toExponential(1) : UI.fdr.toFixed(3);
  $('#lfcSlider').value = UI.lfc; $('#lfcVal').textContent = UI.lfc.toFixed(2);
}
function renderCounts() {
  const info = contrastInfo(UI.contrast); const rows = rowsNow().filter(r => UI.showLod || !r.lod);
  const hits = rows.filter(isHit); const up = hits.filter(r => r.diff > 0).length, down = hits.length - up;
  $('#counts').innerHTML = `<div class="count up">▲ ${up}<small>${esc(info.pos)}</small></div><div class="count down">▼ ${down}<small>${esc(info.neg)}</small></div>`;
  $('#formula').innerHTML = `<code>${info.formula}</code><p>${esc(info.english)}</p>`;
  // false-discovery meter: with BH at level q, expected false discoveries among H hits ≤ q·H
  const exp = UI.fdr * hits.length; const pct = Math.min(100, 100 * UI.fdr / 0.5);
  const p = me(); let verdict;
  if (UI.fdr > 0.25 || UI.lfc < 0.3) { verdict = '🚨 p-hacking zone: many of these hits may be noise.'; p.wasHacker = true; }
  else if (UI.fdr <= 0.1 && UI.lfc >= 1) { verdict = '😇 Responsible cutoffs (course default: FDR < 0.1, |diff| > 1).'; if (p.wasHacker) { p.wasHacker = false; award('reformed'); } }
  else verdict = '🤔 Exploratory cutoffs — say so when you report them.';
  saveStore();
  $('#meter').innerHTML = `<div class="meterbar"><div style="width:${pct}%"></div></div>
    <p><b>Expected false hits ≤ ${exp.toFixed(1)}</b> of ${hits.length} &nbsp;<span class="muted">(FDR cutoff × hits = ${UI.fdr.toFixed(3)} × ${hits.length})</span></p><p>${verdict}</p>`;
}
function onPointClick(id) {
  if (UI.game === 'hunt' && UI.hunt) return huntClick(id);
  const r = DATA.index.get(UI.contrast)?.get(id);
  if (r && r.lod && !(UI.trapDone ||= new Set()).has(id)) return lodTrap(id);
  openCard(id);
}

/* ======================= 6. games ======================= */
/* 6a. LOD Trap — clicking an imputed estimate triggers a question before the card opens */
function lodTrap(id) {
  const e = DATA.entities.get(id);
  const opts = [['A precise fold change — I can trust the number', false],
    ['Presence/absence evidence: at least one group mean was imputed from the detection limit. Check the per-sample values!', true],
    ['It means the protein is a contaminant', false]].sort(() => Math.random() - 0.5);
  showModal(`<h3>🪤 LOD Trap!</h3><p><b>${esc(e.name)}</b> is flagged <code>estimate_type = lod_imputed</code>. What does that tell you?</p>
    ${opts.map(([t, ok], i) => `<button class="opt" data-ok="${ok}" data-i="${i}">${esc(t)}</button>`).join('')}`, box => {
    box.querySelectorAll('.opt').forEach(b => b.onclick = () => {
      UI.trapDone.add(id);
      if (b.dataset.ok === 'true') { const p = me(); p.lodOk++; addXP(15, 'trap spotted'); if (p.lodOk >= 3) award('lod3'); b.classList.add('ok'); }
      else { b.classList.add('bad'); toast('Not quite: imputed means "not measured in a group" — treat it as presence/absence.', 'bad'); }
      setTimeout(() => { closeModal(); openCard(id); }, 900);
    });
  });
}

/* 6b. Sign Sniper — which group has more? */
function newSniperRound() {
  const p = me(); const ints = DATA.contrasts.filter(c => parseContrast(c).type === 'int'); const simples = DATA.contrasts.filter(c => parseContrast(c).type === 'simple');
  const useInt = ints.length && p.streak >= 5 && Math.random() < 0.4;
  const pool = useInt ? ints : simples; if (!pool.length) { $('#sniper').innerHTML = '<p>No suitable contrasts loaded.</p>'; return; }
  let c, r, tries = 0;
  do { c = pick(pool); const cand = (DATA.rows.get(c) || []).filter(x => !x.lod && x.fdr < 0.1 && Math.abs(x.diff) > 1); r = pick(cand); tries++; } while (!r && tries < 20);
  if (!r) { $('#sniper').innerHTML = '<p>No clear hits to quiz on — try another dataset.</p>'; return; }
  const e = DATA.entities.get(r.id), info = contrastInfo(c);
  const answers = [[info.pos, r.diff > 0], [info.neg, r.diff < 0]];
  $('#sniper').innerHTML = `<div class="streak">🔥 streak ${p.streak} · ✅ ${p.sniperOk}</div>
    <p class="q">${useInt ? '<span class="tag hard">HARD: interaction</span>' : ''} <b>${esc(e.name)}</b> <span class="muted">${esc(e.product)}</span></p>
    <p>Contrast <code>${esc(c)}</code> &nbsp; diff = <b>${fmt(r.diff)}</b> &nbsp; FDR = ${fmtP(r.fdr)}</p>
    <p>What does this mean?</p>
    ${answers.map(([t, ok]) => `<button class="opt" data-ok="${ok}">${esc(t)}</button>`).join('')}
    <div id="sniperExplain"></div>`;
  $('#sniper').querySelectorAll('.opt').forEach(b => b.onclick = () => {
    $('#sniper').querySelectorAll('.opt').forEach(x => { x.disabled = true; if (x.dataset.ok === 'true') x.classList.add('ok'); });
    if (b.dataset.ok === 'true') {
      p.streak++; p.sniperOk++; if (useInt) p.intOk++; addXP(10 + Math.min(p.streak, 10) + (useInt ? 10 : 0), `Sign Sniper${p.streak > 1 ? ` ×${p.streak}` : ''}`);
      if (p.streak >= 5) award('sniper5'); if (p.sniperOk >= 25) award('sniper25'); if (p.intOk >= 5) award('int_master');
    } else { b.classList.add('bad'); p.streak = 0; }
    saveStore();
    $('#sniperExplain').innerHTML = `<div class="explain"><code>${info.formula}</code><p>${fmt(r.diff)} is ${r.diff > 0 ? 'positive' : 'negative'} → <b>${esc(r.diff > 0 ? info.pos : info.neg)}</b>.</p>
      <button id="sniperNext">Next ▶</button> <button id="sniperShow" class="ghost">Show on volcano</button></div>`;
    $('#sniperNext').onclick = newSniperRound;
    $('#sniperShow').onclick = () => { setContrast(c); UI.highlight = new Set([r.id]); renderVolcano(); };
  });
}

/* 6c. Mystery Protein — progressive clues */
function newMystery() {
  const named = [...DATA.entities.values()].filter(e => !/^DEMO_/.test(e.id) && e.product && !/hypothetical/i.test(e.product));
  const pool = named.length > 20 ? named : [...DATA.entities.values()];
  let e, best, tries = 0;
  do { e = pick(pool); best = null; for (const c of DATA.contrasts) { const r = DATA.index.get(c).get(e.id); if (r && !r.lod && (!best || r.fdr < best.r.fdr)) best = { c, r }; } tries++; }
  while ((!best || best.r.fdr > 0.05) && tries < 200);
  if (!best) return;
  const info = contrastInfo(best.c);
  const words = e.product.split(/\s+/);
  const clues = [
    `In <code>${esc(best.c)}</code> it is <b>${esc(best.r.diff > 0 ? info.pos : info.neg)}</b> (|diff| ${Math.abs(best.r.diff) > 2 ? '> 2' : '≤ 2'}, FDR ${best.r.fdr < 0.01 ? '< 0.01' : '< 0.05'}).`,
    e.family ? `Protein family: <b>${esc(e.family)}</b>.` : `Its product name has <b>${words.length}</b> word(s).`,
    `Product, masked: <b>${esc(words.map(w => w[0] + w.slice(1).replace(/[A-Za-z0-9]/g, '·')).join(' '))}</b>`,
    `Locus tag starts with <b>${esc((e.locus || e.id).slice(0, -2))}··</b>`,
    `Product: <b>${esc(e.product)}</b> — now name it!`];
  UI.mystery = { e, best, clues, shown: 1, done: false };
  renderMystery();
}
function renderMystery() {
  const M = UI.mystery; if (!M) return;
  const pts = Math.max(10, 60 - 10 * (M.shown - 1));
  $('#mystery').innerHTML = `<ol class="clues">${M.clues.slice(0, M.shown).map(c => `<li>${c}</li>`).join('')}</ol>
    ${M.done ? '' : `<div class="row"><input id="guess" list="entityList" placeholder="Type a gene, locus tag or product…" autocomplete="off">
    <button id="guessBtn">Guess (${pts} XP)</button><button id="clueBtn" class="ghost" ${M.shown >= M.clues.length ? 'disabled' : ''}>Next clue</button>
    <button id="giveUp" class="ghost">Reveal</button></div>`}<div id="mysteryMsg"></div>`;
  if (M.done) return;
  $('#clueBtn').onclick = () => { M.shown++; renderMystery(); };
  $('#giveUp').onclick = () => finishMystery(false);
  const go = () => {
    const g = $('#guess').value.trim().toLowerCase(); if (!g) return; const e = M.e;
    if ([e.name, e.locus, e.product, e.id].some(x => x && x.toLowerCase() === g)) {
      addXP(pts, 'Mystery solved'); if (M.shown <= 2) award('detective'); finishMystery(true);
    } else { $('#mysteryMsg').innerHTML = '<p class="bad">Nope — try again or take a clue.</p>'; }
  };
  $('#guessBtn').onclick = go; $('#guess').onkeydown = ev => { if (ev.key === 'Enter') go(); };
}
function finishMystery(won) {
  const M = UI.mystery; M.done = true; M.shown = M.clues.length; renderMystery();
  $('#mysteryMsg').innerHTML = `<p>${won ? '🎉 Correct!' : 'It was:'} <b>${esc(M.e.name)}</b> — ${esc(M.e.product)}</p>
    <button id="mysteryShow">Show on volcano + open card</button> <button id="mysteryNext" class="ghost">New mystery</button>`;
  $('#mysteryShow').onclick = () => { setContrast(M.best.c); UI.highlight = new Set([M.e.id]); renderVolcano(); openCard(M.e.id); };
  $('#mysteryNext').onclick = newMystery;
}

/* 6d. Family Hunt — find all members of a protein family that pass your cutoffs */
function renderHuntSetup() {
  const fams = FAMILIES.map(f => f[0]).filter(f => [...DATA.entities.values()].some(e => e.family === f));
  $('#hunt').innerHTML = `<p>Pick a family. Every member that passes your current cutoffs in the shown contrast is a target. Hover the volcano (product names!) and click to catch them.</p>
    <div class="row"><select id="huntFam">${fams.map(f => `<option>${esc(f)}</option>`).join('')}</select><button id="huntGo">Start hunt</button></div><div id="huntStatus"></div>`;
  $('#huntGo').onclick = () => startHunt($('#huntFam').value);
}
function startHunt(fam) {
  const targets = new Set(rowsNow().filter(r => isHit(r) && (UI.showLod || !r.lod) && DATA.entities.get(r.id).family === fam).map(r => r.id));
  if (!targets.size) { $('#huntStatus').innerHTML = `<p class="bad">No ${esc(fam)} member passes the cutoffs in <code>${esc(UI.contrast)}</code>. Loosen the cutoffs or switch contrast — or notice that a family can shift as a group without single members passing (that is what GSEA detects).</p>`; return; }
  const goal = Math.min(10, targets.size); // big families: catch any 10
  UI.game = 'hunt'; UI.hunt = { fam, targets, goal, found: new Set(), misses: 0 }; UI.highlight = new Set(); renderVolcano(); renderHuntStatus();
}
function huntClick(id) {
  const H = UI.hunt; const e = DATA.entities.get(id);
  if (H.targets.has(id) && !H.found.has(id)) { H.found.add(id); UI.highlight.add(id); addXP(8, `caught ${e.name}`); renderVolcano(); }
  else if (!H.targets.has(id)) { H.misses++; toast(`${esc(e.name)} is not a ${esc(H.fam)} hit (${esc(e.product).slice(0, 40)})`, 'bad'); }
  renderHuntStatus();
  if (H.found.size >= H.goal) { award('family'); addXP(40, 'family complete'); UI.game = 'explore'; UI.hunt = null; $('#huntStatus').innerHTML += '<p>🎉 Hunt complete! Click points again to open protein cards.</p>'; }
}
function renderHuntStatus() {
  const H = UI.hunt; if (!H) return;
  $('#huntStatus').innerHTML = `<p><b>${esc(H.fam)}</b>: caught ${H.found.size}/${H.goal}${H.targets.size > H.goal ? ` (${H.targets.size} members pass — any ${H.goal} count)` : ''} · misses ${H.misses}</p>
    <button id="huntHint" class="ghost">Hint (−5 XP)</button> <button id="huntStop" class="ghost">Stop hunt</button>`;
  $('#huntHint').onclick = () => { const left = [...H.targets].filter(t => !H.found.has(t)); if (!left.length) return; const p = me(); p.xp = Math.max(0, p.xp - 5); saveStore(); renderHUD(); const e = DATA.entities.get(pick(left)); toast(`Hint: look for “${esc(e.product).slice(0, 50)}” near diff ${fmt(DATA.index.get(UI.contrast).get(e.id).diff, 1)}`); };
  $('#huntStop').onclick = () => { UI.game = 'explore'; UI.hunt = null; UI.highlight = new Set(); renderVolcano(); renderHuntSetup(); };
}

const isLegend = e => !!e && (/^(stk1|stp1|pknB)$/i.test(e.name) || /^B7H15_0635[05]$/.test(e.locus));
/* 6e. Protein Pokédex + protein card (values, per-sample dots, sequence, AlphaFold 3D) */
function collect(id) {
  const p = me(); if (p.dex.includes(id)) return; p.dex.push(id); saveStore(); addXP(5, 'new Pokédex entry');
  if (p.dex.length === 1) award('first_card'); if (p.dex.length >= 10) award('dex10'); if (p.dex.length >= 50) award('dex50');
  const leg = p.dex.map(i => DATA.entities.get(i)).filter(isLegend).map(e => e.locus || e.name.toLowerCase());
  if (new Set(leg).size >= 2) award('legend');
  renderDex();
}
function renderDex() {
  const p = me(); const ents = p.dex.map(id => DATA.entities.get(id)).filter(Boolean);
  $('#dexCount').textContent = `${ents.length} collected in this data set`;
  $('#dexGrid').innerHTML = ents.map(e => `<button class="dex ${p.structs.includes(e.id) ? 'star' : ''} ${isLegend(e) ? 'legend' : ''}" data-id="${esc(e.id)}">
    <b>${esc(e.name)}</b><small>${esc(e.family || e.product).slice(0, 34)}</small></button>`).join('') || '<p class="muted">Open protein cards to fill your Pokédex.</p>';
  $('#dexGrid').querySelectorAll('.dex').forEach(b => b.onclick = () => openCard(b.dataset.id));
}
function openCard(id) {
  const e = DATA.entities.get(id); if (!e) return; collect(id);
  const rows = DATA.contrasts.map(c => [c, DATA.index.get(c).get(id)]).filter(x => x[1]);
  const up = e.uniprot;
  $('#card').innerHTML = `<div class="cardhead"><div><h2>${esc(e.name)} ${isLegend(e) ? '👑' : ''}</h2>
      <p>${esc(e.product)}</p><p class="muted">${esc(e.locus || e.key)}${up ? ` · UniProt ${esc(up)} (NCTC8325 ortholog)` : ''}${e.family ? ` · <span class="tag">${esc(e.family)}</span>` : ''}${DATA.mode === 'demo' ? ' · <span class="tag warn">simulated values</span>' : ''}</p></div>
      <button class="ghost" id="cardClose">✕</button></div>
    <div class="links">${up ? `<a target="_blank" rel="noopener" href="https://www.uniprot.org/uniprotkb/${encodeURIComponent(up)}/entry">UniProt</a>
      <a target="_blank" rel="noopener" href="https://alphafold.ebi.ac.uk/entry/${encodeURIComponent(up)}">AlphaFold DB</a>
      <a target="_blank" rel="noopener" href="https://string-db.org/cgi/network?identifiers=${encodeURIComponent(up)}&species=93061">STRING</a>` : '<span class="muted">No UniProt ortholog mapped.</span>'}</div>
    <div class="cardgrid"><div><h4>Per-sample values <span class="muted">(look before trusting a fold change)</span></h4><div id="cardPlot"></div></div>
      <div><h4>All contrasts</h4><div class="tablewrap"><table><tr><th>contrast</th><th>diff</th><th>FDR</th><th>meaning</th></tr>
      ${rows.map(([c, r]) => { const i = contrastInfo(c); const hit = r.fdr < 0.1 && Math.abs(r.diff) > 1; return `<tr class="${hit ? 'hit' : ''}"><td><code>${esc(c)}</code></td><td>${fmt(r.diff)}${r.lod ? ' <span title="lod_imputed">⚠</span>' : ''}</td><td>${fmtP(r.fdr)}</td><td>${esc(r.diff > 0 ? i.pos : i.neg)}</td></tr>`; }).join('')}
      </table></div></div></div>
    <h4>Sequence & AlphaFold structure</h4>
    <div class="row"><button id="structBtn" ${up ? '' : 'disabled'}>🧊 Load sequence + 3D model</button>
      <input id="motif" placeholder="motif regex, e.g. R..S or [ST]P" ${up || DATA.seqs.has(e.key) ? '' : 'disabled'}><span id="structMsg" class="muted"></span></div>
    <div class="seqgrid"><div id="seq"></div><div id="viewer"></div></div>`;
  $('#cardClose').onclick = () => { $('#card').innerHTML = ''; };
  $('#structBtn').onclick = () => loadStructure(e);
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
async function loadStructure(e) {
  const msg = t => { $('#structMsg').textContent = t; };
  msg('Fetching AlphaFold model…');
  try {
    let entry = STRUCT_CACHE.get(e.uniprot);
    if (!entry) {
      const res = await fetch(`https://alphafold.ebi.ac.uk/api/prediction/${encodeURIComponent(e.uniprot)}`);
      if (!res.ok) throw new Error(`AlphaFold has no model for ${e.uniprot} (${res.status})`);
      entry = (await res.json())[0];
      if (!entry.uniprotSequence) {
        const f = await fetch(`https://rest.uniprot.org/uniprotkb/${encodeURIComponent(e.uniprot)}.fasta`);
        entry.uniprotSequence = (await f.text()).split('\n').slice(1).join('');
      }
      entry.pdbText = await (await fetch(entry.pdbUrl)).text();
      STRUCT_CACHE.set(e.uniprot, entry);
    }
    e.refSeq = entry.uniprotSequence;
    renderSeq(e);
    if (!window.$3Dmol) { msg('3Dmol.js could not load — sequence only.'); return; }
    $('#viewer').innerHTML = '';
    const viewer = $3Dmol.createViewer($('#viewer'), { backgroundColor: 'white' });
    viewer.addModel(entry.pdbText, 'pdb');
    const plddt = b => b > 90 ? '#0053d6' : b > 70 ? '#65cbf3' : b > 50 ? '#ffdb13' : '#ff7d45';
    viewer.setStyle({}, { cartoon: { colorfunc: a => plddt(a.b) } });
    for (const s of sitesFor(e)) if (s.matchesRef) viewer.addStyle({ resi: s.position }, { sphere: { radius: 1.8, color: '#d6007a' } });
    viewer.zoomTo(); viewer.render();
    const p = me(); if (!p.structs.includes(e.id)) { p.structs.push(e.id); saveStore(); addXP(20, 'structure viewed'); if (p.structs.length >= 3) award('struct3'); renderDex(); }
    msg(`AlphaFold ${entry.modelEntityId || ''} · colour = pLDDT confidence (blue high → orange low) · pink spheres = phosphosites`);
  } catch (err) { msg(String(err.message || err)); }
}
function sitesFor(e) {
  const list = DATA.sites.get(e.key) || DATA.sites.get(e.locus) || [];
  const je2 = DATA.seqs.get(e.key) || DATA.seqs.get(e.locus);
  return list.map(s => ({ ...s, matchesJE2: je2 ? je2[s.position - 1] === s.residue : null, matchesRef: e.refSeq ? e.refSeq[s.position - 1] === s.residue : false }));
}
function renderSeq(e) {
  const je2 = DATA.seqs.get(e.key) || DATA.seqs.get(e.locus); const seq = je2 || e.refSeq; if (!seq) return;
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

/* ======================= 7. modal, controls, wiring ======================= */
function showModal(html, wire) { $('#modal').innerHTML = `<div class="modalbox">${html}</div>`; $('#modal').classList.add('open'); wire($('#modal .modalbox')); }
function closeModal() { $('#modal').classList.remove('open'); $('#modal').innerHTML = ''; }

function setContrast(c) { UI.contrast = c; $('#contrastSel').value = c; }
function refreshAll() {
  $('#contrastSel').innerHTML = DATA.contrasts.map(c => `<option>${esc(c)}</option>`).join('');
  setContrast(DATA.contrasts.includes('TSB_vs_SN_at_anc') ? 'TSB_vs_SN_at_anc' : DATA.contrasts[0]);
  $('#entityList').innerHTML = [...DATA.entities.values()].slice(0, 6000).map(e => `<option value="${esc(e.name)}">${esc(e.product)}</option>`).join('');
  $('#dataBadge').className = `databadge ${DATA.mode}`;
  $('#dataBadge').textContent = DATA.mode === 'demo' ? 'DEMO · simulated values' : `REAL · ${DATA.dataset} · local only`;
  UI.highlight = new Set(); UI.game = 'explore'; UI.hunt = null;
  renderVolcano(); newSniperRound(); newMystery(); renderHuntSetup(); renderDex(); renderBadges(); renderHUD(); renderSearch();
  $('#card').innerHTML = '';
}
function renderSearch() {
  const q = ($('#search').value || '').toLowerCase(); const onlyHits = $('#onlyHits').checked;
  let list = rowsNow().filter(r => { const e = DATA.entities.get(r.id); return (!q || `${e.name} ${e.locus} ${e.product} ${e.uniprot}`.toLowerCase().includes(q)) && (!onlyHits || isHit(r)); });
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

function init() {
  loadStore(); me();
  $('#player').value = GAME.current;
  $('#player').onchange = () => { GAME.current = $('#player').value.trim() || 'Player 1'; me(); saveStore(); renderHUD(); renderBadges(); renderDex(); newSniperRound(); };
  $('#contrastSel').onchange = () => { UI.contrast = $('#contrastSel').value; UI.highlight = new Set(); if (UI.hunt) { UI.game = 'explore'; UI.hunt = null; renderHuntSetup(); } renderVolcano(); renderSearch(); };
  $('#fdrSlider').oninput = () => { UI.fdr = +$('#fdrSlider').value; syncSliders(); renderVolcano(); renderSearch(); };
  $('#lfcSlider').oninput = () => { UI.lfc = +$('#lfcSlider').value; syncSliders(); renderVolcano(); renderSearch(); };
  $('#showLod').onchange = () => { UI.showLod = $('#showLod').checked; renderVolcano(); };
  $('#resetCut').onclick = () => { UI.fdr = 0.1; UI.lfc = 1; syncSliders(); renderVolcano(); renderSearch(); };
  $('#exportBtn').onclick = exportHits;
  $('#search').oninput = renderSearch; $('#onlyHits').onchange = renderSearch;
  $('#modal').onclick = ev => { if (ev.target.id === 'modal') closeModal(); };
  document.querySelectorAll('.tabs button').forEach(b => b.onclick = () => {
    document.querySelectorAll('.tabs button').forEach(x => x.classList.toggle('active', x === b));
    document.querySelectorAll('.tabpane').forEach(p => p.hidden = p.id !== b.dataset.tab);
  });
  $('#folder').onchange = async ev => {
    const files = ev.target.files; if (!files.length) return;
    const ok = await loadCourseFolder(files, $('#datasetSel').value);
    if (ok) { award('real'); refreshAll(); }
  };
  $('#datasetSel').onchange = async () => {
    if (!DATA.files) { $('#loadStatus').textContent = 'Choose your course data folder first.'; return; }
    if (await loadCourseFolder(DATA.files, $('#datasetSel').value)) refreshAll();
  };
  $('#demoBtn').onclick = () => { buildDemo(); refreshAll(); $('#loadStatus').textContent = 'Demo data: simulated values.'; };
  buildDemo(); syncSliders(); refreshAll();
}
window.addEventListener('DOMContentLoaded', init);
