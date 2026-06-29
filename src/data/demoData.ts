export type Status = 'complete' | 'running' | 'queued' | 'warning' | 'failed' | 'idle';

export const analysisJobs = [
  { id: 'EXR-2024-031', title: 'Bat CoV BO-11 Panel', detail: 'Pocket Detection · 2024-03-15 09:42', status: 'running' as Status, progress: 67 },
  { id: 'EXR-2024-030', title: 'Murid arenavirus screen', detail: 'Complete · 2024-03-14 14:20', status: 'complete' as Status, progress: 100 },
  { id: 'EXR-2024-029', title: 'Unknown RNA-seq ETV-47', detail: 'Complete · 2024-03-13 08:05', status: 'complete' as Status, progress: 100 },
  { id: 'EXR-2024-028', title: 'Paramyxo candidate screen', detail: 'InterPro annotation · 2024-03-12 16:30', status: 'warning' as Status, progress: 43 },
  { id: 'EXR-2024-027', title: 'Tick-borne virus panel', detail: 'CDD annotation · 2024-03-11 11:00', status: 'failed' as Status, progress: 18 }
];

export const evidenceSummary = [
  { label: 'High Confidence ≥90%', count: 24, tone: 'green' },
  { label: 'Moderate 75–89%', count: 41, tone: 'blue' },
  { label: 'Low 50–74%', count: 57, tone: 'amber' },
  { label: 'Insufficient <50%', count: 18, tone: 'muted' }
];

export const topTargets = [
  { id: 'VP40_EXRC031', score: 97 },
  { id: 'RdRP_EXRC031', score: 94 },
  { id: 'NP_EXRC031', score: 88 },
  { id: 'PROT_EXRC031a', score: 81 },
  { id: 'ENV_EXRC031', score: 76 }
];

export const proteinRows = [
  { proteinId: 'VP40_EXRC031', len: 326, cdd: 'Ebola_VP40', interpro: 'Filovirus matrix', vog: 'VOG0814', function: 'Matrix protein', confidence: 97, priority: 'High' },
  { proteinId: 'RdRP_EXRC031', len: 2142, cdd: 'RdRP_corona', interpro: 'RNA-dep RNA pol', vog: 'VOG0021', function: 'RNA polymerase', confidence: 94, priority: 'High' },
  { proteinId: 'NP_EXRC031', len: 724, cdd: 'Nucprot_param', interpro: 'Nucleoprotein', vog: 'VOG0192', function: 'Nucleoprotein', confidence: 88, priority: 'Medium' },
  { proteinId: 'PROT_EXRC031a', len: 459, cdd: 'Pfam-chymotryp', interpro: 'Serine protease', vog: '—', function: 'Viral protease', confidence: 81, priority: 'Medium' },
  { proteinId: 'ENV_EXRC031', len: 512, cdd: 'Pfam-fusogen', interpro: 'Type-I fusion', vog: 'VOG0311', function: 'Envelope glycoprotein', confidence: 76, priority: 'High' },
  { proteinId: 'UNK_EXRC031b', len: 187, cdd: '—', interpro: 'DUF4699', vog: 'VOG2104', function: 'Unknown', confidence: 34, priority: 'Low' },
  { proteinId: 'UNK_EXRC031c', len: 93, cdd: '—', interpro: '—', vog: '—', function: 'Unannotated', confidence: 12, priority: 'Low' }
];

export const proteinsForViewer = [
  { id: 'VP40_EXRC031', name: 'Matrix protein', structure: true, pockets: 3 },
  { id: 'RdRP_EXRC031', name: 'RNA polymerase', structure: true, pockets: 3 },
  { id: 'NP_EXRC031', name: 'Nucleoprotein', structure: true, pockets: 3 },
  { id: 'PROT_EXRC031a', name: 'Viral protease', structure: true, pockets: 3 },
  { id: 'ENV_EXRC031', name: 'Envelope glycoprotein', structure: true, pockets: 3 }
];

export const pocketRows = [
  { id: 'PKT-001', druggability: 'High', volume: '842 Å³', residues: 24, domain: 'Active site (protease)', priority: 'High' },
  { id: 'PKT-002', druggability: 'Moderate', volume: '521 Å³', residues: 18, domain: 'Allosteric site', priority: 'Medium' },
  { id: 'PKT-003', druggability: 'Low', volume: '294 Å³', residues: 11, domain: 'Surface patch', priority: 'Low' }
];

export const compoundRows = [
  { rank: 1, compound: 'Remdesivir', db: 'ChEMBL', target: 'RdRP_EXRC031', docking: -9.4, rescore: -10.1, literature: 'Strong Evidence', safety: 'Low', resistance: 'Low', score: 94 },
  { rank: 2, compound: 'GS-441524', db: 'PubChem', target: 'RdRP_EXRC031', docking: -8.9, rescore: -9.6, literature: 'Moderate Evidence', safety: 'Low', resistance: 'Low', score: 88 },
  { rank: 3, compound: 'Favipiravir', db: 'ChEMBL', target: 'RdRP_EXRC031', docking: -8.1, rescore: -8.7, literature: 'Moderate Evidence', safety: 'Low', resistance: 'Moderate', score: 79 },
  { rank: 4, compound: 'Camostat mesylate', db: 'DrugBank', target: 'PROT_EXRC031a', docking: -7.8, rescore: -8.2, literature: 'Weak Evidence', safety: 'Low', resistance: 'Low', score: 71 },
  { rank: 5, compound: 'E-64d', db: 'PubChem', target: 'PROT_EXRC031a', docking: -7.3, rescore: -7.9, literature: 'Weak Evidence', safety: 'Moderate', resistance: 'Low', score: 64 },
  { rank: 6, compound: 'Compound X-887', db: 'ZINC', target: 'ENV_EXRC031', docking: -7.1, rescore: -7.4, literature: 'Computational Only', safety: 'Unknown', resistance: 'Unknown', score: 58 },
  { rank: 7, compound: 'BMS-806 analog', db: 'ZINC', target: 'ENV_EXRC031', docking: -6.8, rescore: -7.0, literature: 'Computational Only', safety: 'Unknown', resistance: 'Unknown', score: 52 },
  { rank: 8, compound: 'MG-132', db: 'ChEMBL', target: 'PROT_EXRC031a', docking: -6.4, rescore: -6.9, literature: 'Insufficient Evidence', safety: 'High', resistance: 'Low', score: 38 }
];

export const pipelineSteps = [
  { step: 1, name: 'FASTA Parsing', note: '248 contigs parsed', status: 'complete' as Status, time: '0m 04s' },
  { step: 2, name: 'CDD Annotation', note: '142 domain hits', status: 'complete' as Status, time: '3m 12s' },
  { step: 3, name: 'InterPro Annotation', note: '89 family assignments', status: 'complete' as Status, time: '7m 48s' },
  { step: 4, name: 'VOGDB Annotation', note: '56 viral ORF matches', status: 'complete' as Status, time: '2m 33s' },
  { step: 5, name: 'Structure Prediction', note: 'ColabFold · 31/37 models', status: 'running' as Status, time: '14m 21s+' },
  { step: 6, name: 'Pocket Detection', note: 'Waiting on structures', status: 'queued' as Status, time: '—' },
  { step: 7, name: 'Compound DB Retrieval', note: 'Waiting on pockets', status: 'queued' as Status, time: '—' },
  { step: 8, name: 'Docking & Rescoring', note: 'GNINA scheduled', status: 'queued' as Status, time: '—' },
  { step: 9, name: 'Evidence Grading', note: 'Auto-grade pipeline', status: 'queued' as Status, time: '—' },
  { step: 10, name: 'Report Generation', note: 'Pending all stages', status: 'queued' as Status, time: '—' }
];

export const safeLogs = [
  '2024-03-15 09:42:01 [INFO] Job EXR-2024-031 initialized',
  '2024-03-15 09:42:03 [INFO] FASTA validated: 248 contigs, 3.7 MB',
  '2024-03-15 09:42:04 [INFO] CDD annotation — dispatching batch',
  '2024-03-15 09:45:16 [INFO] CDD annotation — 142 hits returned',
  '2024-03-15 09:45:17 [INFO] InterPro annotation — dispatching batch',
  '2024-03-15 09:53:06 [INFO] InterPro annotation — 89 families assigned',
  '2024-03-15 09:53:06 [INFO] VOGDB annotation — dispatching batch',
  '2024-03-15 09:55:39 [INFO] VOGDB annotation — 56 ORF matches',
  '2024-03-15 09:55:40 [INFO] Structure prediction — ColabFold batch submitted',
  '2024-03-15 10:10:01 [INFO] Structure prediction — 31/37 models complete',
  '2024-03-15 10:10:01 [INFO] Pocket detection — queued, awaiting structures'
];

export const systemCards = [
  { name: 'CPU Server', status: 'complete' as Status, detail: '12-core / 64 GB RAM' },
  { name: 'GPU Worker', status: 'running' as Status, detail: 'RTX 4090 · busy' },
  { name: 'CDD (NCBI)', status: 'complete' as Status, detail: 'v3.21 connected' },
  { name: 'InterPro', status: 'complete' as Status, detail: 'EBI API v5' },
  { name: 'VOGDB', status: 'warning' as Status, detail: 'Rate limit: 80%' },
  { name: 'fpocket', status: 'complete' as Status, detail: 'v4.1 local' },
  { name: 'ColabFold', status: 'running' as Status, detail: 'Batch running' },
  { name: 'GNINA', status: 'queued' as Status, detail: 'Queued' },
  { name: 'Result Storage', status: 'complete' as Status, detail: '2.1 TB / 4 TB' },
  { name: 'Report Generator', status: 'complete' as Status, detail: 'v2.3.1 ready' }
];
