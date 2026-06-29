export type PageKey =
  | 'dashboard'
  | 'new-analysis'
  | 'sequence-review'
  | 'protein-annotation'
  | 'structure-pockets'
  | 'compound-ranking'
  | 'evidence-review'
  | 'pipeline-status'
  | 'system-status'
  | 'settings';

export const pagePaths: Record<PageKey, string> = {
  dashboard: '/',
  'new-analysis': '/runs/new',
  'sequence-review': '/sequence-review',
  'protein-annotation': '/protein-annotation',
  'structure-pockets': '/structure-pockets',
  'compound-ranking': '/compound-ranking',
  'evidence-review': '/evidence-review',
  'pipeline-status': '/pipeline-status',
  'system-status': '/system-status',
  settings: '/settings'
};

export const pageMeta: Record<PageKey, { title: string; subtitle: string }> = {
  dashboard: { title: 'Research Dashboard', subtitle: 'EXORCIST Computational Prioritization Platform' },
  'new-analysis': { title: 'New Analysis Run', subtitle: 'Configure and launch a computational screening job' },
  'sequence-review': { title: 'Sequence Review', subtitle: 'EXR-2024-031 · Bat CoV BO-11 Panel · DEMO DATA' },
  'protein-annotation': { title: 'Protein Annotation Results', subtitle: 'EXR-2024-031 · 7 proteins annotated · DEMO DATA' },
  'structure-pockets': { title: 'Structure & Pocket Analysis', subtitle: 'EXR-2024-031 · AlphaFold2 models · DEMO DATA' },
  'compound-ranking': { title: 'Compound Ranking', subtitle: 'EXR-2024-031 · Candidates for expert review · DEMO DATA' },
  'evidence-review': { title: 'Evidence Review & Reports', subtitle: 'EXR-2024-031 · Research-use only · DEMO DATA' },
  'pipeline-status': { title: 'Pipeline Status', subtitle: 'Real-time job tracking and system log' },
  'system-status': { title: 'System Status', subtitle: 'API connections, workers, and configuration' },
  settings: { title: 'Settings', subtitle: 'Platform preferences and version information' }
};