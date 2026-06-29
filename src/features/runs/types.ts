export type RunStatus =
  | 'queued'
  | 'running'
  | 'completed'
  | 'partial'
  | 'failed'
  | 'cancelled';

export type StageStatus =
  | 'queued'
  | 'running'
  | 'completed'
  | 'warning'
  | 'failed'
  | 'skipped';

export type StageKey =
  | 'fasta_parsing'
  | 'cdd'
  | 'interpro'
  | 'vogdb'
  | 'structure_prediction'
  | 'pocket_detection'
  | 'compound_retrieval'
  | 'docking_rescoring'
  | 'evidence_grading'
  | 'report_generation';

export interface PipelineStage {
  key: StageKey;
  label: string;
  status: StageStatus;
  progress: number;
  message?: string;
  started_at?: string;
  completed_at?: string;
}

export interface AnalysisRun {
  run_id: string;
  project_name: string;
  sample_label: string;
  organism_context?: string;
  suspected_virus?: string;
  sequencing_method?: string;
  metadata_notes?: string;
  fasta_filename: string;
  status: RunStatus;
  current_stage: StageKey | null;
  progress: number;
  created_at: string;
  updated_at: string;
  stages: PipelineStage[];
}