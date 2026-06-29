import { apiRequest } from './client';
import type { AnalysisRun } from './types';

export interface CreateRunInput {
  projectName: string;
  sampleLabel: string;
  organismContext?: string;
  suspectedVirus?: string;
  sequencingMethod?: string;
  metadataNotes?: string;
  fastaFile: File;
}

function appendOptionalField(
  formData: FormData,
  fieldName: string,
  value?: string,
): void {
  const trimmedValue = value?.trim();

  if (trimmedValue) {
    formData.append(fieldName, trimmedValue);
  }
}

export async function createRun(
  input: CreateRunInput,
): Promise<AnalysisRun> {
  const formData = new FormData();

  formData.append(
    'project_name',
    input.projectName.trim(),
  );

  formData.append(
    'sample_label',
    input.sampleLabel.trim(),
  );

  appendOptionalField(
    formData,
    'organism_context',
    input.organismContext,
  );

  appendOptionalField(
    formData,
    'suspected_virus',
    input.suspectedVirus,
  );

  appendOptionalField(
    formData,
    'sequencing_method',
    input.sequencingMethod,
  );

  appendOptionalField(
    formData,
    'metadata_notes',
    input.metadataNotes,
  );

  formData.append(
    'fasta_file',
    input.fastaFile,
    input.fastaFile.name
  );

  return apiRequest<AnalysisRun>('/runs', {
    method: 'POST',
    body: formData,
  });
}

export async function getRun(
  runId: string,
): Promise<AnalysisRun> {
  const encodedRunId = encodeURIComponent(runId);

  return apiRequest<AnalysisRun>(
    `/runs/${encodedRunId}`,
  );
}

export async function listRuns(): Promise<AnalysisRun[]> {
  return apiRequest<AnalysisRun[]>('/runs');
}

export async function cancelRun(
  runId: string,
): Promise<AnalysisRun> {
  const encodedRunId = encodeURIComponent(runId);

  return apiRequest<AnalysisRun>(
    `/runs/${encodedRunId}/cancel`,
    {
      method: 'POST',
    },
  );
}