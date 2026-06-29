function prepareRunId(runId: string): string {
  const trimmedRunId = runId.trim();

  if (!trimmedRunId) {
    throw new Error('A run ID is required to build a run URL.');
  }

  return encodeURIComponent(trimmedRunId);
}

export const runRoutePatterns = {
  pipeline: '/runs/:runId/pipeline',
  sequenceReview: '/runs/:runId/sequences',
  annotations: '/runs/:runId/annotations',
  structures: '/runs/:runId/structures',
  compounds: '/runs/:runId/compounds',
  evidence: '/runs/:runId/evidence',
  report: '/runs/:runId/report',
} as const;

export const runPaths = {
  pipeline: (runId: string) =>
    `/runs/${prepareRunId(runId)}/pipeline`,

  sequenceReview: (runId: string) =>
    `/runs/${prepareRunId(runId)}/sequences`,

  annotations: (runId: string) =>
    `/runs/${prepareRunId(runId)}/annotations`,

  structures: (runId: string) =>
    `/runs/${prepareRunId(runId)}/structures`,

  compounds: (runId: string) =>
    `/runs/${prepareRunId(runId)}/compounds`,

  evidence: (runId: string) =>
    `/runs/${prepareRunId(runId)}/evidence`,

  report: (runId: string) =>
    `/runs/${prepareRunId(runId)}/report`,
} as const;