import {
  useMutation,
  useQuery,
  useQueryClient,
} from '@tanstack/react-query';

import {
  cancelRun,
  createRun,
  getRun,
  listRuns,
} from './runs';

import type {
  AnalysisRun,
  RunStatus,
} from './types';

export const runQueryKeys = {
  all: ['runs'] as const,

  lists: () =>
    [...runQueryKeys.all, 'list'] as const,

  detail: (runId: string) =>
    [...runQueryKeys.all, 'detail', runId] as const,
};

const finishedRunStatuses: RunStatus[] = [
  'completed',
  'partial',
  'failed',
  'cancelled',
];

function isRunFinished(status?: RunStatus): boolean {
  if (!status) {
    return false;
  }

  return finishedRunStatuses.includes(status);
}

export function useRuns() {
  return useQuery({
    queryKey: runQueryKeys.lists(),
    queryFn: listRuns,
    staleTime: 5_000,
  });
}

export function useRun(runId: string) {
  return useQuery<AnalysisRun>({
    queryKey: runQueryKeys.detail(runId),

    queryFn: () => getRun(runId),

    enabled: Boolean(runId),

    refetchInterval: (query) => {
      const run = query.state.data;

      if (isRunFinished(run?.status)) {
        return false;
      }

      return 3_000;
    },
  });
}

export function useCreateRun() {
  const queryClient = useQueryClient();

  return useMutation({
    mutationFn: createRun,

    onSuccess: async (createdRun) => {
      queryClient.setQueryData(
        runQueryKeys.detail(createdRun.run_id),
        createdRun,
      );

      await queryClient.invalidateQueries({
        queryKey: runQueryKeys.lists(),
      });
    },
  });
}

export function useCancelRun() {
  const queryClient = useQueryClient();

  return useMutation({
    mutationFn: cancelRun,

    onSuccess: async (cancelledRun) => {
      queryClient.setQueryData(
        runQueryKeys.detail(cancelledRun.run_id),
        cancelledRun,
      );

      await queryClient.invalidateQueries({
        queryKey: runQueryKeys.lists(),
      });
    },
  });
}