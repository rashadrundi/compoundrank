import type { Status } from '../data/demoData';

export function StatusBadge({ status, label }: { status: Status | string; label?: string }) {
  const normalized = String(status).toLowerCase();
  return <span className={`status-badge ${normalized}`}>{label ?? normalized}</span>;
}

export function EvidenceBadge({ label }: { label: string }) {
  const key = label.toLowerCase().replaceAll(' ', '-');
  return <span className={`evidence-badge ${key}`}>{label}</span>;
}
