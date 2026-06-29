import type { ReactNode } from 'react';

export function Card({ title, children, className = '' }: { title?: string; children: ReactNode; className?: string }) {
  return (
    <section className={`card ${className}`}>
      {title ? <div className="card-title">{title}</div> : null}
      {children}
    </section>
  );
}

export function MetricCard({ label, value, detail, tone = 'blue' }: { label: string; value: string | number; detail?: string; tone?: string }) {
  return (
    <section className="metric-card">
      <div className="metric-label">{label}</div>
      <div className={`metric-value ${tone}`}>{value}</div>
      {detail ? <div className="metric-detail">{detail}</div> : null}
    </section>
  );
}
