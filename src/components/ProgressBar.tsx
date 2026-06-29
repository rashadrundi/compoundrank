export function ProgressBar({ value, tone = 'blue' }: { value: number; tone?: string }) {
  return (
    <div className="progress-wrap" aria-label={`Progress ${value}%`}>
      <div className={`progress-fill ${tone}`} style={{ width: `${Math.min(Math.max(value, 0), 100)}%` }} />
    </div>
  );
}
