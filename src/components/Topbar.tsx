export function Topbar({ title, subtitle }: { title: string; subtitle?: string }) {
  return (
    <header className="topbar">
      <div>
        <h1>{title}</h1>
        {subtitle ? <p>{subtitle}</p> : null}
      </div>
      <div className="topbar-actions">
        <span className="system-dot" />
        <button type="button" className="icon-button" aria-label="Refresh">↻</button>
      </div>
    </header>
  );
}
