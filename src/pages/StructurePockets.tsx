import { Card } from '../components/Card';
import { StatusBadge } from '../components/StatusBadge';
import { pocketRows, proteinsForViewer } from '../data/demoData';

export function StructurePockets() {
  return (
    <div className="page structure-page">
      <div className="structure-shell">
        <Card title="Proteins" className="protein-picker">
          {proteinsForViewer.map((protein, idx) => (
            <div className={`protein-pick ${idx === 0 ? 'selected' : ''}`} key={protein.id}>
              <strong>{protein.id}</strong>
              <span>{protein.name}</span>
              <div>
                <StatusBadge status="complete" label="Structure ✓" />
                <StatusBadge status="queued" label={`${protein.pockets} pockets`} />
              </div>
            </div>
          ))}
        </Card>

        <div className="viewer-column">
          <div className="viewer-toolbar">VP40_EXRC031 · AlphaFold2 · pLDDT avg 84.1</div>
          <div className="molecule-viewer">
            <div className="protein-surface" />
            <div className="orbit o1" />
            <div className="orbit o2" />
            <div className="orbit o3" />
            <div className="orbit o4" />
            <div className="orbit o5" />
            <div className="pocket p1"><span>PKT-001</span></div>
            <div className="pocket p2"><span>PKT-002</span></div>
            <div className="demo-chip">DEMO DATA — Computational model only</div>
          </div>
        </div>
      </div>

      <div className="pocket-grid">
        {pocketRows.map((pocket) => (
          <Card key={pocket.id}>
            <div className="pocket-card-head">
              <strong>{pocket.id}</strong>
              <StatusBadge status={pocket.priority.toLowerCase()} label={pocket.priority} />
            </div>
            <div className="pocket-fields">
              <div><span>Druggability</span><strong>{pocket.druggability}</strong></div>
              <div><span>Volume</span><strong>{pocket.volume}</strong></div>
              <div><span>Residues</span><strong>{pocket.residues}</strong></div>
              <div><span>Domain</span><strong>{pocket.domain}</strong></div>
            </div>
          </Card>
        ))}
      </div>
    </div>
  );
}
