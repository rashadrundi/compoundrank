import { Card, MetricCard } from '../components/Card';

export function SequenceReview() {
  const flow = ['Input Reads / FASTA', 'Quality Summary', 'Candidate Sequence Review', 'Protein Candidates', 'Downstream Annotation'];
  return (
    <div className="page">
      <div className="metric-grid four">
        <MetricCard label="Input Contigs" value="248" detail="From FASTA input" />
        <MetricCard label="Host Filtered" value="31" detail="Background removed" />
        <MetricCard label="Viral Signal" value="187" detail="Candidate contigs" />
        <MetricCard label="Unclassified" value="30" detail="Requires review" tone="amber" />
      </div>

      <Card title="Computational Sequence Flow">
        <div className="flow-strip">
          {flow.map((step, idx) => (
            <div className="flow-node-wrap" key={step}>
              <div className="flow-node">{step}</div>
              {idx < flow.length - 1 ? <span className="flow-arrow">›</span> : null}
            </div>
          ))}
        </div>
      </Card>

      <div className="chart-grid">
        <Card title="Taxonomic Distribution">
          <div className="pie-layout">
            <div className="pie-chart" />
            <div className="legend">
              <div><span className="dot blue" /> Betacoronavirus <strong>34%</strong></div>
              <div><span className="dot violet" /> Filoviridae <strong>27%</strong></div>
              <div><span className="dot cyan" /> Paramyxoviridae <strong>19%</strong></div>
              <div><span className="dot green" /> Reoviridae <strong>11%</strong></div>
              <div><span className="dot muted" /> Unclassified <strong>9%</strong></div>
            </div>
          </div>
        </Card>

        <Card title="Sequence Quality (Q-score)">
          <div className="bar-chart vertical">
            {[12, 38, 86, 142, 96, 44].map((height, idx) => (
              <div className="vbar-wrap" key={idx}>
                <div className="vbar" style={{ height }} />
                <span>Q{idx === 0 ? '10–20' : idx === 1 ? '20–25' : idx === 2 ? '25–30' : idx === 3 ? '30–35' : idx === 4 ? '35–40' : '40+'}</span>
              </div>
            ))}
          </div>
        </Card>

        <Card title="Confidence Distribution">
          <div className="bar-chart horizontal">
            {[
              ['90–100%', 24, 'green'], ['75–89%', 41, 'blue'], ['50–74%', 57, 'amber'], ['25–49%', 39, 'orange'], ['<25%', 18, 'muted']
            ].map(([label, count, tone]) => (
              <div className="hbar-row" key={label as string}>
                <span>{label}</span>
                <div className="hbar-track"><div className={`hbar ${tone}`} style={{ width: `${Number(count)}%` }} /></div>
              </div>
            ))}
          </div>
        </Card>
      </div>

      <div className="two-column">
        <Card title="Known Reference Similarity">
          <div className="reference-list">
            <div><strong>Beta-CoV RdRP homolog</strong><span>Sarbecovirus clade</span><b>91%</b></div>
            <div><strong>Filovirus VP40 homolog</strong><span>Marburg-like matrix</span><b>78%</b></div>
            <div><strong>Paramyxo NP homolog</strong><span>Distant — Respirovirus</span><b>65%</b></div>
            <div><strong>Unknown structural ORF</strong><span>No close reference match</span><b>&lt;30%</b></div>
          </div>
        </Card>

        <Card title="Unclassified Sequence Signal">
          <div className="kv-list">
            <div><span>Total unclassified contigs</span><strong>30</strong></div>
            <div><span>Novel ORF candidates</span><strong className="amber-text">7</strong></div>
            <div><span>Putative structural ORFs</span><strong>3</strong></div>
            <div><span>Low-complexity / repeat</span><strong>11</strong></div>
            <div><span>Host contamination probable</span><strong>9</strong></div>
          </div>
          <div className="inline-note">7 unclassified ORFs flagged for manual review. Expert annotation recommended before target prioritization.</div>
        </Card>
      </div>
    </div>
  );
}
