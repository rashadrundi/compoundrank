import { Card, MetricCard } from '../components/Card';
import { SafeNotice } from '../components/SafeNotice';

const reportSections = [
  'Executive Summary',
  'Input Sequence Summary',
  'Computational Classification',
  'Protein Target Prioritization',
  'Structure & Pocket Findings',
  'Compound Ranking Summary',
  'Evidence Limitations',
  'Recommended Review Pathway'
];

export function EvidenceReview() {
  return (
    <div className="page evidence-page">
      <div className="report-shell">
        <Card title="Report Sections" className="report-nav">
          {reportSections.map((section, idx) => (
            <button className={idx === 0 ? 'active' : ''} key={section}>0{idx + 1} {section}</button>
          ))}
        </Card>

        <div className="report-main">
          <SafeNotice>This report is for research review only and does not provide clinical, diagnostic, therapeutic, or experimental instructions.</SafeNotice>
          <Card className="report-document">
            <div className="report-kicker">Section 1 of 8 — EXR-2024-031 Demo Data</div>
            <h2>Executive Summary</h2>
            <p>
              This report summarizes the computational analysis of run EXR-2024-031, demo dataset. Seven protein targets were identified with high annotation confidence. Eight compounds from existing databases rank as candidates for expert review. All findings require qualified researcher evaluation before any further action.
            </p>
            <div className="metric-grid three compact">
              <MetricCard label="Run ID" value="EXR-2024-031" />
              <MetricCard label="Protein Targets" value="7" />
              <MetricCard label="Top Compounds" value="8" tone="violet" />
            </div>
            <div className="button-row">
              <button className="secondary-action">Export PDF</button>
              <button className="secondary-action">Export CSV</button>
              <button className="secondary-action">Generate Research Handoff Report</button>
              <button className="primary-small">Save Draft Report</button>
            </div>
          </Card>
        </div>
      </div>
    </div>
  );
}
