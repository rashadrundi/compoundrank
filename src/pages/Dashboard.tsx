import { Card, MetricCard } from '../components/Card';
import { ProgressBar } from '../components/ProgressBar';
import { StatusBadge } from '../components/StatusBadge';
import { analysisJobs, evidenceSummary, topTargets } from '../data/demoData';

export function Dashboard() {
  return (
    <div className="page page-dashboard">
      <div className="metric-grid four">
        <MetricCard label="Active Jobs" value="2" detail="Currently running" tone="blue" />
        <MetricCard label="Completed" value="14" detail="Last 30 days" tone="green" />
        <MetricCard label="Warnings / Failed" value="3" detail="Requires attention" tone="amber" />
        <MetricCard label="High-Priority Targets" value="12" detail="From 47 annotated" tone="violet" />
      </div>

      <div className="dashboard-layout">
        <Card title="Analysis Jobs" className="large-card">
          <div className="job-list">
            {analysisJobs.map((job) => (
              <div className="job-row" key={job.id}>
                <div className="job-copy">
                  <div className="job-id">{job.id} <StatusBadge status={job.status} /></div>
                  <div className="job-title">{job.title}</div>
                  <div className="job-detail">{job.detail}</div>
                </div>
                <div className="job-progress">
                  <span>{job.progress}%</span>
                  <ProgressBar value={job.progress} tone={job.status === 'failed' ? 'red' : job.status === 'warning' ? 'amber' : 'blue'} />
                </div>
              </div>
            ))}
          </div>
        </Card>

        <div className="stack">
          <Card title="Evidence Summary">
            <div className="evidence-list">
              {evidenceSummary.map((item) => (
                <div className="evidence-row" key={item.label}>
                  <span className={`dot ${item.tone}`} />
                  <span>{item.label}</span>
                  <strong>{item.count}</strong>
                </div>
              ))}
            </div>
          </Card>

          <Card title="Top Protein Targets">
            <div className="target-list">
              {topTargets.map((target) => (
                <div className="target-row" key={target.id}>
                  <span>{target.id}</span>
                  <strong>{target.score}%</strong>
                </div>
              ))}
            </div>
          </Card>

          <Card title="Recent Exports">
            <div className="export-list">
              <div>▧ EXR-2024-030 Full Report <span>Mar 14</span></div>
              <div>▧ Protein targets EXR-030 <span>Mar 14</span></div>
              <div>▧ EXR-2024-029 Full Report <span>Mar 13</span></div>
            </div>
          </Card>
        </div>
      </div>
    </div>
  );
}
