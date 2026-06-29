import { Card } from '../components/Card';
import { StatusBadge } from '../components/StatusBadge';
import { pipelineSteps, safeLogs } from '../data/demoData';

export function PipelineStatus() {
  return (
    <div className="page pipeline-page">
      <div className="job-toolbar">
        <span>Job:</span>
        <select defaultValue="EXR-2024-031">
          <option>EXR-2024-031 — Bat CoV BO-11 Panel</option>
          <option>EXR-2024-030 — Murid arenavirus screen</option>
        </select>
        <StatusBadge status="running" label="Running — 67% complete" />
      </div>

      <div className="pipeline-layout">
        <Card title="Pipeline Steps">
          <div className="step-list">
            {pipelineSteps.map((step) => (
              <div className={`pipeline-step ${step.status}`} key={step.step}>
                <div className="step-number">{step.step}</div>
                <div className="step-text">
                  <strong>{step.name}</strong>
                  <span>{step.note}</span>
                </div>
                <StatusBadge status={step.status} />
                <div className="step-time">{step.time}</div>
              </div>
            ))}
          </div>
        </Card>

        <Card title="System Log" className="log-card">
          <pre>{safeLogs.join('\n')}\n↳ Structure prediction in progress...</pre>
        </Card>
      </div>
    </div>
  );
}
