import { Card } from '../components/Card';
import { StatusBadge } from '../components/StatusBadge';
import { safeLogs, systemCards } from '../data/demoData';

export function SystemStatus() {
  return (
    <div className="page system-page">
      <div className="system-card-grid">
        {systemCards.map((system) => (
          <Card key={system.name}>
            <div className="system-card-head">
              <strong>{system.name}</strong>
              <StatusBadge status={system.status} />
            </div>
            <span className="system-detail">{system.detail}</span>
          </Card>
        ))}
      </div>

      <div className="two-column status-bottom">
        <Card title="Configuration">
          <div className="config-form">
            <label><span>API Endpoint</span><input value="https://api.exorcist.internal/v2" readOnly /></label>
            <label><span>Job Directory</span><input value="/mnt/data/exorcist/jobs" readOnly /></label>
            <label><span>Result Directory</span><input value="/mnt/data/exorcist/results" readOnly /></label>
            <label><span>Worker Status</span><input value="Active — 2 workers" readOnly /></label>
            <label><span>Last Successful Run</span><input value="2024-03-14 16:02:41" readOnly /></label>
          </div>
        </Card>

        <Card title="System Log" className="log-card">
          <pre>{safeLogs.join('\n')}</pre>
        </Card>
      </div>
    </div>
  );
}
