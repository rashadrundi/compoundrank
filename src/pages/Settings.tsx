import { Card } from '../components/Card';

export function Settings() {
  return (
    <div className="page settings-page">
      <Card title="User Preferences" className="settings-card">
        <div className="config-form">
          <label><span>Default Analysis Preset</span><select defaultValue="standard"><option value="standard">Standard</option><option value="deep">Deep Review</option><option value="quick">Quick Triage</option></select></label>
          <label><span>Report Author Name</span><input placeholder="Researcher name" /></label>
          <label><span>Institution</span><input placeholder="Research institute" /></label>
        </div>
      </Card>

      <Card title="Platform Versions" className="settings-card">
        <div className="version-list">
          <div><span>EXORCIST Platform</span><strong>v2.4.1</strong></div>
          <div><span>CDD Module</span><strong>v3.21</strong></div>
          <div><span>ColabFold</span><strong>v1.5.5</strong></div>
          <div><span>GNINA</span><strong>v1.1</strong></div>
          <div><span>fpocket</span><strong>v4.1</strong></div>
        </div>
      </Card>

      <button className="primary-small">Save Settings</button>
    </div>
  );
}
