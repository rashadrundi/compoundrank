import { Card } from '../components/Card';
import { StatusBadge } from '../components/StatusBadge';
import { ProgressBar } from '../components/ProgressBar';
import { proteinRows } from '../data/demoData';

export function ProteinAnnotation() {
  return (
    <div className="page">
      <div className="filter-row">
        {['All', 'High Confidence', 'High Priority', 'Unknown / Unannotated'].map((filter, index) => (
          <button className={`filter-button ${index === 0 ? 'active' : ''}`} key={filter}>{filter}</button>
        ))}
      </div>

      <Card className="table-card">
        <table className="data-table">
          <thead>
            <tr>
              <th>Protein ID</th>
              <th>Len</th>
              <th>CDD Hit</th>
              <th>InterPro Family</th>
              <th>VOGDB</th>
              <th>Function</th>
              <th>Confidence</th>
              <th>Priority</th>
            </tr>
          </thead>
          <tbody>
            {proteinRows.map((row) => (
              <tr key={row.proteinId}>
                <td><a>{row.proteinId}</a></td>
                <td>{row.len}</td>
                <td>{row.cdd}</td>
                <td>{row.interpro}</td>
                <td>{row.vog}</td>
                <td>{row.function}</td>
                <td className="confidence-cell"><ProgressBar value={row.confidence} tone={row.confidence > 80 ? 'green' : row.confidence > 60 ? 'blue' : 'muted'} /> <span>{row.confidence}%</span></td>
                <td><StatusBadge status={row.priority.toLowerCase()} label={row.priority} /></td>
              </tr>
            ))}
          </tbody>
        </table>
      </Card>
    </div>
  );
}
