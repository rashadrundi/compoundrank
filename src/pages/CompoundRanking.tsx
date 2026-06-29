import { Card } from '../components/Card';
import { EvidenceBadge } from '../components/StatusBadge';
import { ProgressBar } from '../components/ProgressBar';
import { compoundRows } from '../data/demoData';

export function CompoundRanking() {
  return (
    <div className="page">
      <div className="filter-row">
        <span className="filter-label">⌁ Evidence filter:</span>
        {['Strong Evidence', 'Moderate Evidence', 'Weak Evidence', 'Computational Only'].map((filter) => (
          <button className={`filter-button ${filter.toLowerCase().includes('strong') ? 'active green' : ''}`} key={filter}>{filter}</button>
        ))}
      </div>

      <Card className="table-card">
        <table className="data-table compound-table">
          <thead>
            <tr>
              <th>#</th>
              <th>Compound</th>
              <th>DB</th>
              <th>Target</th>
              <th>Docking</th>
              <th>Rescore</th>
              <th>Literature</th>
              <th>Safety</th>
              <th>Resistance</th>
              <th>Score</th>
            </tr>
          </thead>
          <tbody>
            {compoundRows.map((row) => (
              <tr key={row.rank}>
                <td>#{row.rank}</td>
                <td><strong>{row.compound}</strong></td>
                <td>{row.db}</td>
                <td><a>{row.target}</a></td>
                <td>{row.docking}</td>
                <td>{row.rescore}</td>
                <td><EvidenceBadge label={row.literature} /></td>
                <td className={row.safety === 'High' ? 'danger-text' : row.safety === 'Moderate' ? 'amber-text' : 'green-text'}>{row.safety}</td>
                <td className={row.resistance === 'Moderate' ? 'amber-text' : row.resistance === 'Unknown' ? 'muted-text' : 'green-text'}>{row.resistance}</td>
                <td className="score-cell"><ProgressBar value={row.score} /> <span>{row.score}</span></td>
              </tr>
            ))}
          </tbody>
        </table>
      </Card>
    </div>
  );
}
