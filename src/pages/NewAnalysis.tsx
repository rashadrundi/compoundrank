import { Card } from '../components/Card';
import { SafeNotice } from '../components/SafeNotice';

export function NewAnalysis() {
  return (
    <div className="page centered-page">
      <SafeNotice>
        Research-use computational prioritization only. Not for clinical diagnosis, treatment, pathogen engineering, wet-lab execution, or clinical decision-making.
      </SafeNotice>

      <Card title="Run Configuration" className="form-card">
        <div className="form-grid two">
          <label>
            <span>Project Name</span>
            <input placeholder="e.g. Bat CoV Panel 2024" />
          </label>
          <label>
            <span>Sample Label</span>
            <input placeholder="e.g. BO-11-A-2024" />
          </label>
          <label>
            <span>Organism / Context</span>
            <input placeholder="e.g. Rhinolophus affinis" />
          </label>
          <label>
            <span>Suspected Virus</span>
            <input placeholder="e.g. Novel betacoronavirus" />
          </label>
        </div>
        <label className="full-label">
          <span>Sequencing Method</span>
          <select defaultValue="illumina">
            <option value="illumina">Illumina paired-end (150 bp)</option>
            <option value="nanopore">Oxford Nanopore long-read</option>
            <option value="pacbio">PacBio HiFi</option>
            <option value="unknown">Unknown / not provided</option>
          </select>
        </label>
      </Card>

      <Card title="FASTA Input" className="form-card">
        <label className="dropzone">
          <input type="file" accept=".fasta,.fa,.faa,.fna" />
          <div className="drop-icon">⇧</div>
          <strong>Drag and drop a FASTA file, or click to browse</strong>
          <span>Accepts .fa · .fasta · .fna — max 500 MB</span>
        </label>
      </Card>

      <Card title="Metadata Notes (Optional)" className="form-card">
        <textarea placeholder="Collection date, geographic region, sequencing facility, study context, notes..." />
      </Card>

      <button type="button" className="primary-action">▷ Start EXORCIST Analysis</button>
    </div>
  );
}
