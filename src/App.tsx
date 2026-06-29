import { Navigate, Route, Routes, useLocation } from 'react-router';
import { Sidebar } from './components/Sidebar';
import { Topbar } from './components/Topbar';
import { Dashboard } from './pages/Dashboard';
import { NewAnalysis } from './pages/NewAnalysis';
import { SequenceReview } from './pages/SequenceReview';
import { ProteinAnnotation } from './pages/ProteinAnnotation';
import { StructurePockets } from './pages/StructurePockets';
import { CompoundRanking } from './pages/CompoundRanking';
import { EvidenceReview } from './pages/EvidenceReview';
import { PipelineStatus } from './pages/PipelineStatus';
import { SystemStatus } from './pages/SystemStatus';
import { Settings } from './pages/Settings';
import { pageMeta, pagePaths, type PageKey } from './routing';


export default function App() {
  const location = useLocation();
  const matchingPage = (
  Object.entries(pagePaths) as [PageKey, string][]
).find(([, path]) => path === location.pathname);
  const activePage = matchingPage?.[0] ?? 'dashboard';
  const meta = pageMeta[activePage];

  return (
    <div className="app-shell">
      <Sidebar />
      <main className="main-panel">
        <Topbar title={meta.title} subtitle={meta.subtitle} />
        <Routes>
          <Route path={pagePaths.dashboard} element={<Dashboard />} />
          <Route path={pagePaths['new-analysis']} element={<NewAnalysis />} />
          <Route path={pagePaths['sequence-review']} element={<SequenceReview />} />
          <Route path={pagePaths['protein-annotation']} element={<ProteinAnnotation />} />
          <Route path={pagePaths['structure-pockets']} element={<StructurePockets />} />
          <Route path={pagePaths['compound-ranking']} element={<CompoundRanking />} />
          <Route path={pagePaths['evidence-review']} element={<EvidenceReview />} />
          <Route path={pagePaths['pipeline-status']} element={<PipelineStatus />} />
          <Route path={pagePaths['system-status']} element={<SystemStatus />} />
          <Route path={pagePaths.settings} element={<Settings />} />
          <Route path="*" element={<Navigate to={pagePaths.dashboard} replace />} />
        </Routes>  
        <Routes />
      </main>
    </div>
  );
}
