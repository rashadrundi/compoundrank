import { pagePaths, type PageKey } from '../routing';
import { NavLink } from 'react-router';



type NavItem = { key: PageKey; label: string; icon: string };

const navItems: NavItem[] = [
  { key: 'dashboard', label: 'Dashboard', icon: '⌘' },
  { key: 'new-analysis', label: 'New Analysis', icon: '+' },
  { key: 'sequence-review', label: 'Sequence Review', icon: '✣' },
  { key: 'protein-annotation', label: 'Protein Annotation', icon: '◉' },
  { key: 'structure-pockets', label: 'Structure & Pockets', icon: '⬡' },
  { key: 'compound-ranking', label: 'Compound Ranking', icon: '◇' },
  { key: 'evidence-review', label: 'Evidence Review', icon: '▣' },
  { key: 'pipeline-status', label: 'Pipeline Status', icon: '⌁' },
  { key: 'system-status', label: 'System Status', icon: '▤' },
  { key: 'settings', label: 'Settings', icon: '⚙' }
];

export function Sidebar() {
  return (
    <aside className="sidebar">
      <div className="brand">
        <div className="brand-mark">✦</div>
        <div>
          <div className="brand-title">EXORCIST</div>
          <div className="brand-subtitle">Research Platform v2.4</div>
        </div>
      </div>

      <nav className="nav-list" aria-label="Primary navigation">
        {navItems.map((item) => (
          <NavLink
            key={item.key}
            to={pagePaths[item.key]}
            end={item.key === 'dashboard'}
            className={({ isActive }) =>
              `nav-item ${isActive ? 'active' : ''}`
            }
          >
            <span className="nav-icon">{item.icon}</span>
            <span>{item.label}</span>
          </NavLink>
        ))}
      </nav>

      <div className="sidebar-user">
        <div className="avatar">E</div>
        <div>
          <div className="user-name">EXORCIST Team</div>
          <div className="user-role">Research Lead</div>
        </div>
      </div>
    </aside>
  );
}
