export function SafeNotice({ children }: { children: string }) {
  return (
    <div className="safe-notice">
      <span>⚠</span>
      <p>{children}</p>
    </div>
  );
}
