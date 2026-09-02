import React, { useEffect, useState } from 'react';
import { Activity, ShieldCheck, Cpu, Github, ExternalLink } from 'lucide-react';
import { checkHealth } from '../services/api';

export const Navbar: React.FC = () => {
  const [status, setStatus] = useState<string>('checking...');
  const [isHealthy, setIsHealthy] = useState<boolean>(false);

  useEffect(() => {
    let isMounted = true;
    const ping = async () => {
      try {
        const res = await checkHealth();
        if (isMounted) {
          setStatus(res.status);
          setIsHealthy(true);
        }
      } catch (err) {
        if (isMounted) {
          setStatus('cold / waking up');
          setIsHealthy(false);
        }
      }
    };

    ping();
    const interval = setInterval(ping, 30000); // Check every 30s
    return () => {
      isMounted = false;
      clearInterval(interval);
    };
  }, []);

  return (
    <header className="border-b border-dark-800 bg-dark-900/80 backdrop-blur-md sticky top-0 z-50">
      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 h-16 flex items-center justify-between">
        <div className="flex items-center space-x-3">
          <div className="w-10 h-10 rounded-lg bg-genomics-green/10 border border-genomics-green/30 flex items-center justify-center text-genomics-green">
            <Cpu className="w-5 h-5 animate-pulse" />
          </div>
          <div>
            <div className="flex items-center space-x-2">
              <span className="font-bold text-lg tracking-wider text-white">GENOMESIGHT</span>
              <span className="text-xs px-2 py-0.5 rounded bg-genomics-green/20 text-genomics-green font-mono">v2.0 FastAPI</span>
            </div>
            <p className="text-xs text-gray-400">High-Performance Genomic Intelligence & Sequence Lab</p>
          </div>
        </div>

        <div className="flex items-center space-x-4">
          <div className="hidden sm:flex items-center space-x-2 bg-dark-800 px-3 py-1.5 rounded-full border border-dark-700 text-xs font-mono">
            <span className={`w-2 h-2 rounded-full ${isHealthy ? 'bg-genomics-green animate-ping' : 'bg-amber-500'}`} />
            <span className="text-gray-300">Backend:</span>
            <span className={isHealthy ? 'text-genomics-green' : 'text-amber-400'}>{status}</span>
          </div>

          <a
            href="https://github.com/Ardit-Mishra/genomesight"
            target="_blank"
            rel="noopener noreferrer"
            className="flex items-center space-x-1.5 text-sm text-gray-300 hover:text-white bg-dark-800 hover:bg-dark-700 px-3 py-1.5 rounded-lg border border-dark-700 transition"
          >
            <Github className="w-4 h-4" />
            <span>GitHub</span>
            <ExternalLink className="w-3 h-3 ml-0.5 text-gray-500" />
          </a>
        </div>
      </div>
    </header>
  );
};
