import React from 'react';
import { AnalyzeResponse } from '../services/api';
import { PieChart, Cpu } from 'lucide-react';

interface CompositionChartsProps {
  data: AnalyzeResponse;
}

export const CompositionCharts: React.FC<CompositionChartsProps> = ({ data }) => {
  const comp = data.statistics.composition;
  const total = Object.values(comp).reduce((a, b) => a + b, 0) || 1;
  const kmers = data.kmer_analysis.counts;

  return (
    <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
      {/* Nucleotide Composition */}
      <div className="bg-dark-900 border border-dark-800 rounded-xl p-6 shadow-lg">
        <div className="flex items-center space-x-2 mb-4">
          <PieChart className="w-5 h-5 text-genomics-cyan" />
          <h3 className="text-lg font-semibold text-white">Nucleotide Composition</h3>
        </div>

        <div className="space-y-3">
          {Object.entries(comp).map(([nt, count]) => {
            const pct = ((count / total) * 100).toFixed(1);
            const colors: Record<string, string> = {
              A: 'bg-genomics-green',
              T: 'bg-red-500',
              G: 'bg-genomics-amber',
              C: 'bg-genomics-blue',
              N: 'bg-gray-500'
            };

            return (
              <div key={nt} className="flex items-center justify-between">
                <div className="flex items-center space-x-3 w-24">
                  <span className={`w-3 h-3 rounded-full ${colors[nt] || 'bg-gray-400'}`} />
                  <span className="font-mono font-semibold text-gray-200">{nt}</span>
                </div>
                <div className="flex-1 mx-4">
                  <div className="w-full bg-dark-800 h-2.5 rounded-full overflow-hidden">
                    <div
                      className={`h-full rounded-full ${colors[nt] || 'bg-gray-400'}`}
                      style={{ width: `${pct}%` }}
                    />
                  </div>
                </div>
                <div className="w-32 text-right font-mono text-sm text-gray-400">
                  {count.toLocaleString()} <span className="text-xs text-gray-500">({pct}%)</span>
                </div>
              </div>
            );
          })}
        </div>
      </div>

      {/* K-Mer Frequency Analysis */}
      <div className="bg-dark-900 border border-dark-800 rounded-xl p-6 shadow-lg">
        <div className="flex items-center justify-between mb-4">
          <div className="flex items-center space-x-2">
            <Cpu className="w-5 h-5 text-genomics-purple" />
            <h3 className="text-lg font-semibold text-white">K-Mer Frequency (k={data.kmer_analysis.k})</h3>
          </div>
          <span className="text-xs px-2 py-0.5 rounded bg-dark-800 text-gray-400 font-mono">
            Engine: {data.kmer_analysis.engine}
          </span>
        </div>

        <div className="max-h-56 overflow-y-auto pr-2 space-y-2">
          {Object.entries(kmers).length === 0 ? (
            <p className="text-gray-500 text-sm">No k-mer data available.</p>
          ) : (
            <div className="grid grid-cols-2 sm:grid-cols-3 gap-2">
              {Object.entries(kmers).slice(0, 18).map(([kmer, cnt]) => (
                <div key={kmer} className="bg-dark-950 border border-dark-800 p-2 rounded-lg flex items-center justify-between font-mono text-xs">
                  <span className="text-genomics-cyan font-semibold">{kmer}</span>
                  <span className="text-gray-400">{cnt}</span>
                </div>
              ))}
            </div>
          )}
        </div>
      </div>
    </div>
  );
};
