import React from 'react';
import { Dna, BarChart3, Layers, Weight, ShieldAlert } from 'lucide-react';
import { AnalyzeResponse } from '../services/api';

interface SummaryCardsProps {
  data: AnalyzeResponse;
}

export const SummaryCards: React.FC<SummaryCardsProps> = ({ data }) => {
  const stats = data.statistics;
  const qual = data.quality_statistics;

  return (
    <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-4">
      {/* Total Length */}
      <div className="bg-dark-900 border border-dark-800 rounded-xl p-5 shadow-lg">
        <div className="flex items-center justify-between mb-2">
          <span className="text-xs uppercase tracking-wider text-gray-400 font-mono">Total Length</span>
          <Dna className="w-5 h-5 text-genomics-cyan" />
        </div>
        <div className="text-2xl font-bold font-mono text-white">
          {stats.total_length.toLocaleString()} <span className="text-xs font-normal text-gray-400">bp</span>
        </div>
        <p className="text-xs text-gray-500 mt-1">{stats.record_count} sequence record(s)</p>
      </div>

      {/* GC Content */}
      <div className="bg-dark-900 border border-dark-800 rounded-xl p-5 shadow-lg">
        <div className="flex items-center justify-between mb-2">
          <span className="text-xs uppercase tracking-wider text-gray-400 font-mono">GC Content</span>
          <BarChart3 className="w-5 h-5 text-genomics-green" />
        </div>
        <div className="text-2xl font-bold font-mono text-white">
          {stats.gc_content_percentage}%
        </div>
        <div className="w-full bg-dark-800 h-1.5 rounded-full mt-2 overflow-hidden">
          <div
            className="bg-genomics-green h-full rounded-full transition-all duration-500"
            style={{ width: `${Math.min(stats.gc_content_percentage, 100)}%` }}
          />
        </div>
      </div>

      {/* Molecular Weight */}
      <div className="bg-dark-900 border border-dark-800 rounded-xl p-5 shadow-lg">
        <div className="flex items-center justify-between mb-2">
          <span className="text-xs uppercase tracking-wider text-gray-400 font-mono">Approx. MW</span>
          <Weight className="w-5 h-5 text-genomics-purple" />
        </div>
        <div className="text-2xl font-bold font-mono text-white">
          {(stats.molecular_weight_approx / 1000).toFixed(1)} <span className="text-xs font-normal text-gray-400">kDa</span>
        </div>
        <p className="text-xs text-gray-500 mt-1">Single-strand estimate</p>
      </div>

      {/* Quality Score (if FASTQ) */}
      <div className="bg-dark-900 border border-dark-800 rounded-xl p-5 shadow-lg">
        <div className="flex items-center justify-between mb-2">
          <span className="text-xs uppercase tracking-wider text-gray-400 font-mono">Mean Phred Quality</span>
          <ShieldAlert className="w-5 h-5 text-genomics-amber" />
        </div>
        <div className="text-2xl font-bold font-mono text-white">
          {qual ? qual.mean_quality_score : 'N/A'} <span className="text-xs font-normal text-gray-400">{qual ? 'Q-score' : 'FASTA input'}</span>
        </div>
        <p className="text-xs text-gray-500 mt-1">{qual ? `Max Q: ${qual.max_quality}` : 'No quality array'}</p>
      </div>
    </div>
  );
};
