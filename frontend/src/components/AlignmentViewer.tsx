import React, { useState } from 'react';
import { alignSequences, AlignResponse } from '../services/api';
import { Split } from 'lucide-react';

export const AlignmentViewer: React.FC = () => {
  const [seq1, setSeq1] = useState('');
  const [seq2, setSeq2] = useState('');
  const [result, setResult] = useState<AlignResponse | null>(null);
  const [isLoading, setIsLoading] = useState(false);

  const handleAlign = async () => {
    setIsLoading(true);
    try {
      const res = await alignSequences({ seq1, seq2, mode: 'global' });
      setResult(res);
    } catch (err) {
      console.error(err);
    } finally {
      setIsLoading(false);
    }
  };

  return (
    <div className="bg-dark-900 border border-dark-800 rounded-xl p-6 shadow-lg">
      <div className="flex items-center space-x-2 mb-4">
        <Split className="w-5 h-5 text-genomics-blue" />
        <h3 className="text-lg font-semibold text-white">Pairwise Sequence Alignment</h3>
      </div>

      <div className="grid grid-cols-2 gap-4 mb-4">
        <textarea value={seq1} onChange={e => setSeq1(e.target.value)} placeholder="Sequence 1" className="bg-dark-950 border border-dark-700 p-2 text-sm font-mono rounded" rows={3} />
        <textarea value={seq2} onChange={e => setSeq2(e.target.value)} placeholder="Sequence 2" className="bg-dark-950 border border-dark-700 p-2 text-sm font-mono rounded" rows={3} />
      </div>

      <button onClick={handleAlign} disabled={isLoading} className="w-full bg-genomics-blue text-white py-2 rounded-lg text-sm font-semibold hover:bg-blue-600 transition">
        {isLoading ? 'Aligning...' : 'Run Alignment'}
      </button>

      {result && (
        <div className="mt-4 p-4 bg-dark-950 border border-dark-800 rounded-lg">
          <p className="text-xs font-mono text-gray-400 mb-2">Score: {result.score.toFixed(2)} | Identity: {result.identity_percentage.toFixed(1)}%</p>
          <div className="font-mono text-xs overflow-x-auto whitespace-pre p-2 bg-black rounded text-genomics-green">
            {result.aligned_seq1}
            <br />
            {result.match_line}
            <br />
            {result.aligned_seq2}
          </div>
        </div>
      )}
    </div>
  );
};
