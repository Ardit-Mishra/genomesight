import React, { useState } from 'react';
import { calculateCodons, translateSequence, CodonResponse, TranslateResponse } from '../services/api';
import { BookOpen, Dna } from 'lucide-react';

export const CodonTable: React.FC = () => {
  const [seq, setSeq] = useState('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG');
  const [codonRes, setCodonRes] = useState<CodonResponse | null>(null);
  const [transRes, setTransRes] = useState<TranslateResponse | null>(null);

  const handleAnalyze = async () => {
    try {
      const c = await calculateCodons({ sequence: seq });
      const t = await translateSequence({ sequence: seq, frame: 1 });
      setCodonRes(c);
      setTransRes(t);
    } catch (err) {
      console.error(err);
    }
  };

  return (
    <div className="bg-dark-900 border border-dark-800 rounded-xl p-6 shadow-lg">
      <div className="flex items-center space-x-2 mb-4">
        <BookOpen className="w-5 h-5 text-genomics-amber" />
        <h3 className="text-lg font-semibold text-white">Translation & Codon Usage</h3>
      </div>

      <div className="mb-4">
        <input type="text" value={seq} onChange={e => setSeq(e.target.value)} className="w-full bg-dark-950 border border-dark-700 p-2 text-sm font-mono rounded" placeholder="Coding sequence..." />
      </div>

      <button onClick={handleAnalyze} className="w-full bg-genomics-amber text-black font-semibold py-2 rounded-lg text-sm hover:opacity-90 transition mb-4">
        Translate & Calculate RSCU
      </button>

      {transRes && (
        <div className="mb-4 p-3 bg-dark-950 border border-dark-800 rounded">
          <span className="text-xs font-mono text-gray-400">Protein Translation (Frame 1):</span>
          <p className="font-mono text-xs text-genomics-green break-all mt-1">{transRes.protein_sequence}</p>
        </div>
      )}

      {codonRes && (
        <div className="max-h-40 overflow-y-auto">
          <table className="w-full text-left font-mono text-xs">
            <thead>
              <tr className="border-b border-dark-800 text-gray-400">
                <th className="pb-2">Codon</th>
                <th className="pb-2">Count</th>
                <th className="pb-2">RSCU</th>
              </tr>
            </thead>
            <tbody>
              {Object.entries(codonRes.codon_counts).map(([codon, count]) => (
                <tr key={codon} className="border-b border-dark-800/50">
                  <td className="py-1 text-genomics-cyan">{codon}</td>
                  <td className="py-1 text-gray-300">{count}</td>
                  <td className="py-1 text-gray-400">{(codonRes.rscu[codon] || 0).toFixed(2)}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}
    </div>
  );
};
