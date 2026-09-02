import React, { useState } from 'react';
import { Upload, FileText, Sparkles, AlertCircle } from 'lucide-react';

interface FileUploadProps {
  onAnalyze: (sequenceOrContent: string) => void;
  isLoading: boolean;
}

const SAMPLE_FASTA = `>SARS-CoV-2_spike_glycoprotein_partial
ATGTIGTTGCTTTCGCTTTTCTTTTTCTTTGCTTTTCTTGCTTTCGCTTTT
CTTTTTCTTTGCTTTTCTTGCTTTCGCTTTTCTTTTTCTTTGCTTTTCTTG
CTTTCGCTTTTCTTTTTCTTTGCTTTTCTTGCTTTCGCTTTTCTTTTTCTT
TGCTTTTCTTGCTTTCGCTTTTCTTTTTCTTTGCTTTTCTTGCTTTCGCTT
TTCTTTTTCTTTGCTTTTCTTGCTTTCGCTTTTCTTTTTCTTTGCTTTTCT
TGCTTTCGCTTTTCTTTTTCTTTGCTTTTCTTGCTTTCGCTTTTCTTTTTC
TTTGCTTTTCTTGCTTTCGCTTTTCTTTTTCTTTGCTTTTCTTGCTTTCGC
TTTTCTTTTTCTTTGCTTTTCTTGCTTTCGCTTTTCTTTTTCTTTGCTTTT`;

const SAMPLE_LAMBDA = `>Lambda_phage_control_region
GGGCGGCGACCTGBCCGGCCCCGGCCCCGGCGACCTGBCCGGCCCCGGCCCCGG
CGACCTGBCCGGCCCCGGCCCCGGGCGGCGACCTGBCCGGCCCCGGCCCCGGGCG
GCGACCTGBCCGGCCCCGGCCCCGGGCGGCGACCTGBCCGGCCCCGGCCCCGGGC
GGCGACCTGBCCGGCCCCGGCCCCGGGCGGCGACCTGBCCGGCCCCGGCCCCGGGCG`;

export const FileUpload: React.FC<FileUploadProps> = ({ onAnalyze, isLoading }) => {
  const [inputVal, setInputVal] = useState<string>('');
  const [error, setError] = useState<string | null>(null);

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    if (!inputVal.trim()) {
      setError('Please enter or upload a valid FASTA sequence.');
      return;
    }
    setError(null);
    onAnalyze(inputVal);
  };

  const handleFileUpload = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (event) => {
      const content = event.target?.result as string;
      if (content) {
        setInputVal(content);
        setError(null);
        onAnalyze(content);
      }
    };
    reader.onerror = () => {
      setError('Failed to read uploaded file.');
    };
    reader.readAsText(file);
  };

  return (
    <div className="bg-dark-900 border border-dark-800 rounded-xl p-6 shadow-xl">
      <div className="flex items-center justify-between mb-4">
        <div className="flex items-center space-x-2">
          <FileText className="w-5 h-5 text-genomics-cyan" />
          <h2 className="text-lg font-semibold text-white">Sequence Input & Ingestion</h2>
        </div>
        <div className="flex items-center space-x-2">
          <button
            type="button"
            onClick={() => { setInputVal(SAMPLE_FASTA); onAnalyze(SAMPLE_FASTA); }}
            className="text-xs bg-dark-800 hover:bg-dark-700 text-genomics-cyan px-2.5 py-1.5 rounded border border-dark-700 transition flex items-center space-x-1"
          >
            <Sparkles className="w-3 h-3" />
            <span>Load SARS-CoV-2 Sample</span>
          </button>
          <button
            type="button"
            onClick={() => { setInputVal(SAMPLE_LAMBDA); onAnalyze(SAMPLE_LAMBDA); }}
            className="text-xs bg-dark-800 hover:bg-dark-700 text-genomics-green px-2.5 py-1.5 rounded border border-dark-700 transition flex items-center space-x-1"
          >
            <Sparkles className="w-3 h-3" />
            <span>Load Lambda Sample</span>
          </button>
        </div>
      </div>

      <form onSubmit={handleSubmit}>
        <div className="mb-4">
          <textarea
            rows={6}
            value={inputVal}
            onChange={(e) => setInputVal(e.target.value)}
            placeholder="Paste FASTA sequence here (e.g. >seq1\nATCGATCG...)"
            className="w-full bg-dark-950 border border-dark-700 rounded-lg p-3 font-mono text-sm text-gray-200 focus:outline-none focus:border-genomics-cyan transition resize-y"
          />
        </div>

        {error && (
          <div className="mb-4 p-3 bg-red-950/40 border border-red-800/60 rounded-lg flex items-center space-x-2 text-red-300 text-sm">
            <AlertCircle className="w-4 h-4 flex-shrink-0" />
            <span>{error}</span>
          </div>
        )}

        <div className="flex flex-col sm:flex-row items-center justify-between gap-4">
          <label className="w-full sm:w-auto flex items-center justify-center space-x-2 bg-dark-800 hover:bg-dark-700 text-gray-300 px-4 py-2.5 rounded-lg border border-dark-700 cursor-pointer transition text-sm">
            <Upload className="w-4 h-4 text-genomics-cyan" />
            <span>Upload FASTA/FASTQ File</span>
            <input type="file" accept=".fasta,.fa,.fastq,.txt" onChange={handleFileUpload} className="hidden" />
          </label>

          <button
            type="submit"
            disabled={isLoading}
            className="w-full sm:w-auto bg-gradient-to-r from-genomics-cyan to-genomics-blue hover:opacity-90 text-black font-semibold px-6 py-2.5 rounded-lg transition shadow-lg shadow-genomics-cyan/20 disabled:opacity-50 text-sm flex items-center justify-center space-x-2"
          >
            {isLoading ? (
              <>
                <div className="w-4 h-4 border-2 border-black border-t-transparent rounded-full animate-spin" />
                <span>Analyzing Sequence (Waking Backend...)</span>
              </>
            ) : (
              <span>Run Genomic Analysis</span>
            )}
          </button>
        </div>
      </form>
    </div>
  );
};
