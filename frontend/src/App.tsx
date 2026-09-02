import React, { useState } from 'react';
import { Navbar } from './components/Navbar';
import { FileUpload } from './components/FileUpload';
import { SummaryCards } from './components/SummaryCards';
import { CompositionCharts } from './components/CompositionCharts';
import { AlignmentViewer } from './components/AlignmentViewer';
import { CodonTable } from './components/CodonTable';
import { analyzeSequence, AnalyzeResponse } from './services/api';

export function App() {
  const [data, setData] = useState<AnalyzeResponse | null>(null);
  const [isLoading, setIsLoading] = useState<boolean>(false);

  const handleAnalyze = async (content: string) => {
    setIsLoading(true);
    try {
      const res = await analyzeSequence({ sequence: content });
      setData(res);
    } catch (err) {
      console.error('Failed to analyze sequence:', err);
      alert('Analysis failed. Make sure the backend is reachable or check input formatting.');
    } finally {
      setIsLoading(false);
    }
  };

  return (
    <div className="min-h-screen bg-dark-950 text-gray-100 flex flex-col">
      <Navbar />

      <main className="flex-1 max-w-7xl w-full mx-auto px-4 sm:px-6 lg:px-8 py-8 space-y-8">
        {/* Hero section */}
        <div className="text-center sm:text-left flex flex-col sm:flex-row sm:items-center sm:justify-between gap-4 border-b border-dark-800 pb-6">
          <div>
            <h1 className="text-3xl font-extrabold text-white tracking-tight">Genomic Sequence Analysis Studio</h1>
            <p className="text-sm text-gray-400 mt-1">Industrial Utilitarian Bioinformatics — Powered by FastAPI & BioPython</p>
          </div>
          <div className="flex items-center space-x-2">
            <span className="px-3 py-1 bg-dark-900 border border-dark-800 rounded-lg text-xs font-mono text-genomics-green">
              Status: Production Ready
            </span>
          </div>
        </div>

        {/* Input Section */}
        <FileUpload onAnalyze={handleAnalyze} isLoading={isLoading} />

        {/* Results Section */}
        {data ? (
          <div className="space-y-8 animate-fadeIn">
            <SummaryCards data={data} />
            <CompositionCharts data={data} />

            <div className="grid grid-cols-1 lg:grid-cols-2 gap-8">
              <AlignmentViewer />
              <CodonTable />
            </div>
          </div>
        ) : (
          <div className="bg-dark-900/40 border border-dashed border-dark-800 rounded-2xl p-12 text-center text-gray-500">
            <p className="text-base font-mono">No sequence ingested yet. Load a sample or paste FASTA above to begin analysis.</p>
          </div>
        )}
      </main>

      <footer className="border-t border-dark-800 py-6 text-center text-xs text-gray-500 font-mono">
        GenomeSight v2.0 · $0 Cost Architecture (Vercel + Render) · Ardit Mishra Portfolio
      </footer>
    </div>
  );
}

export default App;
