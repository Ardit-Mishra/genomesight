import axios, { AxiosInstance, AxiosError } from 'axios';

const API_BASE_URL = (import.meta as any).env.VITE_API_URL || 'http://localhost:8000';

export interface AnalyzeRequest {
  sequence?: string;
  file_content?: string;
}

export interface AnalyzeResponse {
  success: boolean;
  statistics: {
    record_count: number;
    total_length: number;
    primary_id: string;
    primary_description: string;
    gc_content_percentage: number;
    composition: Record<string, number>;
    molecular_weight_approx: number;
  };
  quality_statistics?: {
    mean_quality_score: number;
    max_quality: number;
    min_quality: number;
  };
  kmer_analysis: {
    k: number;
    engine: string;
    counts: Record<string, number>;
  };
}

export interface AlignRequest {
  seq1: string;
  seq2: string;
  mode: string;
}

export interface AlignResponse {
  success: boolean;
  score: number;
  aligned_seq1: string;
  aligned_seq2: string;
  match_line: string;
  identity_percentage: number;
}

export interface TranslateRequest {
  sequence: string;
  frame: number;
}

export interface TranslateResponse {
  success: boolean;
  protein_sequence: string;
  frame: number;
}

export interface CodonRequest {
  sequence: string;
}

export interface CodonResponse {
  success: boolean;
  codon_counts: Record<string, number>;
  rscu: Record<string, number>;
}

// Create axios instance
const apiClient: AxiosInstance = axios.create({
  baseURL: API_BASE_URL,
  timeout: 30000, // 30s timeout for cold starts
  headers: {
    'Content-Type': 'application/json',
  },
});

// Retry interceptor for Render cold start (handling 5xx or network errors with exponential backoff)
apiClient.interceptors.response.use(
  (response) => response,
  async (error: AxiosError) => {
    const config = error.config as any;
    if (!config) {
      return Promise.reject(error);
    }

    // Initialize retry count if not present
    config.__retryCount = config.__retryCount || 0;
    const maxRetries = 3;
    const delay = Math.pow(2, config.__retryCount) * 1000; // 2s, 4s, 8s

    if (config.__retryCount < maxRetries && (error.code === 'ECONNABORTED' || !error.response || (error.response.status >= 500 && error.response.status <= 599))) {
      config.__retryCount += 1;
      console.warn(`[GenomeSight API] Backend waking up / cold start detected. Retrying request (${config.__retryCount}/${maxRetries}) in ${delay}ms...`);

      await new Promise((resolve) => setTimeout(resolve, delay));
      return apiClient(config);
    }

    return Promise.reject(error);
  }
);

export const analyzeSequence = async (data: AnalyzeRequest): Promise<AnalyzeResponse> => {
  const response = await apiClient.post<AnalyzeResponse>('/api/analyze', data);
  return response.data;
};

export const alignSequences = async (data: AlignRequest): Promise<AlignResponse> => {
  const response = await apiClient.post<AlignResponse>('/api/align', data);
  return response.data;
};

export const translateSequence = async (data: TranslateRequest): Promise<TranslateResponse> => {
  const response = await apiClient.post<TranslateResponse>('/api/translate', data);
  return response.data;
};

export const calculateCodons = async (data: CodonRequest): Promise<CodonResponse> => {
  const response = await apiClient.post<CodonResponse>('/api/codons', data);
  return response.data;
};

export const checkHealth = async (): Promise<{ status: string; service: string; version: string }> => {
  const response = await apiClient.get('/api/health');
  return response.data;
};

// ---------------------------------------------------------------------------
// ORF detection, motif search and restriction sites
// Ported from the Streamlit workbench so the deployed app regains the
// capability its documentation already claimed.
// ---------------------------------------------------------------------------
export interface OrfRequest {
  sequence: string;
  min_length?: number;
  include_reverse?: boolean;
  use_alternative_starts?: boolean;
}

export interface OrfRecord {
  start: number;
  end: number;
  length_nt: number;
  length_aa: number;
  frame: number;
  strand: string;
  start_codon: string;
  stop_codon: string;
  protein: string;
  gc_content: number;
}

export interface OrfResponse {
  success: boolean;
  coordinate_note: string;
  total: number;
  summary: {
    total: number;
    average_length?: number;
    max_length?: number;
    min_length?: number;
    by_frame?: Record<string, number>;
  };
  orfs: OrfRecord[];
}

export interface MotifMatchRecord {
  pattern: string;
  sequence_id: string;
  start_1based: number;
  end: number;
  matched_sequence: string;
  context: string;
}

export interface MotifResponse {
  success: boolean;
  pattern: string;
  regex: string;
  total: number;
  matches: MotifMatchRecord[];
}

export const detectOrfs = async (data: OrfRequest): Promise<OrfResponse> => {
  const response = await apiClient.post<OrfResponse>('/api/orfs', data);
  return response.data;
};

export const findMotifs = async (data: { sequence: string; pattern: string }): Promise<MotifResponse> => {
  const response = await apiClient.post<MotifResponse>('/api/motifs', data);
  return response.data;
};

export const findRestrictionSites = async (data: { sequence: string; enzyme: string }): Promise<MotifResponse> => {
  const response = await apiClient.post<MotifResponse>('/api/restriction-sites', data);
  return response.data;
};

export const listEnzymes = async (): Promise<Record<string, string>> => {
  const response = await apiClient.get<{ success: boolean; enzymes: Record<string, string> }>('/api/enzymes');
  return response.data.enzymes;
};
