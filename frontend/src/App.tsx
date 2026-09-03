import { ChangeEvent, FormEvent, useMemo, useState } from "react"
import {
  Dna,
  FileUp,
  FlaskConical,
  GitCompareArrows,
  LoaderCircle,
  Play,
  ShieldCheck,
} from "lucide-react"

import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert"
import { Badge } from "@/components/ui/badge"
import { Button } from "@/components/ui/button"
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card"
import { Input } from "@/components/ui/input"
import { Label } from "@/components/ui/label"
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@/components/ui/table"
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs"
import { Textarea } from "@/components/ui/textarea"
import {
  AlignResponse,
  AnalyzeResponse,
  CodonResponse,
  alignSequences,
  analyzeSequence,
  calculateCodons,
  translateSequence,
} from "@/services/api"

const SAMPLE = ">example_sequence\nATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
type Operation = "analysis" | "alignment" | "codons" | null
type InputSummary = { format: string; characterCount: number; baseCount: number | null; state: string; detail: string; preview: string; invalid: boolean }

function App() {
  const [sequence, setSequence] = useState("")
  const [analysis, setAnalysis] = useState<AnalyzeResponse | null>(null)
  const [alignment, setAlignment] = useState<AlignResponse | null>(null)
  const [codons, setCodons] = useState<CodonResponse | null>(null)
  const [protein, setProtein] = useState("")
  const [activeOperation, setActiveOperation] = useState<Operation>(null)
  const [error, setError] = useState<string | null>(null)
  const [seq1, setSeq1] = useState("ACGT")
  const [seq2, setSeq2] = useState("AGT")
  const [codingSequence, setCodingSequence] = useState("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")

  const composition = useMemo(() => analysis?.statistics.composition ?? {}, [analysis])
  const inputSummary = useMemo(() => summarizeInput(sequence), [sequence])
  const isBusy = activeOperation !== null

  async function runAnalysis(content = sequence) {
    setActiveOperation("analysis")
    setError(null)
    try {
      setAnalysis(await analyzeSequence({ sequence: content }))
    } catch (cause) {
      setError(readError(cause, "Analysis could not be completed. Check the sequence and try again."))
    } finally {
      setActiveOperation(null)
    }
  }

  async function loadFile(event: ChangeEvent<HTMLInputElement>) {
    const input = event.target
    const file = input.files?.[0]
    if (!file) return
    setActiveOperation("analysis")
    setError(null)
    try {
      const content = await file.text()
      setSequence(content)
      setAnalysis(await analyzeSequence({ file_content: content }))
    } catch (cause) {
      setError(readError(cause, "The selected file could not be analyzed."))
    } finally {
      input.value = ""
      setActiveOperation(null)
    }
  }

  async function runAlignment(event: FormEvent) {
    event.preventDefault()
    setActiveOperation("alignment")
    setError(null)
    try {
      setAlignment(await alignSequences({ seq1, seq2, mode: "global" }))
    } catch (cause) {
      setError(readError(cause, "Alignment could not be completed."))
    } finally {
      setActiveOperation(null)
    }
  }

  async function runCodons(event: FormEvent) {
    event.preventDefault()
    setActiveOperation("codons")
    setError(null)
    try {
      const [codonResponse, translationResponse] = await Promise.all([
        calculateCodons({ sequence: codingSequence }),
        translateSequence({ sequence: codingSequence, frame: 1 }),
      ])
      setCodons(codonResponse)
      setProtein(translationResponse.protein_sequence)
    } catch (cause) {
      setError(readError(cause, "Translation and codon analysis could not be completed."))
    } finally {
      setActiveOperation(null)
    }
  }

  return (
    <div className="min-h-screen">
      <header className="border-b bg-background/90 backdrop-blur">
        <div className="mx-auto flex h-16 max-w-6xl items-center justify-between px-4 sm:px-6">
          <div className="flex items-center gap-3">
            <div className="brand-mark relative grid size-9 place-items-center overflow-hidden rounded-lg border border-primary/40 bg-primary/10 text-primary">
              <Dna className="size-5" aria-hidden="true" />
            </div>
            <div>
              <p className="font-mono text-sm font-semibold tracking-[.16em]">GENOMESIGHT</p>
              <p className="text-xs text-muted-foreground">Sequence analysis workbench</p>
            </div>
          </div>
          <Badge variant="outline" className="hidden gap-1.5 font-mono text-xs sm:flex">
            <span className="status-beacon" aria-hidden="true" />
            <ShieldCheck className="size-3" aria-hidden="true" />
            Validation-first
          </Badge>
        </div>
      </header>

      <main className="mx-auto max-w-7xl px-4 py-6 sm:px-6 sm:py-8">
        <section className="workspace-heading" aria-labelledby="workspace-title">
          <div className="min-w-0">
            <div className="flex items-center gap-3"><p className="font-mono text-xs uppercase tracking-[.14em] text-primary">GenomeSight / Workspace</p><span className="font-mono text-[0.68rem] text-muted-foreground">v2.0</span></div>
            <h1 id="workspace-title" className="mt-3 font-heading text-3xl font-semibold tracking-[-0.02em] sm:text-4xl">Sequence workbench</h1>
            <p className="mt-2 max-w-2xl text-sm leading-6 text-muted-foreground">Validate sequence input, inspect composition, compare aligned reads, and calculate coding bias in one local workflow.</p>
          </div>
          <dl className="workspace-status" aria-label="Workspace status"><div><dt>Pipeline</dt><dd><span className="status-beacon" aria-hidden="true" />Validation active</dd></div><div><dt>Mode</dt><dd>Local operation</dd></div><div><dt>Output</dt><dd>Reproducible readout</dd></div></dl>
        </section>

        {error && <Alert variant="destructive" className="mb-6" aria-live="assertive"><AlertTitle>Analysis needs attention</AlertTitle><AlertDescription>{error}</AlertDescription></Alert>}

        <Tabs defaultValue="analyze" className="min-w-0 space-y-6">
          <TabsList variant="line" className="tool-tabs grid h-auto w-full grid-cols-3 gap-0 overflow-visible rounded-lg border bg-card p-1">
            <TabsTrigger value="analyze" className="min-h-14 min-w-0 flex-col gap-1 whitespace-normal px-2 py-2 text-center text-xs leading-tight data-[state=active]:bg-muted data-[state=active]:text-primary"><FileUp className="size-4" aria-hidden="true" /><span><span className="sm:hidden">Analyze</span><span className="hidden sm:inline">Sequence analysis</span></span></TabsTrigger>
            <TabsTrigger value="align" className="min-h-14 min-w-0 flex-col gap-1 whitespace-normal px-2 py-2 text-center text-xs leading-tight data-[state=active]:bg-muted data-[state=active]:text-primary"><GitCompareArrows className="size-4" aria-hidden="true" /><span><span className="sm:hidden">Align</span><span className="hidden sm:inline">Pairwise alignment</span></span></TabsTrigger>
            <TabsTrigger value="codons" className="min-h-14 min-w-0 flex-col gap-1 whitespace-normal px-2 py-2 text-center text-xs leading-tight data-[state=active]:bg-muted data-[state=active]:text-primary"><FlaskConical className="size-4" aria-hidden="true" /><span><span className="sm:hidden">Codons</span><span className="hidden sm:inline">Translation &amp; RSCU</span></span></TabsTrigger>
          </TabsList>

          <TabsContent value="analyze" className="space-y-6">
            <Card className="workbench-panel rounded-lg shadow-none">
              <CardHeader><CardTitle className="flex items-center gap-2"><FileUp className="size-4 text-primary" aria-hidden="true" />Input</CardTitle><CardDescription>Raw DNA/RNA, FASTA, or FASTQ. Whitespace is normalized; unsupported symbols are rejected.</CardDescription></CardHeader>
              <CardContent className="space-y-4">
                <div className="grid gap-6 lg:grid-cols-[minmax(0,1fr)_17rem]"><div className="space-y-2"><Label htmlFor="sequence">Sequence or FASTA content</Label><Textarea id="sequence" value={sequence} onChange={(event) => setSequence(event.target.value)} placeholder={"ATGCGATC… or >sequence\nATGCGATC…"} className="min-h-56 resize-y font-mono text-sm" /></div><InputInspector summary={inputSummary} /></div>
                <div className="flex flex-col gap-3 sm:flex-row sm:flex-wrap sm:items-center">
                  <Button type="button" variant="outline" className="min-h-11" disabled={isBusy} onClick={() => { setSequence(SAMPLE); void runAnalysis(SAMPLE) }}>Load validated example</Button>
                  <Label className="inline-flex min-h-11 cursor-pointer items-center gap-2 rounded-lg border px-3 text-sm font-medium transition-colors hover:bg-muted focus-within:ring-2 focus-within:ring-ring focus-within:ring-offset-2 focus-within:ring-offset-background"><FileUp className="size-4" aria-hidden="true" />Upload FASTA/FASTQ<Input className="sr-only" type="file" accept=".fasta,.fa,.fastq,.fq,.txt" disabled={isBusy} onChange={loadFile} /></Label>
                  <Button className="min-h-11 sm:ml-auto" onClick={() => void runAnalysis()} disabled={isBusy || !sequence.trim()}>{activeOperation === "analysis" ? <LoaderCircle className="size-4 animate-spin" aria-hidden="true" /> : <Play className="size-4" aria-hidden="true" />}Run analysis</Button>
                </div>
                {activeOperation === "analysis" && <div className="sequence-scan" role="status" aria-live="polite"><span className="font-mono text-xs text-primary">VALIDATING SEQUENCE</span><span className="text-xs text-muted-foreground">Normalizing IUPAC symbols and calculating composition</span></div>}
              </CardContent>
            </Card>
            {analysis && <AnalysisResults data={analysis} composition={composition} />}
          </TabsContent>

          <TabsContent value="align"><Card className="workbench-panel rounded-lg shadow-none"><CardHeader><CardTitle className="flex items-center gap-2"><GitCompareArrows className="size-4 text-primary" aria-hidden="true" />Pairwise alignment</CardTitle><CardDescription>Global Needleman–Wunsch alignment. Identity is calculated over aligned columns, including gaps.</CardDescription></CardHeader><CardContent><form className="space-y-4" onSubmit={runAlignment}><div className="grid gap-4 md:grid-cols-2"><Field id="seq1" label="Sequence 1" value={seq1} onChange={setSeq1} /><Field id="seq2" label="Sequence 2" value={seq2} onChange={setSeq2} /></div><Button className="min-h-11" disabled={isBusy}>{activeOperation === "alignment" && <LoaderCircle className="size-4 animate-spin" aria-hidden="true" />}Run alignment</Button></form>{alignment && <AlignmentResult data={alignment} />}</CardContent></Card></TabsContent>

          <TabsContent value="codons"><Card className="workbench-panel rounded-lg shadow-none"><CardHeader><CardTitle className="flex items-center gap-2"><FlaskConical className="size-4 text-primary" aria-hidden="true" />Translation &amp; RSCU</CardTitle><CardDescription>Frame 1 translation and relative synonymous codon usage for observed sense codons.</CardDescription></CardHeader><CardContent><form className="space-y-4" onSubmit={runCodons}><div className="space-y-2"><Label htmlFor="coding-sequence">Coding sequence</Label><Textarea id="coding-sequence" className="min-h-32 font-mono" value={codingSequence} onChange={(event) => setCodingSequence(event.target.value)} /></div><Button className="min-h-11" disabled={isBusy}>{activeOperation === "codons" && <LoaderCircle className="size-4 animate-spin" aria-hidden="true" />}Translate &amp; calculate</Button></form>{protein && <ProteinResult protein={protein} />}{codons && <CodonResults data={codons} />}</CardContent></Card></TabsContent>
        </Tabs>
      </main>
    </div>
  )
}

function InputInspector({ summary }: { summary: InputSummary }) {
  return <aside className="input-inspector" aria-label="Live input check"><div className="flex items-center justify-between gap-3"><p className="font-mono text-xs uppercase tracking-[.12em] text-muted-foreground">Live input check</p><Badge variant="outline" className={summary.invalid ? "border-destructive/60 text-destructive" : "border-primary/40 text-primary"}>{summary.state}</Badge></div><dl className="mt-5 space-y-4"><div className="inspector-row"><dt>Detected format</dt><dd>{summary.format}</dd></div><div className="inspector-row"><dt>Previewed bases</dt><dd>{summary.baseCount === null ? "Server check" : summary.baseCount.toLocaleString()}</dd></div><div className="inspector-row"><dt>Characters</dt><dd>{summary.characterCount.toLocaleString()}</dd></div></dl><div className="mt-5 border-t pt-4"><p className="text-xs leading-5 text-muted-foreground">{summary.detail}</p>{summary.preview ? <div className="sequence-ribbon mt-3" aria-label="Sequence preview"><span>{summary.preview} · {summary.preview} · {summary.preview}</span></div> : <p className="mt-3 font-mono text-xs text-muted-foreground">Awaiting sequence input</p>}</div></aside>
}
function Field({ id, label, value, onChange }: { id: string; label: string; value: string; onChange: (value: string) => void }) {
  return <div className="space-y-2"><Label htmlFor={id}>{label}</Label><Textarea id={id} className="min-h-28 font-mono" value={value} onChange={(event) => onChange(event.target.value)} /></div>
}

function AnalysisResults({ data, composition }: { data: AnalyzeResponse; composition: Record<string, number> }) {
  const stats = data.statistics
  const metrics = [["Records", stats.record_count], ["Total length", `${Number(stats.total_length).toLocaleString()} bp`], ["GC content", `${stats.gc_content_percentage}%`], ["K-mer engine", data.kmer_analysis.engine]]
  return <section className="motion-enter space-y-6" aria-label="Analysis results" aria-live="polite"><div className="flex items-baseline justify-between gap-4"><h2 className="text-lg font-medium">Analysis results</h2><span className="font-mono text-xs text-primary">Complete</span></div><dl className="grid gap-px overflow-hidden rounded-lg border bg-border sm:grid-cols-2 lg:grid-cols-4">{metrics.map(([label, value]) => <div key={String(label)} className="bg-card p-4"><dt className="text-xs text-muted-foreground">{label}</dt><dd className="mt-2 font-mono text-lg font-medium">{value}</dd></div>)}</dl><div className="grid gap-6 lg:grid-cols-5"><Card className="rounded-lg shadow-none lg:col-span-2"><CardHeader><CardTitle className="text-base">Nucleotide composition</CardTitle><CardDescription>Counts are aggregated across all records.</CardDescription></CardHeader><CardContent className="space-y-3">{Object.entries(composition).map(([base, count], index) => { const percentage = stats.total_length ? (count / stats.total_length) * 100 : 0; return <div key={base}><div className="mb-1 flex justify-between gap-3 font-mono text-xs"><span>{base}</span><span className="text-right">{count.toLocaleString()} · {percentage.toFixed(1)}%</span></div><div className="h-1.5 overflow-hidden rounded-full bg-muted" role="progressbar" aria-label={`${base} composition`} aria-valuemin={0} aria-valuemax={100} aria-valuenow={Number(percentage.toFixed(1))}><div className="metric-bar h-full bg-primary" style={{ width: `${percentage}%`, animationDelay: `${index * 70}ms` }} /></div></div> })}</CardContent></Card><Card className="rounded-lg shadow-none lg:col-span-3"><CardHeader><CardTitle className="text-base">Most frequent 3-mers</CardTitle><CardDescription>Ranked by count, then alphabetically for ties.</CardDescription></CardHeader><CardContent className="overflow-x-auto"><Table><TableHeader><TableRow><TableHead>3-mer</TableHead><TableHead className="text-right">Count</TableHead></TableRow></TableHeader><TableBody>{Object.entries(data.kmer_analysis.counts).slice(0, 10).map(([kmer, count]) => <TableRow key={kmer}><TableCell className="font-mono text-primary">{kmer}</TableCell><TableCell className="text-right font-mono">{count}</TableCell></TableRow>)}</TableBody></Table></CardContent></Card></div></section>
}

function AlignmentResult({ data }: { data: AlignResponse }) {
  return <section className="motion-enter mt-6" aria-label="Alignment result" aria-live="polite"><div className="mb-3 flex items-baseline justify-between gap-4"><h2 className="text-sm font-medium">Alignment result</h2><span className="font-mono text-xs text-primary">Complete</span></div><pre className="overflow-x-auto rounded-lg border bg-muted/30 p-4 font-mono text-sm leading-6">{data.aligned_seq1}{"\n"}{data.match_line}{"\n"}{data.aligned_seq2}{"\n\n"}Score: {data.score.toFixed(2)} · Identity: {data.identity_percentage.toFixed(1)}%</pre></section>
}

function ProteinResult({ protein }: { protein: string }) {
  return <section className="mt-6 rounded-lg border bg-muted/30 p-4" aria-label="Protein translation" aria-live="polite"><p className="text-xs text-muted-foreground">Protein translation · frame 1</p><p className="mt-1 break-all font-mono text-sm">{protein}</p></section>
}

function CodonResults({ data }: { data: CodonResponse }) {
  return <section className="motion-enter mt-6" aria-label="Observed codons" aria-live="polite"><div className="mb-3 flex items-baseline justify-between gap-4"><h2 className="text-sm font-medium">Observed codons</h2><span className="font-mono text-xs text-primary">Complete</span></div><div className="overflow-x-auto rounded-lg border"><Table><TableHeader><TableRow><TableHead>Codon</TableHead><TableHead className="text-right">Count</TableHead><TableHead className="text-right">RSCU</TableHead></TableRow></TableHeader><TableBody>{Object.entries(data.codon_counts).map(([codon, count]) => <TableRow key={codon}><TableCell className="font-mono text-primary">{codon}</TableCell><TableCell className="text-right font-mono">{count}</TableCell><TableCell className="text-right font-mono">{data.rscu[codon]?.toFixed(3) ?? "—"}</TableCell></TableRow>)}</TableBody></Table></div></section>
}

function summarizeInput(input: string): InputSummary {
  const raw = input.trim()
  if (!raw) return { format: "—", characterCount: 0, baseCount: 0, state: "Awaiting input", detail: "Paste a raw sequence or load a FASTA/FASTQ file to inspect it before analysis.", preview: "", invalid: false }
  const format = raw.startsWith(">") ? "FASTA" : raw.startsWith("@") ? "FASTQ" : "Raw sequence"
  if (format === "FASTQ") return { format, characterCount: input.length, baseCount: null, state: "Server check", detail: "FASTQ records are validated with their quality data when the analysis runs.", preview: "", invalid: false }
  const bases = raw.split(/\r?\n/).filter((line) => !line.startsWith(">")).join("").replace(/\s+/g, "").toUpperCase()
  const invalidSymbols = [...new Set(bases.replace(/[ACGTURYSWKMBDHVN]/g, ""))]
  if (!bases) return { format, characterCount: input.length, baseCount: 0, state: "No bases", detail: "A FASTA header was detected, but it does not contain sequence data yet.", preview: "", invalid: false }
  if (invalidSymbols.length) return { format, characterCount: input.length, baseCount: bases.length, state: "Check symbols", detail: `Unsupported symbol${invalidSymbols.length > 1 ? "s" : ""}: ${invalidSymbols.join(", ")}.`, preview: bases.slice(0, 42), invalid: true }
  return { format, characterCount: input.length, baseCount: bases.length, state: "Ready", detail: "Client preview is clean. GenomeSight will run full validation and analysis on submission.", preview: bases.slice(0, 42), invalid: false }
}
function readError(cause: unknown, fallback: string) {
  const detail = (cause as { response?: { data?: { detail?: string } } })?.response?.data?.detail
  return typeof detail === "string" ? detail : fallback
}

export default App
