// GENERATED from design-system/icons.json — do not edit. Run design-system/build.py.
//
// Domain glyphs that no icon library ships: a peptide-binding groove, an ORF
// reading frame, a reliability diagram. Drawn to Lucide's grammar (24px box,
// 1.5 stroke, currentColor, round joins) so the two sets mix without seams —
// keep using lucide-react for generic affordances (upload, search, chevron).

import type { SVGProps } from "react";

export interface IconProps extends SVGProps<SVGSVGElement> {
  /** Rendered edge length in px. The 24-unit geometry scales to it. */
  size?: number;
}

function Frame({ size = 20, children, ...rest }: IconProps) {
  return (
    <svg
      viewBox="0 0 24 24"
      width={size}
      height={size}
      fill="none"
      stroke="currentColor"
      strokeWidth={1.5}
      strokeLinecap="round"
      strokeLinejoin="round"
      aria-hidden={rest["aria-label"] ? undefined : true}
      focusable="false"
      {...rest}
    >
      {children}
    </svg>
  );
}

/** double-stranded nucleic acid. GenomeSight's primary mark and the marker for anything operating on DNA/RNA. */
export function Helix(props: IconProps) {
  return <Frame {...props}><path d="M7 2c0 5 10 5 10 10s-10 5-10 10"/><path d="M17 2c0 5-10 5-10 10s10 5 10 10"/><path d="M9.2 5.4h5.6"/><path d="M8.3 12h7.4"/><path d="M9.2 18.6h5.6"/></Frame>;
}

/** an amino-acid chain — beads on a backbone. Marks anything whose subject is a peptide sequence rather than a nucleotide one. */
export function Peptide(props: IconProps) {
  return <Frame {...props}><circle cx="4" cy="15" r="2"/><circle cx="9.3" cy="10.5" r="2"/><circle cx="14.7" cy="15" r="2"/><circle cx="20" cy="10.5" r="2"/><path d="M5.5 13.7 7.8 11.8"/><path d="M10.8 11.8 13.2 13.7"/><path d="M16.2 13.7 18.5 11.8"/></Frame>;
}

/** PeptideMHC's signature mark. Two alpha-helices (the wavy rails) with a peptide lying in the groove between them — the actual structure the model predicts. The helices are outlined because they are the receptor; the peptide residues are filled because they are what is being measured. */
export function MhcGroove(props: IconProps) {
  return <Frame {...props}><path d="M2 7q2.5-3 5 0t5 0t5 0t5 0"/><path d="M2 17q2.5 3 5 0t5 0t5 0t5 0"/><path d="M6 12h12"/><circle cx="8" cy="12" r="1.3" fill="currentColor"/><circle cx="12" cy="12" r="1.3" fill="currentColor"/><circle cx="16" cy="12" r="1.3" fill="currentColor"/></Frame>;
}

/** aromatic ring. BioStudio's primary mark and the marker for small-molecule / cheminformatics work. */
export function Benzene(props: IconProps) {
  return <Frame {...props}><path d="M12 3.5 19.4 7.75v8.5L12 20.5 4.6 16.25v-8.5z"/><circle cx="12" cy="12" r="4"/></Frame>;
}

/** three atoms and their bonds, one of them a double bond. Marks structure input (SMILES), structure rendering, and descriptor calculation. */
export function Molecule(props: IconProps) {
  return <Frame {...props}><circle cx="5.5" cy="8" r="2.3"/><circle cx="18.5" cy="8" r="2.3"/><circle cx="12" cy="17.5" r="2.3"/><path d="M7.8 6.9h8.4"/><path d="M7.8 9.1h8.4"/><path d="M6.8 9.9 10.7 15.6"/><path d="M17.2 9.9 13.3 15.6"/></Frame>;
}

/** two sequences compared position by position. The vertical pipes are matches; the x is a mismatch — the icon states what alignment produces rather than just showing two lines. */
export function Alignment(props: IconProps) {
  return <Frame {...props}><path d="M3 6h18"/><path d="M3 18h18"/><path d="M7 9v6"/><path d="M11 9v6"/><path d="M19 9v6"/><path d="M14 10.5 16 13.5"/><path d="M16 10.5 14 13.5"/></Frame>;
}

/** an open reading frame: start bracket, directional translation, stop bracket. Marks ORF detection and anything frame-aware. */
export function ReadingFrame(props: IconProps) {
  return <Frame {...props}><path d="M4 6v12"/><path d="M20 6v12"/><path d="M7 12h9"/><path d="M13 9l3 3-3 3"/></Frame>;
}

/** base composition / GC content — four unequal columns over a baseline. The four columns are literal: A, C, G, T. */
export function Composition(props: IconProps) {
  return <Frame {...props}><path d="M3 20h18"/><path d="M6 20v-6"/><path d="M10.5 20V7"/><path d="M15 20v-9"/><path d="M19.5 20v-4"/></Frame>;
}

/** a fixed-width window sliding along a sequence. Marks k-mer frequency analysis; the boxed span is the k. */
export function Kmer(props: IconProps) {
  return <Frame {...props}><path d="M2 12h6"/><path d="M16 12h6"/><path d="M8 7h8v10H8z"/><path d="M10.5 10v4"/><path d="M13.5 10v4"/></Frame>;
}

/** a short pattern recurring at several positions along a sequence. Marks motif search; the diamonds are hits, the line is the sequence. */
export function Motif(props: IconProps) {
  return <Frame {...props}><path d="M3 16h18"/><path d="M6 8l1.8 2L6 12l-1.8-2z"/><path d="M12 8l1.8 2L12 12l-1.8-2z"/><path d="M18 8l1.8 2L18 12l-1.8-2z"/></Frame>;
}

/** a dataset partitioned into a training majority and a held-out remainder. Marks every split disclosure — scaffold split, peptide-grouped split, leave-one-allele-out. The filled block is what the model saw; the outlined block is what it was scored on; the rule between them is the partition, and nothing crosses it. */
export function Split(props: IconProps) {
  return <Frame {...props}><rect x="2.5" y="8.5" width="11.5" height="7" rx="1.5" fill="currentColor"/><rect x="17.5" y="8.5" width="4" height="7" rx="1.5"/><path d="M15.75 4.5v15"/></Frame>;
}

/** a reliability diagram: the stepped line is the ideal y=x reference, the curve is observed frequency. Marks calibration reporting — the one place a probability is allowed to be called trustworthy. */
export function Calibration(props: IconProps) {
  return <Frame {...props}><path d="M4 3v18h17"/><path d="M6.5 18.5 8.5 16.5"/><path d="M11 14 13 12"/><path d="M15.5 9.5 17.5 7.5"/><path d="M5 19.5c5 .5 5.5-6 14-13"/></Frame>;
}

/** a microplate — many specimens measured in one run. Marks batch prediction and any many-at-once operation. */
export function Plate(props: IconProps) {
  return <Frame {...props}><rect x="3" y="5" width="18" height="14" rx="2"/><circle cx="7.5" cy="9.5" r="1.5"/><circle cx="12" cy="9.5" r="1.5"/><circle cx="16.5" cy="9.5" r="1.5"/><circle cx="7.5" cy="14.5" r="1.5"/><circle cx="12" cy="14.5" r="1.5"/><circle cx="16.5" cy="14.5" r="1.5"/></Frame>;
}

/** a lipid bilayer with something crossing it. Marks permeability endpoints — blood-brain barrier, Caco-2, absorption. */
export function Membrane(props: IconProps) {
  return <Frame {...props}><path d="M3 9h18"/><path d="M3 15h18"/><path d="M12 3v18"/><path d="M9 18l3 3 3-3"/></Frame>;
}

/** biotransformation — substrate cycling to product. Marks CYP450 endpoints, clearance and half-life. */
export function Metabolism(props: IconProps) {
  return <Frame {...props}><path d="M20.5 12a8.5 8.5 0 1 1-2.9-6.4"/><path d="M21 4v4.5h-4.5"/><circle cx="12" cy="12" r="2.2"/></Frame>;
}

/** a toxicity or high-risk endpoint. Uses the caveat colour, never the accent — risk is a property of the compound, not a state of the control. */
export function Hazard(props: IconProps) {
  return <Frame {...props}><path d="M12 3.5 21.5 20H2.5z"/><path d="M12 10v4.5"/><circle cx="12" cy="17.3" r="0.9" fill="currentColor"/></Frame>;
}

/** a ligand seated in a binding pocket. Marks target prediction and protein-ligand work. Pocket outlined, ligand filled — same encoding rule as the groove. */
export function Receptor(props: IconProps) {
  return <Frame {...props}><path d="M5 7v6a7 7 0 0 0 14 0V7"/><circle cx="12" cy="10.5" r="2.6" fill="currentColor"/><path d="M12 20v2"/></Frame>;
}

/** a knowledge graph — entities and the relations between them. Marks the drug-target-disease graph and any network view. */
export function Graph(props: IconProps) {
  return <Frame {...props}><circle cx="10" cy="12" r="2.8"/><circle cx="4" cy="4.5" r="1.7"/><circle cx="19.5" cy="6" r="1.7"/><circle cx="20" cy="16.5" r="1.7"/><circle cx="7" cy="20" r="1.7"/><path d="M8.3 9.8 5.1 5.8"/><path d="M12.4 10.5 18.1 6.9"/><path d="M12.6 13.1 18.5 15.8"/><path d="M9 14.6 7.6 18.4"/><path d="M19.7 7.7 19.9 14.8"/></Frame>;
}

/** an assay or wet-lab measurement — the origin of every training label in these apps. Marks datasets, experimental sources and provenance. */
export function Flask(props: IconProps) {
  return <Frame {...props}><path d="M9 3v6.2L3.9 18.6A2 2 0 0 0 5.6 21.5h12.8a2 2 0 0 0 1.7-2.9L15 9.2V3"/><path d="M8 3h8"/><path d="M6.6 15h10.8"/></Frame>;
}

/** one residue substituted in a chain. Marks saturation mutagenesis, variant effects and single-position sensitivity. The substituted residue is filled; its neighbours are not. */
export function Mutation(props: IconProps) {
  return <Frame {...props}><path d="M3 15h18"/><circle cx="6" cy="15" r="2"/><circle cx="12" cy="15" r="2" fill="currentColor"/><circle cx="18" cy="15" r="2"/><path d="M12 3v7"/><path d="M9.5 5.5 12 3l2.5 2.5"/></Frame>;
}

/** a FASTA/FASTQ record — a file whose content is a sequence, not prose. Marks upload, export and file-format disclosure. */
export function SequenceFile(props: IconProps) {
  return <Frame {...props}><path d="M14 2.5H7a2 2 0 0 0-2 2v15a2 2 0 0 0 2 2h10a2 2 0 0 0 2-2V7.5z"/><path d="M14 2.5v5h5"/><path d="M8.5 12.5h4"/><path d="M8.5 15.5h7"/><path d="M8.5 18.5h5"/></Frame>;
}

/** a measurement scale. Marks method/metric disclosure — the sections that say how a number was produced and on what. */
export function Ruler(props: IconProps) {
  return <Frame {...props}><path d="M2.5 9h19v6h-19z"/><path d="M6.5 9v3"/><path d="M10 9v4"/><path d="M13.5 9v3"/><path d="M17 9v4"/></Frame>;
}

export type InstrumentIconName = "helix" | "peptide" | "mhcGroove" | "benzene" | "molecule" | "alignment" | "readingFrame" | "composition" | "kmer" | "motif" | "split" | "calibration" | "plate" | "membrane" | "metabolism" | "hazard" | "receptor" | "graph" | "flask" | "mutation" | "sequenceFile" | "ruler";

/** Lookup for cases where the glyph is chosen at runtime (nav tables, maps). */
export const instrumentIcons: Record<InstrumentIconName, (p: IconProps) => JSX.Element> = {
  helix: Helix,
  peptide: Peptide,
  mhcGroove: MhcGroove,
  benzene: Benzene,
  molecule: Molecule,
  alignment: Alignment,
  readingFrame: ReadingFrame,
  composition: Composition,
  kmer: Kmer,
  motif: Motif,
  split: Split,
  calibration: Calibration,
  plate: Plate,
  membrane: Membrane,
  metabolism: Metabolism,
  hazard: Hazard,
  receptor: Receptor,
  graph: Graph,
  flask: Flask,
  mutation: Mutation,
  sequenceFile: SequenceFile,
  ruler: Ruler,
};

export const instrumentIconNames = ["helix", "peptide", "mhcGroove", "benzene", "molecule", "alignment", "readingFrame", "composition", "kmer", "motif", "split", "calibration", "plate", "membrane", "metabolism", "hazard", "receptor", "graph", "flask", "mutation", "sequenceFile", "ruler"] as const;
