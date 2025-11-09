# Overview

The Genome Sequencing Analyzer is a bioinformatics web application built with Streamlit that provides comprehensive analysis and visualization of genomic sequence data. The application allows users to upload DNA sequence files in FASTA and FASTQ formats, analyze nucleotide composition, calculate GC content, generate quality metrics, and create interactive visualizations. It serves as a powerful tool for researchers, students, and professionals working with genomic data who need quick and accessible sequence analysis capabilities.

# User Preferences

Preferred communication style: Simple, everyday language.

# UI/UX Design System

## Neural Bio-Futurism Theme (Updated November 2025)
The application features a distinctive Neural Bio-Futurism design aesthetic that fuses AI minimalism with biological elegance. This research-grade interface is optimized for extended analysis sessions while maintaining visual sophistication.

### Color Palette
- **Background**: Charcoal Black (#0B0E13) - Deep, neutral base for reduced eye strain
- **Surfaces**: Gunmetal (rgba(22, 27, 34, 0.7)) - Semi-transparent panels with glass-morphism
- **Primary Gradient**: Cyan-Teal (linear-gradient(90deg, #00D8B4, #38BDF8)) - Main accent for interactive elements
- **Typography**: Light Gray (#E2E8F0) for primary text, Slate (#94A3B8) for secondary text
- **Accent**: Magenta-Violet (#C084FC) - Reserved for special highlights and alerts
- **DNA Base Colors**: A=#FF6B6B (red), T=#FFD93D (yellow), C=#6BCB77 (green), G=#4D96FF (blue)

### Visual Effects
- **Glass-morphism**: Semi-transparent surfaces with backdrop-filter blur effects for depth and hierarchy
- **Bioluminescent Glow**: Subtle box-shadow effects using cyan-teal colors on interactive elements
- **Gradient Animations**: Smooth color transitions and hover states for enhanced interactivity
- **Card Elevation**: Multi-layer shadows creating depth and visual separation between content sections

### Design Principles
- **Research-Grade**: Professional interface suitable for academic and clinical use
- **Accessibility**: WCAG-compliant contrast ratios ensuring readability
- **Consistent Theming**: Unified color language throughout all UI components
- **Performance**: CSS-based animations for smooth 60fps interactions

# System Architecture

## Frontend Architecture
The application uses Streamlit as the primary web framework, providing a Python-native approach to building interactive web applications. The frontend follows a single-page application design with custom CSS styling for enhanced user experience. The interface is organized with a wide layout configuration and expandable sidebar for improved usability.

## Core Analysis Engine
The system implements a modular analysis architecture centered around the `SequenceAnalyzer` class, which handles all genomic sequence computations. This design separates concerns by isolating sequence analysis logic from presentation logic, making the codebase maintainable and testable. The analyzer supports both FASTA and FASTQ file formats and provides comprehensive nucleotide analysis capabilities.

## Visualization Layer
Interactive visualizations are handled through a dedicated visualization module using Plotly for creating dynamic, publication-quality charts and graphs. This choice enables rich interactivity, zooming, and data exploration capabilities that static plotting libraries cannot provide. The visualization layer is decoupled from the analysis engine, allowing for easy extension of chart types.

## File Processing Pipeline
The application implements a robust file validation and parsing system that handles multiple sequence file formats. The architecture includes format detection, content validation, and error handling to ensure data integrity throughout the analysis process. This design choice prioritizes data quality and user experience by providing clear feedback on file processing status.

## Data Flow Architecture
The system follows a linear data flow pattern: file upload → validation → parsing → analysis → visualization → report generation. This pipeline approach ensures consistent data processing and makes it easy to add new analysis steps or modify existing ones without affecting other components.

# External Dependencies

## Bioinformatics Libraries
- **Biopython**: Core dependency for biological sequence analysis, file format parsing, and sequence manipulation operations
- **Bio.SeqIO**: Handles reading and writing of sequence files in various formats
- **Bio.SeqUtils**: Provides utility functions for sequence analysis including GC content calculation

## Data Processing and Analysis
- **Pandas**: Used for data manipulation, analysis, and creating structured datasets from sequence analysis results
- **NumPy**: Provides numerical computing capabilities for statistical calculations and array operations on sequence data

## Visualization and Plotting
- **Plotly Express & Graph Objects**: Creates interactive, web-based visualizations for sequence analysis results
- **Plotly Subplots**: Enables creation of complex multi-panel visualizations for comparative analysis

## Web Framework
- **Streamlit**: Primary web application framework that provides the user interface, file upload capabilities, and interactive widgets for parameter configuration

## Standard Library Dependencies
- **io.StringIO & BytesIO**: Handle in-memory file operations for processing uploaded sequence files
- **zipfile**: Manages compressed file archives for batch sequence processing
- **json**: Handles structured data serialization for configuration and report generation
- **datetime**: Provides timestamp functionality for analysis sessions and report generation
- **collections.Counter**: Efficient counting operations for nucleotide composition analysis
- **re**: Regular expression operations for sequence pattern matching and validation