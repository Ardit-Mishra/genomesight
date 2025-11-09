# Overview

The Genome Sequencing Analyzer is a bioinformatics web application built with Streamlit that provides comprehensive analysis and visualization of genomic sequence data. The application allows users to upload DNA sequence files in FASTA and FASTQ formats, analyze nucleotide composition, calculate GC content, generate quality metrics, and create interactive visualizations. It serves as a powerful tool for researchers, students, and professionals working with genomic data who need quick and accessible sequence analysis capabilities.

# User Preferences

Preferred communication style: Simple, everyday language.

# UI/UX Design System

## Classic Scientific Dark Theme (Updated November 2025)
The application features a professional dark theme optimized for extended analysis sessions. The design emphasizes clarity and functionality while providing comprehensive feature information.

### Color Palette
- **Background**: Dark Blue-Gray (#1A1D29) - Reduces eye strain during long sessions
- **Surfaces**: Charcoal (#25293C) - Distinct card backgrounds
- **Primary Accent**: Teal (#00C9A7) - Professional scientific accent
- **Typography**: Light Gray (#E8EAED) for primary text, Medium Gray (#9CA3AF) for secondary text
- **Borders**: Subtle Gray (#373B4D) - Minimal separation
- **DNA Base Colors**: A=#FF6B9D (pink), T=#FFC252 (amber), C=#00C9A7 (teal), G=#4D96FF (blue)

### Visual Style
- **Information-Rich Cards**: Dark backgrounds with comprehensive feature descriptions
- **Clear Typography**: Light text on dark backgrounds for comfortable reading
- **Prominent Upload**: File upload positioned in main content area for easy access
- **Tabbed Information**: Four organized tabs (Overview, Features & Capabilities, File Formats, Use Cases)
- **Feature Showcase**: Detailed explanations of analysis capabilities organized logically

### Design Principles
- **Informative**: Clear descriptions of what each analysis does and why it's useful
- **Professional**: Classic scientific dark theme suitable for research environments
- **User-Centric**: File upload prominently placed; sidebar for analysis options
- **Educational**: Includes use cases (research, clinical, education) to help users understand applications

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