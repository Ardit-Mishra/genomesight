"""
UI Styling Module

This module provides CSS styles and theme configuration for the
Streamlit application. It defines the dark theme color palette
and reusable styling functions.

Color Palette:
    Background: #1A1D29 (dark blue-gray)
    Surface: #25293C (charcoal)
    Primary: #00C9A7 (teal)
    Text Primary: #E8EAED (light gray)
    Text Secondary: #9CA3AF (medium gray)
    Border: #373B4D (subtle gray)
"""


COLORS = {
    'background': '#1A1D29',
    'surface': '#25293C',
    'primary': '#00C9A7',
    'primary_light': '#00FFD1',
    'text_primary': '#E8EAED',
    'text_secondary': '#9CA3AF',
    'border': '#373B4D',
    'error': '#FF6B9D',
    'warning': '#FFC252',
    'success': '#00C9A7'
}

NUCLEOTIDE_COLORS = {
    'A': '#FF6B9D',
    'T': '#FFC252',
    'C': '#00C9A7',
    'G': '#4D96FF',
    'N': '#9CA3AF'
}


def get_dark_theme_css() -> str:
    """
    Get the complete dark theme CSS for the application.
    
    Returns:
        CSS string to be injected via st.markdown
    """
    return """
    <style>
        .stApp {
            background-color: #1A1D29;
            color: #E8EAED;
        }
        
        .main-header {
            font-size: 2rem;
            font-weight: 600;
            color: #FFFFFF;
            margin: 0 0 0.5rem 0;
        }
        
        .subtitle {
            font-size: 1rem;
            color: #9CA3AF;
            margin: 0 0 2rem 0;
        }
        
        .metric-card {
            background: #25293C;
            padding: 1.5rem;
            border-radius: 12px;
            border: 1px solid #373B4D;
            margin: 1rem 0;
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.15);
        }
        
        div[data-testid="stMetric"] {
            background: #25293C;
            border-radius: 8px;
            border: 1px solid #373B4D;
            padding: 1rem;
        }
        
        .stButton > button {
            background: #00C9A7;
            color: #1A1D29;
            border: none;
            border-radius: 6px;
            padding: 0.625rem 1.5rem;
            font-weight: 600;
        }
        
        .stButton > button:hover {
            background: #00FFD1;
        }
        
        button[kind="secondary"] {
            background: transparent;
            color: #E8EAED;
            border: 1px solid #373B4D;
        }
        
        section[data-testid="stSidebar"] {
            background: #1F2230;
            border-right: 1px solid #373B4D;
        }
        
        .stTabs [data-baseweb="tab-list"] {
            gap: 0;
            border-bottom: 1px solid #373B4D;
            background: transparent;
        }
        
        .stTabs [data-baseweb="tab"] {
            background: transparent;
            border: none;
            border-bottom: 2px solid transparent;
            padding: 0.75rem 1.5rem;
            color: #9CA3AF;
            font-weight: 500;
        }
        
        .stTabs [aria-selected="true"] {
            color: #00C9A7;
            border-bottom-color: #00C9A7;
            background: #25293C;
        }
        
        .stSelectbox > div > div,
        .stTextInput > div > div > input,
        .stNumberInput > div > div > input {
            background: #25293C;
            border: 1px solid #373B4D;
            border-radius: 6px;
            color: #E8EAED;
        }
        
        ::-webkit-scrollbar {
            width: 8px;
            height: 8px;
        }
        
        ::-webkit-scrollbar-track {
            background: #1A1D29;
        }
        
        ::-webkit-scrollbar-thumb {
            background: #373B4D;
            border-radius: 4px;
        }
    </style>
    """


def apply_dark_theme(st_module):
    """
    Apply the dark theme to a Streamlit app.
    
    Args:
        st_module: The streamlit module (import streamlit as st)
    """
    st_module.markdown(get_dark_theme_css(), unsafe_allow_html=True)


def get_nucleotide_colors() -> dict:
    """Return nucleotide color mapping."""
    return NUCLEOTIDE_COLORS.copy()


def get_color_palette() -> dict:
    """Return the complete color palette."""
    return COLORS.copy()


def format_sequence_html(sequence: str, max_length: int = 100) -> str:
    """
    Format a sequence with colored nucleotides.
    
    Args:
        sequence: DNA/RNA sequence
        max_length: Maximum length to display
    
    Returns:
        HTML formatted string with colored bases
    """
    if len(sequence) > max_length:
        sequence = sequence[:max_length] + '...'
    
    colored = []
    for base in sequence.upper():
        color = NUCLEOTIDE_COLORS.get(base, COLORS['text_secondary'])
        colored.append(f'<span style="color:{color}">{base}</span>')
    
    return ''.join(colored)
