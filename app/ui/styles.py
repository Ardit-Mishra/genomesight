"""
UI Styling Module

CSS and theme configuration for the Streamlit application.

The palette is NOT defined here — it comes from `tokens.py`, which is generated
from `live-projects/design-system/tokens.json` and shared with BioStudio and
PeptideMHC. Edit the JSON and run `python design-system/build.py`; never edit
tokens.py or re-declare a colour in this file, which is how the charts and the
UI drifted apart in the first place.
"""

from .tokens import (  # generated — see design-system/build.py
    COLORS as _DS,
    NUCLEOTIDE_COLORS,
    FONT_UI,
    FONT_MONO,
)

#: Legacy key names kept so existing call sites (main.py, plots.py, components.py)
#: keep working, mapped onto the shared token names.
COLORS = {
    'background': _DS['ground'],
    'surface': _DS['surface'],
    'primary': _DS['accent'],
    'primary_light': _DS['accentInk'],
    'text_primary': _DS['ink'],
    'text_secondary': _DS['inkMuted'],
    'border': _DS['border'],
    'error': _DS['critical'],
    'warning': _DS['caveat'],
    'success': _DS['safe'],
}


def get_dark_theme_css() -> str:
    """
    Complete "Laboratory Instrument" CSS for the application.

    Derived from COLORS rather than hardcoded, so the charts (which import the
    same dict) and the surrounding UI cannot drift apart — they previously did.

    Returns:
        CSS string to be injected via st.markdown
    """
    return f"""
    <style>
        @import url('https://fonts.googleapis.com/css2?family=Bricolage+Grotesque:opsz,wght@12..96,400;12..96,500;12..96,600;12..96,700&family=IBM+Plex+Mono:wght@400;500;600&display=swap');

        .stApp {{
            background-color: {COLORS['background']};
            color: {COLORS['text_primary']};
        }}

        /* Streamlit ships generated emotion classes whose specificity beats a
           plain element selector, so the UI face has to be declared !important
           to land at all. The data face is re-asserted the same way below. */
        html, body, .stApp, [data-testid="stAppViewContainer"], section[data-testid="stSidebar"],
        h1, h2, h3, h4, h5, h6, p, div, span, label, button, input, select, textarea {{
            font-family: {FONT_UI} !important;
        }}

        /* Left-aligned and restrained; the sequence is the subject, not the
           title. Streamlit's own h1 rule sets 44px, hence !important. */
        .main-header {{
            font-family: {FONT_UI} !important;
            font-size: 1.6rem !important;
            font-weight: 600 !important;
            letter-spacing: -0.02em;
            color: {COLORS['text_primary']};
            margin: 0 0 0.25rem 0 !important;
            padding: 0 !important;
        }}

        .subtitle {{
            font-size: 0.95rem;
            color: {COLORS['text_secondary']};
            margin: 0 0 2rem 0;
        }}

        /* Hairline panels. Elevation is a border, not a drop shadow. */
        .metric-card {{
            background: {COLORS['surface']};
            padding: 1.15rem 1.25rem;
            border-radius: 6px;
            border: 1px solid {COLORS['border']};
            margin: 1rem 0;
        }}

        div[data-testid="stMetric"] {{
            background: {COLORS['surface']};
            border-radius: 6px;
            border: 1px solid {COLORS['border']};
            padding: 0.9rem 1rem;
        }}

        /* Counts, lengths and percentages are read in columns, so they are
           monospace and tabular — digits keep a fixed width as values change. */
        div[data-testid="stMetricValue"] {{
            font-family: {FONT_MONO};
            font-variant-numeric: tabular-nums;
            letter-spacing: -0.02em;
        }}

        div[data-testid="stMetricLabel"] {{
            font-family: {FONT_MONO};
            font-size: 0.68rem;
            letter-spacing: 0.12em;
            text-transform: uppercase;
            color: {COLORS['text_secondary']};
        }}

        .stButton > button {{
            background: {COLORS['primary']};
            color: {COLORS['background']};
            border: 1px solid {COLORS['primary']};
            border-radius: 6px;
            padding: 0.5rem 1.1rem;
            font-family: {FONT_UI};
            font-weight: 600;
            transition: filter .15s ease;
        }}

        .stButton > button:hover {{
            filter: brightness(1.08);
        }}

        button[kind="secondary"] {{
            background: transparent;
            color: {COLORS['text_primary']};
            border: 1px solid {COLORS['border']};
        }}

        button[kind="secondary"]:hover {{
            border-color: {COLORS['primary']};
        }}

        section[data-testid="stSidebar"] {{
            background: {COLORS['surface']};
            border-right: 1px solid {COLORS['border']};
        }}

        /* Sidebar section headings read as instrument labels. */
        section[data-testid="stSidebar"] h4 {{
            font-family: {FONT_MONO};
            font-size: 0.68rem;
            font-weight: 500;
            letter-spacing: 0.12em;
            text-transform: uppercase;
            color: {COLORS['text_secondary']};
        }}

        .stTabs [data-baseweb="tab-list"] {{
            gap: 0;
            border-bottom: 1px solid {COLORS['border']};
            background: transparent;
        }}

        .stTabs [data-baseweb="tab"] {{
            background: transparent;
            border: none;
            border-bottom: 2px solid transparent;
            padding: 0.7rem 1.25rem;
            color: {COLORS['text_secondary']};
            font-weight: 500;
        }}

        .stTabs [aria-selected="true"] {{
            color: {COLORS['primary']};
            border-bottom-color: {COLORS['primary']};
            background: transparent;
        }}

        .stSelectbox > div > div,
        .stTextInput > div > div > input,
        .stNumberInput > div > div > input {{
            background: {COLORS['background']};
            border: 1px solid {COLORS['border']};
            border-radius: 6px;
            color: {COLORS['text_primary']};
        }}

        .stTextInput > div > div > input:focus {{
            border-color: {COLORS['primary']};
        }}

        /* Streamlit paints alerts as a full accent-tinted block with accent body
           text. Neutral panel + one 2px status rail keeps the accent reserved
           for measurement and interaction. */
        [data-testid="stAlert"] {{
            background: {COLORS['surface']} !important;
            border: 1px solid {COLORS['border']} !important;
            border-left: 2px solid {COLORS['primary']} !important;
            border-radius: 0 6px 6px 0 !important;
            box-shadow: none !important;
        }}
        [data-testid="stAlert"] * {{ color: {COLORS['text_primary']} !important; }}
        [data-testid="stAlert"] p {{ max-width: 78ch; }}

        /* Sequences, alignments and any code block: fixed width, always.
           !important because the broad UI-face rule above also matches these. */
        code, pre, .stCode, div[data-testid="stDataFrame"],
        div[data-testid="stMetricValue"], div[data-testid="stMetricLabel"] {{
            font-family: {FONT_MONO} !important;
        }}

        /* An alignment's match line only lines up if nothing wraps or re-spaces. */
        .stCode pre {{
            white-space: pre;
        }}

        ::-webkit-scrollbar {{
            width: 8px;
            height: 8px;
        }}

        ::-webkit-scrollbar-track {{
            background: {COLORS['background']};
        }}

        ::-webkit-scrollbar-thumb {{
            background: {COLORS['border']};
            border-radius: 4px;
        }}

        @media (prefers-reduced-motion: reduce) {{
            * {{ transition-duration: .01ms !important; animation-duration: .01ms !important; }}
        }}

        /* Streamlit's st.columns() keeps a fixed row at every viewport width,
           so a 3- or 4-up metric row (basic stats, ORF summary, alignment
           identity/score/gaps/length) truncates its labels once the sidebar
           (~336px, always present) leaves the main column under ~700px —
           i.e. below roughly a 1100px viewport. Wrapping to a 2-up grid
           there keeps every label and value on-screen without touching the
           desktop layout, where the main column stays well above that. */
        @media (max-width: 1100px) {{
            div[data-testid="stHorizontalBlock"] {{
                flex-wrap: wrap;
                row-gap: 0.75rem;
            }}
            div[data-testid="stHorizontalBlock"] > div[data-testid="stColumn"] {{
                min-width: min(220px, 100%) !important;
                flex: 1 1 45% !important;
                width: auto !important;
            }}
        }}
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
