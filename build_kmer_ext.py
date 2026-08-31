"""
Build script for the optional native (Cython) k-mer accelerator.

This extension is NOT required to run GenomeSight. app/core/kmer_native.py
falls back to a pure-Python k-mer counter with identical output when the
compiled extension isn't present -- and the UI announces plainly which
path is active, every time. Build this only if you want the faster path.

Requires a C compiler on the host (MSVC on Windows / gcc or clang on
Linux/macOS) plus Cython and setuptools, which are NOT part of the app's
normal runtime dependencies -- they're only needed at build time, e.g.:

    uv run --python 3.12 --with cython --with setuptools \\
        python build_kmer_ext.py build_ext --inplace

or, with Cython + a C compiler already on PATH:

    python build_kmer_ext.py build_ext --inplace

On success this produces a compiled module next to
app/core/_kmer_accel.pyx (e.g. app/core/_kmer_accel.cp312-win_amd64.pyd
on Windows, or app/core/_kmer_accel.cpython-312-x86_64-linux-gnu.so on
Linux). app/core/kmer_native.py picks it up automatically on next import --
no other wiring required.

If the build fails (no compiler, no Cython, unsupported platform), that's
fine: the app still runs correctly on the pure-Python path, and says so.
"""
from setuptools import setup
from Cython.Build import cythonize

setup(
    name="genomesight-kmer-accel",
    ext_modules=cythonize(
        "app/core/_kmer_accel.pyx",
        compiler_directives={"language_level": "3"},
    ),
)
