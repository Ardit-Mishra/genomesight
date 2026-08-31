"""
Honest benchmark for the k-mer counting paths.

Times the pure-Python counter (app.core.kmer_native._count_kmers_python)
and, when the compiled extension is available, the native Cython counter
(app.core._kmer_accel.count_kmers), on the SAME generated input for each
(sequence_length, k) combination. Writes results to
benchmarks/kmer_benchmark_results.json.

This script does not assume the native path is faster -- it measures both
paths and reports whatever the numbers say, including "native not
available in this environment" when that's the truth. Do not read
fabricated or hand-edited numbers into the results file; regenerate it by
running this script.

Usage:
    uv run --python 3.12 python benchmarks/benchmark_kmer.py
    uv run --python 3.12 python benchmarks/benchmark_kmer.py --repeats 5
"""
import argparse
import json
import platform
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.core import kmer_native  # noqa: E402

SEQUENCE_LENGTHS = [1_000, 10_000, 100_000, 500_000]
K_VALUES = [4, 8, 12]
ALPHABET = "ATCG"


def make_sequence(length: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice(ALPHABET) for _ in range(length))


def time_call(fn, *args, repeats: int) -> dict:
    times = []
    result = None
    for _ in range(repeats):
        start = time.perf_counter()
        result = fn(*args)
        times.append(time.perf_counter() - start)
    return {
        "seconds_min": min(times),
        "seconds_mean": sum(times) / len(times),
        "repeats": repeats,
        "result_unique_kmers": len(result) if result is not None else None,
    }


def run(repeats: int) -> dict:
    status = kmer_native.backend_status()
    results = {
        "environment": {
            "python_version": sys.version,
            "platform": platform.platform(),
            "native_available": status["native_available"],
            "native_import_error": status["import_error"],
        },
        "runs": [],
    }

    for length in SEQUENCE_LENGTHS:
        sequence = make_sequence(length, seed=length)
        for k in K_VALUES:
            entry = {"sequence_length": length, "k": k}

            entry["python"] = time_call(
                kmer_native._count_kmers_python, [sequence], k, repeats=repeats
            )

            if status["native_available"]:
                entry["native"] = time_call(
                    kmer_native._kmer_accel.count_kmers, [sequence], k, repeats=repeats
                )
                if entry["native"]["result_unique_kmers"] != entry["python"]["result_unique_kmers"]:
                    raise RuntimeError(
                        f"native/python disagreed on unique k-mer count at "
                        f"length={length} k={k} -- refusing to report a benchmark "
                        f"for implementations that don't agree"
                    )
                py_time = entry["python"]["seconds_min"]
                native_time = entry["native"]["seconds_min"]
                entry["speedup_x"] = (py_time / native_time) if native_time > 0 else None
            else:
                entry["native"] = None
                entry["speedup_x"] = None

            results["runs"].append(entry)
            native_note = (
                f"{entry['speedup_x']:.2f}x" if entry.get("speedup_x") else "N/A (native not built)"
            )
            print(
                f"length={length:>7}  k={k:>2}  "
                f"python={entry['python']['seconds_min']*1000:8.2f}ms  "
                f"speedup={native_note}"
            )

    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=3, help="timed repeats per case")
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parent / "kmer_benchmark_results.json",
    )
    args = parser.parse_args()

    results = run(args.repeats)

    args.out.write_text(json.dumps(results, indent=2))
    print(f"\nWrote {args.out}")

    if not results["environment"]["native_available"]:
        print(
            "\nNOTE: native accelerator was NOT built in this environment "
            f"({results['environment']['native_import_error']}). Only the "
            "pure-Python path was timed. Build it (see README.md) and rerun "
            "for a native-vs-Python comparison."
        )


if __name__ == "__main__":
    main()
