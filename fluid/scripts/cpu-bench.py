#!/usr/bin/env python3
"""Alternate two release CPU probes; verify every result hash."""
import argparse
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--repeats", type=int, default=4)
    parser.add_argument("--case", action="append", metavar="NAME:N:COUNT",
                        help="override default cpu_validate workloads; repeat for multiple cases")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.threads < 1 or args.repeats < 2:
        parser.error("threads must be positive and repeats at least two")
    binaries = {"before": args.baseline.resolve(), "after": args.candidate.resolve()}
    for binary in binaries.values():
        if not binary.is_file():
            parser.error(f"missing executable: {binary}")
    environment = dict(os.environ, RAYON_NUM_THREADS=str(args.threads))
    cases = [("smooth", 32, 100), ("smooth", 48, 60),
             ("interface", 24, 60), ("interface", 32, 40),
             ("open", 24, 60), ("open", 32, 40)]
    if args.case:
        cases = []
        for value in args.case:
            try:
                case, n, count = value.split(":")
                n, count = int(n), int(count)
                if not case or n < 1 or count < 1:
                    raise ValueError
            except ValueError:
                parser.error("case must be NAME:N:COUNT with positive N and COUNT")
            cases.append((case, n, count))
    runs = []
    summaries = []
    for case, n, steps in cases:
        times = {v: [] for v in binaries}
        hashes = set()
        for repeat in range(args.repeats):
            order = ("before", "after") if repeat % 2 == 0 else ("after", "before")
            for version in order:
                command = [str(binaries[version]), case, str(n), str(steps)]
                result = subprocess.run(command, env=environment, text=True,
                                        capture_output=True, check=True)
                output = result.stdout.strip()
                values = dict(field.split("=", 1) for field in output.split() if "=" in field)
                times[version].append(float(values["wall_s"]))
                hashes.add(values["hash"])
                row = dict(case=case, n=n, steps=steps, repeat=repeat,
                           version=version, command=command, output=output)
                runs.append(row)
                print(output + f" [{version}]", flush=True)
        if len(hashes) != 1:
            raise RuntimeError(f"result mismatch for {case} {n}: {hashes}")
        median = {v: statistics.median(t) for v, t in times.items()}
        summary = dict(case=case, n=n, steps=steps, median_seconds=median,
                       ranges_seconds={v: [min(t), max(t)] for v, t in times.items()},
                       speedup=median["before"] / median["after"],
                       state_hash=hashes.pop())
        summaries.append(summary)
        print(json.dumps(summary), flush=True)
    data = dict(host=platform.platform(), threads=args.threads,
                repeats=args.repeats, summaries=summaries, runs=runs)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2) + "\n")


if __name__ == "__main__":
    main()
