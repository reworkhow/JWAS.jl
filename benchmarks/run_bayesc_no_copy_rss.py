"""Isolated production BayesC RSS measurements and strictly keyed comparisons.

Default: 36 fresh Julia processes, two workers, 200 outer sweeps/50 burn-in.
Raw exports/logs and complete aligned tables are retained locally and ignored in
Git. Compact metrics, comparisons and inspected joined samples are reviewable.
"""
import argparse
import concurrent.futures
import csv
from decimal import Decimal, InvalidOperation
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def unquote(value):
    value = str(value).strip()
    while len(value) >= 2 and value[0] == value[-1] and value[0] in "\"'":
        value = value[1:-1].strip()
    return value


def marker_key(value, p):
    value = unquote(value)
    # Numeric forms denote explicit original SNP indices, not CSV row order.
    try:
        numeric = Decimal(value)
        if numeric.is_finite() and numeric == numeric.to_integral_value():
            index = int(numeric)
            if 1 <= index <= p:
                return f"rs_{1000 + 7 * index}"
    except InvalidOperation:
        pass
    return value


def keyed(rows, kind, n, p):
    expected = ({("y", f"animal_{i:06d}") for i in range(1, n + 1)} if kind == "animals"
                else {("y", f"rs_{1000 + 7 * i}") for i in range(1, p + 1)})
    assert len(rows) == len(expected), (kind, "unexpected row count", len(rows), len(expected))
    result = {}
    for row in rows:
        identifier = unquote(row["ID"]) if kind == "animals" else marker_key(row["Marker_ID"], p)
        key = (unquote(row["Trait"]), identifier)
        assert key not in result, (kind, "duplicate key", key)
        result[key] = row
    assert set(result) == expected, (kind, "missing or unexpected keys", set(result) ^ expected)
    return result, expected


def align(left, right, kind, n, p):
    a, expected = keyed(left, kind, n, p)
    b, right_expected = keyed(right, kind, n, p)
    assert set(a) == set(b) == expected == right_expected
    columns = ("EBV", "PEV") if kind == "animals" else ("Estimate", "SD", "Model_Frequency")
    identifier = "ID" if kind == "animals" else "Marker_ID"
    rows = []
    for key in sorted(expected):  # Keys have already joined and validated above.
        row = {"Trait": key[0], identifier: key[1]}
        for column in columns:
            u, v = float(a[key][column]), float(b[key][column])
            assert math.isfinite(u) and math.isfinite(v), (kind, key, column)
            row[f"adopted_{column}"] = u
            row[f"comparison_{column}"] = v
        rows.append(row)
    assert len(rows) == len(expected) == len(a) == len(b), "Join lost or duplicated rows"
    summary = {}
    for column in columns:
        deltas = [abs(row[f"adopted_{column}"] - row[f"comparison_{column}"]) for row in rows]
        summary[f"max_abs_{column}"] = max(deltas)
        summary[f"match_{column}"] = all(math.isclose(row[f"adopted_{column}"], row[f"comparison_{column}"],
                                                       abs_tol=1e-4, rel_tol=1e-4) for row in rows)
    return rows, summary


def validate_controls(metric, args):
    expected = {"outer": args.outer, "burnin": args.burnin, "block_size": args.block_size,
                "saved_draws": (args.outer - args.burnin) // 10}
    for key, value in expected.items():
        assert int(metric[key]) == value, ("Stored measurement uses different controls", key, metric[key], value)


def self_test():
    # Reverse order, quoted labels and numeric indices must still align. Duplicate,
    # missing and incorrect marker keys must be rejected before any summary.
    marker = [{"Trait": "y", "Marker_ID": "rs_1007", "Estimate": "0.1", "SD": "0.2", "Model_Frequency": "0.3"},
              {"Trait": "y", "Marker_ID": "rs_1014", "Estimate": "0.4", "SD": "0.5", "Model_Frequency": "0.6"}]
    renamed = [dict(marker[1], Marker_ID='"2.0"'), dict(marker[0], Marker_ID="'1'")]
    rows, summary = align(marker, renamed, "markers", 2, 2)
    assert len(rows) == 2 and summary["max_abs_Estimate"] == 0
    for invalid in (renamed[:1], [renamed[0], renamed[0]], [dict(renamed[0], Marker_ID="other"), renamed[1]]):
        try:
            align(marker, invalid, "markers", 2, 2)
        except AssertionError:
            continue
        raise AssertionError("Invalid comparison was accepted")
    controls = argparse.Namespace(outer=200, burnin=50, block_size=64)
    metric = {"outer": "200", "burnin": "50", "block_size": "64", "saved_draws": "15"}
    validate_controls(metric, controls)
    for key in metric:
        try:
            validate_controls(dict(metric, **{key: "1"}), controls)
        except AssertionError:
            continue
        raise AssertionError("Mismatched measurement controls were accepted")
    print("Comparison checks passed: reordered/numeric/quoted IDs align; dropped, duplicated, wrong keys and mismatched controls fail.")


def source_fingerprint(project):
    digest = hashlib.sha256()
    paths = [project / "Project.toml", *sorted((project / "src").rglob("*.jl")),
             Path(__file__), Path(__file__).with_name("bayesc_no_copy_rss_case.jl")]
    for path in paths:
        digest.update(str(path.relative_to(project)).encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def run_one(config):
    case_dir, args, project, report, fingerprint = config
    metrics_path = case_dir / "metrics.csv"
    fingerprint_path = case_dir / "fingerprint.txt"
    if metrics_path.exists() and fingerprint_path.exists() and fingerprint_path.read_text().strip() == fingerprint:
        previous = read_csv(metrics_path)[0]
        try:
            validate_controls(previous, args)
        except AssertionError:
            pass
        else:
            return previous
    if case_dir.exists():
        shutil.rmtree(case_dir)  # Only this runner's generated raw case directory.
    case_dir.mkdir(parents=True, exist_ok=True)
    log = report / "logs" / (case_dir.name + ".txt")
    log.parent.mkdir(parents=True, exist_ok=True)
    command = [args.julia, f"--project={project}", "--startup-file=no", "--compiled-modules=existing",
               str(Path(__file__).with_name("bayesc_no_copy_rss_case.jl")), str(case_dir)]
    n, p, source, qc, seed = config_values(case_dir.name)
    command += [str(n), str(p), source, str(qc).lower(), str(seed), str(args.outer), str(args.burnin), str(args.block_size)]
    environment = dict(os.environ, JULIA_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", GKSwstype="100")
    start = time.monotonic()
    print("Starting", case_dir.name, flush=True)
    with log.open("w") as handle:
        result = subprocess.run(command, cwd=case_dir, env=environment, stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"Case failed: {case_dir.name}; inspect {log}")
    assert source_fingerprint(project) == fingerprint, "Source changed during measurements; rerun with the final source"
    (case_dir / "fingerprint.txt").write_text(fingerprint + "\n")
    metric = read_csv(metrics_path)[0]
    metric["process_seconds"] = time.monotonic() - start
    write_csv(metrics_path, [metric])
    print("Completed", case_dir.name, f"in {metric['process_seconds']:.1f}s", flush=True)
    return metric


def config_values(name):
    n, p, source, qc, seed = name.split("_")
    return int(n[1:]), int(p[1:]), source, qc == "qc1", int(seed[4:])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, default=Path(__file__).parent / "reports/2026-11-18-bayesc-no-copy-revision")
    parser.add_argument("--sizes", default="128x257,512x1025,2048x4097")
    parser.add_argument("--seeds", default="73,123")
    parser.add_argument("--outer", type=int, default=200)
    parser.add_argument("--burnin", type=int, default=50)
    parser.add_argument("--block-size", type=int, default=64)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--julia", default="julia")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--compare-only", action="store_true")
    args = parser.parse_args()
    self_test()
    if args.self_test:
        return
    project = Path(__file__).resolve().parent.parent
    report = args.report.resolve()
    report.mkdir(parents=True, exist_ok=True)
    sizes = [tuple(map(int, size.split("x"))) for size in args.sizes.split(",")]
    seeds = list(map(int, args.seeds.split(",")))
    fingerprint = source_fingerprint(project)
    configs = []
    for n, p in sizes:
        for qc in (False, True):
            for seed in seeds:
                for source in ("adopted", "copied", "file"):
                    name = f"n{n}_p{p}_{source}_qc{int(qc)}_seed{seed}"
                    configs.append((report / "raw" / name, args, project, report, fingerprint))
    started = time.monotonic()
    if args.compare_only:
        assert all((config[0] / "fingerprint.txt").read_text().strip() == fingerprint for config in configs), "Source fingerprint changed"
        metrics = [read_csv(config[0] / "metrics.csv")[0] for config in configs]
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            metrics = list(pool.map(run_one, configs))
    expected_source = str((project / "src/JWAS.jl").resolve())
    assert all(metric["source_path"] == expected_source for metric in metrics), "Incorrect JWAS package source"
    for metric in metrics:
        validate_controls(metric, args)
    write_csv(report / "metrics.csv", metrics)
    comparisons, samples = [], []
    for n, p in sizes:
        for qc in (False, True):
            for seed in seeds:
                left_name = f"n{n}_p{p}_adopted_qc{int(qc)}_seed{seed}"
                for source in ("copied", "file"):
                    right_name = f"n{n}_p{p}_{source}_qc{int(qc)}_seed{seed}"
                    summary = dict(n=n, p=p, qc=qc, seed=seed, comparison=source)
                    for kind in ("animals", "markers"):
                        aligned, values = align(read_csv(report / "raw" / left_name / f"{kind}.csv"),
                                                read_csv(report / "raw" / right_name / f"{kind}.csv"), kind, n, p)
                        summary.update(values)
                        write_csv(report / "aligned" / f"{left_name}_vs_{source}_{kind}.csv", aligned)
                        for row in aligned[:2]:
                            samples.append({"case": left_name, "comparison": source, "kind": kind,
                                            "key": json.dumps(list(row.values())[:2]), "aligned_values": json.dumps(row, sort_keys=True)})
                    comparisons.append(summary)
    write_csv(report / "comparisons.csv", comparisons)
    write_csv(report / "joined_samples.csv", samples)
    print("Inspected joined-row sample:", json.dumps(samples[:4], indent=2), flush=True)
    metadata = dict(source_fingerprint_sha256=fingerprint, source_path=expected_source,
                    git_head=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=project, text=True).strip(),
                    platform=platform.platform(), python=sys.version.split()[0], sizes=sizes, seeds=seeds,
                    outer=args.outer, burnin=args.burnin, block_size=args.block_size, workers=args.workers,
                    case_count=len(metrics), comparison_count=len(comparisons), elapsed_seconds=time.monotonic() - started,
                    execution_mode="compare_only" if args.compare_only else "measure",
                    rss_units="bytes; Sys.maxrss() process-lifetime maximum, including startup/JIT/fixture generation",
                    all_comparisons_pass=all(value for row in comparisons for name, value in row.items() if name.startswith("match_")))
    (report / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))
    assert metadata["all_comparisons_pass"], "Keyed prediction comparison failed"


if __name__ == "__main__":
    main()
