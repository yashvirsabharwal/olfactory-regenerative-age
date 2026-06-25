#!/usr/bin/env python3
"""Build environment lockfiles and smoke-test reports."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.environment import (
    BASE_PACKAGES,
    LATENT_PACKAGES,
    lockfile_text,
    package_rows,
    pip_freeze,
    r_environment_yml,
    render_environment_markdown,
    runtime_rows,
    smoke_test_rows,
)
from ora.utils import ensure_parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/gateway.yaml")
    parser.add_argument("--report-out", default=None)
    parser.add_argument("--markdown-out", default=None)
    parser.add_argument("--smoke-out", default=None)
    parser.add_argument("--base-lock-out", default=None)
    parser.add_argument("--latent-lock-out", default=None)
    parser.add_argument("--full-freeze-out", default=None)
    parser.add_argument("--r-env-out", default=None)
    args = parser.parse_args()

    outputs = load_config(args.config).get("outputs", {})
    report_out = args.report_out or outputs.get("environment_report_tsv", "results/reports/environment_report.tsv")
    markdown_out = args.markdown_out or outputs.get("environment_report_md", "results/reports/environment_report.md")
    smoke_out = args.smoke_out or outputs.get("environment_smoke_tests_tsv", "results/reports/environment_smoke_tests.tsv")
    base_lock_out = args.base_lock_out or outputs.get("python_base_lock_txt", "requirements/python-base-lock.txt")
    latent_lock_out = args.latent_lock_out or outputs.get("python_latent_lock_txt", "requirements/python-latent-lock.txt")
    full_freeze_out = args.full_freeze_out or outputs.get("python_full_freeze_txt", "requirements/python-full-freeze.txt")
    r_env_out = args.r_env_out or outputs.get("r_bioconductor_environment_yml", "environments/r-bioconductor.yml")

    base_rows = package_rows(BASE_PACKAGES, group="python_base")
    latent_rows = package_rows(LATENT_PACKAGES, group="python_latent")
    package_table = pd.DataFrame([*base_rows, *latent_rows])
    runtime = pd.DataFrame(runtime_rows())
    smoke = pd.DataFrame(smoke_test_rows())

    package_table.to_csv(ensure_parent(report_out), sep="\t", index=False)
    smoke.to_csv(ensure_parent(smoke_out), sep="\t", index=False)
    ensure_parent(markdown_out).write_text(render_environment_markdown(runtime, package_table, smoke), encoding="utf-8")
    ensure_parent(base_lock_out).write_text(lockfile_text(base_rows, title="ORA Python base direct lock"), encoding="utf-8")
    ensure_parent(latent_lock_out).write_text(lockfile_text(latent_rows, title="ORA Python latent direct lock"), encoding="utf-8")
    ensure_parent(full_freeze_out).write_text(pip_freeze(), encoding="utf-8")
    ensure_parent(r_env_out).write_text(r_environment_yml(), encoding="utf-8")

    blocked = int(smoke["status"].eq("blocked").sum()) if "status" in smoke else 0
    failed = int(smoke["status"].eq("fail").sum()) if "status" in smoke else 0
    print(f"Wrote package report: {report_out} ({package_table.shape[0]} rows)")
    print(f"Wrote environment markdown: {markdown_out}")
    print(f"Wrote smoke tests: {smoke_out} ({failed} failed; {blocked} blocked)")
    print(f"Wrote Python base lock: {base_lock_out}")
    print(f"Wrote Python latent lock: {latent_lock_out}")
    print(f"Wrote Python full freeze: {full_freeze_out}")
    print(f"Wrote R/Bioconductor environment: {r_env_out}")


if __name__ == "__main__":
    main()
