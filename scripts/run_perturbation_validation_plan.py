#!/usr/bin/env python3
"""Build perturbation/organoid/ALI validation planning artifacts."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ora.config import load_config
from ora.perturbation_validation import (
    build_minimum_experiment_table,
    build_perturbation_candidate_matrix,
    build_perturbation_search_log,
    render_perturbation_validation_plan,
    write_perturbation_validation_outputs,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/perturbation_validation_candidates.yaml")
    parser.add_argument("--candidate-out", default="results/tables/perturbation_validation_candidate_matrix.tsv")
    parser.add_argument("--search-log-out", default="results/tables/perturbation_validation_search_log.tsv")
    parser.add_argument("--experiment-out", default="results/tables/perturbation_minimum_experiment_design.tsv")
    parser.add_argument("--plan-out", default="docs/perturbation_validation_plan.md")
    args = parser.parse_args()

    config = load_config(args.config)
    candidates = build_perturbation_candidate_matrix(config)
    search_log = build_perturbation_search_log(config)
    minimum_experiment = build_minimum_experiment_table(config)
    plan_md = render_perturbation_validation_plan(
        config=config,
        candidates=candidates,
        search_log=search_log,
        minimum_experiment=minimum_experiment,
    )
    write_perturbation_validation_outputs(
        candidates=candidates,
        search_log=search_log,
        minimum_experiment=minimum_experiment,
        plan_md=plan_md,
        candidate_out=args.candidate_out,
        search_log_out=args.search_log_out,
        experiment_out=args.experiment_out,
        plan_out=args.plan_out,
    )
    high_priority = int(candidates["priority"].le(1).sum()) if "priority" in candidates else 0
    print(
        "Wrote perturbation validation plan: "
        f"{args.plan_out} ({candidates.shape[0]} candidates; {high_priority} high-priority)"
    )


if __name__ == "__main__":
    main()
