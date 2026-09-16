"""Validate or execute the separately frozen MacroMap fixed-program analysis.

Default mode only validates a parent-pinned execution plan. --execute is also
required to read real expression values. No model/source/endpoint search exists.
Outputs always use a fresh ignored directory named in the frozen plan.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

try:
    from . import analyze_macromap_program_context as method
except ImportError:
    import analyze_macromap_program_context as method


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True, type=Path)
    parser.add_argument("--expected-plan-sha256", required=True)
    parser.add_argument("--root", type=Path, default=method.ROOT)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args(argv)
    if method.sha256(args.plan) != args.expected_plan_sha256:
        raise ValueError("Execution plan differs from supplied parent freeze hash")
    if args.execute:
        result = method.run_frozen_analysis(args.plan, expected_plan_sha256=args.expected_plan_sha256, root=args.root)
        print(json.dumps({"stage": result["stage"], "plan_sha256": result["plan_sha256"]}, indent=2))
    else:
        plan = json.loads(args.plan.read_text())
        design = method.validate_execution_plan(plan, root=args.root)
        print(json.dumps({"stage": "execution_plan_validated_no_expression_read",
                          "prepared_summary_sha256": plan["prepared_summary_sha256"],
                          "samples": len(design["samples"]), "real_expression_read": False}, indent=2))


if __name__ == "__main__":
    main()
