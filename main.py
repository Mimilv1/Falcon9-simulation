import argparse
from datetime import datetime
from pathlib import Path

from Events import load_mission
from mission import plot_results, run_mission, terminal_plan


def parse_args():
    parser = argparse.ArgumentParser(
        description="CRS-10 ascent, shutdown, separation and booster coast."
    )
    parser.add_argument("--case", choices=["ascent", "terminal"], default="ascent")
    parser.add_argument("--until", type=float, default=None)
    parser.add_argument("--dt", type=float, default=None)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--no-navigation", action="store_true")
    parser.add_argument("--no-show", action="store_true")
    parser.add_argument("--save", type=Path, default=None)
    parser.add_argument("--plot-only", type=Path, default=None)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.plot_only:
        plot_results(args.plot_only, show=not args.no_show, save=args.save)
        return
    plan = terminal_plan() if args.case == "terminal" else load_mission()
    until = (
        args.until
        if args.until is not None
        else plan["terminal_test"]["duration_s"] if args.case == "terminal" else 160.0
    )
    output = args.output or Path(__file__).parent / "results" / datetime.now().strftime(
        "crs10-" + args.case + "-%Y%m%d-%H%M%S"
    )
    print(f"Running CRS-10 to {until:g} s. Results: {output}", flush=True)
    result = run_mission(
        output,
        until_s=until,
        step_s=args.dt,
        seed=args.seed,
        navigation=not args.no_navigation,
        plan=plan,
    )
    final = result["final_truth"]
    print(
        f"{result['status']}: t={result['actual_end_s']:.3f} s, phase={result['phase']}"
    )
    print(
        f"Altitude={final['cg_altitude_m']/1000:.2f} km, speed={final['speed_mps']:.1f} m/s, mass={final['mass_kg']:.1f} kg"
    )
    print(
        f"Stage-1 propellant={final['stage1_propellant_kg']:.1f} kg; next event={result['next_event']}"
    )
    if result["detail"]:
        print(result["detail"])
    if "navigation_position_rms_m" in result:
        print(
            f"Navigation position RMS={result['navigation_position_rms_m']:.3f} m (truth attitude)"
        )
    print("Full return flight completed:", result["full_mission_completed"])
    plot_results(output, show=not args.no_show, save=args.save)
    if result["status"] != "requested_interval_complete":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
