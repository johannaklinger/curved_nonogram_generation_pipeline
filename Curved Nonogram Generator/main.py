import math
import argparse

import Nonogram_settings as NS
from Nonogram_generation_functions import (calculate_ECEs_from_path, calculate_connection_curves, generate_output,
                                           calculate_bridge_curves)


def do_generation(arguments):
    filled = arguments.fill

    vb_width, vb_height = NS.input_preprocess(arguments.input_file)

    min_dist = (
        arguments.min_dist
        if arguments.min_dist is not None
        else min(vb_width, vb_height) / 20
    )

    min_gap = (
        arguments.min_gap
        if arguments.min_gap is not None
        else min_dist * 0.2
    )

    max_dist = (
        arguments.max_dist
        if arguments.max_dist is not None
        else (vb_width + vb_height) / 5
    )

    NS.set_switches(
        arguments.all_bridge_curves,
        arguments.inner_bridge_curves,
        arguments.bridge_curves_only,
        arguments.find_best_bridge,
        arguments.do_checks,
        arguments.early_checks,
        arguments.check_angles if arguments.w_angle > 0 else False,
        arguments.check_intersections if arguments.w_angle > 0 else False,
        arguments.check_distance if arguments.w_angle > 0 else False,
        arguments.additional_curves
    )

    NS.set_variables(
        arguments.max_loops,
        arguments.min_angle,
        min_dist,
        min_gap,
        arguments.pen_threshold,
        arguments.w_angle,
        arguments.w_intersection,
        arguments.w_distance,
        arguments.alpha,
        arguments.beta,
        max_dist,
        arguments.control_point_dist
    )

    calculate_ECEs_from_path()
    calculate_bridge_curves()
    calculate_connection_curves()
    generate_output(arguments.out_file, filled)


def parse_args():
    parser = argparse.ArgumentParser(description="Generate curved nonogram")

    parser.add_argument(
        "--input-file",
        default="inputs/muffin.svg"
    )

    parser.add_argument(
        "--out-file",
        default="outputs/output.svg"
    )

    parser.add_argument(
        "--fill",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="additionally generate output with filled solution image"
    )

    parser.add_argument(
        "--all-bridge-curves",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="consider all ECEs for bridge curve generation"
    )

    parser.add_argument(
        "--inner-bridge-curves",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="consider only ECEs pointing inwards for bridge curve generation"
    )

    parser.add_argument(
        "--bridge-curves-only",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="generate only bridge curves, no connection curves"
    )

    parser.add_argument(
        "--find-best-bridge",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="find bridges with the lowest penalty"
    )

    parser.add_argument(
        "--do-checks",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="perform optimization"
    )

    parser.add_argument(
        "--early-checks",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="perform optimization immediately after curve generation"
    )

    parser.add_argument(
        "--check-angles",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="consider intersection angles during optimization"
    )

    parser.add_argument(
        "--check-intersections",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="consider distance between intersections during optimization"
    )

    parser.add_argument(
        "--check-distance",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="consider distance between curves during optimization"
    )

    parser.add_argument(
        "--additional-curves",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="generate additional curves"
    )
    parser.add_argument("--max-loops", type=int, default=30,
                        help="maximum amount of iterations per curve during optimization")
    parser.add_argument("--min-angle", type=float, default=20,
                        help="minimum allowed intersection angle")
    parser.add_argument("--min-dist", type=float,
                        help="minimum allowed distance between intersections")
    parser.add_argument("--min-gap", type=float,
                        help="minimum allowed gap between curves")
    parser.add_argument("--pen-threshold", type=float, default=1,
                        help="minimum allowed penalty")
    parser.add_argument("--w-angle", type=float, default=1,
                        help="weight of angle_penalty")
    parser.add_argument("--w-intersection", type=float, default=3,
                        help="weight of intersection_penalty")
    parser.add_argument("--w-distance", type=float, default=1,
                        help="weight of distance_penalty")
    parser.add_argument("--alpha", type=float, default=math.pi / 3,
                        help="outer angle in radians to limit end points on border")
    parser.add_argument("--beta", type=float, default=math.pi / 20,
                        help="inner angle in radians to limit end point on border")
    parser.add_argument("--max-dist", type=float,
                        help="maximum allowed gap between end points on border for additional curve generation")
    parser.add_argument("--control-point-dist", type=float, default=0.3,
                        help="distance between control points on additional curves")
    parser.add_argument(
        "--parameter-help",
        action="store_true",
        help="show explanation of parameters"
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    do_generation(args)
