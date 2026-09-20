import math
import random
from bisect import bisect_left

import numpy as np
import svgpathtools
from svgpathtools import CubicBezier, Path, wsvg, Line

import Nonogram_settings as NS


# remove tuple from list of tuples
def remove_tuple(list_of_tuples, key_tuple):
    x, y = key_tuple
    for s in list_of_tuples:
        a, b = s
        if a == x:
            list_of_tuples.remove(s)


# remove quadruple from list of quadruples
def remove_quadruple(list_of_triples, key_triple):
    a, b, c = key_triple
    for s in list_of_triples:
        d, e, f, _ = s
        if a == d and b == e and c == f:
            list_of_triples.remove(s)


# find index of intersection in list of intersections
def find_index(s1, first_c):
    for i in range(len(first_c.intersections)):
        t, _, _, _ = first_c.intersections[i]
        if s1 == t:
            return i
    raise Exception(f"no intersection at {s1} found on {first_c.curve} with intersections \n{first_c.intersections}")


# transform non-cubic-Bézier-curves into list of cubic Bézier curve
def segment_to_cubics(seg):
    if isinstance(seg, Line) and (seg.start == seg.end):
        return []
    if isinstance(seg, CubicBezier):
        if abs(seg.start - seg.control1) < 1e-5:
            if abs(seg.end - seg.control2) < 1e-5:
                v = seg.end - seg.start
                return [CubicBezier(seg.start, seg.start + v / 3.0, seg.start + 2 * v / 3.0, seg.end)]
            else:
                return [CubicBezier(seg.start, (seg.start + seg.control2) * 0.5, seg.control2, seg.end)]
        elif abs(seg.end - seg.control2) < 1e-5:
            return [CubicBezier(seg.start, seg.control1, (seg.control1 + seg.end) * 0.5, seg.end)]
        return [seg]
    elif hasattr(seg, 'as_cubic_curves'):
        return list(seg.as_cubic_curves())
    else:
        # fallback: Line
        v = seg.end - seg.start
        return [CubicBezier(seg.start, seg.start + v / 3.0, seg.start + 2 * v / 3.0, seg.end)]


# returns True if point pt lies within input figure
def is_enclosed(pt):
    opt = complex(NS.min_x - 10, NS.min_y - 10)
    sum_intersections = 0
    for p in NS.paths:
        sum_intersections += len(Path(Line(pt, opt)).intersect(p))
    if sum_intersections % 2:
        return True
    else:
        return False


# returns True if ece1 and ece2 connect G^1-continuously
def is_g1_continuous_with_ece(ece1, ece2, tol=0.22):
    d1 = ece1.start - ece1.control1
    d2 = ece2.control1 - ece2.start

    # Normalize directions, would otherwise test C^1-continuity
    d1_unit = d1 / abs(d1)
    d2_unit = d2 / abs(d2)

    return abs(d1_unit - d2_unit) < tol


# add new ECEs to list of ECEs if curves prev and curr do not connect to other ECEs G^1-continuously
def add_ECE_if_ness(prev, curr):
    connected_prev = False
    connected_curr = False
    if is_on_border(prev.end):
        connected_prev = True
    if is_on_border(curr.start):
        connected_curr = True
    if connected_prev and connected_curr:
        return

    p2 = prev.control2
    p1 = curr.control1
    if abs(p2 - prev.end) < 1e-5:
        p2 = prev.point(0.9)
    if abs(p1 - curr.control1) < 1e-5:
        p1 = curr.point(0.1)
    ece1 = NS.ECE(prev.end, p2)
    ece2 = NS.ECE(curr.start, p1)

    if abs(prev.end - prev.start) < 1e-5:
        if abs(curr.start - curr.end) < 1e-5:
            return
        connected_prev = True
    elif abs(curr.start - curr.end) < 1e-5:
        connected_curr = True

    # if not continuous at all, search other discontinuous ECEs for possible continuity
    if abs(prev.end - curr.start) > 1e-5 or connected_prev or connected_curr:
        copied_ECEs = NS.open_ECEs.copy()
        for oE in copied_ECEs:
            if not connected_prev and abs(oE.start - prev.end) < 1e-5:
                if not is_g1_continuous_with_ece(oE, ece1):
                    append_ECE(p2, prev.end)
                    append_ECE(oE.control1, oE.start)
                connected_prev = True
                NS.open_ECEs.remove(oE)
                continue
            if not connected_curr and abs(oE.start - curr.start) < 1e-5:
                if not is_g1_continuous_with_ece(oE, ece2):
                    append_ECE(p1, curr.start)
                    append_ECE(oE.control1, oE.start)
                connected_curr = True
                NS.open_ECEs.remove(oE)
                continue
        if not connected_prev:
            NS.open_ECEs.append(ece1)
        if not connected_curr:
            NS.open_ECEs.append(ece2)
        return

    if connected_prev or connected_curr:
        return

    if not is_g1_continuous_with_ece(ece1, ece2):
        append_ECE(p2, prev.end)
        append_ECE(p1, curr.start)


# append ECE to list of ECEs with adjusted distance between start point and second control point
def append_ECE(c1, start):
    v1 = c1 - start
    dist1 = abs(v1)
    v1_norm = v1 / dist1
    dist1 = max(NS.min_dist, dist1)
    c1 = start + v1_norm * dist1

    new_ece = NS.ECE(start, c1, False)
    if any(ece == new_ece for ece in NS.ECEs):
        return

    if NS.inner_bridge_curves:
        check_point = start - 5 * v1_norm
        if is_enclosed(check_point):
            NS.inner_ECEs.append(new_ece)
    NS.ECEs.append(new_ece)


# restrict p1 to border if p1 lies outside the border
def clamp_p1_to_border(p1):
    r = max(min(p1.real, NS.max_x), NS.min_x)
    i = max(min(p1.imag, NS.max_y), NS.min_y)
    return complex(r, i)


# returns True if point p lies on the border
def is_on_border(p) -> bool:
    r = p.real
    i = p.imag
    tol = 1e-8
    if abs(r - NS.min_x) < tol or abs(r - NS.max_x) < tol:
        return True
    if abs(i - NS.min_y) < tol or abs(i - NS.max_y) < tol:
        return True
    return False


# randomly set p2 bounded by p1 and p3
def set_p2_between_p1_and_p3(p1, p3):
    if is_on_border(p1):
        return p1

    lower = int(p1.real)
    upper = int(p3.real)
    if lower < upper:
        p2_real = random.randrange(lower, upper)
    elif lower > upper:
        p2_real = random.randrange(upper, lower)
    else:
        p2_real = lower

    lower = int(p1.imag)
    upper = int(p3.imag)
    if lower < upper:
        p2_imag = random.randrange(lower, upper)
    elif lower > upper:
        p2_imag = random.randrange(upper, lower)
    else:
        p2_imag = lower
    return complex(p2_real, p2_imag)


# set p3 randomly within range limited by alpha and beta, while keeping enough space between other end points,
# then set p2
def p2_and_p3_angled(start, direction, p1):
    if is_on_border(p1):
        d = p3_2dist(p1)
        index = bisect_left([dist for (cords, dist) in NS.end_points], d)
        NS.end_points.insert(index, (p1, d))
        return p1, p1
    p1 = start + direction
    dir_norm = direction / abs(direction)

    sign = 2 * random.getrandbits(1) - 1  # pos: counterclockwise rotation, neg: clockwise rotation
    theta = NS.alpha * sign
    dir1 = dir_norm * np.exp(1j * theta)
    dir1 = dir1 / abs(dir1)

    theta = NS.beta * sign
    dir2 = dir_norm * np.exp(1j * theta)
    dir2 = dir2 / abs(dir2)

    if sign == -1:
        temp = dir1
        dir1 = dir2
        dir2 = temp

    range1 = Line(start, start + dir1 * (NS.vb_height + NS.vb_width))
    range2 = Line(start, start + dir2 * (NS.vb_height + NS.vb_width))

    first = 0
    second = 0

    # find range limited by alpha and beta
    for c in [NS.border_left.curve, NS.border_top.curve, NS.border_bottom.curve, NS.border_right.curve]:
        inters = svgpathtools.bezier_by_line_intersections(c, range1)
        if len(inters) > 0:
            _, s = inters[-1]
            p = range1.point(s)
            first = p3_2dist(p)
            break

    for c in [NS.border_left.curve, NS.border_top.curve, NS.border_bottom.curve, NS.border_right.curve]:
        inters = svgpathtools.bezier_by_line_intersections(c, range2)
        if len(inters) > 0:
            _, s = inters[0]
            p = range2.point(s)
            second = p3_2dist(p)
            break

    p3_extended = False
    if first < second:
        first += 2 * NS.vb_width + 2 * NS.vb_height
        p3_extended = True

    # respect all end point that already lie within the range
    obstacles = [second]
    for i in range(len(NS.end_points)):
        point, distance = NS.end_points[i]
        if p3_extended:
            if second < distance:
                obstacles.append(distance)
                continue
            elif distance < first - 2 * NS.vb_width - 2 * NS.vb_height:
                obstacles.append(distance + 2 * NS.vb_width + 2 * NS.vb_height)
                continue
            if distance > first - 2 * NS.vb_width - 2 * NS.vb_height:
                break
        else:
            if second < distance < first:
                obstacles.append(distance)
                continue
            if distance > first:
                break
    obstacles.append(first)

    ranges = []
    for i in range(len(obstacles) - 1):
        left = obstacles[i]
        right = obstacles[i + 1]
        if not i == 0:
            left += NS.min_dist
        if not i == len(obstacles) - 2:
            right -= NS.min_dist
        if right < 0:
            right += 2 * NS.vb_width + 2 * NS.vb_height
        if int(right) <= int(left):
            continue
        pos_p = random.randrange(int(left), int(right))
        ranges.append(pos_p)

    # assign value to p3 randomly from calculated range
    if len(ranges) == 0:
        p3_dist = (obstacles[0] + obstacles[1]) / 2
    else:
        p3_dist = ranges[random.randint(0, len(ranges) - 1)]

    p3_dist = p3_dist % (NS.vb_width * 2 + NS.vb_height * 2)
    new_p3 = dist2p3(p3_dist)
    index = bisect_left([dist for (cords, dist) in NS.end_points], p3_dist)
    NS.end_points.insert(index, (new_p3, p3_dist))

    # calculate new value for p2
    p2 = set_p2_between_p1_and_p3(p1, new_p3)
    return p2, new_p3


# calculate distance to p3 along the border from coordinates
def p3_2dist(p):
    tol = 1e-4
    r = p.real
    i = p.imag

    if abs(i - NS.min_y) < tol:
        return r - NS.min_x
    if abs(r - NS.max_x) < tol:
        return NS.vb_width + i - NS.min_y
    if abs(i - NS.max_y) < tol:
        return NS.vb_width + NS.vb_height + (NS.max_x - r)
    if abs(r - NS.min_x) < tol:
        return 2 * NS.vb_width + NS.vb_height + (NS.max_y - i)

    raise Exception(f"point {p} not on border")


# calculate coordinates of p3 from distance along the border
def dist2p3(dist):
    while True:
        if dist < NS.vb_width:
            return complex(dist + NS.min_x, NS.min_y)
        dist = dist - NS.vb_width
        if dist < NS.vb_height:
            return complex(NS.max_x, dist + NS.min_y)
        dist = dist - NS.vb_height
        if dist < NS.vb_width:
            return complex(NS.max_x - dist, NS.max_y)
        dist = dist - NS.vb_width
        if dist < NS.vb_height:
            return complex(NS.min_x, NS.max_y - dist)
        dist = dist - NS.vb_height


# calculate distance between two points on the border (given as their distance along the border)
def gap_between(d1, d2):
    if d2 < d1:
        return d2 + (2 * NS.vb_width + 2 * NS.vb_height - d1)
    else:
        return abs(d1 - d2)


# returns the angel between two vectors
def angle_between(v1, v2):
    v1_unit = v1 / abs(v1)
    v2_unit = v2 / abs(v2)

    dot = (v1_unit.real * v2_unit.real + v1_unit.imag * v2_unit.imag)
    dot = max(-1.0, min(1.0, dot))  # numerical safety

    return math.degrees(math.acos(dot))


# calculate angle penalty
def angle_pen(angle, threshold):
    if angle > 90:
        raise Exception("angle > 90")
    if angle < threshold:
        return threshold - angle
    else:
        return 0


# adjust curve during optimization by assigning new parameters
def adjustment_process(nc):
    curve = nc.curve
    delete_intersections(nc)
    if nc.is_additional:
        curve.control2 = set_p2_between_p1_and_p3(nc.curve.control1, nc.curve.end)
    else:
        old_p3 = nc.curve.end
        old_distance = p3_2dist(old_p3)
        remove_tuple(NS.end_points, (old_p3, old_distance))
        curve.control2, curve.end = p2_and_p3_angled(curve.start, curve.unit_tangent(0), curve.control1)
    update_all_intersections(nc)


# calculate intersection penalty
def intersection_pen(curve, t1, t2, threshold):
    p1 = curve.point(t1)
    p2 = curve.point(t2)
    if threshold > abs(p1 - p2):
        if abs(p1 - p2) < threshold / 2:
            return (threshold - abs(p1 - p2)) * 2
        return threshold - abs(p1 - p2)
    return 0


# calculate distance penalty
def distance_pen(cc) -> float:
    penalty = 0
    number_steps = 7
    for bezier in NS.all_bezier_curves:
        if bezier.__eq__(cc):
            continue
        if any(bezier.__eq__(c) for _, _, c, _ in cc.intersections):
            continue
        if cc.curve.start == bezier.curve.start or cc.curve.start == bezier.curve.end:
            continue
        for i in range(number_steps):
            t = i / 10 + 0.3
            point = cc.curve.point(t)
            seg = bezier.curve
            (d_min_np, _), _ = seg.radialrange(point)
            d_min = float(d_min_np)
            if d_min < NS.min_gap:
                if d_min < NS.min_gap * 0.5:
                    penalty += (NS.min_gap - d_min) * 2
                penalty += NS.min_gap - d_min
    return penalty


# calculate intersection penalty of intersections on intersected curve
def neighbors_on_intersecting_curve_too_close(intersection):
    _, s1, first_c, _ = intersection
    index_on_first = find_index(s1, first_c)
    penalty = 0
    if index_on_first > 0:
        i2 = first_c.intersections[index_on_first - 1]
        t2, _, _, _ = i2
        penalty += intersection_pen(first_c.curve, s1, t2, NS.min_dist)
    else:
        penalty += intersection_pen(first_c.curve, 0, s1, NS.min_dist)

    if index_on_first < len(first_c.intersections) - 1:
        i2 = first_c.intersections[index_on_first + 1]
        t2, _, _, _ = i2
        penalty += intersection_pen(first_c.curve, s1, t2, NS.min_dist)
    else:
        penalty += intersection_pen(first_c.curve, s1, 1, NS.min_dist)
    return penalty


# calculate penalty for curve nc and optimize if necessary
def calculate_penalty_and_optimize(nc, distance_threshold, angle_threshold, test_only=False):
    if nc.is_fixed and not NS.find_best_bridge:
        return 0, 0, 0

    new_curve = nc.curve
    cnt = 0
    best_solution = NS.CC(CubicBezier(nc.curve.start, nc.curve.control1, nc.curve.control2, nc.curve.end),
                          [], nc.is_fixed, nc.is_additional, nc.mate, penalties=nc.penalties)
    if best_solution.penalties[0] + best_solution.penalties[1] + best_solution.penalties[2] < NS.pen_threshold:
        return best_solution.penalties[0], best_solution.penalties[1], best_solution.penalties[2]
    angle_p, intersection_p, distance_p = 0, 0, 0
    while cnt < NS.max_loops:
        angle_p, intersection_p, distance_p = 0, 0, 0
        intersection_list = nc.intersections

        if len(intersection_list) > 0:
            i1 = intersection_list[0]
            t1, s1, first_c, angle = i1
            if t1 > 0.1:
                if NS.check_angles:
                    angle_p = angle_pen(angle, angle_threshold) * NS.w_angle
                if NS.check_intersections:
                    t2, _, _, _ = intersection_list[len(intersection_list) - 1]
                    intersection_p = neighbors_on_intersecting_curve_too_close(i1) * NS.w_intersection
                    intersection_p += intersection_pen(new_curve, 0, t1, distance_threshold) * NS.w_intersection
                    intersection_p += intersection_pen(new_curve, t2, 1, distance_threshold) * NS.w_intersection

            for i in range(len(intersection_list) - 1):
                i2 = intersection_list[i + 1]
                t2, s2, second_c, angle = i2
                if t2 > 0.1:
                    if NS.check_angles:
                        angle_p += angle_pen(angle, angle_threshold) * NS.w_angle
                    if NS.check_intersections:
                        i1 = intersection_list[i]
                        t1, s1, first_c, _ = i1
                        intersection_p += neighbors_on_intersecting_curve_too_close(i1) * NS.w_intersection
                        intersection_p += intersection_pen(new_curve, t1, t2,
                                                           distance_threshold) * NS.w_intersection

        if NS.check_distance:
            distance_p = distance_pen(nc) * NS.w_distance
        pen = intersection_p + angle_p + distance_p

        if pen < best_solution.penalties[0] + best_solution.penalties[1] + best_solution.penalties[2]:
            best_solution = NS.CC(CubicBezier(nc.curve.start, nc.curve.control1, nc.curve.control2, nc.curve.end),
                                  [], nc.is_fixed, nc.is_additional, nc.mate,
                                  angle_p, intersection_p, distance_p)

        cnt += 1

        if pen < NS.pen_threshold:
            nc.penalties = [angle_p, intersection_p, distance_p]
            return angle_p, intersection_p, distance_p

        if test_only:
            return angle_p, intersection_p, distance_p
        elif cnt == NS.max_loops:
            d = p3_2dist(nc.curve.end)
            remove_tuple(NS.end_points, (nc.curve.end, d))

            delete_intersections(nc)
            nc.curve = best_solution.curve
            nc.intersections = best_solution.intersections
            nc.penalties = best_solution.penalties

            d = p3_2dist(nc.curve.end)
            index = bisect_left([dist for (cords, dist) in NS.end_points], d)
            NS.end_points.insert(index, (nc.curve.end, d))
            update_all_intersections(nc)
            return nc.penalties[0], nc.penalties[1], nc.penalties[2]
        else:
            adjustment_process(nc)

    nc.penalties = [angle_p, intersection_p, distance_p]
    return angle_p, intersection_p, distance_p


# search best positions to connect sub-curves of a connection curves and the corresponding control points
def search_best_p1(nc1, nc2):
    cnt = 0
    copy_c1 = NS.CC(CubicBezier(nc1.curve.start, nc1.curve.control1, nc1.curve.control2, nc1.curve.end),
                    [], nc1.is_fixed, nc1.is_additional, nc1.mate, penalties=nc1.penalties)
    copy_c2 = NS.CC(CubicBezier(nc2.curve.start, nc2.curve.control1, nc2.curve.control2, nc2.curve.end),
                    [], nc2.is_fixed, nc2.is_additional, nc2.mate, penalties=nc2.penalties)
    best_pen = float('inf')
    pen1 = nc1.penalties[0] + nc1.penalties[1] + nc1.penalties[2]
    pen2 = nc2.penalties[0] + nc2.penalties[1] + nc2.penalties[2]

    while cnt < NS.max_loops:
        pen = pen1 + pen2
        if pen < NS.pen_threshold * 2:
            return

        if pen < best_pen:
            copy_c1 = NS.CC(CubicBezier(nc1.curve.start, nc1.curve.control1, nc1.curve.control2, nc1.curve.end),
                            [], nc1.is_fixed, nc1.is_additional, nc1.mate, penalties=nc1.penalties)
            copy_c2 = NS.CC(CubicBezier(nc2.curve.start, nc2.curve.control1, nc2.curve.control2, nc2.curve.end),
                            [], nc2.is_fixed, nc2.is_additional, nc2.mate, penalties=nc2.penalties)
            best_pen = pen

        cnt += 1
        if cnt == NS.max_loops:
            delete_intersections(nc1)
            delete_intersections(nc2)

            nc1.curve = copy_c1.curve
            nc1.penalties = copy_c1.penalties
            nc2.curve = copy_c2.curve
            nc2.penalties = copy_c2.penalties

            update_all_intersections(nc1)
            update_all_intersections(nc2)
            return

        line = Line(nc1.curve.end, nc2.curve.end)
        line = line.cropped(0.5 - NS.control_point_dist * 0.5, 0.5 + NS.control_point_dist * 0.5)
        middle = line.point(0.5)
        movement = complex(random.randrange(int(NS.min_x - middle.real + NS.min_dist * 2),
                                            int(NS.max_y - middle.real - NS.min_dist * 2)),
                           random.randrange(int(NS.min_y - middle.imag + NS.min_dist * 2),
                                            int(NS.vb_height + NS.min_y - middle.imag - NS.min_dist * 2)))
        line = line.translated(movement)
        start = line.point(0.5)

        delete_intersections(nc1)
        delete_intersections(nc2)
        nc1.curve.start = start
        nc2.curve.start = start
        nc1.curve.control1 = clamp_p1_to_border(line.start)
        nc2.curve.control1 = clamp_p1_to_border(line.end)
        update_all_intersections(nc1)
        update_all_intersections(nc2)

        pen1_a, pen1_i, pen1_d = calculate_penalty_and_optimize(nc1, NS.min_dist, NS.min_angle, test_only=True)
        pen2_a, pen2_i, pen2_d = calculate_penalty_and_optimize(nc2, NS.min_dist, NS.min_angle, test_only=True)
        pen1 = pen1_a + pen1_i + pen1_d
        pen2 = pen2_a + pen2_i + pen2_d


# delete all intersections of curve cc from curves intersected by cc
def delete_intersections(cc):
    for t, s, other_curve, _ in cc.intersections:
        remove_quadruple(other_curve.intersections, (s, t, cc))
    cc.intersections = []


# calculate all intersections of curve cc and insert to list of intersections for cc and intersected curves
def update_all_intersections(cc, tol=1e-2):
    cc.intersections = []
    for bezier in NS.all_bezier_curves:
        if bezier.curve.end == cc.curve.end:
            continue
        try:
            intersections = svgpathtools.bezier_intersections(cc.curve, bezier.curve, 40, tol=1e-1, tol_deC=1e-5)
        except Exception:
            print(f"intersections timed out for {cc.curve} with\n{cc.intersections} and\n{bezier.curve} "
                  f"with\n{bezier.intersections}")
            return
        for t, s in intersections:
            if t < tol or s < tol:
                continue
            index = bisect_left([para for (para, _, intersecting_curve, _) in cc.intersections], t)
            d1 = cc.curve.unit_tangent(t)
            d2 = bezier.curve.unit_tangent(s)
            angle = angle_between(d1, d2)
            if angle > 90:
                angle = 180 - angle
            cc.intersections.insert(index, (t, s, bezier, angle))

            index = bisect_left([para for (para, _, intersecting_curve, _) in bezier.intersections], s)
            bezier.intersections.insert(index, (s, t, cc, angle))


# test for all ECEs if they fulfill the criteria for creating a bridge curves with ece
def find_bridge_partner(ece, open_inner_eces):
    if ece.is_connected():
        raise Exception(f"tried to connect ECE twice:\n\t{ece}")
    dir1 = ece.start - ece.control1
    possible_p1 = ece.start + dir1

    possible_mates = []
    for oie in open_inner_eces:
        if ece.start == oie.start:
            continue
        dir_start2start = oie.start - ece.start
        dir2 = oie.start - oie.control1
        possible_p2 = oie.start + dir2
        angle_to_other_start = angle_between(dir1, dir_start2start)
        if (-100 < angle_to_other_start < 100
                and abs(dir1) < abs(ece.start - possible_p2) - NS.min_dist / 2
                and abs(dir2) < abs(oie.start - possible_p1) - NS.min_dist / 2
                and abs(dir_start2start) > abs(ece.start - possible_p2)
                and abs(dir_start2start) > abs(oie.start - possible_p1)
                and abs(dir_start2start) > NS.min_dist):
            possible_mates.append(oie)

    best_mate_ece = find_best_ece(ece, possible_mates)
    if best_mate_ece:
        open_inner_eces.remove(best_mate_ece)
        open_inner_eces.remove(ece)
        return 2
    return 0


# find ECE among all possible ECEs that creates the lowest penalty when generating a bridge curve to ece
def find_best_ece(ece, possible_mates) -> NS.ECE | None:
    dir1 = ece.start - ece.control1
    length = len(possible_mates)
    if length == 0:
        return None
    best_bridge = NS.CC(None, [], True)
    best_mate_ece = None
    while length > 0:
        e = possible_mates[random.randrange(0, length)]
        dir2 = e.start - e.control1
        p1 = clamp_p1_to_border(ece.start + dir1)
        p2 = clamp_p1_to_border(e.start + dir2)
        new_cc = NS.CC(CubicBezier(ece.start, p1, p2, e.start), [], True)

        if NS.find_best_bridge or NS.do_checks:
            update_all_intersections(new_cc)
        if NS.find_best_bridge:
            angle_p, intersection_p, distance_p = calculate_penalty_and_optimize(
                new_cc, NS.min_dist, NS.min_angle, test_only=True)

            penalty = angle_p + intersection_p + distance_p
            if penalty < best_bridge.penalties[0] + best_bridge.penalties[1] + best_bridge.penalties[2]:
                best_bridge = NS.CC(
                    CubicBezier(new_cc.curve.start, new_cc.curve.control1, new_cc.curve.control2, new_cc.curve.end),
                    [], new_cc.is_fixed, new_cc.is_additional, new_cc.mate,
                    angle_p, intersection_p, distance_p)
                best_mate_ece = e

            if penalty < NS.pen_threshold:
                break
            else:
                length -= 1
                if length == 0:
                    if (best_bridge.penalties[0] + best_bridge.penalties[1] + best_bridge.penalties[2]
                            < NS.pen_threshold):
                        delete_intersections(new_cc)
                        update_all_intersections(best_bridge)
                        break
                    delete_intersections(new_cc)
                    return None
                possible_mates.remove(e)
                delete_intersections(new_cc)
        else:
            best_bridge = NS.CC(
                CubicBezier(new_cc.curve.start, new_cc.curve.control1, new_cc.curve.control2, new_cc.curve.end),
                new_cc.intersections, new_cc.is_fixed, new_cc.is_additional, new_cc.mate)
            best_mate_ece = e
            break
    NS.all_bezier_curves.append(best_bridge)
    NS.all_new_curves.append(best_bridge)

    ece.connected = True
    best_mate_ece.connected = True
    return best_mate_ece


# generate additional curves from every part of the border with large gaps between two end points
def generate_additional_curves():

    number_end_points = len(NS.end_points)
    free_p3s = []

    if number_end_points == 0:
        border_length = 2 * NS.vb_width + 2 * NS.vb_height
        dist_add = border_length / 12
        for i in range(0, 12):
            new_p = dist2p3((0.5 + i) * dist_add)
            free_p3s.append(new_p)
    elif number_end_points == 1:
        raise Exception("there can not be only one point on the border")
    else:
        for i in range(len(NS.end_points)):
            p1, d1 = NS.end_points[i]
            p2, d2 = NS.end_points[(i + 1) % len(NS.end_points)]
            gap = gap_between(d1, d2)
            if gap < NS.max_dist:
                continue
            num_inserts = math.ceil(gap / NS.max_dist)
            new_gap = gap / num_inserts
            for _ in range(num_inserts - 1):
                free_p3s.append(dist2p3(d1 + new_gap))
                d1 += new_gap

    while True:
        if len(free_p3s) < 2:
            break
        rand1 = random.randrange(0, len(free_p3s))
        b1 = free_p3s[rand1]
        free_p3s.remove(b1)
        rand2 = random.randrange(0, len(free_p3s))
        b2 = free_p3s[rand2]
        free_p3s.remove(b2)
        add_additional_curve(b1, b2)


# generate additional curve which starts end ends at border1, border 2
def add_additional_curve(border1, border2):
    line = Line(border1, border2)
    line = line.cropped(0.5 - NS.control_point_dist * 0.5, 0.5 + NS.control_point_dist * 0.5)
    middle = line.point(0.5)
    movement = complex(random.randrange(int(NS.min_x - middle.real + NS.min_dist * 2),
                                        int(NS.max_x - middle.real - NS.min_dist * 2)),
                       random.randrange(int(NS.min_y - middle.imag + NS.min_dist * 2),
                                        int(NS.max_y - middle.imag - NS.min_dist * 2)))
    line = line.translated(movement)
    start = line.point(0.5)

    c1 = line.start
    c2 = line.end
    v = start - c1
    if abs(v) < NS.min_dist * 2.5:
        dir_v_norm = v / abs(v)
        v = dir_v_norm * NS.min_dist * 2.5
        c1 = start - v
        c2 = start + v
    c1 = clamp_p1_to_border(c1)
    c2 = clamp_p1_to_border(c2)

    p2_1 = set_p2_between_p1_and_p3(c1, border1)
    p2_2 = set_p2_between_p1_and_p3(c2, border2)

    nc1 = NS.CC(CubicBezier(start, c1, p2_1, border1), [], additional=True)
    nc2 = NS.CC(CubicBezier(start, c2, p2_2, border2), [], additional=True)

    d = p3_2dist(border1)
    index = bisect_left([dist for (cords, dist) in NS.end_points], d)
    NS.end_points.insert(index, (border1, d))

    d = p3_2dist(border2)
    index = bisect_left([dist for (cords, dist) in NS.end_points], d)
    NS.end_points.insert(index, (border2, d))

    if NS.do_checks:
        update_all_intersections(nc1)
        update_all_intersections(nc2)
        if NS.early_checks:
            search_best_p1(nc1, nc2)
            nc1.curve.control2 = set_p2_between_p1_and_p3(nc1.curve.control1, nc1.curve.end)
            nc2.curve.control2 = set_p2_between_p1_and_p3(nc2.curve.control1, nc2.curve.end)
            calculate_penalty_and_optimize(nc1, NS.min_dist, NS.min_angle)
            calculate_penalty_and_optimize(nc2, NS.min_dist, NS.min_angle)

    for nc in [nc1, nc2]:
        NS.all_new_curves.append(nc)
        NS.all_bezier_curves.append(nc)


# generate new curve from control points p0, p1, p2, p3
def add_curve(p0, p1, p2, p3):
    new_cc = NS.CC(CubicBezier(p0, p1, p2, p3), [])
    if NS.do_checks:
        update_all_intersections(new_cc)
        if NS.early_checks:
            calculate_penalty_and_optimize(new_cc, NS.min_dist, NS.min_angle)
    NS.all_new_curves.append(new_cc)
    NS.all_bezier_curves.append(new_cc)


# connect ece with border by generating new connection curve
def connect_ECE_with_border(ece):
    if ece.is_connected():
        raise Exception(f"tried to connect ECE twice:\n\t{ece}")
    p0 = ece.start
    dir_p1 = p0 - ece.control1
    if dir_p1 == 0:
        raise Exception("dir_p1 == 0")
    if abs(dir_p1) < NS.min_dist:
        dir_p1_norm = dir_p1 / abs(dir_p1)
        dir_p1 = dir_p1_norm * NS.min_dist
    p1 = p0 + dir_p1
    p1 = clamp_p1_to_border(p1)
    p2, p3 = p2_and_p3_angled(ece.start, dir_p1, p1)

    add_curve(p0, p1, p2, p3)
    ece.connected = True


# split input curves at every non-G^1-continuous connection, append resulting ECEs to list
def calculate_ECEs_from_path():
    print(f"processing input file {NS.input_file}...")
    for path in NS.paths:
        first = segment_to_cubics(path[0])[0]
        last = segment_to_cubics(path[-1])[-1]
        add_ECE_if_ness(last, first)
        cc = NS.CC(first, [], True)
        NS.all_bezier_curves.append(cc)
        if NS.do_checks:
            update_all_intersections(cc)

        for i in range(1, len(path)):
            prev = segment_to_cubics(path[i - 1])[-1]
            curr = segment_to_cubics(path[i])[0]
            add_ECE_if_ness(prev, curr)
            cc = NS.CC(curr, [], True)
            NS.all_bezier_curves.append(cc)
            if NS.do_checks:
                update_all_intersections(cc)

    for oE in NS.open_ECEs:
        append_ECE(oE.control1, oE.start)


# try to interconnect ECEs via bridge curves
def calculate_bridge_curves():
    print("generate bridge curves...")
    number_bridge_curves = 0
    if NS.all_bridge_curves:
        # consider all ECEs
        copied_eces = NS.ECEs.copy()
        for ece in NS.ECEs:
            if ece.is_connected():
                continue
            number_bridge_curves += find_bridge_partner(ece, copied_eces)
    elif NS.inner_bridge_curves:
        # consider only inwards-pointing ECEs
        copied_inner_ECEs = NS.inner_ECEs.copy()
        for ece in NS.inner_ECEs:
            if ece.is_connected():
                continue
            number_bridge_curves += find_bridge_partner(ece, copied_inner_ECEs)


# calculate curves for every unconnected ECE that connects ECE with the border
def calculate_connection_curves():
    if NS.early_checks:
        print("generating and optimizing connection curves...")
    else:
        print("generating connection curves...")
    if not NS.bridge_curves_only:
        for i, ece in enumerate(NS.ECEs):
            if ece.is_connected():
                continue
            connect_ECE_with_border(ece)

    if NS.additional_curves:
        if NS.early_checks:
            print("generating and optimizing additional curves...")
        else:
            print("generating additional curves...")
        generate_additional_curves()

    if not NS.early_checks and NS.do_checks:
        print("optimizing generated curves...")
        for c in NS.all_bezier_curves:
            if not c.is_fixed:
                calculate_penalty_and_optimize(c, NS.min_dist, NS.min_angle)

    if not NS.additional_curves:
        for c in [NS.border_left, NS.border_right, NS.border_top, NS.border_bottom]:
            NS.all_bezier_curves.append(c)
            NS.all_new_curves.append(c)


# generate output file from input curves and all newly generated curves
def generate_output(output_file, filled):
    print(f"generating output file {output_file}...")
    if NS.additional_curves:
        for c in [NS.border_left, NS.border_right, NS.border_top, NS.border_bottom]:
            NS.all_bezier_curves.append(c)
            NS.all_new_curves.append(c)

    new_paths = [Path(curve.curve) for curve in NS.all_new_curves]
    new_attr_in = [{'stroke': 'black', 'fill': 'none', 'stroke-width': 0.5} for _ in NS.paths]
    new_attr = [{'stroke': 'black', 'fill': 'none', 'stroke-width': 0.5} for _ in new_paths]

    all_paths = NS.paths + new_paths
    all_attrs = new_attr_in + new_attr

    wsvg(all_paths, attributes=all_attrs, filename=output_file, viewbox=NS.view_box,
         dimensions=(NS.width, NS.height))

    if filled:
        splits = output_file.split(".")
        output_file = splits[0] + "_filled." + splits[1]
        new_attr_in = [{'stroke': 'black', 'fill': 'gray', 'stroke-width': 0.5} for _ in NS.paths]
        all_attrs = new_attr_in + new_attr
        wsvg(all_paths, attributes=all_attrs, filename=output_file, viewbox=NS.view_box,
             dimensions=(NS.width, NS.height))

