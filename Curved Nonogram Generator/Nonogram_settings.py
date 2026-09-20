import xml.etree.ElementTree as ET
from svgpathtools import CubicBezier, svg2paths


# ECE: Extendable Curve End
class ECE:
    def __init__(self, start: complex, control1: complex, connected=False):
        self.start = start
        self.control1 = control1
        self.connected = connected

    def __str__(self):
        return f"start: {self.start}, control1: {self.control1}"

    def __eq__(self, other):
        if not self.start == other.start or not self.control1 == other.control1:
            return False
        return True

    def is_connected(self):
        return self.connected


# CC = connection curve
class CC:
    def __init__(self, curve, intersections, fixed=False, additional=False, mate=None, angle_penalty=float('inf'),
                 intersection_penalty=float('inf'), distance_penalty=float('inf'), penalties=None):
        self.curve = curve
        self.intersections = intersections
        self.is_fixed = fixed
        self.is_additional = additional
        self.mate = mate
        self.penalties = [angle_penalty, intersection_penalty, distance_penalty] if penalties is None else penalties

    def __eq__(self, other):
        if other.curve.__eq__(self.curve):
            return True
        return False


input_file = ""

all_bridge_curves = True
inner_bridge_curves = False
bridge_curves_only = False
find_best_bridge = True

do_checks = False
early_checks = True
check_angles = True
check_intersections = True
check_distance = True

additional_curves = True

paths, attributes = [], []
view_box = None
min_x, min_y, max_x, max_y, vb_width, vb_height, width, height = 0, 0, 0, 0, 0, 0, 0, 0

all_bezier_curves, all_new_curves, ECEs, inner_ECEs, open_ECEs, end_points = [], [], [], [], [], []

border_left, border_top, border_right, border_bottom = (
    CC(CubicBezier(0, 0, 0, 0), []),
    CC(CubicBezier(0, 0, 0, 0), []),
    CC(CubicBezier(0, 0, 0, 0), []),
    CC(CubicBezier(0, 0, 0, 0), []),)

max_loops, min_angle, min_dist, min_gap, pen_threshold, w_angle, w_intersection, w_distance = 0, 0, 0, 0, 0, 0, 0, 0
alpha, beta, max_dist, control_point_dist = 0, 0, 0, 0


def set_switches(ab, ib, bco, fbb, dc, ec, ca, ci, cd, addc):
    global all_bridge_curves  # consider all ECEs for bridge curve generation
    global inner_bridge_curves  # consider only ECEs pointing inwards for bridge curve generation
    global bridge_curves_only  # generate only bridge curves, no connection curves
    global find_best_bridge  # find bridges with the lowest penalty

    global do_checks  # perform optimization
    global early_checks  # perform optimization immediately after curve generation
    global check_angles  # consider intersection angles during optimization
    global check_intersections  # consider distance between intersections during optimization
    global check_distance  # consider distance between curves during optimization

    global additional_curves  # generate additional curves

    all_bridge_curves = ab
    inner_bridge_curves = ib
    bridge_curves_only = bco
    find_best_bridge = fbb

    do_checks = dc
    early_checks = ec
    check_angles = ca
    check_intersections = ci
    check_distance = cd

    additional_curves = addc


def set_variables(ml, ma, md, mg, pt, wa, wi, wd, a, b, max_d, cpd):
    global max_loops  # maximum amount of iterations per curve during optimization
    global min_angle  # minimum allowed intersection angle
    global min_dist  # minimum allowed distance between intersections
    global min_gap  # minimum allowed gap between curves
    global pen_threshold  # minimum allowed penalty
    global w_angle  # weight of angle_penalty
    global w_intersection  # weight of intersection_penalty
    global w_distance  # weight of distance_penalty
    global alpha  # outer angle in radians to limit end points on border
    global beta  # inner angle in radians to limit end points on border
    global max_dist  # maximum allowed gap between end points on border for additional curve generation
    global control_point_dist  # distance between control points on additional curves

    max_loops = ml
    min_angle = ma
    min_dist = md
    min_gap = mg
    pen_threshold = pt

    w_angle = wa
    w_intersection = wi
    w_distance = wd

    alpha = a
    beta = b

    max_dist = max_d
    control_point_dist = cpd


def input_preprocess(input_name):
    global input_file, paths, attributes, view_box, min_x, min_y, max_x, max_y, vb_width, vb_height, width, height

    input_file = input_name
    paths_original, attributes_original = svg2paths(input_name)
    paths = paths_original.copy()
    attributes = attributes_original.copy
    tree = ET.parse(input_name)
    root = tree.getroot()
    try:
        width = float(root.get('width').replace('mm', ''))
    except Exception:
        try:
            width = float(root.get('width').replace('px', ''))
        except Exception:
            width = float(root.get('width').replace('pt', ''))
    try:
        height = float(root.get('height').replace('mm', ''))
    except Exception:
        try:
            height = float(root.get('height').replace('px', ''))
        except Exception:
            height = float(root.get('height').replace('pt', ''))
    view_box = root.get('viewBox')
    min_x, min_y, vb_width, vb_height = map(float, view_box.split())
    min_x = int(min_x)
    min_y = int(min_y)
    max_x = int(vb_width + min_x)
    max_y = int(vb_height + min_y)
    vb_width = int(vb_width)
    vb_height = int(vb_height)

    set_lists()

    return vb_width, vb_height


def set_lists():
    global all_bezier_curves, all_new_curves, ECEs, inner_ECEs, open_ECEs, end_points

    global border_left, border_top, border_right, border_bottom

    all_bezier_curves = []  # input curves, bridge curves, connection curves, additional curves
    all_new_curves = []  # bridge curves, connection curves, additional curves
    ECEs = []  # Extendable Curve Ends
    inner_ECEs = []  # ECEs pointing inwards
    open_ECEs = []  # if path is not continuous
    end_points = []  # end points on border

    # generate border of nonogram:
    top_left = complex(min_x, min_y)
    top_right = complex(max_x, min_y)
    bottom_left = complex(min_x, max_y)
    bottom_right = complex(max_x, max_y)
    border_left = CC(CubicBezier(top_left, top_left, bottom_left, bottom_left),
                     [], True, angle_penalty=0, intersection_penalty=0, distance_penalty=0)
    border_top = CC(CubicBezier(top_left, top_left, top_right, top_right),
                    [], True, angle_penalty=0, intersection_penalty=0, distance_penalty=0)
    border_right = CC(CubicBezier(top_right, top_right, bottom_right, bottom_right),
                      [], True, angle_penalty=0, intersection_penalty=0, distance_penalty=0)
    border_bottom = CC(CubicBezier(bottom_left, bottom_left, bottom_right, bottom_right),
                       [], True, angle_penalty=0, intersection_penalty=0, distance_penalty=0)
