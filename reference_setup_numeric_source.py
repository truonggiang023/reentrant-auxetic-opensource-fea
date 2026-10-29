# -*- coding: utf-8 -*-



from abaqus import *



from abaqusConstants import *



import regionToolset



import mesh



import math



import os







# ============================================================



# USER PARAMETERS



# ============================================================



MODEL_NAME = 'Parametric_ReEntrant_Core'







L1 = 12.0



L2 = L1



L_BEVEL = 5.66



CORE_WALL_T = 0.30







N_ROWS = 6



N_COLS = 5







CORE_HEIGHT = 8.0 * N_ROWS - 0.30



DEPTH = 25







# +1 re-entrant, 0 transition, -1 conventional honeycomb



LAMBDA_M = 1.0







MERGE_CORE = True



DELETE_TEMP_PARTS_AFTER_MERGE = True







SAVE_CAE = False



CAE_FILENAME = 'Parametric_ReEntrant_Core.cae'







# ============================================================



# AlSi10Mg MATERIAL / DAMAGE / ELEMENT SETTINGS



# ============================================================



# Units used here: mm - N - s - tonne - MPa



# Elastic constants below are editable defaults. Replace them if your



# tensile-test source provides different values.



ALSI10MG_E = 69000.0



ALSI10MG_NU = 0.30



ALSI10MG_DENSITY = 2.67e-9







ENABLE_DUCTILE_DAMAGE = True



DUCTILE_FRACTURE_STRAIN = 0.047



DUCTILE_TRIAXIALITY = 0.333



DUCTILE_STRAIN_RATE = 0.0



DAMAGE_FAILURE_DISPLACEMENT = 0.001







ASSIGN_C3D8R = True



GENERATE_CORE_MESH = True

# Global seed used for the core. Paper-reference edge counts below override
# selected long horizontal and inclined edges after this global seeding.
CORE_MESH_SIZE = 1.0

# Paper-reference mesh controls for the baseline validation model.
USE_PAPER_EDGE_SEEDS = True
PAPER_HORIZONTAL_EDGE_ELEMENTS = 19
PAPER_BEVEL_EDGE_ELEMENTS = 10
PAPER_TARGET_TOTAL_CORE_ELEMENTS = 10068
CORE_ELEMENT_COUNT = 0







# ============================================================



# QUASI-STATIC COMPRESSION VALIDATION SETTINGS



# ============================================================



# Compression direction = global Y (U2).

# Core-only model: steel platens contact the upper and lower core surfaces directly.

# For Dong et al. validation, set DEPTH to the paper specimen depth when reproducing it 1:1.



CREATE_STEEL_PLATENS = True







# Elastic steel assumptions. Replace if the source paper gives exact values.



STEEL_E = 210000.0            # MPa



STEEL_NU = 0.30



STEEL_DENSITY = 7.85e-9       # tonne/mm^3







# Steel platen dimensions. The platen contact face is fixed to 76 mm x 25 mm
# to match the cross-sectional area used to convert crosshead force to stress.
PLATEN_WIDTH_X = 76.0          # mm, global X direction
PLATEN_DEPTH_Z = 25.0          # mm, global Z direction
PLATEN_THICKNESS = 25.0        # mm, global Y direction
PLATEN_INITIAL_GAP = 0.0       # mm
REFERENCE_CROSS_SECTION_AREA = PLATEN_WIDTH_X * PLATEN_DEPTH_Z  # mm^2







GENERATE_PLATEN_MESH = True



PLATEN_MESH_SIZE = 2.0        # mm



ASSIGN_PLATEN_C3D8R = True







# General contact



FRICTION_COEFF = 0.35







# Loading target



COMPRESSION_STRAIN = 0.60







# Paper loading speed = 2 mm/min = 0.0333333 mm/s.



PAPER_COMPRESSION_SPEED_MM_MIN = 2.0



PAPER_COMPRESSION_SPEED_MM_S = PAPER_COMPRESSION_SPEED_MM_MIN / 60.0







# Explicit with the true paper speed can require a very long step.



# False = accelerated explicit step; verify quasi-static condition using



# ALLKE / ALLIE. True = use the physical 2 mm/min time.



USE_PAPER_PHYSICAL_SPEED = False



ACCELERATED_STEP_TIME = 1.0   # s; used only when physical speed = False

# Semi-automatic variable mass scaling: Whole Model, Throughout Step,
# target stable time increment = 1e-5 s, scale only elements below target,
# and re-evaluate every 1 increment.
ENABLE_MASS_SCALING = True
MASS_SCALING_TARGET_DT = 1.0e-5
MASS_SCALING_FREQUENCY = 1

# History output sampling interval in step time (s).
# With STEP_TIME = 1.0 s, 0.02 s gives about 50 intervals over the step.
HISTORY_OUTPUT_TIME_INTERVAL = 0.02

STEP_NAME = 'Compression'



AMPLITUDE_NAME = 'Amp_Compression'







# Paper states that out-of-plane displacement of the honeycomb was constrained.



# Here Z is the extrusion/depth direction, so front/back core faces receive U3=0.



CONSTRAIN_CORE_OUT_OF_PLANE = True







# Optional Abaqus job creation. It is not submitted automatically.



CREATE_JOB = True



JOB_NAME = 'Compression_Validation'



SUBMIT_JOB = False







# User-supplied true stress vs true plastic strain data.



# Abaqus *Plastic ordering is: (true stress MPa, true plastic strain).



PLASTIC_DATA = (



    (224.667594, 0.0),



    (226.1780406, 0.000219336),



    (227.78992, 0.000438623),



    (229.4632006, 0.000658854),



    (230.9756288, 0.000878045),



    (232.1861226, 0.001097189),



    (233.4173009, 0.001316284),



    (234.8308537, 0.001536322),



    (236.3759741, 0.001755321),



    (237.8006123, 0.001974272),



    (239.2460455, 0.002193175),



    (240.5812366, 0.00241302),



    (241.9066508, 0.002631827),



    (243.2831397, 0.002850586),



    (244.478304, 0.003069297),



    (245.6840758, 0.003287961),



    (246.9815872, 0.003507566),



    (248.2389497, 0.003726133),



    (249.577755, 0.003944653),



    (250.7451359, 0.004163125),



    (252.0852754, 0.004382538),



    (253.3751222, 0.004600915),



    (254.6958786, 0.004819243),



    (255.8045486, 0.005037525),



    (256.8835483, 0.005256746),



    (258.1450835, 0.005474931),



    (259.4071446, 0.005693069),



    (260.5481176, 0.00591116),



    (261.679685, 0.00613019),



    (262.8824363, 0.006348185),



    (264.1870987, 0.006566133),



    (265.7154586, 0.006784033),



    (267.1534126, 0.007002871),



    (268.2974198, 0.007220676),



    (269.4013, 0.007438434),



    (270.5665487, 0.007656144),



    (271.6002717, 0.007873807),



    (272.563591, 0.008092407),



    (273.5168774, 0.008309975),



    (274.5518412, 0.008527496),



    (275.688862, 0.008744969),



    (276.9587751, 0.008963379),



    (278.1476012, 0.009180758),



    (279.2860684, 0.009398089),



    (280.3639714, 0.009615373),



    (281.4629418, 0.009833593),



    (282.3382032, 0.010050783),



    (283.0916731, 0.010267925),



    (283.9879772, 0.01048502),



    (284.8849234, 0.01070305),



    (285.8430656, 0.010920051),



    (286.8321596, 0.011137005),



    (287.8522252, 0.011353911),



    (288.9137579, 0.011571752),



    (290.0366068, 0.011788564),



    (291.1497188, 0.01200533),



    (292.0695115, 0.012222048),



    (293.0100771, 0.01243872),



    (294.0839513, 0.012656324),



    (295.1579758, 0.012872902),



    (296.252853, 0.013089432),



    (297.3481765, 0.013305916),



    (298.3421193, 0.013523332),



    (299.2953163, 0.013739722),



    (300.2284652, 0.013956064),



    (301.2028641, 0.014172361),



    (302.1370688, 0.014389588),



    (303.0202394, 0.014605791),



    (303.8219619, 0.014821946),



    (304.6137751, 0.015038055),



    (305.3345932, 0.015255095),



    (306.1167875, 0.01547111),



    (306.7969505, 0.015687079),



    (307.4876149, 0.015903001),



    (308.1481345, 0.016119853),



    (308.880295, 0.016335682),



    (309.5512853, 0.016551464),



    (310.2122913, 0.016767199),



    (310.7816308, 0.016983864),



    (311.3918824, 0.017199506),



    (312.0228717, 0.017415102),



    (312.74639, 0.017630652),



    (313.6137778, 0.017846154),



    (314.4510417, 0.018062586),



    (315.1857252, 0.018277995),



    (316.0027985, 0.018493359),



    (316.8304596, 0.018708676),



    (317.6176881, 0.01892492),



    (318.3227616, 0.019140145),



    (319.0281091, 0.019355322),



    (319.7542785, 0.019570454),



    (320.5324231, 0.019786512),



    (321.2488916, 0.020001551),



    (321.8114306, 0.020216544),



    (322.4667264, 0.02043149),



    (323.1740131, 0.020647363),



    (323.8503988, 0.020862216),



    (324.3624151, 0.021077024),



    (324.8334552, 0.021291785),



    (325.4079217, 0.021507472),



    (326.0852494, 0.021722141),



    (326.6804521, 0.021936764),



    (327.0698713, 0.022151341),



    (327.418217, 0.022365872),



    (327.8185244, 0.022581328),



    (328.1671201, 0.022795766),



    (328.5364537, 0.023010159),



    (328.833735, 0.023224506),



    (329.1108052, 0.023439776),



    (329.521752, 0.023654031),



    (330.0463442, 0.02386824),



    (330.6536937, 0.024082403),



    (331.1686907, 0.024297489),



    (331.6732327, 0.02451156),



    (332.2605763, 0.024725585),



    (332.85847, 0.024939564),



    (333.4465805, 0.025154466),



    (333.9415871, 0.025368354),



    (334.3230806, 0.025582196),



    (334.6840307, 0.025795992),



    (335.1178104, 0.02601071),



    (335.5824486, 0.026224414),



    (335.9748418, 0.026438073),



    (336.3156376, 0.026651687),



    (336.7086192, 0.026866221),



    (337.1221154, 0.027079743),



    (337.5461111, 0.027293219),



    (337.9599013, 0.02750665),



    (338.3427634, 0.027720035),



    (338.7157249, 0.02793434),



    (339.1092128, 0.028147634),



    (339.5132022, 0.028360882),



    (339.865499, 0.028574085),



    (340.093815, 0.028788207),



    (340.2907556, 0.029001319),



    (340.5188697, 0.029214386),



    (340.7470501, 0.029427407),



    (340.9756253, 0.029641347),



    (341.2143185, 0.029854277),



    (341.4530824, 0.030067162),



    (341.7438396, 0.030280001),



    (342.0454056, 0.030493759),



    (342.4402403, 0.030706507),



    (342.8456031, 0.030919211),



    (343.2511074, 0.031131869),



    (343.6051057, 0.031345444),



    (343.9069043, 0.031558012),



    (344.1464, 0.031770535),



    (344.3547598, 0.031983012),



    (344.5215595, 0.032195444),



    (344.6887301, 0.032408793),



    (344.8764271, 0.032621134),



    (345.1058169, 0.032833431),



    (345.3456863, 0.033045683),



    (345.5547117, 0.03325885),



    (345.7738807, 0.033471011),



    (346.0139514, 0.033683128),



    (346.2436706, 0.033895199),



    (346.4529399, 0.034108185),



    (346.6827875, 0.034320166),



    (346.8501288, 0.034532103),



    (347.007079, 0.034743994),



    (347.1539641, 0.034956799),



    (347.290112, 0.035168601),



    (347.4158489, 0.035380358),



    (347.6042467, 0.03559207),



    (347.7721422, 0.035804694),



    (347.9710773, 0.036016316),



    (348.1491724, 0.036227894),



    (348.3168631, 0.036439426),



    (348.4218883, 0.036650914),



    (348.5690726, 0.036863314),



    (348.6845885, 0.037074712),



    (348.7896646, 0.037286065),



    (348.8738344, 0.037497374),



    (348.9583423, 0.037709594),



    (349.0215926, 0.037920814),



    (351.151075, 0.038144409),



)











# ============================================================



# DERIVED PARAMETERS / CHECKS



# ============================================================



if N_ROWS < 1 or N_COLS < 1:



    raise ValueError('N_ROWS and N_COLS must be >= 1.')







if L1 <= 0.0 or L2 < 0.0 or L_BEVEL <= 0.0 or CORE_WALL_T <= 0.0:



    raise ValueError('Invalid geometric input.')







if LAMBDA_M < -1.0 or LAMBDA_M > 1.0:



    raise ValueError('LAMBDA_M must be in [-1, 1].')







CELL_H = CORE_HEIGHT / float(N_ROWS)



HALF_H = 0.5 * CELL_H



HALF_T = 0.5 * CORE_WALL_T







if L_BEVEL <= HALF_H:



    raise ValueError('L_BEVEL must be greater than CELL_H/2.')







DX_REF = math.sqrt(L_BEVEL**2 - HALF_H**2)



WAIST_SHIFT = LAMBDA_M * DX_REF







# Same formulation as your CATScript:



# L2 remains the waist-to-waist ligament length.



PITCH_X = L1 + L2 - 2.0 * WAIST_SHIFT







if PITCH_X <= 0.0:



    raise ValueError('Calculated PITCH_X <= 0.')







ACTUAL_INCLINED_LENGTH = math.sqrt(HALF_H**2 + WAIST_SHIFT**2)







if abs(WAIST_SHIFT) < 1.0e-12:



    THETA_DEG = 90.0



else:



    THETA_DEG = math.degrees(math.atan(HALF_H / abs(WAIST_SHIFT)))







THROAT_CLEARANCE = L1 - 2.0 * WAIST_SHIFT - 2.0 * CORE_WALL_T







if LAMBDA_M > 0.0 and THROAT_CLEARANCE <= 0.0:



    raise ValueError('Invalid re-entrant geometry: throat clearance <= 0.')











# ============================================================



# GEOMETRY HELPERS



# ============================================================



def line_intersection(q1, d1, q2, d2):



    det = d1[0] * d2[1] - d1[1] * d2[0]



    if abs(det) < 1.0e-12:



        return None







    rx = q2[0] - q1[0]



    ry = q2[1] - q1[1]



    u = (rx * d2[1] - ry * d2[0]) / det







    return (q1[0] + u*d1[0], q1[1] + u*d1[1])











def offset_closed_loop(points, offset):



    result = []



    n = len(points)







    for i in range(n):



        p_prev = points[(i - 1) % n]



        p_cur = points[i]



        p_next = points[(i + 1) % n]







        d1x = p_cur[0] - p_prev[0]



        d1y = p_cur[1] - p_prev[1]



        d2x = p_next[0] - p_cur[0]



        d2y = p_next[1] - p_cur[1]







        len1 = math.sqrt(d1x*d1x + d1y*d1y)



        len2 = math.sqrt(d2x*d2x + d2y*d2y)







        if len1 < 1.0e-12 or len2 < 1.0e-12:



            raise ValueError('Zero-length segment in centerline.')







        n1x, n1y = -d1y/len1, d1x/len1



        n2x, n2y = -d2y/len2, d2x/len2







        q1 = (p_cur[0] + offset*n1x, p_cur[1] + offset*n1y)



        q2 = (p_cur[0] + offset*n2x, p_cur[1] + offset*n2y)







        inter = line_intersection(q1, (d1x, d1y), q2, (d2x, d2y))







        if inter is None:



            result.append((



                p_cur[0] + offset*0.5*(n1x+n2x),



                p_cur[1] + offset*0.5*(n1y+n2y)



            ))



        else:



            result.append(inter)







    return result











def build_column_centerline(x0):



    x1 = x0 + L1



    pts = [(x0, 0.0)]







    for r in range(N_ROWS):



        z0 = r * CELL_H



        zm = z0 + HALF_H



        z1 = z0 + CELL_H



        pts.append((x0 + WAIST_SHIFT, zm))



        pts.append((x0, z1))







    pts.append((x1, CORE_HEIGHT))







    for r in reversed(range(N_ROWS)):



        z0 = r * CELL_H



        zm = z0 + HALF_H



        pts.append((x1 - WAIST_SHIFT, zm))



        pts.append((x1, z0))







    return pts











def draw_closed_loop(sketch, points):



    for i in range(len(points)):



        p1 = points[i]



        p2 = points[(i + 1) % len(points)]







        if abs(p2[0]-p1[0]) > 1.0e-12 or abs(p2[1]-p1[1]) > 1.0e-12:



            sketch.Line(point1=p1, point2=p2)











def update_bounds(bounds, points):



    for x, y in points:



        bounds[0] = min(bounds[0], x)



        bounds[1] = min(bounds[1], y)



        bounds[2] = max(bounds[2], x)



        bounds[3] = max(bounds[3], y)











def make_ring_part(model, name, centerline, depth, bounds):



    loop_a = offset_closed_loop(centerline, +0.5*CORE_WALL_T)



    loop_b = offset_closed_loop(centerline, -0.5*CORE_WALL_T)







    update_bounds(bounds, loop_a)



    update_bounds(bounds, loop_b)







    sk_name = '__sk_' + name



    sk = model.ConstrainedSketch(name=sk_name, sheetSize=1000.0)







    draw_closed_loop(sk, loop_a)



    draw_closed_loop(sk, loop_b)







    p = model.Part(name=name, dimensionality=THREE_D, type=DEFORMABLE_BODY)



    p.BaseSolidExtrude(sketch=sk, depth=depth)







    del model.sketches[sk_name]



    return p











def make_rect_part(model, name, x0, y0, x1, y1, depth, bounds=None):



    xa, xb = min(x0, x1), max(x0, x1)



    ya, yb = min(y0, y1), max(y0, y1)







    if xb-xa <= 1.0e-12 or yb-ya <= 1.0e-12:



        raise ValueError('Invalid rectangle: %s' % name)







    if bounds is not None:



        update_bounds(bounds, [(xa,ya), (xb,ya), (xb,yb), (xa,yb)])







    sk_name = '__sk_' + name



    sk = model.ConstrainedSketch(name=sk_name, sheetSize=1000.0)



    sk.rectangle(point1=(xa,ya), point2=(xb,yb))







    p = model.Part(name=name, dimensionality=THREE_D, type=DEFORMABLE_BODY)



    p.BaseSolidExtrude(sketch=sk, depth=depth)







    del model.sketches[sk_name]



    return p











# ============================================================
# PAPER-REFERENCE EDGE SEEDING HELPER
# ============================================================

def apply_paper_reference_edge_seeds(part):
    """Apply paper edge counts on the baseline validation geometry.

    Long horizontal XY edges receive 19 elements and inclined/bevel XY edges
    receive 10 elements. The global 1 mm seed remains on the other edges.
    """
    horizontal_edges = []
    bevel_edges = []
    tol = 1.0e-6
    horizontal_min_length = 0.40 * min(L1, max(L2, L1))
    bevel_min_length = 0.40 * ACTUAL_INCLINED_LENGTH
    bevel_max_length = 1.60 * ACTUAL_INCLINED_LENGTH

    for edge in part.edges:
        try:
            vertex_ids = edge.getVertices()
            if len(vertex_ids) != 2:
                continue
            p0 = part.vertices[vertex_ids[0]].pointOn[0]
            p1 = part.vertices[vertex_ids[1]].pointOn[0]
            dx = p1[0] - p0[0]
            dy = p1[1] - p0[1]
            dz = p1[2] - p0[2]
            edge_length = math.sqrt(dx*dx + dy*dy + dz*dz)

            # Only classify edges that lie in an XY end plane. Depth-direction
            # and wall-thickness edges remain controlled by the global seed.
            if abs(dz) > tol:
                continue
            if abs(dy) <= tol and abs(dx) > tol and edge_length >= horizontal_min_length:
                horizontal_edges.append(edge)
            elif abs(dx) > tol and abs(dy) > tol:
                if bevel_min_length <= edge_length <= bevel_max_length:
                    bevel_edges.append(edge)
        except:
            pass

    if len(horizontal_edges) > 0:
        part.seedEdgeByNumber(
            edges=tuple(horizontal_edges),
            number=PAPER_HORIZONTAL_EDGE_ELEMENTS,
            constraint=FIXED
        )
    if len(bevel_edges) > 0:
        part.seedEdgeByNumber(
            edges=tuple(bevel_edges),
            number=PAPER_BEVEL_EDGE_ELEMENTS,
            constraint=FIXED
        )

    print('Paper-reference seeded horizontal edges = %d' % len(horizontal_edges))
    print('Paper-reference seeded bevel edges      = %d' % len(bevel_edges))


# ============================================================



# CREATE MODEL



# ============================================================



if MODEL_NAME in mdb.models.keys():



    del mdb.models[MODEL_NAME]







model = mdb.Model(name=MODEL_NAME)



assembly = model.rootAssembly



assembly.DatumCsysByDefault(CARTESIAN)







core_instances = []



temp_part_names = []







core_bounds = [1.0e30, 1.0e30, -1.0e30, -1.0e30]











# ============================================================



# 1) CONTINUOUS COLUMN RINGS



# ============================================================



for c in range(N_COLS):



    x0 = c * PITCH_X



    centerline = build_column_centerline(x0)







    part_name = 'TMP_Column_%03d' % (c+1)



    p = make_ring_part(model, part_name, centerline, DEPTH, core_bounds)







    inst = assembly.Instance(



        name=part_name + '-1',



        part=p,



        dependent=ON



    )







    core_instances.append(inst)



    temp_part_names.append(part_name)











# ============================================================



# 2) SHARED HORIZONTAL WALLS



# ============================================================



for c in range(N_COLS):



    x0 = c * PITCH_X



    x1 = x0 + L1







    for level in range(1, N_ROWS):



        y = level * CELL_H



        part_name = 'TMP_Shared_C%03d_L%03d' % (c+1, level)







        p = make_rect_part(



            model, part_name,



            x0, y-HALF_T,



            x1, y+HALF_T,



            DEPTH,



            core_bounds



        )







        inst = assembly.Instance(



            name=part_name + '-1',



            part=p,



            dependent=ON



        )







        core_instances.append(inst)



        temp_part_names.append(part_name)











# ============================================================



# 3) LIGAMENTS BETWEEN COLUMNS



# ============================================================



for c in range(N_COLS-1):



    x0 = c * PITCH_X



    x1 = x0 + L1



    next_x0 = (c+1) * PITCH_X







    for r in range(N_ROWS):



        y = r * CELL_H + HALF_H







        right_waist = x1 - WAIST_SHIFT



        next_left_waist = next_x0 + WAIST_SHIFT







        ligament_length = next_left_waist - right_waist







        if abs(ligament_length - L2) > 1.0e-6:



            raise RuntimeError(



                'Ligament length != L2: %.9f vs %.9f'



                % (ligament_length, L2)



            )







        part_name = 'TMP_Lig_C%03d_R%03d' % (c+1, r+1)







        p = make_rect_part(



            model, part_name,



            right_waist, y-HALF_T,



            next_left_waist, y+HALF_T,



            DEPTH,



            core_bounds



        )







        inst = assembly.Instance(



            name=part_name + '-1',



            part=p,



            dependent=ON



        )







        core_instances.append(inst)



        temp_part_names.append(part_name)











assembly.regenerate()











# ============================================================



# BOOLEAN MERGE -> ONE CORE PART



# ============================================================



if MERGE_CORE:



    assembly.InstanceFromBooleanMerge(



        name='Core_AlSi10Mg',



        instances=tuple(core_instances),



        keepIntersections=OFF,



        originalInstances=DELETE,



        domain=GEOMETRY



    )







    core_part = model.parts['Core_AlSi10Mg']







    if len(core_part.cells) > 0:



        core_part.Set(cells=core_part.cells[:], name='ALL_CORE')







    if DELETE_TEMP_PARTS_AFTER_MERGE:



        for part_name in temp_part_names:



            if part_name in model.parts.keys():



                del model.parts[part_name]



else:



    core_part = None











# ============================================================



# AlSi10Mg MATERIAL + PLASTICITY + DUCTILE DAMAGE



# ============================================================



if MERGE_CORE and core_part is not None:







    al_mat = model.Material(name='AlSi10Mg')







    al_mat.Density(



        table=((ALSI10MG_DENSITY,),)



    )







    al_mat.Elastic(



        table=((ALSI10MG_E, ALSI10MG_NU),)



    )







    al_mat.Plastic(



        table=PLASTIC_DATA



    )







    if ENABLE_DUCTILE_DAMAGE:



        # Simplified single-point initiation definition from the value supplied:



        # fracture plastic strain = 0.047 at stress triaxiality = 0.333.



        # If you obtain the full fracture-strain-vs-triaxiality table from the



        # source paper, replace this one-point table for 1:1 validation.



        al_mat.DuctileDamageInitiation(



            table=((



                DUCTILE_FRACTURE_STRAIN,



                DUCTILE_TRIAXIALITY,



                DUCTILE_STRAIN_RATE



            ),)



        )







        al_mat.ductileDamageInitiation.DamageEvolution(



            type=DISPLACEMENT,



            softening=LINEAR,



            table=((DAMAGE_FAILURE_DISPLACEMENT,),)



        )







    model.HomogeneousSolidSection(



        name='AlSi10Mg_Section',



        material='AlSi10Mg',



        thickness=None



    )







    core_region = regionToolset.Region(



        cells=core_part.cells[:]



    )







    core_part.SectionAssignment(



        region=core_region,



        sectionName='AlSi10Mg_Section',



        offset=0.0,



        offsetType=MIDDLE_SURFACE,



        offsetField='',



        thicknessAssignment=FROM_SECTION



    )







    # Assign C3D8R as the intended element formulation.



    # Mesh generation is optional because automatic all-hex meshing can depend



    # on Abaqus version and the final Boolean topology.



    if ASSIGN_C3D8R:



        elem_type = mesh.ElemType(



            elemCode=C3D8R,



            elemLibrary=EXPLICIT



        )







        core_part.setElementType(



            regions=(core_part.cells[:],),



            elemTypes=(elem_type,)



        )







    if GENERATE_CORE_MESH:



        core_part.seedPart(



            size=CORE_MESH_SIZE,



            deviationFactor=0.1,



            minSizeFactor=0.1



        )

        if USE_PAPER_EDGE_SEEDS:
            apply_paper_reference_edge_seeds(core_part)







        try:



            core_part.setMeshControls(



                regions=core_part.cells[:],



                elemShape=HEX,



                technique=SWEEP



            )



        except:



            print('WARNING: automatic HEX/SWEEP controls could not be applied.')



            print('         Partition the core manually before meshing if needed.')

        core_part.generateMesh()
        CORE_ELEMENT_COUNT = len(core_part.elements)
        print('Generated core elements       = %d' % CORE_ELEMENT_COUNT)
        if USE_PAPER_EDGE_SEEDS and CORE_ELEMENT_COUNT != PAPER_TARGET_TOTAL_CORE_ELEMENTS:
            print('WARNING: paper reports %d C3D8R elements; current mesh has %d.'
                  % (PAPER_TARGET_TOTAL_CORE_ELEMENTS, CORE_ELEMENT_COUNT))
            print('         Exact total count depends on partitions/sweep topology, not only edge seeds.')











# ============================================================



# CORE BOUNDS



# ============================================================



CORE_X_MIN = core_bounds[0]



CORE_Y_MIN = core_bounds[1]



CORE_X_MAX = core_bounds[2]



CORE_Y_MAX = core_bounds[3]



assembly.regenerate()



# ============================================================



# STEEL PLATENS + COMPRESSION ASSEMBLY



# ============================================================



#



# Compression is along global Y.



#



# Contact levels for the core-only specimen:

#   bottom = CORE_Y_MIN

#   top    = CORE_Y_MAX

#

# ============================================================



SPECIMEN_BOTTOM_Y = CORE_Y_MIN

SPECIMEN_TOP_Y = CORE_Y_MAX



SPECIMEN_HEIGHT = SPECIMEN_TOP_Y - SPECIMEN_BOTTOM_Y







if SPECIMEN_HEIGHT <= 0.0:



    raise RuntimeError('Calculated specimen compression height <= 0.')







COMPRESSION_DISPLACEMENT = COMPRESSION_STRAIN * SPECIMEN_HEIGHT







if USE_PAPER_PHYSICAL_SPEED:



    STEP_TIME = COMPRESSION_DISPLACEMENT / PAPER_COMPRESSION_SPEED_MM_S



else:



    STEP_TIME = ACCELERATED_STEP_TIME







if STEP_TIME <= 0.0:



    raise ValueError('STEP_TIME must be > 0.')







EFFECTIVE_COMPRESSION_SPEED = COMPRESSION_DISPLACEMENT / STEP_TIME







CORE_X_MID = 0.5 * (CORE_X_MIN + CORE_X_MAX)
CORE_WIDTH_ACTUAL = CORE_X_MAX - CORE_X_MIN

# Fixed 76 mm platen width, centered on the core.
PLATEN_X_MIN = CORE_X_MID - 0.5 * PLATEN_WIDTH_X



PLATEN_X_MAX = CORE_X_MID + 0.5 * PLATEN_WIDTH_X







BOTTOM_INNER_Y = SPECIMEN_BOTTOM_Y - PLATEN_INITIAL_GAP



BOTTOM_OUTER_Y = BOTTOM_INNER_Y - PLATEN_THICKNESS







TOP_INNER_Y = SPECIMEN_TOP_Y + PLATEN_INITIAL_GAP



TOP_OUTER_Y = TOP_INNER_Y + PLATEN_THICKNESS







PLATEN_X_MID = CORE_X_MID



PLATEN_Z_MID = 0.5 * PLATEN_DEPTH_Z

if CORE_WIDTH_ACTUAL > PLATEN_WIDTH_X + 1.0e-6:
    print('WARNING: core width %.6f mm exceeds the 76 mm platen width.' % CORE_WIDTH_ACTUAL)
if abs(DEPTH - PLATEN_DEPTH_Z) > 1.0e-9:
    print('WARNING: core depth %.6f mm differs from platen depth %.6f mm.'
          % (DEPTH, PLATEN_DEPTH_Z))











if CREATE_STEEL_PLATENS:







    # --------------------------------------------------------



    # Steel material and section



    # --------------------------------------------------------



    steel_mat = model.Material(name='Steel_Elastic')



    steel_mat.Density(



        table=((STEEL_DENSITY,),)



    )



    steel_mat.Elastic(



        table=((STEEL_E, STEEL_NU),)



    )







    model.HomogeneousSolidSection(



        name='Steel_Section',



        material='Steel_Elastic',



        thickness=None



    )







    # --------------------------------------------------------



    # Bottom platen



    # --------------------------------------------------------



    bottom_platen = make_rect_part(



        model,



        'Bottom_Platen_Steel',



        PLATEN_X_MIN,



        BOTTOM_OUTER_Y,



        PLATEN_X_MAX,



        BOTTOM_INNER_Y,



        PLATEN_DEPTH_Z



    )







    bottom_platen.SectionAssignment(



        region=regionToolset.Region(cells=bottom_platen.cells[:]),



        sectionName='Steel_Section',



        offset=0.0,



        offsetType=MIDDLE_SURFACE,



        offsetField='',



        thicknessAssignment=FROM_SECTION



    )







    if ASSIGN_PLATEN_C3D8R:



        platen_elem_type = mesh.ElemType(



            elemCode=C3D8R,



            elemLibrary=EXPLICIT



        )







        bottom_platen.setElementType(



            regions=(bottom_platen.cells[:],),



            elemTypes=(platen_elem_type,)



        )







    if GENERATE_PLATEN_MESH:



        bottom_platen.seedPart(



            size=PLATEN_MESH_SIZE,



            deviationFactor=0.1,



            minSizeFactor=0.1



        )







        try:



            bottom_platen.setMeshControls(



                regions=bottom_platen.cells[:],



                elemShape=HEX,



                technique=SWEEP



            )



        except:



            pass







        bottom_platen.generateMesh()







    bottom_platen_inst = assembly.Instance(



        name='Bottom_Platen_Steel-1',



        part=bottom_platen,



        dependent=ON



    )







    # --------------------------------------------------------



    # Top platen



    # --------------------------------------------------------



    top_platen = make_rect_part(



        model,



        'Top_Platen_Steel',



        PLATEN_X_MIN,



        TOP_INNER_Y,



        PLATEN_X_MAX,



        TOP_OUTER_Y,



        PLATEN_DEPTH_Z



    )







    top_platen.SectionAssignment(



        region=regionToolset.Region(cells=top_platen.cells[:]),



        sectionName='Steel_Section',



        offset=0.0,



        offsetType=MIDDLE_SURFACE,



        offsetField='',



        thicknessAssignment=FROM_SECTION



    )







    if ASSIGN_PLATEN_C3D8R:



        platen_elem_type_top = mesh.ElemType(



            elemCode=C3D8R,



            elemLibrary=EXPLICIT



        )







        top_platen.setElementType(



            regions=(top_platen.cells[:],),



            elemTypes=(platen_elem_type_top,)



        )







    if GENERATE_PLATEN_MESH:



        top_platen.seedPart(



            size=PLATEN_MESH_SIZE,



            deviationFactor=0.1,



            minSizeFactor=0.1



        )







        try:



            top_platen.setMeshControls(



                regions=top_platen.cells[:],



                elemShape=HEX,



                technique=SWEEP



            )



        except:



            pass







        top_platen.generateMesh()







    top_platen_inst = assembly.Instance(



        name='Top_Platen_Steel-1',



        part=top_platen,



        dependent=ON



    )







    assembly.regenerate()







    # --------------------------------------------------------



    # Create platen inner / outer surfaces



    # --------------------------------------------------------



    bottom_inner_face = bottom_platen_inst.faces.findAt(



        ((PLATEN_X_MID, BOTTOM_INNER_Y, PLATEN_Z_MID),)



    )



    bottom_outer_face = bottom_platen_inst.faces.findAt(



        ((PLATEN_X_MID, BOTTOM_OUTER_Y, PLATEN_Z_MID),)



    )







    top_inner_face = top_platen_inst.faces.findAt(



        ((PLATEN_X_MID, TOP_INNER_Y, PLATEN_Z_MID),)



    )



    top_outer_face = top_platen_inst.faces.findAt(



        ((PLATEN_X_MID, TOP_OUTER_Y, PLATEN_Z_MID),)



    )







    assembly.Surface(



        name='BOTTOM_PLATEN_INNER',



        side1Faces=(bottom_inner_face,)



    )



    assembly.Surface(



        name='BOTTOM_PLATEN_OUTER',



        side1Faces=(bottom_outer_face,)



    )



    assembly.Surface(



        name='TOP_PLATEN_INNER',



        side1Faces=(top_inner_face,)



    )



    assembly.Surface(



        name='TOP_PLATEN_OUTER',



        side1Faces=(top_outer_face,)



    )







    # --------------------------------------------------------



    # Reference points



    # --------------------------------------------------------



    bottom_rp_feature = assembly.ReferencePoint(



        point=(



            PLATEN_X_MID,



            BOTTOM_OUTER_Y,



            PLATEN_Z_MID



        )



    )



    bottom_rp = assembly.referencePoints[bottom_rp_feature.id]







    top_rp_feature = assembly.ReferencePoint(



        point=(



            PLATEN_X_MID,



            TOP_OUTER_Y,



            PLATEN_Z_MID



        )



    )



    top_rp = assembly.referencePoints[top_rp_feature.id]







    assembly.Set(



        name='RP_BOTTOM_PLATEN',



        referencePoints=(bottom_rp,)



    )



    assembly.Set(



        name='RP_TOP_PLATEN',



        referencePoints=(top_rp,)



    )







    # --------------------------------------------------------



    # Couple each outer platen surface to its reference point.



    #



    # This imposes the loading/support uniformly through the outer



    # surface while leaving the inner steel contact surface deformable.



    # --------------------------------------------------------



    model.Coupling(



        name='Coupling_Bottom_Platen',



        controlPoint=assembly.sets['RP_BOTTOM_PLATEN'],



        surface=assembly.surfaces['BOTTOM_PLATEN_OUTER'],



        influenceRadius=WHOLE_SURFACE,



        couplingType=KINEMATIC,



        localCsys=None,



        u1=ON,



        u2=ON,



        u3=ON,



        ur1=ON,



        ur2=ON,



        ur3=ON



    )







    model.Coupling(



        name='Coupling_Top_Platen',



        controlPoint=assembly.sets['RP_TOP_PLATEN'],



        surface=assembly.surfaces['TOP_PLATEN_OUTER'],



        influenceRadius=WHOLE_SURFACE,



        couplingType=KINEMATIC,



        localCsys=None,



        u1=ON,



        u2=ON,



        u3=ON,



        ur1=ON,



        ur2=ON,



        ur3=ON



    )











# ============================================================



# EXPLICIT COMPRESSION STEP



# ============================================================



model.ExplicitDynamicsStep(



    name=STEP_NAME,



    previous='Initial',



    timePeriod=STEP_TIME,



    improvedDtMethod=ON



)

if ENABLE_MASS_SCALING:
    model.steps[STEP_NAME].setValues(
        massScaling=((
            SEMI_AUTOMATIC,
            MODEL,
            THROUGHOUT_STEP,
            0.0,
            MASS_SCALING_TARGET_DT,
            BELOW_MIN,
            MASS_SCALING_FREQUENCY,
            0,
            0.0,
            0.0,
            0,
            None
        ),)
    )







model.SmoothStepAmplitude(



    name=AMPLITUDE_NAME,



    timeSpan=STEP,



    data=(



        (0.0, 0.0),



        (STEP_TIME, 1.0)



    )



)











# ============================================================



# GENERAL CONTACT



# ============================================================



#



# The paper specifies self-contact and mu = 0.35 between the



# honeycomb and the steel indenter/pedestal.



#



# This implementation assigns mu = 0.35 globally in General Contact,



# so the same tangential property is also used for core self-contact.



# If a separate self-contact friction law is obtained from the paper,



# replace the global assignment with pair-specific assignments.



# ============================================================







contact_prop = model.ContactProperty(



    'GeneralContactProperty'



)







contact_prop.NormalBehavior(



    pressureOverclosure=HARD,



    allowSeparation=ON,



    constraintEnforcementMethod=DEFAULT



)







contact_prop.TangentialBehavior(



    formulation=PENALTY,



    directionality=ISOTROPIC,



    slipRateDependency=OFF,



    pressureDependency=OFF,



    temperatureDependency=OFF,



    dependencies=0,



    table=((FRICTION_COEFF,),),



    shearStressLimit=None,



    maximumElasticSlip=FRACTION,



    fraction=0.005,



    elasticSlipStiffness=None



)







model.ContactExp(



    name='General_Contact',



    createStepName=STEP_NAME



)







model.interactions['General_Contact'].includedPairs.setValuesInStep(



    stepName=STEP_NAME,



    useAllstar=ON



)







model.interactions['General_Contact'].contactPropertyAssignments.appendInStep(



    stepName=STEP_NAME,



    assignments=(



        (GLOBAL, SELF, 'GeneralContactProperty'),



    )



)











# ============================================================



# BOUNDARY CONDITIONS



# ============================================================







if CREATE_STEEL_PLATENS:







    # Fixed lower platen



    model.EncastreBC(



        name='BC_Bottom_Platen_Fixed',



        createStepName='Initial',



        region=assembly.sets['RP_BOTTOM_PLATEN']



    )







    # Top platen:



    # U1 = U3 = 0 and rotations fixed.



    # U2 is prescribed in the compression step.



    model.DisplacementBC(



        name='BC_Top_Platen_Compression',



        createStepName=STEP_NAME,



        region=assembly.sets['RP_TOP_PLATEN'],



        u1=0.0,



        u2=-COMPRESSION_DISPLACEMENT,



        u3=0.0,



        ur1=0.0,



        ur2=0.0,



        ur3=0.0,



        amplitude=AMPLITUDE_NAME,



        fixed=OFF,



        distributionType=UNIFORM,



        fieldName='',



        localCsys=None



    )











# ============================================================



# OUT-OF-PLANE CONSTRAINT FOR CORE



# ============================================================



if CONSTRAIN_CORE_OUT_OF_PLANE and MERGE_CORE:



    core_instance = assembly.instances['Core_AlSi10Mg-1']



    # Select the two extrusion-end faces independently.

    # Do not combine the FaceArray objects into a Python tuple and pass that

    # tuple to assembly.Set(); some Abaqus versions reject that sequence and

    # raise "Feature creation failed".  Apply the U3=0 condition directly

    # to geometry Regions instead.

    z_tol = max(1.0e-5, DEPTH * 1.0e-7)



    front_faces = core_instance.faces.getByBoundingBox(

        zMin=-z_tol,

        zMax=z_tol

    )



    back_faces = core_instance.faces.getByBoundingBox(

        zMin=DEPTH-z_tol,

        zMax=DEPTH+z_tol

    )



    if len(front_faces) > 0:

        front_region = regionToolset.Region(faces=front_faces)

        model.DisplacementBC(

            name='BC_Core_OutOfPlane_Front',

            createStepName='Initial',

            region=front_region,

            u1=UNSET,

            u2=UNSET,

            u3=0.0,

            ur1=UNSET,

            ur2=UNSET,

            ur3=UNSET,

            amplitude=UNSET,

            fixed=OFF,

            distributionType=UNIFORM,

            fieldName='',

            localCsys=None

        )

    else:

        print('WARNING: no front core faces found for U3=0 constraint.')



    if len(back_faces) > 0:

        back_region = regionToolset.Region(faces=back_faces)

        model.DisplacementBC(

            name='BC_Core_OutOfPlane_Back',

            createStepName='Initial',

            region=back_region,

            u1=UNSET,

            u2=UNSET,

            u3=0.0,

            ur1=UNSET,

            ur2=UNSET,

            ur3=UNSET,

            amplitude=UNSET,

            fixed=OFF,

            distributionType=UNIFORM,

            fieldName='',

            localCsys=None

        )

    else:

        print('WARNING: no back core faces found for U3=0 constraint.')





# ============================================================



# OUTPUT REQUESTS



# ============================================================







# Field output for deformation, plasticity and damage.



if 'F-Output-1' in model.fieldOutputRequests.keys():



    model.fieldOutputRequests['F-Output-1'].setValues(



        variables=(



            'S',



            'LE',



            'PEEQ',



            'U',



            'V',



            'A',



            'RF',



            'STATUS',



            'SDEG'



        )



    )







# Platen reaction-force / displacement histories.



if CREATE_STEEL_PLATENS:







    model.HistoryOutputRequest(



        name='H_Top_Platen_RF_U',



        createStepName=STEP_NAME,



        variables=(



            'RF2',



            'U2',



            'V2'



        ),



        timeInterval=HISTORY_OUTPUT_TIME_INTERVAL,



        region=assembly.sets['RP_TOP_PLATEN']



    )







    model.HistoryOutputRequest(



        name='H_Bottom_Platen_RF_U',



        createStepName=STEP_NAME,



        variables=(



            'RF2',



            'U2'



        ),



        timeInterval=HISTORY_OUTPUT_TIME_INTERVAL,



        region=assembly.sets['RP_BOTTOM_PLATEN']



    )







# Whole-model energy outputs.



# Use ALLKE / ALLIE to assess quasi-static quality of the Explicit run.



model.HistoryOutputRequest(



    name='H_Energy',



    createStepName=STEP_NAME,



    variables=(



        'ALLKE',



        'ALLIE',



        'ALLAE',



        'ETOTAL'



    ),



    timeInterval=HISTORY_OUTPUT_TIME_INTERVAL



)











# ============================================================



# OPTIONAL JOB



# ============================================================



if CREATE_JOB:







    if JOB_NAME in mdb.jobs.keys():



        del mdb.jobs[JOB_NAME]







    mdb.Job(



        name=JOB_NAME,



        model=MODEL_NAME,



        description='Parametric AlSi10Mg re-entrant quasi-static compression validation',



        type=ANALYSIS,



        explicitPrecision=SINGLE,



        nodalOutputPrecision=SINGLE,



        multiprocessingMode=DEFAULT,



        numCpus=1,



        numDomains=1



    )







    if SUBMIT_JOB:



        mdb.jobs[JOB_NAME].submit(consistencyChecking=OFF)











# ============================================================



# DISPLAY / SUMMARY



# ============================================================



try:



    vp = session.viewports[session.currentViewportName]



    if MERGE_CORE:



        vp.setValues(displayedObject=core_part)



    else:



        vp.setValues(displayedObject=assembly)



    vp.view.fitView()



except:



    pass







print('')



print('============================================================')



print('PARAMETRIC ABAQUS GEOMETRY CREATED')



print('============================================================')



print('L1                     = %.6f mm' % L1)



print('L2                     = %.6f mm' % L2)



print('l_ref                  = %.6f mm' % L_BEVEL)



print('t_core                 = %.6f mm' % CORE_WALL_T)



print('N_ROWS                 = %d' % N_ROWS)



print('N_COLS                 = %d' % N_COLS)



print('CORE_HEIGHT            = %.6f mm' % CORE_HEIGHT)



print('CELL_H                 = %.6f mm' % CELL_H)



print('DEPTH                  = %.6f mm' % DEPTH)



print('lambda_m               = %.6f' % LAMBDA_M)



print('dx_ref                 = %.6f mm' % DX_REF)



print('WaistShift             = %.6f mm' % WAIST_SHIFT)



print('Actual inclined length = %.6f mm' % ACTUAL_INCLINED_LENGTH)



print('Theta                  = %.6f deg' % THETA_DEG)



print('PitchX                 = %.6f mm' % PITCH_X)



print('Throat clearance       = %.6f mm' % THROAT_CLEARANCE)



print('Core X                 = %.6f -> %.6f mm' % (CORE_X_MIN, CORE_X_MAX))



print('Core Y                 = %.6f -> %.6f mm' % (CORE_Y_MIN, CORE_Y_MAX))







if MERGE_CORE:



    print('Final core Part         = Core_AlSi10Mg')



    print('Core cells              = %d' % len(core_part.cells))



    print('Material                = AlSi10Mg')



    print('Elastic E               = %.3f MPa' % ALSI10MG_E)



    print('Poisson ratio           = %.4f' % ALSI10MG_NU)



    print('Density                 = %.6e tonne/mm^3' % ALSI10MG_DENSITY)



    print('Plastic data points     = %d' % len(PLASTIC_DATA))



    print('Ductile damage enabled  = %s' % str(ENABLE_DUCTILE_DAMAGE))



    if ENABLE_DUCTILE_DAMAGE:



        print('Fracture strain         = %.6f' % DUCTILE_FRACTURE_STRAIN)



        print('Stress triaxiality      = %.6f' % DUCTILE_TRIAXIALITY)



        print('Failure displacement    = %.6f mm' % DAMAGE_FAILURE_DISPLACEMENT)



    print('Element formulation     = C3D8R' if ASSIGN_C3D8R else 'Element formulation     = not assigned')



    print('Core mesh generated     = %s' % str(GENERATE_CORE_MESH))
    print('Core global mesh size   = %.6f mm' % CORE_MESH_SIZE)
    if GENERATE_CORE_MESH:
        print('Core element count      = %d' % CORE_ELEMENT_COUNT)
    print('Paper horiz. divisions  = %d' % PAPER_HORIZONTAL_EDGE_ELEMENTS)
    print('Paper bevel divisions   = %d' % PAPER_BEVEL_EDGE_ELEMENTS)
    print('Paper target elements   = %d' % PAPER_TARGET_TOTAL_CORE_ELEMENTS)







if CREATE_STEEL_PLATENS:



    print('Steel platen Parts      = Bottom_Platen_Steel / Top_Platen_Steel')



    print('Steel E                 = %.3f MPa' % STEEL_E)



    print('Steel nu                = %.4f' % STEEL_NU)
    print('Platen width X          = %.6f mm' % PLATEN_WIDTH_X)
    print('Platen depth Z          = %.6f mm' % PLATEN_DEPTH_Z)
    print('Platen thickness Y      = %.6f mm' % PLATEN_THICKNESS)
    print('Reference area          = %.6f mm^2' % REFERENCE_CROSS_SECTION_AREA)



    print('Friction coefficient    = %.4f' % FRICTION_COEFF)



    print('Compression strain      = %.4f' % COMPRESSION_STRAIN)



    print('Compression disp.       = %.6f mm' % COMPRESSION_DISPLACEMENT)



    print('Step time               = %.6f s' % STEP_TIME)



    print('Effective speed         = %.6f mm/s' % EFFECTIVE_COMPRESSION_SPEED)



    print('Paper speed             = %.6f mm/s' % PAPER_COMPRESSION_SPEED_MM_S)



    print('Physical paper speed    = %s' % str(USE_PAPER_PHYSICAL_SPEED))



    print('General Contact         = ON')



    print('Core U3 constraint      = %s' % str(CONSTRAIN_CORE_OUT_OF_PLANE))



    print('Field output            = S/LE/PEEQ/U/V/A/RF/STATUS/SDEG')

    print('History output          = RF2/U2/V2 at top RP')

    print('History time interval   = %.6f s' % HISTORY_OUTPUT_TIME_INTERVAL)



    print('Energy output           = ALLKE/ALLIE/ALLAE/ETOTAL')
    print('Mass scaling enabled    = %s' % str(ENABLE_MASS_SCALING))
    if ENABLE_MASS_SCALING:
        print('Mass scaling target dt  = %.6e s' % MASS_SCALING_TARGET_DT)
        print('Mass scaling frequency  = every %d increment(s)' % MASS_SCALING_FREQUENCY)
        print('Mass scaling mode       = SEMI_AUTOMATIC / MODEL / THROUGHOUT_STEP / BELOW_MIN')



    if not USE_PAPER_PHYSICAL_SPEED:



        print('IMPORTANT               = Check ALLKE/ALLIE for quasi-static validity')







print('============================================================')



print('')







if SAVE_CAE:



    cae_path = os.path.abspath(CAE_FILENAME)



    mdb.saveAs(pathName=cae_path)



    print('CAE saved to: %s' % cae_path)

# Reuse the exact model dimensions when exporting the coupled-RP ODB histories.
# Run export_abaqus_curve.py in Abaqus Python after the job has completed.
if CREATE_STEEL_PLATENS:
    import json
    with open('rp_curve_metadata.json', 'w') as curve_metadata_file:
        json.dump({
            'initial_height_mm': SPECIMEN_HEIGHT,
            'reference_area_mm2': REFERENCE_CROSS_SECTION_AREA,
            'step': STEP_NAME,
            'top_set': 'RP_TOP_PLATEN',
            'bottom_set': 'RP_BOTTOM_PLATEN',
            'strain_convention': '-(U2_top-U2_bottom)/H0',
            'stress_convention': '-RF2_top/A0',
            'coupling': 'kinematic outer platen surface to reference point'
        }, curve_metadata_file, indent=2)
