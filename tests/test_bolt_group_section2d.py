"""Verification tests for BoltGroupSection2d -- a whole fastener group as ONE
order-3 {P, Vy, Mz} section, hosted in a single zeroLengthSection element.

The oracle throughout is the hand-built model these tests' companion file
test_hankinson_dowel_group.py exercises: one zeroLengthSection per fastener,
all tied to a shared master node with rigidLink('beam', ...) (see
EXAMPLES/ExamplePython/HankinsonDowelGroup.py). The two are two spellings of
the same mechanics and must agree to round-off, because rigidLink builds an
MP_Constraint with a CONSTANT transformation matrix -- i.e. it is itself
linear in the rotation, exactly like BoltGroupSection2d's default
'-rot linear' kinematics. Anything worse than round-off between them is a
bug in the section, which is why these tests compare absolutely rather than
with a loose engineering tolerance.

Geometry and calibration reproduce the 3x3 dowel connection of Schweigler et
al.'s worked example (schweigler_implementation/MainAnalysis.py).
"""

import os
import sys
import math
from math import isclose

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..',
                                 'EXAMPLES', 'ExamplePython'))
import HankinsonDowelGroup as hdg  # noqa: E402

ops = hdg.ops

COORDS = hdg.dowel_coordinates_rp(hdg.N_ROWS, hdg.N_COLUMNS,
                                   hdg.A1, hdg.A2, hdg.A3, hdg.A4)
NB = len(COORDS)
PROTO_TAG, GROUP_TAG = 10, 1
GND, MASTER = 1, hdg.MASTER_NODE


# ---------------------------------------------------------------------
# model builders
# ---------------------------------------------------------------------

def _prototype():
    """One HankinsonSection2d = one fastener, the group's prototype."""
    hdg.saws_dowel_type(hdg.DOWEL_MAT, hdg.F0_0, hdg.K0_0, hdg.K1_0)
    ops.section('HankinsonSection2d', PROTO_TAG, '-paramHankinson', hdg.DOWEL_MAT,
                '-param', 'F0', hdg.F0_0, hdg.F0_90,
                '-param', 'S0', hdg.K0_0, hdg.K0_90,
                '-param', 'K1', hdg.K1_0, hdg.K1_90)


def build_group(rot='linear', explicit=False, contact=None, ndm=2):
    """The whole connection as ONE element. Returns its tag."""
    ops.wipe()
    if ndm == 2:
        ops.model('basic', '-ndm', 2, '-ndf', 3)
    else:
        ops.model('basic', '-ndm', 3, '-ndf', 6)
    _prototype()

    args = ['BoltGroupSection2d', GROUP_TAG, '-section', PROTO_TAG]
    if explicit:
        args += ['-bolts', NB] + [c for xy in COORDS for c in xy]
    else:
        args += ['-grid', hdg.N_ROWS, hdg.N_COLUMNS, hdg.A1, hdg.A2]
    if contact is not None:
        mat_tag, ys = contact
        args += ['-contact', mat_tag, len(ys)] + list(ys)
    args += ['-rot', rot]
    ops.section(*args)

    if ndm == 2:
        ops.node(GND, 0.0, 0.0)
        ops.fix(GND, 1, 1, 1)
        ops.node(MASTER, 0.0, 0.0)
        ops.element('zeroLengthSection', 1, GND, MASTER, GROUP_TAG,
                    '-orient', 1, 0, 0)
    else:
        ops.node(GND, 0.0, 0.0, 0.0)
        ops.fix(GND, 1, 1, 1, 1, 1, 1)
        ops.node(MASTER, 0.0, 0.0, 0.0)
        ops.element('zeroLengthSection', 1, GND, MASTER, GROUP_TAG,
                    '-orient', 1, 0, 0, 0, 1, 0)
        # the section resists P, Vy, Mz only -- the out-of-plane DOFs get no
        # stiffness from anywhere and must be restrained, exactly as in the
        # hand-built 3D model (HankinsonDowelGroup3d.build_model)
        ops.fix(MASTER, 0, 0, 1, 1, 1, 0)
    return 1


def _analyze(load, n_steps, alg=('Newton',), tol=1.0e-11, maxit=100, ndm=2):
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(MASTER, *load)
    ops.wipeAnalysis()
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('UmfPack')
    ops.test('NormDispIncr', tol, maxit)
    ops.algorithm(*alg)
    ops.integrator('LoadControl', 1.0 / n_steps)
    ops.analysis('Static')
    for i in range(n_steps):
        assert ops.analyze(1) == 0, "failed to converge at step %d" % i
    return ops.nodeDisp(MASTER)


def push_m(m, n_steps, **kw):
    return _analyze((0.0, 0.0, m), n_steps, **kw)


def push_m3d(m, n_steps, **kw):
    return _analyze((0.0, 0.0, 0.0, 0.0, 0.0, m), n_steps, ndm=3, **kw)


def bolt_forces(ele=1):
    return ops.eleResponse(ele, 'section', 'boltForces')


def bolt_slips(ele=1):
    return ops.eleResponse(ele, 'section', 'boltSlips')


# Beware the singular/plural distinction, which is inherited from the
# prototype and not something BoltGroupSection2d invents: HankinsonSection2d
# intercepts 'force' with its own SCALAR slip-law force (a fastener's law is
# one-dimensional along its slip direction), so only the plural 'forces'
# falls through to SectionForceDeformation's 2-component stress resultant.
# 'boltForces' is the concatenation of the plural one.


# ---------------------------------------------------------------------
# 1. Equivalence with the hand-built 9-element + rigidLink model
# ---------------------------------------------------------------------

def test_equivalence_with_multi_element_model_2d():
    """The single-element group and the 9-element + rigidLink model must give
    the same master displacement, the same per-fastener force, and the same
    per-fastener slip, at every load level including well into yield. Both
    are linear in the rotation, so this is an identity, not an approximation:
    the tolerances below are round-off, not engineering tolerances."""
    for n_pre, m_target in [(0.0, 8.0), (40.0, 8.0), (80.0, 10.0)]:
        hdg.build_model()
        if n_pre:
            hdg.apply_normal_force(n_pre)
        ref_u = push_m(m_target, 80)
        ref_f = [ops.eleResponse(i + 1, 'section', 'forces') for i in range(NB)]
        ref_d = [ops.eleResponse(i + 1, 'section', 'deformation') for i in range(NB)]

        build_group()
        if n_pre:
            hdg.apply_normal_force(n_pre)
        new_u = push_m(m_target, 80)
        new_f, new_d = bolt_forces(), bolt_slips()

        f_scale = max(abs(v) for row in ref_f for v in row)
        d_scale = max(abs(v) for row in ref_d for v in row)

        for a, b in zip(ref_u, new_u):
            assert abs(a - b) < 1.0e-12
        for i in range(NB):
            for c in (0, 1):
                assert abs(ref_f[i][c] - new_f[2 * i + c]) < 1.0e-9 * f_scale
                assert abs(ref_d[i][c] - new_d[2 * i + c]) < 1.0e-9 * d_scale


def test_equivalence_with_multi_element_model_3d():
    """Same identity in a 3D model (ndm=3, ndf=6), where the group plane is
    placed by the element's '-orient' rather than by anything in the section."""
    import HankinsonDowelGroup3d as hdg3

    hdg3.build_model()
    hdg3.push_moment(8.0, 60)
    ref_phi = ops.nodeDisp(hdg3.MASTER_NODE, 6)
    ref_f = [ops.eleResponse(i + 1, 'section', 'forces') for i in range(NB)]

    build_group(ndm=3)
    push_m3d(8.0, 60)
    new_phi = ops.nodeDisp(MASTER, 6)
    new_f = bolt_forces()

    f_scale = max(abs(v) for row in ref_f for v in row)
    assert abs(ref_phi - new_phi) < 1.0e-12
    for i in range(NB):
        for c in (0, 1):
            assert abs(ref_f[i][c] - new_f[2 * i + c]) < 1.0e-9 * f_scale


# ---------------------------------------------------------------------
# 2. Geometry input
# ---------------------------------------------------------------------

def test_grid_and_explicit_bolts_agree():
    """'-grid nRows nCols a1 a2' must reproduce dowelCoordinatesRP's grid --
    same values AND same fastener ordering, which is what makes per-fastener
    responses line up with the reference implementation's indices. The two
    are not bit-identical because -grid forms a1*(col-(nc-1)/2) while the
    python helper forms a3+a1*col-centroid, which round differently in the
    last bit."""
    build_group(explicit=False)
    push_m(8.0, 40)
    u_grid, f_grid = ops.nodeDisp(MASTER), bolt_forces()

    build_group(explicit=True)
    push_m(8.0, 40)
    u_expl, f_expl = ops.nodeDisp(MASTER), bolt_forces()

    for a, b in zip(u_grid, u_expl):
        assert abs(a - b) < 1.0e-14
    for a, b in zip(f_grid, f_expl):
        assert abs(a - b) < 1.0e-9


def test_initial_rotational_stiffness_matches_lever_arm_theory():
    """Small-rotation group stiffness must equal sum_i k(alpha_i) * r_i^2 with
    each fastener's own Hankinson-blended k -- the same closed form
    test_hankinson_dowel_group.py checks for the hand-built model."""
    def hankinson(x0, x90, w):
        return x0 * x90 / (x0 * w + x90 * (1.0 - w))

    theory = 0.0
    for x, y in COORDS:
        r2 = x * x + y * y
        if r2 < 1e-14:
            continue
        theory += hankinson(hdg.K0_0, hdg.K0_90, (x * x) / r2) * r2

    build_group()
    m_hist, phi_hist = [], []
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(MASTER, 0.0, 0.0, 0.01)
    ops.wipeAnalysis()
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('UmfPack')
    ops.test('NormDispIncr', 1.0e-12, 50)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0 / 60)
    ops.analysis('Static')
    for i in range(60):
        assert ops.analyze(1) == 0
        m_hist.append(0.01 * (i + 1) / 60)
        phi_hist.append(ops.nodeDisp(MASTER, 3))

    k_fe = (m_hist[-1] - m_hist[-2]) / (phi_hist[-1] - phi_hist[-2])
    assert isclose(k_fe, theory, rel_tol=5e-3)


# ---------------------------------------------------------------------
# 3. Copy semantics
# ---------------------------------------------------------------------

def test_elements_sharing_a_section_tag_stay_independent():
    """Two elements built from the same section tag must not share fastener
    state. Each element gets its own getCopy(), which must deep-copy every
    fastener (and hence every sub-material and its bound Parameters) --
    otherwise both elements would drive one shared set of laws."""
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    _prototype()
    ops.section('BoltGroupSection2d', GROUP_TAG, '-section', PROTO_TAG,
                '-grid', hdg.N_ROWS, hdg.N_COLUMNS, hdg.A1, hdg.A2)
    for k in (0, 1):
        ops.node(100 + k, 0.0, 0.0)
        ops.fix(100 + k, 1, 1, 1)
        ops.node(200 + k, 0.0, 0.0)
        ops.element('zeroLengthSection', 1 + k, 100 + k, 200 + k, GROUP_TAG,
                    '-orient', 1, 0, 0)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(200, 0.0, 0.0, 8.0)        # loaded hard, deep into yield
    ops.load(201, 0.0, 0.0, 0.8)        # loaded 10x less, stays near elastic
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('UmfPack')
    ops.test('NormDispIncr', 1.0e-11, 50)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0 / 60)
    ops.analysis('Static')
    for _ in range(60):
        assert ops.analyze(1) == 0
    together = ops.nodeDisp(201, 3)

    build_group()
    alone = push_m(0.8, 60)[2]
    assert abs(together - alone) < 1.0e-14


# ---------------------------------------------------------------------
# 4. Group-level responses
# ---------------------------------------------------------------------

def test_bolt_response_vectors_match_per_bolt_responses():
    """'boltForces'/'boltSlips' must be exactly the concatenation of the
    per-fastener 'bolt i forces'/'bolt i deformation' responses -- one
    recorder in place of NB of them."""
    build_group()
    push_m(8.0, 60)
    vf, vd = bolt_forces(), bolt_slips()
    assert len(vf) == 2 * NB and len(vd) == 2 * NB
    for i in range(NB):
        f = ops.eleResponse(1, 'section', 'bolt', i, 'forces')
        d = ops.eleResponse(1, 'section', 'bolt', i, 'deformation')
        assert vf[2 * i] == f[0] and vf[2 * i + 1] == f[1]
        assert vd[2 * i] == d[0] and vd[2 * i + 1] == d[1]


def test_bolt_scalar_force_is_the_slip_direction_magnitude():
    """The singular 'bolt i force' reaches HankinsonSection2d's own scalar
    slip-law force rather than the base class's stress resultant. Since a
    fastener's law acts along its slip direction, that scalar's magnitude
    must equal the norm of its (P, Vy) pair -- pinned here because the
    singular/plural distinction is an easy trap when reading these out."""
    build_group()
    push_m(8.0, 60)
    vf = bolt_forces()
    for i in range(NB):
        scalar = ops.eleResponse(1, 'section', 'bolt', i, 'force')
        assert len(scalar) == 1
        assert isclose(abs(scalar[0]),
                       math.hypot(vf[2 * i], vf[2 * i + 1]), abs_tol=1e-12)


def test_center_of_rotation():
    """The centre of rotation is where the rigid-body slip field vanishes,
    (x, y) = (-u_y/theta, u_x/theta) -- the point ModelFunctions.findCR
    locates geometrically. Under a pure moment it sits at the reference
    point; an axial preload pushes it off along y."""
    build_group()
    push_m(8.0, 60)
    cr = ops.eleResponse(1, 'section', 'centerOfRotation')
    assert abs(cr[0]) < 1.0e-9 and abs(cr[1]) < 1.0e-9

    build_group()
    hdg.apply_normal_force(60.0)
    u = push_m(8.0, 60)
    cr = ops.eleResponse(1, 'section', 'centerOfRotation')
    assert isclose(cr[0], -u[1] / u[2], rel_tol=1e-12, abs_tol=1e-12)
    assert isclose(cr[1], u[0] / u[2], rel_tol=1e-12)
    assert abs(cr[1]) > 1.0e-3          # genuinely shifted off the RP


# ---------------------------------------------------------------------
# 5. Against the reference implementation's own formulas
# ---------------------------------------------------------------------

def _slip_curve(disp, alpha):
    """ModelFunctions.slipCurve, verbatim: the Foschi envelope whose F, ki
    and kf are each Hankinson-interpolated by the load-to-grain angle."""
    def hank(x0, x90):
        return (x0 * x90) / (x0 * math.sin(alpha) ** 2 + x90 * math.cos(alpha) ** 2)

    f_a = hank(hdg.F0_0, hdg.F0_90)
    ki_a = hank(hdg.K0_0, hdg.K0_90)
    kf_a = hank(hdg.K1_0, hdg.K1_90)
    return (f_a + kf_a * disp) * (1.0 - math.exp(-(ki_a * disp) / f_a))


def test_pure_translation_matches_reference_slip_curve():
    """Translate the whole group at a fixed 30-degree load-to-grain angle: no
    rotation, so every fastener sees the same slip at the same constant
    alpha, the section's committed-alpha lag is irrelevant after the first
    step, and the group force must be exactly NB times the reference
    implementation's own monotonic slip curve."""
    alpha = math.radians(30.0)
    n_target = 120.0
    build_group()
    u = _analyze((n_target * math.cos(alpha), n_target * math.sin(alpha), 0.0),
                 60, tol=1.0e-12)

    disp = math.hypot(u[0], u[1])
    assert isclose(math.atan2(u[1], u[0]), alpha, rel_tol=1e-9)

    q = ops.eleResponse(1, 'section', 'force')
    f_group = math.hypot(q[0], q[1])
    assert isclose(f_group, NB * _slip_curve(disp, alpha), rel_tol=2e-3)
    assert isclose(f_group, n_target, rel_tol=1e-9)      # equilibrium


# ---------------------------------------------------------------------
# 6. -rot tan
# ---------------------------------------------------------------------

def test_rot_tan_slips_use_the_tangent_map():
    """'-rot tan' must distribute the rotation as tan(theta)*lever, exactly
    as ModelFunctions.dowelDisp does."""
    build_group(rot='tan')
    u = push_m(10.0, 80)
    t = math.tan(u[2])
    slips = bolt_slips()
    for i, (x, y) in enumerate(COORDS):
        assert abs(slips[2 * i] - (u[0] - t * y)) < 1.0e-15
        assert abs(slips[2 * i + 1] - (u[1] + t * x)) < 1.0e-15


def test_rot_tan_moment_is_energy_consistent():
    """Under the tan map the work conjugate of theta is sec^2(theta) times
    the raw lever-arm sum, so the reported M must carry that factor. The
    factor is ~1% at the rotations this pushover reaches, so this genuinely
    discriminates -- note the reference implementation omits it, which is
    why -rot tan reproduces the reference's slips exactly but its moment
    only to O(theta^2)."""
    build_group(rot='tan')
    u = push_m(10.0, 80)
    f = bolt_forces()
    raw = sum(x * f[2 * i + 1] - y * f[2 * i] for i, (x, y) in enumerate(COORDS))
    c = 1.0 + math.tan(u[2]) ** 2

    q = ops.eleResponse(1, 'section', 'force')
    assert isclose(q[2], c * raw, rel_tol=1e-12)
    assert abs(c - 1.0) > 5.0e-3        # the factor is not negligible here


def test_rot_tan_residual_is_independent_of_the_tangent():
    """Initial-stiffness iteration never uses the tangent this section
    assembles, so landing on the same answer as full Newton proves the
    residual (including the sec^2 moment factor) is right on its own, and
    that the geometric term in k(2,2) is a genuine tangent term rather than
    something propping up a wrong force."""
    build_group(rot='tan')
    newton = push_m(10.0, 80)[2]
    build_group(rot='tan')
    initial = push_m(10.0, 80, alg=('ModifiedNewton', '-initial'),
                     tol=1.0e-9, maxit=3000)[2]
    assert isclose(newton, initial, rel_tol=1e-5)


def test_rot_tan_agrees_with_linear_to_second_order():
    """The two rotation modes must differ only at O(theta^2) -- the whole
    justification for shipping 'linear' as the default."""
    build_group(rot='linear')
    lin = push_m(10.0, 80)[2]
    build_group(rot='tan')
    tan = push_m(10.0, 80)[2]

    rel = abs(lin - tan) / abs(lin)
    assert rel < 20.0 * lin ** 2        # O(theta^2), with room for the
    assert rel > 1.0e-4                 # plateau amplification; but not zero


# ---------------------------------------------------------------------
# 7. -contact rows
# ---------------------------------------------------------------------

def test_contact_rows_add_the_closed_form_elastic_stiffness():
    """Two elastic contact rows at +-y act along local x only, so they must
    add exactly 2E to k(0,0) and 2*E*y^2 to k(2,2), and nothing anywhere
    else."""
    e_contact, y_contact = 5.0e4, 0.15
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    _prototype()
    ops.uniaxialMaterial('Elastic', 77, e_contact)
    ops.section('BoltGroupSection2d', GROUP_TAG, '-section', PROTO_TAG,
                '-grid', hdg.N_ROWS, hdg.N_COLUMNS, hdg.A1, hdg.A2,
                '-contact', 77, 2, -y_contact, y_contact)
    ops.node(GND, 0.0, 0.0)
    ops.fix(GND, 1, 1, 1)
    ops.node(MASTER, 0.0, 0.0)
    ops.element('zeroLengthSection', 1, GND, MASTER, GROUP_TAG, '-orient', 1, 0, 0)
    with_contact = ops.eleResponse(1, 'basicStiffness')

    build_group()
    without = ops.eleResponse(1, 'basicStiffness')

    d = [a - b for a, b in zip(with_contact, without)]
    assert isclose(d[0], 2.0 * e_contact, rel_tol=1e-12)          # k(0,0)
    assert isclose(d[8], 2.0 * e_contact * y_contact ** 2, rel_tol=1e-12)  # k(2,2)
    for idx in (1, 2, 3, 4, 5, 6, 7):                             # everything else
        assert abs(d[idx]) < 1.0e-9 * max(abs(with_contact[idx]), 1.0)


def test_contact_rows_carry_load_and_report_forces():
    """A compression-only gap material must stay inert in tension and pick up
    force in compression, and 'contactForces' must report it."""
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    _prototype()
    # ElasticPPGap: closes at -1e-4, then E = 1e6 in compression only
    ops.uniaxialMaterial('ElasticPPGap', 78, 1.0e6, -400.0, -1.0e-4)
    ops.section('BoltGroupSection2d', GROUP_TAG, '-section', PROTO_TAG,
                '-grid', hdg.N_ROWS, hdg.N_COLUMNS, hdg.A1, hdg.A2,
                '-contact', 78, 2, -0.15, 0.15)
    ops.node(GND, 0.0, 0.0)
    ops.fix(GND, 1, 1, 1)
    ops.node(MASTER, 0.0, 0.0)
    ops.element('zeroLengthSection', 1, GND, MASTER, GROUP_TAG, '-orient', 1, 0, 0)

    # pull the group apart: both rows stay open, so they carry nothing
    _analyze((60.0, 0.0, 0.0), 40)
    assert all(abs(v) < 1.0e-12 for v in ops.eleResponse(1, 'section', 'contactForces'))

    # push it together: both rows close and carry compression
    build_group()   # reset
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    _prototype()
    ops.uniaxialMaterial('ElasticPPGap', 78, 1.0e6, -400.0, -1.0e-4)
    ops.section('BoltGroupSection2d', GROUP_TAG, '-section', PROTO_TAG,
                '-grid', hdg.N_ROWS, hdg.N_COLUMNS, hdg.A1, hdg.A2,
                '-contact', 78, 2, -0.15, 0.15)
    ops.node(GND, 0.0, 0.0)
    ops.fix(GND, 1, 1, 1)
    ops.node(MASTER, 0.0, 0.0)
    ops.element('zeroLengthSection', 1, GND, MASTER, GROUP_TAG, '-orient', 1, 0, 0)
    _analyze((-60.0, 0.0, 0.0), 40)
    fc = ops.eleResponse(1, 'section', 'contactForces')
    assert all(v < -1.0e-6 for v in fc)
    assert isclose(fc[0], fc[1], rel_tol=1e-12)     # symmetric rows, pure axial


# ---------------------------------------------------------------------
# 8. Input validation -- must warn and fail, never take the process down
# ---------------------------------------------------------------------

def test_invalid_input_is_rejected_cleanly():
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    _prototype()
    ops.uniaxialMaterial('Elastic', 99, 1000.0)
    ops.section('Aggregator', 21, 99, 'P', 99, 'Vy', 99, 'Mz')   # order 3

    bad = [
        ('order-3 prototype',
         lambda: ops.section('BoltGroupSection2d', 1, '-section', 21,
                             '-grid', 3, 3, 0.07, 0.03)),
        ('missing prototype',
         lambda: ops.section('BoltGroupSection2d', 1, '-section', 999,
                             '-grid', 3, 3, 0.07, 0.03)),
        ('no fastener positions',
         lambda: ops.section('BoltGroupSection2d', 1, '-section', PROTO_TAG)),
        ('short coordinate list',
         lambda: ops.section('BoltGroupSection2d', 1, '-section', PROTO_TAG,
                             '-bolts', 3, 0.0, 0.0)),
        ('unknown flag',
         lambda: ops.section('BoltGroupSection2d', 1, '-section', PROTO_TAG,
                             '-nonsense', 1)),
        ('bad -rot value',
         lambda: ops.section('BoltGroupSection2d', 1, '-section', PROTO_TAG,
                             '-grid', 3, 3, 0.07, 0.03, '-rot', 'quadratic')),
        ('missing contact material',
         lambda: ops.section('BoltGroupSection2d', 1, '-section', PROTO_TAG,
                             '-grid', 3, 3, 0.07, 0.03,
                             '-contact', 998, 1, 0.0)),
    ]
    for why, thunk in bad:
        try:
            thunk()
        except Exception:
            pass
        else:
            raise AssertionError("%s was accepted" % why)

    # ...and the interpreter is still usable afterwards
    assert ops.getNDM()[0] == 2
