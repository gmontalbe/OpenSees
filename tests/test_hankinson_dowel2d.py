"""Verification tests for the HankinsonDowel2d coupled section.

HankinsonDowel2d implements the Schweigler et al. modelling approach for the
M-N-V interaction of timber dowel-type connections: a dowel's 2-D in-plane
slip is Hankinson-interpolated (n = 2) between a 0-degree and a 90-degree
material calibration, using an evolving load-to-grain angle alpha that is a
state variable of the current slip direction (see SRC/material/section/
HankinsonDowel2d.h for the full design notes). Two interpolation modes are
tested: response-level (-responseHankinson, two independent UniaxialMaterial
copies blended) and parameter-level (-paramHankinson, one material whose
named parameters are re-interpolated at every commit -- the literal
Schweigler formulation, exercised here via DowelType's Foschi envelope,
which was extended with setParameter/updateParameter for exactly this
purpose).

All rigs use a 2D model (ndm=2, ndf=3) with a zeroLengthSection element:
node 1 fully fixed, node 2 free in-plane (dof 3, the unused rotation, is
fixed since HankinsonDowel2d's default response codes are {P, Vy}).

Scope note: several tests intentionally stay in the monotonic / mildly
cyclic regime. Driving a post-yield hysteretic material through a reversal
under load control is a well known Newton-robustness edge case in OpenSees
generally (confirmed here to fail identically for a plain zeroLength+Steel01
rig with no HankinsonDowel2d involved at all) -- exercising that robustness
is not what these tests are for. The regularization used to keep the Mode B
Hankinson blend finite when the two calibrations' forces have opposite sign
is instead checked directly against the closed-form formula in pure Python
(no FE model needed).
"""

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

import math
from math import isclose, radians, sin, cos, hypot


# ---------------------------------------------------------------------
# Shared rig helpers
# ---------------------------------------------------------------------

def _fresh_2node():
    """Fully-fixed node 1, node 2 free in-plane (dof3/rotation fixed)."""
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    ops.node(1, 0, 0)
    ops.fix(1, 1, 1, 1)
    ops.node(2, 0, 0)
    ops.fix(2, 0, 0, 1)


def _ramp_load(ang_deg, n_steps, fx=1.0):
    """Ramp a load at a fixed angle over n_steps LoadControl increments."""
    ang = radians(ang_deg)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, fx * cos(ang), fx * sin(ang), 0.0)
    ops.analysis('Static', '-noWarnings')
    for i in range(n_steps):
        ok = ops.analyze(1)
        assert ok == 0, "analyze failed at step %d (angle %s)" % (i, ang_deg)


# ---------------------------------------------------------------------
# 1. Closed-form Hankinson interpolation (Mode A, elastic material)
# ---------------------------------------------------------------------

def test_hankinson_closed_form_elastic():
    """Ramping a fixed-direction load lets the lagged alpha converge to the
    true load angle; the resulting secant stiffness must then match the
    closed-form Hankinson interpolation of the two calibration stiffnesses
    to essentially machine precision (the material is linear, so once alpha
    has converged there is nothing left to approximate)."""
    E0, E90 = 1000.0, 500.0
    for ang in (0, 15, 30, 45, 60, 75, 90):
        _fresh_2node()
        ops.uniaxialMaterial('Elastic', 1, E0)
        ops.section('HankinsonDowel2d', 1, '-paramHankinson', 1,
                     '-param', 'E', E0, E90)
        ops.element('zeroLengthSection', 1, 1, 2, 1, '-orient', 1, 0, 0)
        _ramp_load(ang, 150)

        d = ops.nodeDisp(2)
        f = ops.eleResponse(1, 'basicForce')
        delta = hypot(d[0], d[1])
        fmag = hypot(f[0], f[1])

        w = sin(radians(ang)) ** 2
        c2 = cos(radians(ang)) ** 2
        E_expect = E0 * E90 / (E0 * w + E90 * c2)
        assert isclose(fmag / delta, E_expect, rel_tol=1e-8), ang


# ---------------------------------------------------------------------
# 2. Consistent tangent (Mode B): analytic kb must equal a true FD probe
# ---------------------------------------------------------------------

def test_tangent_matches_finite_difference():
    """The section's reported 2x2 tangent must be the true derivative of its
    force w.r.t. deformation. Ramp a hysteretic (post-yield) Steel01 pair to
    a converged, committed state, freeze that load with loadConst, read the
    analytic tangent, then apply one tiny additional load increment and
    verify K_analytic @ delta_e reproduces the applied delta_P. This is
    exactly the check that would catch CoupledZeroLength's known defect: an
    incomplete tangent (missing the radial g*(I-n(x)n) term, or the wrong
    sign of g on a reversed branch)."""
    _fresh_2node()
    ops.uniaxialMaterial('Steel01', 1, 15.0, 24460.0, 0.02)
    ops.uniaxialMaterial('Steel01', 2, 10.0, 19640.0, 0.05)
    ops.section('HankinsonDowel2d', 1, '-responseHankinson', 1, 2)
    ops.element('zeroLengthSection', 1, 1, 2, 1, '-orient', 1, 0, 0)

    ang = 35.0
    _ramp_load(ang, 60, fx=20.0)
    ops.loadConst('-time', 0.0)

    d0 = ops.nodeDisp(2)
    K = ops.eleResponse(1, 'basicStiffness')  # [Kxx, Kxy, Kyx, Kyy]

    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    dPx, dPy = 1.0e-4, -0.5e-4
    ops.load(2, dPx, dPy, 0.0)
    ok = ops.analyze(1)
    assert ok == 0

    d1 = ops.nodeDisp(2)
    de = (d1[0] - d0[0], d1[1] - d0[1])

    dPx_pred = K[0] * de[0] + K[1] * de[1]
    dPy_pred = K[2] * de[0] + K[3] * de[1]

    assert isclose(dPx_pred, dPx, rel_tol=1e-3, abs_tol=1e-9)
    assert isclose(dPy_pred, dPy, rel_tol=1e-3, abs_tol=1e-9)
    # tangent must be symmetric (committed-alpha scheme => symmetric kb)
    assert isclose(K[1], K[2], rel_tol=1e-9, abs_tol=1e-12)


def test_radial_term_present():
    """Regression guard against silently reverting to a diagonal (radial-
    term-free) tangent, as CoupledZeroLength ships: at a yielded state the
    transverse stiffness must equal the secant modulus F/s, and must not be
    zero."""
    _fresh_2node()
    ops.uniaxialMaterial('Steel01', 1, 15.0, 24460.0, 0.02)
    ops.uniaxialMaterial('Steel01', 2, 15.0, 24460.0, 0.02)
    ops.section('HankinsonDowel2d', 1, '-responseHankinson', 1, 2)
    ops.element('zeroLengthSection', 1, 1, 2, 1, '-orient', 1, 0, 0)

    # Pure axial (alpha = 0 throughout -- no ramp needed to "discover" it).
    _ramp_load(0.0, 40, fx=25.0)

    F = ops.eleResponse(1, 'basicForce')
    d = ops.nodeDisp(2)
    K = ops.eleResponse(1, 'basicStiffness')

    s = d[0]
    g_expect = F[0] / s
    assert abs(K[3]) > 1.0e-6
    assert isclose(K[3], g_expect, rel_tol=1e-9)


# ---------------------------------------------------------------------
# 3. Sign rule through a reversal (elastic regime, so Newton is trivial)
# ---------------------------------------------------------------------

def test_sign_rule_elastic_reversal():
    """The rotating sign rule must reduce to an ordinary 1-D spring: through
    a reversal the reported force must stay an odd function of slip, exactly
    tracking E0*disp (elastic -- no hysteresis to obscure a sign bug), and
    the 'sign' response must flip exactly at each zero-force crossing."""
    _fresh_2node()
    E0 = 24460.0
    ops.uniaxialMaterial('Elastic', 1, E0)
    ops.uniaxialMaterial('Elastic', 2, E0)
    ops.section('HankinsonDowel2d', 1, '-responseHankinson', 1, 2)
    ops.element('zeroLengthSection', 1, 1, 2, 1, '-orient', 1, 0, 0)

    def frange(a, b, n):
        return [a + (b - a) * i / n for i in range(0, n + 1)]

    factors = [0.0] + frange(0, 10, 20) + frange(10, -10, 40)[1:] + frange(-10, 10, 40)[1:]
    ops.timeSeries('Path', 1, '-dt', 1.0, '-values', *factors)
    ops.pattern('Plain', 1, 1)
    ops.load(2, 1.0, 0.0, 0.0)
    ops.analysis('Static', '-noWarnings')

    signs = []
    for i in range(len(factors) - 1):
        ok = ops.analyze(1)
        assert ok == 0, i
        d = ops.nodeDisp(2)
        f = ops.eleResponse(1, 'basicForce')
        assert isclose(f[0], E0 * d[0], rel_tol=1e-9, abs_tol=1e-9)
        signs.append(ops.eleResponse(1, 'section', 'sign')[0])

    flips = sum(1 for a, b in zip(signs, signs[1:]) if a != b)
    assert flips == 2  # one flip per zero-crossing (down-up-down => 2)


# ---------------------------------------------------------------------
# 4. Alpha converges to the true load angle (Mode A, hysteretic DowelType)
# ---------------------------------------------------------------------

def test_alpha_converges_to_load_angle():
    _fresh_2node()
    ops.uniaxialMaterial('DowelType', 10, 5.0, 200.0, 0.8, 0.2, 0.6, 0.9,
                          0.5, 2.0, 0.8, 0.9, 0.9,
                          '-exponential', 24460.0, 462.0 / 24460.0, 19.1, 0.05, 5000.0)
    ops.section('HankinsonDowel2d', 1, '-paramHankinson', 10,
                '-param', 'F0p', 19.1, 14.18,
                '-param', 'K0p', 24460.0, 19640.0,
                '-param', 'K1p', 462.0, 739.0)
    ops.element('zeroLengthSection', 1, 1, 2, 1, '-orient', 1, 0, 0)

    _ramp_load(45.0, 25, fx=1.0)

    alpha_deg = ops.eleResponse(1, 'section', 'alphaDeg')[0]
    assert isclose(alpha_deg, 45.0, abs_tol=0.05)


# ---------------------------------------------------------------------
# 5. Mode A reproduces the exact Schweigler Foschi + Hankinson formula
# ---------------------------------------------------------------------

def _foschi(disp, f0, kf, ki):
    return (f0 + kf * disp) * (1.0 - math.exp(-ki * disp / f0))


def _hankinson(x0, x90, w):
    c2 = 1.0 - w
    return x0 * x90 / (x0 * w + x90 * c2)


def test_foschi_hankinson_matches_closed_form():
    """Mode A + DowelType's -exponential envelope is exactly the Python
    reference's slipCurve (F0/K0/K1 <-> F0/ki/kf). Reproduce the reference
    formula's own Hankinson interpolation at several angles and check the
    element's force matches the closed-form Foschi-at-X_alpha value."""
    F0_0, KF_0, KI_0 = 19.1, 462.0, 24460.0
    F0_9, KF_9, KI_9 = 14.18, 739.0, 19640.0

    for ang in (0, 15, 30, 45, 60, 75, 90):
        _fresh_2node()
        ops.uniaxialMaterial('DowelType', 10, 5.0, 200.0, 0.8, 0.2, 0.6, 0.9,
                              0.5, 2.0, 0.8, 0.9, 0.9,
                              '-exponential', KI_0, KF_0 / KI_0, F0_0, 0.05, 5000.0)
        ops.section('HankinsonDowel2d', 1, '-paramHankinson', 10,
                    '-param', 'F0p', F0_0, F0_9,
                    '-param', 'K0p', KI_0, KI_9,
                    '-param', 'K1p', KF_0, KF_9)
        ops.element('zeroLengthSection', 1, 1, 2, 1, '-orient', 1, 0, 0)

        # Stay comfortably below the envelope's cap displacement (~0.05 m)
        # so the response is monotonic-hardening throughout (load control
        # cannot trace the post-cap softening branch -- an FE technique
        # limitation, not something this test needs to exercise).
        _ramp_load(ang, 25, fx=1.0)

        d = ops.nodeDisp(2)
        f = ops.eleResponse(1, 'basicForce')
        delta = hypot(d[0], d[1])
        fmag = hypot(f[0], f[1])

        w = sin(radians(ang)) ** 2
        f0 = _hankinson(F0_0, F0_9, w)
        kf = _hankinson(KF_0, KF_9, w)
        ki = _hankinson(KI_0, KI_9, w)
        f_expect = _foschi(delta, f0, kf, ki)

        assert isclose(fmag, f_expect, rel_tol=1e-8), ang


# ---------------------------------------------------------------------
# 6. Calibration angles are reproduced exactly by the Mode B blend
# ---------------------------------------------------------------------

def test_calibration_angles_exact():
    """At alpha = 0 or alpha = 90 the regularized Hankinson blend must
    collapse to exactly mat0 or mat90's own response (the blend's defining
    property, per the closed-form check: at w=0, H = F0*F90/F90 = F0 and
    L = F0, for ANY confidence beta)."""
    for pure_angle, tag_expected in ((0.0, 1), (90.0, 2)):
        _fresh_2node()
        ops.uniaxialMaterial('Steel01', 1, 15.0, 24460.0, 0.02)
        ops.uniaxialMaterial('Steel01', 2, 10.0, 19640.0, 0.05)
        ops.section('HankinsonDowel2d', 1, '-responseHankinson', 1, 2)
        ops.element('zeroLengthSection', 1, 1, 2, 1, '-orient', 1, 0, 0)
        _ramp_load(pure_angle, 40, fx=25.0)

        d = ops.nodeDisp(2)
        f = ops.eleResponse(1, 'basicForce')
        delta = hypot(d[0], d[1])
        fmag = hypot(f[0], f[1])

        # Independent reference: drive a bare copy of the expected material.
        ops.wipe()
        Fy, E, b = ((15.0, 24460.0, 0.02) if tag_expected == 1
                    else (10.0, 19640.0, 0.05))
        ops.uniaxialMaterial('Steel01', 99, Fy, E, b)
        ops.testUniaxialMaterial(99)
        ops.setStrain(delta)
        f_ref = ops.getStress()

        assert isclose(fmag, f_ref, rel_tol=1e-9, abs_tol=1e-9)


# ---------------------------------------------------------------------
# 7. Mode B regularization: bounded and continuous, verified against the
#    closed-form blend formula directly (no FE model needed)
# ---------------------------------------------------------------------

def _smoothstep(x):
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    return x * x * (3.0 - 2.0 * x)


def _mode_b_blend(F0, F90, w, ftol):
    """Pure-Python mirror of HankinsonDowel2d::evaluateLaw's Mode B force
    blend (see HankinsonDowel2d.h/.cpp), used to check the regularization's
    boundedness/continuity across the full (F0, F90, w) input space --
    including the opposite-sign case that only arises transiently on a
    hysteretic unloading branch and is impractical to reach reliably via a
    load-controlled FE reversal (a Newton-robustness question, not a
    property of this formula)."""
    c2 = 1.0 - w
    L = c2 * F0 + w * F90
    D = c2 * F90 + w * F0
    u0, u90 = F0 / ftol, F90 / ftol
    beta = (_smoothstep(u0) * _smoothstep(u90)
            + _smoothstep(-u0) * _smoothstep(-u90))
    if beta <= 0.0 or D == 0.0:
        return L
    H = F0 * F90 / D
    return beta * H + (1.0 - beta) * L


def test_modeB_regularization_bounded_and_exact_at_calibration():
    ftol = 0.1
    import random
    random.seed(0)
    for _ in range(2000):
        F0 = random.uniform(-20.0, 20.0)
        F90 = random.uniform(-20.0, 20.0)
        w = random.uniform(0.0, 1.0)
        F = _mode_b_blend(F0, F90, w, ftol)
        assert math.isfinite(F)
        assert abs(F) <= max(abs(F0), abs(F90)) + 1.0e-9

    # Exact at the calibration angles, for a range of ftol and even when
    # F0/F90 disagree in sign (only F0 -- or only F90 -- can matter at w=0/1).
    for ftol in (1.0e-6, 1.0e-2, 1.0, 100.0):
        for F0, F90 in ((5.0, -3.0), (-5.0, 3.0), (0.0, -3.0), (5.0, 0.0)):
            assert isclose(_mode_b_blend(F0, F90, 0.0, ftol), F0, abs_tol=1e-12)
            assert isclose(_mode_b_blend(F0, F90, 1.0, ftol), F90, abs_tol=1e-12)


# ---------------------------------------------------------------------
# 8. getCopy() independence
# ---------------------------------------------------------------------

def test_getCopy_independence():
    """ZeroLengthSection copies the prototype section per element; loading
    one element into yield must not affect a second, unloaded element that
    shares the same section tag (guards the Parameter-rebinding trap noted
    in HankinsonDowel2d.h Sec. getCopy())."""
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    ops.uniaxialMaterial('Steel01', 1, 15.0, 24460.0, 0.02)
    ops.uniaxialMaterial('Steel01', 2, 10.0, 19640.0, 0.05)
    ops.section('HankinsonDowel2d', 1, '-responseHankinson', 1, 2)

    for i in (1, 2):
        ops.node(10 * i, 0, 0)
        ops.fix(10 * i, 1, 1, 1)
        ops.node(10 * i + 1, 0, 0)
        ops.fix(10 * i + 1, 0, 0, 1)
        ops.element('zeroLengthSection', i, 10 * i, 10 * i + 1, 1, '-orient', 1, 0, 0)

    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(11, 20.0, 0.0, 0.0)  # only element 1 is loaded
    ops.analysis('Static', '-noWarnings')
    ok = ops.analyze(1)
    assert ok == 0

    d1 = ops.nodeDisp(11)
    d2 = ops.nodeDisp(21)
    assert abs(d1[0]) > 1.0e-4          # element 1 yielded
    assert d2[0] == 0.0                 # element 2 untouched


# ---------------------------------------------------------------------
# 9. DowelType setParameter/updateParameter
# ---------------------------------------------------------------------

def test_dowel_type_setparameter():
    """Mutating F0p via the generic Parameter mechanism must change the
    envelope force at a fixed probe displacement to exactly the analytic
    Foschi value with the new F0 (and the recomputed derived quantities)."""
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    ops.node(1, 0, 0)
    ops.fix(1, 1, 1, 1)
    ops.node(2, 0, 0)
    ops.fix(2, 0, 1, 1)

    K0, R1, F0, Dc, Kd = 24460.0, 462.0 / 24460.0, 19.1, 0.05, 5000.0
    ops.uniaxialMaterial('DowelType', 10, 5.0, 200.0, 0.8, 0.2, 0.6, 0.9,
                          0.5, 2.0, 0.8, 0.9, 0.9,
                          '-exponential', K0, R1, F0, Dc, Kd)
    ops.element('zeroLength', 1, 1, 2, '-mat', 10, '-dir', 1)

    ops.parameter(1, 'element', 1, 'F0p')
    F0_new = 30.0
    ops.updateParameter(1, F0_new)

    probe = 0.001
    K1 = R1 * K0
    f_expect = (F0_new + K1 * probe) * (1.0 - math.exp(-K0 * probe / F0_new))

    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, f_expect, 0.0, 0.0)
    ops.analysis('Static', '-noWarnings')
    ok = ops.analyze(1)
    assert ok == 0

    d = ops.nodeDisp(2)
    assert isclose(d[0], probe, rel_tol=2e-2)

    # Recompute the analytic Foschi force at the actual converged
    # displacement (not the target probe) for an exact, solver-tolerance
    # independent check.
    f_at_converged = (F0_new + K1 * d[0]) * (1.0 - math.exp(-K0 * d[0] / F0_new))
    f = ops.eleResponse(1, 'force')
    assert isclose(f[3], f_at_converged, rel_tol=1e-9)
    assert isclose(f[0], -f_at_converged, rel_tol=1e-9)
