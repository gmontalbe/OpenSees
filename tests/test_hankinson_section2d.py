"""Verification tests for the HankinsonSection2d coupled section.

HankinsonSection2d implements the Schweigler et al. modelling approach for the
M-N-V interaction of timber dowel-type connections: a dowel's 2-D in-plane
slip is Hankinson-interpolated (n = 2) between a 0-degree and a 90-degree
material calibration, using an evolving load-to-grain angle alpha that is a
state variable of the current slip direction (see SRC/material/section/
HankinsonSection2d.h for the full design notes).

This is parameter-level interpolation only (-paramHankinson): one material
whose named parameters are re-interpolated between 0-deg and 90-deg values
and pushed into it (via setParameter/updateParameter) at every commit -- the
literal Schweigler formulation. Exercised here via DowelType's Foschi
envelope, SAWSMaterial's (identical) Foschi envelope, Steel01 and
ElasticMaterial -- all extended with setParameter/updateParameter for
exactly this purpose. An earlier response-level blending mode
(-responseHankinson, two independent materials whose scalar force responses
were blended) was removed: it had no clear physical basis (see
HankinsonSection2d.h for the reasoning) and is not part of this class anymore.

All rigs use a 2D model (ndm=2, ndf=3) with a zeroLengthSection element:
node 1 fully fixed, node 2 free in-plane (dof 3, the unused rotation, is
fixed since HankinsonSection2d's default response codes are {P, Vy}).

Scope note: several tests intentionally stay in the monotonic / mildly
cyclic regime. Driving a post-yield hysteretic material through a reversal
under load control is a well known Newton-robustness edge case in OpenSees
generally (confirmed here to fail identically for a plain zeroLength+Steel01
rig with no HankinsonSection2d involved at all) -- exercising that robustness
is not what these tests are for.
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


def _dowel_type_exponential(tag):
    """A DowelType material with the Schweigler 0-deg Foschi envelope, used
    as the Mode A driven material in several rigs below (its F0p/K0p/K1p
    parameters are what HankinsonSection2d re-interpolates)."""
    ops.uniaxialMaterial('DowelType', tag, 5.0, 200.0, 0.8, 0.2, 0.6, 0.9,
                          0.5, 2.0, 0.8, 0.9, 0.9,
                          '-exponential', 24460.0, 462.0 / 24460.0, 19.1, 0.05, 5000.0)


# ---------------------------------------------------------------------
# 1. Closed-form Hankinson interpolation (elastic material)
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
        ops.section('HankinsonSection2d', 1, '-paramHankinson', 1,
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
# 2. Consistent tangent: analytic kb must equal a true FD probe
# ---------------------------------------------------------------------

def test_tangent_matches_finite_difference():
    """The section's reported 2x2 tangent must be the true derivative of its
    force w.r.t. deformation. Ramp a hysteretic (post-yield) DowelType rig to
    a converged, committed state, freeze that load with loadConst, read the
    analytic tangent, then apply one tiny additional load increment and
    verify K_analytic @ delta_e reproduces the applied delta_P. This is
    exactly the check that would catch an incomplete tangent (missing the
    radial g*(I-n(x)n) term, or the wrong sign of g on a reversed branch) --
    the known defect CoupledZeroLength ships with."""
    _fresh_2node()
    _dowel_type_exponential(10)
    ops.section('HankinsonSection2d', 1, '-paramHankinson', 10,
                '-param', 'F0p', 19.1, 14.18,
                '-param', 'K0p', 24460.0, 19640.0,
                '-param', 'K1p', 462.0, 739.0)
    ops.element('zeroLengthSection', 1, 1, 2, 1, '-orient', 1, 0, 0)

    # _ramp_load has no explicit fractional LoadControl increment, so each
    # analyze(1) call re-applies the full pattern (cumulative load after
    # n_steps is fx*n_steps, not fx) -- fx here is chosen so that product
    # lands comfortably between yield onset (~16.5 kN blended at 35 deg) and
    # the envelope cap (~35-40 kN blended), not literally "20 kN".
    ang = 35.0
    _ramp_load(ang, 60, fx=0.5)
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
    _dowel_type_exponential(10)
    ops.section('HankinsonSection2d', 1, '-paramHankinson', 10,
                '-param', 'F0p', 19.1, 14.18,
                '-param', 'K0p', 24460.0, 19640.0,
                '-param', 'K1p', 462.0, 739.0)
    ops.element('zeroLengthSection', 1, 1, 2, 1, '-orient', 1, 0, 0)

    # Pure axial (alpha = 0 throughout -- no ramp needed to "discover" it),
    # well past F0 = 19.1 so the material is on the yielded branch. (See the
    # comment in test_tangent_matches_finite_difference on why fx*n_steps,
    # not fx, is the actual cumulative applied force here.)
    _ramp_load(0.0, 40, fx=0.6)

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
    ops.section('HankinsonSection2d', 1, '-paramHankinson', 1,
                '-param', 'E', E0, E0)
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
# 4. Alpha converges to the true load angle (hysteretic DowelType)
# ---------------------------------------------------------------------

def test_alpha_converges_to_load_angle():
    _fresh_2node()
    _dowel_type_exponential(10)
    ops.section('HankinsonSection2d', 1, '-paramHankinson', 10,
                '-param', 'F0p', 19.1, 14.18,
                '-param', 'K0p', 24460.0, 19640.0,
                '-param', 'K1p', 462.0, 739.0)
    ops.element('zeroLengthSection', 1, 1, 2, 1, '-orient', 1, 0, 0)

    _ramp_load(45.0, 25, fx=1.0)

    alpha_deg = ops.eleResponse(1, 'section', 'alphaDeg')[0]
    assert isclose(alpha_deg, 45.0, abs_tol=0.05)


# ---------------------------------------------------------------------
# 5. Reproduces the exact Schweigler Foschi + Hankinson formula
# ---------------------------------------------------------------------

def _foschi(disp, f0, kf, ki):
    return (f0 + kf * disp) * (1.0 - math.exp(-ki * disp / f0))


def _hankinson(x0, x90, w):
    c2 = 1.0 - w
    return x0 * x90 / (x0 * w + x90 * c2)


def test_foschi_hankinson_matches_closed_form():
    """DowelType's -exponential envelope is exactly the Python reference's
    slipCurve (F0/K0/K1 <-> F0/ki/kf). Reproduce the reference formula's own
    Hankinson interpolation at several angles and check the element's force
    matches the closed-form Foschi-at-X_alpha value."""
    F0_0, KF_0, KI_0 = 19.1, 462.0, 24460.0
    F0_9, KF_9, KI_9 = 14.18, 739.0, 19640.0

    for ang in (0, 15, 30, 45, 60, 75, 90):
        _fresh_2node()
        ops.uniaxialMaterial('DowelType', 10, 5.0, 200.0, 0.8, 0.2, 0.6, 0.9,
                              0.5, 2.0, 0.8, 0.9, 0.9,
                              '-exponential', KI_0, KF_0 / KI_0, F0_0, 0.05, 5000.0)
        ops.section('HankinsonSection2d', 1, '-paramHankinson', 10,
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


def test_saws_hankinson_matches_closed_form():
    """SAWSMaterial's monotonic envelope (LPATH==1, first loading) is the
    identical Foschi curve, F0/K1<->F0/kf and S0<->ki -- confirmed directly
    against SAWSMaterial.cpp's setTrialStrain(). Same check as
    test_foschi_hankinson_matches_closed_form but driving SAWSMaterial's own
    setParameter/updateParameter (F0/S0/K1) instead of DowelType's, i.e. the
    worked example's actual material (EXAMPLES/ExamplePython/
    HankinsonDowelGroup.py)."""
    F0_0, KF_0, KI_0 = 19.1, 462.0, 24460.0
    F0_9, KF_9, KI_9 = 14.18, 739.0, 19640.0

    for ang in (0, 15, 30, 45, 60, 75, 90):
        _fresh_2node()
        # DU large enough that the failure displacement DF never bites at
        # the small probe displacements this test stays within.
        ops.uniaxialMaterial('SAWS', 10, F0_0, 0.1 * F0_0, 0.05, KI_0,
                              KF_0 / KI_0, -0.05, 1.0, 0.0, 0.8, 1.05)
        ops.section('HankinsonSection2d', 1, '-paramHankinson', 10,
                    '-param', 'F0', F0_0, F0_9,
                    '-param', 'S0', KI_0, KI_9,
                    '-param', 'K1', KF_0, KF_9)
        ops.element('zeroLengthSection', 1, 1, 2, 1, '-orient', 1, 0, 0)

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
# 6. getCopy() independence
# ---------------------------------------------------------------------

def test_getCopy_independence():
    """ZeroLengthSection copies the prototype section per element; loading
    one element into yield must not affect a second, unloaded element that
    shares the same section tag. This is exactly the Parameter-rebinding
    trap noted in HankinsonSection2d.h/getCopy(): a memberwise copy would
    leave every element's bound Parameters pointing at the PROTOTYPE
    material, so all elements sharing this section tag would silently
    mutate one shared law."""
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    _dowel_type_exponential(10)
    ops.section('HankinsonSection2d', 1, '-paramHankinson', 10,
                '-param', 'F0p', 19.1, 14.18,
                '-param', 'K0p', 24460.0, 19640.0,
                '-param', 'K1p', 462.0, 739.0)

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
# 7. DowelType / SAWSMaterial setParameter/updateParameter
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


def test_saws_setparameter():
    """Same check as test_dowel_type_setparameter, but for SAWSMaterial's
    F0/S0/K1 (mutating K1 directly, not R1 -- see SAWSMaterial::setParameter
    for why): the monotonic envelope force at a fixed probe displacement
    must match the analytic Foschi value with the new F0."""
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)
    ops.node(1, 0, 0)
    ops.fix(1, 1, 1, 1)
    ops.node(2, 0, 0)
    ops.fix(2, 0, 1, 1)

    S0, K1, F0 = 24460.0, 462.0, 19.1
    ops.uniaxialMaterial('SAWS', 10, F0, 0.1 * F0, 0.05, S0,
                          K1 / S0, -0.05, 1.0, 0.0, 0.8, 1.05)
    ops.element('zeroLength', 1, 1, 2, '-mat', 10, '-dir', 1)

    ops.parameter(1, 'element', 1, 'F0')
    F0_new = 30.0
    ops.updateParameter(1, F0_new)

    probe = 0.001
    f_expect = (F0_new + K1 * probe) * (1.0 - math.exp(-S0 * probe / F0_new))

    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, f_expect, 0.0, 0.0)
    ops.analysis('Static', '-noWarnings')
    ok = ops.analyze(1)
    assert ok == 0

    d = ops.nodeDisp(2)
    assert isclose(d[0], probe, rel_tol=2e-2)

    f_at_converged = (F0_new + K1 * d[0]) * (1.0 - math.exp(-S0 * d[0] / F0_new))
    f = ops.eleResponse(1, 'force')
    assert isclose(f[3], f_at_converged, rel_tol=1e-9)
    assert isclose(f[0], -f_at_converged, rel_tol=1e-9)
