"""Verification tests for a timber dowel GROUP built out of HankinsonDowel2d
sections -- as opposed to test_hankinson_dowel2d.py, which only ever
exercises a single dowel (one zeroLengthSection at the origin).

HankinsonDowel2d has no notion of a dowel group: it only ever sees one
dowel's 2-D in-plane slip. The group's geometry (how many dowels, their
spacing, and therefore their lever arms about the connection's reference
point) is assembled entirely with ordinary OpenSees node coordinates plus
rigidLink('beam', ...) tying each dowel's free node to a shared master node
-- see EXAMPLES/ExamplePython/HankinsonDowelGroup.py, which this file
imports its model-building helpers from. The M-N-V interaction of the whole
connection then falls out of nodal equilibrium; there is no group-level
bookkeeping inside HankinsonDowel2d/DowelType to test here.

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


# ---------------------------------------------------------------------
# 1. Geometry helper reproduces the Schweigler reference grid
# ---------------------------------------------------------------------

def test_dowel_coordinates_match_reference_grid():
    """ModelFunctions.dowelCoordinatesRP(3, 3, a1, a2, a3, a4), evaluated by
    hand for the Schweigler 3x3 connection, gives a symmetric grid centred
    on the group's geometric centroid: corners at (+-a1, +-a2), edge
    midpoints at (0, +-a2)/(+-a1, 0), one dowel at the centroid."""
    coords = hdg.dowel_coordinates_rp(3, 3, hdg.A1, hdg.A2, hdg.A3, hdg.A4)
    assert len(coords) == 9

    expected = set()
    for cx in (-hdg.A1, 0.0, hdg.A1):
        for cy in (-hdg.A2, 0.0, hdg.A2):
            expected.add((round(cx, 9), round(cy, 9)))

    got = {(round(x, 9), round(y, 9)) for x, y in coords}
    assert got == expected

    # Centroid of the 9 points is exactly the reference point (0, 0).
    xg = sum(x for x, _ in coords) / len(coords)
    yg = sum(y for _, y in coords) / len(coords)
    assert isclose(xg, 0.0, abs_tol=1e-12)
    assert isclose(yg, 0.0, abs_tol=1e-12)


# ---------------------------------------------------------------------
# 2. Initial (small-rotation) rotational stiffness comes from lever arms
#    and each dowel's own Hankinson-blended stiffness -- exactly the
#    "how does geometry enter" question this file exists to answer.
# ---------------------------------------------------------------------

def _hankinson(x0, x90, w):
    return x0 * x90 / (x0 * w + x90 * (1.0 - w))


def _analytic_K_phi(coords, K0_par, K0_perp):
    """Small-rotation rigid-body kinematics: a dowel at radius r sees a pure
    tangential slip of magnitude theta*r, at a load-to-grain angle whose
    sin^2 = x^2/r^2 (slip direction (-y, x)/r is perpendicular to the
    dowel's own position vector). Its secant stiffness there is the
    Hankinson blend of the two calibrations; its moment contribution about
    the centroid is k(alpha)*r^2. Standard rigid bolt-group elastic theory,
    generalized from an isotropic k to Hankinson-anisotropic k(alpha)."""
    total = 0.0
    for x, y in coords:
        r2 = x * x + y * y
        if r2 < 1e-14:
            continue
        w = (x * x) / r2
        k = _hankinson(K0_par, K0_perp, w)
        total += k * r2
    return total


def test_group_initial_rotational_stiffness_matches_lever_arm_theory():
    """A tiny-moment ramp (deep in the pre-yield K0 regime) lets each
    dowel's lagged/committed Hankinson weight settle to its true
    rotation-consistent direction after a couple of steps (the same
    committed-alpha lag documented in HankinsonDowel2d.h and exercised by
    test_hankinson_dowel2d.test_hankinson_closed_form_elastic's ramp). The
    resulting group rotational stiffness must then match the closed-form
    sum of each dowel's own Hankinson-blended K0 times its lever arm
    squared -- i.e. the group's stiffness is exactly the geometry (lever
    arms) combined with each dowel's own anisotropic calibration, with
    nothing else involved.
    """
    coords = hdg.dowel_coordinates_rp(hdg.N_ROWS, hdg.N_COLUMNS,
                                       hdg.A1, hdg.A2, hdg.A3, hdg.A4)
    K_phi_theory = _analytic_K_phi(coords, hdg.K0_0, hdg.K0_90)

    hdg.build_model()
    m_hist, phi_hist = hdg.push_moment(0.01, 60)
    K_fe = (m_hist[-1] - m_hist[-2]) / (phi_hist[-1] - phi_hist[-2])

    assert isclose(K_fe, K_phi_theory, rel_tol=5e-3)


# ---------------------------------------------------------------------
# 3. Global equilibrium: applied (N, V, M) at the master node must be
#    exactly balanced by the grounded dowels' reactions, at every step
#    (including well into the yielded regime) -- this is what "the
#    geometry/lever-arm/moment-assembly comes from the constraint plus
#    equilibrium" (rather than anything inside HankinsonDowel2d) means
#    concretely.
# ---------------------------------------------------------------------

def test_group_pure_moment_equilibrium_through_yield():
    dowels = hdg.build_model()
    m_target, n_steps = 8.0, 80

    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(hdg.MASTER_NODE, 0.0, 0.0, m_target)
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('UmfPack')
    ops.test('NormDispIncr', 1.0e-8, 25)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0 / n_steps)
    ops.analysis('Static')

    for i in range(n_steps):
        ok = ops.analyze(1)
        assert ok == 0, "pushover failed to converge at step %d" % i

        if i % 20 != 19:
            continue  # only check equilibrium periodically (cheap enough,
                       # but no need for all 80 steps)

        ops.reactions()
        Fx = Fy = Msum = 0.0
        for etag, x, y in dowels:
            gnode = 1000 + (etag - 1)
            Rx, Ry, _ = ops.nodeReaction(gnode)
            Fx += Rx
            Fy += Ry
            Msum += x * Ry - y * Rx

        m_applied = m_target * (i + 1) / n_steps
        assert isclose(Fx, 0.0, abs_tol=1e-8)
        assert isclose(Fy, 0.0, abs_tol=1e-8)
        assert isclose(-Msum, m_applied, rel_tol=1e-6, abs_tol=1e-8)


# ---------------------------------------------------------------------
# 4. Geometry actually matters: a wider dowel spacing (bigger lever arms)
#    must give a stiffer and stronger connection in rotation. This is the
#    property that would be completely invisible from
#    test_hankinson_dowel2d.py's single-dowel tests.
# ---------------------------------------------------------------------

def test_wider_spacing_increases_rotational_stiffness():
    def build_and_probe(a1, a2, a3, a4):
        ops.wipe()
        ops.model('basic', '-ndm', 2, '-ndf', 3)
        hdg.saws_dowel_type(hdg.DOWEL_MAT_0, hdg.F0_0, hdg.K0_0, hdg.K1_0)
        hdg.saws_dowel_type(hdg.DOWEL_MAT_90, hdg.F0_90, hdg.K0_90, hdg.K1_90)
        ops.section('HankinsonDowel2d', hdg.SEC_TAG, '-responseHankinson',
                    hdg.DOWEL_MAT_0, hdg.DOWEL_MAT_90)
        coords = hdg.dowel_coordinates_rp(hdg.N_ROWS, hdg.N_COLUMNS, a1, a2, a3, a4)
        dowels = hdg.build_dowel_group(coords, hdg.SEC_TAG, hdg.MASTER_NODE)
        m_hist, phi_hist = hdg.push_moment(0.01, 60)
        return (m_hist[-1] - m_hist[-2]) / (phi_hist[-1] - phi_hist[-2])

    K_reference = build_and_probe(hdg.A1, hdg.A2, hdg.A3, hdg.A4)
    K_wider = build_and_probe(2.0 * hdg.A1, 2.0 * hdg.A2, 2.0 * hdg.A3, 2.0 * hdg.A4)

    # Doubling every spacing quadruples every lever arm squared, so the
    # small-rotation stiffness (~sum of k*r^2) should be ~4x larger.
    assert K_wider > 3.5 * K_reference
