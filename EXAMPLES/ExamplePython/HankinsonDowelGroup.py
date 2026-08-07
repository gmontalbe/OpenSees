# Worked example: assembling a timber dowel GROUP (not a single dowel) out of
# HankinsonDowel2d sections, following the "Building a dowel group" recipe --
# one zeroLengthSection per physical dowel location, all tied to a shared
# rigid master node with rigidLink('beam', ...).
#
# HankinsonDowel2d itself only ever sees a single dowel's 2-D in-plane slip;
# it has no notion of a "group". The group geometry (how many dowels, their
# spacing, and therefore their lever arms about the connection's reference
# point) lives entirely in ordinary node coordinates plus the rigidLink
# constraint -- the M-N-V interaction of the whole connection then falls out
# of nodal equilibrium (sum of F_i x r_i) with no additional bookkeeping.
#
# Each dowel's 0-deg/90-deg calibration is Hankinson-blended by
# HankinsonDowel2d in RESPONSE mode (-responseHankinson mat0 mat90): two
# SAWSMaterial instances (the Folz & Filiatrault CUREE pinching hysteretic
# model), one calibrated parallel to grain and one perpendicular, whose
# response is blended per dowel according to its own current slip direction.
# SAWSMaterial has no setParameter/updateParameter (checked directly against
# SRC/material/uniaxial/SAWSMaterial.cpp), so it cannot be driven through
# HankinsonDowel2d's PARAMETER mode (-paramHankinson, used with e.g.
# DowelType in tests/test_hankinson_dowel2d.py) -- response-level blending is
# the only option, which is exactly the "works with any UniaxialMaterial"
# case Mode B exists for.
#
# The default group below reproduces the geometry of Schweigler et al.'s
# worked example (a 3x3 dowel connection): see /schweigler_implementation/
# MainAnalysis.py, function dowelCoordinatesRP in ModelFunctions.py. Units
# are consistent with the calibration values (force in kN, length in m) --
# the same scale used throughout tests/test_hankinson_dowel2d.py.
#
# Run standalone:
#   cd EXAMPLES/ExamplePython && cp ../../build/Release/OpenSeesPy.so ./opensees.so
#   python HankinsonDowelGroup.py

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


def dowel_coordinates_rp(n_rows, n_columns, a1, a2, a3, a4):
    """Grid of dowel (x, y) coordinates referred to the group's geometric
    centre (the connection's reference point RP). x runs parallel to grain
    at spacing a1 (column pitch), y runs perpendicular to grain at spacing
    a2 (row pitch); a3/a4 are the edge distances of the first row/column.
    Mirrors ModelFunctions.dowelCoordinatesRP in the Python reference."""
    coords = []
    for row in range(n_rows):
        for col in range(n_columns):
            coords.append((a3 + a1 * col, a4 + a2 * row))
    x_cg = a3 + (0.5 * n_columns - 0.5) * a1
    y_cg = a4 + (0.5 * n_rows - 0.5) * a2
    return [(x - x_cg, y - y_cg) for x, y in coords]


def build_dowel_group(dowel_coords, sec_tag, master_node,
                       ground_node_start=1000, free_node_start=2000,
                       ele_tag_start=1):
    """Place one zeroLengthSection element (hosting section `sec_tag`) at
    each (x, y) in dowel_coords, and rigid-link its free node to
    `master_node`. Requires ops.model('basic', '-ndm', 2, '-ndf', 3) and
    `sec_tag` to already be defined (typically a HankinsonDowel2d section --
    every element below shares the same prototype tag; ZeroLengthSection
    calls the section's getCopy() per element, so each dowel gets its own
    independent state/history).

    The grain direction is assumed uniform across the whole connection (a
    single timber member), so every element uses the same '-orient 1 0 0'
    (local x = global x = parallel-to-grain); only the dowels' positions --
    not their orientation -- vary across the group.

    Returns a list of (ele_tag, x, y) in the same order as dowel_coords.
    """
    ops.node(master_node, 0.0, 0.0)
    info = []
    for i, (x, y) in enumerate(dowel_coords):
        gnode = ground_node_start + i
        fnode = free_node_start + i
        etag = ele_tag_start + i
        ops.node(gnode, x, y)
        ops.fix(gnode, 1, 1, 1)
        ops.node(fnode, x, y)
        ops.element('zeroLengthSection', etag, gnode, fnode, sec_tag,
                    '-orient', 1, 0, 0)
        ops.rigidLink('beam', master_node, fnode)
        info.append((etag, x, y))
    return info


def saws_dowel_type(tag, F0, K0, K1, DU=0.05, FI_frac=0.10,
                     R2=-0.05, R3=1.0, R4=0.0, ALPHA=0.8, BETA=1.05):
    """A SAWSMaterial (Folz & Filiatrault CUREE pinching hysteretic model)
    whose monotonic backbone is the identical Foschi exponential envelope
    used by DowelType's '-exponential' option:
        F = (F0 + R1*S0*|d|)(1 - exp(-S0*|d|/F0)),   R1 = K1/K0, S0 = K0
    (confirmed against the HYSTR envelope logic in SAWSMaterial.cpp, which
    is a direct port of the same Foschi curve). Only F0/K0/K1 are
    Hankinson-relevant (they set the 0-deg/90-deg calibration); DU (cap
    displacement) and the R2/R3/R4/ALPHA/BETA pinching/degradation shape
    parameters are representative CUREE-model values, not fit to any
    particular test data. DU=0.05 m keeps the failure displacement DF
    (computed internally from FI/R2/R4) comfortably beyond DU for both
    calibrations here (DF ~ 0.085 m at 0 deg, ~0.10 m at 90 deg).
    """
    R1 = K1 / K0
    FI = FI_frac * F0
    ops.uniaxialMaterial('SAWS', tag, F0, FI, DU, K0,
                          R1, R2, R3, R4, ALPHA, BETA)


# 0-degree (parallel-to-grain) / 90-degree (perpendicular-to-grain)
# calibration -- the Foschi envelope values from Schweigler et al., already
# validated bit-exact against ModelFunctions.slipCurve in
# test_foschi_hankinson_matches_closed_form (there, driving DowelType).
F0_0, K0_0, K1_0 = 19.1, 24460.0, 462.0
F0_90, K0_90, K1_90 = 14.18, 19640.0, 739.0

# 3x3 group geometry (Schweigler's MainAnalysis.py): spacing parallel to
# grain a1, perpendicular a2, edge distances a3/a4, all in metres.
A1, A2, A3, A4 = 0.070, 0.030, 0.080, 0.030
N_ROWS, N_COLUMNS = 3, 3

DOWEL_MAT_0, DOWEL_MAT_90, SEC_TAG, MASTER_NODE = 20, 21, 1, 9000


def build_model():
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    saws_dowel_type(DOWEL_MAT_0, F0_0, K0_0, K1_0)     # 0-deg calibration
    saws_dowel_type(DOWEL_MAT_90, F0_90, K0_90, K1_90)  # 90-deg calibration
    ops.section('HankinsonDowel2d', SEC_TAG, '-responseHankinson',
                DOWEL_MAT_0, DOWEL_MAT_90)

    coords = dowel_coordinates_rp(N_ROWS, N_COLUMNS, A1, A2, A3, A4)
    dowels = build_dowel_group(coords, SEC_TAG, MASTER_NODE)
    return dowels


def apply_normal_force(n_target, pattern_tag=5, n_steps=10):
    """Ramp a constant axial preload N (global/local x = parallel-to-grain,
    per build_dowel_group's '-orient 1 0 0') onto the master node, then
    loadConst() it so it stays in place under a later push_moment() pattern.

    N acts as a pure rigid-body translation of the whole group (no lever
    arm, unlike M), so every dowel picks up the *same* extra parallel-to-
    grain slip. The group's axial stiffness is the 9 dowels' K0_0 acting in
    parallel (~9x24460 = 220,000 kN/m here), so small N (order 10-20 kN)
    only produces ~1e-4 m of preload -- negligible next to the 1e-2-1e-1 m
    of rotation-induced slip a moment pushover produces at the outer
    dowels, and the M-phi curve will look unchanged. N needs to be a
    sizeable fraction of the group's axial capacity (~9*F0_0 = 172 kN here)
    before its effect on the moment response is visible.
    """
    ops.timeSeries('Linear', pattern_tag)
    ops.pattern('Plain', pattern_tag, pattern_tag)
    ops.load(MASTER_NODE, n_target, 0.0, 0.0)
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('UmfPack')
    ops.test('NormDispIncr', 1.0e-8, 25)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0 / n_steps)
    ops.analysis('Static')

    for _ in range(n_steps):
        ok = ops.analyze(1)
        assert ok == 0, "normal-force preload failed to converge"
    ops.loadConst('-time', 0.0)


def push_moment(m_target, n_steps, pattern_tag=1):
    ops.timeSeries('Linear', pattern_tag)
    ops.pattern('Plain', pattern_tag, pattern_tag)
    ops.load(MASTER_NODE, 0.0, 0.0, m_target)
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('UmfPack')
    ops.test('NormDispIncr', 1.0e-8, 25)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0 / n_steps)
    ops.analysis('Static')

    m_hist, phi_hist = [0.0], [0.0]
    for _ in range(n_steps):
        ok = ops.analyze(1)
        assert ok == 0, "pushover failed to converge"
        m_hist.append(m_target * len(phi_hist) / n_steps)
        phi_hist.append(ops.nodeDisp(MASTER_NODE, 3))
    return m_hist, phi_hist


def plot_dowel_group(dowels, ax=None, force_scale=None, title=None):
    """Plot the dowel group's current geometry together with each dowel's
    force vector, read from the CURRENT analysis state (call ops.reactions()
    -- done here internally -- right after the ops.analyze() step you want
    to visualize; no time-stepping happens in this function).

    Each ground node's reaction is the force the timber exerts on the
    dowel; the force the dowel (connection) exerts back on the timber, at
    that dowel's location, is its negative -- that is what gets drawn, so
    arrows point in the direction the connection is "pushing" the timber at
    each dowel, which is the physically intuitive picture of load transfer
    through the group.
    """
    import matplotlib.pyplot as plt
    if ax is None:
        _, ax = plt.subplots(figsize=(6.5, 6.5))

    ops.reactions()
    xs = [x for _, x, y in dowels]
    ys = [y for _, x, y in dowels]
    fx, fy = [], []
    for etag, x, y in dowels:
        gnode = 1000 + (etag - 1)   # ground_node_start=1000, ele_tag_start=1
        Rx, Ry, _ = ops.nodeReaction(gnode)
        fx.append(-Rx)
        fy.append(-Ry)

    fmag = [(u * u + v * v) ** 0.5 for u, v in zip(fx, fy)]
    fmax = max(fmag) if fmag else 0.0
    if force_scale is None:
        span = max(max(xs) - min(xs), max(ys) - min(ys), 1.0e-6)
        force_scale = 0.4 * span / fmax if fmax > 0 else 1.0

    mx, my = ops.nodeCoord(MASTER_NODE)

    ax.scatter(xs, ys, s=70, c='saddlebrown', zorder=3, label='dowels')
    ax.scatter([mx], [my], marker='+', s=140, c='black', zorder=3,
               linewidths=2, label='rigid reference node')
    for x, y in zip(xs, ys):
        ax.plot([mx, x], [my, y], color='0.75', linewidth=0.8,
                linestyle='--', zorder=1)

    ax.quiver(xs, ys, [u * force_scale for u in fx],
              [v * force_scale for v in fy],
              angles='xy', scale_units='xy', scale=1,
              color='crimson', width=0.006, zorder=2,
              label='dowel force on timber')

    for (etag, x, y), m in zip(dowels, fmag):
        ax.annotate("%d\n%.1f kN" % (etag, m), (x, y),
                    textcoords="offset points", xytext=(7, 7), fontsize=7)

    ax.set_aspect('equal', adjustable='datalim')
    ax.set_xlabel('x (m, parallel to grain)')
    ax.set_ylabel('y (m, perpendicular to grain)')
    ax.grid(True, linestyle=':')
    ax.legend(loc='best', fontsize=8)
    if title:
        ax.set_title(title)
    return ax


# Axial preload levels (kN) to sweep in the __main__ demo below. 0/40/80 kN
# span from no preload up to ~46% of the group's axial capacity
# (9*F0_0 = 172 kN) -- large enough for the N-M interaction to be visible,
# while M_TARGET stays inside the moment capacity at every preload level (a
# larger N eats into each dowel's capacity and can push the higher-preload
# curves into non-convergence if M_TARGET is raised much further).
NORMAL_FORCE_LEVELS = [0.0, 40.0, 80.0]
M_TARGET = 5

if __name__ == '__main__':
    try:
        import matplotlib.pyplot as plt
        have_plt = True
        # Dowel-group force-vector panels: one per preload level N, all at
        # the SAME final moment M_TARGET, so the panels isolate the effect
        # of axial preload on the force pattern rather than of rotation.
        # ops holds only ONE live domain at a time (build_model() wipes it),
        # so each panel must be drawn immediately after its own analysis --
        # not after the whole N-sweep, by which point only the last N's
        # converged state would still be queryable.
        fig, axes = plt.subplots(1, len(NORMAL_FORCE_LEVELS), figsize=(16, 5.5))
        if len(NORMAL_FORCE_LEVELS) == 1:
            axes = [axes]
    except ImportError:
        have_plt = False
        print("\nmatplotlib not available -- skipping plots.\n")

    curves = {}
    for i, n_val in enumerate(NORMAL_FORCE_LEVELS):
        dowels = build_model()
        if n_val:
            apply_normal_force(n_val)
        m_hist, phi_hist = push_moment(m_target=M_TARGET, n_steps=80)
        curves[n_val] = (phi_hist, m_hist)

        # Report both an early-stage secant (where the preload's effect on
        # each dowel's Hankinson-blended stiffness is clearest) and the
        # final one (where the response has flattened into a plastic
        # plateau and relative differences shrink even though they're real).
        i_mid = len(m_hist) // 3
        k_mid = m_hist[i_mid] / phi_hist[i_mid]
        k_final = m_hist[-1] / phi_hist[-1]
        print("N = %5.1f kN   K_secant(M=%.1f) = %7.2f   "
              "K_secant(M=%.1f) = %6.2f kN.m/rad"
              % (n_val, m_hist[i_mid], k_mid, m_hist[-1], k_final))

        if have_plt:
            plot_dowel_group(dowels, ax=axes[i], title="N = %.0f kN" % n_val)

    if have_plt:
        fig.suptitle("3x3 dowel group: force transfer under varying axial "
                      "preload, at M = %.1f kN.m"
                      "\n(HankinsonDowel2d + SAWSMaterial, "
                      "-responseHankinson)" % M_TARGET)
        fig.tight_layout()

        fig2, ax2 = plt.subplots(figsize=(6.5, 5))
        for n_val, (phi_hist, m_hist) in curves.items():
            ax2.plot(phi_hist, m_hist, label="N = %.0f kN" % n_val)
        ax2.set_xlabel("Rotation phi (rad)")
        ax2.set_ylabel("Moment M (kN.m)")
        ax2.set_title("3x3 dowel group: M-phi under axial preload")
        ax2.grid(True)
        ax2.legend()

        plt.show()
