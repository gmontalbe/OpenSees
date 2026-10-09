# Worked example: assembling a timber dowel GROUP (not a single dowel) out of
# HankinsonSection2d sections, following the "Building a dowel group" recipe --
# one zeroLengthSection per physical dowel location, all tied to a shared
# rigid master node with rigidLink('beam', ...).
#
# HankinsonSection2d itself only ever sees a single dowel's 2-D in-plane slip;
# it has no notion of a "group". The group geometry (how many dowels, their
# spacing, and therefore their lever arms about the connection's reference
# point) lives entirely in ordinary node coordinates plus the rigidLink
# constraint -- the M-N-V interaction of the whole connection then falls out
# of nodal equilibrium (sum of F_i x r_i) with no additional bookkeeping.
#
# Each dowel's 0-deg/90-deg calibration is Hankinson-blended by
# HankinsonSection2d's parameter-level interpolation (-paramHankinson mat):
# a single SAWSMaterial (the Folz & Filiatrault CUREE pinching hysteretic
# model) has its F0/S0/K1 parameters re-interpolated between the 0-deg and
# 90-deg calibration values -- via SAWSMaterial's setParameter/
# updateParameter -- and pushed back into the SAME material at every commit.
# This is the literal Schweigler formulation: there is one physical dowel
# law whose *parameters* vary with grain angle, not two independent
# materials whose responses get blended (an earlier "response-level" mode
# that did the latter was removed from HankinsonSection2d -- see
# SRC/material/section/HankinsonSection2d.h for why it isn't physical).
#
# The demo below builds the SAME 3x3 dowel-group connection two ways and
# checks they agree to round-off:
#
#   1. by hand: one zeroLengthSection per dowel, all rigidLink('beam', ...)ed
#      to a shared master node -- 9 elements, 9 constraints;
#   2. as ONE element: the same 9 HankinsonSection2d dowels, but wrapped in a
#      single BoltGroupSection2d, which takes the single-dowel section as a
#      PROTOTYPE, copies it once per fastener ('-grid nRows nCols a1 a2'
#      lays out the grid), and internally does what the rigidLink model does
#      with nodes and constraints: map the reference point's (u_x, u_y,
#      theta) onto each fastener by its lever arm and sum the forces back
#      into (N, V, M).
#
# The 3x3 group geometry below reproduces Schweigler et al.'s worked example
# (see /schweigler_implementation/MainAnalysis.py, function dowelCoordinatesRP
# in ModelFunctions.py). Units are consistent with the calibration values
# (force in kN, length in m) -- the same scale used throughout
# tests/test_hankinson_dowel2d.py.
#
# Run standalone:
#   cd EXAMPLES/ExamplePython && cp ../../build/Release/OpenSeesPy.so ./opensees.so
#   python HankinsonDowelGroup.py

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

# 0-degree (parallel-to-grain) / 90-degree (perpendicular-to-grain)
# calibration -- the Foschi envelope values from Schweigler et al., already
# validated bit-exact against ModelFunctions.slipCurve in
# test_foschi_hankinson_matches_closed_form (there, driving DowelType):
#     F = (F0 + K1*|d|)(1 - exp(-S0*|d|/F0)),   S0 = K0
F0_0, K0_0, K1_0 = 19.1, 24460.0, 462.0
F0_90, K0_90, K1_90 = 14.18, 19640.0, 739.0

# SAWSMaterial (Folz & Filiatrault CUREE pinching hysteretic model) shape
# parameters: DU (cap displacement) and R2/R3/R4/ALPHA/BETA are
# representative CUREE-model values, not fit to any particular test data.
# DU=0.05 m keeps the failure displacement DF (computed internally from
# FI/R2/R4) comfortably beyond DU for both calibrations here (DF ~ 0.085 m
# at 0 deg, ~0.10 m at 90 deg). F0/K0(S0)/K1 are the three Hankinson-relevant
# parameters -- SAWSMaterial exposes them directly via setParameter/
# updateParameter as "F0"/"S0"/"K1" (K1 rather than only the ratio R1, since
# Hankinson-interpolating R1 and S0 separately does not give the same result
# as interpolating K1 itself).
SAWS_DU, SAWS_FI_FRAC = 0.05, 0.10
SAWS_R2, SAWS_R3, SAWS_R4, SAWS_ALPHA, SAWS_BETA = -0.05, 1.0, 0.0, 0.8, 1.05

# 3x3 group geometry (Schweigler's MainAnalysis.py): spacing parallel to
# grain a1, perpendicular a2, edge distances a3/a4, all in metres.
A1, A2, A3, A4 = 0.070, 0.030, 0.080, 0.030
N_ROWS, N_COLUMNS = 3, 3

DOWEL_MAT, SEC_TAG, MASTER_NODE = 20, 1, 9000
GROUP_SEC_TAG, GROUP_ELE_TAG, GROUP_GROUND_NODE = 2, 500, 500

# Axial preload levels (kN) to sweep below. 0/40/80 kN span from no preload
# up to ~46% of the group's axial capacity (9*F0_0 = 172 kN) -- large enough
# for the N-M interaction to be visible, while M_TARGET stays inside the
# moment capacity at every preload level (a larger N eats into each dowel's
# capacity and can push the higher-preload curves into non-convergence if
# M_TARGET is raised much further).
NORMAL_FORCE_LEVELS = [0.0, 40.0, 80.0]
M_TARGET = 10

# Grid of dowel (x, y) coordinates referred to the group's geometric centre
# (the connection's reference point RP). x runs parallel to grain at spacing
# A1 (column pitch), y runs perpendicular to grain at spacing A2 (row
# pitch); A3/A4 are the edge distances of the first row/column. Mirrors
# ModelFunctions.dowelCoordinatesRP in the Python reference. This is fixed
# geometry (the constants above never change through the file), so it is
# computed once and reused everywhere below.
raw_coords = [(A3 + A1 * col, A4 + A2 * row)
              for row in range(N_ROWS) for col in range(N_COLUMNS)]
x_cg = A3 + (0.5 * N_COLUMNS - 0.5) * A1
y_cg = A4 + (0.5 * N_ROWS - 0.5) * A2
coords = [(x - x_cg, y - y_cg) for x, y in raw_coords]

try:
    import matplotlib.pyplot as plt
    have_plt = True
    # Dowel-group force-vector panels: one per preload level N, all at the
    # SAME final moment M_TARGET, so the panels isolate the effect of axial
    # preload on the force pattern rather than of rotation.
    fig, axes = plt.subplots(1, len(NORMAL_FORCE_LEVELS), figsize=(16, 5.5))
    if len(NORMAL_FORCE_LEVELS) == 1:
        axes = [axes]
except ImportError:
    have_plt = False
    print("\nmatplotlib not available -- skipping plots.\n")

# ---------------------------------------------------------------------
# Cross-check, at the hardest preload level: the hand-built 9-element +
# rigidLink model must agree with the single BoltGroupSection2d element to
# round-off (both are linear in the rotation) -- the demo itself
# demonstrates that the two model forms are the same mechanics.
# ---------------------------------------------------------------------

# --- 9 elements + rigidLink -----------------------------------------------
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

# One SAWSMaterial, started at the 0-deg calibration; HankinsonSection2d
# re-interpolates its F0/S0/K1 between the 0-deg and 90-deg values below
# (via setParameter/updateParameter) at every commit.
R1 = K1_0 / K0_0
FI = SAWS_FI_FRAC * F0_0
ops.uniaxialMaterial('SAWS', DOWEL_MAT, F0_0, FI, SAWS_DU, K0_0,
                      R1, SAWS_R2, SAWS_R3, SAWS_R4, SAWS_ALPHA, SAWS_BETA)
ops.section('HankinsonSection2d', SEC_TAG, '-paramHankinson', DOWEL_MAT,
            '-param', 'F0', F0_0, F0_90,
            '-param', 'S0', K0_0, K0_90,
            '-param', 'K1', K1_0, K1_90)

# One zeroLengthSection element (hosting SEC_TAG) per dowel, its free node
# rigid-linked to the master node. The grain direction is assumed uniform
# across the whole connection (a single timber member), so every element
# uses the same '-orient 1 0 0' (local x = global x = parallel-to-grain);
# only the dowels' positions -- not their orientation -- vary across the
# group. Requires ops.model('basic', '-ndm', 2, '-ndf', 3).
ops.node(MASTER_NODE, 0.0, 0.0)
for i, (x, y) in enumerate(coords):
    gnode = 1000 + i
    fnode = 2000 + i
    etag = 1 + i
    ops.node(gnode, x, y)
    ops.fix(gnode, 1, 1, 1)
    ops.node(fnode, x, y)
    ops.element('zeroLengthSection', etag, gnode, fnode, SEC_TAG,
                '-orient', 1, 0, 0)
    ops.rigidLink('beam', MASTER_NODE, fnode)

# Ramp a constant axial preload N (global/local x = parallel-to-grain) onto
# the master node, then loadConst() it so it stays in place under the
# push_moment pattern below. N acts as a pure rigid-body translation of the
# whole group (no lever arm, unlike M), so every dowel picks up the *same*
# extra parallel-to-grain slip. The group's axial stiffness is the 9 dowels'
# K0_0 acting in parallel (~9x24460 = 220,000 kN/m here), so N needs to be a
# sizeable fraction of the group's axial capacity (~9*F0_0 = 172 kN here)
# before its effect on the moment response is visible -- hence 80 kN here.
ops.timeSeries('Linear', 5)
ops.pattern('Plain', 5, 5)
ops.load(MASTER_NODE, 80.0, 0.0, 0.0)
ops.wipeAnalysis()
ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('UmfPack')
ops.test('NormDispIncr', 1.0e-8, 25)
ops.algorithm('Newton')
ops.integrator('LoadControl', 1.0 / 10)
ops.analysis('Static')
for _ in range(10):
    ok = ops.analyze(1)
    assert ok == 0, "normal-force preload failed to converge"
ops.loadConst('-time', 0.0)

# Push the in-plane moment M (about local z), 80 LoadControl steps.
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(MASTER_NODE, 0.0, 0.0, M_TARGET)
ops.wipeAnalysis()
ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('UmfPack')
ops.test('NormDispIncr', 1.0e-8, 25)
ops.algorithm('Newton')
ops.integrator('LoadControl', 1.0 / 80)
ops.analysis('Static')
phi_hist = [0.0]
for _ in range(80):
    ok = ops.analyze(1)
    assert ok == 0, "pushover failed to converge"
    phi_hist.append(ops.nodeDisp(MASTER_NODE, 3))
phi_multi = phi_hist[-1]

# --- the same connection as ONE BoltGroupSection2d element ----------------
# This replaces 9 elements plus 9 rigidLink constraints with 1 element and
# 3 free DOF, and -- because there are no multi-point constraints left -- it
# also sidesteps the TransformationConstraintHandler problem documented in
# HankinsonDowelGroup3d.py's apply_normal_force step.
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

R1 = K1_0 / K0_0
FI = SAWS_FI_FRAC * F0_0
ops.uniaxialMaterial('SAWS', DOWEL_MAT, F0_0, FI, SAWS_DU, K0_0,
                      R1, SAWS_R2, SAWS_R3, SAWS_R4, SAWS_ALPHA, SAWS_BETA)
ops.section('HankinsonSection2d', SEC_TAG, '-paramHankinson', DOWEL_MAT,
            '-param', 'F0', F0_0, F0_90,
            '-param', 'S0', K0_0, K0_90,
            '-param', 'K1', K1_0, K1_90)
ops.section('BoltGroupSection2d', GROUP_SEC_TAG, '-section', SEC_TAG,
            '-grid', N_ROWS, N_COLUMNS, A1, A2)

ops.node(GROUP_GROUND_NODE, 0.0, 0.0)
ops.fix(GROUP_GROUND_NODE, 1, 1, 1)
ops.node(MASTER_NODE, 0.0, 0.0)
ops.element('zeroLengthSection', GROUP_ELE_TAG,
            GROUP_GROUND_NODE, MASTER_NODE, GROUP_SEC_TAG,
            '-orient', 1, 0, 0)

ops.timeSeries('Linear', 5)
ops.pattern('Plain', 5, 5)
ops.load(MASTER_NODE, 80.0, 0.0, 0.0)
ops.wipeAnalysis()
ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('UmfPack')
ops.test('NormDispIncr', 1.0e-8, 25)
ops.algorithm('Newton')
ops.integrator('LoadControl', 1.0 / 10)
ops.analysis('Static')
for _ in range(10):
    ok = ops.analyze(1)
    assert ok == 0, "normal-force preload failed to converge"
ops.loadConst('-time', 0.0)

ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(MASTER_NODE, 0.0, 0.0, M_TARGET)
ops.wipeAnalysis()
ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('UmfPack')
ops.test('NormDispIncr', 1.0e-8, 25)
ops.algorithm('Newton')
ops.integrator('LoadControl', 1.0 / 80)
ops.analysis('Static')
phi_hist = [0.0]
for _ in range(80):
    ok = ops.analyze(1)
    assert ok == 0, "pushover failed to converge"
    phi_hist.append(ops.nodeDisp(MASTER_NODE, 3))
phi_single = phi_hist[-1]

print("cross-check at N = %.0f kN: phi(9 elements + rigidLink) = %.12e"
      % (NORMAL_FORCE_LEVELS[-1], phi_multi))
print("                            phi(1 BoltGroupSection2d)  = %.12e"
      % phi_single)
print("                            difference                 = %.3e\n"
      % abs(phi_multi - phi_single))

# ---------------------------------------------------------------------
# The BoltGroupSection2d model swept over axial preload, at the same
# M_TARGET each time -- so the panels below isolate the effect of preload
# on the force pattern rather than of rotation.
# ---------------------------------------------------------------------

curves = {}
for i, n_val in enumerate(NORMAL_FORCE_LEVELS):
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    R1 = K1_0 / K0_0
    FI = SAWS_FI_FRAC * F0_0
    ops.uniaxialMaterial('SAWS', DOWEL_MAT, F0_0, FI, SAWS_DU, K0_0,
                          R1, SAWS_R2, SAWS_R3, SAWS_R4, SAWS_ALPHA, SAWS_BETA)
    ops.section('HankinsonSection2d', SEC_TAG, '-paramHankinson', DOWEL_MAT,
                '-param', 'F0', F0_0, F0_90,
                '-param', 'S0', K0_0, K0_90,
                '-param', 'K1', K1_0, K1_90)
    ops.section('BoltGroupSection2d', GROUP_SEC_TAG, '-section', SEC_TAG,
                '-grid', N_ROWS, N_COLUMNS, A1, A2)

    ops.node(GROUP_GROUND_NODE, 0.0, 0.0)
    ops.fix(GROUP_GROUND_NODE, 1, 1, 1)
    ops.node(MASTER_NODE, 0.0, 0.0)
    ops.element('zeroLengthSection', GROUP_ELE_TAG,
                GROUP_GROUND_NODE, MASTER_NODE, GROUP_SEC_TAG,
                '-orient', 1, 0, 0)

    if n_val:
        ops.timeSeries('Linear', 5)
        ops.pattern('Plain', 5, 5)
        ops.load(MASTER_NODE, n_val, 0.0, 0.0)
        ops.wipeAnalysis()
        ops.constraints('Transformation')
        ops.numberer('RCM')
        ops.system('UmfPack')
        ops.test('NormDispIncr', 1.0e-8, 25)
        ops.algorithm('Newton')
        ops.integrator('LoadControl', 1.0 / 10)
        ops.analysis('Static')
        for _ in range(10):
            ok = ops.analyze(1)
            assert ok == 0, "normal-force preload failed to converge"
        ops.loadConst('-time', 0.0)

    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(MASTER_NODE, 0.0, 0.0, M_TARGET)
    ops.wipeAnalysis()
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('UmfPack')
    ops.test('NormDispIncr', 1.0e-8, 25)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0 / 80)
    ops.analysis('Static')
    m_hist, phi_hist = [0.0], [0.0]
    for _ in range(80):
        ok = ops.analyze(1)
        assert ok == 0, "pushover failed to converge"
        m_hist.append(M_TARGET * len(phi_hist) / 80)
        phi_hist.append(ops.nodeDisp(MASTER_NODE, 3))
    curves[n_val] = (phi_hist, m_hist)

    # Report both an early-stage secant (where the preload's effect on each
    # dowel's Hankinson-blended stiffness is clearest) and the final one
    # (where the response has flattened into a plastic plateau and relative
    # differences shrink even though they're real).
    i_mid = len(m_hist) // 3
    k_mid = m_hist[i_mid] / phi_hist[i_mid]
    k_final = m_hist[-1] / phi_hist[-1]
    print("N = %5.1f kN   K_secant(M=%.1f) = %7.2f   "
          "K_secant(M=%.1f) = %6.2f kN.m/rad"
          % (n_val, m_hist[i_mid], k_mid, m_hist[-1], k_final))

    if have_plt:
        # The group's per-fastener force pattern at the end of this preload
        # level's pushover, read from the group section's own 'boltForces'
        # response (the concatenated (Fx, Fy) pairs -- note the plural:
        # HankinsonSection2d intercepts the SINGULAR 'force' with its own
        # scalar slip-law force, so only 'forces' gives the (P, Vy) pair).
        ax = axes[i]
        f = ops.eleResponse(GROUP_ELE_TAG, 'section', 'boltForces')
        xs = [x for x, _ in coords]
        ys = [y for _, y in coords]
        fx = [f[2 * k] for k in range(len(coords))]
        fy = [f[2 * k + 1] for k in range(len(coords))]

        fmag = [(u * u + v * v) ** 0.5 for u, v in zip(fx, fy)]
        fmax = max(fmag) if fmag else 0.0
        span = max(max(xs) - min(xs), max(ys) - min(ys), 1.0e-6)
        force_scale = 0.4 * span / fmax if fmax > 0 else 1.0

        mx, my = ops.nodeCoord(MASTER_NODE)
        cr = ops.eleResponse(GROUP_ELE_TAG, 'section', 'centerOfRotation')

        ax.scatter(xs, ys, s=70, c='saddlebrown', zorder=3, label='dowels')
        ax.scatter([mx], [my], marker='+', s=140, c='black', zorder=3,
                   linewidths=2, label='reference point')
        if cr and abs(cr[0]) < 10 * max(abs(v) for v in xs + ys):
            ax.scatter([cr[0]], [cr[1]], marker='x', s=110, c='tab:blue',
                       zorder=3, linewidths=2, label='centre of rotation')
        for x, y in zip(xs, ys):
            ax.plot([mx, x], [my, y], color='0.75', linewidth=0.8,
                    linestyle='--', zorder=1)

        ax.quiver(xs, ys, [u * force_scale for u in fx],
                  [v * force_scale for v in fy],
                  angles='xy', scale_units='xy', scale=1,
                  color='crimson', width=0.006, zorder=2,
                  label='dowel force')

        for j, ((x, y), m) in enumerate(zip(coords, fmag)):
            ax.annotate("%d\n%.1f kN" % (j, m), (x, y),
                        textcoords="offset points", xytext=(7, 7), fontsize=7)

        ax.set_aspect('equal', adjustable='datalim')
        ax.set_xlabel('x (m, parallel to grain)')
        ax.set_ylabel('y (m, perpendicular to grain)')
        ax.grid(True, linestyle=':')
        ax.legend(loc='best', fontsize=8)
        ax.set_title("N = %.0f kN" % n_val)

if have_plt:
    fig.suptitle("3x3 dowel group as ONE element: force transfer under "
                  "varying axial preload, at M = %.1f kN.m"
                  "\n(BoltGroupSection2d over a HankinsonSection2d "
                  "prototype + SAWSMaterial)" % M_TARGET)
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
