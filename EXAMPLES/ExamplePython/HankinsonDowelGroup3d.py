# Worked example: the same 3x3 timber dowel GROUP as HankinsonDowelGroup.py,
# but embedded in a full 3D model (ndm=3, ndf=6) instead of the planar
# ndm=2, ndf=3 model used there.
#
# HankinsonSection2d is inherently a 2-DOF coupled section (getOrder()==2,
# codes {P, Vy} by default): hosted in a zeroLengthSection element with
# ndf=6, it resists ONLY the two in-plane translational relative DOFs it
# has codes for -- grain-parallel P (local x) and in-plane-perpendicular Vy
# (local y') -- between the element's two nodes. It contributes exactly
# ZERO stiffness to the third (out-of-plane, Vz) translation and to all
# three relative rotations (My, Mz, T): ZeroLengthSection::setTransformation
# (SRC/element/zeroLength/ZeroLengthSection.cpp) only ever fills the rows of
# its element-deformation transformation matrix A that correspond to the
# section's own getType() codes, so directions the section doesn't report
# never get any force or stiffness at all.
#
# That is physically appropriate for a real dowel-type fastener -- it
# transfers load by bearing against the timber in its own shear plane and
# has no meaningful resistance to being pulled or rotated out of that
# plane -- but it does mean this model, unlike the 2D one, must explicitly
# restrain the group's out-of-plane DOFs itself. Nothing in
# HankinsonSection2d or zeroLengthSection does that for you: leaving those
# DOFs free gives a singular stiffness matrix (the analysis fails to
# factorize, it does not silently do the wrong thing). See the master
# node's ops.fix() call below for where and why.
#
# Run standalone:
#   cd EXAMPLES/ExamplePython && cp ../../build/opensees.so .
#   python HankinsonDowelGroup3d.py

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

# 0-degree (parallel-to-grain) / 90-degree (perpendicular-to-grain)
# calibration -- the Foschi envelope values from Schweigler et al.:
#     F = (F0 + K1*|d|)(1 - exp(-S0*|d|/F0)),   S0 = K0
F0_0, K0_0, K1_0 = 19.1, 24460.0, 462.0
F0_90, K0_90, K1_90 = 14.18, 19640.0, 739.0

# SAWSMaterial (Folz & Filiatrault CUREE pinching hysteretic model) shape
# parameters -- see HankinsonDowelGroup.py for the rationale of these
# values. F0/K0(S0)/K1 are the three Hankinson-relevant parameters,
# re-interpolated between the 0-deg and 90-deg values below (via
# HankinsonSection2d's '-paramHankinson') at every commit.
SAWS_DU, SAWS_FI_FRAC = 0.05, 0.10
SAWS_R2, SAWS_R3, SAWS_R4, SAWS_ALPHA, SAWS_BETA = -0.05, 1.0, 0.0, 0.8, 1.05

# 3x3 group geometry (Schweigler's MainAnalysis.py): spacing parallel to
# grain a1, perpendicular a2, edge distances a3/a4, all in metres. The
# connection plane is the global X-Y plane: grain parallel to X, in-plane
# perpendicular to Y, out-of-plane normal Z.
A1, A2, A3, A4 = 0.070, 0.030, 0.080, 0.030
N_ROWS, N_COLUMNS = 3, 3

DOWEL_MAT, SEC_TAG, MASTER_NODE = 20, 1, 9000
GROUP_SEC_TAG, GROUP_ELE_TAG, GROUP_GROUND_NODE = 2, 500, 500

# Axial preload levels (kN) to sweep below -- identical choice to
# HankinsonDowelGroup.py, see its NORMAL_FORCE_LEVELS comment for the
# sizing rationale (same calibration, same 3x3 group).
NORMAL_FORCE_LEVELS = [0.0, 40.0, 80.0]
M_TARGET = 10

# Grid of dowel (x, y) coordinates referred to the group's geometric centre
# (the connection's reference point RP), computed once since the geometry
# constants above never change through the file. Mirrors
# ModelFunctions.dowelCoordinatesRP in the Python reference.
raw_coords = [(A3 + A1 * col, A4 + A2 * row)
              for row in range(N_ROWS) for col in range(N_COLUMNS)]
x_cg = A3 + (0.5 * N_COLUMNS - 0.5) * A1
y_cg = A4 + (0.5 * N_ROWS - 0.5) * A2
coords = [(x - x_cg, y - y_cg) for x, y in raw_coords]

try:
    import matplotlib.pyplot as plt
    have_plt = True
    fig, axes = plt.subplots(1, len(NORMAL_FORCE_LEVELS), figsize=(16, 5.5))
    if len(NORMAL_FORCE_LEVELS) == 1:
        axes = [axes]
except ImportError:
    have_plt = False
    print("\nmatplotlib not available -- skipping plots.\n")

# ---------------------------------------------------------------------
# Cross-check, at the hardest preload level: the hand-built 9-element +
# rigidLink model must agree with the single BoltGroupSection2d element to
# round-off (both are linear in the rotation).
# ---------------------------------------------------------------------

# --- 9 elements + rigidLink, 3D -------------------------------------------
ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 6)

R1 = K1_0 / K0_0
FI = SAWS_FI_FRAC * F0_0
ops.uniaxialMaterial('SAWS', DOWEL_MAT, F0_0, FI, SAWS_DU, K0_0,
                      R1, SAWS_R2, SAWS_R3, SAWS_R4, SAWS_ALPHA, SAWS_BETA)
ops.section('HankinsonSection2d', SEC_TAG, '-paramHankinson', DOWEL_MAT,
            '-param', 'F0', F0_0, F0_90,
            '-param', 'S0', K0_0, K0_90,
            '-param', 'K1', K1_0, K1_90)

# One zeroLengthSection element per dowel, all rigid-linked to a shared
# master node -- the 3D counterpart of the 2D layout, but every node now
# carries 6 DOF. '-orient' must give both the local x and y' vectors
# explicitly here -- 3D zeroLengthSection has no auto y = z cross x
# shortcut the way the 2D case does (see ZeroLengthSection's
# OPS_ZeroLengthSection parser). Each ground (timber-side) node is fixed in
# all 6 DOF; the free (fastener-side) node's out-of-plane translation and
# all three rotations are left unfixed here, but since HankinsonSection2d
# never resists them, they end up entirely governed by the rigidLink
# constraint -- i.e. slaved to the master node's own out-of-plane DOFs, not
# resisted by anything. Restraining those AT the master node (below) is
# what actually gives them stiffness.
ops.node(MASTER_NODE, 0.0, 0.0, 0.0)
for i, (x, y) in enumerate(coords):
    gnode = 1000 + i
    fnode = 2000 + i
    etag = 1 + i
    ops.node(gnode, x, y, 0.0)
    ops.fix(gnode, 1, 1, 1, 1, 1, 1)
    ops.node(fnode, x, y, 0.0)
    ops.element('zeroLengthSection', etag, gnode, fnode, SEC_TAG,
                '-orient', 1, 0, 0, 0, 1, 0)
    ops.rigidLink('beam', MASTER_NODE, fnode)

# No dowel resists the master node's out-of-plane translation (Uz) or any
# rotation (Rx, Ry, Rz) directly. The in-plane rotation Rz DOES pick up
# stiffness despite that, but only indirectly: it comes from the lever-arm
# coupling of each dowel's in-plane (Fx, Fy) force through the rigidLink
# constraint, exactly as in the 2D example, not from any element resisting
# Rz itself. Uz, Rx and Ry get NO stiffness from anywhere in this model and
# must be restrained explicitly or the analysis is singular in those
# directions. Physically this fixation stands in for whatever restrains the
# connection out of its own plane in the real structure (the timber
# member's own flexural continuity, other fasteners, a test rig).
ops.fix(MASTER_NODE, 0, 0, 1, 1, 1, 0)

# Ramp a constant axial preload N onto the master node, then loadConst() it
# so it stays in place under the moment pushover below -- see
# HankinsonDowelGroup.py for the sizing rationale (same calibration, same
# 3x3 group).
#
# 'Lagrange', NOT 'Transformation': in this build, TransformationConstraintHandler
# mishandles a node that is BOTH (a) partially ops.fix()ed and (b) the
# RETAINED node of an MP_Constraint with other, still-free DOFs -- which is
# exactly MASTER_NODE here (its Uz/Rx/Ry are fixed above, every dowel's
# rigidLink('beam', MASTER_NODE, ...) retains it, and Ux/Uy/Rz stay free).
# Verified in total isolation (one zeroLengthSection + one rigidLink + a
# Bidirectional section, no Hankinson code involved at all, in BOTH 2D and
# 3D): across n repeated LoadControl analyze(1) calls to the same
# cumulative load, the master's OWN reported displacement comes out
# (n+1)/2 times too large, while the slave node's displacement (and every
# reaction force) stays exactly correct throughout -- so the bug is
# invisible in equilibrium/force checks and only shows up when you read the
# master's own displacement after a multi-step analysis. A master node with
# NO fixed DOFs (as in the 2D example, which never fixes MASTER_NODE) does
# not trigger it. Only 'Transformation' is affected; 'Lagrange' reproduces
# the theoretical/2D-model rotational stiffness exactly for any number of
# steps.
ops.timeSeries('Linear', 5)
ops.pattern('Plain', 5, 5)
ops.load(MASTER_NODE, 80.0, 0.0, 0.0, 0.0, 0.0, 0.0)
ops.wipeAnalysis()
ops.constraints('Lagrange')
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

# Push the in-plane moment Mz (about the global/local Z axis, out of the
# dowel group's X-Y plane) -- DOF 6 of 6 here rather than DOF 3 of 3 in the
# 2D example.
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(MASTER_NODE, 0.0, 0.0, 0.0, 0.0, 0.0, M_TARGET)
ops.wipeAnalysis()
ops.constraints('Lagrange')
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
    phi_hist.append(ops.nodeDisp(MASTER_NODE, 6))
phi_multi = phi_hist[-1]

# --- the same connection as ONE BoltGroupSection2d element, 3D ------------
# The group plane is placed entirely by the element's '-orient' (local x
# along the grain, local y' in-plane perpendicular); the fastener offsets
# inside the section are always expressed in that local frame. Nothing else
# about the section is 3D-aware: it still reports {P, Vy, Mz}, so Uz/Rx/Ry
# still get no stiffness from anywhere and must be restrained, exactly as
# above. Two things the hand-built model needs and this one does not: there
# are no rigidLink constraints, so 'Transformation' is safe again (see the
# note above on why it is not, in general, for a partially-fixed retained
# node), and the model has 3 free DOF instead of 9 nodes' worth plus 9
# constraints.
ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 6)

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

ops.node(GROUP_GROUND_NODE, 0.0, 0.0, 0.0)
ops.fix(GROUP_GROUND_NODE, 1, 1, 1, 1, 1, 1)
ops.node(MASTER_NODE, 0.0, 0.0, 0.0)
ops.element('zeroLengthSection', GROUP_ELE_TAG,
            GROUP_GROUND_NODE, MASTER_NODE, GROUP_SEC_TAG,
            '-orient', 1, 0, 0, 0, 1, 0)
ops.fix(MASTER_NODE, 0, 0, 1, 1, 1, 0)

ops.timeSeries('Linear', 5)
ops.pattern('Plain', 5, 5)
ops.load(MASTER_NODE, 80.0, 0.0, 0.0, 0.0, 0.0, 0.0)
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
ops.load(MASTER_NODE, 0.0, 0.0, 0.0, 0.0, 0.0, M_TARGET)
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
    phi_hist.append(ops.nodeDisp(MASTER_NODE, 6))
phi_single = phi_hist[-1]

print("cross-check at N = %.0f kN: phi(9 elements + rigidLink) = %.12e"
      % (NORMAL_FORCE_LEVELS[-1], phi_multi))
print("                            phi(1 BoltGroupSection2d)  = %.12e"
      % phi_single)
print("                            difference                 = %.3e\n"
      % abs(phi_multi - phi_single))

# ---------------------------------------------------------------------
# The BoltGroupSection2d model swept over axial preload, at the same
# M_TARGET each time.
# ---------------------------------------------------------------------

curves = {}
for i, n_val in enumerate(NORMAL_FORCE_LEVELS):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)

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

    ops.node(GROUP_GROUND_NODE, 0.0, 0.0, 0.0)
    ops.fix(GROUP_GROUND_NODE, 1, 1, 1, 1, 1, 1)
    ops.node(MASTER_NODE, 0.0, 0.0, 0.0)
    ops.element('zeroLengthSection', GROUP_ELE_TAG,
                GROUP_GROUND_NODE, MASTER_NODE, GROUP_SEC_TAG,
                '-orient', 1, 0, 0, 0, 1, 0)
    ops.fix(MASTER_NODE, 0, 0, 1, 1, 1, 0)

    if n_val:
        ops.timeSeries('Linear', 5)
        ops.pattern('Plain', 5, 5)
        ops.load(MASTER_NODE, n_val, 0.0, 0.0, 0.0, 0.0, 0.0)
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
    ops.load(MASTER_NODE, 0.0, 0.0, 0.0, 0.0, 0.0, M_TARGET)
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
        phi_hist.append(ops.nodeDisp(MASTER_NODE, 6))
    curves[n_val] = (phi_hist, m_hist)

    i_mid = len(m_hist) // 3
    k_mid = m_hist[i_mid] / phi_hist[i_mid]
    k_final = m_hist[-1] / phi_hist[-1]
    print("N = %5.1f kN   K_secant(M=%.1f) = %7.2f   "
          "K_secant(M=%.1f) = %6.2f kN.m/rad"
          % (n_val, m_hist[i_mid], k_mid, m_hist[-1], k_final))

    if have_plt:
        # The 2D group section's own 'boltForces' response is dimension-
        # agnostic (only Fx/Fy per fastener are ever nonzero -- see the
        # module docstring), so the same force-vector picture as the 2D
        # example applies unchanged here.
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

        mx, my, _ = ops.nodeCoord(MASTER_NODE)
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
        ax.set_ylabel('y (m, perpendicular to grain, in-plane)')
        ax.grid(True, linestyle=':')
        ax.legend(loc='best', fontsize=8)
        ax.set_title("N = %.0f kN" % n_val)

if have_plt:
    fig.suptitle("3x3 dowel group as ONE element in a 3D model (ndm=3, "
                  "ndf=6): force transfer under varying axial preload, at "
                  "M = %.1f kN.m\n(BoltGroupSection2d over a "
                  "HankinsonSection2d prototype; out-of-plane DOFs fixed "
                  "at the master node)" % M_TARGET)
    fig.tight_layout()

    fig2, ax2 = plt.subplots(figsize=(6.5, 5))
    for n_val, (phi_hist, m_hist) in curves.items():
        ax2.plot(phi_hist, m_hist, label="N = %.0f kN" % n_val)
    ax2.set_xlabel("Rotation phi (rad, about out-of-plane Z axis)")
    ax2.set_ylabel("Moment M (kN.m)")
    ax2.set_title("3x3 dowel group (3D model): M-phi under axial preload")
    ax2.grid(True)
    ax2.legend()

    plt.show()
