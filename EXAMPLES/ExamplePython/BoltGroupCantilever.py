# Worked example: a timber cantilever whose BASE CONNECTION is a dowel group,
# modelled as one BoltGroupSection2d element.
#
# This is the smallest model in which the group does something structurally
# interesting: the connection carries the full (N, V, M) the cantilever
# delivers to it, and because it is flexible the tip deflection is no longer
# just the beam's own flexure. Three things add up at the tip:
#
#     delta_tip  =  u_y(connection)  +  theta(connection) * L  +  delta_flexure
#
# and for a realistic dowel group the middle term dominates -- which is the
# whole reason connection flexibility is worth modelling at all. The demo
# below reports that split explicitly, and compares against the same beam
# with its base fully fixed.
#
# Model (all one line of nodes, the connection is a zero-length element):
#
#     node 1 ==[ BoltGroupSection2d ]== node 2 -----[ beam ]----- node 3
#      (fixed)     zeroLengthSection                              (tip, loaded)
#     x = 0            x = 0                                      x = L
#
# The fastener law is the Schweigler et al. calibration used throughout
# HankinsonDowelGroup.py (SAWSMaterial + HankinsonSection2d, parameters
# Hankinson-interpolated between a 0-deg and a 90-deg calibration by each
# fastener's own evolving load-to-grain angle). The group geometry here is
# NOT Schweigler's 3x3 test specimen but a 5x3 pattern sized to the beam
# depth, since a cantilever connection wants its fasteners spread through the
# depth to pick up moment.
#
# Units: kN, m (consistent with the calibration in HankinsonDowelGroup.py).
#
# Run standalone:
#   cd EXAMPLES/ExamplePython && cp ../../build/OpenSeesPy.so ./opensees.so
#   python BoltGroupCantilever.py

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

# ---------------------------------------------------------------------
# Geometry and properties
# ---------------------------------------------------------------------

L = 1.0                      # cantilever span (m)
B, H = 0.200, 0.320          # beam width and depth (m)
E_TIMBER = 11.5e6            # kN/m2 (11500 MPa, glulam order of magnitude)

A_BEAM = B * H
IZ_BEAM = B * H ** 3 / 12.0

# Connection: 5 rows x 3 columns. a1 is the spacing along the beam axis
# (parallel to grain), a2 the spacing through the depth. The outermost
# fasteners sit at y = +-0.12 m, comfortably inside the 0.32 m depth.
N_ROWS, N_COLS = 5, 3
A1, A2 = 0.070, 0.060

# The fastener law: the Schweigler dowel calibration, Hankinson-interpolated
# between a 0-deg (parallel-to-grain) and 90-deg (perpendicular-to-grain)
# Foschi envelope by each fastener's own evolving load-to-grain angle. Values
# already validated bit-exact against the closed-form Foschi curve in
# test_foschi_hankinson_matches_closed_form. SAWSMaterial (Folz &
# Filiatrault CUREE pinching hysteretic model) is used for the monotonic
# backbone: F = (F0 + K1*|d|)(1 - exp(-S0*|d|/F0)), S0 = K0. DU/R2/R3/R4/
# ALPHA/BETA are representative CUREE-model shape values, not fit to any
# particular test data; DU=0.05 m keeps the failure displacement comfortably
# beyond it for both calibrations.
F0_0, K0_0, K1_0 = 19.1, 24460.0, 462.0
F0_90, K0_90, K1_90 = 14.18, 19640.0, 739.0
SAWS_DU, SAWS_FI_FRAC = 0.05, 0.10
SAWS_R2, SAWS_R3, SAWS_R4, SAWS_ALPHA, SAWS_BETA = -0.05, 1.0, 0.0, 0.8, 1.05

MAT_DOWEL, SEC_BOLT, SEC_GROUP = 20, 10, 11
TRANSF = 1
N_FIX, N_CONN, N_TIP = 1, 2, 3
ELE_CONN, ELE_BEAM = 1, 2

P_MAX = 30.0                 # tip load (kN), enough to yield the group
N_STEPS = 120

# Fastener (x, y) offsets referred to the group's geometric centre, for
# labelling output -- edge distances (a3, a4) are zero here since only the
# centroid-referred grid ('-grid nRows nCols a1 a2') matters to the section.
raw_coords = [(A1 * col, A2 * row) for row in range(N_ROWS) for col in range(N_COLS)]
x_cg = (0.5 * N_COLS - 0.5) * A1
y_cg = (0.5 * N_ROWS - 0.5) * A2
coords = [(x - x_cg, y - y_cg) for x, y in raw_coords]

# ---------------------------------------------------------------------
# Flexible base: cantilever on the dowel-group connection
# ---------------------------------------------------------------------

ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

ops.node(N_FIX, 0.0, 0.0)
ops.fix(N_FIX, 1, 1, 1)
ops.node(N_TIP, L, 0.0)

# One fastener: the Schweigler dowel law, Hankinson-interpolated.
R1 = K1_0 / K0_0
FI = SAWS_FI_FRAC * F0_0
ops.uniaxialMaterial('SAWS', MAT_DOWEL, F0_0, FI, SAWS_DU, K0_0,
                      R1, SAWS_R2, SAWS_R3, SAWS_R4, SAWS_ALPHA, SAWS_BETA)
ops.section('HankinsonSection2d', SEC_BOLT, '-paramHankinson', MAT_DOWEL,
            '-param', 'F0', F0_0, F0_90,
            '-param', 'S0', K0_0, K0_90,
            '-param', 'K1', K1_0, K1_90)
# The whole group: that prototype copied onto the grid below. Local x (the
# element's '-orient' below) is the beam axis and the grain direction;
# local y runs through the beam depth.
ops.section('BoltGroupSection2d', SEC_GROUP, '-section', SEC_BOLT,
            '-grid', N_ROWS, N_COLS, A1, A2)

ops.node(N_CONN, 0.0, 0.0)
ops.element('zeroLengthSection', ELE_CONN, N_FIX, N_CONN, SEC_GROUP,
            '-orient', 1, 0, 0)

ops.geomTransf('Linear', TRANSF)
ops.element('elasticBeamColumn', ELE_BEAM, N_CONN, N_TIP,
            A_BEAM, E_TIMBER, IZ_BEAM, TRANSF)

# Pushover: ramp a downward tip load to P_MAX over N_STEPS.
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(N_TIP, 0.0, -P_MAX, 0.0)
ops.wipeAnalysis()
ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('UmfPack')
ops.test('NormDispIncr', 1.0e-10, 50)
ops.algorithm('Newton')
ops.integrator('LoadControl', 1.0 / N_STEPS)
ops.analysis('Static')

p_flex, d_flex = [0.0], [0.0]
for i in range(N_STEPS):
    if ops.analyze(1) != 0:
        print("   stopped: no convergence at step %d of %d" % (i + 1, N_STEPS))
        break
    p_flex.append(P_MAX * (i + 1) / N_STEPS)
    d_flex.append(-ops.nodeDisp(N_TIP, 2))

# Decompose the current tip deflection into its three contributions. The
# identity delta = u_y + theta*L + flexure is exact for small rotations,
# which is the regime an elasticBeamColumn with a Linear transformation
# assumes anyway.
uy = -ops.nodeDisp(N_CONN, 2)
rot = -ops.nodeDisp(N_CONN, 3) * L
flex = d_flex[-1] - uy - rot
tip = d_flex[-1]

# (N, V, M) carried by the connection, in its local frame.
nvm = ops.eleResponse(ELE_CONN, 'section', 'forces')

# ---------------------------------------------------------------------
# Rigid base, same beam: the reference case that isolates how much of the
# tip deflection is the connection rather than the beam.
# ---------------------------------------------------------------------

ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

ops.node(N_FIX, 0.0, 0.0)
ops.fix(N_FIX, 1, 1, 1)
ops.node(N_TIP, L, 0.0)

ops.geomTransf('Linear', TRANSF)
ops.element('elasticBeamColumn', ELE_BEAM, N_FIX, N_TIP,
            A_BEAM, E_TIMBER, IZ_BEAM, TRANSF)

ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(N_TIP, 0.0, -P_MAX, 0.0)
ops.wipeAnalysis()
ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('UmfPack')
ops.test('NormDispIncr', 1.0e-10, 50)
ops.algorithm('Newton')
ops.integrator('LoadControl', 1.0 / N_STEPS)
ops.analysis('Static')

p_rig, d_rig = [0.0], [0.0]
for i in range(N_STEPS):
    if ops.analyze(1) != 0:
        print("   stopped: no convergence at step %d of %d" % (i + 1, N_STEPS))
        break
    p_rig.append(P_MAX * (i + 1) / N_STEPS)
    d_rig.append(-ops.nodeDisp(N_TIP, 2))

print("Cantilever on a %dx%d dowel-group connection "
      "(%d fasteners, L = %.2f m)\n" % (N_ROWS, N_COLS, N_ROWS * N_COLS, L))

print("At P = %.1f kN:" % p_flex[-1])
print("   tip deflection, flexible base : %8.3f mm" % (1000 * d_flex[-1]))
print("   tip deflection, rigid base    : %8.3f mm" % (1000 * d_rig[-1]))
print("   connection accounts for       : %8.1f %% of the total\n"
      % (100.0 * (d_flex[-1] - d_rig[-1]) / d_flex[-1]))

print("   deflection split (flexible base):")
print("      connection slip  u_y        : %8.3f mm  (%4.1f %%)"
      % (1000 * uy, 100 * uy / tip))
print("      connection rotation * L     : %8.3f mm  (%4.1f %%)"
      % (1000 * rot, 100 * rot / tip))
print("      beam flexure                : %8.3f mm  (%4.1f %%)"
      % (1000 * flex, 100 * flex / tip))
print("      total                       : %8.3f mm\n" % (1000 * tip))

print("   connection forces: N = %8.3f kN   V = %8.3f kN   M = %8.3f kN.m"
      % (nvm[0], nvm[1], nvm[2]))
print("   (M should equal P*L = %.3f kN.m)\n" % (p_flex[-1] * L))

# ---------------------------------------------------------------------
# Axial preload: N-M interaction at the connection. A tension preload eats
# into each fastener's capacity, so the same tip load costs more deflection
# -- the interaction the group section exists to capture, and which no
# uncoupled spring-per-DOF connection model can reproduce.
# ---------------------------------------------------------------------

print("Axial preload sweep (tip axial force held constant, then pushed):")
curves = {}
for n_val in (0.0, 60.0, 120.0):
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    ops.node(N_FIX, 0.0, 0.0)
    ops.fix(N_FIX, 1, 1, 1)
    ops.node(N_TIP, L, 0.0)

    R1 = K1_0 / K0_0
    FI = SAWS_FI_FRAC * F0_0
    ops.uniaxialMaterial('SAWS', MAT_DOWEL, F0_0, FI, SAWS_DU, K0_0,
                          R1, SAWS_R2, SAWS_R3, SAWS_R4, SAWS_ALPHA, SAWS_BETA)
    ops.section('HankinsonSection2d', SEC_BOLT, '-paramHankinson', MAT_DOWEL,
                '-param', 'F0', F0_0, F0_90,
                '-param', 'S0', K0_0, K0_90,
                '-param', 'K1', K1_0, K1_90)
    ops.section('BoltGroupSection2d', SEC_GROUP, '-section', SEC_BOLT,
                '-grid', N_ROWS, N_COLS, A1, A2)

    ops.node(N_CONN, 0.0, 0.0)
    ops.element('zeroLengthSection', ELE_CONN, N_FIX, N_CONN, SEC_GROUP,
                '-orient', 1, 0, 0)

    ops.geomTransf('Linear', TRANSF)
    ops.element('elasticBeamColumn', ELE_BEAM, N_CONN, N_TIP,
                A_BEAM, E_TIMBER, IZ_BEAM, TRANSF)

    if n_val:
        ops.timeSeries('Linear', 9)
        ops.pattern('Plain', 9, 9)
        ops.load(N_TIP, n_val, 0.0, 0.0)
        ops.wipeAnalysis()
        ops.constraints('Transformation')
        ops.numberer('RCM')
        ops.system('UmfPack')
        ops.test('NormDispIncr', 1.0e-10, 40)
        ops.algorithm('Newton')
        ops.integrator('LoadControl', 1.0 / 20)
        ops.analysis('Static')
        for _ in range(20):
            assert ops.analyze(1) == 0, "axial preload failed to converge"
        ops.loadConst('-time', 0.0)

    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(N_TIP, 0.0, -P_MAX, 0.0)
    # wipeAnalysis() first, otherwise a previously created analysis object is
    # silently kept and the settings below are ignored with a warning.
    ops.wipeAnalysis()
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('UmfPack')
    ops.test('NormDispIncr', 1.0e-10, 50)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0 / N_STEPS)
    ops.analysis('Static')

    p_h, d_h = [0.0], [0.0]
    for i in range(N_STEPS):
        if ops.analyze(1) != 0:
            print("   stopped: no convergence at step %d of %d" % (i + 1, N_STEPS))
            break
        p_h.append(P_MAX * (i + 1) / N_STEPS)
        d_h.append(-ops.nodeDisp(N_TIP, 2))
    curves[n_val] = (d_h, p_h)

    nvm_n = ops.eleResponse(ELE_CONN, 'section', 'forces')
    print("   N = %6.1f kN -> P %5.2f kN, tip %7.3f mm, "
          "connection (N, M) = (%7.3f kN, %7.3f kN.m)"
          % (n_val, p_h[-1], 1000 * d_h[-1], nvm_n[0], nvm_n[2]))

# ---------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("\nmatplotlib not available -- skipping plots.")
    raise SystemExit

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

ax1.plot([1000 * d for d in d_flex], p_flex, label='dowel-group base')
ax1.plot([1000 * d for d in d_rig], p_rig, '--', label='rigid base')
for n_val, (d_h, p_h) in curves.items():
    if n_val:
        ax1.plot([1000 * d for d in d_h], p_h, ':',
                 label='dowel-group base, N = %.0f kN' % n_val)
ax1.set_xlabel('Tip deflection (mm)')
ax1.set_ylabel('Tip load P (kN)')
ax1.set_title('Cantilever response: the connection dominates')
ax1.grid(True, linestyle=':')
ax1.legend(fontsize=8)

# Fastener force pattern at the end of the unpreloaded pushover: rebuild and
# resolve the plain flexible-base case, since the axial-preload sweep above
# left the domain in the N = 120 kN state.
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)

ops.node(N_FIX, 0.0, 0.0)
ops.fix(N_FIX, 1, 1, 1)
ops.node(N_TIP, L, 0.0)

R1 = K1_0 / K0_0
FI = SAWS_FI_FRAC * F0_0
ops.uniaxialMaterial('SAWS', MAT_DOWEL, F0_0, FI, SAWS_DU, K0_0,
                      R1, SAWS_R2, SAWS_R3, SAWS_R4, SAWS_ALPHA, SAWS_BETA)
ops.section('HankinsonSection2d', SEC_BOLT, '-paramHankinson', MAT_DOWEL,
            '-param', 'F0', F0_0, F0_90,
            '-param', 'S0', K0_0, K0_90,
            '-param', 'K1', K1_0, K1_90)
ops.section('BoltGroupSection2d', SEC_GROUP, '-section', SEC_BOLT,
            '-grid', N_ROWS, N_COLS, A1, A2)

ops.node(N_CONN, 0.0, 0.0)
ops.element('zeroLengthSection', ELE_CONN, N_FIX, N_CONN, SEC_GROUP,
            '-orient', 1, 0, 0)

ops.geomTransf('Linear', TRANSF)
ops.element('elasticBeamColumn', ELE_BEAM, N_CONN, N_TIP,
            A_BEAM, E_TIMBER, IZ_BEAM, TRANSF)

ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(N_TIP, 0.0, -P_MAX, 0.0)
ops.wipeAnalysis()
ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('UmfPack')
ops.test('NormDispIncr', 1.0e-10, 50)
ops.algorithm('Newton')
ops.integrator('LoadControl', 1.0 / N_STEPS)
ops.analysis('Static')
for i in range(N_STEPS):
    if ops.analyze(1) != 0:
        break

# The group's per-fastener force pattern, read from the group section's own
# 'boltForces' response (the concatenated (Fx, Fy) pairs -- note the plural:
# HankinsonSection2d intercepts the SINGULAR 'force' with its own scalar
# slip-law force, so only 'forces' gives the (P, Vy) pair).
f = ops.eleResponse(ELE_CONN, 'section', 'boltForces')
xs = [x for x, _ in coords]
ys = [y for _, y in coords]
fx = [f[2 * k] for k in range(len(coords))]
fy = [f[2 * k + 1] for k in range(len(coords))]

fmag = [(u * u + v * v) ** 0.5 for u, v in zip(fx, fy)]
fmax = max(fmag) if fmag else 0.0
span = max(max(xs) - min(xs), max(ys) - min(ys), 1.0e-6)
force_scale = 0.4 * span / fmax if fmax > 0 else 1.0

mx, my = ops.nodeCoord(N_CONN)
cr = ops.eleResponse(ELE_CONN, 'section', 'centerOfRotation')

ax2.scatter(xs, ys, s=70, c='saddlebrown', zorder=3, label='dowels')
ax2.scatter([mx], [my], marker='+', s=140, c='black', zorder=3,
            linewidths=2, label='reference point')
if cr and abs(cr[0]) < 10 * max(abs(v) for v in xs + ys):
    ax2.scatter([cr[0]], [cr[1]], marker='x', s=110, c='tab:blue',
                zorder=3, linewidths=2, label='centre of rotation')
for x, y in zip(xs, ys):
    ax2.plot([mx, x], [my, y], color='0.75', linewidth=0.8,
             linestyle='--', zorder=1)

ax2.quiver(xs, ys, [u * force_scale for u in fx],
           [v * force_scale for v in fy],
           angles='xy', scale_units='xy', scale=1,
           color='crimson', width=0.006, zorder=2,
           label='dowel force')

for j, ((x, y), m) in enumerate(zip(coords, fmag)):
    ax2.annotate("%d\n%.1f kN" % (j, m), (x, y),
                 textcoords="offset points", xytext=(7, 7), fontsize=7)

ax2.set_aspect('equal', adjustable='datalim')
ax2.set_xlabel('x (m, parallel to grain)')
ax2.set_ylabel('y (m, perpendicular to grain)')
ax2.grid(True, linestyle=':')
ax2.legend(loc='best', fontsize=8)
ax2.set_title('Fastener forces at P = %.0f kN' % P_MAX)

fig.suptitle('Timber cantilever on a %dx%d dowel-group connection '
             '(BoltGroupSection2d, one element)' % (N_ROWS, N_COLS))
fig.tight_layout()
plt.show()
