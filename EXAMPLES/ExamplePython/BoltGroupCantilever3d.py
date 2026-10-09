# Worked example: the same dowel-group cantilever as BoltGroupCantilever.py,
# but in a full 3D model (ndm=3, ndf=6) -- which is where the one real
# limitation of BoltGroupSection2d has to be dealt with head-on.
#
# BoltGroupSection2d reports codes {P, Vy, Mz}: it resists the two in-plane
# translations and the in-plane rotation of the connection, and NOTHING else.
# That is physically right for dowel-type fasteners, which transfer load by
# bearing in their own shear plane and have no meaningful resistance to being
# pulled or twisted out of it. But a 3D cantilever hanging off such a
# connection is then free to translate in Z and rotate about X and Y with no
# stiffness anywhere in the model -- the stiffness matrix is singular, and the
# analysis fails to factorize rather than silently doing something wrong.
#
# There are two honest ways to deal with that, and this file runs BOTH and
# shows they give the same in-plane answer:
#
#   'aggregate'  Wrap the group in a SectionAggregator that supplies elastic
#                Vz / My / T stiffness. This is the better option when
#                something real does restrain the connection out of plane
#                (the timber member's own continuity, a steel plate, other
#                fasteners) and you want to represent it with a finite
#                stiffness. SectionAggregator is block-diagonal, so the
#                added terms cannot contaminate the group's own {P, Vy, Mz}.
#
#   'fix'        Simply restrain Uz / Rx / Ry at the connection node. This is
#                the option HankinsonDowelGroup3d.py uses; it is the right
#                one when the out-of-plane restraint is effectively rigid.
#
# The model is otherwise identical to the 2D example, and the in-plane
# response must (and does, see the cross-check the demo prints) reproduce it
# exactly: the group plane is placed by the element's '-orient', and the
# fastener offsets live in that local frame.
#
# Units: kN, m.
#
# Run standalone:
#   cd EXAMPLES/ExamplePython && cp ../../build/OpenSeesPy.so ./opensees.so
#   python BoltGroupCantilever3d.py

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

# ---------------------------------------------------------------------
# Geometry and properties -- same beam and connection as BoltGroupCantilever.py
# ---------------------------------------------------------------------

L = 1.0                      # cantilever span (m)
B, H = 0.200, 0.320          # beam width and depth (m)
E_TIMBER = 11.5e6            # kN/m2 (11500 MPa, glulam order of magnitude)

N_ROWS, N_COLS = 5, 3
A1, A2 = 0.070, 0.060

A_BEAM = B * H
IZ_BEAM = B * H ** 3 / 12.0        # in-plane bending (about local z)
IY_BEAM = H * B ** 3 / 12.0        # out-of-plane bending (about local y)
J_BEAM = IY_BEAM + IZ_BEAM         # crude polar approximation, fine here
G_TIMBER = E_TIMBER / (2.0 * (1.0 + 0.3))

# Out-of-plane stand-in stiffnesses for the 'aggregate' option. These are NOT
# fastener behaviour -- they represent whatever else restrains the connection
# out of its own plane, and are deliberately stiff enough that the
# out-of-plane response is governed by the beam rather than by them.
K_VZ, K_MY, K_T = 1.0e5, 1.0e4, 1.0e4

# The fastener law: the Schweigler dowel calibration, Hankinson-interpolated
# between a 0-deg (parallel-to-grain) and 90-deg (perpendicular-to-grain)
# Foschi envelope by each fastener's own evolving load-to-grain angle,
# driven through a SAWSMaterial (Folz & Filiatrault CUREE pinching
# hysteretic model) backbone.
F0_0, K0_0, K1_0 = 19.1, 24460.0, 462.0
F0_90, K0_90, K1_90 = 14.18, 19640.0, 739.0
SAWS_DU, SAWS_FI_FRAC = 0.05, 0.10
SAWS_R2, SAWS_R3, SAWS_R4, SAWS_ALPHA, SAWS_BETA = -0.05, 1.0, 0.0, 0.8, 1.05

MAT_DOWEL, SEC_BOLT, SEC_GROUP, SEC_AGG = 20, 10, 11, 12
MAT_VZ, MAT_MY, MAT_T = 31, 32, 33
TRANSF = 1
N_FIX, N_CONN, N_TIP = 1, 2, 3
ELE_CONN, ELE_BEAM = 1, 2

P_MAX, N_STEPS = 30.0, 120

# Fastener (x, y) offsets referred to the group's geometric centre, for
# labelling output.
raw_coords = [(A1 * col, A2 * row) for row in range(N_ROWS) for col in range(N_COLS)]
x_cg = (0.5 * N_COLS - 0.5) * A1
y_cg = (0.5 * N_ROWS - 0.5) * A2
coords = [(x - x_cg, y - y_cg) for x, y in raw_coords]

# ---------------------------------------------------------------------
# In-plane pushover under both out-of-plane treatments: the two must agree.
# ---------------------------------------------------------------------

print("3D cantilever on a %dx%d dowel-group connection "
      "(%d fasteners, L = %.2f m)\n" % (N_ROWS, N_COLS, N_ROWS * N_COLS, L))

results = {}
for restrain in ('aggregate', 'fix'):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)

    ops.node(N_FIX, 0.0, 0.0, 0.0)
    ops.fix(N_FIX, 1, 1, 1, 1, 1, 1)
    ops.node(N_CONN, 0.0, 0.0, 0.0)
    ops.node(N_TIP, L, 0.0, 0.0)

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

    if restrain == 'aggregate':
        ops.uniaxialMaterial('Elastic', MAT_VZ, K_VZ)
        ops.uniaxialMaterial('Elastic', MAT_MY, K_MY)
        ops.uniaxialMaterial('Elastic', MAT_T, K_T)
        # Aggregator order: the added uniaxial DOFs first, then '-section'
        # splices in the group's own {P, Vy, Mz}. The result is order 6.
        ops.section('Aggregator', SEC_AGG,
                    MAT_VZ, 'Vz', MAT_MY, 'My', MAT_T, 'T',
                    '-section', SEC_GROUP)
        sec_tag = SEC_AGG
    else:
        sec_tag = SEC_GROUP

    # '-orient' takes BOTH local axes in 3D: local x along the beam and the
    # grain, local y through the beam depth. That -- not anything inside the
    # section -- is what places the group plane.
    ops.element('zeroLengthSection', ELE_CONN, N_FIX, N_CONN, sec_tag,
                '-orient', 1, 0, 0, 0, 1, 0)

    if restrain == 'fix':
        # Uz, Rx, Ry get no stiffness from the bare group; restrain them at
        # the connection node (the beam carries them onward to the tip).
        ops.fix(N_CONN, 0, 0, 1, 1, 1, 0)

    # vecxz = (0, 0, 1) puts local y along global Y, so local-z bending (Iz)
    # is the in-plane bending that matches the 2D example.
    ops.geomTransf('Linear', TRANSF, 0, 0, 1)
    ops.element('elasticBeamColumn', ELE_BEAM, N_CONN, N_TIP,
                A_BEAM, E_TIMBER, G_TIMBER, J_BEAM, IY_BEAM, IZ_BEAM, TRANSF)

    # In-plane pushover: a tip load along -Y, exercising the group's Vy/Mz.
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(N_TIP, 0.0, -P_MAX, 0.0, 0.0, 0.0, 0.0)
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

    nvm = ops.eleResponse(ELE_CONN, 'section', 'forces')
    results[restrain] = (p_h, d_h, nvm)
    print("   restrain='%s': tip %8.3f mm at P = %.1f kN" %
          (restrain, 1000 * d_h[-1], p_h[-1]))

d_agg, d_fix = results['aggregate'][1][-1], results['fix'][1][-1]
print("   difference between the two: %.3e mm\n"
      % (1000 * abs(d_agg - d_fix)))

# ---------------------------------------------------------------------
# Cross-check against the 2D model (BoltGroupCantilever.py): both must
# reproduce it exactly.
# ---------------------------------------------------------------------

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

d2d = [0.0]
for i in range(N_STEPS):
    if ops.analyze(1) != 0:
        break
    d2d.append(-ops.nodeDisp(N_TIP, 2))

print("cross-check against the 2D model (BoltGroupCantilever.py):")
print("   2D tip deflection : %.12e m" % d2d[-1])
print("   3D tip deflection : %.12e m" % d_agg)
print("   difference        : %.3e m\n" % abs(d2d[-1] - d_agg))

# ---------------------------------------------------------------------
# Out-of-plane load case: the group is inert, the aggregate carries it.
# ---------------------------------------------------------------------

ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 6)

ops.node(N_FIX, 0.0, 0.0, 0.0)
ops.fix(N_FIX, 1, 1, 1, 1, 1, 1)
ops.node(N_CONN, 0.0, 0.0, 0.0)
ops.node(N_TIP, L, 0.0, 0.0)

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

ops.uniaxialMaterial('Elastic', MAT_VZ, K_VZ)
ops.uniaxialMaterial('Elastic', MAT_MY, K_MY)
ops.uniaxialMaterial('Elastic', MAT_T, K_T)
ops.section('Aggregator', SEC_AGG,
            MAT_VZ, 'Vz', MAT_MY, 'My', MAT_T, 'T',
            '-section', SEC_GROUP)

ops.element('zeroLengthSection', ELE_CONN, N_FIX, N_CONN, SEC_AGG,
            '-orient', 1, 0, 0, 0, 1, 0)

ops.geomTransf('Linear', TRANSF, 0, 0, 1)
ops.element('elasticBeamColumn', ELE_BEAM, N_CONN, N_TIP,
            A_BEAM, E_TIMBER, G_TIMBER, J_BEAM, IY_BEAM, IZ_BEAM, TRANSF)

# Out-of-plane pushover: a tip load along -Z, which the group cannot resist
# at all -- everything here is carried by the aggregated Vz/My/T stiffness.
P_OUT = 5.0
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(N_TIP, 0.0, 0.0, -P_OUT, 0.0, 0.0, 0.0)
ops.wipeAnalysis()
ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('UmfPack')
ops.test('NormDispIncr', 1.0e-10, 50)
ops.algorithm('Newton')
ops.integrator('LoadControl', 1.0 / N_STEPS)
ops.analysis('Static')

p_out, d_out = [0.0], [0.0]
for i in range(N_STEPS):
    if ops.analyze(1) != 0:
        print("   stopped: no convergence at step %d of %d" % (i + 1, N_STEPS))
        break
    p_out.append(P_OUT * (i + 1) / N_STEPS)
    d_out.append(-ops.nodeDisp(N_TIP, 3))

# The group's per-fastener forces. This reaches the group even through the
# SectionAggregator: the aggregator does not recognise 'boltForces', and its
# setResponse falls through to the wrapped section.
f_out = ops.eleResponse(ELE_CONN, 'section', 'boltForces')

print("Out-of-plane load case (the group resists nothing here):")
print("   tip deflection out of plane : %8.3f mm at P = %.1f kN"
      % (1000 * d_out[-1], p_out[-1]))
print("   largest fastener force       : %.3e kN  (exactly zero: the "
      "group\n                                  sees no deformation in "
      "Vz/My/T)" % max(abs(v) for v in f_out))

# ---------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("\nmatplotlib not available -- skipping plots.")
    raise SystemExit

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

for restrain, style in (('aggregate', '-'), ('fix', '--')):
    p_h, d_h, _ = results[restrain]
    ax1.plot([1000 * d for d in d_h], p_h, style,
             label="in-plane, restrain='%s'" % restrain)
ax1.plot([1000 * d for d in d_out], p_out, ':',
         label='out-of-plane (aggregated stiffness only)')
ax1.set_xlabel('Tip deflection (mm)')
ax1.set_ylabel('Tip load P (kN)')
ax1.set_title('In-plane vs out-of-plane response')
ax1.grid(True, linestyle=':')
ax1.legend(fontsize=8)

# Fastener force pattern at the end of the in-plane pushover, 'aggregate'
# treatment: rebuild and resolve, since the out-of-plane case above left the
# domain in the out-of-plane-loaded state.
ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 6)

ops.node(N_FIX, 0.0, 0.0, 0.0)
ops.fix(N_FIX, 1, 1, 1, 1, 1, 1)
ops.node(N_CONN, 0.0, 0.0, 0.0)
ops.node(N_TIP, L, 0.0, 0.0)

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

ops.uniaxialMaterial('Elastic', MAT_VZ, K_VZ)
ops.uniaxialMaterial('Elastic', MAT_MY, K_MY)
ops.uniaxialMaterial('Elastic', MAT_T, K_T)
ops.section('Aggregator', SEC_AGG,
            MAT_VZ, 'Vz', MAT_MY, 'My', MAT_T, 'T',
            '-section', SEC_GROUP)

ops.element('zeroLengthSection', ELE_CONN, N_FIX, N_CONN, SEC_AGG,
            '-orient', 1, 0, 0, 0, 1, 0)

ops.geomTransf('Linear', TRANSF, 0, 0, 1)
ops.element('elasticBeamColumn', ELE_BEAM, N_CONN, N_TIP,
            A_BEAM, E_TIMBER, G_TIMBER, J_BEAM, IY_BEAM, IZ_BEAM, TRANSF)

ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(N_TIP, 0.0, -P_MAX, 0.0, 0.0, 0.0, 0.0)
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

f = ops.eleResponse(ELE_CONN, 'section', 'boltForces')
xs = [x for x, _ in coords]
ys = [y for _, y in coords]
fx = [f[2 * k] for k in range(len(coords))]
fy = [f[2 * k + 1] for k in range(len(coords))]

fmag = [(u * u + v * v) ** 0.5 for u, v in zip(fx, fy)]
fmax = max(fmag) if fmag else 0.0
span = max(max(xs) - min(xs), max(ys) - min(ys), 1.0e-6)
force_scale = 0.4 * span / fmax if fmax > 0 else 1.0

mx, my, _ = ops.nodeCoord(N_CONN)
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
ax2.set_title('Fastener forces at P = %.0f kN (in-plane)' % P_MAX)

fig.suptitle('Timber cantilever on a %dx%d dowel-group connection, 3D '
             'model\n(BoltGroupSection2d + SectionAggregator for the '
             'out-of-plane DOFs)' % (N_ROWS, N_COLS))
fig.tight_layout()
plt.show()
