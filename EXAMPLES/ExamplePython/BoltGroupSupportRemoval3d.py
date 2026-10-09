# Worked example: SUDDEN SUPPORT REMOVAL on a timber beam whose base
# connection is a 3x2 dowel group modelled as one BoltGroupSection2d element.
#
# This is the standard alternate-load-path / progressive-collapse procedure,
# and it is here for two reasons: it is the load case in which a semi-rigid
# timber connection matters most (the connection has to absorb, dynamically,
# a moment it never saw in service), and it exercises BoltGroupSection2d in
# a TRANSIENT analysis -- inertia, damping and path-dependent unloading --
# rather than the monotonic pushovers of BoltGroupCantilever{,3d}.py.
#
# The model is a propped cantilever: the dowel group at one end, a vertical
# prop at the other, a uniformly distributed gravity load in between.
#
#     node 1 ==[ BoltGroupSection2d ]== 2 --[beam]--[beam]--...--[beam]-- 8
#      (fixed)     zeroLengthSection                                     /_\
#                                       <------ w, uniformly distributed -->
#
# The four phases:
#
#   1. STATIC     Ramp w on with the prop in place, then read the prop's
#                 reaction R.
#   2. FREEZE     loadConst() makes w a constant load at pseudo-time 0, so it
#                 keeps acting for the whole transient.
#   3. RELEASE    remove('sp', tip, 2) deletes the prop, and the SAME instant
#                 an external load R is applied in its place. The structure
#                 cannot tell the difference: it is in exactly the equilibrium
#                 it was already in, now at rest with no support. A Path time
#                 series then ramps R to zero over T_RAMP, which is the
#                 "removal". T_HOLD of quiet time before the ramp is the
#                 check: if R had the wrong sign or the wrong magnitude, the
#                 beam would move during it. The demo measures that drift.
#   4. TRANSIENT  Newmark, ordinary implicit time stepping, ~7 cycles.
#
# The one thing that must be got right for this to mean anything: THE MASS
# MUST BE THE MASS OF THE LOAD. The distributed load w is the weight of a
# tributary mass, so the beam carries rho = w/g per unit length. If you load
# a massless beam and then release it, the released load has nothing to
# accelerate and the "dynamic" answer is meaningless (and the analysis is
# ill-conditioned besides). Here RHO is derived from W_LINE, not chosen.
#
# What comes out is the dynamic amplification the connection actually sees.
# The demo brackets it by running the identical release QUASI-STATICALLY
# (same load path, no inertia), so the peak dynamic response can be reported
# against the true static answer for the released structure rather than
# against the textbook factor of 2. That factor does not survive contact
# with a real connection, and the two amplifications the demo prints
# straddle it in opposite directions for a reason worth knowing:
#
#   * on the connection MOMENT the factor comes out below 2, because the
#     group yields and cannot deliver a moment proportional to the
#     overshoot, and because damping bleeds off the first peak;
#   * on DEFLECTION it comes out above 2 anyway, because that same yielding
#     means the extra moment is bought at a much lower tangent stiffness.
#
# Reporting one and assuming the other is how connections get under-designed
# for this load case.
#
# Units: kN, m, s -- so mass is in kN*s^2/m = tonne, and g = 9.81 m/s^2.
#
# Run standalone:
#   cd EXAMPLES/ExamplePython && cp ../../build/OpenSeesPy.so ./opensees.so
#   python BoltGroupSupportRemoval3d.py

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

import math

# ---------------------------------------------------------------------
# Geometry, properties, load
# ---------------------------------------------------------------------

L = 0.50                      # span, prop to connection (m)
N_ELE = 6                     # beam elements along the span

B, H = 0.200, 0.320           # glulam section (m)
E_TIMBER = 11.5e6             # kN/m2 (11500 MPa, glulam order of magnitude)
G_TIMBER = E_TIMBER / (2.0 * (1.0 + 0.3))

A_BEAM = B * H
IZ_BEAM = B * H ** 3 / 12.0   # in-plane bending, about local z
IY_BEAM = H * B ** 3 / 12.0
J_BEAM = IY_BEAM + IZ_BEAM    # crude polar approximation, adequate here

# The fastener law: the Schweigler dowel calibration, Hankinson-interpolated
# between a 0-deg (parallel-to-grain) and 90-deg (perpendicular-to-grain)
# Foschi envelope by each fastener's own evolving load-to-grain angle,
# driven through a SAWSMaterial (Folz & Filiatrault CUREE pinching
# hysteretic model) backbone.
F0_0, K0_0, K1_0 = 19.1, 24460.0, 462.0
F0_90, K0_90, K1_90 = 14.18, 19640.0, 739.0
SAWS_DU, SAWS_FI_FRAC = 0.05, 0.10
SAWS_R2, SAWS_R3, SAWS_R4, SAWS_ALPHA, SAWS_BETA = -0.05, 1.0, 0.0, 0.8, 1.05

# Connection: 3 rows x 2 columns. A1 is the column pitch (along the beam
# axis and the grain), A2 the row pitch through the depth, so the fasteners
# sit at x = +-0.035 and y = -0.10, 0, +0.10 -- spread through the 0.32 m
# depth with a 60 mm edge distance, which is what gives the group its
# moment arm.
N_ROWS, N_COLS = 3, 3
A1, A2 = 0.070, 0.100

W_LINE = 30.0                 # distributed gravity load (kN/m)
GRAV = 9.81                   # m/s2
RHO = W_LINE / GRAV           # tonne/m -- NOT independent of W_LINE

ZETA = 0.03                   # modal damping ratio, applied to N_MODES modes
N_MODES = 4

DT = 2.0e-4                   # transient time step (s)
T_HOLD = 0.010                # quiet time before the release (s)
T_RAMP = 0.001                # duration of the removal itself (s)
T_END = 1.00                  # end of the transient (s), ~13 cycles

MAT_DOWEL, SEC_BOLT, SEC_GROUP = 20, 10, 11
TRANSF = 1
N_FIX, N_CONN = 1, 2
N_TIP = N_CONN + N_ELE        # the propped end
ELE_CONN = 1
BEAM_ELES = [101 + i for i in range(N_ELE)]

PAT_GRAV, TS_GRAV = 1, 1      # distributed load
PAT_REL, TS_REL = 2, 2        # the released prop reaction

DOF_PROP = 2                  # the propped DOF: global Y (in the group plane)

# Fastener (x, y) offsets referred to the group's geometric centre.
raw_coords = [(A1 * col, A2 * row) for row in range(N_ROWS) for col in range(N_COLS)]
x_cg = (0.5 * N_COLS - 0.5) * A1
y_cg = (0.5 * N_ROWS - 0.5) * A2
coords = [(x - x_cg, y - y_cg) for x, y in raw_coords]

print("Sudden support removal, %dx%d dowel-group connection "
      "(%d fasteners)" % (N_ROWS, N_COLS, N_ROWS * N_COLS))
print("   span %.2f m, w = %.1f kN/m (total %.1f kN), "
      "tributary mass %.4f t\n" % (L, W_LINE, W_LINE * L, RHO * L))

# ---------------------------------------------------------------------
# 1. static, propped
# ---------------------------------------------------------------------

ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 6)

ops.node(N_FIX, 0.0, 0.0, 0.0)
ops.fix(N_FIX, 1, 1, 1, 1, 1, 1)
for i in range(N_ELE + 1):
    ops.node(N_CONN + i, L * i / N_ELE, 0.0, 0.0)

# One fastener: the Schweigler dowel law, its F0/S0/K1 Hankinson-
# interpolated by that fastener's own evolving load-to-grain angle.
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

# '-orient' places the group plane: local x along the beam and the grain,
# local y through the depth. The fastener offsets live in that frame.
ops.element('zeroLengthSection', ELE_CONN, N_FIX, N_CONN, SEC_GROUP,
            '-orient', 1, 0, 0, 0, 1, 0)

# The group reports {P, Vy, Mz} and nothing else, so Uz/Rx/Ry get no
# stiffness from it. Restrain them at the connection node -- the
# 'aggregate' alternative is shown in BoltGroupCantilever3d.py; here the
# rigid choice is the right one and it keeps the transient's mass matrix
# from having free, massless, unrestrained DOFs.
ops.fix(N_CONN, 0, 0, 1, 1, 1, 0)

# vecxz = (0, 0, 1) puts local y along global Y, so local-z bending is the
# in-plane bending the group resists.
ops.geomTransf('Linear', TRANSF, 0, 0, 1)
for i, etag in enumerate(BEAM_ELES):
    ops.element('elasticBeamColumn', etag, N_CONN + i, N_CONN + i + 1,
                A_BEAM, E_TIMBER, G_TIMBER, J_BEAM, IY_BEAM, IZ_BEAM,
                TRANSF, '-mass', RHO)

# The prop: a vertical roller at the far end. This is the support that gets
# removed.
ops.fix(N_TIP, 0, 1, 0, 0, 0, 0)

# Ramp the distributed gravity load on with the prop in place.
ops.timeSeries('Linear', TS_GRAV)
ops.pattern('Plain', PAT_GRAV, TS_GRAV)
# 3D -beamUniform takes (Wy, Wz) in local axes; local y is global Y.
ops.eleLoad('-ele', *BEAM_ELES, '-type', '-beamUniform', -W_LINE, 0.0)

ops.wipeAnalysis()
ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('UmfPack')
ops.test('NormDispIncr', 1.0e-10, 50)
ops.algorithm('Newton')
ops.integrator('LoadControl', 1.0 / 20)
ops.analysis('Static')
for i in range(20):
    if ops.analyze(1) != 0:
        raise RuntimeError("static phase failed at step %d" % (i + 1))

ops.reactions()
reaction = ops.nodeReaction(N_TIP, DOF_PROP)

ops.reactions()
r_base = ops.nodeReaction(N_FIX, DOF_PROP)
nvm_propped = ops.eleResponse(ELE_CONN, 'section', 'forces')
# Identically zero -- that is what the prop does. It is kept as the
# reference for the hold-phase drift check below, which is therefore a
# check against exact zero.
d_propped = -ops.nodeDisp(N_TIP, DOF_PROP)

print("Phase 1 -- static, prop in place:")
print("   prop reaction        : %8.3f kN" % reaction)
print("   base reaction        : %8.3f kN" % r_base)
print("   sum vs applied load  : %8.3f kN vs %.3f kN  (residual %.2e)"
      % (reaction + r_base, W_LINE * L,
         abs(reaction + r_base - W_LINE * L)))
print("   connection moment    : %8.3f kN.m" % nvm_propped[2])
print("   connection rotation  : %8.3e rad" % ops.nodeDisp(N_CONN, 6))
print("   midspan deflection   : %8.4f mm"
      % (-1000 * ops.nodeDisp(N_CONN + N_ELE // 2, DOF_PROP)))
print("   (the propped end does not move: %.1e mm)\n"
      % (1000 * d_propped))

# ---------------------------------------------------------------------
# 2/3. freeze, release, and damp
# ---------------------------------------------------------------------

# Freeze the gravity load, delete the prop, put its reaction back as an
# external load, and schedule that load to vanish over T_RAMP seconds.
#
# On the sign: nodeReaction() returns the force the SUPPORT applies to the
# structure (upward, positive, under downward load). Replacing the support
# by that same force leaves the structure in identical equilibrium, so it
# is applied with no sign change. The hold phase is what proves this: any
# error in sign or magnitude shows up immediately as drift, and the demo
# reports it.
#
# Constant, at pseudo-time 0: the gravity load stays on for the whole
# transient rather than scaling with time as a Linear series would.
ops.loadConst('-time', 0.0)
ops.remove('sp', N_TIP, DOF_PROP)
ops.timeSeries('Path', TS_REL,
               '-time', [0.0, T_HOLD, T_HOLD + T_RAMP, T_END],
               '-values', [1.0, 1.0, 0.0, 0.0])
ops.pattern('Plain', PAT_REL, TS_REL)
ops.load(N_TIP, 0.0, reaction, 0.0, 0.0, 0.0, 0.0)

# Incremental modal damping: a flat ratio ZETA on the first N_MODES modes of
# the RELEASED structure.
#
# Modal rather than Rayleigh, because Rayleigh damping is anchored to a
# stiffness this structure does not keep. The group softens by roughly a
# factor of two as it yields, so a beta*K term calibrated on the initial
# tangent fades away exactly when the response is largest, while the
# alpha*M term over-damps whatever is left. Modal damping fixes zeta per
# mode and is indifferent to that.
#
# 'Incremental' is ops.modalDamping(): the modal damping matrix
# 2*zeta_n*w_n*(M phi_n)(M phi_n)^T is assembled into the TANGENT every
# increment (TransientIntegrator::formTangent calls addModalDampingMatrix)
# as well as contributing its force to the residual, so Newton's Jacobian
# is consistent with the damping it is being asked to balance. The
# alternative, ops.modalDampingQ(), adds only the force, leaving the
# tangent ignorant of it; that costs iterations at each reversal and is
# what you fall back to if forming the modal matrix is too expensive.
#
# Two things this buys and one it does not: it is exact for the modes it
# covers and adds nothing outside them, but the modes are frozen at the
# state where eigen() was called. The connection yields after that, so the
# true modes drift; damping stays anchored to the released-but-elastic
# shapes. That is the usual and accepted approximation, and the reason
# eigen() is called AFTER the static phase rather than on the pristine
# model.
#
# Drop the static analysis first: eigen() builds its own, and reusing a
# UmfPack SOE that was sized for the propped structure fails outright once
# the prop has been removed.
ops.wipeAnalysis()
lam = ops.eigen('-genBandArpack', N_MODES)
ops.modalDamping(ZETA)
# The eigen data lives on the domain and the nodes, not on the analysis, so
# the wipeAnalysis() that sets up the transient below does not lose it --
# but modalDamping() itself must be called while the eigen SOE is still
# alive, i.e. here and not after.
w = [math.sqrt(v) for v in lam]
w1 = w[0]

print("Phase 2/3 -- prop removed, reaction re-applied then ramped out:")
print("   released periods     : %s s"
      % ", ".join("%.4f" % (2 * math.pi / wi) for wi in w))
print("   damping              : incremental modal, %.1f %% on %d modes"
      % (100 * ZETA, N_MODES))
print("   time step            : %.1e s  (T1/%.0f)"
      % (DT, (2 * math.pi / w1) / DT))
print("   removal ramp         : %.1e s  (T1/%.1f)\n"
      % (T_RAMP, (2 * math.pi / w1) / T_RAMP))

# ---------------------------------------------------------------------
# 4. transient -- Newmark average acceleration to T_END.
# ---------------------------------------------------------------------

# Beside the time series plotted below, the sampling records what the two
# families of hysteresis loops need:
#
#   'th', 'm'    connection rotation and moment  -- the GROUP's loop
#   'uy', 'v'    connection slip and shear       -- likewise
#   'bs', 'bf'   per-fastener signed slip and force along that fastener's
#                own slip direction -- the loops of the individual dowel
#                laws that the group's loop is the sum of
#
# The per-fastener pair comes from HankinsonSection2d's 'state' response
# (alpha, delta, slip, sign, force, k, secant) rather than from the group's
# 'boltForces'/'boltSlips' vectors, and deliberately so: 'state' reports the
# SCALAR slip-law quantities -- the actual 1-D hysteretic law being
# exercised, signed along the current slip direction -- whereas the group's
# vectors are the 2-D components, whose x-y decomposition rotates as the
# load-to-grain angle evolves and does not draw a clean loop.
ops.wipeAnalysis()
ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('UmfPack')
ops.test('NormDispIncr', 1.0e-10, 50)
ops.algorithm('Newton')
ops.integrator('Newmark', 0.5, 0.25)
ops.analysis('Transient')

nb = N_ROWS * N_COLS
h = {'t': [], 'd': [], 'm': [], 'v': [], 'th': [], 'uy': [], 'fmax': [],
     'bs': [[] for _ in range(nb)], 'bf': [[] for _ in range(nb)]}

h['t'].append(ops.getTime())
h['d'].append(-ops.nodeDisp(N_TIP, DOF_PROP))
q = ops.eleResponse(ELE_CONN, 'section', 'forces')
e = ops.eleResponse(ELE_CONN, 'section', 'deformations')
h['v'].append(q[1])
h['m'].append(q[2])
h['uy'].append(e[1])
h['th'].append(e[2])
f = ops.eleResponse(ELE_CONN, 'section', 'boltForces')
h['fmax'].append(max(math.hypot(f[2 * i], f[2 * i + 1]) for i in range(nb)))
for i in range(nb):
    st = ops.eleResponse(ELE_CONN, 'section', 'bolt', i, 'state')
    h['bs'][i].append(st[2])
    h['bf'][i].append(st[4])

n_steps = int(round(T_END / DT))
for i in range(n_steps):
    if ops.analyze(1, DT) != 0:
        # Ordinary sub-stepping fallback: the release itself is the only
        # place this model has ever needed it.
        ok = -1
        for _ in range(4):
            ok = ops.analyze(2, DT / 2.0)
            if ok == 0:
                break
        if ok != 0:
            print("   transient stopped at t = %.4f s (step %d of %d)"
                  % (ops.getTime(), i + 1, n_steps))
            break

    h['t'].append(ops.getTime())
    h['d'].append(-ops.nodeDisp(N_TIP, DOF_PROP))
    q = ops.eleResponse(ELE_CONN, 'section', 'forces')
    e = ops.eleResponse(ELE_CONN, 'section', 'deformations')
    h['v'].append(q[1])
    h['m'].append(q[2])
    h['uy'].append(e[1])
    h['th'].append(e[2])
    f = ops.eleResponse(ELE_CONN, 'section', 'boltForces')
    h['fmax'].append(max(math.hypot(f[2 * i], f[2 * i + 1]) for i in range(nb)))
    for i in range(nb):
        st = ops.eleResponse(ELE_CONN, 'section', 'bolt', i, 'state')
        h['bs'][i].append(st[2])
        h['bf'][i].append(st[4])

t, d, m, f = h['t'], h['d'], h['m'], h['fmax']

# The hold phase is the equilibrium check: with the prop replaced by its own
# reaction the structure must not move at all before the ramp starts.
hold_idx = [i for i, ti in enumerate(t) if ti <= T_HOLD]
drift = max(abs(d[i] - d_propped) for i in hold_idx)

i_peak = max(range(len(d)), key=lambda i: d[i])
i_mpeak = max(range(len(m)), key=lambda i: abs(m[i]))

# The beam is still ringing at t_end, so the settled value is the mean of
# the last cycle rather than the last sample.
tail = [i for i in range(len(t)) if t[i] >= t[-1] - 2 * math.pi / w1]
d_settled = sum(d[i] for i in tail) / len(tail)
m_settled = sum(m[i] for i in tail) / len(tail)

# Modal damping is easy to ask for and get nothing (it silently does
# nothing if eigen() was not called, or if the eigen data was lost), so
# measure it back out of the response: the logarithmic decrement of the
# free vibration about the settled position.
peaks = [i for i in range(1, len(d) - 1)
         if d[i] > d[i - 1] and d[i] >= d[i + 1] and t[i] > T_HOLD + T_RAMP]
amp = [d[i] - d_settled for i in peaks]
n_dec = min(6, len(amp) - 1)
zeta_eff = float('nan')
if n_dec > 0 and amp[0] > 0 and amp[n_dec] > 0:
    delta = math.log(amp[0] / amp[n_dec]) / n_dec
    zeta_eff = delta / math.sqrt(4 * math.pi ** 2 + delta ** 2)

print("Phase 4 -- transient (%d steps to t = %.2f s):" % (len(t) - 1, t[-1]))
print("   drift during the hold : %.3e mm  "
      "(the equilibrium check: the\n                           structure "
      "must not move before the ramp)" % (1000 * drift))
# Close to ZETA is the check that modal damping took effect at all; it is
# not expected to land on it. The log decrement is a linear, constant-mean
# measure being applied to a response that is neither: the group dissipates
# hysteretically on top of the modal damping, and its mean position is
# still migrating while it yields, which biases amplitudes measured about
# the final one. Undamped, this would read 0.
print("   damping, measured back: %.2f %% over %d cycles  "
      "(asked for %.1f %% -- see\n                           the note in "
      "the source for why it does not land exactly)"
      % (100 * zeta_eff, n_dec, 100 * ZETA))

# ---------------------------------------------------------------------
# The quasi-static bracket
# ---------------------------------------------------------------------

ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 6)

ops.node(N_FIX, 0.0, 0.0, 0.0)
ops.fix(N_FIX, 1, 1, 1, 1, 1, 1)
for i in range(N_ELE + 1):
    ops.node(N_CONN + i, L * i / N_ELE, 0.0, 0.0)

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
ops.element('zeroLengthSection', ELE_CONN, N_FIX, N_CONN, SEC_GROUP,
            '-orient', 1, 0, 0, 0, 1, 0)
ops.fix(N_CONN, 0, 0, 1, 1, 1, 0)
ops.geomTransf('Linear', TRANSF, 0, 0, 1)
for i, etag in enumerate(BEAM_ELES):
    ops.element('elasticBeamColumn', etag, N_CONN + i, N_CONN + i + 1,
                A_BEAM, E_TIMBER, G_TIMBER, J_BEAM, IY_BEAM, IZ_BEAM,
                TRANSF, '-mass', RHO)
ops.fix(N_TIP, 0, 1, 0, 0, 0, 0)

ops.timeSeries('Linear', TS_GRAV)
ops.pattern('Plain', PAT_GRAV, TS_GRAV)
ops.eleLoad('-ele', *BEAM_ELES, '-type', '-beamUniform', -W_LINE, 0.0)
ops.wipeAnalysis()
ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('UmfPack')
ops.test('NormDispIncr', 1.0e-10, 50)
ops.algorithm('Newton')
ops.integrator('LoadControl', 1.0 / 20)
ops.analysis('Static')
for i in range(20):
    if ops.analyze(1) != 0:
        raise RuntimeError("static phase failed at step %d" % (i + 1))
ops.reactions()
reaction_qs = ops.nodeReaction(N_TIP, DOF_PROP)

# The identical release with the inertia left out: ramp the reaction to zero
# under LoadControl. This is the static answer for the released structure
# -- the honest denominator for a dynamic amplification factor, because the
# connection yields on the way down and the released system is not the
# propped one.
#
# Same substitution as the transient release above, but parked at pseudo-
# time 1 so that a Linear series reproduces the reaction exactly: the
# structure starts in the propped equilibrium. LoadControl then walks the
# time BACKWARDS to 0, which retires the reaction load in n_steps
# increments.
ops.loadConst('-time', 1.0)
ops.remove('sp', N_TIP, DOF_PROP)
ops.timeSeries('Linear', TS_REL)
ops.pattern('Plain', PAT_REL, TS_REL)
ops.load(N_TIP, 0.0, reaction_qs, 0.0, 0.0, 0.0, 0.0)

ops.wipeAnalysis()
ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('UmfPack')
ops.test('NormDispIncr', 1.0e-10, 100)
ops.algorithm('Newton')
ops.integrator('LoadControl', -1.0 / 100)
ops.analysis('Static')
for i in range(100):
    if ops.analyze(1) != 0:
        raise RuntimeError("quasi-static release failed at step %d" % (i + 1))

d_static_rel = -ops.nodeDisp(N_TIP, DOF_PROP)
m_static_rel = ops.eleResponse(ELE_CONN, 'section', 'forces')[2]

print("\n   released structure, no inertia (quasi-static):")
print("      tip deflection     : %8.3f mm" % (1000 * d_static_rel))
print("      connection moment  : %8.3f kN.m" % m_static_rel)
print("   released structure, dynamic:")
print("      peak deflection    : %8.3f mm at t = %.4f s"
      % (1000 * d[i_peak], t[i_peak]))
print("      peak conn. moment  : %8.3f kN.m at t = %.4f s"
      % (m[i_mpeak], t[i_mpeak]))
# That mean sits well below the peak and well ABOVE the quasi-static
# release: the group yielded on the way down and does not come back.
last = [d[i] for i in tail]
print("      settled defl.      : %8.3f mm  (mean of the last cycle, "
      "still swinging +-%.3f mm)"
      % (1000 * d_settled, 1000 * 0.5 * (max(last) - min(last))))
print("      permanent set      : %8.3f mm  (settled - quasi-static: "
      "the group yielded)" % (1000 * (d_settled - d_static_rel)))
# Statics has the last word once the motion stops, and it does not care how
# much the group yielded on the way: the connection must carry w*L^2/2, the
# same moment the quasi-static release ends at. That the ringing settles
# onto it -- while the deflection settles well below the quasi-static one
# -- is the cleanest check in this example, and the clearest statement of
# what "damage" means here.
print("      settled conn. M    : %8.3f kN.m  (statics demands "
      "-w*L^2/2 = %.3f, error %.2e)"
      % (m_settled, -0.5 * W_LINE * L ** 2,
         abs(m_settled + 0.5 * W_LINE * L ** 2)))
v_settled = sum(h['v'][i] for i in tail) / len(tail)
print("      settled conn. V    : %8.3f kN    (statics demands "
      "-w*L   = %.3f, error %.2e)"
      % (v_settled, -W_LINE * L, abs(v_settled + W_LINE * L)))
print("   dynamic amplification:")
print("      on deflection      : %8.3f" % (d[i_peak] / d_static_rel))
print("      on connection M    : %8.3f"
      % (abs(m[i_mpeak]) / abs(m_static_rel)))
print("   peak fastener force   : %8.3f kN  (one fastener's capacity "
      "F0 = %.1f kN)" % (max(f), F0_0))

# ---------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("\nmatplotlib not available -- skipping plots.")
    raise SystemExit

fig, axes = plt.subplots(2, 3, figsize=(16, 9.5))
(ax1, ax2, ax3), (ax4, ax5, ax6) = axes

ax1.plot([1000 * ti for ti in t], [1000 * di for di in d], lw=1.2)
ax1.axhline(1000 * d_propped, ls='--', c='0.4', lw=1,
            label='propped (static)')
ax1.axhline(1000 * d_static_rel, ls='-.', c='C3', lw=1,
            label='released (quasi-static)')
ax1.axhline(1000 * d_settled, ls=':', c='C2', lw=1,
            label='settled (permanent set)')
ax1.axvspan(1000 * T_HOLD, 1000 * (T_HOLD + T_RAMP), color='C1',
            alpha=0.3, label='removal ramp')
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Tip deflection, downward (mm)')
ax1.set_title('Sudden support removal')
ax1.grid(True, linestyle=':')
ax1.legend(fontsize=8)

ax2.plot([1000 * ti for ti in t], m, lw=1.2)
ax2.axhline(nvm_propped[2], ls='--', c='0.4', lw=1, label='propped')
ax2.axhline(m_static_rel, ls='-.', c='C3', lw=1,
            label='released (quasi-static)')
ax2.set_xlabel('Time (ms)')
ax2.set_ylabel('Connection moment Mz (kN.m)')
ax2.set_title('What the connection has to absorb')
ax2.grid(True, linestyle=':')
ax2.legend(fontsize=8)

# --- hysteresis of the whole connection ---------------------------
# Both loaded group DOFs, because they tell different stories. The moment
# loop is fat and centred near zero: the group is driven hard in both
# directions and dissipates through it. The shear loop swings by about
# +-50 % of the weight the connection holds up but never reverses SIGN --
# it oscillates about -w*L instead of about zero -- so the fasteners never
# unload through zero in shear, and the loop stays much thinner than the
# moment one for the same reason.
ax4.plot([1000 * x for x in h['th']], h['m'], lw=0.9)
ax4.plot([1000 * h['th'][0]], [h['m'][0]], 'o', ms=5, c='0.4',
         label='propped (start)')
ax4.plot([1000 * h['th'][i_mpeak]], [h['m'][i_mpeak]], 'x', ms=7, c='C3',
         label='peak M')
ax4.set_xlabel('Connection rotation (mrad)')
ax4.set_ylabel('Connection moment Mz (kN.m)')
ax4.set_title('Connection hysteresis, M-theta')
ax4.grid(True, linestyle=':')
ax4.legend(fontsize=8)

ax5.plot([1000 * x for x in h['uy']], h['v'], lw=0.9, c='C4')
ax5.set_xlabel('Connection slip uy (mm)')
ax5.set_ylabel('Connection shear Vy (kN)')
ax5.axhline(-W_LINE * L, ls='-.', c='C3', lw=1,
            label='-w*L (what it must hold)')
ax5.set_title('Connection hysteresis, V-uy (swings, never reverses)')
ax5.grid(True, linestyle=':')
ax5.legend(fontsize=8)

# --- hysteresis of every dowel, overlaid ---------------------------
# The whole point of the group section in one picture: six copies of the
# SAME slip law, driven to six different amplitudes by their lever arms and
# running at six different load-to-grain angles, so no two loops coincide.
# The group loop in ax4 is their sum. Ordered by (y, x) so the legend reads
# bottom row upwards, but drawn in flatly distinct colours rather than a
# colour ramp: several of these loops overlap closely, and neighbouring
# shades of one ramp make them impossible to tell apart where they do.
order = sorted(range(len(coords)), key=lambda i: (coords[i][1], coords[i][0]))
for rank, i in enumerate(order):
    x_i, y_i = coords[i]
    ax6.plot([1000 * s for s in h['bs'][i]], h['bf'][i], lw=1.1,
             color='C%d' % rank, alpha=0.85,
             label='(%+.0f, %+.0f) mm' % (1000 * x_i, 1000 * y_i))
ax6.set_xlabel('Fastener slip along its own slip direction (mm)')
ax6.set_ylabel('Fastener force (kN)')
ax6.set_title('The six dowel laws that add up to it')
ax6.grid(True, linestyle=':')
ax6.legend(fontsize=7, title='fastener (x, y)', title_fontsize=7)

# Fastener pattern at the instant of peak connection moment: re-run and stop
# there, since the group's state is path dependent and cannot be recovered
# after the fact.
ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 6)

ops.node(N_FIX, 0.0, 0.0, 0.0)
ops.fix(N_FIX, 1, 1, 1, 1, 1, 1)
for i in range(N_ELE + 1):
    ops.node(N_CONN + i, L * i / N_ELE, 0.0, 0.0)

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
ops.element('zeroLengthSection', ELE_CONN, N_FIX, N_CONN, SEC_GROUP,
            '-orient', 1, 0, 0, 0, 1, 0)
ops.fix(N_CONN, 0, 0, 1, 1, 1, 0)
ops.geomTransf('Linear', TRANSF, 0, 0, 1)
for i, etag in enumerate(BEAM_ELES):
    ops.element('elasticBeamColumn', etag, N_CONN + i, N_CONN + i + 1,
                A_BEAM, E_TIMBER, G_TIMBER, J_BEAM, IY_BEAM, IZ_BEAM,
                TRANSF, '-mass', RHO)
ops.fix(N_TIP, 0, 1, 0, 0, 0, 0)

ops.timeSeries('Linear', TS_GRAV)
ops.pattern('Plain', PAT_GRAV, TS_GRAV)
ops.eleLoad('-ele', *BEAM_ELES, '-type', '-beamUniform', -W_LINE, 0.0)
ops.wipeAnalysis()
ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('UmfPack')
ops.test('NormDispIncr', 1.0e-10, 50)
ops.algorithm('Newton')
ops.integrator('LoadControl', 1.0 / 20)
ops.analysis('Static')
for i in range(20):
    if ops.analyze(1) != 0:
        raise RuntimeError("static phase failed at step %d" % (i + 1))
ops.reactions()
r2 = ops.nodeReaction(N_TIP, DOF_PROP)

ops.loadConst('-time', 0.0)
ops.remove('sp', N_TIP, DOF_PROP)
ops.timeSeries('Path', TS_REL,
               '-time', [0.0, T_HOLD, T_HOLD + T_RAMP, T_END],
               '-values', [1.0, 1.0, 0.0, 0.0])
ops.pattern('Plain', PAT_REL, TS_REL)
ops.load(N_TIP, 0.0, r2, 0.0, 0.0, 0.0, 0.0)

ops.wipeAnalysis()
lam = ops.eigen('-genBandArpack', N_MODES)
ops.modalDamping(ZETA)

ops.wipeAnalysis()
ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('UmfPack')
ops.test('NormDispIncr', 1.0e-10, 50)
ops.algorithm('Newton')
ops.integrator('Newmark', 0.5, 0.25)
ops.analysis('Transient')
while ops.getTime() < t[i_mpeak] - 0.5 * DT:
    if ops.analyze(1, DT) != 0:
        break

f_peak = ops.eleResponse(ELE_CONN, 'section', 'boltForces')
xs = [x for x, _ in coords]
ys = [y for _, y in coords]
fx = [f_peak[2 * k] for k in range(len(coords))]
fy = [f_peak[2 * k + 1] for k in range(len(coords))]

fmag = [(u * u + v * v) ** 0.5 for u, v in zip(fx, fy)]
fmax = max(fmag) if fmag else 0.0
span = max(max(xs) - min(xs), max(ys) - min(ys), 1.0e-6)
force_scale = 0.4 * span / fmax if fmax > 0 else 1.0

mx, my, _ = ops.nodeCoord(N_CONN)
cr = ops.eleResponse(ELE_CONN, 'section', 'centerOfRotation')

ax3.scatter(xs, ys, s=70, c='saddlebrown', zorder=3, label='dowels')
ax3.scatter([mx], [my], marker='+', s=140, c='black', zorder=3,
            linewidths=2, label='reference point')
if cr and abs(cr[0]) < 10 * max(abs(v) for v in xs + ys):
    ax3.scatter([cr[0]], [cr[1]], marker='x', s=110, c='tab:blue',
                zorder=3, linewidths=2, label='centre of rotation')
for x, y in zip(xs, ys):
    ax3.plot([mx, x], [my, y], color='0.75', linewidth=0.8,
             linestyle='--', zorder=1)

ax3.quiver(xs, ys, [u * force_scale for u in fx],
           [v * force_scale for v in fy],
           angles='xy', scale_units='xy', scale=1,
           color='crimson', width=0.006, zorder=2,
           label='dowel force')

for j, ((x, y), mm) in enumerate(zip(coords, fmag)):
    ax3.annotate("%d\n%.1f kN" % (j, mm), (x, y),
                 textcoords="offset points", xytext=(7, 7), fontsize=7)

ax3.set_aspect('equal', adjustable='datalim')
ax3.set_xlabel('x (m, parallel to grain)')
ax3.set_ylabel('y (m, perpendicular to grain)')
ax3.grid(True, linestyle=':')
ax3.legend(loc='best', fontsize=8)
ax3.set_title('Fastener forces at peak M (t = %.1f ms)' % (1000 * t[i_mpeak]))

fig.suptitle('Sudden support removal on a %dx%d dowel-group connection '
             '(BoltGroupSection2d, 3D transient, incremental modal '
             'damping)' % (N_ROWS, N_COLS))
fig.tight_layout()
plt.show()
