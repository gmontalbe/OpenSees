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

import HankinsonDowelGroup as hdg
import BoltGroupCantilever as bgc     # geometry, properties, calibration

L = bgc.L
B, H = bgc.B, bgc.H
E_TIMBER = bgc.E_TIMBER
N_ROWS, N_COLS, A1, A2 = bgc.N_ROWS, bgc.N_COLS, bgc.A1, bgc.A2

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

MAT_DOWEL, SEC_BOLT, SEC_GROUP, SEC_AGG = 20, 10, 11, 12
MAT_VZ, MAT_MY, MAT_T = 31, 32, 33
TRANSF = 1
N_FIX, N_CONN, N_TIP = 1, 2, 3
ELE_CONN, ELE_BEAM = 1, 2

P_MAX, N_STEPS = bgc.P_MAX, bgc.N_STEPS


def build_cantilever_3d(restrain='aggregate'):
    """3D cantilever on a BoltGroupSection2d base connection.

    restrain='aggregate' wraps the group in a SectionAggregator carrying
    elastic Vz/My/T; restrain='fix' leaves the group bare and restrains the
    connection node's out-of-plane DOFs instead. Returns the fastener
    offsets.
    """
    if restrain not in ('aggregate', 'fix'):
        raise ValueError("restrain must be 'aggregate' or 'fix'")

    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)

    ops.node(N_FIX, 0.0, 0.0, 0.0)
    ops.fix(N_FIX, 1, 1, 1, 1, 1, 1)
    ops.node(N_CONN, 0.0, 0.0, 0.0)
    ops.node(N_TIP, L, 0.0, 0.0)

    hdg.saws_dowel_type(MAT_DOWEL, hdg.F0_0, hdg.K0_0, hdg.K1_0)
    ops.section('HankinsonSection2d', SEC_BOLT, '-paramHankinson', MAT_DOWEL,
                '-param', 'F0', hdg.F0_0, hdg.F0_90,
                '-param', 'S0', hdg.K0_0, hdg.K0_90,
                '-param', 'K1', hdg.K1_0, hdg.K1_90)
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

    return hdg.dowel_coordinates_rp(N_ROWS, N_COLS, A1, A2, 0.0, 0.0)


def pushover_3d(p_max=P_MAX, n_steps=N_STEPS, direction='in-plane'):
    """Ramp a tip load either in-plane (-Y, exercising the group's Vy/Mz) or
    out-of-plane (-Z, which the group cannot resist at all). Returns
    (P history, tip deflection history along the loaded axis)."""
    dof = 2 if direction == 'in-plane' else 3
    load = [0.0] * 6
    load[dof - 1] = -p_max

    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(N_TIP, *load)
    ops.wipeAnalysis()
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('UmfPack')
    ops.test('NormDispIncr', 1.0e-10, 50)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0 / n_steps)
    ops.analysis('Static')

    p_hist, d_hist = [0.0], [0.0]
    for i in range(n_steps):
        if ops.analyze(1) != 0:
            print("   stopped: no convergence at step %d of %d" % (i + 1, n_steps))
            break
        p_hist.append(p_max * (i + 1) / n_steps)
        d_hist.append(-ops.nodeDisp(N_TIP, dof))
    return p_hist, d_hist


def bolt_forces():
    """The group's per-fastener forces. This reaches the group even through
    the SectionAggregator: the aggregator does not recognise 'boltForces',
    and its setResponse falls through to the wrapped section."""
    return ops.eleResponse(ELE_CONN, 'section', 'boltForces')


if __name__ == '__main__':
    print("3D cantilever on a %dx%d dowel-group connection "
          "(%d fasteners, L = %.2f m)\n" % (N_ROWS, N_COLS, N_ROWS * N_COLS, L))

    # --- the two out-of-plane treatments must agree in plane ----------
    results = {}
    for how in ('aggregate', 'fix'):
        coords = build_cantilever_3d(how)
        p_h, d_h = pushover_3d()
        nvm = ops.eleResponse(ELE_CONN, 'section', 'forces')
        results[how] = (p_h, d_h, nvm)
        print("   restrain='%s': tip %8.3f mm at P = %.1f kN" %
              (how, 1000 * d_h[-1], p_h[-1]))

    d_agg, d_fix = results['aggregate'][1][-1], results['fix'][1][-1]
    print("   difference between the two: %.3e mm\n"
          % (1000 * abs(d_agg - d_fix)))

    # --- and both must reproduce the 2D model exactly -----------------
    bgc.build_cantilever()
    _, d2d = bgc.pushover()
    print("cross-check against the 2D model (BoltGroupCantilever.py):")
    print("   2D tip deflection : %.12e m" % d2d[-1])
    print("   3D tip deflection : %.12e m" % d_agg)
    print("   difference        : %.3e m\n" % abs(d2d[-1] - d_agg))

    # --- out-of-plane: the group is inert, the aggregate carries it ---
    print("Out-of-plane load case (the group resists nothing here):")
    build_cantilever_3d('aggregate')
    p_out, d_out = pushover_3d(p_max=5.0, direction='out-of-plane')
    f_out = bolt_forces()
    print("   tip deflection out of plane : %8.3f mm at P = %.1f kN"
          % (1000 * d_out[-1], p_out[-1]))
    print("   largest fastener force       : %.3e kN  (exactly zero: the "
          "group\n                                  sees no deformation in "
          "Vz/My/T)" % max(abs(v) for v in f_out))

    # --- plots ---------------------------------------------------------
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nmatplotlib not available -- skipping plots.")
        raise SystemExit

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

    for how, style in (('aggregate', '-'), ('fix', '--')):
        p_h, d_h, _ = results[how]
        ax1.plot([1000 * d for d in d_h], p_h, style,
                 label="in-plane, restrain='%s'" % how)
    ax1.plot([1000 * d for d in d_out], p_out, ':',
             label='out-of-plane (aggregated stiffness only)')
    ax1.set_xlabel('Tip deflection (mm)')
    ax1.set_ylabel('Tip load P (kN)')
    ax1.set_title('In-plane vs out-of-plane response')
    ax1.grid(True, linestyle=':')
    ax1.legend(fontsize=8)

    build_cantilever_3d('aggregate')
    pushover_3d()
    hdg.plot_bolt_group(ELE_CONN, coords, ax=ax2, ref_node=N_CONN,
                        title='Fastener forces at P = %.0f kN (in-plane)' % P_MAX)

    fig.suptitle('Timber cantilever on a %dx%d dowel-group connection, 3D '
                 'model\n(BoltGroupSection2d + SectionAggregator for the '
                 'out-of-plane DOFs)' % (N_ROWS, N_COLS))
    fig.tight_layout()
    plt.show()
