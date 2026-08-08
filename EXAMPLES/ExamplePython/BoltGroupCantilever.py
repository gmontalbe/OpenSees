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

# The fastener calibration, the SAWSMaterial helper and the grid helper are
# all reused as-is -- none of that is specific to the dowel-group test model.
import HankinsonDowelGroup as hdg


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

MAT_DOWEL, SEC_BOLT, SEC_GROUP = 20, 10, 11
TRANSF = 1
N_FIX, N_CONN, N_TIP = 1, 2, 3
ELE_CONN, ELE_BEAM = 1, 2

P_MAX = 30.0                 # tip load (kN), enough to yield the group
N_STEPS = 120


def build_cantilever(rigid_base=False, axial_preload=False):
    """Cantilever on a BoltGroupSection2d base connection.

    With rigid_base=True the connection is omitted entirely (the beam is
    built straight off the fixed node), which is the reference case the demo
    compares against -- it isolates how much of the tip deflection is the
    connection rather than the beam.

    Returns the list of fastener (x, y) offsets, for labelling output.
    """
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    ops.node(N_FIX, 0.0, 0.0)
    ops.fix(N_FIX, 1, 1, 1)
    ops.node(N_TIP, L, 0.0)

    coords = hdg.dowel_coordinates_rp(N_ROWS, N_COLS, A1, A2, 0.0, 0.0)

    if rigid_base:
        beam_start = N_FIX
    else:
        # One fastener: the Schweigler dowel law, Hankinson-interpolated.
        hdg.saws_dowel_type(MAT_DOWEL, hdg.F0_0, hdg.K0_0, hdg.K1_0)
        ops.section('HankinsonSection2d', SEC_BOLT, '-paramHankinson', MAT_DOWEL,
                    '-param', 'F0', hdg.F0_0, hdg.F0_90,
                    '-param', 'S0', hdg.K0_0, hdg.K0_90,
                    '-param', 'K1', hdg.K1_0, hdg.K1_90)
        # The whole group: that prototype copied onto the grid below. Local x
        # (the element's '-orient' below) is the beam axis and the grain
        # direction; local y runs through the beam depth.
        ops.section('BoltGroupSection2d', SEC_GROUP, '-section', SEC_BOLT,
                    '-grid', N_ROWS, N_COLS, A1, A2)

        ops.node(N_CONN, 0.0, 0.0)
        ops.element('zeroLengthSection', ELE_CONN, N_FIX, N_CONN, SEC_GROUP,
                    '-orient', 1, 0, 0)
        beam_start = N_CONN

    ops.geomTransf('Linear', TRANSF)
    ops.element('elasticBeamColumn', ELE_BEAM, beam_start, N_TIP,
                A_BEAM, E_TIMBER, IZ_BEAM, TRANSF)

    return coords


def pushover(p_max=P_MAX, n_steps=N_STEPS, axial=0.0):
    """Ramp a downward tip load (and optionally a constant axial preload
    applied first and held). Returns (P history, tip deflection history)."""
    if axial:
        ops.timeSeries('Linear', 9)
        ops.pattern('Plain', 9, 9)
        ops.load(N_TIP, axial, 0.0, 0.0)
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
    ops.load(N_TIP, 0.0, -p_max, 0.0)
    # wipeAnalysis() first, otherwise a previously created analysis object is
    # silently kept and the settings below are ignored with a warning.
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
        d_hist.append(-ops.nodeDisp(N_TIP, 2))
    return p_hist, d_hist


def deflection_split():
    """Decompose the current tip deflection into the three contributions.

    Only valid for the flexible-base model (it reads the connection node).
    The identity delta = u_y + theta*L + flexure is exact for small
    rotations, which is the regime an elasticBeamColumn with a Linear
    transformation assumes anyway.
    """
    uy = -ops.nodeDisp(N_CONN, 2)
    th = -ops.nodeDisp(N_CONN, 3)
    tip = -ops.nodeDisp(N_TIP, 2)
    return uy, th * L, tip - uy - th * L, tip


def connection_forces():
    """(N, V, M) carried by the connection, in its local frame."""
    return ops.eleResponse(ELE_CONN, 'section', 'forces')


if __name__ == '__main__':
    print(__doc__ if __doc__ else "")
    print("Cantilever on a %dx%d dowel-group connection "
          "(%d fasteners, L = %.2f m)\n" % (N_ROWS, N_COLS, N_ROWS * N_COLS, L))

    # --- flexible base ------------------------------------------------
    coords = build_cantilever()
    p_flex, d_flex = pushover()
    uy, rot, flex, tip = deflection_split()
    nvm = connection_forces()

    # --- rigid base, same beam ---------------------------------------
    build_cantilever(rigid_base=True)
    p_rig, d_rig = pushover()

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

    # --- axial preload: N-M interaction at the connection -------------
    # A tension preload eats into each fastener's capacity, so the same tip
    # load costs more deflection -- the interaction the group section exists
    # to capture, and which no uncoupled spring-per-DOF connection model can
    # reproduce.
    print("Axial preload sweep (tip axial force held constant, then pushed):")
    curves = {}
    for n_val in (0.0, 60.0, 120.0):
        build_cantilever()
        p_h, d_h = pushover(axial=n_val)
        curves[n_val] = (d_h, p_h)
        nvm_n = connection_forces()
        print("   N = %6.1f kN -> P %5.2f kN, tip %7.3f mm, "
              "connection (N, M) = (%7.3f kN, %7.3f kN.m)"
              % (n_val, p_h[-1], 1000 * d_h[-1], nvm_n[0], nvm_n[2]))

    # --- plots ---------------------------------------------------------
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

    # fastener force pattern at the end of the unpreloaded pushover
    build_cantilever()
    pushover()
    hdg.plot_bolt_group(ELE_CONN, coords, ax=ax2, ref_node=N_CONN,
                        title='Fastener forces at P = %.0f kN' % P_MAX)

    fig.suptitle('Timber cantilever on a %dx%d dowel-group connection '
                 '(BoltGroupSection2d, one element)' % (N_ROWS, N_COLS))
    fig.tight_layout()
    plt.show()
