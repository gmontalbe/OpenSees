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
# plane -- but it does mean this file, unlike the 2D example, must
# explicitly restrain the group's out-of-plane DOFs itself. Nothing in
# HankinsonSection2d or zeroLengthSection does that for you: leaving those
# DOFs free gives a singular stiffness matrix (the analysis fails to
# factorize, it does not silently do the wrong thing). See build_model()
# below for where and why.
#
# Run standalone:
#   cd EXAMPLES/ExamplePython && cp ../../build/opensees.so .
#   python HankinsonDowelGroup3d.py

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

# Reuse the 2D example's geometry helper, SAWSMaterial calibration helper,
# and the Schweigler-derived calibration/geometry constants -- none of that
# is 2D-specific, it's just node/material bookkeeping. Importing the module
# does not run its __main__ demo.
import HankinsonDowelGroup as hdg2d


def build_dowel_group_3d(dowel_coords, sec_tag, master_node,
                          ground_node_start=1000, free_node_start=2000,
                          ele_tag_start=1):
    """3D counterpart of HankinsonDowelGroup.build_dowel_group(): same
    layout (one zeroLengthSection per dowel, all rigidLink('beam', ...)ed to
    a shared master node), but every node now carries 6 DOF and the
    connection plane is the global X-Y plane (grain parallel to X, in-plane
    perpendicular to Y, out-of-plane normal Z). Requires
    ops.model('basic', '-ndm', 3, '-ndf', 6).

    '-orient' must give both the local x and y' vectors explicitly here --
    3D zeroLengthSection has no auto y = z cross x shortcut the way the 2D
    case does (see ZeroLengthSection's OPS_ZeroLengthSection parser).

    Each ground (timber-side) node is fixed in all 6 DOF, same as the 2D
    example. The free (fastener-side) node's out-of-plane translation and
    all three rotations are left unfixed, but since HankinsonSection2d
    never resists them, they end up entirely governed by the rigidLink
    constraint -- i.e. slaved to the master node's own out-of-plane DOFs,
    not resisted by anything. Restraining those AT the master node is
    build_model()'s job, not this function's.

    Returns a list of (ele_tag, x, y) in the same order as dowel_coords.
    """
    ops.node(master_node, 0.0, 0.0, 0.0)
    info = []
    for i, (x, y) in enumerate(dowel_coords):
        gnode = ground_node_start + i
        fnode = free_node_start + i
        etag = ele_tag_start + i
        ops.node(gnode, x, y, 0.0)
        ops.fix(gnode, 1, 1, 1, 1, 1, 1)
        ops.node(fnode, x, y, 0.0)
        ops.element('zeroLengthSection', etag, gnode, fnode, sec_tag,
                    '-orient', 1, 0, 0, 0, 1, 0)
        ops.rigidLink('beam', master_node, fnode)
        info.append((etag, x, y))
    return info


DOWEL_MAT, SEC_TAG, MASTER_NODE = 20, 1, 9000


def build_model():
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)

    # One SAWSMaterial, started at the 0-deg calibration; HankinsonSection2d
    # re-interpolates its F0/S0/K1 between the 0-deg and 90-deg values
    # below at every commit -- identical to the 2D example, see
    # HankinsonDowelGroup.saws_dowel_type() for why K1 (not R1) is exposed.
    hdg2d.saws_dowel_type(DOWEL_MAT, hdg2d.F0_0, hdg2d.K0_0, hdg2d.K1_0)
    ops.section('HankinsonSection2d', SEC_TAG, '-paramHankinson', DOWEL_MAT,
                '-param', 'F0', hdg2d.F0_0, hdg2d.F0_90,
                '-param', 'S0', hdg2d.K0_0, hdg2d.K0_90,
                '-param', 'K1', hdg2d.K1_0, hdg2d.K1_90)

    coords = hdg2d.dowel_coordinates_rp(hdg2d.N_ROWS, hdg2d.N_COLUMNS,
                                         hdg2d.A1, hdg2d.A2, hdg2d.A3, hdg2d.A4)
    dowels = build_dowel_group_3d(coords, SEC_TAG, MASTER_NODE)

    # No dowel resists the master node's out-of-plane translation (Uz) or
    # any rotation (Rx, Ry, Rz) directly -- see the module docstring. The
    # in-plane rotation Rz DOES pick up stiffness despite that, but only
    # indirectly: it comes from the lever-arm coupling of each dowel's
    # in-plane (Fx, Fy) force through the rigidLink constraint, exactly as
    # in the 2D example, not from any element resisting Rz itself. Uz, Rx
    # and Ry get NO stiffness from anywhere in this model and must be
    # restrained explicitly or the analysis is singular in those
    # directions. Physically this fixation stands in for whatever restrains
    # the connection out of its own plane in the real structure (the timber
    # member's own flexural continuity, other fasteners, a test rig) --
    # exactly the "handle the missing DOFs yourself" requirement discussed
    # for HankinsonSection2d in 3D.
    ops.fix(MASTER_NODE, 0, 0, 1, 1, 1, 0)
    return dowels


GROUP_SEC_TAG, GROUP_ELE_TAG, GROUP_GROUND_NODE = 2, 500, 500


def build_group_model():
    """The same 3D connection as ONE element, via BoltGroupSection2d -- the
    3D counterpart of HankinsonDowelGroup.build_group_model().

    The group plane is placed entirely by the element's '-orient' (local x
    along the grain, local y' in-plane perpendicular); the fastener offsets
    inside the section are always expressed in that local frame. Nothing
    else about the section is 3D-aware: it still reports {P, Vy, Mz}, so
    Uz/Rx/Ry still get no stiffness from anywhere and must be restrained,
    exactly as in build_model() above.

    Two things the hand-built model needs and this one does not: there are
    no rigidLink constraints, so 'Transformation' is safe again (see
    apply_normal_force()'s note on why it is not, above), and the model has
    3 free DOF instead of 9 nodes' worth plus 9 constraints.

    Returns (ele_tag, coords).
    """
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)

    hdg2d.saws_dowel_type(DOWEL_MAT, hdg2d.F0_0, hdg2d.K0_0, hdg2d.K1_0)
    ops.section('HankinsonSection2d', SEC_TAG, '-paramHankinson', DOWEL_MAT,
                '-param', 'F0', hdg2d.F0_0, hdg2d.F0_90,
                '-param', 'S0', hdg2d.K0_0, hdg2d.K0_90,
                '-param', 'K1', hdg2d.K1_0, hdg2d.K1_90)
    ops.section('BoltGroupSection2d', GROUP_SEC_TAG, '-section', SEC_TAG,
                '-grid', hdg2d.N_ROWS, hdg2d.N_COLUMNS, hdg2d.A1, hdg2d.A2)

    ops.node(GROUP_GROUND_NODE, 0.0, 0.0, 0.0)
    ops.fix(GROUP_GROUND_NODE, 1, 1, 1, 1, 1, 1)
    ops.node(MASTER_NODE, 0.0, 0.0, 0.0)
    ops.element('zeroLengthSection', GROUP_ELE_TAG,
                GROUP_GROUND_NODE, MASTER_NODE, GROUP_SEC_TAG,
                '-orient', 1, 0, 0, 0, 1, 0)
    ops.fix(MASTER_NODE, 0, 0, 1, 1, 1, 0)

    coords = hdg2d.dowel_coordinates_rp(hdg2d.N_ROWS, hdg2d.N_COLUMNS,
                                         hdg2d.A1, hdg2d.A2, hdg2d.A3, hdg2d.A4)
    return GROUP_ELE_TAG, coords


def apply_normal_force(n_target, pattern_tag=5, n_steps=10):
    """Ramp a constant axial preload N (global/local x = parallel-to-grain,
    per build_dowel_group_3d's '-orient 1 0 0 0 1 0') onto the master node,
    then loadConst() it so it stays in place under a later push_moment()
    pattern. See HankinsonDowelGroup.apply_normal_force() for the sizing
    rationale (unchanged: same calibration, same 3x3 group).
    """
    ops.timeSeries('Linear', pattern_tag)
    ops.pattern('Plain', pattern_tag, pattern_tag)
    ops.load(MASTER_NODE, n_target, 0.0, 0.0, 0.0, 0.0, 0.0)
    # wipeAnalysis() first: if a previous analysis object still exists (e.g.
    # push_moment() ran before this), the analysis commands below silently
    # keep the OLD constraint handler and print "WARNING can't set handler
    # after analysis is created" instead of actually reconfiguring it.
    ops.wipeAnalysis()
    # 'Lagrange', NOT 'Transformation': in this build, TransformationConstraintHandler
    # mishandles a node that is BOTH (a) partially ops.fix()ed and (b) the
    # RETAINED node of an MP_Constraint with other, still-free DOFs -- which
    # is exactly MASTER_NODE here (build_model() fixes its Uz/Rx/Ry, every
    # dowel's rigidLink('beam', MASTER_NODE, ...) retains it, and Ux/Uy/Rz
    # stay free). Verified in total isolation (one zeroLengthSection + one
    # rigidLink + a Bidirectional section, no Hankinson code involved at
    # all, in BOTH 2D and 3D): across n repeated LoadControl analyze(1)
    # calls to the same cumulative load, the master's OWN reported
    # displacement comes out (n+1)/2 times too large, while the slave
    # node's displacement (and every reaction force) stays exactly correct
    # throughout -- so the bug is invisible in equilibrium/force checks and
    # only shows up when you read the master's own displacement after a
    # multi-step analysis. A master node with NO fixed DOFs (as in the 2D
    # example, which never fixes MASTER_NODE) does not trigger it. Only
    # 'Transformation' is affected; 'Lagrange' reproduces the theoretical/
    # 2D-model rotational stiffness exactly for any number of steps.
    ops.constraints('Lagrange')
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
    """Push the in-plane moment Mz (about the global/local Z axis, out of
    the dowel group's X-Y plane) -- the 3D counterpart of the 2D example's
    third DOF (theta_z), now DOF 6 of 6 rather than DOF 3 of 3.
    """
    ops.timeSeries('Linear', pattern_tag)
    ops.pattern('Plain', pattern_tag, pattern_tag)
    ops.load(MASTER_NODE, 0.0, 0.0, 0.0, 0.0, 0.0, m_target)
    # See the matching comments in apply_normal_force(): without
    # wipeAnalysis() first, this reconfiguration is silently ignored (with a
    # warning) whenever apply_normal_force() ran first and left an analysis
    # object in place; and 'Lagrange' (not 'Transformation') is required
    # because MASTER_NODE is both partially fixed and the retained node of
    # every dowel's rigidLink.
    ops.wipeAnalysis()
    ops.constraints('Lagrange')
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
        phi_hist.append(ops.nodeDisp(MASTER_NODE, 6))
    return m_hist, phi_hist


def plot_dowel_group_3d(dowels, ax=None, force_scale=None, title=None):
    """3D counterpart of HankinsonDowelGroup.plot_dowel_group(): the same
    in-plane force-vector picture, adapted only for the extra reaction/
    coordinate components 3D nodes carry (only the Fx/Fy components are
    ever nonzero -- see the module docstring -- so there is nothing new to
    plot, just more numbers to unpack per node).
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
        Rx, Ry, Rz, Mx, My, Mz = ops.nodeReaction(gnode)
        fx.append(-Rx)
        fy.append(-Ry)

    fmag = [(u * u + v * v) ** 0.5 for u, v in zip(fx, fy)]
    fmax = max(fmag) if fmag else 0.0
    if force_scale is None:
        span = max(max(xs) - min(xs), max(ys) - min(ys), 1.0e-6)
        force_scale = 0.4 * span / fmax if fmax > 0 else 1.0

    mx, my, _ = ops.nodeCoord(MASTER_NODE)

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
    ax.set_ylabel('y (m, perpendicular to grain, in-plane)')
    ax.grid(True, linestyle=':')
    ax.legend(loc='best', fontsize=8)
    if title:
        ax.set_title(title)
    return ax


# Axial preload levels (kN) to sweep in the __main__ demo below -- identical
# choice to HankinsonDowelGroup.py, see its NORMAL_FORCE_LEVELS comment for
# the sizing rationale (same calibration, same 3x3 group).
NORMAL_FORCE_LEVELS = [0.0, 40.0, 80.0]
M_TARGET = 10

if __name__ == '__main__':
    try:
        import matplotlib.pyplot as plt
        have_plt = True
        fig, axes = plt.subplots(1, len(NORMAL_FORCE_LEVELS), figsize=(16, 5.5))
        if len(NORMAL_FORCE_LEVELS) == 1:
            axes = [axes]
    except ImportError:
        have_plt = False
        print("\nmatplotlib not available -- skipping plots.\n")

    # Same cross-check as the 2D example: the single-element group and the
    # 9-element + rigidLink model are the same mechanics.
    build_model()
    apply_normal_force(NORMAL_FORCE_LEVELS[-1])
    phi_multi = push_moment(m_target=M_TARGET, n_steps=80)[1][-1]
    build_group_model()
    apply_normal_force(NORMAL_FORCE_LEVELS[-1])
    phi_single = push_moment(m_target=M_TARGET, n_steps=80)[1][-1]
    print("cross-check at N = %.0f kN: phi(9 elements + rigidLink) = %.12e"
          % (NORMAL_FORCE_LEVELS[-1], phi_multi))
    print("                            phi(1 BoltGroupSection2d)  = %.12e"
          % phi_single)
    print("                            difference                 = %.3e\n"
          % abs(phi_multi - phi_single))

    curves = {}
    for i, n_val in enumerate(NORMAL_FORCE_LEVELS):
        ele_tag, coords = build_group_model()
        if n_val:
            apply_normal_force(n_val)
        m_hist, phi_hist = push_moment(m_target=M_TARGET, n_steps=80)
        curves[n_val] = (phi_hist, m_hist)

        i_mid = len(m_hist) // 3
        k_mid = m_hist[i_mid] / phi_hist[i_mid]
        k_final = m_hist[-1] / phi_hist[-1]
        print("N = %5.1f kN   K_secant(M=%.1f) = %7.2f   "
              "K_secant(M=%.1f) = %6.2f kN.m/rad"
              % (n_val, m_hist[i_mid], k_mid, m_hist[-1], k_final))

        if have_plt:
            # the 2D plotter reads the group section's own 'boltForces', which
            # is dimension-agnostic -- unlike plot_dowel_group_3d(), which has
            # to unpack 6-component nodal reactions
            hdg2d.plot_bolt_group(ele_tag, coords, ax=axes[i],
                                  title="N = %.0f kN" % n_val)

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
