.. _bolt-group-section-2d:

``BoltGroupSection2d``: a fastener group as one element
=======================================================

.. note::

   The canonical command reference lives in the separate documentation
   repository `OpenSees/OpenSeesDocumentation
   <https://github.com/OpenSees/OpenSeesDocumentation>`_ (rendered at
   https://OpenSees.github.io/OpenSeesDocumentation); this file is kept in the
   source tree as the authoring copy and can be contributed upstream verbatim.

``BoltGroupSection2d`` turns a whole group of mechanical fasteners (dowels,
bolts, screws) into a **single 3-DOF section** ``{P, Vy, Mz}``, hosted in an
ordinary ``zeroLengthSection`` element.  It takes a **prototype** order-2
section describing *one* fastener, copies it once per fastener, distributes
the reference point's :math:`(u_x, u_y, \theta)` to each copy by that
fastener's lever arm, and sums the results back into the connection's
:math:`(N, V, M)`:

.. math::

   \mathbf{B}_i = \begin{bmatrix} 1 & 0 & -y_i \\ 0 & 1 & x_i \end{bmatrix},
   \qquad
   \mathbf{e}_i = \mathbf{B}_i\,\mathbf{e},
   \qquad
   \mathbf{q} = \sum_i \mathbf{B}_i^{T} \mathbf{q}_i,
   \qquad
   \mathbf{k} = \sum_i \mathbf{B}_i^{T} \mathbf{k}_i \mathbf{B}_i

This is the assembly a ``FiberSection`` performs, with a 2-DOF coupled
"fiber" — a whole fastener — in place of a uniaxial one.

Why one element instead of many
-------------------------------

The same connection can always be built by hand: one ``zeroLengthSection``
per fastener, all tied to a shared master node with ``rigidLink('beam', ...)``
(see ``EXAMPLES/ExamplePython/HankinsonDowelGroup.py``).  The two are the same
mechanics and agree to round-off, because ``rigidLink`` builds an
``MP_Constraint`` with a *constant* transformation matrix — it is itself
linear in the rotation, exactly like :math:`\mathbf{B}_i` above.

What the single-element form buys:

* **3 free DOF instead of ~27 plus 9 multi-point constraints** for a 3×3 group.
* **No constraint handler involved at all**, which sidesteps the
  ``TransformationConstraintHandler`` problem documented in
  ``HankinsonDowelGroup3d.apply_normal_force()`` (a master node that is both
  partially fixed and constraint-retained reports its own displacement
  incorrectly).  ``Transformation`` is safe again.
* **Group-level output** — per-fastener force/slip vectors and the centre of
  rotation from one recorder.

Command
-------

.. code-block:: python

   ops.section('BoltGroupSection2d', secTag, '-section', protoTag,
               # exactly one of:
               '-grid', nRows, nCols, a1, a2,
               '-bolts', n, x1, y1, x2, y2, ...,
               # optional:
               '-contact', matTag, nRows, y1, y2, ...,
               '-rot', 'linear')      # or 'tan'

``-section protoTag``
    The prototype: any **order-2** section, copied once per fastener so every
    fastener carries fully independent state and history.  Component 0 is
    taken as the section-local *x* slip (for timber, parallel to grain) and
    component 1 as local *y* — positionally, not by response code.  A
    prototype not reporting ``{P, Vy}`` warns but is still used.

``-grid nRows nCols a1 a2``
    A rectangular grid referred to its own centroid, i.e. the connection's
    reference point, row-major.  This is
    ``ModelFunctions.dowelCoordinatesRP`` from the Schweigler reference
    implementation, including fastener ordering.  That function's edge
    distances ``a3``/``a4`` are deliberately absent: they shift every
    fastener *and* the centroid by the same amount, so they cancel exactly
    out of centroid-referred coordinates.

``-bolts n x1 y1 ...``
    Explicit offsets from the reference point, in section-local axes.
    ``-dowels`` is accepted as an alias.

``-contact matTag n y1 ...``
    Optional compression/bearing rows at heights :math:`y_j`, acting along
    **local x only** — the timber-to-timber contact that closes on one side
    of a connection in bending, mirroring ``ModelFunctions.contactDisp`` /
    ``contactForce``.  Each row is an independent copy of one
    ``UniaxialMaterial`` (typically ``ElasticPPGap``) driven by
    :math:`u_x - y_j\theta`, contributing
    :math:`\mathbf{B}_j = [1, 0, -y_j]`.  The material is a
    force-vs-displacement law: **any tributary area or contact length must
    already be baked into it.**

``-rot linear|tan``
    Rotation kinematics, default ``linear``.  See below.

Orientation
-----------

The section always reports ``{P, Vy, Mz}``; the group plane is placed
entirely by the host element's ``-orient``, and the fastener offsets are
always in that local frame.  ``ZeroLengthSection`` only fills the ``Mz`` row
of its transformation matrix for ``numDOF`` 6 and 12, so use ``ndf=3`` in 2D,
and ``ndf=6`` in 3D with **both** local axes given:

.. code-block:: python

   ops.element('zeroLengthSection', 1, iNode, jNode, secTag, '-orient', 1,0,0)          # 2D
   ops.element('zeroLengthSection', 1, iNode, jNode, secTag, '-orient', 1,0,0, 0,1,0)   # 3D

.. note::

   In **Tcl**, ``zeroLengthSection`` requires all six ``-orient`` components
   even in 2D — its parser differs from the Python one here, which accepts
   three.  Use ``-orient 1 0 0 0 1 0`` in both dimensions from Tcl.

In 3D the section resists ``P``, ``Vy`` and ``Mz`` only — ``Uz``, ``Rx`` and
``Ry`` get no stiffness from anywhere and must be restrained explicitly, or
the analysis is singular in those directions.

Rotation kinematics
-------------------

``linear`` (default)
    :math:`\theta` is used directly; :math:`\mathbf{B}_i` is constant.

``tan``
    :math:`\tan\theta` is used, as ``ModelFunctions.dowelDisp`` does, with
    :math:`c = \sec^2\theta`.  The work conjugate of :math:`\theta` under
    that map carries the same factor, so
    :math:`M = c \sum_i (x_i q_{iy} - y_i q_{ix})`, and the tangent picks up
    the geometric term :math:`2\tan\theta\,c\,m` in :math:`k_{22}`.

    Note that the reference implementation applies the ``tan`` map to the
    slips but **not** the matching factor to its moment, so it is not energy
    consistent; ``-rot tan`` therefore reproduces the reference's slips
    exactly and its moment only to :math:`O(\theta^2)`.

The two modes agree to :math:`O(\theta^2)`, which at the rotations dowel-type
connections reach is below the staggering error the sub-section's own lagged
Hankinson angle already introduces.  ``linear`` is the default for that
reason, and because it is the mode that matches a ``rigidLink``-built model
exactly.  ``-rot tan`` errors out beyond :math:`|\theta| = 1.5` rad, where the
map breaks down.

Output
------

.. code-block:: python

   ops.eleResponse(ele, 'section', 'boltForces')        # 2*n: (P, Vy) per fastener
   ops.eleResponse(ele, 'section', 'boltSlips')         # 2*n: slips per fastener
   ops.eleResponse(ele, 'section', 'centerOfRotation')  # (x, y)
   ops.eleResponse(ele, 'section', 'contactForces')     # n_contact
   ops.eleResponse(ele, 'section', 'bolt', i, ...)      # forward to fastener i
   ops.eleResponse(ele, 'section', 'contact', j, ...)   # forward to contact row j

The centre of rotation is where the rigid-body slip field vanishes,
:math:`(x, y) = (-u_y/\theta,\; u_x/\theta)` — the point
``ModelFunctions.findCR`` locates geometrically.  It is undefined for a pure
translation and reported as ``(0, 0)`` there; the honest answer is "at
infinity".

.. warning::

   Mind the singular/plural distinction when reading one fastener at a time.
   It is inherited from the prototype, not invented here:
   ``HankinsonSection2d`` intercepts the **singular** ``'force'`` with its own
   *scalar* slip-law force, so only the **plural** ``'bolt i forces'`` falls
   through to the 2-component stress resultant.  ``'boltForces'`` is the
   concatenation of the plural one.

``setParameter`` follows the same shape: ``'bolt' i ...`` and
``'contact' j ...`` address one sub-object, and anything else is offered to
every fastener, so one ``Parameter`` can drive the whole group.

Worked example
--------------

The 3×3 dowel connection of Schweigler et al.'s worked example, with each
fastener a Schweigler-style dowel whose slip law has its ``F0``/``S0``/``K1``
Hankinson-interpolated between a 0° and a 90° calibration by that fastener's
own evolving load-to-grain angle:

.. code-block:: python

   ops.uniaxialMaterial('SAWS', 20, 19.1, 1.91, 0.05, 24460.,
                        0.0189, -0.05, 1.0, 0.0, 0.8, 1.05)
   ops.section('HankinsonSection2d', 10, '-paramHankinson', 20,
               '-param', 'F0', 19.1,   14.18,
               '-param', 'S0', 24460., 19640.,
               '-param', 'K1', 462.,   739.)
   ops.section('BoltGroupSection2d', 1, '-section', 10, '-grid', 3, 3, 0.070, 0.030)

   ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1, 1)
   ops.node(2, 0.0, 0.0)
   ops.element('zeroLengthSection', 1, 1, 2, 1, '-orient', 1, 0, 0)

Because every fastener in a loaded group sits at a *different* load-to-grain
angle, the group's response is not reducible to any single section — which is
what this class exists for.

Runnable versions of this model, in 2D and 3D, are in
``EXAMPLES/ExamplePython/HankinsonDowelGroup.py`` and
``HankinsonDowelGroup3d.py`` (both also build the multi-element form and print
the agreement between the two).  The verification suite is
``tests/test_bolt_group_section2d.py``.

Three further examples put the group in a structure rather than on its own:

``BoltGroupCantilever.py``
    A cantilever whose base connection is the group, in 2D.  Splits the tip
    deflection into connection slip, connection rotation and beam flexure,
    and sweeps an axial preload to show the N–M interaction.

``BoltGroupCantilever3d.py``
    The same model at ``ndm=3``/``ndf=6``, running both ways of supplying the
    missing ``Vz``/``My``/``T`` stiffness (``SectionAggregator`` versus fixing
    those DOFs) and showing they agree.

``BoltGroupSupportRemoval3d.py``
    A 3D **transient** analysis: a propped cantilever on the group, loaded
    statically, then the prop is replaced by its own reaction and that
    reaction ramped to zero over a fraction of a period — the standard
    sudden-support-removal procedure, damped with incremental modal damping
    (``ops.modalDamping``).  It brackets the peak dynamic response with a
    quasi-static release of the identical load path, and checks that the
    ringing settles onto the moment statics demands, :math:`wL^2/2`, while
    the deflection keeps a permanent set.  It also plots the group's own
    :math:`M`–:math:`\theta` and :math:`V`–:math:`u_y` hysteresis against
    the overlaid loops of the individual fasteners it is the sum of, read
    from ``'bolt' i 'state'``.

Caveats
-------

* A group whose fasteners all sit at the reference point has **no rotational
  stiffness**: the ``Mz`` row and column of the tangent are identically zero.
  That is legal but makes the assembled system singular unless the host
  element's relative rotation is restrained elsewhere.  The section warns at
  construction.
* Each fastener holds its own copy of the prototype's sub-material, so memory
  and per-step cost scale with the fastener count.
* The prototype's own limitations carry over unchanged — in particular
  ``HankinsonSection2d``'s Hankinson angle is lagged to the last committed
  state, which is a staggering (explicit) approximation.  Its
  ``'explicitError'`` response monitors that, per fastener.
