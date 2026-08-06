.. _batch-input-commands:

Batch input: ``nodes`` and ``elements``
=======================================

.. note::

   These are OpenSeesPy (Python) commands.  The canonical command reference lives
   in the separate documentation repository
   `OpenSees/OpenSeesDocumentation <https://github.com/OpenSees/OpenSeesDocumentation>`_
   (rendered at https://OpenSees.github.io/OpenSeesDocumentation); this file is
   kept in the source tree as the authoring copy and can be contributed upstream
   verbatim.

``nodes`` and ``elements`` create many nodes or elements in a single call using
**columnar (broadcast) input**.  They are batch forms of the single-item
:func:`node` and :func:`element` commands and produce exactly the same model —
every item is dispatched through the identical machinery, so all node/element
types, options and error checks behave the same.

Broadcasting rules
------------------

For both commands the arguments are interpreted item-by-item:

* Any argument that is a **list or tuple** is a *per-item column*: its ``i``-th
  entry is used for item ``i``.
* Any argument that is a **scalar** (int, float, or string flag) is *broadcast*
  unchanged to every item.
* All list/tuple arguments must share the **same length** ``N`` (the number of
  items).  Mismatched lengths raise ``ValueError``.
* If a per-item entry is *itself* a list/tuple, it is **spliced** in as several
  tokens.  This lets a list-of-lists supply a distinct vector (e.g. a mass row)
  per item.
* If no argument is a list/tuple, a single item is created.

``nodes``
---------

.. function:: nodes(idList, xInput, yInput[, zInput], *args)

   Create a batch of nodes.  Positionally this mirrors :func:`node`
   (``node(tag, x, y[, z], *options)``) with each field promoted to a column:

   * ``idList`` -- list of node tags.
   * ``xInput``, ``yInput``, ``zInput`` -- coordinate columns (lists) or scalars
     to broadcast.  The number of coordinates required follows the model
     ``ndm`` set by :func:`model`.
   * ``*args`` -- any options accepted by :func:`node` (``-mass``, ``-disp``,
     ``-vel``, ``-ndf`` ...), broadcast or columnar under the rules above.

**Examples**

.. code-block:: python

   ops.model('basic', '-ndm', 2, '-ndf', 3)

   # three 2D nodes; y broadcast to 0.0
   ops.nodes([1, 2, 3], [0.0, 1.0, 2.0], 0.0)

   # optional flag + values broadcast to every node
   ops.nodes([1, 2, 3], [0.0, 1.0, 2.0], 0.0, '-mass', mx, my, mz)

   # per-node mass rows via a list-of-lists (spliced per item)
   ops.nodes([1, 2], [0.0, 1.0], 0.0,
             '-mass', [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])

``elements``
------------

.. function:: elements(eleType, idList, *args)

   Create a batch of elements of a single ``eleType``.  Positionally this mirrors
   :func:`element` (``element(eleType, tag, *eleArgs)``):

   * ``eleType`` -- element type string (a scalar; broadcast to every item).
   * ``idList`` -- list of element tags.
   * ``*args`` -- the remaining :func:`element` arguments for ``eleType``
     (connected node tags, then the type-specific properties), each given as a
     per-item column (list) or a broadcast scalar.

**Example**

.. code-block:: python

   # two elasticBeamColumn elements sharing A, E, Iz and the transform tag;
   # tags and end-node connectivity vary per element
   ops.geomTransf('Linear', 1)
   ops.elements('elasticBeamColumn',
                [1, 2],        # element tags
                [1, 2],        # i-node of each element
                [2, 3],        # j-node of each element
                A, E, Iz, 1)   # broadcast properties + transform tag

Notes
-----

* ``elements`` builds one element *type* per call (the type string is a scalar).
  Mixing types requires one ``elements`` call per type.
* Because each item is handed to the existing :func:`node`/:func:`element`
  handler, any element added in the future is supported automatically, and all
  option flags parse exactly as in the single-item commands.
