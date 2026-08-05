"""Tests for the batch-input commands ``nodes`` and ``elements``.

These commands create many nodes/elements in a single call using columnar
(broadcast) input: any list/tuple argument is a per-item column, any scalar
argument is broadcast to every item.  Each item is dispatched through the same
machinery as the single-item ``node``/``element`` commands, so behaviour must be
identical to building the model one item at a time.
"""
try:
   import opensees as ops
except ModuleNotFoundError:
   import openseespy.opensees as ops

from math import isclose
import pytest

# elasticBeamColumn cantilever properties (two equal segments of length L)
L = 48.0
E = 29000.0
A = 20.0
I = 800.0
P = 10.0


def test_nodes_equivalent_to_singles():
   """nodes(ids, xs, ys) must place every node exactly like node(id, x, y)."""
   ops.wipe()
   ops.model('basic', '-ndm', 2, '-ndf', 3)

   ids = [1, 2, 3, 4]
   xs = [0.0, 1.5, 3.0, 4.5]
   ys = [0.0, 10.0, 20.0, 30.0]
   ops.nodes(ids, xs, ys)

   for tag, x, y in zip(ids, xs, ys):
      cx, cy = ops.nodeCoord(tag)
      assert isclose(cx, x)
      assert isclose(cy, y)


def test_nodes_scalar_broadcast():
   """A scalar coordinate is broadcast to every node."""
   ops.wipe()
   ops.model('basic', '-ndm', 2, '-ndf', 3)

   ids = [1, 2, 3]
   xs = [0.0, 5.0, 10.0]
   ops.nodes(ids, xs, 0.0)  # y broadcast to 0.0

   for tag, x in zip(ids, xs):
      cx, cy = ops.nodeCoord(tag)
      assert isclose(cx, x)
      assert isclose(cy, 0.0)


def test_nodes_3d():
   """Batch nodes honour the model ndm (3 coordinates in 3D)."""
   ops.wipe()
   ops.model('basic', '-ndm', 3, '-ndf', 6)

   ids = [1, 2]
   ops.nodes(ids, [0.0, 1.0], [0.0, 2.0], [0.0, 3.0])
   assert ops.nodeCoord(2) == [1.0, 2.0, 3.0]


def _cantilever_singles():
   ops.wipe()
   ops.model('basic', '-ndm', 2, '-ndf', 3)
   ops.node(1, 0.0, 0.0)
   ops.node(2, L, 0.0)
   ops.node(3, 2 * L, 0.0)
   ops.fix(1, 1, 1, 1)
   ops.geomTransf('Linear', 1)
   ops.element('elasticBeamColumn', 1, 1, 2, A, E, I, 1)
   ops.element('elasticBeamColumn', 2, 2, 3, A, E, I, 1)
   ops.timeSeries('Constant', 1)
   ops.pattern('Plain', 1, 1)
   ops.load(3, 0.0, P, 0.0)
   ops.analysis('Static', '-noWarnings')
   ops.analyze(1)
   ops.reactions()


def _cantilever_batch():
   ops.wipe()
   ops.model('basic', '-ndm', 2, '-ndf', 3)
   ops.nodes([1, 2, 3], [0.0, L, 2 * L], 0.0)
   ops.fix(1, 1, 1, 1)
   ops.geomTransf('Linear', 1)
   # per-element tag + connectivity as lists; shared section props broadcast
   ops.elements('elasticBeamColumn', [1, 2], [1, 2], [2, 3], A, E, I, 1)
   ops.timeSeries('Constant', 1)
   ops.pattern('Plain', 1, 1)
   ops.load(3, 0.0, P, 0.0)
   ops.analysis('Static', '-noWarnings')
   ops.analyze(1)
   ops.reactions()


def test_elements_batch_analysis():
   """A model built with nodes/elements must analyse identically to node/element."""
   _cantilever_singles()
   disp_single = ops.nodeDisp(3, 2)
   reac_single = ops.nodeReaction(1, 2)
   moment_single = ops.nodeReaction(1, 3)
   tags_single = ops.getEleTags()

   _cantilever_batch()
   assert ops.getEleTags() == tags_single
   assert isclose(ops.nodeDisp(3, 2), disp_single)
   assert isclose(ops.nodeReaction(1, 2), reac_single)
   assert isclose(ops.nodeReaction(1, 3), moment_single)


def test_mass_flag_broadcast():
   """An optional flag and its scalar values broadcast to every node."""
   ops.wipe()
   ops.model('basic', '-ndm', 2, '-ndf', 3)

   ids = [1, 2, 3]
   ops.nodes(ids, [0.0, 1.0, 2.0], 0.0, '-mass', 4.0, 5.0, 6.0)
   for tag in ids:
      assert isclose(ops.nodeMass(tag, 1), 4.0)
      assert isclose(ops.nodeMass(tag, 2), 5.0)
      assert isclose(ops.nodeMass(tag, 3), 6.0)


def test_nodes_per_item_vector_splice():
   """A list-of-lists supplies a distinct mass row per node (splice path)."""
   ops.wipe()
   ops.model('basic', '-ndm', 2, '-ndf', 3)

   ids = [1, 2]
   masses = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
   ops.nodes(ids, [0.0, 1.0], 0.0, '-mass', masses)

   assert isclose(ops.nodeMass(1, 1), 1.0)
   assert isclose(ops.nodeMass(1, 3), 3.0)
   assert isclose(ops.nodeMass(2, 1), 4.0)
   assert isclose(ops.nodeMass(2, 3), 6.0)


def test_length_mismatch_raises():
   """List arguments of differing lengths are rejected."""
   ops.wipe()
   ops.model('basic', '-ndm', 2, '-ndf', 3)
   with pytest.raises(Exception):
      ops.nodes([1, 2, 3], [0.0, 1.0], [0.0, 0.0, 0.0])
