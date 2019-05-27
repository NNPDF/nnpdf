# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 13:19:18 2015

@author: zah
"""
from collections import deque
#import weakref


class Node:
    def __init__(self, value, inputs=None, outputs=None):
        if inputs is None:
            inputs = set()
        if outputs is None:
            outputs = set()

        self.inputs = inputs
        self.outputs = outputs
        self.value = value

    def __str__(self):
        return "Node({!r})".format(self.value)

    def __repr__(self):
        return "Node({!r})".format(self.value)


class DAGError(Exception): pass

class CycleError(DAGError):
    def __init__(self, node ,cycle):
        self.node = node
        self.cyle = cycle
        msg = "%s introduces a cycle: %s" % (node, cycle)
        super().__init__(msg)

class DAG:
    """A direct Acyclic graph where every node has a different value.
    Retrieving any node can be done in constant time from its value."""
    def __init__(self, ):
        self._head_nodes = set()
        self._leaf_nodes = set()
        #Maybe we can do weakrefs in the future, but is not nice to test them
        #without support for basic types such as str and int.
        #We mostly care about function keys, which there isn't too much sense to
        #weakref anyway.
        self._node_refs = {} #weakref.WeakKeyDictionary()

    def to_node(self, value):
        return value if isinstance(value, Node) else self[value]

    def to_nodes(self, values):
        if values is None:
            return set()
        return {self.to_node(val) for val in values}

    def add_node(self, value, inputs=None, outputs=None):
        if value in self:
            raise ValueError("Value already included in graph: %s" % value)

        inputs, outputs = self.to_nodes(inputs), self.to_nodes(outputs)
        n = Node(value, inputs, outputs)
        self._node_refs[value] = n
        try:
            self._wire_node(n)
        except CycleError:
            self.delete_node(n)
            raise


    def _wire_node(self, n):
        """Update data structures taking into account the new state of
        node ``n``."""

        if not n.inputs:
            self._head_nodes.add(n)
        else:
            for parent in n.inputs:
                parent.outputs.add(n)
            self._head_nodes.discard(n)

        if not n.outputs:
            self._leaf_nodes.add(n)
        else:
            for child in n.outputs:
                child.inputs.add(n)
            self._leaf_nodes.discard(n)

        self._leaf_nodes -=  n.inputs
        self._head_nodes -=  n.outputs
        #Check for cycles
        #if we have no inputs or no outputs we cannot create a cycle.
        if n.inputs and n.outputs:
            visited = set()
            for o in n.outputs:
                now =[]
                for u in self.deepfirst_iter(o, visited):
                    now.append(u)
                    if u == n:
                        raise CycleError(n, now)
                now = []


    def add_or_update_node(self, value, inputs=None, outputs=None):
        """Create the node if it doesn't exists. Else add the new inputs and
        outputs to the existing node."""
        if value not in self._node_refs:
            self.add_node(value, inputs, outputs)
        else:
            n = self._node_refs[value]
            inputs, outputs = self.to_nodes(inputs), self.to_nodes(outputs)
            newinputs = inputs - n.inputs
            newoutputs = outputs - n.outputs
            #It's much easier to fail here than to restore the state after
            #_wire_node
            if n in inputs or n in outputs:
                raise CycleError(n,[n])
            n.inputs |= newinputs
            n.outputs |= newoutputs
            try:
                self._wire_node(n)
            except CycleError:
                self.delete_node(n)
                n.inputs -= newinputs
                n.outputs -= newoutputs
                self.add_node(n.value, inputs=n.inputs, outputs=n.outputs)
                raise

    def delete_node(self, n):
        del self._node_refs[n.value]
        self._head_nodes -= {n}
        self._leaf_nodes -= {n}
        for parent in n.inputs:
            parent.outputs.remove(n)
            if not parent.outputs:
                self._leaf_nodes.add(parent)

        for child in n.outputs:
           child.inputs.remove(n)
           if not child.inputs:
               self._head_nodes.add(child)



    def dependency_resolver(self):
        """Yield the nodes that have all dependencies satisfied. Send the next
        completed task."""

        can_run = {n.value for n in self._head_nodes}

        blocked = {output: len(output.inputs) for node in self
                   for output in node.outputs if output.inputs}

        pending = set()
        while True:
            pending |= can_run
            next_completed = yield can_run
            try:
                pending.remove(next_completed)
            except KeyError:
                raise ValueError("Sent value must be pending")
            next_completed = self.to_node(next_completed)

            if not blocked:
                break

            can_run = set()

            for output in next_completed.outputs:
                blocked[output] -= 1
                if blocked[output] == 0:
                    blocked.pop(output)
                    can_run.add(output.value)

    def topological_iter(self):
       """Simplified version of dependency resolver. Yield nodes in such an
       order than dependencies are resolved when actions are executed
       sequentially."""
       can_run = deque(self._head_nodes)

       blocked = {output: len(output.inputs) for node in self.deepfirst_iter() for output in node.outputs}

       while can_run:
           next_node = can_run.popleft()
           yield next_node
           for output in next_node.outputs:
               blocked[output] -= 1
               if blocked[output] == 0:
                   blocked.pop(output)
                   can_run.append(output)

    def deepfirst_iter(self, heads=None, visited=None):
        if heads is None:
            heads = self._head_nodes
        elif isinstance(heads, Node):
            heads = {heads}

        if visited is None:
            visited = set()
        for head in heads:
            if head not in visited:
                yield head
                visited.add(head)
                yield from self.deepfirst_iter(heads=head.outputs,
                                            visited=visited)

    def deepfirst_iter_back(self, leafs=None, visited=None):
        if leafs is None:
            leafs = self._leaf_nodes
        elif isinstance(leafs, Node):
            leafs = {leafs}

        if visited is None:
            visited = set()
        for leaf in leafs:
            if leaf not in visited:
                yield leaf
                visited.add(leaf)
                yield from self.deepfirst_iter_back(leafs=leaf.inputs,
                                                 visited=visited)

    #While this is not recursive, we keep the same interface
    def breadthfirst_iter(self, heads=None, visited=None):

        if heads is None:
            heads = self._head_nodes
        elif isinstance(heads, Node):
            heads = {heads}

        if visited is None:
            visited = set()

        l = deque(heads)

        while l:
            node = l.popleft()
            if not node in visited:
                yield node
                l.extend(node.outputs)
                visited.add(node)

    def breadthfirst_iter_back(self, leafs=None, visited=None):

        if leafs is None:
            leafs = self._leaf_nodes
        elif isinstance(leafs, Node):
            leafs = {leafs}

        if visited is None:
            visited = set()

        l = deque(leafs)

        while l:
            node = l.popleft()
            if not node in visited:
                yield node
                l.extend(node.inputs)
                visited.add(node)

    def __getitem__(self, value):
        return self._node_refs[value]

    def __len__(self):
        return len(self._node_refs)


    def __iter__(self):
        yield from self.topological_iter()

    def __contains__(self, value):
        return value in self._node_refs
