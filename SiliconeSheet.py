"""
Classes to model a deformable sheet
"""

from typing import Tuple, Union, List  # used for type hints

import pychrono as chrono
import pychrono.fea as fea
import numpy as np


class Sheet:
    """
    Deformable parallelepiped for the simulation.
    This setup a mesh made up of tetrahedical elements with the given
    physical properties

    Args:
        pos: position of the sheet. This is also the position of the first corner node generated
        size: dimensions of the sheet as (width, length, height)
        subdiv: the number of times that the original parallelepiped is subdivided on each dimension \
            e.g. [1, 1, 2] means that on the z (vertical) axis the parallelepiped is split in two
        mass: mass of the sheet, assumed to be uniform
    """

    def __init__(
        self,
        pos: Tuple[float, float, float],
        size: Tuple[float, float, float],
        subdiv: Union[int, Tuple[int, int, int]],
        mass: float = 0,
    ) -> None:

        self._nodes = None  # store nodes in a 3d list, indeces _nodes[x][y][z]
        self._elements = None
        self.material = None  # TODO
        self.contact_surface = None  # TODO
        self._mesh = fea.ChMesh()

        # use the same subdivision on all axis
        if not isinstance(subdiv, tuple):
            subdiv = (subdiv,) * 3

        self._make_nodes(pos, size, subdiv, mass)
        self._make_elements()

    @property
    def mesh(self) -> fea.ChMesh:
        """
        Returns:
            pychrono.fea.ChMesh used for this object

        """
        return self._mesh

    def sef_fix(self, faces: List[str] = [], edges: List[str] = [], fix=True) -> None:
        """
        Fix the nodes of some faces and edges of the mesh, i.e. they cannot be moves.

        The faces are identified by the following strings

        * x+: left face
        * x-: right face
        * y+: back face
        * y-: front face
        * z+: top face
        * z-: bottom face

        The edges are identified by the intersection of two faces:
            x+y+ is the edge between the top and the left face
        Faces with no intersections are ignored

        Args:
            faces: list of faces to fix
            edges: list of edges to fix
            fix: whether to fix or unfix the selected nodes
        Raises:
            ValueError: onincorrect strings

        """
        # TODO find a better way to identify faces and edges
        pass

    def set_visualization(self, *args) -> None:
        """
        Set how the mesh will be displayed in the simulation.

        Multiple calls to this methods does not overwrite the previous visualization, but
        but displays all the visualizations set

        Args:

        Returns:
            None

        """
        # TODO

    def _make_nodes(
        self,
        pos: Tuple[float, float, float],
        size: Tuple[float, float, float],
        subdiv: Tuple[int, int, int],
        mass,
    ) -> None:
        """
        Generate the list of nodes

        Args:
            pos: position of the mesh w.r.t the left most bottom front node
            size: size of the sheet
            subdiv: subdivisions on each axis
            mass: total mass
        """
        pos = chrono.ChVectorD(*pos)
        mass /= (subdiv[0] + 1) * (subdiv[1] + 1) * (subdiv[2] + 1)
        # for each dimension generate the list of positions where to generate the nodes
        # the index of the node identify the posiion in the grid
        grid = [np.linspace(0, dim, sub + 1) for dim, sub in zip(size, subdiv)]

        # 3d array containing the nodes
        self._nodes = np.ndarray(
            [i + 1 for i in subdiv], dtype=np.dtype(fea.ChNodeFEAxyz)
        )

        # iterator for self._nodes
        with np.nditer(
            self._nodes, flags=["multi_index", "refs_ok"], op_flags=["readwrite"]
        ) as it:
            for i in it:
                ind = it.multi_index  # index of the current node
                node_pos = chrono.ChVectorD(
                    grid[0][ind[0]], grid[1][ind[1]], grid[2][ind[2]]
                )  # local position of the node
                tmp = fea.ChNodeFEAxyz(pos + node_pos)
                tmp.SetMass(mass)
                i[...] = tmp  # add the node to self._nodes
                self._mesh.AddNode(tmp)  # add the node to the mesh

    def _make_elements(self) -> None:
        """
        Generate the tetrahedrons from the list of nodes

        Idea:
            * find a "cube"
            * generate 5 fea.ChElementTetra_4
            * repeat until it traversed the entire grid

            Nodes in the cube:

              8-----7
             /|    /|
            5-----6 |
            | 4---|-3
            |/    |/
            1-----2

            7 has index (i, j, k), 1 has (i-1, j-1, k-1)
        Returns:
            None.

        """
        # number of nodes per dimension
        x, y, z = self._nodes.shape

        # iterators for x, y, z
        i = j = k = 1

        e = 0  # index of the element being created -> easier than retrieving it with i,j,k
        self._elements = [None] * (x * y * z * 5)

        # TODO find a better way to represent the cube
        n1 = n2 = n3 = n4 = n5 = n6 = n7 = n8 = None

        while i < x:
            while j < y:
                while k < z:
                    # update cube
                    n1 = self._nodes[i - 1, j - 1, k - 1]
                    n2 = self._nodes[i, j - 1, k - 1]
                    n3 = self._nodes[i, j, k - 1]
                    n4 = self._nodes[i - 1, j, k - 1]
                    n5 = self._nodes[i - 1, j - 1, k]
                    n6 = self._nodes[i, j - 1, k]
                    n7 = self._nodes[i, j, k]
                    n8 = self._nodes[i - 1, j, k]

                    # tetra 1
                    tmp = self._make_tetra(n1, n2, n3, n6)
                    self._elements[e] = tmp
                    e += 1

                    # tetra 2
                    tmp = self._make_tetra(n1, n3, n4, n8)
                    self._elements[e] = tmp
                    e += 1

                    # tetra 3
                    tmp = self._make_tetra(n1, n3, n6, n8)
                    self._elements[e] = tmp
                    e += 1

                    # tetra 4
                    tmp = self._make_tetra(n1, n5, n6, n8)
                    self._elements[e] = tmp
                    e += 1

                    # tetra 5
                    tmp = self._make_tetra(n3, n6, n7, n8)
                    self._elements[e] = tmp
                    e += 1

                    k += 1
                j += 1
            i += 1

    def _make_tetra(self, nodeA, nodeB, nodeC, nodeD):
        """
        Make single tetrahedron
        """
        tetra = fea.ChElementTetra_4()
        self._mesh.AddElement(tetra)
        tetra.SetNodes(nodeA, nodeB, nodeC, nodeD)
        return tetra
