"""
Classes to model a deformable sheet
"""
# Notes: without material the simulation chrashes

from typing import Tuple, Union, List, Optional  # used for type hints

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
        sphere_swept_thickness: radius used during the collision detection
    """

    def __init__(
        self,
        pos: Tuple[float, float, float],
        size: Tuple[float, float, float],
        subdiv: Union[int, Tuple[int, int, int]],
        mass: float = 0,
        sphere_swept_thickness: Optional[float] = 0.2,
    ) -> None:

        # use the same subdivision on all axis
        if not isinstance(subdiv, tuple):
            subdiv = (subdiv,) * 3

        self._nodes = None  # store nodes in a 3d list, indeces _nodes[x][y][z]
        self._elements = None
        self._sphere_swept_thickness = sphere_swept_thickness
        self._pos = pos
        self._size = size
        self._subdiv = subdiv
        self._mass = mass
        self._material = fea.ChContinuumElastic()  # elastic: no permanet deformation
        self._contact_surface = None  # fea.ChContactSurfaceMesh()
        self._contact_material = chrono.ChMaterialSurfaceSMC()
        self._mesh = fea.ChMesh()

    @property
    def mesh(self) -> fea.ChMesh:
        """
        Returns:
            pychrono.fea.ChMesh used for this object

        """
        return self._mesh

    def build(self) -> None:
        """
        Assemble the different components together.
        Call this method after setting the material properties and before using
        this object.

        Returns:
            None
        """

        # this method is necessary because the contact surface and the elements
        # need a material

        self._make_nodes()
        self._make_elements()

        self._contact_surface = fea.ChContactSurfaceMesh(self._contact_material)
        self._contact_surface.SetMesh(self._mesh)
        self._mesh.AddContactSurface(self._contact_surface)
        self._contact_surface.AddFacesFromBoundary(self._sphere_swept_thickness)

    def set_elastic_properies(
        self,
        young: Optional[float] = None,
        poisson: Optional[float] = None,
        shear_modulus: Optional[float] = None,
        density: Optional[float] = None,
        rayleight_damping_m: Optional[float] = None,
        rayleight_damping_k: Optional[float] = None,
    ) -> None:
        """
        The material is assumed to be elastic (all deformations are reversible)
        The properties of the material are related to the stress (internal force
        generated when a force is applied to the material) and the strain
        (deformation due to the stress).

        Given a planar section of the material, the normal stress is
        :σ = F_normal / Area:
        The average shear stress (parallel to the section) is :τ = Force / Area:

        The linear strain is the ratio of the change in length and the original length.
        :ε = ΔL / L:
        The shear strain  γ is the change of the angle.

        Then the Young modulus E is the ralation between the normal stress and
        the linear strain (Hooke's law).
        :σ = E ε:

        The shear modulus G is the relation between the shear stress and the
        shear strain.
        :τ = G γ:

        "Poisson's ratio ν (nu) is a measure of the Poisson effect,
        that describes the expansion or contraction of a material in
        directions perpendicular to the direction of loading"
        (https://en.wikipedia.org/wiki/Poisson%27s_ratio).

        It is defined as :ν = ε lat / ε axial:, where ε axial is the strain
        in the direction of the applied force and ε lat is the strain perpendiculat to it.
        The value is between 0 and 0.5, where 0.5 is an incopressible material.

        The density is the ratio between the mass and the volume of the material.
        :d = m / V:

        The damping effect of the material is described by the mass-proportional
        damping factor and by the stiffness-proportional damping factor. (# TODO document better)



        If a parameter is set to None then the default Chrono::Engine is used.

        Args:
            young: Young modulus E (in Pa, Pascal)
            poisson: Poisson ratio ν (betwenn 0 and 0.5 included)
            shear_modulus: Shear modulus G (in Pa, Pascal)
            density: density of the material (in kg/m^2)
            rayleight_damping_m: mass-proportional damping factor
            rayleight_damping_k: stiffness-proportional damping factor

        Returns:
            None
        """
        # list of method used to set the values and the values to be set
        setter = {
            self._material.Set_E: young,
            self._material.Set_v: poisson,
            self._material.Set_G: shear_modulus,
            self._material.Set_density: density,
            self._material.Set_RayleighDampingK: rayleight_damping_k,
            self._material.Set_RayleighDamping: rayleight_damping_m,
        }

        # set the value if it is not None
        for s, d in setter.items():
            if d is not None:
                s(d)

    def set_contact_properties(
        self,
        young: Optional[float] = None,
        poisson: Optional[float] = None,
        static_friction: Optional[float] = None,
        kinetic__friction: Optional[float] = None,
        rolling_friction: Optional[float] = None,
        spinning_friction: Optional[float] = None,
        restitution: Optional[float] = None,
        adhesion: Optional[float] = None,
        adhesoinMultDMT: Optional[float] = None,
        adhesionSPerko: Optional[float] = None,
        kn: Optional[float] = None,
        kt: Optional[float] = None,
        gn: Optional[float] = None,
        gt: Optional[float] = None,
    ) -> None:
        """
        Set the parameter of the material of contact surface, i.e. the behaviour
        during a collision.

        If a value is None then the default Chrono::Engine is used

        Args:
            young: see :set_elastic_properies:
            poisson: set_elastic_properies
            static_friction: friction coefficient when the object is not moving yet
            kinetic__friction: friction coeffiction when the object is moving
            rolling_friction: friction coefficient when the object is rolling
            spinning_friction: friction coefficient when the object is spinning
            restitution: ratio between the velocity after and before the collision (between 0 and 1)
            adhesion: adesion force
            adhesoinMultDMT:
            adhesionSPerko:
            kn: normal stiffness coefficient
            kt: tangential stiffness coefficient
            gn: normal damping coefficient
            gt: tangential damping coefficient

        Returns:
            None

        """

        # pairs of values to set and setter methods
        setter = {
            self._contact_material.SetYoungModulus: young,
            self._contact_material.SetPoissonRatio: poisson,
            self._contact_material.SetSfriction: static_friction,
            self._contact_material.SetKfriction: kinetic__friction,
            self._contact_material.SetRollingFriction: rolling_friction,
            self._contact_material.SetSpinningFriction: spinning_friction,
            self._contact_material.SetRestitution: restitution,
            self._contact_material.SetAdhesion: adhesion,
            self._contact_material.SetAdhesionMultDMT: adhesoinMultDMT,
            self._contact_material.SetAdhesionSPerko: adhesionSPerko,
            self._contact_material.SetKn: kn,
            self._contact_material.SetKt: kt,
            self._contact_material.SetGn: gn,
            self._contact_material.SetGt: gt,
        }

        # set the value only if it s not None
        for s, d in setter.items():
            if d is not None:
                s(d)

    def sef_fix(
        self,
        x: Optional[int] = None,
        y: Optional[int] = None,
        z: Optional[int] = None,
        fix=True,
    ) -> None:
        """
        Set all the nodes in a slice as fixed, so that they cannot move

        Args:
            x: x coord to fix, if None use all
            y: y coord to fix, if None use all
            z: z coord to fix, if None use all
            fix: whether to fix or unfix

        Returns:
            None.

        """
        with np.nditer(
            self._nodes[x, y, z],
            flags=["multi_index", "refs_ok"],
            op_flags=["readwrite"],
        ) as it:
            for i in it:
                i.item().SetFixed(fix)

    def set_visualization(
        self,
        original_wireframe: bool = False,
        display_nodes: bool = False,
        node_thikness: float = 0.015,
        display_surf: bool = False,
        smooth_surf: bool = True,
        display_contact: bool = False,
    ) -> None:
        """
        Set how the mesh will be displayed in the simulation.

        Multiple calls to this methods does not overwrite the previous visualization,
        but displays all the visualizations set

        Based on demo_FEA_brick.py by Simone Benatti

        Args:
            original_wireframe: display the original structure as wireframe
            display_nodes: display a cube on each node
            node_thikness: size of the (visual) node
            display_surf: display the surface of the mesh
            smooth_surf: use smooth shading for the surface
            display_contact: display the contact mesh
        Returns:
            None

        """
        # TODO show the wireframe alone

        # Original ref
        if original_wireframe:
            visualization = fea.ChVisualizationFEAmesh(self._mesh)
            visualization.SetFEMdataType(fea.ChVisualizationFEAmesh.E_PLOT_SURFACE)
            visualization.SetWireframe(True)
            visualization.SetDrawInUndeformedReference(True)
            self._mesh.AddAsset(visualization)

        # Node Position
        if display_nodes:
            visualization = fea.ChVisualizationFEAmesh(self._mesh)
            visualization.SetFEMglyphType(
                fea.ChVisualizationFEAmesh.E_GLYPH_NODE_DOT_POS
            )
            visualization.SetFEMdataType(fea.ChVisualizationFEAmesh.E_PLOT_SURFACE)
            visualization.SetSymbolsThickness(node_thikness)
            self._mesh.AddAsset(visualization)

        # Surface
        if display_surf:
            visualization = fea.ChVisualizationFEAmesh(self._mesh)
            visualization.SetFEMdataType(fea.ChVisualizationFEAmesh.E_PLOT_NODE_P)
            visualization.SetSmoothFaces(smooth_surf)
            self._mesh.AddAsset(visualization)

        # Contact
        if display_contact:
            visualization = fea.ChVisualizationFEAmesh(self._mesh)
            visualization.SetFEMdataType(
                fea.ChVisualizationFEAmesh.E_PLOT_CONTACTSURFACES
            )
            visualization.SetWireframe(True)
            visualization.SetDefaultMeshColor(chrono.ChColor(1, 0.5, 0))
            self._mesh.AddAsset(visualization)

    def _make_nodes(self) -> None:
        """
        Generate the list of nodes

        """
        pos = chrono.ChVectorD(*self._pos)
        mass = (
            self._mass
            / (self._subdiv[0] + 1)
            * (self._subdiv[1] + 1)
            * (self._subdiv[2] + 1)
        )
        # for each dimension generate the list of positions where to generate the nodes
        # the index of the node identify the posiion in the grid
        grid = [
            np.linspace(0, dim, sub + 1) for dim, sub in zip(self._size, self._subdiv)
        ]

        # 3d array containing the nodes
        self._nodes = np.ndarray(
            [i + 1 for i in self._subdiv], dtype=np.dtype(fea.ChNodeFEAxyz)
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

              8------7
             /|     /|
            5------6 |
            | 4----|-3
            |/     |/
            1------2

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
        tetra.SetMaterial(self._material)
        return tetra
