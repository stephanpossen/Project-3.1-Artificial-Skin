"""
Class representing the entire simulation.
It contains the physical system, the sheets, the sensors and those objectes
needed to manage them.
Moreover it uncludes the built-in visualization
"""
# Parameters
# General
#    Effective curvature radius = 1 (no idea what it is)
#    Suggested margin = 0.006 (no idea)
#    time step = 0.004  (how much time is simulated at each step)
#    length = 10 (how long is the simulation, in seconds)

# Built-in visualization
#    Data path = "..." (path to the folder with the textures)
#    title = "..." (title of the window)
#    window_width = 1024 (dimensions of the window)
#    window_height = 768
#    fullscreen = False
#    shadows = True
#    antialias = True
#    camera_pos_x = some float (position and rotation of the camera)
#    camera_pos_y =
#    camera_pos_z =
#    camera_rot_x =
#    camera_rot_y =
#    camera_rot_z =

# Solver -> solve differential equations
#    iterations = 1000 (int)
#    tolerance = 1e-12 (float)
#    diagonal preconditioner = True

# Stepper -> time integration
#    alpha = -0.2
#    max iterations = 100
#    absolute tolerance = 1e-5
#    scaling = True

# Silicone sheets   -> 2 sheets
#    pos_x1 = float
#    pos_y1
#    pos_z1
#    size_x1 = float
#    size_y1
#    size_z1
#    subdiv_x1 = int
#    subdiv_y1
#    subdiv_z1
#    pos_x2 = float -> can we assume that they have the same size, subdivisions and mass?
#    pos_y2
#    pos_z2
#    size_x2 = float
#    size_y2
#    size_z2
#    subdiv_x2 = int
#    subdiv_y2
#    subdiv_z2
#    mass1
#    mass2
#    sphere_swept_thickness = float (something for the collisions, but I don't know the details)
#   Elastic
#    young
#    poisson
#    density
#    shear  -> optional
#    rayleight damping m -> optional
#    rayleight damping k -> optional
#   Contact
#    young
#    poisson
#    static friction
#    kinetic friction
#    rolling friction -> optional
#    spinning friction -> optional
#    restitution
#    adhesion
#    adhesionMultDMt, adhesionSPerko, kn, kt, gn, gt (some other parameters probably not used)
#   Visualization
#    original wireframe = False (show initial position of the edges)
#    display nodes = False (show a small cube on each node)
#    node thikness = 0.015 (size of those cube)
#    display surface = False (disaply the surface of the sheets)
#    smooth_surf = True
#    display contact = False

# Pressure sensors
#    n = int (how many sensors)
#    pos_x  (we need the position of each sensor)
#    pos_y
#    pos_z
#    Optional: diameter, height, mass, density, min detectable force, friction, young modulus, adhesions, etc, ...

from typing import Dict, Any, Optional  # used for type hints

from SiliconeSheet import Sheet
from PressureSensor import PressureSensor
import SimulationSettings

import pychrono as chrono
# import pychrono.fea as fea
import pychrono.irrlicht as irr
# import numpy as np


class Simulation:
    """
    This class simulates the artificial skin and export the result.
    It constructs the physical system given a configuation object and runs
    the simulation.

    Multiple simulation can be run one after the other.

    The configuarion contains general information about the simulation itself
    (e.g. the gravity value, the time step, ...) and information of the objects
    of the physical system (e.g. the positions, the sizes, ...).

    Args:
        configuation: configuation object used
        visualization: whether the buil-in simulation is used or not
        output: object receiving the output of the simulation
    Returns:
        None
    """

    def __init__(
        self, configuration, visualization: bool = False, output: Optional = None
    ) -> None:
        # TODO a copy of the configuration is better
        self._conf = configuration
        self._visualization = visualization
        self._output = output
        self._app = None

        # TODO not sure about the meaning of theese two parameters, take them from the configuation
        chrono.ChCollisionInfo.SetDefaultEffectiveCurvatureRadius(1)
        chrono.ChCollisionModel.SetDefaultSuggestedMargin(0.006)

        # TODO look into the parameters of the system
        # e.g. MinBounceSpeed, Gravity
        self._system = chrono.ChSystemSMC()

        if self._visualization:
            # TODO take the data path, the title, the dimension and do_[something] from the configuration
            chrono.SetChronoDataPath(
                "/home/gianluca/anaconda3/envs/chrono/share/chrono/data/"
            )
            self._app = irr.ChIrrApp(
                self._system,
                "Artificial Skin",
                irr.dimension2du(1024, 768),
                do_fullscreen=False,
                do_shadows=True,
                do_antialias=True,
            )
            self._app.AddTypicalSky()
            self._app.AddTypicalLights()
            # TODO take the info for the came and lighs from the configuation
            self._app.AddTypicalCamera(irr.vector3df(-1, -1, 0), irr.vector3df(0, 0, 0))
            # TODO this has too many parameters and it is not that important
            self._app.AddLightWithShadow(
                irr.vector3df(1.5, 5.5, -2.5),
                irr.vector3df(0, 0, 0),
                3,
                2.2,
                7.2,
                40,
                512,
                irr.SColorf(1, 1, 1),
            )
            self._app.AddShadowAll()
            self._app.SetTimestep(0.004)

        self._make_sheets()
        self._make_sensors()
        # TODO make the load class -> it takes a node / pair of indexes and time and return the force
        self._add_forces()

        # TODO look at all the parameters of the solver and the stepper
        # TODO consider also mkl.ChSolverMKL
        # TODO take the parameters from the configuration
        self._solver = chrono.ChSolverMINRES()
        self._system.SetSolver(self._solver)
        self._solver.SetMaxIterations(1000)
        self._solver.SetTolerance(1e-12)
        self._solver.EnableDiagonalPreconditioner(True)
        self._solver.SetVerbose(False)  # don't take this from the configuration

        # HHT implicit integrator for II order systems, adaptive
        # TODO take the parameters from the configuaration
        self._stepper = chrono.ChTimestepperHHT(self._system)
        self._system.SetTimestepper(self._stepper)
        self._stepper.SetAlpha(-0.2)
        self._stepper.SetMaxiters(100)
        self._stepper.SetAbsTolerances(1e-5)
        self._stepper.SetMode(chrono.ChTimestepperHHT.POSITION)
        self._stepper.SetScaling(True)

        # TODO from configuarion
        # length of the simulation
        self._life = 5  # seconds

    def run(self):
        if self._visualization:
            while self._app.GetDevice().run() or self._system.GetChTime() > self._life:
                self._app.BeginScene()
                self._app.DrawAll()
                self.computeData()  # it also export the data if self._output is not none
                self._app.DoStep()
                self._app.EndScene()
        else:
            while self._system.GetChTime() > self._life:
                self.computeData()
                # TODO parameter from the configuration
                self._system.DoStepDynamics(0.04)
        self.computeData()  # final state

    @classmethod
    def run_series(
        cls,
        configurations: Dict[Any, int],
        visualization: bool = False,
        output: Optional = None,
    ) -> None:
        """
        Run multiple simulations one after the other

        Args:
            configurations: dictionary of configuations and the number of simulations for each of them
            visualization: whether the buil-in simulation is used or not
            output: object receiving the output of the simulation
        Returns:
            None
        """
        for conf, n in configurations.items():
            for i in range(n):
                cls(conf, visualization, output())

    def _make_sheets(self) -> None:
        """
        Generate 2 sheets with the parameters given by the configuration

        Returns:
            None.
        """

        self._sheets = [None, None]
        # TODO parameters from the configuration
        pos1 = (0, 0, 0)
        pos2 = (0, 0, 0)
        size = (0, 0, 0)
        subdiv = (0, 0, 0)
        mass = 0
        sphere_swept_thickness = 0

        young = 0
        poisson = 0
        shear = 0
        density = 0
        static_friction = 0
        kinetic_friction = 0
        rolling_friction = 0
        spinning_friction = 0
        restitution = 0
        adhesion = 0

        sheet = Sheet(pos1, size, subdiv, mass, sphere_swept_thickness)
        sheet.set_elastic_properies(young, poisson, shear, density)
        sheet.set_contact_properties(
            young,
            poisson,
            static_friction,
            kinetic_friction,
            rolling_friction,
            spinning_friction,
            restitution,
            adhesion)
        sheet.build()
        self._sheets[0] = sheet

        sheet = Sheet(pos2, size, subdiv, mass, sphere_swept_thickness)
        sheet.set_elastic_properies(young, poisson, shear, density)
        sheet.set_contact_properties(
            young,
            poisson,
            static_friction,
            kinetic_friction,
            rolling_friction,
            spinning_friction,
            restitution,
            adhesion)
        sheet.build()
        self._sheets[1] = sheet

        self._system.AddMesh(self._sheets[0].mesh)
        self._system.AddMesh(self._sheets[1].mesh)

    def _make_sensors(self) -> None:
        """
        Generate a list of sensors in the given positions

        Returns:
            None.

        """
        n = 0
        # TODO consider using a numpy array
        self._sensors = [None] * n
        tmp = None
        for i in range(n):
            # TODO from the configuaration
            pos = (0, 0, 0)
            tmp = PressureSensor(pos)
            # TODO make material
            tmp.build()
            self._sensors[i] = tmp
            self._system.AddBody(tmp)

    def _add_force(self):
        pass

    def computeData(self):
        data = [0] * len(self._sensors)
        if self._output is not None:
            for i, sens in enumerate(self._sensors):
                # TODO use correct method
                data[i] = sens.getAppliedForce()
            # TODO use correct method
            self._output.setData(data)
