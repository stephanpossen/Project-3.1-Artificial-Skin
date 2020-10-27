"""
Class representing the entire simulation.
It contains the physical system, the sheets, the sensors and those objectes
needed to manage them.
Moreover it uncludes the built-in visualization
"""
from typing import Tuple, Union, List, Dict, Any, Optional  # used for type hints

from SiliconeSheet import Sheet
from PressureSensor import PressureSensor
import SimulationSettings

import pychrono as chrono
import pychrono.fea as fea
import pychrono.irrlicht as irr
import numpy as np


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
                self._system.DoStepDynamics(0.02)
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

    def _make_sheets(self):
        pass

    def _make_sensors(self):
        pass

    def _add_force(self):
        pass

    def computeData(self):
        pass
