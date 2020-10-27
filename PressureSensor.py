import pychrono as chrono

import numpy as np

from typing import Tuple, Union, List, Optional


class PressureSensor:
    """
    Sensor that records the pressure that is applied to the skin at its position.
    It is a single object, representing the SingleTact sensor and its physical characteristics
    Reference: www.SingleTact.com
    Material: polyimide
    Weight sensor: 0.23g
    For this project, the tail of the sensor is ignored.

    Args:
        Size: dimensions of the sensor as (diameter, height) in mm
        Mass: mass of the sensor
        Position: position of the sensor in the sheet
        Min_detectable_force (grams): The minimum amount of force that is needed for the sensor to be able to detect it

    """

    def __init__(self, position: Tuple[float, float, float]):
        # The mass and size of the sensor are known and don't need to be changed at any time
        self._position = position
        self._diameter = 8  # in mm
        self._height = 0.35  # in mm
        self._mass = 0.23  # in grams
        self._density = 1.23  # in g/cm^-3, density of polyimide
        self._contact_material = chrono.ChMaterialSurfaceSMC
        self._min_detectable_force = 20  # in grams

        self._body = chrono.ChBodyEasyCylinder(self._diameter, self._height, self._density, False, True, self._contact_material)

    @property
    def body(self) -> chrono.ChBodyEasyCylinder:
        """
        Returns:
             pychrono.ChBodyEasyCilinder used for this project

        """
        return self._body

    def build(self) -> None:
        """

        :param self:
        :return:
        """
        self._body.SetPos(chrono.ChVectorD(*self._position))
        self._body.SetMass(self._mass)

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

    def getAppliedForce(self):
        return self._body.GetAppliedForce()

    def getContactForce(self):
        return self._body.GetContactForce()



    

