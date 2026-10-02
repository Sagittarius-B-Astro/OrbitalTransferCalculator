"""Orbit description and perifocal -> inertial (ECI) conversion."""
from dataclasses import dataclass
import numpy as np


def Rz(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])


def Rx(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[1.0, 0, 0], [0, c, -s], [0, s, c]])


@dataclass(frozen=True)
class Orbit:
    """Keplerian orbit from apoapsis/periapsis radii and orientation angles (radians)."""
    ra: float
    rp: float
    inc: float = 0.0
    raan: float = 0.0
    argp: float = 0.0

    @property
    def a(self):
        return (self.ra + self.rp) / 2

    @property
    def e(self):
        return (self.ra - self.rp) / (self.ra + self.rp)

    @property
    def p(self):
        return self.a * (1 - self.e ** 2)

    @property
    def dcm(self):
        """Rotation matrix perifocal -> inertial."""
        return Rz(self.raan) @ Rx(self.inc) @ Rz(self.argp)

    def normal(self):
        return self.dcm @ np.array([0.0, 0.0, 1.0])


def state_at(orbit, nu, mu):
    """Inertial position/velocity on `orbit` at true anomaly `nu` (radians)."""
    p, e = orbit.p, orbit.e
    r = p / (1 + e * np.cos(nu))
    rpf = np.array([r * np.cos(nu), r * np.sin(nu), 0.0])
    vpf = np.sqrt(mu / p) * np.array([-np.sin(nu), e + np.cos(nu), 0.0])
    Q = orbit.dcm
    return Q @ rpf, Q @ vpf
