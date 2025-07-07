from dataclasses import dataclass
from enum import Enum

class FlareProfile(str, Enum):
    """Enumeration for the horn flare profile types."""
    EXPONENTIAL = "exponential"
    HYPERBOLIC = "hyperbolic"
    CONICAL = "conical"

@dataclass
class HornParameters:
    """
    A structured representation of a horn's geometric parameters.
    All dimensions are in meters.
    """
    flare_profile: FlareProfile
    throat_radius: float
    mouth_radius: float
    length: float
    
    # Flare constant for exponential or hyperbolic horns
    m: float = 1.0 