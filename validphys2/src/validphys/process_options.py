"""
    Module to hold process dependent options

    Only variables included in the `_Vars` enum and processes included in the ``Processes`` enum are allowed
"""
import dataclasses
from enum import Enum
from typing import Callable, Optional, Tuple

from validobj.custom import Parser


class _Vars(Enum):
    x = "x"
    Q2 = "Q2"
    Q = "Q"
    y = "y"


@dataclasses.dataclass(frozen=True)
class _Process:
    name: str
    description: str
    accepted_variables: Tuple[str]
    xq2map_function: Optional[Callable] = None

    def __hash__(self):
        return hash(self.name)

    def xq2map(self, kin_df, metadata):
        """Transform the kinematics dataframe into a x q dataframe
        For double hadronic processes the number of output point will be 2x ninput
        These functions should have access to both the kinematic dataframe and the
        metadata of the commondata
        """
        if self.xq2map_function is None:
            raise NotImplementedError(f"xq2map is not implemented for {self.name}")
        # check that all variables in the dataframe are accepted by this process
        return self.xq2map_function(kin_df, metadata)

    def __str__(self):
        return self.name


def _dis_xq2map(kin_df, metadata):
    """In the old style commondata, the variables in the dataframe were ``x, Q2, y``
    but due to the transformations that happen inside validphys they become ``x, Q, y``
    """
    if all(kin_df.columns == ["k1", "k2", "k3"]):
        # old style format: (x, Q, y)
        x, q, _ = kin_df.values.T
        return x, q * q

    raise NotImplementedError("No new style DIS dataset implemented anyway")


DIS = _Process(
    "DIS",
    "Deep Inelastic Scattering",
    accepted_variables=(_Vars.x, _Vars.Q2, _Vars.y, _Vars.Q),
    xq2map_function=_dis_xq2map,
)


PROCESSES = {
    "DIS": DIS,
    "DIS_NC": dataclasses.replace(DIS, name="DIS_NC"),
    "DIS_CC": dataclasses.replace(DIS, name="DIS_CC"),
    "DIS_NCE": dataclasses.replace(DIS, name="DIS_NCE"),
}


@Parser
def ValidProcess(process_name) -> _Process:
    return PROCESSES.get(process_name.upper(), process_name.upper())
