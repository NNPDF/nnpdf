"""
    Module to hold process dependent options

    Only variables included in the `_Vars` enum and processes included in the ``Processes`` enum are allowed
"""
import dataclasses
from enum import Enum
from typing import Callable, Optional, Tuple


class _Vars(Enum):
    x = "x"
    Q2 = "Q2"
    y = "y"


@dataclasses.dataclass(frozen=True)
class _Process:
    name: str
    description: str
    accepted_variables: Tuple[str]
    xq2map_function: Optional[Callable] = None

    def __hash__(self):
        return hash(self.name)

    def xq2map(self, kin_df):
        """Transform the kinematics dataframe into a x q dataframe
        For double hadronic processes the number of output point will be 2x ninput
        """
        if self.xq2map_function is None:
            raise NotImplementedError(f"xq2map is not implemented for {self.name}")
        # check that all variables in the dataframe are accepted by this process
        return self.xq2map(kin_df)


def _dis_xq2map(kin_df):
    if all(kin_df.columns == ["kin1", "kin2", "kin3"]):
        # old style format: (x, Q2, y)
        return kin_df.drop("kin3", axis=1)
    raise NotImplementedError("No new style DIS dataset implemented anyway")


DIS = _Process(
    "DIS",
    "Deep Inelastic Scattering",
    accepted_variables=(_Vars.x, _Vars.Q, _Vars.y),
    xq2map_function=_dis_xq2map,
)


class Processes:
    dis = DIS
    dis_nc = dataclasses.replace(DIS, name="DIS NC")
    dis_cc = dataclasses.replace(DIS, name="DIS CC")
