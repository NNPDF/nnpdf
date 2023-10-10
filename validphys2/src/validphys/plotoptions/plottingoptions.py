import dataclasses
import enum
import typing

from validobj import ValidationError

from reportengine.utils import get_functions
from validphys.plotoptions import kintransforms, labelers, resulttransforms
from validphys.plotoptions.utils import get_subclasses

default_labels = ('idat', 'k1', 'k2', 'k3')


labeler_functions = get_functions(labelers)
transform_functions = get_subclasses(kintransforms, kintransforms.Kintransform)
result_functions = get_functions(resulttransforms)


ResultTransformations = enum.Enum('ResultTransformations', list(result_functions.keys()))
TransformFunctions = enum.Enum('TransformFunctions', list(transform_functions.keys()))


class Scale(enum.Enum):
    linear = enum.auto()
    log = enum.auto()
    symlog = enum.auto()


@dataclasses.dataclass
class PlottingOptions:
    func_labels: dict = dataclasses.field(default_factory=dict)
    dataset_label: typing.Optional[str] = None
    experiment: typing.Optional[str] = None
    nnpdf31_process: typing.Optional[str] = None
    data_reference: typing.Optional[str] = None
    theory_reference: typing.Optional[str] = None
    process_description: typing.Optional[str] = None
    y_label: typing.Optional[str] = None
    x_label: typing.Optional[str] = None

    kinematics_override: typing.Optional[TransformFunctions] = None

    result_transform: typing.Optional[ResultTransformations] = None

    # TODO: change this to x: typing.Optional[KinLabel] = None
    # but this currently fails CI because some datasets have
    # a kinlabel of $x_1$ or " "!!
    x: typing.Optional[str] = None

    # TODO: the old commondata uses x, the new uses plot_x
    # the old commondata only allowed the x to be k1, k2, k3 or idat
    plot_x: typing.Optional[str] = None

    x_scale: typing.Optional[Scale] = None
    y_scale: typing.Optional[Scale] = None

    line_by: typing.Optional[list] = None
    figure_by: typing.Optional[list] = None

    extra_labels: typing.Optional[typing.Mapping[str, typing.List]] = None

    def parse_figure_by(self):
        if self.figure_by is not None:
            for el in self.figure_by:
                if el in labeler_functions:
                    self.func_labels[el] = labeler_functions[el]

    def parse_line_by(self):
        if self.line_by is not None:
            for el in self.line_by:
                if el in labeler_functions:
                    self.func_labels[el] = labeler_functions[el]

    def parse_x(self):
        if self.x is not None and self.x not in self.all_labels:
            raise ValidationError(
                f"The label {self.x} is not in the set of known labels {self.all_labels}"
            )

    @property
    def all_labels(self):
        if self.extra_labels is None:
            return set(default_labels)
        return set(self.extra_labels.keys()).union(set(default_labels))

    def __post_init__(self):
        if self.kinematics_override is not None:
            self.kinematics_override = transform_functions[self.kinematics_override.name]()
        if self.result_transform is not None:
            self.result_transform = result_functions[self.result_transform.name]

        # TODO:
        # add a deprecation warning in case anyone start using `x` in the new commondata

        self.parse_figure_by()
        self.parse_line_by()
        self.parse_x()
