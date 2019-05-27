import sklearn.datasets
import sklearn.ensemble
import sklearn.svm
import sklearn.gaussian_process

from reportengine.app import App
from reportengine.report import Config
from reportengine.configparser import ConfigError, element_of
from reportengine.figure import figure
from reportengine.table import table


ALGORITHMS = {
    'svc': sklearn.svm.SVC,
    'gp': sklearn.gaussian_process.GaussianProcessClassifier,
    'random_forest': sklearn.ensemble.RandomForestClassifier
}

class FlowersConfig(Config):
    def produce_dataset(self):
        """Produce the IRIS dataset"""
        return sklearn.datasets.load_iris()

    def parse_xaxis(self, feature:str):
        """Feature in the X axis"""
        return feature

    def parse_yaxis(self, feature:str):
        """Feature in the Y axis"""
        return feature

    def produce_scatterdata(self, dataset, xaxis, yaxis):
        """Check that features exist and return the data values"""
        feat = dataset.feature_names
        def ind(name):
            try:
                 return feat.index(xaxis)
            except ValueError as e:
                raise ConfigError(f'No such feature', bad_item=xaxis,
                    alternatives=feat)
        i,j = ind(xaxis), ind(yaxis)
        return dataset.data[:, i], dataset.data[:,j], dataset.target

    @element_of('algorithms')
    def parse_algorithm(self, alg:str):
        """The name of the classifier algorithm"""
        try:
            return ALGORITHMS[alg]
        except KeyError as e:
            raise ConfigError("Unknown algorithm", bad_item=alg,
                    alternatives=ALGORITHMS) from e

class FlowersApp(App):
    config_class = FlowersConfig


def main():
    a = FlowersApp(name="Flowers",
        default_providers=['flowers.actions', 'reportengine.report'])
    a.main()

if __name__ == '__main__':
    main()

