"""
The HyperScanner class is basically a dictionary containing all parameters,
the functions defined as hp_ (from hyperspace)

The goal of this module is to read all parameters in the `hyperopt` section of the runcard
and modify the parameter dictionary so that now it is filled with the `hyperop` sampler objects

The idea behind the wrappers if that if you ever want to use another hyperoptimization library, assuming
that it also takes just
    - a function
    - a dictionary of spaces of parameters
you can do so by simply modifying the wrappers to point somewhere else
(and, of course the function in the fitting action that calls the miniimization).
"""
# TODO: make this part of Report Engine
# in principle Report Engine would be reading everything from the runcard and then providing
# basically a list of strings or similar to ModelTrainer, which would in turn make them into 
# the keras/tensorflow/whatever objects it needs.
# If the --hpyeropt flag is active, then report engine should return a sampler from the hyperopt dict
# basically returning the wrapper functions defined below (but instead of simply calling the hyperopt ones
# having done some checks on them)
# even more, depending on the syntax of the fitting part of the runcard, it should do one thing or another
# but the design of this needs someone who is actually good at UX

import hyperopt
import numpy as np
from n3fit.backends import MetaModel, MetaLayer

# These are just wrapper around some hyperopt's sampling expresions defined in here
# https://github.com/hyperopt/hyperopt/wiki/FMin#21-parameter-expressions
# with a bit of extra documentation for the ones that are not obvious
def hp_uniform(key, lower_end, higher_end):
    """ Sample uniformly between lower_end and higher_end """
    return hyperopt.hp.uniform(key, lower_end, higher_end)
def hp_quniform(key, lower_end, higher_end, step_size):
    """ Like uniform but admits a step_size """
    return hyperopt.hp.quniform(key, lower_end, higher_end, step_size)
def hp_loguniform(key, lower_end, higher_end):
    """
    Sample from lower_end to higher_end lograithmically.
    Note that it is different from numpy's logspace in that it takes the
    lower and higher boundaries, not the value of the exponent
    """
    return hyperopt.hp.loguniform(key, lower_end, higher_end)
def hp_choice(key, choices):
    """ Sample from the list `choices` """
    return hyperopt.hp.choice(key, choices)


class HyperScanner:
    """
    The HyperScanner generates a dictionary of parameters for scanning
    It takes cares of known correlation between parameters by tying them together
    It also provides methods for updating the parameter dictionaries after using hyperopt

    It manages:
        stopping
        optimizer
        positivity

        NN architecture (WIP)
    """

    def __init__(self, parameters):
        """
        Takes as input a dictionary of parameters
        and saves it
        """
        self.parameters = parameters
        self.keys = parameters.keys()
        self.hyper_keys = set([])
        self.dict_keys = set([])

    def dict(self):
        return self.parameters

    def update_dict(self, newdict):
        self.parameters.update(newdict)
        # Since some hyperparameters were made in the form of dictionaries, we need to "sviluppare" them
        for key in self.dict_keys:
            if key in newdict.keys():
                value = newdict[key]
                # But this might be a conditional thing and not all of the possible options be dictionaries
                if isinstance(value, dict):
                    self.update_dict(value)

    def _update_param(self, key, value, hyperopt=True, dkeys=None):
        if dkeys is None:
            dkeys = []
        if key not in self.keys:
            raise ValueError(
                "Trying to update a parameter that was not declared in the dictionary: {0}. HyperScanner @ _update_param".format(
                    key
                )
            )
        if hyperopt:
            self.hyper_keys.add(key)
            print("Adding the key {0} with the following value: {1}".format(key, value))
            # Add possible extra keys
            _ = [self.dict_keys.add(nk) for nk in dkeys]
        self.parameters[key] = value

    # Calling the functions below turn parameters of the dictionary into hyperparameter scans
    def stopping(self, min_epochs=5e3, max_epochs=30e3, min_patience=0.10, max_patience=0.3):
        """
        Modifies the following entries of the parameter dictionary:
            - epochs
            - stopping_patience
        """
        epochs_key = "epochs"
        stopping_key = "stopping_patience"

        epoch_step = (max_epochs - min_epochs) / 4
        if min_epochs < epoch_step:
            epoch_step = min_epochs

        epochs = hp_quniform(epochs_key, min_epochs, max_epochs, min_epochs)
        stopping_patience = hp_quniform(stopping_key, min_patience, max_patience, 0.05)

        self._update_param(epochs_key, epochs)
        self._update_param(stopping_key, stopping_patience)

    def optimizer(self, names=None, min_lr=0.0005, max_lr=0.5):
        # But this might be a conditional thing and not all of the possible options be dictionaries
        """
        This function look at the optimizers implemented in MetaModel and adds the learning rate to (only)
        those who use it. The special keyword "ALL" will make it use all optimizers implemented
            - optimizer
            - learning rate
        """
        if names is None:
            names = ["RMSprop"]

        opt_key = "optimizer"
        lr_key = "learning_rate"

        optimizer_dict = MetaModel.optimizers
        choices = []

        min_lr_exp = np.log(min_lr)
        max_lr_exp = np.log(max_lr)
        lr_choice = hp_loguniform(lr_key, min_lr_exp, max_lr_exp)

        if names == "ALL":
            names = optimizer_dict.keys()

        for opt_name in names:
            if opt_name not in optimizer_dict.keys():
                # Do we want an early crash or just drop the optimizer silently? We'll see...
                raise NotImplementedError("HyperScanner: Optimizer {0} not implemented in MetaModel.py".format(opt_name))

            args = optimizer_dict[opt_name][1]
            if "lr" in args.keys():
                choices.append({opt_key: opt_name, lr_key: lr_choice})
            else:
                choices.append(opt_name)

        opt_val = hp_choice(opt_key, choices)
        # Tell the HyperScanner this key might contain a dictionary to store it separately
        self._update_param(opt_key, opt_val, dkeys=[opt_key])

    def positivity(self, min_multiplier=1.01, max_multiplier=1.3, min_initial=0.5, max_initial=100):
        """
        Modifies
            - pos_multiplier
            - pos_initial
        """
        mul_key = "pos_multiplier"
        if not min_initial and not max_initial:
            ini_key = None
        else:
            ini_key = "pos_initial"

        mul_val = hp_uniform(mul_key, min_multiplier, max_multiplier)
        self._update_param(mul_key, mul_val)

        if ini_key:
            min_initial = np.log(min_initial)
            max_initial = np.log(max_initial)
            ini_val = hp_loguniform(ini_key, min_initial, max_initial)
            self._update_param(ini_key, ini_val)

    def NN_architecture(
        self,
        n_layers=None,
        max_units=50,
        min_units=5,
        activations=None,
        initializers=None,
        layer_types=None,
        max_drop=0.0,
    ):
        """
        Uses all the given information to generate the parameters for the NN
            - nodes_per_layer
            - activation_per_layer
            - initializer
            - layer_type
            - dropout
        """
        if activations is None:
            activations = ["sigmoid", "tanh"]
        if initializers is None:
            initializers = ["glorot_normal"]
        if n_layers is None:
            n_layers = [1, 2, 5]
        if layer_types is None:
            layer_types = ["dense"]
        # Do we want to fix or generate dinamically the units? for now let's fix
        # there will be time for more sofisticated functions

        act_key = "activation_per_layer"
        nodes_key = "nodes_per_layer"
        # the information contained in the variable is indeed the number of nodes per layer e.g. [5,10,30]
        # but the information the user will be interested in is how many layers are there...

        # Generate the possible activation choices
        act_choices = []
        # Generate a function that will returns as many copies of the str as the number of layers
        # and a linear one at the end
        for afun in activations:
            def activation_str(nla):
                cop = [afun]*(nla-1)
                cop.append("linear")
                return cop
            act_choices.append(activation_str)

        nodes_choices = []
        for n in n_layers:
            units = []
            for i in range(n):
                # We can even play games as lowering the maximum as the number of layers grows
                units_label = "nl{1}:-{0}/{1}".format(i, n)
                units.append(hp_quniform(units_label, min_units, max_units, 5))

            # And then the last one is a dense with 8 units
            units.append(8)
            nodes_choices.append(units)

        act_functions = hp_choice(act_key, act_choices)
        nodes = hp_choice(nodes_key, nodes_choices)

        # Now let's select the initializers looking at the ones implemented in MetaLayer
        ini_key = "initializer"
        imp_inits = MetaLayer.initializers
        imp_init_names = imp_inits.keys()
        if initializers == "ALL":
            initializers = imp_init_names

        ini_choices = []
        for ini_name in initializers:
            if ini_name not in imp_init_names:
                # Do we want an early crash or just drop the optimizer silently? We'll see...
                raise NotImplementedError("HyperScanner: Initializer {0} not implemented in MetaLayer.py".format(ini_name))
            # For now we are going to use always all initializers and with default values
            ini_choices.append(ini_name)

        ini_choice = hp_choice(ini_key, ini_choices)

        # Finally select the layer types
        if layer_types:
            layer_key = "layer_type"
            layer_choices = hp_choice(layer_key, layer_types)
            self._update_param(layer_key, layer_choices)

        # And add the dropout parameter
        drop_key = "dropout"
        n_drops = 3
        drop_step = max_drop / n_drops
        drop_val = hp_quniform(drop_key, 0.0, max_drop, drop_step)

        self._update_param(act_key, act_functions)
        self._update_param(nodes_key, nodes)
        self._update_param(ini_key, ini_choice)
        self._update_param(drop_key, drop_val)
