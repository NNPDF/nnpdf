from reportengine.figure import figure
import numpy as np
import matplotlib.pyplot as plt
from reportengine.table  import table
import seaborn as sns


# 0 = normal scatter plot, 1 = violin, 2 = log
plotting_styles = {
    "iteration": 0,
    "optimizer": 1,
    "learning_rate": 2,
    "initializer": 1,
    "epochs": 0,
    "stopping_epochs": 0,
    "stopping_patience": 0,
    "multiplier": 0,
    "number_of_layers": 1,
    "activation_per_layer": 1,
    "dropout": 0,
    "clipnorm": 0
}


@table
def best_setup(hyperopt_dataframe):
    _, best_trial = hyperopt_dataframe
    return best_trial


@table
def hyperopt_table(hyperopt_dataframe):
    dataframe, _ = hyperopt_dataframe
    return dataframe


@figure
def plot_iterations(hyperopt_dataframe):
    dataframe, best_trial = hyperopt_dataframe
    fig = plot_scans(dataframe, best_trial, "iteration")
    return fig


@figure
def plot_optimizers(hyperopt_dataframe):
    dataframe, best_trial = hyperopt_dataframe
    fig = plot_scans(dataframe, best_trial, "optimizer")
    return fig


@figure
def plot_clipnorm(hyperopt_dataframe, optimizer_name):
    dataframe, best_trial = hyperopt_dataframe
    filtered_dataframe = dataframe[dataframe.optimizer==optimizer_name]
    best_filtered_idx = filtered_dataframe.loss.idxmin()
    best_idx = best_trial.iteration.iloc[0]
    if best_filtered_idx == best_idx:
        include_best = True
    else: 
        include_best = False
    fig = plot_scans(filtered_dataframe, best_trial, "clipnorm", include_best=include_best)
    return fig


@figure
def plot_learning_rate(hyperopt_dataframe, optimizer_name):
    dataframe, best_trial = hyperopt_dataframe
    filtered_dataframe = dataframe[dataframe.optimizer==optimizer_name]
    best_filtered_idx = filtered_dataframe.loss.idxmin()
    best_idx = best_trial.iteration.iloc[0]
    if best_filtered_idx == best_idx:
        include_best = True
    else: 
        include_best = False
    fig = plot_scans(filtered_dataframe, best_trial, "learning_rate", include_best=include_best)
    return fig


@figure
def plot_initializer(hyperopt_dataframe):
    dataframe, best_trial = hyperopt_dataframe
    fig = plot_scans(dataframe, best_trial, "initializer")
    return fig


@figure
def plot_epochs(hyperopt_dataframe):
    dataframe, best_trial = hyperopt_dataframe
    fig = plot_scans(dataframe, best_trial, "epochs")
    return fig


@figure
def plot_number_of_layers(hyperopt_dataframe):
    dataframe, best_trial = hyperopt_dataframe
    fig = plot_scans(dataframe, best_trial, "number_of_layers")
    return fig


@figure
def plot_activation_per_layer(hyperopt_dataframe):
    dataframe, best_trial = hyperopt_dataframe
    fig = plot_scans(dataframe, best_trial, "activation_per_layer")
    return fig


def order_axis(df, bestdf, key):
    """
    Helper function for ordering the axis and make sure the best is always first
    """
    best_x_lst = bestdf.get(key).tolist()
    ordering = set(df.get(key).tolist())
    ordering.remove(best_x_lst[0])
    ordering_true = best_x_lst + list(ordering)
    best_x = np.array([str(best_x_lst[0])])
    return ordering_true, best_x


def plot_scans(df, best_df, plotting_parameter, include_best=True):
    """
    This function plots all trials
    """
    figs, ax = plt.subplots()

    # Set the quantity we will be plotting in the y axis
    loss_k = "loss"

    key = plotting_parameter
    mode = plotting_styles[plotting_parameter]

    if mode == 0 or mode == 2:  # normal scatter plot
        ax = sns.scatterplot(x=key, y=loss_k, data=df, ax=ax)
        best_x = best_df.get(key)
        if mode == 2:
            ax.set_xscale("log")
    elif mode == 1:
        ordering_true, best_x = order_axis(df, best_df, key=key)
        ax = sns.violinplot(
            x=key,
            y=loss_k,
            data=df,
            ax=ax,
            palette="Set2",
            cut=0.0,
            order=ordering_true,
        )
        ax = sns.stripplot(
            x=key,
            y=loss_k,
            data=df,
            ax=ax,
            color="gray",
            order=ordering_true,
            alpha=0.4,
        )

    # Finally plot the "best" one, which will be first
    if include_best:
        ax = sns.scatterplot(
            x=best_x, y=best_df.get(loss_k), ax=ax, color="orange", marker="s"
        )
    ax.set_ylabel("Loss")
    ax.set_xlabel(key)

    return figs
