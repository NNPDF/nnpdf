% NNPDF hyperplot report

## Best Setup
{@ best_setup @}

## Optimizers
{@ plot_optimizers @}

{@ with hyperscan::optimizer @}
### {@ optimizer_name @}
{@ plot_learning_rate @}
{@ plot_clipnorm @}
{@ endwith @}

## Initializer
{@ plot_initializer @}

## Epochs
{@ plot_epochs @}

## Number of Layers
{@ plot_number_of_layers @}

## Activation function
{@ plot_activation_per_layer @}

## Results Table
[Detailed hyperopt results]({@ results_table report @})
