import tensorflow as tf
from tensorflow.keras.metrics import Metric

import n3fit.backends.keras_backend.operations as op


class LossMetric(Metric):
    """
    Implementation of the (validation) loss as a metric.
    Keeps track of per replica loss internally, aggregates just for logging.

    Parameters
    ----------
        loss_layer : tf.keras.layers.Layer
            The loss layer to use for the metric.
        agg : str
            Aggregation method to use for the replicas. Can be 'sum' or 'mean'.
    """

    def __init__(self, loss_layer, agg='sum', name='val_loss', **kwargs):
        super().__init__(name=name, **kwargs)
        self.loss_layer = loss_layer
        if agg == 'sum':
            self.agg = op.sum
        elif agg == 'mean':
            self.agg = op.mean
        else:
            raise ValueError(f'agg must be sum or mean, got {agg}')
        num_replicas = loss_layer.output.shape[0]
        self.per_replica_losses = self.add_weight(
            name="per_replica_losses", shape=(num_replicas,), initializer="zeros"
        )

    def update_state(self, y_true, y_pred, sample_weight=None):
        self.per_replica_losses.assign(self.loss_layer(y_pred))

    def result(self):
        return self.agg(self.per_replica_losses)

    def reset_state(self):
        self.per_replica_losses.assign(tf.zeros_like(self.per_replica_losses))
