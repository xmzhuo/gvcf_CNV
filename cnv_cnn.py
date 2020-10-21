import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.layers.experimental.preprocessing import CenterCrop
from tensorflow.keras.layers.experimental.preprocessing import Rescaling

#build a CNN 
dense = keras.layers.Dense(units=16)
# Let's say we expect our inputs to be RGB images of arbitrary size
inputs = keras.Input(shape=(None, None, 3))

from tensorflow.keras import layers
# Center-crop images to 101 x 200 to fit sample size
x = CenterCrop(height=101, width=200)(inputs)
# Rescale images to [0, 1]
x = Rescaling(scale=1./255)(x)
# Apply some convolution and pooling layers
x = layers.Conv2D(filters=32, kernel_size=(3, 3), padding='SAME', activation='relu')(x)
x = layers.MaxPooling2D(pool_size=(3, 3))(x)
x = layers.Conv2D(filters=32, kernel_size=(3, 3), padding='SAME', activation='relu')(x)
# Apply global average pooling to get flat feature vectors
x = layers.GlobalAveragePooling2D()(x)
# add a dense layer
x = layers.Dense(20, activation='relu')(x)
# Add a dense classifier on top
num_classes = 10
outputs = layers.Dense(num_classes, activation='softmax')(x)
model = keras.Model(inputs=inputs, outputs=outputs)
model.summary()
#compile and keep metrics
model.compile(optimizer='adam', loss='sparse_categorical_crossentropy',metrics=[keras.metrics.SparseCategoricalAccuracy(name='acc')])
#model.fit(x_train, y_train, batch_size=32, epochs=10)
batch_size = 20
dataset = tf.data.Dataset.from_tensor_slices((x_train, y_train)).batch(batch_size)
print('Fit on Dataset')
#history = model.fit(dataset, epochs=1)
#validation
val_dataset = tf.data.Dataset.from_tensor_slices((x_test, y_test)).batch(batch_size)
history = model.fit(dataset, epochs=1, validation_data=val_dataset)

callbacks = [
keras.callbacks.ModelCheckpoint(
filepath='/home/jovyan/work/CNV/mlmodel/model_{epoch}',
save_freq='epoch')
]
history = model.fit(dataset, epochs=2, callbacks=callbacks)

# Load the TensorBoard notebook extension
%load_ext tensorboard
# Create some TensorBoard logs
callbacks = [
keras.callbacks.TensorBoard(log_dir='./logs'),
keras.callbacks.ModelCheckpoint(
filepath='/home/jovyan/work/CNV/mlmodel/model_{epoch}',
save_freq='epoch')
]
model.fit(dataset, epochs=2, callbacks=callbacks)
# Launch in-line TensorBoard
%tensorboard --logdir logs --host localhost

loss, acc = model.evaluate(val_dataset) # returns loss and metrics
print('loss: %.2f' % loss)
print('acc: %.2f' % acc)

predictions = model.predict(val_dataset)
print(predictions.shape)

import kerastuner
tuner = kerastuner.tuners.Hyperband(
build_model,
objective='val_loss',
max_epochs=100,
max_trials=200,
executions_per_trial=2,
directory='my_dir')
tuner.search(dataset, validation_data=val_dataset)
tuner.results_summary()
models = tuner.get_best_models(num_models=2)

model.save('path/to/location')
model = keras.models.load_model('path/to/location')