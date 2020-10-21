# make 2D CNN model
#adjust to generate imgegenerator to allow multi-labels
import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.layers.experimental.preprocessing import CenterCrop
from tensorflow.keras.layers.experimental.preprocessing import Rescaling
!pip install -q -U keras-tuner
import kerastuner


def model_builder(hp):
    #build a CNN 
    dense = keras.layers.Dense(units=16)
    # Let's say we expect our inputs to be RGB images of arbitrary size
    inputs = keras.Input(shape=(None, None, 3))

    from tensorflow.keras import layers
    # Center-crop images to 101 x 200 to fit sample size
    x = CenterCrop(height=201, width=200)(inputs)
    # Rescale images to [0, 1]
    x = Rescaling(scale=1./255)(x)
    # Apply some convolution and pooling layers
    x = layers.Conv2D(filters=32, kernel_size=(3, 3), padding='SAME', activation='relu')(x)
    x = layers.MaxPooling2D(pool_size=(3, 3))(x)
    x = layers.Conv2D(filters=32, kernel_size=(3, 3), padding='SAME', activation='relu')(x)
    # Apply global average pooling to get flat feature vectors
    x = layers.GlobalAveragePooling2D()(x)
    # add a dense layer
    x = layers.Dense(16, activation='relu')(x)
    # Add a dense classifier on top
    num_classes = 3
    outputs = layers.Dense(num_classes, activation='sigmoid')(x)
    model = keras.Model(inputs=inputs, outputs=outputs)
    model.summary()
    # Tune the learning rate for the optimizer 
    # Choose an optimal value from 0.01, 0.001, or 0.0001
    hp_learning_rate = hp.Choice('learning_rate', values = [1e-2, 1e-3, 1e-4]) 
    #compile and keep metrics
    #model.compile(optimizer=keras.optimizers.Adam(learning_rate = hp_learning_rate), 
    #                loss=keras.losses.SparseCategoricalCrossentropy(from_logits = True),
    #                metrics=[keras.metrics.SparseCategoricalAccuracy(name='acc')])
    model.compile(optimizer=keras.optimizers.Adam(learning_rate = hp_learning_rate), 
                    loss=keras.losses.BinaryCrossentropy()(from_logits = True),
                    metrics=['accuracy'])
    #model.fit(x_train, y_train, batch_size=32, epochs=10)
    return model


batch_size=20
dataset = keras.preprocessing.image_dataset_from_directory(
'/home/jovyan/work/CNV/img/train', batch_size=64, image_size=(201, 200), class_names=['nor', 'del','dup'])
val_dataset = keras.preprocessing.image_dataset_from_directory(
'/home/jovyan/work/CNV/img/val', batch_size=64, image_size=(201, 200), class_names=['nor', 'del','dup'])

import pandas as pd
from keras_preprocessing.image import ImageDataGenerator
df = pd.read_csv("/home/jovyan/work/CNV/img/test_unpack/test.csv")
df["labels"]=df["labels"].apply(lambda x:x.split(","))

datagen=ImageDataGenerator(rescale=1./255.)
test_datagen=ImageDataGenerator(rescale=1./255.)
train_generator=datagen.flow_from_dataframe(
dataframe=df,
directory="/home/jovyan/work/CNV/img/test_unpack",
x_col="filename",
y_col="labels",
batch_size=32,
seed=42,
shuffle=True,
class_mode="categorical",
classes=["nor","del","dup"],
target_size=(201,200))
valid_generator=test_datagen.flow_from_dataframe(
dataframe=df,
directory="./miml_dataset/images",
x_col="filename",
y_col="labels",
batch_size=32,
seed=42,
shuffle=True,
class_mode="categorical",
classes=["nor","del","dup"],
target_size=(201,200))
test_generator=test_datagen.flow_from_dataframe(
dataframe=df[1900:],
directory="./miml_dataset/images",
x_col="Filenames",
batch_size=1,
seed=42,
shuffle=False,
class_mode=None,
target_size=(100,100))


# Let's say we expect our inputs to be RGB images of arbitrary size
inputs = keras.Input(shape=(None,None,3))
targetsize=512
from tensorflow.keras import layers
import math
# 1D cropping to fit sample size
x = CenterCrop(height=1,width=targetsize)(inputs)
#print(x.shape)
# Rescale images to [0, 1]
x = Rescaling(scale=1./255)(inputs)
# Apply some convolution and pooling layers
x = layers.Conv2D(filters=32, kernel_size=3, strides=2, padding='SAME', activation='relu')(x)
x = layers.Conv2D(filters=32, kernel_size=3, strides=2, padding='SAME', activation='relu')(x)
x = layers.Conv2D(filters=32, kernel_size=3, strides=2, padding='SAME', activation='relu')(x)
# Apply global average pooling to get flat feature vectors
x = layers.GlobalAveragePooling2D()(x)
# add a dense layer
x = layers.Dense(16, activation='relu')(x)
# Add a dense classifier on top
num_classes = 3
outputs = layers.Dense(num_classes, activation='sigmoid')(x)
model = keras.Model(inputs=inputs, outputs=outputs)
model.summary()
model.compile(optimizer=keras.optimizers.Adam(), 
                    loss=keras.losses.BinaryCrossentropy(),
                    metrics=['accuracy'])

STEP_SIZE_TRAIN=train_generator.n//train_generator.batch_size
STEP_SIZE_VALID=valid_generator.n//valid_generator.batch_size
STEP_SIZE_TEST=test_generator.n//test_generator.batch_size

callbacks = [
keras.callbacks.ModelCheckpoint(
filepath='/home/jovyan/work/CNV/CNNmlmodel/model_{epoch}',
save_freq='epoch')
]
model.fit_generator(generator=train_generator,
                    steps_per_epoch=STEP_SIZE_TRAIN,
                    validation_data=valid_generator,
                    validation_steps=STEP_SIZE_VALID,
                    epochs=10,
                    callbacks=callbacks
                    
)

'''
datagen = tf.keras.preprocessing.image.ImageDataGenerator(rescale=1./255)
train_data = datagen.flow_from_directory(
        '/home/jovyan/work/CNV/img/train',
        target_size=(201, 200),
        batch_size=batch_size,
        class_mode='binary',
        classes=['nor', 'del','dup'])

val_data = datagen.flow_from_directory(
        '/home/jovyan/work/CNV/img/val',
        target_size=(201, 200),
        batch_size=batch_size,
        class_mode='binary',
        classes=['nor', 'del','dup'])
'''
'''
class ClearTrainingOutput(tf.keras.callbacks.Callback):
  def on_train_end(*args, **kwargs):
    IPython.display.clear_output(wait = True)

callbacks = [
keras.callbacks.ModelCheckpoint(
filepath='/home/jovyan/work/CNV/mlmodel/model_{epoch}',
save_freq='epoch'),
ClearTrainingOutput()
'''
callbacks = [
keras.callbacks.ModelCheckpoint(
filepath='/home/jovyan/work/CNV/mlmodel/model_{epoch}',
save_freq='epoch')
]

tuner = kerastuner.Hyperband(model_builder,
                     objective = 'accuracy', 
                     max_epochs = 10,
                     factor = 3,
                     directory = '/home/jovyan/work/tuner',
                     project_name = 'intro_to_kt')  
tuner.search(dataset, validation_data=val_dataset, epochs=10, callbacks=callbacks)
tuner.results_summary()

models = tuner.get_best_models(num_models=2)

models.save('/home/jovyan/work/bestmodel')