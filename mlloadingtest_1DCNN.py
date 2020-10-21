# make 2D CNN model
import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.layers.experimental.preprocessing import CenterCrop
from tensorflow.keras.layers.experimental.preprocessing import Rescaling
!pip install -q -U keras-tuner
import kerastuner


def model_builder(hp):
    #build a CNN 
    # Let's say we expect our inputs to be RGB images of arbitrary size
    inputs = keras.Input(shape=(None, 3))
    targetsize=512
    from tensorflow.keras import layers
    import math
    # 1D cropping to fit sample size
    x = layers.Cropping1D(cropping=[math.floor((inputs.shape[1]-targetsize)/2),math.ceil((inputs.shape[1]-targetsize)/2)])(inputs)
    #print(x.shape)
    # Rescale images to [0, 1]
    x = Rescaling(scale=1./255)(x)
    # Apply some convolution and pooling layers
    x = layers.Conv1D(filters=32, kernel_size=3, strides=2, padding='SAME', activation='relu')(x)
    x = layers.Conv1D(filters=32, kernel_size=3, strides=2, padding='SAME', activation='relu')(x)
    x = layers.Conv1D(filters=32, kernel_size=3, strides=2, padding='SAME', activation='relu')(x)
    # Apply global average pooling to get flat feature vectors
    x = layers.GlobalAveragePooling1D()(x)
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
    model.compile(optimizer=keras.optimizers.Adam(learning_rate = hp_learning_rate), 
                    loss=keras.losses.BinaryCrossentropy(from_logits = True),
                    metrics=[keras.metrics.SparseCategoricalAccuracy(name='acc')])
    #model.fit(x_train, y_train, batch_size=32, epochs=10)
    return model

'''
batch_size=64
dataset = keras.preprocessing.image_dataset_from_directory(
'/home/jovyan/work/CNV/img/train', batch_size=batch_size, image_size=(201, 200), class_names=['nor', 'del','dup'])
val_dataset = keras.preprocessing.image_dataset_from_directory(
'/home/jovyan/work/CNV/img/val', batch_size=batch_size, image_size=(201, 200), class_names=['nor', 'del','dup'])
'''

import pandas as pd
import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.layers.experimental.preprocessing import CenterCrop
from tensorflow.keras.layers.experimental.preprocessing import Rescaling
import kerastuner
from keras_preprocessing.image import ImageDataGenerator

df = pd.read_csv("/home/jovyan/work/CNV/1Dimg/train/train.csv")
df["labels"]=df["labels"].apply(lambda x:x.split(","))

datagen=ImageDataGenerator(rescale=1./255.)
train_generator=datagen.flow_from_dataframe(
dataframe=df,
directory="/home/jovyan/work/CNV/1Dimg/train",
x_col="filename",
y_col="labels",
batch_size=32,
seed=42,
shuffle=True,
class_mode="categorical",
classes=["del","dup","nor"],
target_size=(1,512))

test_datagen=ImageDataGenerator(rescale=1./255.)
df = pd.read_csv("/home/jovyan/work/CNV/1Dimg/val/val.csv")
df["labels"]=df["labels"].apply(lambda x:x.split(","))
valid_generator=test_datagen.flow_from_dataframe(
dataframe=df,
directory="/home/jovyan/work/CNV/1Dimg/val",
x_col="filename",
y_col="labels",
batch_size=32,
seed=42,
shuffle=True,
class_mode="categorical",
classes=["del","dup","nor"],
target_size=(1,512))

df = pd.read_csv("/home/jovyan/work/CNV/1Dimg/test/test.csv")
df["labels"]=df["labels"].apply(lambda x:x.split(","))
test_generator=test_datagen.flow_from_dataframe(
dataframe=df,
directory="/home/jovyan/work/CNV/1Dimg/test",
x_col="filename",
batch_size=1,
seed=42,
shuffle=False,
class_mode=None,
target_size=(1,512))

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


test_generator.reset()
pred=model.predict_generator(test_generator,
steps=STEP_SIZE_TEST,
verbose=1)

#df[[ True if 'dup' in x else False for x in df.labels]]

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
import IPython
class ClearTrainingOutput(tf.keras.callbacks.Callback):
  def on_train_end(*args, **kwargs):
    IPython.display.clear_output(wait = True)
callbacks = [
keras.callbacks.ModelCheckpoint(
filepath='/home/jovyan/work/CNV/mlmodel/model_{epoch}',
save_freq='epoch'),
ClearTrainingOutput()
]

tuner = kerastuner.Hyperband(model_builder,
                     objective = 'val_acc', 
                     max_epochs = 10,
                     factor = 3,
                     directory = '/home/jovyan/work/CNN_tuner',
                     project_name = '1D')  
tuner.search(train_generator, validation_data=valid_generator, epochs=10, callbacks=callbacks)
tuner.results_summary()

models = tuner.get_best_models(num_models=2)

models.save('/home/jovyan/work/bestmodel')