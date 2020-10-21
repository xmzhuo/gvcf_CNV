# make 2D CNN model
# 1D CNN with weight
import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.layers.experimental.preprocessing import CenterCrop
from tensorflow.keras.layers.experimental.preprocessing import Rescaling
#!pip install -q -U keras-tuner
#import kerastuner
import pandas as pd
from keras_preprocessing.image import ImageDataGenerator

#https://www.tensorflow.org/tutorials/structured_data/imbalanced_data

METRICS = [
      keras.metrics.TruePositives(name='tp'),
      keras.metrics.FalsePositives(name='fp'),
      keras.metrics.TrueNegatives(name='tn'),
      keras.metrics.FalseNegatives(name='fn'), 
      keras.metrics.BinaryAccuracy(name='accuracy'),
      keras.metrics.Precision(name='precision'),
      keras.metrics.Recall(name='recall'),
      keras.metrics.AUC(name='auc'),
]

def make_model(metrics = METRICS, output_bias=None):
    if output_bias is not None:
        output_bias = tf.keras.initializers.Constant(output_bias)
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
    outputs = layers.Dense(num_classes, activation='sigmoid',bias_initializer=output_bias)(x)
    model = keras.Model(inputs=inputs, outputs=outputs)
    model.summary()
    model.compile(optimizer=keras.optimizers.Adam(lr=1e-3), 
                        loss=keras.losses.BinaryCrossentropy(),
                        metrics=metrics)
    return model


df = pd.read_csv("/home/jovyan/work/CNV/1Dimg/train/train.csv")
df["labels"]=df["labels"].apply(lambda x:x.split(","))
ct_del = df[[ True if 'del' in x else False for x in df.labels]].shape[0]
ct_dup = df[[ True if 'dup' in x else False for x in df.labels]].shape[0]
ct_nor = df[[ True if 'nor' in x else False for x in df.labels]].shape[0]
wt_del = (ct_del + ct_dup + ct_nor) / ct_del / 3.0
wt_dup = (ct_del + ct_dup + ct_nor) / ct_dup / 3.0
wt_nor = (ct_del + ct_dup + ct_nor) / ct_nor / 3.0
print('del', wt_del, 'dup' , wt_dup, 'nor', wt_nor)
#class_weight = {'del':wt_del, 'dup':wt_dup, 'nor':wt_nor}
class_weight = {0:wt_del, 1:wt_dup, 2:wt_nor}
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

early_stopping = tf.keras.callbacks.EarlyStopping(
    monitor='val_auc', 
    verbose=1,
    patience=10,
    mode='max',
    restore_best_weights=True)

callbacks = [
keras.callbacks.ModelCheckpoint(
filepath='/home/jovyan/work/CNV/CNNmlmodel/model_wt_{epoch}',
save_freq='epoch'),
early_stopping
]
'''
model.fit_generator(generator=train_generator,
                    steps_per_epoch=STEP_SIZE_TRAIN,
                    validation_data=valid_generator,
                    validation_steps=STEP_SIZE_VALID,
                    epochs=10,
                    callbacks=callbacks
                    
)
'''


model = make_model()

import tempfile
initial_weights = os.path.join(tempfile.mkdtemp(),'initial_weights')
model.save_weights(initial_weights)
weighted_model.load_weights(initial_weights)

history = model.fit(x=train_generator,
            steps_per_epoch=STEP_SIZE_TRAIN,
            validation_data=valid_generator,
            validation_steps=STEP_SIZE_VALID,
            epochs=5,
            callbacks=callbacks,
            class_weight=class_weight              
)

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
mpl.rcParams['figure.figsize'] = (12, 10)
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
def plot_metrics(history):
  metrics =  ['loss', 'auc', 'precision', 'recall']
  for n, metric in enumerate(metrics):
    name = metric.replace("_"," ").capitalize()
    plt.subplot(2,2,n+1)
    plt.plot(history.epoch,  history.history[metric], color=colors[0], label='Train')
    plt.plot(history.epoch, history.history['val_'+metric],
             color=colors[0], linestyle="--", label='Val')
    plt.xlabel('Epoch')
    plt.ylabel(name)
    if metric == 'loss':
      plt.ylim([0, plt.ylim()[1]])
    elif metric == 'auc':
      plt.ylim([0.6,1])
    else:
      plt.ylim([0,1])

    plt.legend()

plot_metrics(history)

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