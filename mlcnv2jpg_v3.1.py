'''
dx run app-cloud_workstation --instance-type mem3_ssd2_x8 --ssh -imax_session_length=3d

##initialized the work station
    unset DX_WORKSPACE_ID
    dx cd $DX_PROJECT_CONTEXT_ID:
    echo $PYTHONPATH

    sudo -E bash
    export PYTHONPATH="/usr/share/dnanexus/lib/python2.7/site-packages:/usr/share/dnanexus/lib/python2.7/site-packages:/usr/share/dnanexus/lib/python2.7/site-packages:/usr/share/dnanexus/lib/python2.7/site-packages:"
    mkdir /mydata
    cd /mydata

docker run --rm -p 32775:8888 -e JUPYTER_ENABLE_LAB=yes -e GRANT_SUDO=yes --user root -v /mydata:/home/jovyan/work jupyter/tensorflow-notebook

hostname -I
dig +short myip.opendns.com @resolver1.opendns.com
18.208.156.157

lsof -ti:32775 | xargs kill -9

myIP=18.208.156.157
server_name=$(echo $myIP | sed 's/\./-/g')

ssh -i ~/.dnanexus_config/ssh_id -o "StrictHostKeyChecking no" -f -L 32775:localhost:32775 -N dnanexus@ec2-$server_name.compute-1.amazonaws.com

'''

#convert csv to jpg , 
# input: the csv has structure of (3 chanel x 5 section (-2000,-1000,0,1000,2000)) x N row
# output: a jpeg file with the same name
import numpy as np
import pandas as pd
import sys
from PIL import Image
import re
#import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.layers.experimental.preprocessing import CenterCrop
from tensorflow.keras.layers.experimental.preprocessing import Rescaling
!pip install -q -U keras-tuner
import kerastuner


fname = sys.argv[1]
#split ratio for train, validation, test
split_data = [6,2,2]

def model_builder(hp):
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
    # Tune the learning rate for the optimizer 
    # Choose an optimal value from 0.01, 0.001, or 0.0001
    hp_learning_rate = hp.Choice('learning_rate', values = [1e-2, 1e-3, 1e-4]) 
    #compile and keep metrics
    model.compile(optimizer=keras.optimizers.Adam(learning_rate = hp_learning_rate), 
                    loss=keras.losses.SparseCategoricalCrossentropy(from_logits = True),
                    metrics=[keras.metrics.SparseCategoricalAccuracy(name='acc')])
    #model.fit(x_train, y_train, batch_size=32, epochs=10)
    return model
batch_size = 20

def movinglist(newlist,steps,strides):
    #movinglist([1]*200,100,50) > 2d array of [100*2+1,200+2*50*100]
    #input: newlist: the list for creating a 2darray of moving
    #       steps: create the width of the expect 2d array, number of moving in either direction
    #       strides: each steps how many value were moved 
    #output: a 2d array of designated size [2*steps+1,len(newlist)+2*strides*steps]
    for i in range(steps*2+1):
        templist = [0]*(2*steps-i)*strides + newlist + [0]*i*strides
        #print(templist)
        if i == 0 :
            newlist2 = np.array([templist])
        else :
            #print(newlist2.shape,len(templist))
            newlist2 = np.vstack((newlist2,np.array([templist])))
    return newlist2


file1 = open(fname, 'r') 
Lines = file1.readlines() 
#tally the samples used for train, val and test set
split_count = [round(len(Lines) * u / sum(split_data)) for u in split_data]
print("Prepare train set :{0}, val set :{1}, test set :{2}".format(*split_count))  
#loop over to add all samples into array for training
n_s = 0
for line in Lines: 
    n_s = n_s + 1
    print("Prepare File{0}: {1}".format(n_s, line.strip())) 
    if n_s <= split_count[0]:
        print("sample {0} is for train set".format(n_s))
    elif n_s <= split_count[0]+split_count[1]:
        print(n_s-1, Xtrain.shape,Ytrain.shape)
        print("sample {0} is for val set".format(n_s))
    else :
        print(n_s-1, Xval.shape,Yval.shape)
        print("sample {0} is for test set".format(n_s))

    f1name = line.strip()
    print(n_s, f1name)
    df = pd.read_csv(f1name, sep=',', names=["CHROM","POS","DP","GQ","AF","CNV"], dtype={"CHROM":str,"POS":str,"DP":float,"GQ":float,"AF":float,"CNV":str})
    #df = pd.read_csv(f1name, sep=',', names=["CHROM","POS","DP","GQ","AF"], dtype={"CHROM":str,"POS":str,"DP":float,"GQ":float,"AF":float})
    #transfrom relative DP to log2 scale around 100 and amplify it by 50, cap by 0-250 for image, 0.25 > 0; 8 > 250 
    df.loc[df.CNV == '.', 'CNV'] = 0
    #df.loc[df.CNV == 'NOR', 'CNV'] = 0
    #df.loc[df.CNV == 'DEL', 'CNV'] = 1
    #df.loc[df.CNV == 'DUP', 'CNV'] = 2
    newCNV = np.array([0] * df.shape[0])
    newCNV[df.CNV == '.'] = 0
    newCNV[df.CNV == 'NOR'] = 0
    newCNV[df.CNV == 'DEL'] = 1
    newCNV[df.CNV == 'DUP'] = 2

    logDP = np.log2(df.DP + 0.0000001)*50+100
    logDP[ logDP < 0 ] = 0
    logDP[ logDP > 250 ] = 250
    logDPa = np.array(list(map(int, logDP)))
    ##transfrom relative GQ to log2 scale around 100 and amplify it by 50, cap by 0-250 for image, 0.25 > 0; 8 > 250 
    logGQ = np.log2(df.GQ + 0.0000001)*50+100
    logGQ[ logGQ < 0 ] = 0
    logGQ[ logGQ > 250 ] = 250
    logGQa = np.array(list(map(int, logGQ)))
    #make AF to be 0-250, 250 for 0.5, 0 for 0 or 1 loss of het
    newAFa = np.array(list(map(int,(0.5-abs(df.AF - 0.5))*500)))

    #check how many chromsome
    uniqlist = []
    for i in df.CHROM:
        if i not in uniqlist:
            uniqlist.append(i)

    n_i=0
    #prepare for each chromosome
    for i in uniqlist:
        n_i = n_i+1
        print(n_i,i)
        logDPt = list(logDPa[df.CHROM == i])
        logGQt = list(logGQa[df.CHROM == i])
        newAFt = list(newAFa[df.CHROM == i])
        CNVt = list(newCNV[df.CHROM == i])   

        imgsteps=50
        imgstrides=50
        #101 width x length + 2* imgsteps*imgstrides
        my_img1 = movinglist(logDPt,imgsteps,imgstrides)
        my_img2 = movinglist(logGQt,imgsteps,imgstrides)
        my_img3 = movinglist(newAFt,imgsteps,imgstrides)
        myCNV = movinglist(CNVt,imgsteps,imgstrides)
        #my_img1 = np.array([[0]*400+logDPt,[0]*300+logDPt+[0]*100,[0]*200+logDPt+[0]*200,[0]*100+logDPt+[0]*300,logDPt+[0]*400])
        #my_img2 = np.array([[0]*400+logGQt,[0]*300+logGQt+[0]*100,[0]*200+logGQt+[0]*200,[0]*100+logGQt+[0]*300,logGQt+[0]*400])
        #my_img3 = np.array([[0]*400+newAFt,[0]*300+newAFt+[0]*100,[0]*200+newAFt+[0]*200,[0]*100+newAFt+[0]*300,newAFt+[0]*400])
        #f1newname = re.sub('.csv','.'+i+'.jpeg',f1name)
        #my_data = np.genfromtxt(f1name, delimiter=',', dtype=int)
        #my_img1 = np.array(my_data[:,[0,3,6,9,12]])
        #my_img2 = np.array(my_data[:,[1,4,7,10,13]])
        #my_img3 = np.array(my_data[:,[2,5,8,11,14]])

        my_img = np.array([my_img1,my_img2,my_img3])
        #my_img.shape (3, 100, 5)
        img2=np.rollaxis(my_img, 0, 3) 
        img3=img2.astype(np.uint8)
        #img2.shape (5,100,3)
        
        #calculate how many DUP or DEL events in each section
        myDEL = np.array([ myCNV[myCNV[:,x] == 1,x].shape[0] for x in range(myCNV.shape[1])])
        myDUP = np.array([ myCNV[myCNV[:,x] == 2,x].shape[0] for x in range(myCNV.shape[1])])

        #for j in range(0,img3.shape[1],int((imgsteps*imgstrides)/2)):
        #    f1newname = re.sub('.csv','.'+str(i)+'_'+str(j)+'-'+str(j+int(2*imgsteps*imgstrides))+'.jpeg',f1name)
        #    im = Image.fromarray(img3[:,j:j+imgstrides*2*2,:], mode='RGB')
        #    im.save(f1newname)
        
        #prepare img like structure array for each range
        imloc = []
        imstart = []
        imend = []

        for j in range(0,img3.shape[1]-int((imgsteps*imgstrides)/2),int((imgsteps*imgstrides)/2)):
            imloc = imloc + [int(j + 1)]
            imstart = imstart + [int(j - 2*imgstrides*imgsteps)]
            imend = imend + [int(j + imgstrides*2*2)]
            imtemp = img3[:,j:j+imgstrides*2*2,:]
            imtempdel = sum(myDEL[j:j+imgstrides*2*2])
            imtempdup = sum(myDUP[j:j+imgstrides*2*2])
            if j == 0 :
                im = np.array([imtemp])
                im_del = [imtempdel]
                im_dup = [imtempdup]
            else :
                #print(newlist2.shape,len(templist))
                im = np.vstack((im,np.array([imtemp])))
                im_del = im_del + [imtempdel]
                im_dup = im_dup + [imtempdup]
            #print(imstart[-1],imend[-1],im.shape,imtempdel,imtempdup)
        im_del = np.array(im_del)
        im_dup = np.array(im_dup)

        imydel=np.array([0]*im_del.shape[0])
        imydel[im_del > 4 ] = 1
        imydup=np.array([0]*im_dup.shape[0])
        imydup[im_dup > 4 ] = 2
        imy = np.array(imydel + imydup) 

        if n_i == 1 :
            X1 = im
            Y1 = imy
        else:
            X1 = np.vstack((X1,im))
            Y1 = np.append(Y1,imy)      
        print(im.shape,X1.shape,imy.shape,Y1.shape)
    
    if n_s == 1 :
        Xtrain = X1
        Ytrain = Y1
    elif n_s <= split_count[0]:
        Xtrain = np.vstack((Xtrain,X1))
        Ytrain = np.append(Ytrain,Y1)
    elif n_s == split_count[0] + 1:
        Xval= X1
        Yval = Y1
    elif n_s <= split_count[0]+split_count[1]:
        Xval = np.vstack((Xval,X1))
        Yval = np.append(Yval,Y1)
    elif n_s == split_count[0] + split_count[1] + 1:
        Xtest= X1
        Ytest = Y1
    else :
        Xtest = np.vstack((Xtest,X1))
        Ytest = np.append(Ytest,Y1)
    
print(n_s, Xtrain.shape,Ytrain.shape, Xval.shape,Yval.shape, Xtest.shape,Ytest.shape)

    #any hits will be 1 : 2*steps+1 (101 section show positive)
 #docker pull jupyter/tensorflow-notebook
 #docker run --rm -p 10000:8888 -e JUPYTER_ENABLE_LAB=yes -v "$PWD":/home/jovyan/work jupyter/tensorflow-notebook
dataset = tf.data.Dataset.from_tensor_slices((X, Y)).batch(batch_size)

class ClearTrainingOutput(tf.keras.callbacks.Callback):
  def on_train_end(*args, **kwargs):
    IPython.display.clear_output(wait = True)

callbacks = [
keras.callbacks.ModelCheckpoint(
filepath='/home/jovyan/work/CNV/mlmodel/model_{epoch}',
save_freq='epoch'),
ClearTrainingOutput()
]
"""
history = model.fit(dataset, epochs=3, callbacks=callbacks)

###
loss, acc = model.evaluate(val_dataset) # returns loss and metrics
print('loss: %.2f' % loss)
print('acc: %.2f' % acc)
"""

tuner = kt.Hyperband(model_builder,
                     objective = 'val_accuracy', 
                     max_epochs = 10,
                     factor = 3,
                     directory = '/home/jovyan/work/CNV/tuner',
                     project_name = 'intro_to_kt')  
tuner.search(dataset, validation_data=val_dataset. epochs=10, callbacks=callbacks)
tuner.results_summary()
# Get the optimal hyperparameters
best_hps = tuner.get_best_hyperparameters(num_trials = 1)[0]

print(f"""
The hyperparameter search is complete. The optimal number of units in the first densely-connected
layer is {best_hps.get('units')} and the optimal learning rate for the optimizer
is {best_hps.get('learning_rate')}.
""")
models = tuner.get_best_models(num_models=2)

model.save('path/to/location')
model = keras.models.load_model('path/to/location')