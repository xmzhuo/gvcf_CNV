#convert csv to jpg , 
# input: the csv has structure of (3 chanel x 5 section (-2000,-1000,0,1000,2000)) x N row
# output: a jpeg file with the same name
import numpy as np
import pandas as pd
import sys
from PIL import Image
import re


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


f1name = sys.argv[1]

df = pd.read_csv(f1name, sep=',', names=["CHROM","POS","DP","GQ","AF"], dtype={"CHROM":str,"POS":str,"DP":float,"GQ":float,"AF":float})
#transfrom relative DP to log2 scale around 100 and amplify it by 50, cap by 0-250 for image, 0.25 > 0; 8 > 250 
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

uniqlist = []
for i in df.CHROM:
    if i not in uniqlist:
        uniqlist.append(i)

#prepare for each chromosome
for i in uniqlist:
    print(i)
    logDPt = list(logDPa[df.CHROM == i])
    logGQt = list(logGQa[df.CHROM == i])
    newAFt = list(newAFa[df.CHROM == i])

    imgsteps=50
    imgstrides=50
    #101 width x length + 2* imgsteps*imgstrides
    my_img1 = movinglist(logDPt,imgsteps,imgstrides)
    my_img2 = movinglist(logGQt,imgsteps,imgstrides)
    my_img3 = movinglist(newAFt,imgsteps,imgstrides)
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

    #for j in range(0,img3.shape[1],int((imgsteps*imgstrides)/2)):
    #    f1newname = re.sub('.csv','.'+str(i)+'_'+str(j)+'-'+str(j+int(2*imgsteps*imgstrides))+'.jpeg',f1name)
    #    im = Image.fromarray(img3[:,j:j+imgstrides*2*2,:], mode='RGB')
    #    im.save(f1newname)
    
    #prepare img like structure array for each range
    imloc = []
    imstart = []
    imend = []
    #type(imloc)
    for j in range(0,img3.shape[1],int((imgsteps*imgstrides)/2)):
        imloc = imloc + [int(j + 1)]
        imstart = imstart + [int(j - 2*imgstrides*imgsteps)]
        imend = imend + [int(j + imgstrides*2*2)]
        imtemp = img3[:,j:j+imgstrides*2*2,:]
        if j == 0 :
            im = np.array([imtemp])
        else :
            #print(newlist2.shape,len(templist))
            im = np.vstack((im,np.array([imtemp])))
    #any hits will be 1 : 2*steps+1 (101 section show positive)
 #docker pull jupyter/tensorflow-notebook
 #docker run --rm -p 10000:8888 -e JUPYTER_ENABLE_LAB=yes -v "$PWD":/home/jovyan/work jupyter/tensorflow-notebook