#calculate the coverage and gq for each interval in a tsv
#normalize the result per sample

import pandas as pd
import numpy as np
import math
import time
import re
import sys
#bed filename
f1name=sys.argv[1]
f1name_re=re.sub('.tsv','',f1name)

#try:
#    f1hand=open(f1name)
#    print(f1name)
#except:
#    print('Input file bed no exist')

df = pd.read_csv(f1name, sep ='\t', names=["CHROM","POS","DP","AF","AD","GQ"], dtype={"CHROM":str,"POS":str,"DP":str,"AF":str,"AD":str,"GQ":str})
#df = pd.read_csv(f1name, sep ='\t', header=0)
print(df.shape)
print(df.head)

maxAD = [ 0 if max(i.split(','))=='.' else int(max(i.split(','))) for i in df.AD]
sumDP = [ 0 if i=='.' else int(i) for i in df.DP ]
maxAF = [ 0.00 if sumDP[i]==0 else round(maxAD[i]/sumDP[i],2) for i in range(len(df.AD)) ]
avgDP = sum(sumDP)/len(sumDP)
relDP = [ int(100*i/avgDP) for i in sumDP]
#newGQ = [ 0 if i=='.' else int(i) for i in df.GQ ]
newGQ = [ 0 if max(i.split(','))=='.' else int(max(i.split(','))) for i in df.GQ]
newCHROM = [ str(i) for i in df.CHROM ]
newPOS = [ int(i) for i in df.POS ]
#print(df.head())
#print((df['end'] - df['start']).head())
#print(df)
#df_length = df['end'] + 1 - df['start']
#total_DP = np.dot(df_length, df['DP'])
#total_GQ = np.dot(df_length, df['GQ'])
f2 = {'CHROM':newCHROM,'POS':newPOS,'relDP':relDP, 'maxAF':maxAF, 'GQ':newGQ}
df2 = pd.DataFrame(f2)
df2 = df2[['CHROM','POS','relDP','maxAF','GQ']]
df2.to_csv(f1name_re+'.norm.tsv', sep='\t', index=False)

