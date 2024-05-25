# Open each contrast csv file from Contrasts directory
# Copy >0 fold change and <0 fold change into different csv files

import os
import pandas as pd

dir = '/Users/u1501359/Documents/Lab work:protocols/mTSCs_ATAC_2022/Contrasts'
os.mkdir(dir+'/updownDARs')

os.chdir(dir)

def UpDown(contrasts):
    df = pd.read_csv(contrasts)
    f = 0
    upDARs = df[df['Fold'] >= f]
    downDARs = df[df['Fold'] < f]
    upDARs.to_csv(dir+'/updownDARs/upDARs_'+contrasts, sep=',', header = True, index=False)
    downDARs.to_csv(dir+'/updownDARs/downDARs_'+contrasts, sep=',', header = True, index=False)

# 'upDARs' = more accessible regions
# 'downDARs' = less accessible regions

files = [x for x in os.listdir(dir) if '.csv' in x]

for file in files:
    UpDown(file)
