'''Removes part of the GW170814 data because there is some kind of
 glitch in Livingston before the event'''
import pandas as pd
df = pd.read_csv('GW170814_old.dat', sep='\s+')
N = len(df)
df[3*N//8 : 6*N//8+1].to_csv('GW170814.dat', sep='\t', index=False)