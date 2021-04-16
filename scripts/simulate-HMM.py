#!/usr/bin/env python
"""
Refer to https://notebook.community/jmschrei/pomegranate/tutorials/B_Model_Tutorial_3_Hidden_Markov_Models
"""
import numpy as np 
import pandas as pd

emission = pd.read_csv("params/CpG.emission.txt",sep="\t",index_col=0)
transition = pd.read_csv("params/CpG.transition.txt",sep="\t",index_col=0)
emission.index = emission.index.astype(str)

def simulate(transition,emission):
    idx = np.where(np.random.multinomial(1,transition.loc["s",:].values.reshape(-1))==1)[0][0]
    state = transition.columns[idx]
    states = ""
    signal = ""
    while(True):
        idx = np.where(np.random.multinomial(1,transition.loc[state,:].values.reshape(-1))==1)[0][0]
        state = transition.columns[idx]
        if state == "e":
            break
        states += state
        idx = np.where(np.random.multinomial(1,emission.loc[state,:].values.reshape(-1))==1)[0][0]
        emit = emission.columns[idx]
        signal += emit
    return states,signal

for i in range(100):
    states,signal = simulate(transition,emission)
    print(f">{i}")
    print(signal)
    print(states)


