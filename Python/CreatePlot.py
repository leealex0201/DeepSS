#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 11 15:21:46 2018

@author: Alex
"""

import numpy as np
import matplotlib.pyplot as plt
import os

def FiltData(Data, Cutoff):
    import scipy.signal
    N = 4 # order of butterworth filter
    SamplingRate = 30000 # 30 kHz sampling rate 
    Wn = 250/(SamplingRate*0.5)
    b, a = scipy.signal.butter(N, Wn, 'highpass')
    FilteredData = scipy.signal.filtfilt(b, a, Data)
    
    return FilteredData

def GetFilteredData(Data):
    ThreshThresh = 0.002 # threshold of threshold
    n, bins = np.histogram(Data, bins = np.linspace(np.amin(Data), np.amax(Data), num=500))
    
    n = np.float64(np.append(n, n[-1])) # seems like there is a size mismatch
    n = n/np.amax(n) # normalize it

    MaxInd = np.argmin(np.absolute(n-1))    
    MinInd = np.argmin(np.absolute(n-ThreshThresh))
    if MaxInd > MinInd:
        # peak is on the right side
        AnotherMinInd = MaxInd+(MaxInd-MinInd)
        if AnotherMinInd > n.size:
            AnotherMinInd = n.size
    elif MaxInd < MinInd:
        # peak is on the left side
        AnotherMinInd = MaxInd-(MinInd-MaxInd)
        if AnotherMinInd < 0:
            AnotherMinInd = 0
    UpBnd = np.max([bins[MinInd], bins[AnotherMinInd]])
    LowBnd = np.min([bins[MinInd], bins[AnotherMinInd]])
    
    return Data[(Data<UpBnd) & (Data>LowBnd)]
    
def GetThreshold(Data):
    ThresholdOffset = 10
    ThreshThresh = 0.02 # threshold of threshold
    n, bins = np.histogram(Data, bins = np.linspace(np.amin(Data), np.amax(Data), num=500))
    
    n = np.float64(np.append(n, n[-1])) # seems like there is a size mismatch
    n = n/np.amax(n) # normalize it
    
    MinInd = np.argmin(np.absolute(n-ThreshThresh))
    
    return -(np.abs(bins[MinInd])+ThresholdOffset), np.abs(bins[MinInd])+ThresholdOffset

def get_filtered_values(dist, seq):

    prev_val = seq[0]
    compare_to = prev_val + dist
    filtered = [prev_val]

    for elt in seq[1:]:
        if elt <= compare_to:           # <-- change to `<` to match desired results; 
                                        # this matches the results of your implementation 
            continue
        else:
            compare_to = elt + dist
            filtered.append(elt)
    return filtered
    
def GetWfAndTimestamp(ContData, Cutoff):
    # Cutoff: cut off frequency for filtering. Default = 250 Hz

    preThreshold = 10
    spikelen = 48 # searched spike waveform has length of 1.6 ms (48 points considering 30 kHz sampling rate of continuous data)
    
    # filter data
    ContData = FiltData(ContData, Cutoff)
#    ContData = GetFilteredData(ContData)
    
    # thresholding
    Threshold, ThresholdUp = GetThreshold(ContData) # compute threshold value first
    Timestamp = np.where(ContData<Threshold)[0]
    Timestamp = Timestamp - preThreshold
    Timestamp = np.delete(Timestamp, np.where(Timestamp<1)[0]) # remove elements that are less than 1
    Timestamp = np.delete(Timestamp, np.where(Timestamp+spikelen>ContData.size)[0])
    SelectedTimestamp = get_filtered_values(spikelen, Timestamp) # make sure all spikes are one spike apart (spikelen)
    SelectedTimestamp = np.asanyarray(SelectedTimestamp)
    
    # collect the waveforms. Make sure each waveforms are one spike away
    Wf = np.zeros((SelectedTimestamp.size, spikelen))
    for i in range(SelectedTimestamp.size):
        Wf[i,:] = ContData[SelectedTimestamp[i]:SelectedTimestamp[i]+48]
    
    return Wf, SelectedTimestamp
    
def ProcessData(FilePath, ind, Cutoff=250):
    import csv
    with open("ContDataExample.csv") as f:
        reader = csv.reader(f)
        # next(reader) # skip header
        FullData = [r for r in reader]
        
    # we will process only one row of data
    Data = FullData[ind]
    
    Data = np.array(Data).astype(np.float) # str -> float
    
    Wf, Timestamp = GetWfAndTimestamp(Data, Cutoff)
    
    return Wf, Timestamp

def plot_pca(x,y,T,S,N):    
    plt.scatter(x[T], y[T], alpha=0.2, c=(.5, .5, .5))
    plt.scatter(x[S], y[S], c=np.arange(S.size), cmap='gist_rainbow', s=100)
    plt.xticks([])
    plt.yticks([])
    plt.box(on=None)
    plt.savefig(N+'PCA.png')
    plt.show()
    
def plot_PT(PT,y,T,S,N):
    plt.scatter(PT[T], y[T], alpha=0.2, c=(.5, .5, .5))
    plt.scatter(PT[S], y[S], c=np.arange(S.size), cmap='gist_rainbow', s=100)
    plt.xticks([])
    plt.yticks([])
    plt.box(on=None)
    plt.savefig(N+'PT.png')
    plt.show()

def GetPT(Wf):
    PT = np.zeros((Wf.shape[0],))
    for i in range(Wf.shape[0]):
        PT[i] = np.max(Wf[i,:])-np.min(Wf[i,:])
    return PT

def plot_Wf(Wf,T,S,N):
    plt.plot(np.transpose(Wf[T,:]), alpha=0.1, c=(.5, .5, .5))
    color_idx = np.linspace(0, 1, S.size)
    for i in range(color_idx.size):
        plt.plot(Wf[S[i],:], color=plt.cm.gist_rainbow(color_idx[i]))
    plt.xticks([])
    plt.yticks([])
    plt.box(on=None)
    plt.savefig(N+'Wf.png')
    plt.show()

def CombinePCAAndPT(NameOfThis):
    # In this function, I try to combine two or three image files
    import PIL

    # PCA and PT first
    PCAImg = NameOfThis+'PCA.png'
    PTImg = NameOfThis+'PT.png'
    ListofImgs = [PCAImg, PTImg]
    Imgs = [PIL.Image.open(i) for i in ListofImgs]
    ImgsComb = np.hstack((np.asarray(i) for i in Imgs))
    ImgsComb = PIL.Image.fromarray(ImgsComb)
    SizeImgsComb = ImgsComb.size
    
    # Add Wf
    WfImgName = NameOfThis+'Wf.png'
    WfImg = PIL.Image.open(WfImgName)
    
    AllImgs = [ImgsComb,WfImg]
    AllCombinedImages = np.vstack( (np.asarray( i.resize(SizeImgsComb) ) for i in AllImgs ) )
    AllCombinedImages = PIL.Image.fromarray(AllCombinedImages)
    AllCombinedImages.save(NameOfThis+'Combined.png')

def GetImportantPCElements(F,S):
    '''
    In this function, I will attempt to get the index of timestamp that corresponds
    to the coordinate (FPC or SPC space) that has high histogram value in PC
    
    I have first PC (F) and the second PC (S). Now, what I want to do is sort
    of clustering, but not really. I would like to first, normalize it so that
    the maximum value of histogram is 1. Then, would like to get all indicies 
    of dots (or samples) that exceeds that threshold. I will have a few or seve
    ral clusters. I will color FROM those clusters. EVENLY DISTRIBUTED
    ''' 

    # get 2D histogram
    SIZE = 100    
    H, xedges, yedges = np.histogram2d(F, S, bins=(np.linspace(np.amin(F), np.amax(F), num=SIZE), np.linspace(np.amin(S), np.amax(S), num=SIZE)))
    TempH = np.append(H, np.zeros((SIZE-1,1)), axis=1)
    TempTempH = np.append(TempH, np.zeros((1,SIZE)), axis=0)
    H = TempTempH
    del TempH
    del TempTempH
    
    H = H/np.amax(H) # normalize it
    Threshold = 0.15
    
    IndX, IndY = np.where(H>Threshold) # indicies of edges
    FInds = xedges[IndX]
    SInds = yedges[IndY]

    return FInds, SInds

def GetMinDistInds(F,S,FPC,SPC):
    Dist = np.zeros((FPC.shape[0],))
    Dist = np.sqrt(np.square(FPC-F)+np.square(SPC-S))
    
    return np.argmin(Dist)

def PlotData(Wf, TotalInds, NameOfThis):
    
    # 1. PCA
    from sklearn.decomposition import PCA
    
    pca = PCA(n_components=2).fit_transform(Wf)
    FPC = pca[:,0]
    SPC = pca[:,1]
    
    FInds, SInds = GetImportantPCElements(FPC,SPC) 
    NumbOfSelectedWf = 25
    SelectedInds = np.random.choice(FInds.shape[0], NumbOfSelectedWf, replace=False)
    WfSelectedInds = np.zeros((NumbOfSelectedWf,), dtype=int)
    for i in range(NumbOfSelectedWf):
        WfSelectedInds[i] = GetMinDistInds(FInds[SelectedInds[i]],SInds[SelectedInds[i]],FPC,SPC)
    
    plot_pca(FPC,SPC,TotalInds,WfSelectedInds,NameOfThis)
    
    # 2. Peak-trough
    PT = GetPT(Wf)
    
    plot_PT(PT,Timestamp,TotalInds,WfSelectedInds,NameOfThis)
    
    # 3. Waveform plot
    plot_Wf(Wf,TotalInds,WfSelectedInds,NameOfThis)
    
    # combine all
    CombinePCAAndPT(NameOfThis)
    
# # # MAIN # # #
ContFilePath = '/Users/Alex/Google Drive/N_Hatsopoulos/SpikeSorter/ContDataExample.csv'
NameOfThis = 'Test'
NumbOfNeurons = 3
ProcessingIndex = 2
Cutoff = 250 # cut off frequency

Wf, Timestamp = ProcessData(ContFilePath, ProcessingIndex, Cutoff)

NumbTotalInds = 700 # number of total wf to plot (grey dots)
TotalInds = np.random.choice(Wf.shape[0], NumbTotalInds)

directory = os.getcwd() # current directory
directory = directory + '/' + str(NumbOfNeurons)
if not os.path.exists(directory):
    os.makedirs(directory)
Title = directory + '/' + NameOfThis

PlotData(Wf, TotalInds, Title)
