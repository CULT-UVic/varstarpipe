#!/usr/bin/python

import numpy as np
import pyfits
import os
from pyraf import iraf
import matplotlib.pyplot as plt
from scipy.ndimage.fourier import fourier_shift as imshift

#pipeline for doing variable star astronomy at UVic observatory

def imdisp(data):

	plt.imshow(data, cmap='gray', vmin=0.9*np.median(data),vmax=1.1*np.median(data))
	plt.colorbar()
	plt.show()


def fitsopen(filename):

	hdu = pyfits.open(filename)
	data = hdu[0].data
	head = hdu[0].header
	hdu.close()

	return data,head

def masterdark(filename):

    f = open(filename,'r')
    l = f.readlines()
    for x,i in enumerate(l): l[x]=i.replace('\n',"")

    darkcube = []
    for i in l:
        d,h=fitsopen(i)
        darkcube.append(d)

    mdark = np.array(np.median(darkcube, axis=0))
    return mdark

def masterflat(filename, mdark):

    f = open(filename,'r')
    l = f.readlines()
    for x,i in enumerate(l): l[x]=i.replace('\n',"")

    flatcube = []
    for i in l:
        d,h=fitsopen(i)
        flatcube.append(list(np.array(d)-np.array(mdark)))

    mflat = np.mean(flatcube, axis=0)
    normflat=np.median(mflat)
    mflat=mflat/normflat
    return mflat

def calibscience(filename, mflat, mdark):

    f = open(filename,'r')
    l = f.readlines()
    for x,i in enumerate(l): l[x]=i.replace('\n',"")

    sciencecube = []
    for i in l:
        d,h=fitsopen(i)
        sciencecalib=(list((np.array(d)-np.array(mdark))/np.array(mflat)))
        pyfits.writeto(i.replace('.fits','_calib.fits'),sciencecalib,h)

def padwithzeros(vector, pad_width, iaxis, kwargs):
    vector[:pad_width[0]] = 0
    vector[-pad_width[1]:] = 0
    return vector

def imgxcorr(im1,im2,posx,posy,boxsize,size):

    arr = range((-1)*size,size+1)

    xcorrmtx = np.zeros((2*size+1,2*size+1))
    shifts = np.zeros((2*size+1,2*size+1,2))

    im1 = im1[posx-boxsize:posx+boxsize,posy-boxsize:posy+boxsize]
    im2 = im2[posx-boxsize:posx+boxsize,posy-boxsize:posy+boxsize]

    im1 = np.lib.pad(im1,(size/2),mode='median')
    im2 = np.lib.pad(im2,(size/2),mode='median')

    #imdisp(im1)
    #imdisp(im2)

    im2 = np.fft.fft2(im2)
    sz = np.shape(im2)
    for i,x in enumerate(arr):
        for j,y in enumerate(arr):
	    im2_tst = imshift(im2,[x,y])
	    im2_tst = np.abs(np.real(np.fft.ifft2(im2_tst))) 
	    xcorrmtx[i-1,j-1] = np.sum(im2_tst*im1)
	    shifts[i-1,j-1] = [x,y]

    corr = np.where(xcorrmtx == np.max(xcorrmtx))
    print shifts[corr]
    #plt.imshow(xcorrmtx)
    #plt.show()
    return (shifts[corr])[0]
		
		

def imagestack(filename,posx,posy,boxsize,size):

    f = open(filename,'r')
    l = f.readlines()
    for x,i in enumerate(l): l[x]=i.replace('\n',"")

    img,h = fitsopen(l[0])
    
    for i in l[1:]:
	img2,h = fitsopen(i) 
	shft = imgxcorr(img,img2,posx,posy,boxsize,size)
	print shft
	img2 = np.fft.fft2(img2)
	img2 = imshift(img2,shft)
	img2 = np.abs(np.real(np.fft.ifft2(img2)))
	img = img+img2

    pyfits.writeto(i.replace('.fits','_stacked.fits'),img,h)


def vstarphot(filename):
    f = open(filename,'r')
    l = f.readlines()
    for x,i in enumerate(l): l[x]=i.replace('\n',"")

    iraf.noao()
    iraf.digiphot()
    iraf.apphot()
    iraf.images()

    os.system('rm *.coo*')

    for i in l:

        x=iraf.imstat(i, field="image, mean,stddev, midpt, min, max", Stdout=1)
        string=x[1]
        arrstring = ' '.join(string.split())
        arrstring=arrstring.split(" ")
        print arrstring
        midpt=arrstring[3]
        min=arrstring[4]
        max=arrstring[5]

        datapars=iraf.datapars.getParList()
        iraf.datapars.setParam('sigma',midpt)
        iraf.datapars.setParam('datamin',min)
        iraf.datapars.setParam('datamax',max)
        iraf.datapars.setParam('fwhmpsf','3.0')
        iraf.datapars.saveParList(filename='datapars.par')

        centerpars=iraf.centerpars.getParList()
        iraf.centerpars.setParam('cbox','40.0')
        iraf.centerpars.saveParList(filename='centerpars.par')

        findpars=iraf.findpars.getParList()
        iraf.findpars.setParam('threshold',2.0)
        iraf.findpars.saveParList(filename='findpars.par')

        daofind=iraf.daofind.getParList()
        iraf.daofind.setParam('image', i)
        iraf.daofind.setParam('verify','no')
        iraf.daofind.saveParList(filename='daofind.par')
        iraf.daofind(i)
