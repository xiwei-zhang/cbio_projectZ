
import pdb
# import Image
import os
import re
import sys
import math

import numpy as np
from numpy import linalg

###############################
## Function sizeEstimate
## Get size information of the image. An estimation of the size of OpticDisc, vessel, and lesion from the size of ROI. (on pixels)
## Algo by Xiwei ZHANG
## code by Xiwei ZHANG
##
## input: imROI
## output: [Diam_ROI, Diam_OD, Diam_Vessel, Diam_MA] 
##
## variable:	L: python list, image of ROI
##		size: [width, height] image size
##		x0: First non-zero pixel x coordonner from left side of image
## 		f: flag, to stop loop
##		i,j: coordonner
##		Diam_ROI, Diam_OD, Diam_Vessel, Diam_MA: Diameter of ROI, estimated diameters of OD,Vessel and MA
def sizeEstimate(imROI):

	L = ImageToArray(imROI) ## transfer MM image to Python list
	size = imROI.getSize()

	x0 = size[0]  ## find first non-zero pixel x coordonner from left side of image
	f = 0
	i = 0
	while i<size[0]:
		j = 0
		while j<size[1]:
			if L[i+size[0]*j]!=0:
				x0 = i%size[0]
				f = 1
				break
			j = j+1
		if f == 1:
			break
		i = i+1

	if x0>size[0]/4: ## in case of totally black image
		return [2,2,2,2]

	Diam_ROI = size[0]-2*x0 ## estimate parameters
	Diam_OD = int(round(float(Diam_ROI)/5))
	Diam_Vessel = int(round(float(Diam_ROI)/76.6)) 
	Diam_MA = int(round(float(Diam_ROI)/910*9)) ## old ver: int(round(Diam_ROI/87)) 
	# print "Disque:",Diam_ROI, "Optic Disc:", Diam_OD, "Vessel:", Diam_Vessel, "MA", Diam_MA
	return [Diam_ROI, Diam_OD, Diam_Vessel, Diam_MA]
	
## end of function sizeEstimate
###############################



###############################
## Function: resize
## resize the images using Python Image Library. For the detection of OD, large size is not necessary
##
## Input:	imin: image to be resized
##		size: the width of new image
## Output:	imout: resized image
#def resize(imin,size,size2=0):
#    osType = getOsType() ## 1 for windows, 2 for linux
#    fileWrite(imin,'temp1.png')
#    im = Image.open('temp1.png')
#
#    if size2!=0:
#        size_new = (size,size2)
#    else:
#        size_orig = im.size
#        r0 = float(size_orig[0])/size
#        r1 = float(size_orig[0])/size_orig[1]
#        size_new = (size,int(size/r1))
#
#    im_resize = im.resize(size_new)
#    im_resize.save('temp.png','PNG')
#    imout = fileRead('temp.png')
#
#    if osType == 1:
#        os.system("del temp.png") ## windows
#    elif osType == 2:
#        os.system("rm -f ./temp.png")  ## linux
#
#    return imout
### End of function resize
################################



###############################
## Function MeanFilterFFT
## Using FFT to do large scale mean filter
##
## Input:	imin: input image
##		l: kernal size
## Output:	imtempD3: filtered image (IMP note, it's a double image)
def MeanFilterFFT(imin,l):
	imPad1 = FFT.ImPadding(imin,0)
	imtempC1 = ImCreateSame(imPad1,"complex<F_DOUBLE>")
	imtempC2 = ImCreateSame(imPad1,"complex<F_DOUBLE>")
	imtempC3 = ImCreateSame(imPad1,"complex<F_DOUBLE>")
	imtempD4 = ImCreateSame(imPad1,"F_DOUBLE")
	imtempD3 = ImCreateSame(imin,"F_DOUBLE")
	listKernel = [1/float(l)**2]*l**2
	imKernel = ArrayToImage(listKernel,l,l,1)
	imKernelPad = FFT.ImPadding(imKernel,0,imPad1.getWxSize(),imPad1.getWySize(),1)

	FFT.ImFreeFFT(imKernelPad,imtempC1)	# kernel fft
	FFT.ImFreeFFT(imPad1,imtempC2)		# image orig fft	

	arithMultImage(imtempC1,imtempC2,imtempC2)
	FFT.ImFreeIFFT(imtempC2,imtempD4)
	FFT.ImFreeFFTShift(imtempD4)
	FFT.ImUnPadding(imtempD4,imtempD3)
	return imtempD3
## End of MeanFilterFFT
###############################


###############################
## Function ColorToS
## Extract saturation channel from color image
##
## Input:	image: address of color image
## Output:	imout: saturation image
def ColorToS(image):
	order = "./ColorToS "+image
	os.system(order)
	imout = fileRead("saturation.png")
	return imout	
## End of RGB2H
###############################


	
### Creat a structure element ###	
def CreatSE12(size,n):
	""" Creat a structure element in 12 directions. (used in detection of vessel)
	Input: size(int)
	Output: list of 12 vectors. (need to convert to SE form)
	"""
	
	l = size
	print "taille opening: ",l
	pi = 3.1415
	t = pi/n
	angle = []

	for i in range(1,n+1):
		angle.append((i-1)*t)

	x_l = range(-l/2,l/2+1)
	x = []
	for theta in angle:
		temp = []
		for xx in x_l:
			temp.append(xx*math.cos(theta))
		x.append(temp)

	y = []
	for i in range(0,n):
		temp = []
		for xx in x[i]:
			temp.append(xx*math.tan(angle[i]))
		y.append(temp)

	ne = []
	for i in range(0,n):
		temp = []
		for index in range(len(x[i])):
			temp.append(int(round(x[i][index])))
			temp.append(int(round(y[i][index])))
		ne.append(temp)
	return ne


### Creat a structure element ###	
def CreatSESquare(radius):
	""" Creat a square structure element. (used in filters, eg: median filter)
	Input: radius(int)
	Output: SE
	"""
	i = 0-radius
	pointList = []
	while i <= radius:
		j = 0-radius
		while j <= radius:
			pointList.append(i)
			pointList.append(j)
			j = j+1
		i = i+1
	SE = NeighborListFactory(ConnexityType.Square, 2, pointList)
	return SE


### Creat a structure element ###
def CreatSECross(size):
	pointList=[]
	i = 0-size
	while i <= size:
		j = 0-size
		while j<=size:
			if i==0 or j==0:
				pointList.append(i)
				pointList.append(j)
			j = j+1
		i = i+1
	SE = NeighborListFactory(ConnexityType.Square, 2, pointList)
	return SE


### Creat a structure element ###
def CreatSEH(size):
	pointList = []
	i = 0-size
	while i<=size:
		pointList.append(i)
		pointList.append(0)
		i = i+1
	SE = NeighborListFactory(ConnexityType.Square, 2, pointList)
	return SE

### Creat a structure element ###
def CreatSEV(size):
	pointList = []
	i = 0-size
	while i<=size:
		pointList.append(0)
		pointList.append(i)
		i = i+1
	SE = NeighborListFactory(ConnexityType.Square, 2, pointList)
	return SE


### Creat structure element ###
def CreatSEX(size):
	pointList=[]
	i = 0-size
	while i <= size:
		j = 0-size
		while j<=size:
			if abs(i) == abs(j):
				pointList.append(i)
				pointList.append(j)
			j = j+1
		i = i+1
	SE = NeighborListFactory(ConnexityType.Square, 2, pointList)
	return SE

    

### Get files ###
def getFiles(path, extension, suffix = None):
	""" Get all the files with 'extension' in a directory.
	Input: path (directory of images), extension (jpg, png..., ignor case)
	Output: a list of files
	"""	
	listFile = []
	for root, dirs, files in os.walk(path):
		dirs.sort()
		files.sort()
		for item in files:
			if os.path.splitext(item)[1][1:].lower() == extension.lower():	
				if suffix!=None:
					m = re.search(suffix,item)
					if m!=None:
						listFile.append(os.path.join(root,item))
				else:
					listFile.append(os.path.join(root,item))
	return listFile
	
    
### Write files ###
def writeFiles(fileName,input_path,current_path, output_path, name, extension):

    order = 'mv -f ./'+fileName+' '+ output_path + "/" + name  ## linux
    os.system(order)



#    osType = getOsType()
#    if osType == 1:
#        order = 'move '+nameOut+' '+rep   ## windows
#    elif osType == 2:
#        order = 'mv -f ./'+nameOut+' '+ rep   ## linux
#    l = len(input_path)
#    dossier = current_path[l+1:]
#    i = 0 
#    rep = output_path
#    p1 = 0
#    while i<len(dossier):
#        if(dossier[i]=='/' or dossier[i]=='\\'):
#            if osType == 1:
#                rep = rep+"\\"+dossier[p1:i]
#            elif osType == 2:
#                rep = rep + "/" + dossier[p1:i]  ## linux
#            if  ~os.path.isdir(rep):
#                os.system("mkdir "+rep)
#            p1 = i+1
#        if dossier[i]=='.':
#            nameOut = dossier[p1:i]+name+'.'+extension
#            break
#        i = i+1
#    pdb.set_trace()
#    fileWrite(imout,nameOut)
	
def writeText(fileName,input_path,current_path, output_path, name, extension):
    osType = getOsType()
    l = len(input_path)
    dossier = current_path[l+1:]

    i = 0
    rep = output_path
    p1 = 0
    while i<len(dossier):
        if(dossier[i]=='/' or dossier[i]=='\\'):
            if osType == 1:
                rep = rep+"\\"+dossier[p1:i]
            elif osType == 2:
                rep = rep + "/" + dossier[p1:i]  ## linux
            if  ~os.path.isdir(rep):
                os.system("mkdir "+rep)
            p1 = i+1
        if dossier[i]=='.':
            nameOut = dossier[p1:i]+name+'.'+extension
            break
        i = i+1

    if osType == 1:
        order = 'move '+fileName+ ' '+rep+'\\'+nameOut
    elif osType == 2:
        order = 'mv -f ./'+fileName+' '+ rep+'/'+nameOut
    os.system(order)

### Combine the results of the segmentation ###
def combine(imin,imColor):
	imout = getSame(imColor)
	im3 = Initialisation(imColor)
	for im in im3:
		arithSupImage(imin,im,im)
	imout = combineChannels(im3[0],im3[1],im3[2])
	return imout
	
### find path ###
def findPath(origname,path,form,suffixe):
    osType = getOsType() ## 1 for windows, 2 for linux
    files = getFiles(path,form)
    m = re.search('jpg', origname, re.I)
    end = m.start()-1
    slash = []
    if osType == 1:
        for m in re.finditer('\\\\',origname):  ## for windows
            slash.append(m.start())
            start = slash[len(slash)-2]
            pattern = origname[start:slash[len(slash)-1]]+"\\\\"+origname[slash[len(slash)-1]+1:end]
    elif osType == 2:
        for m in re.finditer('/',origname):  ## for linux
            slash.append(m.start())
            start = slash[len(slash)-2]
            pattern = origname[start:end]

    for rep in files:
        m = re.search(pattern,rep)
        if m!= None:
            for m2 in re.finditer('/',rep):
                start = m2.start()+1
            temp = rep[start:len(rep)]
            m2 = re.search(suffixe,temp)
            if m2!=None:
                print rep
                return rep
    if m==None:
        print "File not found"
        return m
        pdb.set_trace()

### draw cross in MA ###
def drawCross(image,nl,SE):
	# SECross = CreatSECross(5)
	imtemp1 = getSame(image)
	imtemp2 = getSame(image)
	imtemp3 = getSame(image)
	imtemp4 = getSame(image)
	imtemp16_1 = ImCreateSame(image,'UINT16')
	imtemp16_2 = ImCreateSame(image,'UINT16')
	size = image.getSize()
	
	ImDilate(image,nl,imtemp2)
	Distance(imtemp2,nl,imtemp1)
	ImMaxima(imtemp1,nl,imtemp2)
	ImLabel(image,nl,imtemp16_1)
	ImCompare(imtemp2,">",0,imtemp16_1,0,imtemp16_2)
	hist = measHistogram(imtemp16_2)
	keys = hist.keys()
	del keys[0]
	p = []
	L = ImageToArray(imtemp16_2)
	for k in keys:
		i = 0
		while i< len(L):
			if L[i]==k:
				p.append(i)
				break

			i = i+1
	ImSetConstant(imtemp4,0)
	for co in p:
		imtemp4.setPixel(co,255)
	ImDilate(imtemp4,SE,imtemp1)
	return imtemp1

##### Get OS type #####
def getOsType():
    if (re.search('win',sys.platform)!=None):
        osType = 1
    elif (re.search('linux',sys.platform)!=None):
        osType = 2
    return osType


###### get file name #####
def getFileName(input_path, suffix):
    input_path = input_path.replace(' ', '_')
    m = re.search(suffix, input_path)
    if m == None:
        print "Wrong suffix"
    end = m.start()
    slash = []
    for m in re.finditer('/', input_path):  ## for linux
        slash.append(m.start())
        start = slash[len(slash)-1]+1
    pattern = input_path[start:end-1]

    return pattern


###### color deconvolution #####
def colorDeconv(imin,stainMatrix=0):
    M_h_e_dab_meas = np.array([[0.650, 0.072, 0.268],\
                               [0.704, 0.990, 0.570],\
                               [0.286, 0.105, 0.776]])
    M_h_e_meas = np.array([[0.644211, 0.092789],\
                           [0.716556, 0.954111],\
                           [0.266844, 0.283111]])
    
    M_h_e_meas_2 = np.array([[0.49015734, 0.04615336],\
                         [0.76897085, 0.8420684],\
                         [0.41040173, 0.5373925]])
 

    if stainMatrix == 0:
        M = M_h_e_dab_meas 
    elif stainMatrix == 1:
        M = M_h_e_meas 
    elif stainMatrix == 2:
        M = M_h_e_meas_2

    M_inv =  np.dot(linalg.inv(np.dot(M.T, M)), M.T)
    imout = np.dot(imin, M_inv.T)

    for i in range(imout.shape[2]):
        vmax = np.max(imout[:,:,i])
        vmin = np.min(imout[:,:,i])
        if vmax == vmin:
            imout[:,:,i] = 0
        else:
            imout[:,:,i] = (imout[:,:,i] - vmin) / (vmax - vmin) * 255

    return imout

           
 
