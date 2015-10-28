import ipdb
import math
import pdb
import re
# import Image
import os

import numpy as np
import vigra

import basicOperations


## A window centered at GT coordinate
winSize = 231
winSizeR = winSize/2
winSizeR2 = winSize/2+1
cut_orig = True
cut_result = True

## directories
path_result = "/Volumes/Xiwei_MacExt/output/test2"
path_imin = "/Volumes/Xiwei_MacExt/cecog_data/train_40"
path_csv = "/Volumes/Xiwei_MacExt/databases/mitos_atypia"
path_out_mitosis  = "/Volumes/Xiwei_MacExt/output/test2_crops"
# path_out_not_mitosis  = "/Volumes/Xiwei_MacExt/output/test1_crop"
# path_result = "/home/seawave/work/output/test2_div2"
# path_imin = "/home/seawave/work/database/train_40"
# path_csv = "/home/seawave/work/database/mitos_atypia"
# path_out_mitosis  = "/home/seawave/work/output/test2_crops"



## select mitosis count GTs
GT_list_all = basicOperations.getFiles(path_csv,'csv')
orig_image_all = basicOperations.getFiles(path_imin,'tiff')
results = basicOperations.getFiles(path_result,'png')

mitosis_set = []
not_mitosis_set = []
L1 = len(path_csv)
for filetemp1 in GT_list_all:
    filetemp2 = filetemp1[L1+1:]
    m = re.search("mitosis", filetemp2)
    if m!=None:
        m2 = re.search("not", filetemp2)
        if m2!=None:
            not_mitosis_set.append(filetemp1)
        else:
            mitosis_set.append(filetemp1)
    
## cut images
if (1):  ## for mitosis image 
    
    pattern1 = ".+train_40/(?P<slide>[A-Za-z0-9]+)/(?P<filename>[A-Za-z0-9_]+).tiff"
    toto = re.compile(pattern1)

    count = 0
    for image in orig_image_all:
        print count, image
        count += 1

        tata = toto.match(image)
        slide = tata.groupdict()['slide']
        filename = tata.groupdict()['filename']
        csv_file = os.path.join(path_csv, slide)
        csv_file = os.path.join(csv_file, "mitosis")
        csv_file = os.path.join(csv_file, filename+"_mitosis.csv")
        result_file = os.path.join(path_result, filename + ".png")

        if not os.path.isfile(csv_file):
            print "ERROR1 !!", image
        if not os.path.isfile(result_file):
            print "ERROR1 !!", image
   
        if cut_orig: ### cut original image
            imOrig = vigra.readImage(image) 
            imSize = imOrig.shape

            imOrigLarge = vigra.RGBImage((imOrig.shape[0]+winSizeR2*2,imOrig.shape[1]+winSizeR2*2))
            imOrigLarge[:,:,:] = 255
            imOrigLarge[winSizeR2:winSizeR2+imOrig.shape[0], winSizeR2:winSizeR2+imOrig.shape[1], :] = imOrig[:,:,:]

        if cut_result: ### cut result image
            imResult = vigra.readImage(result_file) 
            imSize = imResult.shape

            imResultLarge = vigra.Image((imResult.shape[0]+winSizeR2*2,imResult.shape[1]+winSizeR2*2))
            imResultLarge[:,:,:] = 255
            imResultLarge[winSizeR2:winSizeR2+imResult.shape[0], winSizeR2:winSizeR2+imResult.shape[1], :] = imResult[:,:,:]
   
        ## get mitosis GT coordinates
        textTemp = open(csv_file,'r')
        lines =  textTemp.readlines()

        if len(lines)==0:
            continue

        i = 0
        for line in lines:
            words = line.split(',')
            x = int(words[0]) + winSizeR2
            y = int(words[1]) + winSizeR2
   
            if cut_orig:
                imCrop = imOrigLarge[x-winSizeR:x+winSizeR+1, y-winSizeR:y+winSizeR+1]
                outName = filename+"_"+str(i)+".png"
                outfile = os.path.join(path_out_mitosis, outName)
                vigra.impex.writeImage(imCrop,outfile)

            if cut_result:
                imCrop = imResultLarge[x-winSizeR:x+winSizeR+1, y-winSizeR:y+winSizeR+1]
                outName = filename+"_"+str(i)+"_test.png"
                outfile = os.path.join(path_out_mitosis, outName)
                vigra.impex.writeImage(imCrop,outfile)


            i+=1

ipdb.set_trace()

if (1):  ## for mitosis image 
    
    pattern1 = ".+train_40/(?P<slide>[A-Za-z0-9]+)/(?P<filename>[A-Za-z0-9_]+).tiff"
    toto = re.compile(pattern1)

    count = 0
    for image in orig_image_all:
        print count, image
        count += 1

        tata = toto.match(image)
        slide = tata.groupdict()['slide']
        filename = tata.groupdict()['filename']
        csv_file = os.path.join(path_csv, slide)
        csv_file = os.path.join(csv_file, "mitosis")
        csv_file = os.path.join(csv_file, filename+"_mitosis.csv")
        result_file = os.path.join(path_result, filename + ".png")

        if not os.path.isfile(csv_file):
            print "ERROR1 !!", image
        if not os.path.isfile(result_file):
            print "ERROR1 !!", image
   
        imOrig = vigra.readImage(image) 
        imSize = imOrig.shape

        imOrigLarge = vigra.RGBImage((imOrig.shape[0]+winSizeR2*2,imOrig.shape[1]+winSizeR2*2))
        imOrigLarge[:,:,:] = 255
        imOrigLarge[winSizeR2:winSizeR2+imOrig.shape[0], winSizeR2:winSizeR2+imOrig.shape[1], :] = imOrig[:,:,:]
    
        ## get mitosis GT coordinates
        textTemp = open(csv_file,'r')
        lines =  textTemp.readlines()

        if len(lines)==0:
            continue

        i = 0
        for line in lines:
            words = line.split(',')
            x = int(words[0]) + winSizeR2
            y = int(words[1]) + winSizeR2
   
            imCrop = imOrigLarge[x-winSizeR:x+winSizeR+1, y-winSizeR:y+winSizeR+1]
            outName = filename+"_"+str(i)+".png"
            outfile = os.path.join(path_out_mitosis, outName)
            vigra.impex.writeImage(imCrop,outfile)

            i+=1
 



if (0):  ## for not mitosis image 
    count = 0
    for filetemp in not_mitosis_set:
        print count, filetemp
        count += 1
    
        ## find original image
        words = filetemp.split("/")
        m = re.search("not_mitosis", words[-1])
        fileName = words[-1][:m.start()-1]
    
        for image in orig_image_all:
            m1 = re.search(fileName, image)
            m2 = re.search("x40", image)
            if m1 != None and m2 != None:
                break
    
        imOrig = vigra.readImage(image) 
        imSize = imOrig.shape

        imOrigLarge = vigra.RGBImage((imOrig.shape[0]+winSizeR2*2,imOrig.shape[1]+winSizeR2*2))
        imOrigLarge[:,:,:] = 255
        imOrigLarge[winSizeR2:winSizeR2+imOrig.shape[0], winSizeR2:winSizeR2+imOrig.shape[1], :] = imOrig[:,:,:]
    
        ## get mitosis GT coordinates
        textTemp = open(filetemp,'r')
        lines =  textTemp.readlines()

        if len(lines)==0:
            continue

        i = 0
        for line in lines:
            words = line.split(',')
            x = int(words[0]) + winSizeR2
            y = int(words[1]) + winSizeR2
   
            imCrop = imOrigLarge[x-winSizeR:x+winSizeR+1, y-winSizeR:y+winSizeR+1]
            outName = fileName+"_"+str(i)+".png"
            vigra.impex.writeImage(imCrop,outName)
            basicOperations.writeFiles(outName, path, image, path_out_not_mitosis, "", "png")

            i+=1
 



    ## Find original image

#    osType = getOsType() ## 1 for windows, 2 for linux
#    files = getFiles(path,form)
#    m = re.search('jpg', origname, re.I)
#    end = m.start()-1
#    slash = []
#    if osType == 1:
#        for m in re.finditer('\\\\',origname):  ## for windows
#            slash.append(m.start())
#            start = slash[len(slash)-2]
#            pattern = origname[start:slash[len(slash)-1]]+"\\\\"+origname[slash[len(slash)-1]+1:end]
#    elif osType == 2:
#        for m in re.finditer('/',origname):  ## for linux
#            slash.append(m.start())
#            start = slash[len(slash)-2]
#            pattern = origname[start:end]
#
#    for rep in files:
#        m = re.search(pattern,rep)
#        if m!= None:
#            for m2 in re.finditer('/',rep):
#                start = m2.start()+1
#            temp = rep[start:len(rep)]
#            m2 = re.search(suffixe,temp)
#            if m2!=None:
#                print rep
#                return rep
#    if m==None:
#        print "File not found"
#        return m
#        pdb.set_trace()




pdb.set_trace()

i = 0
for image in im_list:
    tt1 = time.time()
    print image

    # order = "./TeleOphta_EX "+image ## linux
    # order = order.replace('(','\(')
    # order = order.replace(')','\)')
    order = "TeleOphta.exe "+image ## windows
    os.system(order)

    imROI = fileRead("imROI.png")
    imVessel = fileRead("imVessel.png")
    imOD = fileRead("imOD.png")

#     file_OD = basicOperations.findPath(image,path_OD,'png','')
    file_GT = basicOperations.findPath(image,path_GT,'png','')
    
    ## Read image
    im = Image.open(image)
    im.save("temp.png","PNG")
    imin = fileRead("temp.png")
    im3 = extractChannels(imin)
    
    ## Detect ROI (visible region in retinal image)
#     imROI = detectROI.detectROI(im3,10)
    # imROI = getSame(im3[1])
    # ImThreshold(im3[0],1,255,255,0,imROI)
#     fileWrite(imROI,'ROI.png')
    
    ## Estimated sizes of the image
    imInfo = basicOperations.sizeEstimate(imROI)
    print "Disque:",imInfo[0], "Optic Disc:", imInfo[1], "Vessel:", imInfo[2], "MA",imInfo[3]
    
    ## Segment vessels
#     imVessel = detectVessel.detectVessel(im3[1],imROI,imInfo)
    # imVesselC = fileRead(file_vessel)  ## Color image
    # imVessel3 = extractChannels(imVesselC)
    # imVessel = imVessel3[1]
#     fileWrite(imVessel,'Vessel.png')
    
    ## Detect Optic Disc
    # imOD = detectOD.detectOpticDisc(im3,imROI,imVessel)
#     imOD = fileRead(file_OD)
#     fileWrite(imOD,"OpticDisc.png")
    
    ## Detect Exudate
    ## To choose the characts you want to calculate 
    ## For details, see the beginning of ccAnalyse.py
#    imGT = fileRead(file_GT)
    charact = [	1,1,1,1,1,1,1,1,1,1,\
    		    1,1,1,1,1,1,1,1,1,1,\
    		    1,1,1,1,1,1,1,1,1	]
    imExudate = detectExudat.detectExudat(im3,imVessel,imROI,imOD,imInfo,image,0,charact)
    
    pdb.set_trace()
    if imExudate[1]==0:
        f = open("EXcandis2.txt","w")
        f.write(' ')
        f.close()

    basicOperations.writeFiles(imExudate[0], path, image, output_path1, "","png")
    basicOperations.writeFiles(imOD, path, image, output_path2, "","png")
    basicOperations.writeFiles(imVessel, path, image, output_path3, "","png")
    
    basicOperations.writeText("EXcandis2.txt", path, image, output_path1,"","txt")
    
    # basicOperations.writeFiles_b(imExudate[2], path, image, output_path3, "","png")

   
    # if imExudate[1]!=0:
        # f = open("ExdList.txt","w")
        # for x in imExudate[1]:
            # for y in x:
                # f.write(y+' ')
            # f.write('\n')
        # f.close()
    # else:
        # f = open("ExdList.txt","w")
        # f.write(' ')
        # f.close()


    




pdb.set_trace()



