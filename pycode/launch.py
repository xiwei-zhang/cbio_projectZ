"""
This is the version used to run over entire database.
"""
import ipdb

import pdb
import os
import sys
import re
import time
import math
import numpy as np
import random
import cPickle
# from openopt import QP
import thread

# sys.path.append("../../")
import basicOperations

tbegin = time.time()

def shuffleSeq(N):
    temp = []
    for i in range(N):
        temp.append(i)
    random.shuffle(temp)
    return temp

# def myGaussianKernal(x,y):
def getRBFParameter(trSet, teSet):
    allSet = np.vstack((trSet, teSet))
    # dist1 = metric.pairwise_distances(trSet)
    dist2 = metric.pairwise_distances(allSet)
    gamma = 1 / (np.median(dist2))
    return gamma


def loadAllDataset(pickleFile, train_folders, N_featPF):
   
    pklFiles = basicOperations.getFiles(input_path, "cpkl")

    # pattern1 = ".+(?P<slide>slide[0-9]).+cpkl"
    pattern1 = ".+pickles/(?P<slide>[A-Za-z0-9]+)/.+cpkl"
    toto1 = re.compile(pattern1)
    # pattern2 = ".+/(?P<filename>[A-Za-z0-9]+_[A-Za-z0-9]+_[0-9]+).+cpkl"
    pattern2 = ".+/(?P<filename>[A-Za-z0-9]+_[A-Za-z0-9]+).+cpkl"
    toto2 = re.compile(pattern2)
    
    feat = {}
    count = 1
    Nobjects = [0,0,0] # number of mitos, notmitos and others
    for folder in train_folders:
        mito = []
        notmito = []
        notannot = []
        others = []
        for pkl in pklFiles:
            m = re.search(folder, pkl)
            if m == None:
                continue
            pkltemp = open(pkl, 'rb')
            feats = cPickle.load(pkltemp)
            for i in range(len(feats)):
                if feats[i][-1]==1:
                    mito.append(feats[i][:-3])
                elif feats[i][-1]==2:
                    notmito.append(feats[i][:-3])
                if feats[i][-1]==3:
                    notannot.append(feats[i][:-3])
            # centers.append(feats[i][-3:-1])
            # slides.append(slide)
            pkltemp.close()
        # idx_shuffle = shuffleSeq(len(others))
        # for i in range(N_featPF):
        #     notannot.append(others[idx_shuffle[i]])
        # for x in mito:
        #     mitos.append(x)
        # for x in notmito:
        #     notmitos.append(x)


        ######################################
        ### check is there invalide (nan or inf) element
        mito = np.array(mito)
        notmito = np.array(notmito)
        notannot = np.array(notannot)
        toto = [mito, notmito, notannot]
        for i in range(3):
            if np.isnan(toto[i]).any():
                nanElem = np.isnan(toto[i])
                nanElemRow = np.unique(np.where(nanElem == True)[0])
                infElemRow = np.unique(np.where(toto[i] >= np.finfo(np.float32).max)[0])
                nanElemRow = np.unique(np.hstack((nanElemRow, infElemRow)))
                arBool = np.ones(len(toto[i]))
                for x in nanElemRow:
                    arBool[x] = 0
                arBool = arBool.astype(np.bool)
                toto[i] = toto[i][arBool]
                ## centers = centers[arBool]
                ## slides = slides[arBool]
                ## filenames = filenames[arBool]
        
        mito = toto[0]
        notmito = toto[1]
        notannot = toto[2]

        for i in range(3):
            for j in range(len(toto[i])):
                Nobjects[i] += 1
                feat[count] = {}
                feat[count]['folder'] = folder
                feat[count]['features'] = toto[i][j]
                if i==0:
                    feat[count]['annotation'] = 'mito'
                elif i==1:
                    feat[count]['annotation'] = 'notmito'
                else:
                    feat[count]['annotation'] = 'others'
                count += 1


    # f = open(pickleFile,'r')
    # feat = pickle.load(f)
    # f.close()

    ## get feature_name ##
    # featureName = feat[feat.keys()[0]]['features'].keys()
    print Nobjects
    featureName = 0
    return feat, featureName



def loadDataset(pickleFile, train_folders, N_featPF):
  
    pklFiles = basicOperations.getFiles(input_path, "cpkl")

    # pattern1 = ".+(?P<slide>slide[0-9]).+cpkl"
    pattern1 = ".+pickles/(?P<slide>[A-Za-z0-9]+)/.+cpkl"
    toto1 = re.compile(pattern1)
    # pattern2 = ".+/(?P<filename>[A-Za-z0-9]+_[A-Za-z0-9]+_[0-9]+).+cpkl"
    pattern2 = ".+/(?P<filename>[A-Za-z0-9]+_[A-Za-z0-9]+).+cpkl"
    toto2 = re.compile(pattern2)
    
    feat = {}
    count = 1
    Nobjects = [0,0,0] # number of mitos, notmitos and others
    for folder in train_folders:
        mito = []
        notmito = []
        notannot = []
        others = []
        for pkl in pklFiles:
            m = re.search(folder, pkl)
            if m == None:
                continue
            pkltemp = open(pkl, 'rb')
            feats = cPickle.load(pkltemp)
            for i in range(len(feats)):
                if feats[i][-1]==1:
                    mito.append(feats[i][:-3])
                elif feats[i][-1]==2:
                    notmito.append(feats[i][:-3])
                if feats[i][-1]==3:
                    others.append(feats[i][:-3])
            # centers.append(feats[i][-3:-1])
            # slides.append(slide)
            pkltemp.close()
        idx_shuffle = shuffleSeq(len(others))
        for i in range(N_featPF):
            notannot.append(others[idx_shuffle[i]])
        # for x in mito:
        #     mitos.append(x)
        # for x in notmito:
        #     notmitos.append(x)


        ######################################
        ### check is there invalide (nan or inf) element
        mito = np.array(mito)
        notmito = np.array(notmito)
        notannot = np.array(notannot)
        toto = [mito, notmito, notannot]
        for i in range(3):
            if np.isnan(toto[i]).any():
                nanElem = np.isnan(toto[i])
                nanElemRow = np.unique(np.where(nanElem == True)[0])
                infElemRow = np.unique(np.where(toto[i] >= np.finfo(np.float32).max)[0])
                nanElemRow = np.unique(np.hstack((nanElemRow, infElemRow)))
                arBool = np.ones(len(toto[i]))
                for x in nanElemRow:
                    arBool[x] = 0
                arBool = arBool.astype(np.bool)
                toto[i] = toto[i][arBool]
                ## centers = centers[arBool]
                ## slides = slides[arBool]
                ## filenames = filenames[arBool]
        
        mito = toto[0]
        notmito = toto[1]
        notannot = toto[2]

        for i in range(3):
            for j in range(len(toto[i])):
                Nobjects[i] += 1
                feat[count] = {}
                feat[count]['folder'] = folder
                feat[count]['features'] = toto[i][j]
                if i==0:
                    feat[count]['annotation'] = 'mito'
                elif i==1:
                    feat[count]['annotation'] = 'notmito'
                else:
                    feat[count]['annotation'] = 'others'
                count += 1


    # f = open(pickleFile,'r')
    # feat = pickle.load(f)
    # f.close()

    ## get feature_name ##
    # featureName = feat[feat.keys()[0]]['features'].keys()
    print Nobjects
    featureName = 0
    return feat, featureName

    
def getDataset(feat, featureName, folderLeftOut, N_featPF): 
    folderNames = ["A03","A04","A05","A07","A10","A11","A12","A14","A15",\
            "A17","A18", "H03","H04","H05","H07","H10","H11","H12","H14","H15",\
            "H17","H18"]

    train_set = []
    train_class = []
    test_set = []
    test_class = []

    #######
    ### clustering to extract notannotated samples
    for f in folderNames:
        print " ",f
        # if f!="A10":
        #     continue
        data_set = []
        temp = []
        if f!=folderLeftOut:
            for k in feat.keys():
                if feat[k]['folder'] != f:
                    continue
                if feat[k]['annotation'] == 'others':
                    data_set.append(feat[k]['features'])
            data_set_, test_set_ = featNormalization(data_set, data_set[0][:])
            clusters = clustering(data_set_, 200)

            idx_select = []
            for i in range(200):
                idx_select.append([])

            if (0): ## rank according to the histogram of clusters
                hist = np.histogram(clusters, bins=200)[0]
                ranklist = np.argsort(hist)[::-1]
                ranklist_map = np.zeros(200).astype(np.int)
                for i in range(len(ranklist)):
                    ranklist_map[ranklist[i]] = i

                for i in range(len(clusters)): 
                    idx_select[ranklist_map[clusters[i]]].append(i)

            else: ## don't rank
                for i in range(len(clusters)): 
                    idx_select[clusters[i]].append(i)

            ## shuffle
            for i in range(len(idx_select)):
                random.shuffle(idx_select[i])

            c = 0
            p = 0
            count = 0
            while count < N_featPF and p < len(data_set):
                if p<len(idx_select[c]):
                    train_set.append(data_set[idx_select[c][p]])
                    train_class.append(0)
                    temp.append(0) 
                    count += 1
                c+=1
                if c>=200:
                    c = 0
                    p += 1
 
            # print len(temp), len(data_set), count

    for k in feat.keys():
        if feat[k]['folder'] == folderLeftOut: ## test set
            continue
        if feat[k]['annotation'] == 'mito':
            train_set.append(feat[k]['features'])
            train_class.append(1)
        elif feat[k]['annotation'] == 'notmito':
            train_set.append(feat[k]['features'])
            train_class.append(0) 

    return train_set, train_class, test_set, test_class

def featNormalization(learn_set, test_set, method = 3):
    ml = len(learn_set)
    mt = len(test_set)
    if method == 1:
        Xl = preprocessing.scale(learn_set)
        Xt = preprocessing.scale(test_set)
    elif method == 2:
        all_set = learn_set + test_set
        XX = preprocessing.scale(all_set)
        Xl = XX[:ml]
        Xt = XX[ml:]
    elif method == 3:
        mm = np.mean(learn_set, axis=0)
        std = np.std(learn_set, axis=0)
        std[ np.where(mm==0) ] = 1 ## prevent dividing by zero
        Xl = (learn_set - mm) / std
        Xt = (test_set - mm) / std
    return Xl, Xt


def kernelMeanMatching(Xl, Xt, gamma=0.04, B=1.2, e=0.001):
## list of potential parameters:
## gamma: 0.03, B: 1000, e = 0.1, 1.026, 0.0001, 0.9, 0.99
## gamma: 0.02, B: 1000, e = 0.1, 1.853, 0.0000, 0.9, 1.001
    ##for gm in (10.0**np.arange(-2,2)):
    ml = len(Xl)
    mt = len(Xt)
    GramK1 = pairwise.rbf_kernel(Xl, gamma = gamma)
    GramK_ = pairwise.rbf_kernel(Xl, Xt, gamma = gamma)
    # GramK1 = pairwise.linear_kernel(Xl)
    # GramK_ = pairwise.linear_kernel(Xl, Xt)
    GramK2 = 2  * np.sum(GramK_, axis=1) / (mt * ml) # axis=1 means sum of each row
    GramK1 = GramK1 / (ml * ml)
    
    ## QP optimization
    ## IMP!! Assign values to B and e
    ## default form of QP solver:
    ## p = QP(np.diag([1, 2, 3]), [15, 8, 80], A = np.matrix('1 2 3; 8 15 80'), b = [150, 800], Aeq = [0, 1, -1], beq = 25.5, ub = [15,np.inf,np.inf])
    A = np.vstack((np.ones(ml), -1*np.ones(ml)))
    b = [ml*(e+1), ml*(e-1)]
    p = QP(GramK1, GramK2, A = A, b = b, lb = np.zeros(ml), ub = np.ones(ml)*B)
    # p = QP(GramK1, GramK2, lb = np.zeros(ml), ub = np.ones(ml)*B, Aeq = np.ones(ml), beq = ml)
    r = p._solve('cvxopt_qp',iprint=1)
    f_opt, x_opt = r.ff, r.xf
    print "gamma:",gamma, "B:", B, "e:",e
    print np.max(x_opt), np.min(x_opt), np.mean(x_opt), np.median(x_opt)

    return x_opt




def svm(train_set, train_class, test_set, C=1, gamma=0.001, sample_weight = False, beta=-1, class_weight='auto'):

    # clf = SVC(C = C, gamma = gamma, cache_size = 2000, class_weight = 'auto',\
    #         probability=True)
    clf = LinearSVC(C = C, penalty='l1', dual=False, class_weight = 'auto')

    if sample_weight:
        clf.fit(train_set, train_class, sample_weight=beta )
    else:
        clf.fit(train_set, train_class)

    # result = clf.predict(test_set)
    dist = clf.decision_function(test_set)
    # dist = clf.predict_proba(test_set)


    return dist


def randomForest(train_set, train_class, test_set, N_est=500, n_jobs=2, beta=-1):

    rf = RandomForestClassifier(n_estimators = N_est,  n_jobs = n_jobs)
    # rf.fit(train_set, train_class, sample_weight=beta )
    rf.fit(train_set, train_class )
    dist = rf.predict_proba(test_set)
    return dist


def clustering(data_set, n_cluster, batch_size = 100 ):
    k_means = MiniBatchKMeans(n_clusters = n_cluster, batch_size = batch_size)

    t0 = time.time()
    k_means.fit(data_set)
    t1 = time.time()
    print "Mini batch kmeans finished. Time used:", t1 - t0

    toto = k_means.labels_
    return toto

    #pickleFile = open("clusters.cpkl", "wb")
    #cPickle.dump(toto, pickleFile, -1)
    #pickleFile.close()


def clustering_DBSC(data_set, eps=0.5, min_samples=10, leaf_size=30):
    dbsc = DBSCAN(eps = eps, min_samples=min_samples, leaf_size = leaf_size)

    t0 = time.time()
    db = dbsc.fit(data_set)
    t1 = time.time()
    print "DBSCAN finished. Time used:", t1 - t0

    toto = k_means.labels_
    return toto



def resultAnalyse(result, test_class):
    auc = metric.roc_auc_score(test_class,result)
    fpr, tpr, th = metric.roc_curve(test_class, result)

    f1 = []
    for t in th:
        aa = (result > t).astype(int)[:,0]
        f1.append(metric.f1_score(test_class, aa))

    Fs = max(f1)

    ### plot ###
    xx = np.linspace(0.0, 1.0, len(fpr))
    fig, ax = plt.subplots()
    plt.plot(fpr, tpr, label="ROC curve (area = %0.3f)" % auc)
    plt.plot(xx, f1, 'r-', label="F-score (max = %0.3f)" % Fs)
    plt.plot([0,1],[0,1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC of Kenel mean matching')
    plt.legend(loc="lower right")
    plt.savefig("svm_tf_0.037.png", dpi=200)
    plt.show()

def featSelection(learn_set, test_set, learn_class, k=10, method = 0):
    if 0: ## reserved
        ml = len(learn_set)
        mt = len(test_set)
        all_set = learn_set + test_set
        XX = SelectKBest(chi2, k=k).scale(all_set)
        Xl = XX[:ml]
        Xt = XX[ml:]
    if 0:
        ml = len(learn_set)
        xl = (learn_set - np.min(learn_set)) / np.max(learn_set)
        XX  = SelectKBest(chi2, k=k)
        kb = XX.fit(xl, learn_class)
        feat = kb.get_support()
        return feat

    if 1: ## RFE
        estimator = SVC(kernel = 'linear', cache_size = 2000, class_weight = 'auto')
        # estimator = SVC(C = 1, gamma = 0.001, cache_size = 2000, class_weight = 'auto')

        selector = RFE(estimator,  200, 50)
        selector.fit(learn_set, learn_class)

        f = open("temp2.txt","a")
        f.write("\n")
        selector.get_support(True).tofile(f, sep=" ")
        f.write("\n")
        f.close()

        return selector.support_

        

def getTestset(pkl_file, testFolder):
    test_set = []
    test_class = []
    f = open(pkl_file, "rb")
    feat = cPickle.load(f)
    f.close()

    for k in feat.keys():
        if feat[k]['folder'] == testFolder:
            test_set.append(feat[k]['features'])
            if feat[k]['annotation'] == 'mito':
                test_class.append(1)
            else:
                test_class.append(0)

    return test_set, test_class


def writeIntoFile(result, test_class, filename):
    f = open(filename, 'w')
    for x,y in zip(result, test_class):
        f.write(str(x)+" "+str(y)+"\n")
    f.close()


def multi_proc(images, n_job):
    time.sleep(0.3*n_job)
    print "thread:", n_job, "starting"

    pattern1 = ".+train_40/(?P<slide>[A-Za-z0-9]+)/(?P<filename>[A-Za-z0-9_]+).tiff"
    toto = re.compile(pattern1)

    for image in images:
        tata = toto.match(image)
        slide = tata.groupdict()["slide"]
        filename = tata.groupdict()["filename"]
        
        imout = filename + ".png"
        output_file = os.path.join(output_path, imout)

        csv = filename + "_mitosis.csv"
        csvfile = os.path.join(csv_path, slide)
        csvfile = os.path.join(csvfile, "mitosis")
        csvfile = os.path.join(csvfile, csv)

        csv2 = filename + "_not_mitosis.csv"
        csvfile2 = os.path.join(csv_path, slide)
        csvfile2 = os.path.join(csvfile2, "mitosis")
        csvfile2 = os.path.join(csvfile2, csv2)

        featout = filename + ".txt"
        featout_file = os.path.join(output_path, featout)

        order = "./projectZ " + image + " " + output_file + " " \
                + csvfile + " " + csvfile2 + " " + featout_file

        os.system(order)
    




if __name__ == "__main__":
    ######
    # parameter
    n_job = 4
    
    input_path = "/home/seawave/work/database/train_40"
    csv_path = "/home/seawave/work/database/mitos_atypia"
    output_path = "/home/seawave/work/output/test3"
    
    
    train_folders = ["A03","A04","A05","A07","A10","A11","A12","A14","A15",\
            "A17","A18", "H03","H04","H05","H07","H10","H11","H12","H14","H15",\
            "H17","H18"]
    # train_folders = ["A03","A04","A05"]
   
    images = basicOperations.getFiles(input_path, "tiff")

    if 0: ### check not analysed images
        images_not_analysed = []
        pattern1 = ".+train_40/(?P<slide>[A-Za-z0-9]+)/(?P<filename>[A-Za-z0-9_]+).tiff"
        toto = re.compile(pattern1)
        for image in images:
            tata = toto.match(image)
            filename = tata.groupdict()['filename']
            imout = filename+".png"
            outfile = os.path.join(output_path, imout)
            if not os.path.isfile(outfile):
                images_not_analysed.append(image)

        if len(images_not_analysed) == 0:
            print "finished"
            ipdb.set_trace()
        else:
            images = images_not_analysed

    imageGroup = []
    interV = len(images)/4
    for i in range(n_job):
        subGroup = []
        for j in range(i*interV, (i+1)*interV):
            subGroup.append(images[j])
        imageGroup.append(subGroup)
        


    # multi_proc(imageGroup[0], 0)
    # ipdb.set_trace()

    try:
        thread.start_new_thread( multi_proc, (imageGroup[0], 0) )
        thread.start_new_thread( multi_proc, (imageGroup[1], 1) )
        thread.start_new_thread( multi_proc, (imageGroup[2], 2) )
        thread.start_new_thread( multi_proc, (imageGroup[3], 3) )


    except:
        print "Error: unable to start thread"



    ipdb.set_trace()

    for f in train_folders:
    
        print f
        t_begin = time.time()
    
        ## Get data sets ##
        train_set, train_class, test_set, test_class = getDataset(feat, featureNames, f, N_featPF)
        test_set, test_class = getTestset("dataset.cpkl", f)
    
        ## Features normalizations ##
        train_set, test_set = featNormalization(train_set, test_set)
    
    
        ## feature selection
        if (0):
            featSupport = featSelection(train_set, test_set, train_class, 100)
            continue
            train_set = train_set[:,featSupport]
            test_set = test_set[:,featSupport]
    
        ## Estimatae the RBF kernel feature ##
        # GM = getRBFParameter(train_set, test_set)
        # print " estimated RBF kernel parameter:", GM
    
        ## Kerner mean matching ## 
        # beta = kernelMeanMatching(train_set, test_set, GM, 1.2, 0.001)
    
        ## SVM ##
        # dist = svm(train_set, train_class, test_set, 1, 0.001, True, beta)
        # dist = svm(train_set, train_class, test_set, 1, 0.001)
        dist = randomForest(train_set, train_class, test_set, 500, 4)
    
        ## Write into file ##
        # writeIntoFile(dist[:,0], test_class, "svm_result_"+f+".txt")
        writeIntoFile(dist, test_class, "svm_result_"+f+".txt")
    
        t_end = time.time()
        print "Time used:", t_end - t_begin, "s"
    
    pdb.set_trace()



##########
## grid search
# C_range = 10.0 ** np.arange(-1,2)
# gamma_range = 10.0 ** np.arange(-3,0)
