EXAMPLE = '''
Example run:

    python bootstrappingError.py -179.00 ii-nh3

where:
    -179.00 --- is the dihedral angle for which the bootstrapping will be carried out
    ii-nh3  --- is the molecule name for which the bootstrapping will be ran (auxiliary data will be read from "ii-nh3.conf")
'''


import numpy
import cPickle as pkl
import gzip
from random import randrange, choice
import glob
import os
import time
import sys
import smallFunctions


INTERACTIONS = ["ele","vdw","bonds","angle","tors"]
N_OF_BINS = 1000
N_OF_RUNS = 100
CALC_CLASSIC_BOOTSTRAP = False # classic bootstrap is VERY time consuming
PATH_TO_LARGE_OUTPUT = "/home/ponadto/Desktop/largeOutput2/"
PATH_TO_BOOTSTRAP_OUTPUT = "/home/ponadto/molecular-cogs/bootstrapOutput/"
   


'''
fileNames = glob.glob("ge*gz")
angles = []
for file in fileNames:
	angles.append(float(file.split("=")[1].split("_")[0]))
print "len(angles): ",len(angles)
for angle in angles: print "python bootstrappingError.py ",angle
'''


def blockBootstrap(statistic,nOfBins,sequence,nOfRuns):
    statistics = []
    lengthOfBins = [int(len(sequence)/nOfBins) for i in xrange(nOfBins-1)]
    lengthOfBins.append(len(sequence)-(nOfBins-1)*int(len(sequence)/nOfBins))
    
    for run in xrange(nOfRuns):
        bins = [ randrange(nOfBins) for i in xrange(nOfBins) ]
        tmp = []
        for b in bins:
            tmp += sequence[b*nOfBins:(b+1)*nOfBins]
        statistics.append(statistic(tmp))
    return (numpy.mean(statistics),numpy.var(statistics))
    

    
def classicBootstrap(statistic,sequence,nOfRuns):
    ''' 
        SLOOOOOOOOOOOW, but not much slower than when written in C++
    '''
    sequence_np = numpy.array(sequence)
    statistics = [ statistic(numpy.random.choice(sequence_np,size=len(sequence),replace=True)) for run in xrange(nOfRuns)]
    return ( numpy.mean(statistics), numpy.var(statistics) )

    
def getSequences(fileName,interaction,angle,zksis,conf):
    index = 0
    nOfAtoms = conf["nOfAtoms"]
    prefix = conf["prefix"]
    streams = [ [ open(PATH_TO_BOOTSTRAP_OUTPUT+"tmp_%s_%.2f_%s_%d_%d.dat" % (interaction,angle,prefix,j,i),"w") for i in xrange(j+1,nOfAtoms)] for j in xrange(nOfAtoms-1)]
    try:
	with gzip.open(fileName) as f:
    	    row = 0
    	    for line in f:
    	        words = line.split(" ")[:-1]
    	        for column in xrange(len(words)):
    	            if column>row:
    	                streams[row][column-row-1].write(str(float(words[column])*zksis[index])+" ")
    	                #streams[row][column-row-1].write(str(float(words[column]))+" ")
    	        row += 1
    	        if row>nOfAtoms-1:
    	           row = 0
    	           index += 1
    except IndexError:
	print "index = ",index
    for i in xrange(nOfAtoms-1):
        for j in xrange(i+1,nOfAtoms):
            streams[i][j-i-1].close()











if __name__=="__main__":
    
    try:
        angle = sys.argv[1]
        prefix = sys.argv[2]
    except IndexError:
        print "Wrong arguments"
        print EXAMPLE
        sys.exit(1)

    conf = smallFunctions.parseConf(prefix)

    for key in sorted(conf.keys()):
        print "%s : %s" % (key,str(conf[key]))

    prefix = conf["prefix"]
    timeStep = conf["timeStep"]
    twistPeriod = conf["twistPeriod"]
    nOfSamples = conf["nOfSamples"]
    binWidth = conf["binWidth"]
    idleSteps = conf["idleSteps"]
    nOfAtoms = conf["nOfAtoms"]

    angle = float(angle)
    zksis = []
    dAdKsis = []
    entropic = []
    fileName = PATH_TO_LARGE_OUTPUT+"generalOutput_prefix=%s_angle=%.2f_timestep=%s_twistPeriod=%d_nOfSamples=%d_binWidth=%s_idleSteps=%d.gz" % (prefix,angle,str(timeStep),twistPeriod,nOfSamples,str(binWidth),idleSteps)
    print fileName
    with gzip.open(fileName,"r") as wput:
        for line in wput:
            words = line.split()
            if "dihedral"==words[0]:
                continue
            zksi = float(words[5])
            zksis.append(zksi**(-0.5))
            dAdKsi = float(words[3])
            dAdKsis.append(dAdKsi*zksi**(-0.5))
            #dAdKsis.append(dAdKsi)
            gradUgradKsi = float(words[4])
            entropic.append( (dAdKsi-gradUgradKsi)*zksi**(-0.5) )
    
    wyputDoEnergiiSwobodnej = open(PATH_TO_BOOTSTRAP_OUTPUT+"dAdKsiAndZksi=%.2f.dat"%angle,"w")
    (mean,var) = blockBootstrap(numpy.mean,N_OF_BINS,dAdKsis,N_OF_RUNS)
    wyputDoEnergiiSwobodnej.write( "dAdKsi = %f +/- %f\n" % (mean,numpy.sqrt(var)) )
    (mean,var) = blockBootstrap(numpy.mean,N_OF_BINS,zksis,N_OF_RUNS)
    wyputDoEnergiiSwobodnej.write( "<Zksi**(-.5)> = %f +/- %f\n" % (mean,numpy.sqrt(var)) )
    (mean,var) = blockBootstrap(numpy.mean,N_OF_BINS,entropic,N_OF_RUNS)
    wyputDoEnergiiSwobodnej.write( "<entropic> = %f +/- %f\n" % (mean,numpy.sqrt(var)) )
    wyputDoEnergiiSwobodnej.close()
    
    tot = 0
    
    for interaction in INTERACTIONS: 
        start = time.clock()   
        fileName = PATH_TO_LARGE_OUTPUT+"matrices_prefix=%s_angle=%.2f_timestep=%s_twistPeriod=%d_nOfSamples=%d_binWidth=%s_idleSteps=%d_interaction=%s.gz" % (prefix,angle,str(timeStep),twistPeriod,nOfSamples,str(binWidth),idleSteps,interaction)
        print fileName
        getSequences(fileName,interaction,angle,zksis,conf)    
        print "sequences loaded! (it took %.2f sec)" % (time.clock()-start)
            
        blockBootstrapVarianceWyput = open(PATH_TO_BOOTSTRAP_OUTPUT+"blockBootstrapVariance_%s_%.2f.dat" % (interaction,angle),"w")
        if CALC_CLASSIC_BOOTSTRAP:
            classicBootstrapVarianceWyput = open(PATH_TO_BOOTSTRAP_OUTPUT+"classicBootstrapVariance_%s_%.2f.dat" % (interaction,angle),"w")
        dumbVarianceWyput = open(PATH_TO_BOOTSTRAP_OUTPUT+"dumbVariance_%s_%.2f.dat" % (interaction,angle),"w")
        blockBootstrapMeanWyput = open(PATH_TO_BOOTSTRAP_OUTPUT+"blockBootstrapMean_%s_%.2f.dat" % (interaction,angle),"w")
        if CALC_CLASSIC_BOOTSTRAP:
            classicBootstrapMeanWyput = open(PATH_TO_BOOTSTRAP_OUTPUT+"classicBootstrapMean_%s_%.2f.dat" % (interaction,angle),"w")
        dumbBootstrapMeanWyput = open(PATH_TO_BOOTSTRAP_OUTPUT+"dumbMean_%s_%.2f.dat" % (interaction,angle),"w")
     
        for row in xrange(nOfAtoms-1):
            for column in xrange(row+1,nOfAtoms):
                start = time.clock() 
                with open(PATH_TO_BOOTSTRAP_OUTPUT+"tmp_%s_%.2f_%s_%d_%d.dat" % (interaction,angle,prefix,row,column)) as wput:
                    print "opening "+wput.name
                    print "calculating..."
                    start = time.clock() 
                    seq = map(float,wput.readline().split(" ")[:-1])
                    seq_np = numpy.array(seq)
                    print "...blockBootstrap..."
                    (tmp_mean,tmp_var) = blockBootstrap(numpy.mean,N_OF_BINS,seq,N_OF_RUNS)
                    if not numpy.isnan(tmp_mean): tot += tmp_mean
                    else: print "NO NIE!!"
                    print "tot = ",tot
                    blockBootstrapVarianceWyput.write("%.15f " % tmp_var )
                    blockBootstrapMeanWyput.write("%.15f " % tmp_mean )
                    if CALC_CLASSIC_BOOTSTRAP:
                        print "...classicBootstrap..."
                        (tmp_mean,tmp_var) = classicBootstrap(numpy.mean,seq_np,N_OF_RUNS)
                        classicBootstrapVarianceWyput.write("%.15f " % tmp_var)
                        classicBootstrapMeanWyput.write("%.15f " % tmp_mean)
                    dumbVarianceWyput.write("%.15f " % (numpy.var(seq_np)/len(seq))) # TODO
                    dumbBootstrapMeanWyput.write("%.15f " % numpy.mean(seq_np)) # TODO
                    print "... whole took %.2f sec\n" % (time.clock()-start )
            blockBootstrapVarianceWyput.write("\n")
            if CALC_CLASSIC_BOOTSTRAP:
                classicBootstrapVarianceWyput.write("\n")
            dumbVarianceWyput.write("\n")
            blockBootstrapMeanWyput.write("\n")
            if CALC_CLASSIC_BOOTSTRAP:
                classicBootstrapMeanWyput.write("\n")
            dumbBootstrapMeanWyput.write("\n")
        
           
        for fileName in glob.glob("tmp_%s_%.2f_%s_*dat" % (interaction,angle,prefix)):
            os.remove(fileName)
                
        blockBootstrapVarianceWyput.close()
        if CALC_CLASSIC_BOOTSTRAP:
            classicBootstrapVarianceWyput.close()
        dumbVarianceWyput.close()
        blockBootstrapMeanWyput.close()
        if CALC_CLASSIC_BOOTSTRAP:
            classicBootstrapMeanWyput.close()
        dumbBootstrapMeanWyput.close()
    
    print "tot = ",tot
    (mean,var) = blockBootstrap(numpy.mean,N_OF_BINS,dAdKsis,N_OF_RUNS)
    print "dAdKsi = %f +/- %f\n" % (mean,numpy.sqrt(var))
 

