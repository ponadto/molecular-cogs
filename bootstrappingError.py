import numpy
import cPickle as pkl
import gzip
from random import randrange, choice
import glob
import os
import time
import sys


angles = sys.argv[1:]
'''
fileNames = glob.glob("ge*gz")
angles = []
for file in fileNames:
	angles.append(float(file.split("=")[1].split("_")[0]))
print "len(angles): ",len(angles)
for angle in angles: print "python bootstrappingError.py ",angle
'''
interactions = ["ele","vdw","bonds","angle","tors"]

N_OF_BINS = 1000
N_OF_RUNS = 100

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
        WOOOOOOOOOOOLNE, ale w C++ bylo rownie wolne
    '''
    sequence_np = numpy.array(sequence)
    statistics = [ statistic(numpy.random.choice(sequence_np,size=len(sequence),replace=True)) for run in xrange(nOfRuns)]
    return ( numpy.mean(statistics), numpy.var(statistics) )

    
def getSequences(fileName,interaction,angle,zksis):
    index = 0
    streams = [ [ open("tmp_%s_%.2f_%d_%d.dat" % (interaction,angle,j,i),"w") for i in xrange(j+1,11)] for j in xrange(10)]
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
    	        if row>10:
    	           row = 0
    	           index += 1
    except IndexError:
	print "index = ",index
    for i in xrange(10):
        for j in xrange(i+1,11):
            streams[i][j-i-1].close()




for angle in angles:
    angle = float(angle)
    zksis = []
    dAdKsis = []
    entropic = []
    fileName = "generalOutput_angle=%.2f_timestep=0.002_twistPeriod=1000_nOfSamples=1000000_binWidth=0.05_idleSteps=50.gz" % angle
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
    
    wyputDoEnergiiSwobodnej = open("dAdKsiAndZksi=%.2f.dat"%angle,"w")
    (mean,var) = blockBootstrap(numpy.mean,N_OF_BINS,dAdKsis,N_OF_RUNS)
    wyputDoEnergiiSwobodnej.write( "dAdKsi = %f +/- %f\n" % (mean,numpy.sqrt(var)) )
    (mean,var) = blockBootstrap(numpy.mean,N_OF_BINS,zksis,N_OF_RUNS)
    wyputDoEnergiiSwobodnej.write( "<Zksi**(-.5)> = %f +/- %f\n" % (mean,numpy.sqrt(var)) )
    (mean,var) = blockBootstrap(numpy.mean,N_OF_BINS,entropic,N_OF_RUNS)
    wyputDoEnergiiSwobodnej.write( "<entropic> = %f +/- %f\n" % (mean,numpy.sqrt(var)) )
    wyputDoEnergiiSwobodnej.close()
    
    tot = 0
    
    for interaction in interactions: 
        start = time.clock()   
        fileName = "matrices_angle=%.2f_timestep=0.002_twistPeriod=1000_nOfSamples=1000000_binWidth=0.05_idleSteps=50_interaction=%s.gz" % (angle,interaction)
        print fileName
        getSequences(fileName,interaction,angle,zksis)    
        print "sequences loaded! (it took %.2f sec)" % (time.clock()-start)
        
        blockBootstrapVarianceWyput = open("blockBootstrapVariance_%s_%.2f.dat" % (interaction,angle),"w")
        classicBootstrapVarianceWyput = open("classicBootstrapVariance_%s_%.2f.dat" % (interaction,angle),"w")
        dumbVarianceWyput = open("dumbVariance_%s_%.2f.dat" % (interaction,angle),"w")
        blockBootstrapMeanWyput = open("blockBootstrapMean_%s_%.2f.dat" % (interaction,angle),"w")
        classicBootstrapMeanWyput = open("classicBootstrapMean_%s_%.2f.dat" % (interaction,angle),"w")
        dumbBootstrapMeanWyput = open("dumbMean_%s_%.2f.dat" % (interaction,angle),"w")
        
        for row in xrange(10):
            for column in xrange(row+1,11):
                start = time.clock() 
                with open("tmp_%s_%.2f_%d_%d.dat" % (interaction,angle,row,column)) as wput:
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
                    print "...classicBootstrap..."
                    (tmp_mean,tmp_var) = classicBootstrap(numpy.mean,seq_np,N_OF_RUNS)
                    classicBootstrapVarianceWyput.write("%.15f " % tmp_var)
                    classicBootstrapMeanWyput.write("%.15f " % tmp_mean)
                    dumbVarianceWyput.write("%.15f " % (numpy.var(seq_np)/len(seq))) # TODO
                    dumbBootstrapMeanWyput.write("%.15f " % numpy.mean(seq_np)) # TODO
                    print "... whole took %.2f sec\n" % (time.clock()-start )
            blockBootstrapVarianceWyput.write("\n")
            classicBootstrapVarianceWyput.write("\n")
            dumbVarianceWyput.write("\n")
            blockBootstrapMeanWyput.write("\n")
            classicBootstrapMeanWyput.write("\n")
            dumbBootstrapMeanWyput.write("\n")
            
            
        for fileName in glob.glob("tmp_%s_%.2f_*dat" % (interaction,angle)):
            os.remove(fileName)
                    
        blockBootstrapVarianceWyput.close()
        classicBootstrapVarianceWyput.close()
        dumbVarianceWyput.close()
        blockBootstrapMeanWyput.close()
        classicBootstrapMeanWyput.close()
        dumbBootstrapMeanWyput.close()
        
    print "tot = ",tot
    (mean,var) = blockBootstrap(numpy.mean,N_OF_BINS,dAdKsis,N_OF_RUNS)
    print "dAdKsi = %f +/- %f\n" % (mean,numpy.sqrt(var))
    
