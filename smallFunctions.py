def parseConf(prefix):
    conf = {}
    conf["prefix"] = prefix
    with open(prefix+".conf") as inp:
        for line in inp:
            if line[0]=="#":
                continue
            words = line.split()
            if words[0] in conf:
                print "two specs for one thing in configuration file: %s" % fileName
                sys.exit(1)
            else:
                conf[words[0]] = words[1]

    # to float
    conf["timeStep"] = float(conf["timeStep"])
    conf["binWidth"] = float(conf["binWidth"])
    conf["nu"] = float(conf["nu"])

    # timeStep stuff
    conf["timeStep2"] = conf["timeStep"]**2
    conf["invTimeStep"] = 1./conf["timeStep"]
    conf["nuTimesTimeStep"] = conf["nu"]*conf["timeStep"]

    # to tuple
    conf["whichAtomsForDih"] = tuple([int(i) for i in conf["whichAtomsForDih"].split(",")])
    if "whichAtomsForAux" in conf:
        conf["whichAtomsForAux"] = tuple([int(i) for i in conf["whichAtomsForAux"].split(",")])

    # to int
    conf["nOfAtoms"] = int(conf["nOfAtoms"])
    conf["nOfSamples"] = int(conf["nOfSamples"])
    conf["twistPeriod"] = int(conf["twistPeriod"])
    conf["nOfStepsInBetween"] = int(conf["nOfStepsInBetween"])
    conf["idleSteps"] = int(conf["idleSteps"])

    return conf


def matrixForm(productsMatrices):
    for key in productsMatrices.keys():
        print "%s: " % key
        for line in productsMatrices[key]:
            print " ".join(["%f\t" % number for number in line])
        print

def drange(start, stop, step):
     r = start
     while r < stop:
     	yield r
     	r += step





if __name__=="__main__":
    print "TESTING"
