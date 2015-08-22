import openbabel
import numpy
import math
import random
import time
import operator
import pylab
from itertools import izip
import gzip
import os
import glob

# TODO:
# popraw: 
# mksi_gradU_invMasses_gradKsis[index] = tmp_mksi_gradU_invMasses_gradKsi / ZksiToMinusHalf  (jest obecnie tmp_dAdKsi)


LARGE_OUTPUT = True
PREFIX = "ii-nh3_forGAFF"
WHICH_ATOMS = (1,5,6,7) # used to define the dihedral angle
RELEVANT_COORDINATES = []
for atomNumber in WHICH_ATOMS:
    for i in xrange(3):
        RELEVANT_COORDINATES.append(3*(atomNumber-1)+i)
# FOR TESTING PURPOSES:
#RELEVANT_COORDINATES = [i for i in xrange(33)]
MOL2_FILE_NAME = PREFIX+".mol2"
FF_NAME = "GAFF" # GAFF MMFF94 UFF
BOLTZMANN_CONSTANT = 0.001987 # [kcal/K/mol]
BIN_WIDTH = 0.05
TEMPERATURE = 300.
TIMESTEP = 0.002
TIMESTEP2 = TIMESTEP**2
INV_TIMESTEP = 1./TIMESTEP
IDLE_STEPS = 50
MEMORY_BUFFER_SIZE = 20000
PATH_TO_LARGE_OUTPUT = "/home/ponadto/Desktop/largeOutput/"
NUxTIMESTEP = 0.01*TIMESTEP # for AndersenIntegrator
    
    
    
    
class DataAndCalculations:
    """
        A class for holding TEMPORARY data produced during simulation, and for processing that data
    """

    def __init__(self,coords,masses):
        self.nOfAtomsX3 = len(coords)
        self.nOfAtoms = self.nOfAtomsX3/3
        self.setMasses(masses)
        self.reset(coords)
        
    def reset(self,coords):
        self.coords = [ coords[:] for i in xrange(3) ] # for times t-dt, t, t+dt
        
        # initialize velocities
        self.velocity = [self.nOfAtomsX3*[0] for i in xrange(4) ] # for times t-3dt/2, t-dt, t-dt/2, t
        for i in xrange(self.nOfAtoms):
                for j in xrange(3): self.velocity[-1][3*i+j] = random.gauss(0,self.AndersenSigmas[3*i+j])
                
        # set initial momentum to 0:
        totMomentum = [0,0,0]
        totMass = sum(self.masses)/3
        for i in xrange(self.nOfAtoms):
                for j in xrange(3): 
                    totMomentum[j] += self.masses[3*i+j]*self.velocity[-1][3*i+j]
        totMomentum = map(lambda x: x/totMass,totMomentum)
        for i in xrange(self.nOfAtoms):
                for j in xrange(3): self.velocity[-1][3*i+j] -= totMomentum[j]
        
        # set initial kinetic energy to nDim*KbT/2
        nDim = 3
        kinEnergy = 0
        for i in xrange(self.nOfAtoms):
                for j in xrange(3): kinEnergy += 0.5*self.masses[3*i+j]*self.velocity[-1][3*i+j]**2
        velocityScalingFactor = numpy.sqrt(nDim*BOLTZMANN_CONSTANT*TEMPERATURE/kinEnergy/2)
        for i in xrange(self.nOfAtoms):
            for j in xrange(3): self.velocity[-1][3*i+j] *= velocityScalingFactor
            
        '''kinEnergy = 0
        for i in xrange(self.nOfAtoms):
                for j in xrange(3): kinEnergy += 0.5*self.masses[3*i+j]*self.velocity[-1][3*i+j]**2
        print "kinEnergy = ",kinEnergy
        print "nDim*KbT/2 = ",nDim*BOLTZMANN_CONSTANT*TEMPERATURE/2'''
        
        
        self.acceleration = [ self.nOfAtomsX3*[0] for i in xrange(3) ] # for times t-dt, t, t+dt
        self.gradKsi = 3*[None] # for times t-dt, t, t+dt # TODO: nie potrzebuje az trzech, jak sadze
        self.gradU = 3*[None] # for times t-dt, t, t+dt # TODO: nie potrzebuje az trzech, jak sadze
        self.Zksi = 3*[None] # for times t-dt, t, t+dt  # TODO: nie potrzebuje az trzech, jak sadze
        self.ksiMomentum = 3*[None]
        self.hessian = 3*[None]
        
    def setNonVelocityAtribute(self,atribute,newValue):
        assert(len(atribute)==3)
        for i in xrange(2): atribute[i]=atribute[i+1]
        atribute[2] = newValue
        
    def setMasses(self,masses):
        self.masses = masses[:]
        self.invMasses = self.nOfAtomsX3*[0]
        for i in xrange(self.nOfAtomsX3): self.invMasses[i] = 1./self.masses[i]
        #dampingCoefficient = 1.0 # TODO: na razie przyjmuje, ze jest rowne 1, bez mozliwosci modyfikacji
        self.AndersenSigmas = [ numpy.sqrt(BOLTZMANN_CONSTANT*TEMPERATURE/mass) for mass in masses ] # TODO: sprawdz, czy faktycznie
        
    def calculateZksi(self,gradKsi):
        self.setNonVelocityAtribute(self.gradKsi,gradKsi[:])
        Zksi = scalarTrippleProduct(self.invMasses,gradKsi,gradKsi) # it's faster than numpy
        self.setNonVelocityAtribute(self.Zksi,Zksi)
        
    def getZksi(self):
        return self.Zksi[-1]
        
    def getGradKsi(self):
        return self.gradKsi[-1]
        
    def getGradU(self):
        return self.gradU[-1]
        
    def getCoords(self):
        return self.coords[-1]
        
    def setCoords(self,coords):
        self.setNonVelocityAtribute(self,coords,coords[:])
                
    def setCoords(self,coords):
        self.setNonVelocityAtribute(self.coords,coords[:])
        #self.setNonVelocityAtribute(self.nCoords,numpy.array(coords))
        
    def calculateGradU(self,ff,mol): # TODO: przepisz to tak, zeby korzystac z RELEVANT_COORDINATES
        ff.Setup(mol)
        energy = ff.Energy(True) # True -> gradients will be calculated and held 
        gradU = [ 0 ] * self.nOfAtomsX3
        for i in xrange(self.nOfAtoms):
            coordinate = 3*i 
            tmpGradient = ff.GetGradientPublically(mol.GetAtom(i+1))
            gradU[coordinate]   = -tmpGradient.GetX()
            gradU[coordinate+1] = -tmpGradient.GetY()
            gradU[coordinate+2] = -tmpGradient.GetZ()
                        
        self.setNonVelocityAtribute(self.gradU,gradU)
        return energy
            
            
    def AndersenIntegrator(self,mol,ff,h,reactionCoordinate):
        currCoords = self.coords[-1]
        currVelocity = self.velocity[-1]
        currAcc = self.acceleration[-1]
        
        # calculating next position:  
        nextCoords = self.nOfAtomsX3*[0]
        for i in xrange(self.nOfAtomsX3):
            nextCoords[i] = currCoords[i] + TIMESTEP*currVelocity[i] + currAcc[i]*TIMESTEP2/2 # nie jestem pewien tego currAcc
        
        for i in xrange(self.nOfAtomsX3):
            currVelocity[i] = currVelocity[i] + 0.5*TIMESTEP*currAcc[i] # that's the ACTUAL currentVelocity
            
        mol.SetCoordinates( openbabel.double_array(nextCoords) )   
        energy = self.calculateGradU(ff,mol)
        gradKsi = calcGrad(mol,reactionCoordinate,h)
        self.calculateZksi( gradKsi )
        
        Zksi = self.Zksi[-1]
        gradU = self.gradU[-1]
  
        if Zksi==0: constraintLambda = 0
        else: constraintLambda = (scalarTrippleProduct(self.invMasses,gradKsi,gradU)-calcVectorHessianVectorProduct(mol,reactionCoordinate,h,currVelocity)) / Zksi
        
        ###########################
        #constraintLambda = 0 # if set to 0 it will be a simple Andersen dynamics simulation
        ###########################

        currAcc = self.nOfAtomsX3*[0]                                
        for i in xrange(self.nOfAtomsX3):
            currAcc[i] = (-gradU[i] + constraintLambda*gradKsi[i])/self.masses[i] # TODO: upewnij sie, ze tu jest dzielenie przez mase
            
        # second velocity half-step
        for i in xrange(self.nOfAtomsX3):
            currVelocity[i] = currVelocity[i] + 0.5*TIMESTEP*currAcc[i] # that's the ACTUAL currentVelocity

        # Andersen thermostat
        for i in xrange(self.nOfAtoms):
            if random.random()<NUxTIMESTEP:
                for j in xrange(3): currVelocity[3*i+j] = random.gauss(0,self.AndersenSigmas[3*i+j])
                       
        self.setNonVelocityAtribute(self.coords,nextCoords)
        self.setNonVelocityAtribute(self.acceleration,currAcc)
        self.setNonVelocityAtribute(self.gradU,gradU)
        self.velocity[3] = currVelocity  
        
        return energy
         
                      
    def idleKSteps(self,K,mol,ff,reactionCoordinate,h):      
        for i in xrange(K): 
            self.AndersenIntegrator(mol,ff,h,reactionCoordinate)

            
        
    def calculateKsiMomentum(self):
        Zksi = self.Zksi[-1] 
        gradKsi = self.gradKsi[-1]
        velocity = self.velocity[-1] # v[t]
              
        if not (Zksi and gradKsi and velocity):
            self.setNonVelocityAtribute(self.ksiMomentum,None)
        else:
            ksiMomentum = 0
            for i in xrange(self.nOfAtomsX3):
                ksiMomentum += gradKsi[i] * velocity[i]
            ksiMomentum /= Zksi 
            
            self.setNonVelocityAtribute(self.ksiMomentum,ksiMomentum)
    
        
    def getTimeDerivativeOfKsiMomentum(self):
        '''
            Calculated at t-dt
        '''
        if all(self.ksiMomentum)==False: # all need to be != None
            return None
        return (self.ksiMomentum[-1]-self.ksiMomentum[-3])/2./TIMESTEP
            
            
    def calculateHessianDiagonal(self,mol,reactionCoordinate,h):        
        hessian =  [ [ 0 for i in xrange(self.nOfAtomsX3)  ] for j in xrange(self.nOfAtomsX3) ]
        currentReactionCoordinateValue = reactionCoordinate(mol)
        nOfAtoms = mol.NumAtoms()
        tmp = openbabel.doubleArray_frompointer(mol.GetCoordinates())

        for i in RELEVANT_COORDINATES:
            hessian[i][i] = calcSecondDerivativeInOneCoordinate(mol,nOfAtoms,reactionCoordinate,currentReactionCoordinateValue,h,i,tmp)
        self.setNonVelocityAtribute( self.hessian, hessian )
            
    def calculateHessian(self,mol,reactionCoordinate,h): # TODO: check if it's calculated at the right time
        self.setNonVelocityAtribute( self.hessian, calcHessian(mol,reactionCoordinate,h) )
        
        
    def getDivergenceMksiGradKsi(self):
        '''
            Needs to be further multiplied by - k_B * T ;
            Calculated at t-dt
        '''
        
        if all(self.ksiMomentum)==False:
            return None
        
        Zksi = self.Zksi[-2] # Zksi[t-dt]
        gradKsi = self.gradKsi[-2] # gradKsi[t-dt]
        hessian = self.hessian[-2] # hessian[t-dt]

        result = 0
        for i in xrange(self.nOfAtomsX3):
            tmp = 0
            for j in xrange(self.nOfAtomsX3):
                tmp -= gradKsi[j] * hessian[i][j] / self.masses[j]
            tmp *= 2 * gradKsi[i] / Zksi / Zksi
            tmp += hessian[i][i] / Zksi 
            result += tmp / self.masses[i]
            
        return result
        
        
    def getVelocity_gradMKsiGradKsi_Velocity(self):
        '''
            Calculated at t-dt
        '''
        
        if all(self.ksiMomentum)==False:
            return None
            
        Zksi = self.Zksi[-2] # Zksi[t-dt]
        gradKsi = self.gradKsi[-2] # gradKsi[t-dt]
        velocity = self.velocity[-3] # v[t-dt]
        hessian = self.hessian[-2] # hessian[t-dt]

        result = 0
        for i in xrange(self.nOfAtomsX3):
            for k in xrange(self.nOfAtomsX3):
                tmp = 0
                for j in xrange(self.nOfAtomsX3):
                    tmp -= gradKsi[j] * hessian[k][j] / self.masses[j]
                tmp *= 2 * gradKsi[i] / Zksi / Zksi
                tmp += hessian[i][k] / Zksi 
                result += velocity[i] * tmp * velocity[k] 
        
        return result
        
    def getMksi_gradU_invMasses_gradKsi(self):
        '''
            Calculated at t
        '''
                
        Zksi = self.Zksi[-1] # Zksi[t-dt]
        gradKsi = self.gradKsi[-1] # gradKsi[t-dt]
        gradU = self.gradU[-1]
        
        if Zksi==None: return None
        
        return ( scalarTrippleProduct(gradU,self.invMasses,gradKsi) * (Zksi**(-1.5)) ) # TODO: czy faktycznie **(-1.5) ?
                        
    def getKsiDerivativeOfFreeEnergy(self):
        '''
            Calculated at t
        '''
        
        Zksi = self.Zksi[-1] # Zksi[t-dt]
        gradKsi = self.gradKsi[-1] # gradKsi[t-dt]
        gradU = self.gradU[-1]
        hessian = self.hessian[-1]

        if Zksi==None: return None
                
        laplacianKsi = 0
        for i in xrange(len(hessian)): laplacianKsi += hessian[i][i]
        
        gradKsiLength = vecLength(gradKsi)
        
        dAdKsi = scalarProduct(gradKsi,gradU) - BOLTZMANN_CONSTANT*TEMPERATURE*(laplacianKsi - 2*multVecMatVec(gradKsi,hessian,gradKsi)/gradKsiLength/gradKsiLength)
        dAdKsi /= (gradKsiLength*gradKsiLength)
        dAdKsi *= Zksi**(-0.5)
        
        return dAdKsi
        
        
# END class DataAndCalculations



def vecLength(vector):
    return numpy.linalg.norm(vector)

def addVec(v,w):
    #assert(len(v)==len(w)) # slows down by 25%
    #return map(operator.add,v,w) 
    return [i+j for i,j in izip(v,w)] # this solution proved to be the fastest

def scalarProduct(x,y):
    #assert(len(x)==len(y))
    total = 0
    for i in xrange(len(x)):
        total += x[i] * y[i]
    return total

def scalarTrippleProduct(x,y,z):
    #assert(len(x)==len(y)==len(z))
    total = 0
    for i in xrange(len(x)):
        total += x[i] * y[i] * z[i]
    return total
    
def multVec(alpha,w):
    return [alpha*c for c in w]
    
def multVecMatVec(vec1,mat,vec2):
    assert len(vec1)==len(mat)==len(vec2)
    total = 0
    for i in xrange(len(vec1)):
        for j in xrange(len(vec2)):
            total += vec1[i] * mat[i][j] * vec2[j]
    return total

def calcDih(mol,WHICH_ATOMS):
    return mol.GetTorsion(*WHICH_ATOMS)
        
def calcDihWithCenterOfMass(mol,WHICH_ATOMS,WHICH_ATOMSBelongToGroup,masses):
    centerOfMass = [0.0,0.0,0.0]
    for atomIndex in WHICH_ATOMSBelongToGroup:
        mass = masses[atomIndex-1]
        vec = mol.GetAtom(atomIndex).GetVector()
        centerOfMass[0] += vec.GetX()*mass
        centerOfMass[1] += vec.GetY()*mass
        centerOfMass[2] += vec.GetZ()*mass
    centerOfMass = multVec(1.0/len(WHICH_ATOMSBelongToGroup),centerOfMass)
    theActualPosition = mol.GetAtom(WHICH_ATOMS[0]).GetVector()
    theActualPosition = [theActualPosition.GetX(),theActualPosition.GetY(),theActualPosition.GetZ()]
    mol.GetAtom(WHICH_ATOMS[0]).SetVector(centerOfMass[0],centerOfMass[1],centerOfMass[2])
    tor = mol.GetTorsion(*WHICH_ATOMS)
    mol.GetAtom(WHICH_ATOMS[0]).SetVector(theActualPosition[0],theActualPosition[1],theActualPosition[2])
    return tor  
            
'''def calcDirectionalGrad(mol,reactionCoordinate,versor,h,currentReactionCoordinate):
    molCopy = openbabel.OBMol(mol)
    tmp = openbabel.doubleArray_frompointer(mol.GetCoordinates())
    coords = [ tmp[i] for i in xrange(3*mol.NumAtoms()) ]
    coords = addVec(coords,multVec(h,versor))
    molCopy.SetCoordinates(openbabel.double_array(coords))
    newReactionCoordinate = reactionCoordinate(molCopy)
    if currentReactionCoordinate<-170 and newReactionCoordinate>170:
        return (currentReactionCoordinate-newReactionCoordinate+360)/h
    if newReactionCoordinate<-170 and currentReactionCoordinate>170:
        return (newReactionCoordinate-currentReactionCoordinate+360)/h
    return (newReactionCoordinate - currentReactionCoordinate)/h
'''
    
def shiftAndCalculateReactionCoordinate(mol,nOfAtoms,reactionCoordinate,vector,h,coords=None): # this function is called many times
    if coords==None:
        tmp = openbabel.doubleArray_frompointer(mol.GetCoordinates())
        coords = [ tmp[i] for i in xrange(3*nOfAtoms) ]
    coords = addVec(coords,multVec(h,vector))
    mol.SetCoordinates(openbabel.double_array(coords))
    return reactionCoordinate(mol)
        
    
def differenceQuotientForDihedralAngle(leftReactionCoordinate,rightReactionCoordinate,h):
    if leftReactionCoordinate<-170 and rightReactionCoordinate>170:
        return (leftReactionCoordinate-rightReactionCoordinate+360)/h/2
    if rightReactionCoordinate<-170 and leftReactionCoordinate>170:
        return (rightReactionCoordinate-leftReactionCoordinate+360)/h/2
    return (rightReactionCoordinate - leftReactionCoordinate)/h/2
    
def calcDirectionalGrad2(mol,reactionCoordinate,versor,h,currentReactionCoordinate):
    tmp = openbabel.doubleArray_frompointer(mol.GetCoordinates())
    for i in xrange(len(versor)): tmp[i] -= (h*versor[i])
    mol.SetCoordinates(tmp)
    leftReactionCoordinate = reactionCoordinate(mol)
    
    for i in xrange(len(versor)): tmp[i] += (2*h*versor[i])
    mol.SetCoordinates(tmp)
    rightReactionCoordinate = reactionCoordinate(mol)
    
    return differenceQuotientForDihedralAngle(leftReactionCoordinate,rightReactionCoordinate,h)
    
      
def calcGrad(mol,reactionCoordinate,h):
    numOfAtomsX3 = 3*mol.NumAtoms()
    handleToOpenBabelDoubleArray = openbabel.doubleArray_frompointer(mol.GetCoordinates())
    molCopy = openbabel.OBMol(mol)
    grad = numOfAtomsX3 * [ 0.0 ]
    for position in RELEVANT_COORDINATES:
        tmp = handleToOpenBabelDoubleArray[position]
        handleToOpenBabelDoubleArray[position] -= h 
        molCopy.SetCoordinates(handleToOpenBabelDoubleArray) # WATCH OUT!
        leftReactionCoordinate = reactionCoordinate(molCopy)
        handleToOpenBabelDoubleArray[position] += (2*h) 
        molCopy.SetCoordinates(handleToOpenBabelDoubleArray)
        rightReactionCoordinate = reactionCoordinate(molCopy)
        grad[position] = differenceQuotientForDihedralAngle(leftReactionCoordinate,rightReactionCoordinate,h)
        handleToOpenBabelDoubleArray[position] = tmp
    return grad


def calcVectorHessianVectorProduct(mol,reactionCoordinate,h,vector):
    '''
        According to:
            
           v^T * H(f(x)) * v := ( directionalGrad_f(x+h*v) - directionalGrad_f(x-h*v) ) / 2 / h
            
        with h->0
    '''
    
    molLeftCopy = openbabel.OBMol(mol) 
    molRightCopy = openbabel.OBMol(mol)
    nOfAtoms = mol.NumAtoms()
    tmp = openbabel.doubleArray_frompointer(mol.GetCoordinates())
    coords = [ tmp[i] for i in xrange(3*nOfAtoms) ]
    
    reactionCoordinateValueLeft = shiftAndCalculateReactionCoordinate(molLeftCopy,nOfAtoms,reactionCoordinate,vector,-h,coords) # Left => (-h)
    reactionCoordinateValueRight = shiftAndCalculateReactionCoordinate(molRightCopy,nOfAtoms,reactionCoordinate,vector,h,coords) # Right => (+h)
            
    vectorLength = vecLength(vector)
    if vectorLength==0: return 0
    versor = multVec(1./vectorLength,vector)
        
    return vectorLength * ( calcDirectionalGrad2(molRightCopy,reactionCoordinate,versor,h,reactionCoordinateValueRight) - calcDirectionalGrad2(molLeftCopy,reactionCoordinate,versor,h,reactionCoordinateValueLeft) ) / 2 / h # WATCH OUT!
    
    
def calcVectorHessianVectorProduct2(mol,reactionCoordinate,h,vector1,vector2):
    '''
        According to:
            
           v1^T * H(f(x)) * v2 := ( directionalGrad_f(x+h*v2) - directionalGrad_f(x-h*v2) ) / 2 / h
            
        with h->0
    '''
    
    molLeftCopy = openbabel.OBMol(mol) 
    molRightCopy = openbabel.OBMol(mol)
    nOfAtoms = mol.NumAtoms()
    tmp = openbabel.doubleArray_frompointer(mol.GetCoordinates())
    coords = [ tmp[i] for i in xrange(3*nOfAtoms) ]
    reactionCoordinateValueLeft = shiftAndCalculateReactionCoordinate(molLeftCopy,nOfAtoms,reactionCoordinate,vector2,-h,coords) # WATCH OUT!
    reactionCoordinateValueRight = shiftAndCalculateReactionCoordinate(molRightCopy,nOfAtoms,reactionCoordinate,vector2,h,coords) # WATCH OUT!    
            
    vector1Length = vecLength(vector1)
    vector2Length = vecLength(vector2)
    if vector1Length==0 or vector2Length==0: return 0
    versor1 = multVec(1./vector1Length,vector1)
        
    return vector1Length * ( calcDirectionalGrad2(molRightCopy,reactionCoordinate,versor1,h,reactionCoordinateValueRight) - calcDirectionalGrad2(molLeftCopy,reactionCoordinate,versor1,h,reactionCoordinateValueLeft) ) / 2 / h # WATCH OUT!
    

def calcSecondDerivativeInOneCoordinate(molCopy,nOfAtoms,reactionCoordinate,currentReactionCoordinateValue,h,coordinateNumber,handleToOpenBabelDoubleArray):  
    
    handleToOpenBabelDoubleArray[coordinateNumber] -= h
    molCopy.SetCoordinates(handleToOpenBabelDoubleArray)
    reactionCoordinateValueLeft = reactionCoordinate(molCopy)
    
    handleToOpenBabelDoubleArray[coordinateNumber] += (2*h)
    molCopy.SetCoordinates(handleToOpenBabelDoubleArray)    
    reactionCoordinateValueRight = reactionCoordinate(molCopy)
    
    handleToOpenBabelDoubleArray[coordinateNumber] -= h
    molCopy.SetCoordinates(handleToOpenBabelDoubleArray)    
    
    return (reactionCoordinateValueLeft + reactionCoordinateValueRight - 2*currentReactionCoordinateValue) / h / h
    
    
def calcSecondDerivativeMixedCoordinates(molCopy,nOfAtoms,reactionCoordinate,h,coordinateNumberA,coordinateNumberB,handleToOpenBabelDoubleArray):
    
    handleToOpenBabelDoubleArray[coordinateNumberA] += h
    handleToOpenBabelDoubleArray[coordinateNumberB] += h
    molCopy.SetCoordinates(handleToOpenBabelDoubleArray)
    reactionCoordinateValueRightRight = reactionCoordinate(molCopy)
    
    handleToOpenBabelDoubleArray[coordinateNumberA] -= (2*h)
    handleToOpenBabelDoubleArray[coordinateNumberB] -= (2*h)
    molCopy.SetCoordinates(handleToOpenBabelDoubleArray)
    reactionCoordinateValueLeftLeft = reactionCoordinate(molCopy)
    
    handleToOpenBabelDoubleArray[coordinateNumberA] += (2*h)
    molCopy.SetCoordinates(handleToOpenBabelDoubleArray)
    reactionCoordinateValueRightLeft = reactionCoordinate(molCopy)
    
    handleToOpenBabelDoubleArray[coordinateNumberA] -= (2*h)
    handleToOpenBabelDoubleArray[coordinateNumberB] += (2*h)
    molCopy.SetCoordinates(handleToOpenBabelDoubleArray)
    reactionCoordinateValueLeftRight = reactionCoordinate(molCopy)
    
    handleToOpenBabelDoubleArray[coordinateNumberA] += h
    handleToOpenBabelDoubleArray[coordinateNumberB] -= h
    molCopy.SetCoordinates(handleToOpenBabelDoubleArray)
        
    return (reactionCoordinateValueRightRight + reactionCoordinateValueLeftLeft - reactionCoordinateValueRightLeft - reactionCoordinateValueLeftRight) / h / h / 4
    
def calcHessian(mol,reactionCoordinate,h):
    nOfAtoms = mol.NumAtoms()
    hessian =  [ [ 0 for i in xrange(3*nOfAtoms)  ] for j in xrange(nOfAtoms*3) ]
    handleToOpenBabelDoubleArray = openbabel.doubleArray_frompointer(mol.GetCoordinates())
    molCopy = openbabel.OBMol(mol)
    
    currentReactionCoordinateValue = reactionCoordinate(mol)
    for i in xrange(len(RELEVANT_COORDINATES)):
        coordinateA = RELEVANT_COORDINATES[i]
        hessian[coordinateA][coordinateA] = calcSecondDerivativeInOneCoordinate(molCopy,nOfAtoms,reactionCoordinate,currentReactionCoordinateValue,h,coordinateA,handleToOpenBabelDoubleArray)
        for j in xrange(i+1,len(RELEVANT_COORDINATES)):
            coordinateB = RELEVANT_COORDINATES[j]
            hessian[coordinateA][coordinateB] = calcSecondDerivativeMixedCoordinates(molCopy,nOfAtoms,reactionCoordinate,h,coordinateA,coordinateB,handleToOpenBabelDoubleArray)
            hessian[coordinateB][coordinateA] = hessian[coordinateA][coordinateB]
                    
    return hessian 


def assignToSet(setA,setB,interactionAtoms):
    if sum([ atomId in setA for atomId in interactionAtoms])/len(interactionAtoms)==1:
        return 1
    if sum([ atomId in setB for atomId in interactionAtoms])/len(interactionAtoms)==1:
        return 3
    return 2


def parseVariousForces2(variousForces,gradKsi,Zksi,masses):
    products = [] # this will be the output of this method
    lines = variousForces.split("\n")

    for line in lines:
        words = line.split()
        if len(words)==0: continue
        
        products.append([])

        whichInteraction = words[0]
        products[-1].append(whichInteraction)
        
        if whichInteraction in ("oop","strbnd"):
            pass
        
        if whichInteraction in ("ele","vdw","bonds"):
            which = 3
        elif whichInteraction=="angle":
            which = 4
        elif whichInteraction=="tors":
            which = 5
        else: raise Exception("Nieznany typ oddzialywan")

        atoms = [int(word) for word in words[1:which]]
        products[-1].append(atoms)
        contributingForces = [ float(word) for word in words[which:] ]
        for i in xrange(len(atoms)):
            atomNumberCoordinate = 3*(atoms[i]-1) # I assume here that gradKsi is relevant starting from the first (zeroth) position
            forceCoordinate = 3*i
            product = 0
            for coo in xrange(3):
                # < grad_U . m_ksi . M^(-1) . grad_ksi > ~= Z_ksi^(-1) . grad_U . M^(-1) . grad_ksi
                # because 
                #    m_ksi = Z_ksi^(-1)
                #
                product -= gradKsi[atomNumberCoordinate+coo] * contributingForces[forceCoordinate+coo] / masses[atomNumberCoordinate+coo] / Zksi
            products[-1].append(product) 
            
            
    return products
            
            
def newInteractionTypeDictionary(initializer):
    return {"ele":initializer[:],"vdw":initializer[:],"bonds":initializer[:],"angle":initializer[:],"tors":initializer[:]}#,"oop":initializer[:],"strbnd":initializer[:]}
     
     
def fillOutProductsMatrices(products,productsMatrices):
    for product in products:
        whichInteraction = product[0]
        atoms = product[1]
        actualProducts = product[2:]
        assert(len(atoms)==len(actualProducts))
        nOfAtoms = len(atoms)
        #print "\ninteraction type = %s" % whichInteraction
        for i in xrange(nOfAtoms-1):
            for j in xrange(i+1,nOfAtoms):
                #if atoms[i]==1 and atoms[j]==5 and whichInteraction=="ele":
                #    print "atoms: %d <---> %d, products = %s" % (atoms[i]-1,atoms[j]-1,str(["%.2f" % k for k in actualProducts]))
                tmp = (actualProducts[i]+actualProducts[j])/(nOfAtoms-1)
                productsMatrices[whichInteraction][atoms[i]-1][atoms[j]-1] += tmp # ISTOTNE: atomy numerowane sa od 1, ale macierze sa od 0
                productsMatrices[whichInteraction][atoms[j]-1][atoms[i]-1] += tmp
 
           
def printOutProductsMatrices(productsMatrices,ticks,step=1,dihedral=0.0,fileName=None,limit=50):
    pylab.clf()
    fig, axes = plt.subplots(nrows=2, ncols=3)
    
    plt.setp(axes)
    for interaction, ax in  zip(["ele","vdw","bonds","angle","tors"],axes.flat):
        im = ax.matshow(productsMatrices[interaction],interpolation='none',cmap=pylab.cm.seismic,vmin=-limit,vmax=limit)
        plt.sca(ax)
        plt.tick_params(labelsize=8)
        plt.xticks(range(len(ticks)),ticks)
        plt.yticks(range(len(ticks)),ticks)
        plt.xlabel("%s, %.3f" % (interaction,dihedral))
    #plt.colorbar()
    if fileName!=None: 
        print fileName
        fig.savefig(fileName) 
    else: fig.savefig("productsMatrices_%05d.png" % step)


def printOutAngles(angles):
    nOfAngles = len(angles)
    h = 0.001
    reactionCoordinate = lambda x: calcDih(x,WHICH_ATOMS)
    ff = openbabel.OBForceField.FindForceField(FF_NAME)
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol2", "mol2")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, MOL2_FILE_NAME)
    nOfAtoms = mol.NumAtoms()
    
    ticks = []
    for i in xrange(nOfAtoms):
        ticks.append(mol.GetAtom(i+1).GetType()[0]+str(i))
    
    
    for angle in angles:
        mol.SetTorsion(WHICH_ATOMS[0],WHICH_ATOMS[1],WHICH_ATOMS[2],WHICH_ATOMS[3],angle*math.pi/180.)
        ff.Setup(mol)
        energy = ff.Energy(True)
        dihedral = mol.GetTorsion(*WHICH_ATOMS)
        print "E(%.2f) = %.2f" % (dihedral,energy)
        obConversion.WriteFile(mol,"nh2-ii_ANGLE_%d.mol2" % int(angle))
        
    
def initProductsMatrices(nOfAtoms):
    productsMatrices = {}
    for key in ["ele","vdw","bonds","angle","tors"]:
        productsMatrices[key] = [ [0.0 for j in xrange(nOfAtoms) ] for i in xrange(nOfAtoms) ] # graph prototypes shall be storred in this data structure
    return productsMatrices
    
def addToAllProductsMatrices(allProductsMatrices_ACCUMULATED,productsMatrices,nOfAtoms,ZksiToMinusHalf):
    for key in productsMatrices.keys():
        for i in xrange(nOfAtoms):
            for j in xrange(nOfAtoms):
                allProductsMatrices_ACCUMULATED[key][i][j] += productsMatrices[key][i][j]*ZksiToMinusHalf

    
       
def checkHessians(mol,ff,angle,nOfSamples,nOfStepsInBetween,twistPeriod):
    
    assert mol != None and ff != None
    mol.SetTorsion(WHICH_ATOMS[0],WHICH_ATOMS[1],WHICH_ATOMS[2],WHICH_ATOMS[3], angle*math.pi/180.)
    
    # constants:
    h = 0.0001 # 10^-4 seemed best
    reactionCoordinate = lambda x: calcDih(x,WHICH_ATOMS)
    
    vector = [ random.random() for i in xrange(3*mol.NumAtoms()) ]
    print "my estimate: ",calcVectorHessianVectorProduct(mol,reactionCoordinate,h,vector)
    start = time.clock() 
    for i in xrange(1000): calcVectorHessianVectorProduct(mol,reactionCoordinate,h,vector)
    print "which takes %f sec" % (time.clock()-start)
    
    hessian = calcHessian(mol,reactionCoordinate,h)
    total = 0
    for i in xrange(len(hessian)):
        for j in xrange(len(hessian)):
            total += vector[i] * vector[j] * hessian[i][j]
    print "true estimate: ",total
    start = time.clock()
    for k in xrange(10):
        hessian = calcHessian(mol,reactionCoordinate,h)
        total = 0
        for i in xrange(len(hessian)):
            for j in xrange(len(hessian)):
                total += vector[i] * vector[j] * hessian[i][j]
    print "which takes %f sec" % (100*(time.clock()-start))

                                                                                        
def samplingWithConstantReactionCoordinate(mol,ff,angle,nOfSamples,nOfStepsInBetween,twistPeriod):
    
    assert mol != None and ff != None
        
    random.seed(13) 
    
    # constants:
    h = 0.0001 # 10^-4 seemed best
    def reactionCoordinate(x): 
        return calcDih(x,WHICH_ATOMS)

    mol.SetTorsion(WHICH_ATOMS[0],WHICH_ATOMS[1],WHICH_ATOMS[2],WHICH_ATOMS[3], angle*math.pi/180.)
    nOfAtoms = mol.NumAtoms()
    tmp = openbabel.doubleArray_frompointer(mol.GetCoordinates())
    coords = 3*nOfAtoms*[0]
    for i in xrange(3*nOfAtoms): coords[i] = tmp[i]     
            
    masses = [ 0 ] * 3*nOfAtoms
    for i in xrange(nOfAtoms): masses[3*i] = masses[3*i+1] = masses[3*i+2] = 1000 * mol.GetAtom(i+1).GetAtomicMass() # TODO: tak jest w openbabelu
    
    # data storage and calculation
    dataCalculator = DataAndCalculations(coords,masses)    
    
    
    # data extracted from the simulation will be storred in the following data structures:
    index = 0
    lowestEnergy = 100000000
    configurationWithLowestEnergy = None
    dihedrals = MEMORY_BUFFER_SIZE * [ None ]
    energies = MEMORY_BUFFER_SIZE * [ None ]
    Zksis = MEMORY_BUFFER_SIZE * [ None ]
    mksi_gradU_invMasses_gradKsis = MEMORY_BUFFER_SIZE * [ None ]
    dAdKsis = MEMORY_BUFFER_SIZE * [ None ]
    #TODO: allProductsMatrices = MEMORY_BUFFER_SIZE * [ None ]
    ZksiToMinusHalf_ACCUMULATED = 0
    mksi_gradU_invMasses_gradKsi_ACCUMULATED = 0
    dAdKsi_ACCUMULATED = 0
    allProductsMatrices_ACCUMULATED = initProductsMatrices(nOfAtoms)
    
    # quality / speed testing:
    noShiftTimes = []
    dihedrals_HNCC = []
    shiftTime = 0
    nOfAccepts = 0
    
    # tmp for simulation purposes:
    coords = [ 0 ] * 3*nOfAtoms
        
    convergenceOutput = open("convOutput_angle=%.2f_timestep=%s_twistPeriod=%d_nOfSamples=%d_binWidth=%s_idleSteps=%d.dat" % (angle,str(TIMESTEP),twistPeriod,nOfSamples,str(BIN_WIDTH),IDLE_STEPS),"w")
    convergenceOutput.write( "ZksiToMinusHalf_ACCUMULATED\tdAdKsi-<m_ksi gradU.M^(-1).gradKsi>\tdAdKsi\t\t<m_ksi gradU.M^(-1).gradKsi>\t-<v.grad(m_ksi gradKsi).v>\t-kT<m_ksi div(M^(-1).gradKsi)>\n" )

    gatherStatisticsFlag = False
    start = time.clock() 
    
    if LARGE_OUTPUT:
        generalOutputStream = gzip.open(PATH_TO_LARGE_OUTPUT+"generalOutput_angle=%.2f_timestep=%s_twistPeriod=%d_nOfSamples=%d_binWidth=%s_idleSteps=%d.gz" % (angle,str(TIMESTEP),twistPeriod,nOfSamples,str(BIN_WIDTH),IDLE_STEPS),"wb")
        generalOutputStream.write("dihedral\t\t\tdihedral_HNCC\t\t\tenergy\t\t\tdAdKsi\t\t\t<...gradU...>\t\t\tZksi\n")
        fullMatricesOutputStreams = {}
        for key in allProductsMatrices_ACCUMULATED.keys():
            fullMatricesOutputStreams[key] = gzip.open(PATH_TO_LARGE_OUTPUT+"matrices_angle=%.2f_timestep=%s_twistPeriod=%d_nOfSamples=%d_binWidth=%s_idleSteps=%d_interaction=%s.gz" % (angle,str(TIMESTEP),twistPeriod,nOfSamples,str(BIN_WIDTH),IDLE_STEPS,key),"w")
         
    allProductsMatrices = MEMORY_BUFFER_SIZE * [ None ]
         
    for iteration in xrange(nOfSamples*nOfStepsInBetween):
        
        if gatherStatisticsFlag==False:
            dataCalculator.idleKSteps(IDLE_STEPS,mol,ff,reactionCoordinate,h)
            gatherStatisticsFlag = True
        else:
            #gradKsi = calcGrad(mol,reactionCoordinate,h)
            #dataCalculator.calculateZksi( gradKsi )
            #energy = dataCalculator.calculateGradU(ff,mol) this is now in the AndersenIntegrator method
                          
            # ONE STEP OF SIMULATION:
            energy = dataCalculator.AndersenIntegrator(mol,ff,h,reactionCoordinate)
            
            if numpy.abs(reactionCoordinate(mol)-angle)>BIN_WIDTH:
                print "dihedral = ",reactionCoordinate(mol)
                mol.SetTorsion(WHICH_ATOMS[0],WHICH_ATOMS[1],WHICH_ATOMS[2],WHICH_ATOMS[3],angle*math.pi/180.)
                tmp = openbabel.doubleArray_frompointer(mol.GetCoordinates())
                for i in xrange(3*nOfAtoms): coords[i] = tmp[i]     
                dataCalculator.reset( coords )
                print "shift!!!!"      
                newAngle = reactionCoordinate(mol)
                print "now it's dihedral = ",newAngle
                if numpy.abs(newAngle-angle)>BIN_WIDTH:
                    print "WARNING!!!"
                    obConversion.WriteFile(mol,"ERROR_ANGLE=%.3f_ENERGY=%.5f.mol2" % (angle,energy))
                    os.exit()
                noShiftTimes.append(shiftTime)
                shiftTime = 0
                gatherStatisticsFlag = False
                continue
            else:
                shiftTime += 1
            
            # MEMORIZE LOWEST-ENERGY CONFIGURATION:
            if energy<lowestEnergy:
                lowestEnergy = energy
                configurationWithLowestEnergy = openbabel.OBMol(mol)

            # CALCULATE SECOND DERIVATIVES OF KSI:
            dataCalculator.calculateHessianDiagonal(mol,reactionCoordinate,h)
            Zksi = dataCalculator.getZksi()
            gradKsi = dataCalculator.getGradKsi()              
                  
            if Zksi==0: print gradKsi
            
            tmp_mksi_gradU_invMasses_gradKsi = dataCalculator.getMksi_gradU_invMasses_gradKsi() # should be multiplied by Zksi^(-0.5) in this method
            tmp_dAdKsi = dataCalculator.getKsiDerivativeOfFreeEnergy() # is multiplied by Zksi^(-0.5) in this method
            if tmp_mksi_gradU_invMasses_gradKsi!=None and tmp_dAdKsi!=None:
                mksi_gradU_invMasses_gradKsi_ACCUMULATED += tmp_mksi_gradU_invMasses_gradKsi
                dAdKsi_ACCUMULATED += tmp_dAdKsi
                ZksiToMinusHalf = Zksi**(-0.5)
                ZksiToMinusHalf_ACCUMULATED += ZksiToMinusHalf
                dAdKsis[index] = tmp_dAdKsi / ZksiToMinusHalf 
                dihedrals_HNCC.append( mol.GetTorsion(2,1,5,6) ) 
                dihedral = reactionCoordinate(mol)
                energies[index] = energy  
                dihedrals[index] = dihedral    
                Zksis[index] = Zksi 
                #mksi_gradU_invMasses_gradKsis[index] = tmp_dAdKsi / ZksiToMinusHalf # CO TO BYLO?
                mksi_gradU_invMasses_gradKsis[index] = tmp_mksi_gradU_invMasses_gradKsi / ZksiToMinusHalf
                productsMatrices = initProductsMatrices(nOfAtoms)
                variousForces = ff.GetVariousForces()
                products = parseVariousForces2(variousForces,gradKsi,Zksi,masses)
                fillOutProductsMatrices(products,productsMatrices)
                addToAllProductsMatrices(allProductsMatrices_ACCUMULATED,productsMatrices,nOfAtoms,ZksiToMinusHalf) # it's where productsMatrices are multiplied by ZksiToMinusHalf
                allProductsMatrices[index] = productsMatrices
                if numpy.abs(dihedrals[index]-angle)>BIN_WIDTH:
			print "gathering WRONG stats for dihedral = ",dihedrals[-1]
			sys.exit(1)
                if energies[-1] != None:
                    if LARGE_OUTPUT:
                        for i in xrange(MEMORY_BUFFER_SIZE):
				generalOutputStream.write( 
					"%f\t\t\t%f\t\t\t%f\t\t\t%f\t\t\t%f\t\t\t%f\n" % (dihedrals[i],
											dihedrals_HNCC[i],
											energies[i],
											dAdKsis[i],
											mksi_gradU_invMasses_gradKsis[i],
											Zksis[i]) )
                    index = 0
                    dihedrals = MEMORY_BUFFER_SIZE * [ None ]
                    energies = MEMORY_BUFFER_SIZE * [ None ]
                    Zksis = MEMORY_BUFFER_SIZE * [ None ]
                    mksi_gradU_invMasses_gradKsis = MEMORY_BUFFER_SIZE * [ None ]
                    dAdKsis = MEMORY_BUFFER_SIZE * [ None ]
                    for i in xrange(MEMORY_BUFFER_SIZE):
                        if LARGE_OUTPUT:
                            for key in fullMatricesOutputStreams.keys():
                                for row in allProductsMatrices[i][key]:
                                    for elementInRow in row:
                                        if abs(elementInRow)>0.000001: fullMatricesOutputStreams[key].write( "%f " % elementInRow )
                                        else: fullMatricesOutputStreams[key].write( "0 " )
                                    fullMatricesOutputStreams[key].write( "\n" )
                    allProductsMatrices = MEMORY_BUFFER_SIZE * [ None ]
                else:
                    index += 1
                
            if numpy.mod(iteration+1,int(0.001*nOfStepsInBetween*nOfSamples))==0:
                gradKsi = dataCalculator.getGradKsi()   
		percent = int(1000.0*iteration/nOfSamples/nOfStepsInBetween)/10.
                elapsed = (time.clock() - start)
                start = time.clock()
                outTmp = "\t%.8f\t\t\t%.8f\t\t\t%.8f\t\t%.8f\t\t\t%.8f\t\t\t%.8f\t\t(%s percents, one took %s secs)\n" % (ZksiToMinusHalf_ACCUMULATED,(dAdKsi_ACCUMULATED-mksi_gradU_invMasses_gradKsi_ACCUMULATED)/ZksiToMinusHalf_ACCUMULATED,dAdKsi_ACCUMULATED/ZksiToMinusHalf_ACCUMULATED,mksi_gradU_invMasses_gradKsi_ACCUMULATED/ZksiToMinusHalf_ACCUMULATED,0,0, str(percent),str(elapsed))
                convergenceOutput.write( outTmp )
                print outTmp
                                    
            # SETTING NEW COORDINATES:
            #mol.SetCoordinates( openbabel.double_array(nextCoords) )   now in the AndersenIntegrator
            
            # CHECKING THE H-N-C-C DIHEDRAL ANGLE:
            if numpy.mod(iteration+1,twistPeriod)==0:
                oldAngle = mol.GetTorsion(2,1,5,6)
                newAngle = random.random()*359.999-179.999
                mol.SetTorsion(6,5,1,2,newAngle*math.pi/180.)
                ff.Setup(mol)
                newEnergy = ff.Energy()
                if newEnergy <= energy or math.exp((energy-newEnergy)/BOLTZMANN_CONSTANT/TEMPERATURE) > random.random():
                    tmp = openbabel.doubleArray_frompointer(mol.GetCoordinates())
                    coords = 3*nOfAtoms*[0]
                    for i in xrange(3*nOfAtoms): coords[i] = tmp[i]     
                    dataCalculator.reset( coords )
                    nOfAccepts += 1
                    gatherStatisticsFlag = False
                else:
                    mol.SetTorsion(6,5,1,2,oldAngle*math.pi/180.)
                    
                    

      
    if LARGE_OUTPUT:
        for i in xrange(index): generalOutputStream.write( "%f\t\t\t%f\t\t\t%f\t\t\t%f\t\t\t%f\t\t\t%f\n" % (dihedrals[i],dihedrals_HNCC[i],energies[i],dAdKsis[i],mksi_gradU_invMasses_gradKsis[i],Zksis[i]) )
        generalOutputStream.close()        
        for key in fullMatricesOutputStreams.keys():
            for i in xrange(index):
                for row in allProductsMatrices[i][key]:
                    for elementInRow in row:
                        if abs(elementInRow)>0.000001: fullMatricesOutputStreams[key].write( "%f " % elementInRow )
                        else: fullMatricesOutputStreams[key].write( "0 " )
                    fullMatricesOutputStreams[key].write( "\n" )
            fullMatricesOutputStreams[key].close()
    
    check_ACCUMULATED = 0
    for key in allProductsMatrices_ACCUMULATED.keys():
        for i in xrange(nOfAtoms):
            for j in xrange(nOfAtoms):
                allProductsMatrices_ACCUMULATED[key][i][j] /= ZksiToMinusHalf_ACCUMULATED
                check_ACCUMULATED += allProductsMatrices_ACCUMULATED[key][i][j]
                
                
    fileNames = glob.glob("./lowestEnergyConfigurations/ii-nh3*ANGLE=%.3f*" % angle)
    if len(fileNames)>1:
        print "what?"
    elif len(fileNames)==1:
        lowestEnergySoFar = float(fileNames[0].split("ENERGY=")[1].split(".mol2")[0])
        if lowestEnergySoFar>lowestEnergy:
            obConversion.WriteFile(configurationWithLowestEnergy,"./lowestEnergyConfigurations/ii-nh3_forGAFF_ANGLE=%.3f_ENERGY=%.5f.mol2" % (angle,lowestEnergy))
            os.remove(fileNames[0])
    elif len(fileNames)==0:
        obConversion.WriteFile(configurationWithLowestEnergy,"./lowestEnergyConfigurations/ii-nh3_forGAFF_ANGLE=%.3f_ENERGY=%.5f.mol2" % (angle,lowestEnergy))
    
    print "check_ACCUMULATED = ",check_ACCUMULATED/2
    print "mksi_gradU_invMasses_gradKsi_ACCUMULATED = ",mksi_gradU_invMasses_gradKsi_ACCUMULATED/ZksiToMinusHalf_ACCUMULATED

    
    print "P(acc) = ",  float(nOfAccepts) / nOfSamples / nOfStepsInBetween * twistPeriod
    print "mean shift time = ",numpy.mean(noShiftTimes)
    
    convergenceOutput.close()
    
    return allProductsMatrices_ACCUMULATED, noShiftTimes, ZksiToMinusHalf_ACCUMULATED
    
    
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
    import sys
    angles = sys.argv[1:]
    
    # OPENBABEL STUFF:
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol2", "mol2")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, MOL2_FILE_NAME)
    ff = openbabel.OBForceField.FindForceField(FF_NAME)    
    #printOutAngles([-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180])

    
    # SIMULATION STUFF:    
    nOfSamples = 10**5
    twistPeriod = 1000 # TODO: zastanow sie
    nOfStepsInBetween = 1

    
    import pickle
    for angle in angles:
        angle = int(100*float(angle))/100.
        print "angle = ",angle
        #checkHessians(mol,ff,angle,nOfSamples,nOfStepsInBetween,twistPeriod)
        allProductsMatrices_ACCUMULATED, noShiftTimes, ZksiToMinusHalf_ACCUMULATED = samplingWithConstantReactionCoordinate(mol,ff,angle,nOfSamples,nOfStepsInBetween,twistPeriod)
        
        
    for interaction in allProductsMatrices_ACCUMULATED.keys():
        stream = open("productsMatrices/%s/angle=%.2f_interaction=%s_timestep=%s_twistPeriod=%d_nOfSamples=%d_binWidth=%s_idleSteps=%d.dat" % (PREFIX,angle,interaction,str(TIMESTEP),twistPeriod,nOfSamples,str(BIN_WIDTH),IDLE_STEPS),"w")
        for line in allProductsMatrices_ACCUMULATED[interaction]:
            stream.write(" ".join(["%f\t" % number for number in line])+"\n")
        stream.close()

    
