#!/usr/bin/python

import sys
import string
import re
import itertools

def findOccurences(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]
    
def missStrMatch(str1,str2):
    '''compares a noMiss str witha at-least-one-missStr
        returns 1 if matches 0 otherwise'''
    #print str1
    #print str2    
    if (str1.find('_') != -1):
        r = 0
        return r
    r = int(str1 == str2)
    if(r==0):
        s2 = findOccurences(str2,'_')
        #print 's2 = ',s2
        if(s2):
            L = list(str1)
            for i in s2:
                L[i] = '_'
            str1 = ''.join(L)
            #print 'str2 = ',str2, 'str1 = ',str1
            r = int(str1==str2)
    return r
    
def explainData(baseSet,data,nPE):
    ''' can explained by only 2 Haplotypes '''
    m = 0
    r = 1
    for j in range(nPE): 
        m1 = missStrMatch(str(baseSet[0]),str(data[j]))
        m2 = missStrMatch(str(baseSet[1]),str(data[j]))
        #print str(baseSet[0]),' , ',str(baseSet[1]),' <=> ',str(data[j])
        #print 'm1 = ',m1,' m2 = ',m2
        if( (m1 == 0) & (m2 == 0) ):
            r = 0
            return r
            
    return r
            
#######################################################
########### main function #############################
#print 'Number of arguments:', len(sys.argv)-1, 'arguments.'
#print 'Argument List:', str(sys.argv)
if len(sys.argv) != 3:
    print "Usage: python hcf_filter_file.py <input file> <output file>\n"
    sys.exit()
else:
    print "Running hcf filter on ", str(sys.argv[1])

TOT_FLD = 9
LHV_MIN_CNT = 3
SIG_HAP_CNT_POS = 3
SIG_HAP_POS = 4
REF_POS = 2
PE_POS = 7
DATA_POS = 8
MIN_FREQ = 1

fName_in = str(sys.argv[1])
fName_out = str(sys.argv[2])

fid2=open(fName_out, 'w')
        
with open(fName_in, 'r') as fid1:
    for tStr1 in fid1:
        if((tStr1[0] != '#') & (len(tStr1) > 1)):
            #print tStr1
            T1 = tStr1.split('\t')
            if(len(T1) == TOT_FLD):
                nSig = int(T1[SIG_HAP_CNT_POS])
                if(nSig >= LHV_MIN_CNT):
                    nRG = len(T1[REF_POS]); 
                    T2 = re.split('[(|)|,]',T1[SIG_HAP_POS])
                    
                    calledHap = []
                    calledHapStr = ""
                    for i in range(nSig):
                        k = 3*i
                        #print str(T2[k])
                        calledHap.append(str(T2[k])) 
                        calledHapStr = "%s,%s" % (calledHapStr, calledHap[i])
                    L=list(calledHapStr)
                    L[0]=''
                    calledHapStr = ''.join(L)
                    #print calledHapStr
                    tStr2 = T1[PE_POS]
                    T3 = tStr2.split('=')
                    nPE = int(T3[1])
                    T4 = re.split('[;|=]',T1[DATA_POS])
                    j = 0
                    data = []
                    freq = []
                    for i in range(1,nPE+1):
                        id1 = 2*(i-1)
                        id2 = id1+1
                        f = int(T4[id2])
                        d = str(T4[id1])
                        if( f > MIN_FREQ):
                            j = j+1
                            data.append(d)
                            freq.append(f)
                    nPE = j        
                    #print 'nPE = ',nPE
                    #print data
                    #print freq    
                    #print 'All necessaty parsing done!!'
                    M = list(itertools.combinations(range(0,nSig),2))
                    L = len(M)
                    bFlagLHV = 0
                    baseSet = []
                    for i in range(L):
                        baseSet = []
                        baseSet.append(calledHap[M[i][0]])
                        baseSet.append(calledHap[M[i][1]])
                        #print baseSet
                        #print data
                        r = explainData(baseSet,data,nPE)
                        
                        if(r == 1):
                            tmpStr1 = "%s\tSHRINK=1\t(%s,%s)\n" % (tStr1.rstrip(), baseSet[0],baseSet[1])
                            fid2.write(tmpStr1)
                            break
                        else:
                            bFlagLHV = bFlagLHV + r
                            
                    if((r != 1) & (bFlagLHV == 0)):
                         tmpStr2 = "%s\tSHRINK=0\t(%s)\n" % (tStr1.rstrip(), calledHapStr)
                         fid2.write(tmpStr2)
                else:
                    fid2.write(tStr1)
            else:                    
                print 'Not correct type of hcf file !!';
                break
        else:
            fid2.write(tStr1)
    
    print 'input file ',fName_in,' is done!!'    
    fid2.close()                  
                                  
 #################################    