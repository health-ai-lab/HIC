import math
import sys
from statsmodels.stats.proportion import proportions_ztest
import numpy as np
import pickle

""" 

Checks if a is a decendant of b

ddepth: the depth of each node like ddepth['root']=0,ddepth['250']=1, and ddepth['250.0']=2,etc
p: dictionary where the key is a node x and the value is the parent of node x (e.g. p['250.00']='250.0')

"""

def isDes(a,b,p,ddepth):
        #is b a des of a (a is ancestor)
        if ddepth[b]<ddepth[a]:
                return False
        for i in range(ddepth[b]-ddepth[a]):
                b=p[b] # get the parent of b until the depths are the same
        if a==b:
                return True
        return False

""" 

Gets all leaf nodes in tree d based on icd9 f

f: current icd9 code
d: dictionary of all icd9 codes
dans: the dictionary that stores all leaf nodes

"""

def getleaf(f,d,dans):
        if f in d:
                for q in d[f]:
                        if q is not None:
                                currans = getleaf2(q,d,dans)
                                if currans  is not None:
                                        dans[currans]=1
        else:
                if f is not None:
                        return f
        return dans

""" 

getleaf helper function to iterate over all children

f: current icd9 code
d: dictionary of all icd9 codes
dans: the dictionary that stores all leaf nodes


"""

def getleaf2(f,d,dans):
        if f in d:
                for q in d[f]:
                        if q is not None:
                                currans = getleaf2(q,d,p,dans)
                                if currans  is not None:
                                        dans[currans]=1
        else:
                if f is not None:
                        return f

""" 

Fills the depth of each icd9 code

f: current icd9 code
d: dictionary of all icd9 codes
ddepth: the depth of each node like ddepth['root']=0,ddepth['250']=1, and ddepth['250.0']=2,etc [will be filled in]

"""

def fillDepth(f,d,ddepth):
        if f in d:
                for q in d[f]:
                        ddepth[q]=ddepth[f]+1
                        fillDepth(q,d,ddepth)

""" 

Imports the weights based on equations 6 and 7 (for HIC) 
The file should have 2 rows:
branch,X
tree,Y

weightsfn: path to file with the branch and tree weight

"""

def importWeights(weightsfn):
        f=open(weightsfn,"r")
        f1=f.readlines()
        for x in f1:
                a=x.split(",")
                if a[0].lower()=="branch":
                        weightbranch=float(a[1])
                elif a[0].lower()=="tree"
                        weighttree=float(a[1])
                else:
                        continue
        return weightbranch, weighttree

""" 

Imports an ontology in the form of a 2 column seperated list
first,second means that second is a decendant of first ex: 250.0,250.00
all level 1 icd9 codes must have 'root' as a parent ex: root,250
all comments beggining with # are ignored

ontofn: path to ontology (icd9 hierarchy)

"""

def importOnto(ontofn):
        f=open(ontofn,"r")
        f1=f.readlines()
        d={}
        p={}
        for x in f1:

                if len(x)<=1 or x[0]=="#":
                        continue
                a=x.split(",")
                if a[1][-1]=='\n':
                        a[1]=a[1][:-1] #take out the \n
                if a[0] in d:
                        d[a[0]].append(a[1])
                else:
                        d[a[0]]=[a[1]]
                p[a[1]]=a[0]
        ddepth={}
        ddepth['root']=0
        fillDepth('root',d,ddepth)
        highestLevel=0 #init to find the longest depth
        for kdep in ddepth:
                if ddepth[kdep]>highestLevel:
                        highestLevel=ddepth[kdep]
        return d,p,ddepth,highestLevel

""" 

Imports an ontology in the form of a dictionary num with keys 'all' and 'dead'
num[x]['all']: signifying the amount of people who have icd9 code x
num[x]['dead']: signifying the amount of people dead at that icd9 code x

ontonumfn: path to file with the format icd9,amount of people with that icd9 code,amount of people with that icd9 code who have died
p: dictionary where the key is a node x and the value is the parent of node x (e.g. p['250.00']='250.0')
ddepth: the depth of each node in the icd9 hierarchy
highestLevel: length of the largest path in the tree

"""

def importOntoNum(ontonumfn,p,ddepth,highestLevel):
        f=open(ontonumfn,"r")
        f1=f.readlines()
        num={}
        num['root']={}
        num['root']['all']=0 #the parent is root
        num['root']['dead']=0 #the parent is root
        for x in f1:
                if len(x)<=1 or x[0]=="#":
                        continue
                a=x.split(",")
                if a[0] not in num:
                        num[a[0]]={}
                num[a[0]]['all']=int(a[1])
                num[a[0]]['dead']=int(a[2])
        currlevel=highestLevel
        bdone=True
        while bdone:
                for x in ddepth.keys():
                        if ddepth[x]==currlevel:
                                if p[x] not in num:
                                        num[p[x]]={}
                                        num[p[x]]['all']=0
                                        num[p[x]]['dead']=0
                                if x not in num:
                                        num[x]={}
                                        num[x]['all']=0
                                        num[x]['dead']=0
                                #num[p[x]]['all']+=num[x]['all'] #add if your ontonumfn doesn't propogate up
                                #num[p[x]]['dead']+=num[x]['dead'] #add if your ontonumfn doesn't propogate up
                currlevel-=1
                if currlevel==0:
                        bdone=False
        return num

""" 

Returns the Mutual Information of x based on the probabilities in the dictionary pr with the base of the log baseMI

""" 

def MI(x,pr,baseMI=2):
        try:
                sum1=1.0*pr[x+'dead']*math.log( 1.0*pr[x+'dead']/(pr[x]*pr['dead']) , baseMI) 
        except:
                sum1=0
        try:
                sum2=1.0*pr['not'+x+'dead']*math.log( 1.0*pr['not'+x+'dead']/(pr['not'+x]*pr['dead']) , baseMI) 
        except:
                sum2=0
        try:
                sum3=1.0*pr[x+'notdead']*math.log( 1.0*pr[x+'notdead']/(pr[x]*pr['notdead']) , baseMI)
        except:
                sum3=0
        try:
                sum4=1.0*pr['not'+x+'notdead']*math.log( 1.0*pr['not'+x+'notdead']/(pr['not'+x]*pr['notdead']) , baseMI) 
        except:
                sum4=0
        totalsum=sum1+sum2+sum3+sum4
        return totalsum

""" 

Returns the sigmoid of x

""" 

def sigmoid(x):
        try:
                return 1.0/(1+math.exp(-x))
        except:
                return 0 

""" 

Returns the maximum statistical signifincace based on a ztest for the whole tree

num: dictionary of all icd9 codes with keys 'all' and 'dead'
num[x]['all']: signifying the amount of people who have icd9 code x
num[x]['dead']: signifying the amount of people dead at that icd9 code x
d: dictionary of all icd9 codes
p: dictionary where the key is a node x and the value is the parent of node x (e.g. p['250.00']='250.0')
ddepth: the depth of each node in the icd9 hierarchy
highestLevel: length of the largest path in the tree
toprint: boolean whether to print the output

"""

def ztest_tree(num,d,p,ddepth,highestLevel,toprint=False):
        pdict={}
        if toprint:
                print ("Node 1    Node 2    P-value tree")
        for x in sorted(ddepth.keys()):
                for y in sorted(ddepth.keys()):
                        if x==y or x=='root' or y=='root': 
                                continue
                        if not isDes(x,y,p,ddepth) and not isDes(y,x,p,ddepth):
                                continue #it is on the wrong branch
                        stat, pval = proportions_ztest([num[x]['dead'],num[y]['dead']], [num[x]['all'],num[y]['all']])
                        if num[x]['all']==0 or math.isnan(pval):
                                pval=1
                        else:
                                pval/=num[x]['all']
                        if toprint:
                                print (x,"       ",y,"       ",'{0:0.3f}'.format(pval))
                        if x in pdict:
                                pdict[x].append(pval)
                        else:
                                pdict[x]=[pval]
                        if y in pdict:
                                pdict[y].append(pval)
                        else:
                                pdict[y]=[pval]
        pdict2={} #maximum z-value for each pair
        for x in pdict:
                pdict2[x]=max(pdict[x])
                if math.isnan(pdict2[x]):
                        pdict2[x]=1
        return pdict2


""" 

Returns the maximum statistical signifincace based on a ztest for each node in branch

num: dictionary of all icd9 codes with keys 'all' and 'dead'
num[x]['all']: signifying the amount of people who have icd9 code x
num[x]['dead']: signifying the amount of people dead at that icd9 code x
d: dictionary of all icd9 codes
p: dictionary where the key is a node x and the value is the parent of node x (e.g. p['250.00']='250.0')
ddepth: the depth of each node in the icd9 hierarchy
highestLevel: length of the largest path in the tree
pathcomb: dictionary of all possible combinations
toprint: boolean whether to print the output

"""

def ztest_branch(num,d,p,ddepth,highestLevel,pathcomb,toprint=False):
        pdict={}
        if toprint:
                print ("Node 1    Node 2    P-value branch")
        for x2 in sorted(pathcomb.keys()):
                for x in sorted(pathcomb[x2]):
                        for y in sorted(pathcomb[x2]):
                                if x==y or x=='root' or y=='root':
                                        continue
                                if not isDes(x,y,p,ddepth) and not isDes(y,x,p,ddepth):
                                        continue #it is on the wrong branch
                                stat, pval = proportions_ztest([num[x]['dead'],num[y]['dead']], [num[x]['all'],num[y]['all']])
                                if num[x]['all']==0 or math.isnan(pval):
                                        pval=1
                                else:
                                        pval/=num[x]['all']
                                if toprint:
                                        print (x,"       ",y,"       ",'{0:0.3f}'.format(pval))
                                if (x2,x) in pdict:
                                        pdict[(x2,x)].append(pval)
                                else:
                                        pdict[(x2,x)]=[pval]
                                if (x2,y) in pdict:
                                        pdict[(x2,y)].append(pval)
                                else:
                                        pdict[(x2,y)]=[pval]
        pdict2={} #maximum z-value for each pair
        for x in pdict:
                pdict2[x]=max(pdict[x])
                if math.isnan(pdict2[x]):
                        pdict2[x]=1
        return pdict2

""" 

Calculates the HIC value based on the HIC formula

num: dictionary of all icd9 codes with keys 'all' and 'dead'
num[x]['all']: signifying the amount of people who have icd9 code x
num[x]['dead']: signifying the amount of people dead at that icd9 code x
d: dictionary of all icd9 codes
p: dictionary where the key is a node x and the value is the parent of node x (e.g. p['250.00']='250.0')
ddepth: the depth of each node in the icd9 hierarchy
weightbranch: weight of branch statistical significance as defined by equation 6
weighttree: weight of tree statistical significance as defined by equation 7
highestLevel: length of the largest path in the tree
toprint: boolean whether to print the output

"""

def calcHIC(num,d,p,ddepth,weightbranch,weighttree,highestLevel,toprint=False):
        pr={}
        pathcomb={}
        allleaf=getleaf('root',d,p,{})
        for qq in allleaf.keys():
                pathcomb[qq]=[]
        for q in ddepth.keys():
                for qq in allleaf.keys():
                        if isDes(q,qq,p,ddepth) or isDes(qq,q,p,ddepth):
                                if q!='root' and qq!='root':
                                        pathcomb[qq].append(q)
        pdict_tree=ztest_tree(num,d,p,ddepth,highestLevel,False)
        pdict_branch=ztest_branch(num,d,p,ddepth,highestLevel,pathcomb,False)
        if toprint:
                print ("BranchDepth(BDepth) depends on the branch,  s is the sigmoid function")
                print ("Node      Level BDepth Num       p(Dead)    s(MI) s(Log(lvl+max)) max(P-val branch) max(P-val tree) HIC")
        for x in sorted(ddepth.keys()):
                if x=='root':
                        continue
                currlevel=ddepth[x]
                if x not in pr:
                        pr[x]=1.0*num[x]['all']/num['root']['all']
                if 'not'+x not in pr:
                        pr['not'+x]=1.0-pr[x]
                if 'dead' not in pr:
                        pr['dead']=1.0*num['root']['dead']/num['root']['all']
                if 'notdead' not in pr:
                        pr['notdead']=1.0-pr['dead']
                if x+'dead' not in pr:
                        if num[x]['all']==0:
                                pr[x+'dead']=0 #div by 0
                        else:
                                pr[x+'dead']=1.0*num[x]['dead']/num[x]['all']
                if x+'notdead' not in pr:
                        if num[x]['all']==0:
                                pr[x+'notdead']=0 #div by 0
                        else:
                                pr[x+'notdead']=1.0*(num[x]['all']-num[x]['dead'])/num[x]['all']
                if 'not'+x+'dead' not in pr:
                        pr['not'+x+'dead']=0 #no division by 0
                        if num['root']['all'] != num[x]['all']:
                                pr['not'+x+'dead']=1.0*(num['root']['dead']-num[x]['dead'])/(num['root']['all']-num[x]['all'])
                if 'not'+x+'notdead' not in pr:
                        pr['not'+x+'notdead']=1-pr['not'+x+'dead'] #must sum to 1
                        if pr['not'+x+'notdead']<=0:
                                pr['not'+x+'notdead']=0
        tmphictree={}
        for x2 in sorted(pathcomb.keys()):
                for x in sorted(pathcomb[x2]):
                        maxlevel = len(pathcomb[x2])
                        if maxlevel == 1: #only one in the branch, skip
                                break
                        currlevel=ddepth[x]
                        HIC=sigmoid(MI(x,pr))-sigmoid(math.log(ddepth[x]+maxlevel,maxlevel))-weightbranch*pdict_branch[x]-weighttree*pdict_tree[x]+math.log(num[x]['all'],num['root']['all'])
                        if x not in tmphictree:
                                tmphictree[x]=[HIC]
                        else:
                                tmphictree[x].append(HIC)
                        if toprint:
                                if x==sorted(pathcomb[x2])[0]:
                                        print ("_"*82)
                                print ("{:6}".format(x),"{:6}".format(ddepth[x]),"{:6}".format(maxlevel),"{:6}".format(num[x]['all']),"{:10.3f}".format(pr[x+'dead']),"{:10.3f}".format(sigmoid(MI(x,pr))),"{:10.3f}".format(sigmoid(math.log(ddepth[x]+maxlevel,maxlevel))),"{:10.3f}".format(weightbranch*pdict_branch[x]),"{:10.3f}".format(weighttree*pdict_tree[x]),"{:10.3f}".format(HIC))
        hictree={}                        
        for icd9 in tmphictree:
                hictree[icd9]=max(tmphictree[icd9]) #get HIC based on maxlevel
        return hictree

"""

ontofn: path to ontology (icd9 hierarchy)
ontonumfn: path to file with the format icd9,amount of people with that icd9 code,amount of people with that icd9 code who have died
weightsfn: path to file with the branch and tree weight
picklefn: path to file which you want to save all the HIC values to

"""

if __name__ == '__main__':

        if len(sys.argv)==5:
                ontofn=sys.argv[1]
                ontonumfn=sys.argv[2]
                weightsfn=sys.argv[3]
                picklefn=sys.argv[4]
                dd,p,ddepth,highestLevel=importOnto(ontofn)
                num=importOntoNum(ontonumfn,p,ddepth,highestLevel)
                weightbranch, weighttree = importWeights(weightsfn)
                hictree=calcHIC(num,dd,p,ddepth,weightbranch,weighttree,highestLevel)
                picklefile = open(picklefn,"wb")
                savetree={}
                savetree['hictree']=hictree
                pickle.dump(savetree,picklefile)
                picklefile.close()
        else:
                print ("Usage: python HIC.py ontologyfn ontologynumfn weightsfn picklefn")