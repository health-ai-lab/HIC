# Hierarchical Information Criterion for Variable Abstraction

## How to cite

M. Mirtchouk, B. Srikishan, and S. Kleinberg. Hierarchical Information Criterion for Variable Abstraction. In: *MLHC*, 2021

## Overview:

This code implements the HIC formula described in the paper above. This approach takes a hierarchy (e.g. ICD-9 codes), data about the hierarchy in regards to the outcome (e.g. for each ICD-9 code, the amount of people with that code, and how many of those had died) 
The primary application for the method is feature ranking based on hierarcical data.
We further assume an ontology is provided, which contains paths such as: root->250->250.0->250.00. In our ontology format:
root,250  \
250,250.0  \
250.0,250.00  \

## Using the code

### Preparation:

1. ontofn format: ensure ontofn is a single CSV file with each row being a parent child relationship (such as a ICD-9 codes) of the format: A,B where B is a child of A.

2. ontonumfn format: ensure ontonumfn matches the nodes presented in ontofn. Also, each row should heve 3 values (e.g. ICD-9 code, number of patients, number of patients who have died)

3. weightsfn format: ensure there are 2 rows: 1 for branch statistical significance and 1 for tree statistical significance which sum to 1. 
e.g.  \
branch,X  \
tree,Y  \

4. Prepare a directory to save the output file picklefn


### Usage with detailed parameters:

Usage: python HIC.py ontofn ontonumfn weightsfn picklefn

1. ontofn: path to ontology (icd9 hierarchy)
2. ontonumfn: path to file with the format icd9,amount of people with that icd9 code,amount of people with that icd9 code who have died
3. weightsfn: path to file with the branch and tree weight
4. picklefn: path to file which you want to save all the HIC values to

### To run the code:

python HIC.py /data/ontofn.csv /data/ontonumfn.csv /data/weightsfn.csv /data/picklefn.pkl

### File format examples

ontofn.csv example

A,B   \
A,C   \
A,D   \
B,E   \
B,F   \
F,G   \
C,H

Explanation: The parent is A who has 3 children B,C, and D. B has 2 children E and F. F has a child G. C has a child H. Therefore: D,E,G, and H are leaf nodes. 

ontonumfn.csv example  

A,1000,500 \
B,700,300 \
C,100,40 \
D,200,160 \
E,400,150 \
F,300,150 \
G,75,30 \
H,50,10 \

weightsfn.csv example \

branch,0.452 \
tree,0.548
