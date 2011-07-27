"""

generate inputs for the weighting-cdim code, run it, and return the results.
Note this is only for the "calcweights" code

"""

from __future__ import print_function

import os

import numpy
from numpy import where

import tempfile

import recfile

class WeightCalculator(dict):
    def __init__(self, data1, data2):
        self.data1=data1
        self.data2=data2

    def calc(self, n_near1, n_near2, cleanup=True):
        """

        Generate a set of weights such that the distribution of observables
        for dataset 1 are matched to dataset 2.  Use "nnear1" nearest
        neighbors on the first pass, nnear2 for the second.

        Dataset 1 is called the "training" set for photozs and dataset 2 is
        called the "photometric" sample. 

        The inputs for the calcweights code are wierd for the dataset 1, 
        consisting of 

            zspec extra weight n-dim-point

        For the weights calculation, the first three are not necessary at all.

        Inputs for dataset 2 is

            id n-dim-point

        Limitations:
            for now only support 1-d
        """

        
        data1=self.data1
        data2=self.data2
        
        ndim=len(data1.shape)
        if len(data1.shape) > 1 or len(data2.shape) > 1:
            raise ValueError("only support 1-d now")

        data1_dtype = [('junk1','f8'),('junk2','f8'),('junk3','f8'),
                       ('data','f8')]
        data2_dtype = [('id','i8'), ('data','f8')]

        input1 = numpy.zeros(data1.size, dtype=data1_dtype)
        input2 = numpy.zeros(data2.size, dtype=data2_dtype)

        input1['data'] = data1

        input2['id'] = numpy.arange(data2.size)
        input2['data'] = data2

        front='wts{ndim}'.format(ndim=ndim)
        f1=tempfile.mktemp(prefix=front+'-input-data1-', suffix='.dat')
        f2=tempfile.mktemp(prefix=front+'-input-data2-', suffix='.dat')

        wfile_pass1        = tempfile.mktemp(prefix=front+'-weights-pass1-', suffix='.dat')
        wfile_nozero_pass1 = tempfile.mktemp(prefix=front+'-weights-nozero-pass1-', suffix='.dat')
        numfile_pass1      = tempfile.mktemp(prefix=front+'-num-pass1-', suffix='.dat')

        wfile_pass2        = tempfile.mktemp(prefix=front+'-weights-pass2-', suffix='.dat')
        numfile_pass2      = tempfile.mktemp(prefix=front+'-num-pass2-', suffix='.dat')


        script_file = tempfile.mktemp(prefix=front+'-script-', suffix='.sh')

        try:
            print("writing temp file for data1:",f1)
            with recfile.Open(f1, mode='w', delim=' ') as fobj1:
                fobj1.write(input1)
            print("writing temp file for data2:",f2)
            with recfile.Open(f2, mode='w', delim=' ') as fobj2:
                fobj2.write(input2)

            script = _script.format(ndim=ndim,
                                    train_file           = f1,
                                    photo_file           = f2,
                                    n_near1              = n_near1,
                                    n_near2              = n_near2,
                                    weights_file1        = wfile_pass1,
                                    num_file1            = numfile_pass1,
                                    weights_file_nozero1 = wfile_nozero_pass1,
                                    weights_file2        = wfile_pass2,
                                    num_file2            = numfile_pass2)

            print("writing script file:",script_file)
            with open(script_file,'w') as sobj:
                sobj.write(script)

        finally:
            if cleanup:
                if os.path.exists(f1):
                    os.remove(f1)
                if os.path.exists(f2):
                    os.remove(f2)
                if os.path.exists(script_file):
                    os.remove(script_file)

_script = """#!/bin/bash
photo_file="{photo_file}"
train_file="{train_file}"

n_near1={n_near1}
n_near2={n_near2}

weights_file1="{weights_file1}"
weights_file_nozero1="{weights_file_nozero1}"
weights_file2="{weights_file2}"

num_file1="{num_file1}"
num_file2="{num_file2}"

# first run

echo "
Running calcweights

First run with n_near=$n_near1
"

calcweights{ndim}         \\
    "$train_file"    \\
    "$photo_file"    \\
    "$n_near1"       \\
    "$weights_file1" \\
    "$num_file1"

if [ "$?" != "0" ]; then
    echo Halting
    exit 45
fi


# remove training set objects with zero weight at first
# n_near
echo "
removing zero weight training objects to
    "$weights_file_nozero1"
"
awk '$3>0' < "$weights_file1" 1> "$weights_file_nozero1"

if [ "$?" != "0" ]; then
    echo Error running awk.  Halting
    exit 45
fi

# second run
echo "
Second run with n_near=$n_near2
"

calcweights{ndim}                \\
    "$weights_file_nozero1" \\
    "$photo_file"           \\
    "$n_near2"              \\
    "$weights_file2"        \\
    "$num_file2"

if [ "$?" != "0" ]; then
    echo Halting
    exit 45
fi
    """
