#!/usr/bin/env python3
import sys
#import java.io.FileReader as FileReader
import weka.core.Instances as Instances
import weka.classifiers.functions.MultilayerPerceptron as MP
import weka.core.converters.ConverterUtils.DataSource;
import weka.core.Instances;

#first paramter is arff file that is to be used

# check commandline parameters
if (not (len(sys.argv) == 2)):
    print ("Usage: UsingJ48.py <ARFF-file>")
    sys.exit()

# load data file

data1 = DataSource.read(sys.argv[1]);

# set the class Index - the index of the dependent variable
data1.setClassIndex(data1.numAttributes() - 1)

# create the model
print ("Training Multilayer Perceptron...")
multilayerperceptron = MP()
multilayerperceptron.buildClassifier(data)

# print out the built model
print ("Generated model:\n")
print (multilayerperceptron)