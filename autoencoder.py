#!/usr/bin/env python
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tensorflow as tf

from sklearn.datasets import make_classification
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input
from tensorflow.keras.layers import Dense
from tensorflow.keras.layers import LeakyReLU
from tensorflow.keras.layers import BatchNormalization
from tensorflow.keras.utils import plot_model
from matplotlib import pyplot

from sklearn.datasets import make_classification
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from tensorflow.keras.models import load_model
from sklearn.metrics import classification_report,confusion_matrix
# Random Forest
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import matthews_corrcoef


def autoencoder_generation(data):
	y=data["local_CA_nom"]
	data.pop("local_CA_nom")
	#data.pop("local_CA")
	data.pop("1")
	data.pop("2")
	data.pop("3")
	data.pop("4")
	data.pop("5")
	print(data.shape, y.shape)
	X_train, X_test, y_train, y_test = train_test_split(data, y, test_size=0.3, random_state=1)
	n_inputs = data.shape[1]
	print(data.shape)
	# scale data
	t = MinMaxScaler()
	t.fit(X_train)
	X_train = t.transform(X_train)
	X_test = t.transform(X_test)
	# define encoder
	visible = Input(shape=(n_inputs,))
	# encoder level 1
	e = Dense(n_inputs*2)(visible)
	e = BatchNormalization()(e)
	e = LeakyReLU()(e)
	# encoder level 2
	e = Dense(n_inputs)(e)
	e = BatchNormalization()(e)
	e = LeakyReLU()(e)
	# bottleneck
	n_bottleneck = round(float(n_inputs) / 5.0)
	bottleneck = Dense(n_bottleneck)(e)
	# define decoder, level 1
	d = Dense(n_inputs)(bottleneck)
	d = BatchNormalization()(d)
	d = LeakyReLU()(d)
	# decoder level 2
	d = Dense(n_inputs*2)(d)
	d = BatchNormalization()(d)
	d = LeakyReLU()(d)
	# output layer
	output = Dense(n_inputs, activation='linear')(d)
	# define autoencoder model
	model = Model(inputs=visible, outputs=output)
	# compile autoencoder model
	model.compile(optimizer='adam', loss='mse')
	# plot the autoencoder
	plot_model(model, 'autoencoder_no_compress.png', show_shapes=True)
	# fit the autoencoder model to reconstruct input
	history = model.fit(X_train, X_train, epochs=200, batch_size=16, verbose=2, validation_data=(X_test,X_test))
	# plot loss
	pyplot.plot(history.history['loss'], label='train')
	pyplot.plot(history.history['val_loss'], label='test')
	plt.title("Autoencoder Learning Curve", fontsize=14, fontweight='bold')
	plt.xlabel('Epoch')
	plt.ylabel('Loss')
	pyplot.legend()
	pyplot.show()
	# define an encoder model (without the decoder)
	encoder = Model(inputs=visible, outputs=bottleneck)
	plot_model(encoder, 'encoder_no_compress.png', show_shapes=True)
	# save the encoder to file
	encoder.save('encoder.h5')


def model(data):
	y=data["local_CA_nom"]
	data.pop("local_CA_nom")
	#data.pop("local_CA")
	data.pop("1")
	data.pop("2")
	data.pop("3")
	data.pop("4")
	data.pop("5")
	print(data.shape, y.shape)
	X_train, X_test, y_train, y_test = train_test_split(data, y, test_size=0.3, random_state=1)
	n_inputs = data.shape[1]
	print(data.shape)
	# scale data
	t = MinMaxScaler()
	t.fit(X_train)
	X_train = t.transform(X_train)
	X_test = t.transform(X_test)
	encoder = load_model('encoder.h5')
	# encode the train data
	X_train_encode = encoder.predict(X_train)
	# encode the test data
	X_test_encode = encoder.predict(X_test)
	# define the model
	model = RandomForestClassifier(n_estimators=200)
	model.fit(X_train_encode, y_train)
	RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini',max_depth=None, max_features='auto', max_leaf_nodes=None,
		min_impurity_decrease=0.0, min_impurity_split=None,
		min_samples_leaf=1, min_samples_split=2,
		min_weight_fraction_leaf=0.0, n_estimators=200, n_jobs=None,
		oob_score=False, random_state=None, verbose=0,
		warm_start=False)
	# fit the model on the training set
	model.fit(X_train_encode, y_train)
	rf_prediction = model.predict(X_test_encode)# Evaluations
	print('Classification Report: \n')
	MCC=matthews_corrcoef(y_test,rf_prediction, sample_weight=None )
	print("MCC",MCC)
	print(classification_report(y_test,rf_prediction))
	print('\nConfusion Matrix: \n')
	print(confusion_matrix(y_test,rf_prediction))
	scores = model.score(X_test_encode, y_test)
	print(scores)

	model.fit(X_train, y_train)
	rf_prediction = model.predict(X_test)# Evaluations
	print('Classification Report: \n')
	MCC=matthews_corrcoef(y_test,rf_prediction, sample_weight=None )
	print("MCC",MCC)
	print(classification_report(y_test,rf_prediction))
	print('\nConfusion Matrix: \n')
	print(confusion_matrix(y_test,rf_prediction))
	scores = model.score(X_test, y_test)
	print(scores)

data = pd.read_csv(sys.argv[1], header=0)
print(data.shape)
print(data.columns)
autoencoder_generation(data)
model(data)