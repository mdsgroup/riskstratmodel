# -*- coding: utf-8 -*-
"""
After preprocess_ver2_nsfilter.R
Copyright Choi H

"""

import pandas as pd
import numpy as np

from keras.models import Sequential
from keras.layers.core import Dense, Activation, Dropout, Flatten
from keras.optimizers import SGD, RMSprop
from keras.regularizers import l2, activity_l2
from keras.layers import Convolution1D
import theano.tensor as T

from lifelines.utils import concordance_index, k_fold_cross_validation
from lifelines import CoxPHFitter

from sklearn.model_selection import ShuffleSplit

from pandas import DataFrame

edat1=pd.read_csv('./ModuleGenes/GSEred.csv')
edat2=pd.read_csv('./ModuleGenes/GSEturquoise.csv')
Clindata=pd.read_csv('./ModuleGenes/GSEpheno.csv')

X=np.dstack([np.asarray(edat1),np.asarray(edat2)])

E=np.asarray(Clindata["status"]=="dead", dtype='float64')
Y=np.asarray(Clindata["time"])

#Loss Function
def negative_log_likelihood(E):
	def loss(y_true,y_pred):
		hazard_ratio = T.exp(y_pred)
		log_risk = T.log(T.extra_ops.cumsum(hazard_ratio))
		uncensored_likelihood = y_pred.T - log_risk
		censored_likelihood = uncensored_likelihood * E
		neg_likelihood = -T.sum(censored_likelihood)
		return neg_likelihood
	return loss
	
cv=ShuffleSplit(n_splits=5, test_size=0.2, random_state=0)
iter=1
cv_result_val=[]
cv_result_train=[]
for train_index, val_index in cv.split(X):
	X_train=X[train_index]
	X_val=X[val_index]
	Y_train=Y[train_index]
	Y_val=Y[val_index]
	E_train=E[train_index]
	E_val=E[val_index]

	#Standardize_Scaling... 
	mean_train=np.mean(X_train,axis=0)
	std_train=np.std(X_train,axis=0)
	X_train_=(X_train-mean_train)/std_train
	X_val_=(X_val-mean_train)/std_train

	#Sorting for NNL!
	sort_idx = np.argsort(Y_train)[::-1]
	X_train_=X_train_[sort_idx]
	Y_train=Y_train[sort_idx]
	E_train=E_train[sort_idx]

	#Keras model
	np.random.seed(1337)
	model = Sequential()
	model.add(Convolution1D(24, 10, border_mode='valid', input_shape=(10, 2)))
	model.add(Activation('relu'))
	model.add(Dropout(0.5))
	model.add(Flatten())
	model.add(Dense(24)) # shape= length, dimension  init='glorot_uniform'
	model.add(Activation('relu'))
	model.add(Dropout(0.5))
	model.add(Dense(24))
	model.add(Activation('relu'))
	model.add(Dropout(0.5))
	model.add(Dense(1,activation="linear")) ##W_regularizer=l2(0.01), activity_regularizer=activity_l2(0.01)
	#
	
	sgd = SGD(lr=1e-5, decay=0, momentum=0.9, nesterov=True)
	rmsprop=RMSprop(lr=1e-4, rho=0.9, epsilon=1e-8)
	model.compile(loss=negative_log_likelihood(E_train), optimizer=rmsprop)

	print('Training...', iter, 'CVset')
	iteration=50
	ci_train=[]
	ci_val=[]
	for i in range(iteration):
		model.fit(X_train_, Y_train, batch_size=np.size(X_train,axis=0), 
				nb_epoch=10, shuffle=False, verbose=0)  # Shuffle False --> Important!!
		hr_pred=model.predict(X_train_)
		hr_pred=np.exp(hr_pred)
		ci=concordance_index(Y_train,-hr_pred,E_train)
	
		hr_pred2=model.predict(X_val_)
		hr_pred2=np.exp(hr_pred2)
		ci2=concordance_index(Y_val,-hr_pred2,E_val)
		print 'Iteration number', i+1
		print 'Concordance Index for training dataset:', ci
		print 'Concordance Index for test dataset:', ci2
		ci_train.append(ci)
		ci_val.append(ci2)

	#plt.figure()
	#plt.plot(ci_train, label='Training')
	#plt.plot(ci_val, label='Validation')
	#plt.xlabel('Epoch')
	#plt.ylabel('Concordance Index')
	#plt.legend()
	#plt.show()
	
	cv_result_train.append(ci_train[-1])
	cv_result_val.append(ci_val[-1])	
	
	iter+=1

print '-'*60
print 'CV result of deep neural network model: ', np.mean(cv_result_val)


######################################
#Complete Model Using same Parameters#
######################################
X_train=X
Y_train=Y
E_train=E

mean_train=np.mean(X_train,axis=0)
std_train=np.std(X_train,axis=0)
X_train_=(X_train-mean_train)/std_train

sort_idx = np.argsort(Y_train)[::-1]
X_train_=X_train_[sort_idx]
Y_train=Y_train[sort_idx]
E_train=E_train[sort_idx]

#Keras model
np.random.seed(1337)
model = Sequential()
model.add(Convolution1D(24, 10, border_mode='valid', input_shape=(10, 2)))
model.add(Activation('relu'))
model.add(Dropout(0.5))
model.add(Flatten())
model.add(Dense(24)) # shape= length, dimension  init='glorot_uniform'
model.add(Activation('relu'))
model.add(Dropout(0.5))
model.add(Dense(24))
model.add(Activation('relu'))
model.add(Dropout(0.5))
model.add(Dense(1,activation="linear")) ##W_regularizer=l2(0.01), activity_regularizer=activity_l2(0.01)
#
sgd = SGD(lr=1e-5, decay=0, momentum=0.9, nesterov=True)
rmsprop=RMSprop(lr=1e-4, rho=0.9, epsilon=1e-8)
model.compile(loss=negative_log_likelihood(E_train), optimizer=rmsprop)
model.fit(X_train_, Y_train, batch_size=np.size(X_train,axis=0), 
				nb_epoch=500, shuffle=False, verbose=0)  # Shuffle False --> Important!!
hr_pred=model.predict(X_train_)
hr_pred=np.exp(hr_pred)
ci=concordance_index(Y_train,-hr_pred,E_train)
print 'Concordance Index for all dataset:', ci


###########
#Test set##
###########

edattest1=pd.read_csv('./ModuleGenes_validation/GSEred.csv')
edattest2=pd.read_csv('./ModuleGenes_validation/GSEturquoise.csv')
Clindatatest=pd.read_csv('./ModuleGenes_validation/GSEpheno.csv')

X_test=np.dstack([np.asarray(edattest1),np.asarray(edattest2)])

E_test=np.asarray(Clindatatest["status"]=="dead", dtype='float64')
Y_test=np.asarray(Clindatatest["time"])

X_test_=(X_test-mean_train)/std_train

hr_test=model.predict(X_test_)
#hr_test=np.exp(hr_test)
ci_test=concordance_index(Y_test,-hr_test,E_test)

print '-'*50
print 'Testset CI:', ci_test

hr_df=DataFrame(hr_test, columns=['HR'])
hr_df['Time']=Y_test
hr_df['event']=E_test

cf=CoxPHFitter()
cf.fit(hr_df,'Time',event_col='event')
cf.print_summary()

############
#Test set2##
############

edattest1=pd.read_csv('./ModuleGenes_validation2/GSEred.csv')
edattest2=pd.read_csv('./ModuleGenes_validation2/GSEturquoise.csv')
Clindatatest=pd.read_csv('./ModuleGenes_validation2/GSEpheno.csv')

X_test=np.dstack([np.asarray(edattest1),np.asarray(edattest2)])

E_test=np.asarray(Clindatatest["status"]=="dead", dtype='float64')
Y_test=np.asarray(Clindatatest["time"])

X_test_=(X_test-mean_train)/std_train

hr_test2=model.predict(X_test_)
#hr_test=np.exp(hr_test)
ci_test2=concordance_index(Y_test,-hr_test2,E_test)

print '-'*50
print 'Testset2 CI:', ci_test2

hr_df2=DataFrame(hr_test2, columns=['HR'])
hr_df2['Time']=Y_test
hr_df2['event']=E_test

cf=CoxPHFitter()
cf.fit(hr_df2,'Time',event_col='event')
cf.print_summary()
