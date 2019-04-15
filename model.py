#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import sklearn.preprocessing
from sklearn.metrics import cohen_kappa_score,hamming_loss,jaccard_similarity_score,hinge_loss
from datetime import datetime

# N = 100 # ÿ�����е�������
# D = 12 # ά��
# K = 3 # ������
# X = np.zeros((N*K,D)) # ����input
# y = np.zeros(N*K, dtype='uint8') # ����ǩ
# for j in xrange(K):
#   ix = range(N*j,N*(j+1))
#   r = np.linspace(0.0,1,N) # radius
#   t = np.linspace(j*4,(j+1)*4,N) + np.random.randn(N)*0.2 # theta
#   X[ix] = np.c_[r*np.sin(t), r*np.cos(t)]
#   y[ix] = j
# # ���ӻ�һ�����ǵ�������
# plt.scatter(X[:, 0], X[:, 1], c=y, s=40, cmap=plt.cm.Spectral)
# plt.show()

# 数据预处理
data = np.load('feature.npz')
pos=data['arr_0']
feature = data['arr_1']
label = data['arr_2'].astype(int)
gnd = data['arr_3'].astype(int).reshape(data['arr_3'].shape[0])

x = feature[np.where(label >-1)[0]]
X = sklearn.preprocessing.scale(x)
y = label[label >-1]
D = X.shape[1]
print 'num_examples',X.shape[0]
K = 4  # 类别

# 初次训练
h = 100
W = 0.01 * np.random.randn(D, h)
b = np.zeros((1, h))
W2 = 0.01 * np.random.randn(h, K)
b2 = np.zeros((1, K))

step_size = 1e-0
reg = 1e-3

num_examples = X.shape[0]
for i in xrange(5000):

    hidden_layer = np.maximum(0, np.dot(X, W) + b)
    scores = np.dot(hidden_layer, W2) + b2

    exp_scores = np.exp(scores)
    probs = exp_scores / np.sum(exp_scores, axis=1, keepdims=True)  # [N x K]

    corect_logprobs = -np.log(probs[range(num_examples), y])
    data_loss = np.sum(corect_logprobs) / num_examples
    reg_loss = 0.5 * reg * np.sum(W * W) + 0.5 * reg * np.sum(W2 * W2)
    loss = data_loss + reg_loss
    if i % 1000 == 0:
        print "iteration %d: loss %f" % (i, loss)

    dscores = probs
    dscores[range(num_examples), y] -= 1
    dscores /= num_examples

    dW2 = np.dot(hidden_layer.T, dscores)
    db2 = np.sum(dscores, axis=0, keepdims=True)

    dhidden = np.dot(dscores, W2.T)

    dhidden[hidden_layer <= 0] = 0

    dW = np.dot(X.T, dhidden)
    db = np.sum(dhidden, axis=0, keepdims=True)

    dW2 += reg * W2
    dW += reg * W

    W += -step_size * dW
    b += -step_size * db
    W2 += -step_size * dW2
    b2 += -step_size * db2

hidden_layer = np.maximum(0, np.dot(X, W) + b)
scores = np.dot(hidden_layer, W2) + b2
predicted_class = np.argmax(scores, axis=1)
print 'first training accuracy: %.2f' % (np.mean(predicted_class == y))

# 标签传播后的大模型训练
x = feature
X = sklearn.preprocessing.scale(x)
hidden_layer = np.maximum(0, np.dot(X, W) + b)
scores = np.dot(hidden_layer, W2) + b2
predicted_class = np.argmax(scores, axis=1)
print 'first accuracy: %.2f' % (np.mean(predicted_class == gnd))
predicted_class[np.where(label==1)[0]]=1
predicted_class[np.where(label==2)[0]]=2
predicted_class[np.where(label==3)[0]]=3
y = predicted_class
D = X.shape[1]
K = 4  # 类别
print 'num_examples',X.shape[0]

h = 100
W = 0.01 * np.random.randn(D, h)
b = np.zeros((1, h))
W2 = 0.01 * np.random.randn(h, K)
b2 = np.zeros((1, K))

step_size = 1e-0
reg = 1e-3

num_examples = X.shape[0]
for i in xrange(5000):

    hidden_layer = np.maximum(0, np.dot(X, W) + b)
    scores = np.dot(hidden_layer, W2) + b2

    exp_scores = np.exp(scores)
    probs = exp_scores / np.sum(exp_scores, axis=1, keepdims=True)  # [N x K]

    corect_logprobs = -np.log(probs[range(num_examples), y])
    data_loss = np.sum(corect_logprobs) / num_examples
    reg_loss = 0.5 * reg * np.sum(W * W) + 0.5 * reg * np.sum(W2 * W2)
    loss = data_loss + reg_loss
    if i % 1000 == 0:
        print "iteration %d: loss %f" % (i, loss)

    dscores = probs
    dscores[range(num_examples), y] -= 1
    dscores /= num_examples

    dW2 = np.dot(hidden_layer.T, dscores)
    db2 = np.sum(dscores, axis=0, keepdims=True)

    dhidden = np.dot(dscores, W2.T)

    dhidden[hidden_layer <= 0] = 0

    dW = np.dot(X.T, dhidden)
    db = np.sum(dhidden, axis=0, keepdims=True)

    dW2 += reg * W2
    dW += reg * W

    W += -step_size * dW
    b += -step_size * db
    W2 += -step_size * dW2
    b2 += -step_size * db2

hidden_layer = np.maximum(0, np.dot(X, W) + b)
scores = np.dot(hidden_layer, W2) + b2
predicted_class = np.argmax(scores, axis=1)
print 'second training accuracy: %.2f' % (np.mean(predicted_class == y))
print 'final accuracy: %.2f' % (np.mean(predicted_class == gnd))


np.savez('result.npz',predicted_class,gnd)
kappa=cohen_kappa_score(gnd,predicted_class)
hamming=hamming_loss(gnd,predicted_class)
jaccard=jaccard_similarity_score(gnd,predicted_class)
#hinge=hinge_loss(gnd,predicted_class)
print 'cohen_kappa_score:',kappa
print 'hamming_loss:',hamming
print 'jaccard_similarity_score:',jaccard
#print 'hinge_loss:',hinge

np.savetxt('metrics.txt',[kappa,hamming,jaccard],fmt='%.3e')
np.savetxt('result.txt',np.c_[pos,feature,predicted_class],fmt='%s')

print "完成"
print  datetime.now().strftime('%Y-%m-%d %H:%M:%S')