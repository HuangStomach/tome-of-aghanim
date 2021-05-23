import os
import sys
import random as rn
import numpy as np
import pandas as pd

import tensorflow as tf
import tensorflow.keras as keras
from tensorflow.keras import backend as K

# 固定种子，创造可复现结果
os.environ['PYTHONHASHSEED'] = 0
session_conf = tf.compat.v1.ConfigProto(
    intra_op_parallelism_threads=1,
    inter_op_parallelism_threads=1
)
sess = tf.compat.v1.Session(graph=tf.compat.v1.get_default_graph() ,config=session_conf)
tf.compat.v1.keras.backend.set_session(session_conf)

from tensorflow.keras import metrics, optimizers
from tensorflow.keras.layers import Input, Dense, Lambda, Layer, Activation, Dropout, BatchNormalization
from tensorflow.keras.layers import BatchNormalization
from tensorflow.keras.models import Model
from tensorflow.keras.callbacks import Callback

output_dir = './Output'
train_file_path = './Data/V15.CCLE.4VAE.ZS.tsv'
val_file_1_path = './Data/V15.TCGA.4VAE.ZS.tsv'
train_latent_file = 'CCLE_latent.tsv'
train_weight_file = 'CCLE_weight.tsv'

predict_file = 'PANCAN_prediction.tsv'
encoder_file = 'CCLE_encoder_onehidden_vae.hdf5'
decoder_file = 'CCLE_decoder_onehidden_vae.hdf5'
print("output_dir: {}".output_dir)

epsilon_std = 1.0

def sampling(args):
    #import tensorflow as tf
    z_mean, z_log_var = args
    epsilon = K.random_normal(shape=tf.shape(input=z_mean), mean=0.,stddev=epsilon_std)
    z = z_mean + K.exp(z_log_var / 2) * epsilon
    return z