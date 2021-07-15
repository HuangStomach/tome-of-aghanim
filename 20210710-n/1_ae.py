import os
import pandas as pd
# import numpy as np

import tensorflow as tf
from tensorflow.keras import metrics, optimizers, Input
from tensorflow.keras.layers import Dense, Lambda, Layer, Activation, Dropout, BatchNormalization
from tensorflow.keras import Model

import tensorflow.python.util.deprecation as deprecation
deprecation._PRINT_DEPRECATION_WARNINGS = False
tf.get_logger().setLevel('ERROR')

output_dir = './Output/1/'
train_ccle_file_path = './Output/0/V15.CCLE.4VAE.RANK.tsv'
train_tcga_file_path = './Output/0/V15.TCGA.4VAE.RANK.tsv'
print("output_dir: {}".format(output_dir))

rnaseq_df_ccle = pd.read_table(train_ccle_file_path, index_col = 0)
rnaseq_df_tcga = pd.read_table(train_tcga_file_path, index_col = 0)

original_dim = rnaseq_df_ccle.shape[1] # 原始癌症种类 6203
learning_rate = 0.00001
units = 100 # 目标维度
batch_size = 100
epochs = 100

def train(i):
    i = str(_i)

    train_ccle_latent_file = i + '.ccle_latent.tsv'
    # train_ccle_weight_file = i + '.CCLE_weight.tsv'
    train_tcga_latent_file = i + '.tcga_latent.tsv'
    # encoder_file = i + '.CCLE_encoder_onehidden_vae.hdf5'
    # decoder_file = i + '.CCLE_decoder_onehidden_vae.hdf5'

    input_layer = Input(shape=(original_dim,))
    encoded = Dense(units, activation='relu')(input_layer)
    encoded = BatchNormalization(name='bottleneck')(encoded)

    decoded = Dense(original_dim, activation='sigmoid')(encoded)
    
    autoencoder = Model(input_layer, decoded)
    adam = optimizers.Adam(learning_rate=learning_rate)
    autoencoder.compile(optimizer=adam, loss='binary_crossentropy')
    autoencoder.fit(
        rnaseq_df_ccle, rnaseq_df_ccle,
        epochs=epochs,
        batch_size=batch_size,
        shuffle=True,
        validation_split=.1,
        verbose=2
    )
    
    encoder_model = Model(autoencoder.input, autoencoder.get_layer(name='bottleneck').output)

    encoded_rnaseq_df_ccle = encoder_model.predict_on_batch(rnaseq_df_ccle)
    encoded_rnaseq_df_ccle = pd.DataFrame(encoded_rnaseq_df_ccle, index=rnaseq_df_ccle.index)

    encoded_rnaseq_df_ccle.columns.name = 'sample_id'
    encoded_rnaseq_df_ccle.columns = encoded_rnaseq_df_ccle.columns + 1
    encoded_file = os.path.join(output_dir, train_ccle_latent_file)
    encoded_rnaseq_df_ccle.to_csv(encoded_file, sep='\t')

    encoded_rnaseq_df_tcga = encoder_model.predict_on_batch(rnaseq_df_tcga)
    encoded_rnaseq_df_tcga = pd.DataFrame(encoded_rnaseq_df_tcga, index=rnaseq_df_tcga.index)

    encoded_rnaseq_df_tcga.columns.name = 'sample_id'
    encoded_rnaseq_df_tcga.columns = encoded_rnaseq_df_tcga.columns + 1
    encoded_file = os.path.join(output_dir, train_tcga_latent_file)
    encoded_rnaseq_df_tcga.to_csv(encoded_file, sep='\t')

    print(encoded_rnaseq_df_ccle.head(2))

for _i in range(1, 3):
    tf.function(train(_i))
