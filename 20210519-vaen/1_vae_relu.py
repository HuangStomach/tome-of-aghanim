import os
import pandas as pd

import tensorflow as tf
from tensorflow.keras import backend as K
from tensorflow.python.keras import activations

# session_conf = tf.compat.v1.ConfigProto(
#     intra_op_parallelism_threads=1,
#     inter_op_parallelism_threads=1
# )
# sess = tf.compat.v1.Session(graph=tf.compat.v1.get_default_graph() ,config=session_conf)
# tf.compat.v1.keras.backend.set_session(session_conf)

from tensorflow.keras import metrics, optimizers, Input
from tensorflow.keras.layers import Dense, Lambda, Layer, Activation, Dropout, BatchNormalization
from tensorflow.keras import Model
from tensorflow.keras.callbacks import Callback

import tensorflow.python.util.deprecation as deprecation
deprecation._PRINT_DEPRECATION_WARNINGS = False
tf.get_logger().setLevel('ERROR')

output_dir = './Output/1/'
train_file_path = './Output/0/V15.CCLE.4VAE.RANK.tsv'
val_file_1_path = './Output/0/V15.TCGA.4VAE.RANK.tsv'
print("output_dir: {}".format(output_dir))

def sampling(args):
    z_mean, z_log_var = args
    epsilon = K.random_normal(shape = tf.shape(input = z_mean), mean = 0., stddev = epsilon_std)
    z = z_mean + K.exp(z_log_var / 2) * epsilon
    return z

# 自定义变分层
class CustomVariationalLayer(Layer):
    def __init__(self, **kwargs):
        self.is_placeholder = True
        super(CustomVariationalLayer, self).__init__(**kwargs)

    def vae_loss(self, x_input, x_decoded):
        reconstruction_loss = original_dim * metrics.binary_crossentropy(x_input, x_decoded)
        kl_loss = - 0.5 * K.sum(
            1 + z_log_var_encoded - K.square(z_mean_encoded) - K.exp(z_log_var_encoded)
            , axis=-1)
        return K.mean(reconstruction_loss + (K.get_value(beta) * kl_loss))

    def call(self, inputs):
        x = inputs[0]
        x_decoded = inputs[1]
        loss = self.vae_loss(x, x_decoded)
        self.add_loss(loss, inputs=inputs)
        return x

class WarmUpCallback(Callback):
    def __init__(self, beta, kappa):
        self.beta = beta
        self.kappa = kappa
    def on_epoch_end(self, epoch, logs={}):
        if K.get_value(self.beta) <= 1:
            K.set_value(self.beta, K.get_value(self.beta) + self.kappa)

rnaseq_df = pd.read_table(train_file_path, index_col = 0)
val_df_1 = pd.read_table(val_file_1_path, index_col = 0)

test_set_percent = 0.1
rnaseq_test_df = rnaseq_df.sample(frac=test_set_percent) # 使用10%作为测试样本
rnaseq_train_df = rnaseq_df.drop(rnaseq_test_df.index) # 使用90%作为训练样本

original_dim = rnaseq_df.shape[1] # 原始癌症种类
units = 100 # 隐藏层单元数量

batch_size = 100
epochs = 100
learning_rate = 0.0005 # 学习步长

epsilon_std = 1.0
beta = K.variable(0)
kappa = 1

def train(i):
    i = str(_i)
    #tf.random.set_seed(_i)

    train_latent_file = i + '.CCLE_latent.tsv'
    train_weight_file = i + '.CCLE_weight.tsv'
    predict_file = i + '.TCGA_latent.tsv'
    encoder_file = i + '.CCLE_encoder_onehidden_vae.hdf5'
    decoder_file = i + '.CCLE_decoder_onehidden_vae.hdf5'

    rnaseq_input = Input(shape=(original_dim, )) # 作为神经网络入口的输入层

    z_mean_dense_linear = Dense(units, kernel_initializer='glorot_uniform')(rnaseq_input) # 有规律的密集连接的NN层
    z_mean_dense_batchnorm = BatchNormalization()(z_mean_dense_linear)
    z_mean_encoded = Activation('relu')(z_mean_dense_batchnorm) # 激活函数

    z_log_var_dense_linear = Dense(units, kernel_initializer='glorot_uniform')(rnaseq_input) # 有规律的密集连接的NN层
    z_log_var_dense_batchnorm = BatchNormalization()(z_log_var_dense_linear)
    z_log_var_encoded = Activation('relu')(z_log_var_dense_batchnorm) # 激活函数


    z = Lambda(sampling, output_shape=(units, ))([z_mean_encoded, z_log_var_encoded])

    drop_layer = Dropout(rate = 0.2, noise_shape = None)(z) # 将Dropout应用于输入。
    decoder_to_reconstruct = Dense(original_dim, kernel_initializer='glorot_uniform', activation='relu')
    rnaseq_reconstruct = decoder_to_reconstruct(drop_layer)

    """ 
    # 训练VAE模型
    adam = optimizers.Adam(learning_rate=learning_rate)
    vae_layer = CustomVariationalLayer()([rnaseq_input, rnaseq_reconstruct])
    vae = Model(rnaseq_input, vae_layer)
    vae.compile(optimizer=adam, loss=None, loss_weights=[beta])

    hist = vae.fit(np.array(rnaseq_train_df),
        shuffle=True,
        epochs=epochs,
        verbose=0,
        batch_size=batch_size,
        validation_data=(np.array(rnaseq_test_df), None),
        callbacks=[WarmUpCallback(beta, kappa)])
        #TQDMNotebookCallback(leave_inner=True, leave_outer=True)])

    history_df = pd.DataFrame(hist.history)
    loss_log_file = os.path.join(output_dir + "sampling.K.loss.log.txt")
    history_df.to_csv(loss_log_file, sep='\t') 
    """

    encoder = Model(rnaseq_input, z_mean_encoded)
    encoded_rnaseq_df = encoder.predict_on_batch(rnaseq_df)
    encoded_rnaseq_df = pd.DataFrame(encoded_rnaseq_df, index=rnaseq_df.index)

    encoded_rnaseq_df.columns.name = 'sample_id'
    encoded_rnaseq_df.columns = encoded_rnaseq_df.columns + 1
    encoded_file = os.path.join(output_dir, train_latent_file)
    encoded_rnaseq_df.to_csv(encoded_file, sep='\t')

    decoder_input = Input(shape=(units, ))  # can generate from any sampled z vector
    _x_decoded_mean = decoder_to_reconstruct(decoder_input)
    decoder = Model(decoder_input, _x_decoded_mean)

    encoder_model_file = os.path.join(output_dir, encoder_file)
    encoder.save(encoder_model_file)
    decoder_model_file = os.path.join(output_dir, decoder_file)
    decoder.save(decoder_model_file)
    # sum_node_activity = encoded_rnaseq_df.sum(axis=0).sort_values(ascending=False)

    weights = []
    for layer in decoder.layers:
        weights.append(layer.get_weights())

    weight_layer_df = pd.DataFrame(weights[1][0], columns=rnaseq_df.columns, index=range(1, 101))
    weight_layer_df.index.name = 'encodings'

    weight_file = os.path.join(output_dir, train_weight_file)
    weight_layer_df.to_csv(weight_file, sep='\t')

    ###################################
    encoded_val_df = encoder.predict_on_batch(val_df_1)
    encoded_val_df = pd.DataFrame(encoded_val_df, index=val_df_1.index)

    encoded_val_df.columns.name = 'sample_id'
    encoded_val_df.columns = encoded_val_df.columns + 1
    encoded_file = os.path.join(output_dir, predict_file)
    encoded_val_df.to_csv(encoded_file, sep='\t')
    print(encoded_file)
    print(encoded_val_df.head(2))

for _i in range(92, 101):
    tf.function(train(_i))
