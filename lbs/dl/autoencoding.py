import keras
from keras.models import Model, Input
from keras.layers import Dense
import numpy as np
def encoder_correlation_check(data_to_pca, n_fet=6, labels=nm, prev_best=2):
    for i in range(1, n_fet):
        inp = Input(shape=(data_to_pca.shape[1], data_to_pca.shape[2]))
        l1 = Dense(i, activation = 'tanh')(inp)
        l2 = Dense(data_to_pca.shape[2])(l1)
        autoencoder = Model(inp, l2)
        encoder = Model(inp, l1)
        autoencoder.compile(optimizer='adam', loss='mse', metrics=['accuracy'])
        autoencoder.fit(data_to_pca, data_to_pca,
                        epochs=100,
                        batch_size=256,
                        shuffle=False)
        decoded = autoencoder.predict(data_to_pca)
        orginal = data_to_pca.mean(axis=1)
        decoded_check = decoded.mean(axis=1)
        decoded_check = pd.DataFrame(decoded_check)
        decoded_check[n_fet] = labels
        wq.append(decoded_check)
        result = decoded_check[prev_best].corr(decoded_check[n_fet])
        result_dev = (result - orginal_corr) ** 2
        ts.append(result_dev)
        return ts, wq
