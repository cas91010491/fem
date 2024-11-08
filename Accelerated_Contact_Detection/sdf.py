import pandas as pd
import tensorflow as tf
from tensorflow.keras import layers, Model

# Load data
data = pd.read_csv('sampled_data_1_percent.csv')
x_train = data[['x', 'y', 'z']].values
y_classification = tf.keras.utils.to_categorical(data['p_id'].values, num_classes=96)

# Define test model for classification branch only
input_layer = layers.Input(shape=(3,))
shared_layer = layers.Dense(64, activation='relu')(input_layer)
classification_branch = layers.Dense(32, activation='relu')(shared_layer)
classification_output = layers.Dense(96, activation='softmax', name='classification_output')(classification_branch)

# Compile model with just classification output
classification_model = Model(inputs=input_layer, outputs=classification_output)
classification_model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])

# Test training to verify shape compatibility
history = classification_model.fit(x_train, y_classification, epochs=10, batch_size=32, validation_split=0.2)
