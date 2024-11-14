# %%
import tensorflow as tf
from tensorflow.keras import layers, models, Model
from pdb import set_trace



def model_architecture(input_shape=(None,3)):
    # Define Model 1
    # Define the input layer
    inputs = layers.Input(shape=input_shape)

    # Shared Layer 1
    layer1 = layers.Dense(128, activation='relu')(inputs)

    # Branch 1 - Patch classification and (xi1, xi2) regression
    # Layer 2.1 for classification and segmented regression
    layer2_1 = layers.Dense(64, activation='relu')(layer1)

    # Classification output for patches (96 classes)
    classification_output = layers.Dense(96, activation='softmax', name="classification_output")(layer2_1)

    # Layer 3 for (xi1, xi2) regression
    layer3 = layers.Dense(128, activation='relu')(layer2_1)

    # Regression output for (xi1, xi2) with 192 outputs
    xi_output = layers.Dense(192, activation=None, name="xi_output")(layer3)

    # Branch 2 - Single output regression for gn
    # Layer 2.2 for gn regression
    layer2_2 = layers.Dense(64, activation='relu')(layer1)

    # Regression output for gn
    gn_output = layers.Dense(1, activation=None, name="gn_output")(layer2_2)

    # Combine branches into a model
    model = Model(inputs=inputs, outputs=[classification_output, xi_output, gn_output])

    # Display the model summary
    model.summary()
    return model



# # Custom loss for xi_output to calculate only for relevant patches
# @tf.function
# def xi_loss(y_true, y_pred):
#     # Extract one-hot encoded p_id, xi1, and xi2 from y_true
#     p_id_mask = y_true[:, :96]    # One-hot encoded patch ID mask
#     xi1_true = y_true[:, 96:192]  # True xi1 values for all patches
#     xi2_true = y_true[:, 192:288] # True xi2 values for all patches

#     # Split y_pred into xi1 and xi2 predictions
#     xi1_pred = y_pred[:, :96]
#     xi2_pred = y_pred[:, 96:192]

#     # Calculate masked loss
#     xi1_loss = tf.reduce_sum(p_id_mask * tf.square(xi1_true - xi1_pred), axis=1)
#     xi2_loss = tf.reduce_sum(p_id_mask * tf.square(xi2_true - xi2_pred), axis=1)

#     return tf.reduce_mean(xi1_loss + xi2_loss)



@tf.function
def xi_loss(y_true, y_pred):
    tf.print("xi_loss - y_true shape:", tf.shape(y_true))  # Debug print to verify shape
    tf.print("xi_loss - y_pred shape:", tf.shape(y_pred))  # Debug print to verify shape

    # Ensure y_true has 288 columns, which includes both xi1 and xi2 values for all 96 patches
    if tf.shape(y_true)[1] != 289:
        raise ValueError("Expected y_true to have 288 columns for xi_loss")

    # Extract the components for xi1 and xi2 from y_true
    p_id_mask = y_true[:, :96]        # One-hot encoded patch ID mask (first 96 columns)
    xi1_true = y_true[:, 96:192]      # True xi1 values for all patches (next 96 columns)
    xi2_true = y_true[:, 192:288]     # True xi2 values for all patches (last 96 columns)

    # Now slice y_pred accordingly for xi1 and xi2 predictions
    xi1_pred = y_pred[:, :96]
    xi2_pred = y_pred[:, 96:192]

    # Calculate the loss using p_id_mask to focus only on relevant patches
    xi1_loss = tf.reduce_sum(p_id_mask * tf.square(xi1_true - xi1_pred), axis=1)
    xi2_loss = tf.reduce_sum(p_id_mask * tf.square(xi2_true - xi2_pred), axis=1)

    return tf.reduce_mean(xi1_loss + xi2_loss)





# %%
import pandas as pd
import numpy as np
import tensorflow as tf

# Load data
sampled_data = pd.read_csv('sampled_data_1_percent.csv')

# Split into features and labels
features = sampled_data[['x', 'y', 'z']].values
labels = sampled_data[['p_id', 'xi1', 'xi2', 'gn']].values


set_trace()



# Prepare classification labels (one-hot encoding of p_id)
p_id = labels[:, 0].astype(int)
y_class = tf.keras.utils.to_categorical(p_id, num_classes=96)  # Shape: (num_samples, 96)

# Prepare xi labels for segmented regression (96 patches * 2 for xi1, xi2)
y_xi = np.zeros((len(labels), 192))
for i, idx in enumerate(p_id):
    y_xi[i, idx] = labels[i, 1]      # xi1 at patch idx
    y_xi[i, 96 + idx] = labels[i, 2]  # xi2 at patch idx + 96

# Prepare gn labels
y_gn = labels[:, 3].reshape(-1, 1)  # Shape: (num_samples, 1)

# Instead of combining, use a dictionary for separate outputs in the dataset
dataset = tf.data.Dataset.from_tensor_slices(
    (features, {'classification_output': y_class, 'xi_output': y_xi, 'gn_output': y_gn})
)

# Split into training and testing datasets
train_size = int(0.8 * len(dataset))
train_data = dataset.take(train_size).batch(32).prefetch(tf.data.AUTOTUNE)
test_data = dataset.skip(train_size).batch(32).prefetch(tf.data.AUTOTUNE)


# %%
# Define loss functions for each output
losses = {
    'classification_output': 'categorical_crossentropy',
    'xi_output': xi_loss,
    'gn_output': 'mean_squared_error'
}
loss_weights = {
    'classification_output': 1.0,
    'xi_output': 1.5,
    'gn_output': 1.5
}

# Compile the model with the specified losses and loss weights
model.compile(optimizer='adam', loss=losses, loss_weights=loss_weights,
              metrics={'classification_output': 'accuracy', 'xi_output': 'mae', 'gn_output': 'mae'})

# Train and save model
train_and_save_model(model, train_data, test_data, model_name="Simple_Dense_Model", epochs=5)


# %%


# %%
train_and_save_model(model, train_data, test_data, model_name="Simple_Dense_Model", epochs=5)

# Plot loss and validation loss
import matplotlib.pyplot as plt
plt.plot(history.history['loss'], label='Training Loss')
plt.plot(history.history['val_loss'], label='Validation Loss')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.legend()
plt.show()

sample_comparisons(model, test_data, num_samples=10)

# %%



