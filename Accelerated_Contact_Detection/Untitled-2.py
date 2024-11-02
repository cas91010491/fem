# %%
import pandas as pd
import glob
import tensorflow as tf
from pdb import set_trace
set_trace()



# Define the path where your CSV files are located
csv_files_path = '../csv_files/*.csv'

# Uncomment this section if you need to recreate a random 1% sample
# data_frames = []
# for file in glob.glob(csv_files_path):
#     df = pd.read_csv(file, header=None)  # Load without headers
#     if df.shape[1] == 7:                 # Ensure it has exactly 7 columns
#         data_frames.append(df)           # Append if structure is correct
#     else:
#         print(f"Skipping file {file} due to unexpected number of columns: {df.shape[1]}")
# 
# # Concatenate all valid DataFrames
# all_data = pd.concat(data_frames, ignore_index=True)
# 
# # Assign column names
# all_data.columns = ['x', 'y', 'z', 'p_id', 'xi1', 'xi2', 'gn']
# 
# # Drop any rows with NaN values
# all_data.dropna(inplace=True)
# 
# # Randomly sample 1% of the data
# sampled_data = all_data.sample(frac=0.01, random_state=42)
# 
# # Save the sample
# sampled_data.to_csv('sampled_data_1_percent.csv', index=False)

# Load the pre-saved 1% sampled data
sampled_data = pd.read_csv('Accelerated\ Contact\ Detection/sampled_data_1_percent.csv')

# Split into features (x, y, z) and labels (p_id, gn, xi1, xi2)
features = sampled_data[['x', 'y', 'z']].values
labels = sampled_data[['p_id', 'xi1', 'xi2', 'gn']].values

# Convert to TensorFlow dataset with (features, labels) tuples
sampled_data = tf.data.Dataset.from_tensor_slices((features, labels))

# Define the split ratio (80% train, 20% test)
train_size = int(0.8 * len(sampled_data))
train_data = sampled_data.take(train_size)
test_data = sampled_data.skip(train_size)

# Print dataset sizes to confirm
print(f"Total data size (1% sample): {len(sampled_data)}")
print(f"Training set size: {len(train_data)}")
print(f"Test set size: {len(test_data)}")


# %%
import numpy as np

# Adjust labels to match the output structure of the model (4 values per patch * 96 patches)
labels_expanded = np.zeros((len(labels), 4 * 96))

# Populate the expanded labels array for each sample
for idx, (p_id, xi1, xi2, gn) in enumerate(labels):
    # Assuming p_id is the patch index (0 to 95)
    patch_idx = int(p_id)
    labels_expanded[idx, patch_idx] = 1  # One-hot encoding for p_id
    labels_expanded[idx, 96 + patch_idx] = xi1
    labels_expanded[idx, 192 + patch_idx] = xi2
    labels_expanded[idx, 288 + patch_idx] = gn

# Update the TensorFlow dataset with expanded labels
sampled_data = tf.data.Dataset.from_tensor_slices((features, labels_expanded))

# Define the split ratio (80% train, 20% test)
train_size = int(0.8 * len(sampled_data))
train_data = sampled_data.take(train_size)
test_data = sampled_data.skip(train_size)




import tensorflow as tf

def custom_loss(y_true, y_pred):
    # Classification loss for the first 96 outputs (p_id prediction)
    p_id_true = y_true[:, :96]  # True one-hot encoded patch IDs
    p_id_pred = y_pred[:, :96]  # Predicted patch probabilities

    # Use categorical cross-entropy for patch classification
    classification_loss = tf.reduce_mean(
        tf.keras.losses.categorical_crossentropy(p_id_true, p_id_pred, from_logits=True)
    )


    # Regression loss for the xi1, xi2, and gn predictions
    xi1_true = y_true[:, 96:192]
    xi1_pred = y_pred[:, 96:192]
    xi2_true = y_true[:, 192:288]
    xi2_pred = y_pred[:, 192:288]
    gn_true = y_true[:, 288:]
    gn_pred = y_pred[:, 288:]

    # Identify the target patch for each sample
    target_patch = tf.argmax(p_id_true, axis=1)  # Indices of the target patches

    # Gather true and predicted values only for the target patch
    xi1_true_target = tf.gather(xi1_true, target_patch, axis=1, batch_dims=0)
    xi1_pred_target = tf.gather(xi1_pred, target_patch, axis=1, batch_dims=0)
    xi2_true_target = tf.gather(xi2_true, target_patch, axis=1, batch_dims=0)
    xi2_pred_target = tf.gather(xi2_pred, target_patch, axis=1, batch_dims=0)
    gn_true_target = tf.gather(gn_true, target_patch, axis=1, batch_dims=0)
    gn_pred_target = tf.gather(gn_pred, target_patch, axis=1, batch_dims=0)

    # Calculate MSE loss for xi1, xi2, and gn only for the target patch
    regression_loss = (
        tf.reduce_mean(tf.square(xi1_true_target - xi1_pred_target)) +
        tf.reduce_mean(tf.square(xi2_true_target - xi2_pred_target)) +
        tf.reduce_mean(tf.square(gn_true_target - gn_pred_target))
    )

    # Total loss: Combine classification and regression losses
    total_loss = classification_loss + regression_loss
    return total_loss




# %%
# Verify the shapes of features and labels
print("Features shape:", features.shape)
print("Labels shape:", labels.shape)

# Configure training and test data for batching and shuffling
batch_size = 32
train_data = train_data.batch(batch_size).shuffle(1000)
test_data = test_data.batch(batch_size)

# Display one batch to confirm
for batch_features, batch_labels in train_data.take(1):
    print("Batch features:", batch_features.numpy())
    print("Batch labels:", batch_labels.numpy())


# %%
from tensorflow.keras import layers, models

# Define the model architecture
model = models.Sequential()

# Input layer (assuming (x, y, z) as input)
model.add(layers.Input(shape=(3,)))

# Dense layers to increase dimensionality
model.add(layers.Dense(128, activation='relu'))
model.add(layers.Dense(512, activation='relu'))

# Reshape layer for 3D convolution
model.add(layers.Reshape((8, 8, 8, 1)))  # Adjust dimensions if necessary

# 3D convolutional layer to capture spatial dependencies
model.add(layers.Conv3D(96, kernel_size=(3, 3, 3), activation='relu', padding='same'))

# Pooling layer to reduce dimensionality
model.add(layers.MaxPooling3D(pool_size=(2, 2, 2)))

# Flatten the output from the convolutional layers
model.add(layers.Flatten())

# Dense layer for further feature processing
model.add(layers.Dense(256, activation='relu'))
model.add(layers.Dropout(0.3))  # Regularization to prevent overfitting

# Output layer with 384 outputs: (p_id, xi1, xi2, gn) for each patch
model.add(layers.Dense(4 * 96, activation=None))

# Compile the model with the custom loss function
model.compile(optimizer='adam', loss=custom_loss, metrics=['mae'])

# Display the model summary
model.summary()


# %%
# Set training parameters
epochs = 5

# Train the model
history = model.fit(train_data, 
                    validation_data=test_data, 
                    epochs=epochs)

# Optional: Plotting training history
import matplotlib.pyplot as plt

# Plot loss and validation loss
plt.plot(history.history['loss'], label='Training Loss')
plt.plot(history.history['val_loss'], label='Validation Loss')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.legend()
plt.show()


# %%
import numpy as np

# Number of samples to display for comparison
num_samples = 5

# Take a batch from the test set
for test_features, test_labels in test_data.take(1):
    # Make predictions on the batch
    predictions = model.predict(test_features)

    # Randomly select a few samples for comparison
    indices = np.random.choice(range(len(test_features)), num_samples, replace=False)
    print("Sample Comparisons (Predicted vs True):\n")
    
    for i in indices:
        print(f"Sample {i + 1}:")
        print("Predicted ->", predictions[i])
        print("True      ->", test_labels[i].numpy())
        print("-" * 30)

# %%
# Save the model in the .keras format
model_name = 'CLASSIF_REGRESS_128_512_CONV3D_256_20EPOCHS'

model.save(model_name+'.keras')

import json

# Save the training history
with open(model_name+'.json', 'w') as f:
    json.dump(history.history, f)


# To later load doing:
# from tensorflow.keras.models import load_model
# model = load_model('my_model.keras')

# with open('training_history.json', 'r') as f:
#     history_data = json.load(f)

# %%
model.name

# %%



