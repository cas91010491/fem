import pandas as pd
import tensorflow as tf
import glob
from tensorflow.keras import layers, Model
import matplotlib.pyplot as plt
from tensorflow.keras import backend as K
import numpy as np
from pdb import set_trace
import os
from datetime import datetime
from tensorflow.keras.utils import plot_model
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Clear any lingering session state
K.clear_session()

# Define the path where your CSV files are located
csv_files_path = '../csv_files/*.csv'

# # Uncomment this section if you need to recreate a random 1% sample
# data_frames = []
# for file in glob.glob(csv_files_path):
#     df = pd.read_csv(file, header=None)  # Load without headers
#     if df.shape[1] == 7:                 # Ensure it has exactly 7 columns
#         data_frames.append(df)           # Append if structure is correct
#     else:
#         print(f"Skipping file {file} due to unexpected number of columns: {df.shape[1]}")

# # Concatenate all valid DataFrames
# all_data = pd.concat(data_frames, ignore_index=True)

# # Assign column names
# all_data.columns = ['x', 'y', 'z', 'p_id', 'xi1', 'xi2', 'gn']

# # Drop any rows with NaN values
# all_data.dropna(inplace=True)

# # Randomly sample 1% of the data
# sampled_data = all_data.sample(frac=1.0, random_state=42)

# # Save the sample
# sampled_data.to_csv('../sampled_data_100_percent.csv', index=False)




percent = 1
epochs = 200
Drop_factor = 0.2


model_name = 'multitask_DropOut'+str(Drop_factor)+'_'+str(epochs)+'epochs_'+str(percent)+'percent_BatchNorm'



# Load the data
data = pd.read_csv('../sampled_data_'+str(percent)+'_percent.csv')
# Filter rows where gn is within [-0.5, 1.5]
filtered_data = data[(data['gn'] >= -0.5) & (data['gn'] <= 1.5)]
n_data = filtered_data.shape[0]


set_trace()

# Get the current time
current_time = datetime.now().strftime('%Y%m%d%H%M')

# Create the output directory name
output_dir = f'OUTPUT_{current_time}_{model_name}'

# Create the directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)











# Separate features and labels
x_train = filtered_data[['x', 'y', 'z']].values
y_distance = filtered_data['gn'].values.reshape(-1, 1)
y_classification = filtered_data['p_id'].values.reshape(-1, 1)  # Integer labels for classification
y_projection = filtered_data[['xi1', 'xi2']].values


# Now `x_train`, `y_distance`, `y_classification`, and `y_projection` contain only filtered samples
y_projection_one_hot = np.zeros((n_data,2*96+1))        # the '+1' is to store the patch number to be used in the custom loss function
y_projection_one_hot[:,2*96] = y_classification.ravel()

for i in range(n_data):
    p_i = int(y_classification[i])
    y_projection_one_hot[i,[p_i,p_i+96]] = y_projection[i]


# Define the shared base model
input_layer = layers.Input(shape=(3,))
shared_layer = layers.Dense(64, activation='relu')(input_layer)


# Branch for Signed Distance Regression
distance_branch = layers.Dense(32, activation='relu')(input_layer)
concat_for_dist_layer = layers.Concatenate()([input_layer,distance_branch])
distance_layer = layers.Dense(64, activation='relu')(concat_for_dist_layer)
concat_for_dist_layer2 = layers.Concatenate()([distance_branch,distance_layer])
distance_layer2 = layers.Dense(32, activation='relu')(concat_for_dist_layer2)
signed_distance_output = layers.Dense(1, name='signed_distance_output')(distance_layer2)

# Shared Layer for Classification and Regression
shared_layer2 = layers.Dense(128, activation='relu')(shared_layer)
shared_layer2 = layers.BatchNormalization()(shared_layer2)

# Branch for Patch Classification
concat_for_class_branch = layers.Concatenate()([shared_layer2,input_layer])
classification_branch = layers.Dense(64, activation='relu')(concat_for_class_branch)
classification_branch = layers.Dropout(Drop_factor)(classification_branch)
classification_layer = layers.Dense(128, activation='relu')(classification_branch)
classification_output = layers.Dense(96, activation='softmax', name='classification_output')(classification_layer)

# Branch for Projection Regression
projection_branch = layers.Dense(256, activation='relu')(shared_layer2)
concat_for_proj_layer = layers.Concatenate()([shared_layer,projection_branch,classification_branch])
projection_layer = layers.Dense(128, activation='relu')(concat_for_proj_layer)
projection_layer = layers.Dropout(Drop_factor)(projection_layer)
concat_for_proj_layer2 = layers.Concatenate()([projection_layer, input_layer, classification_output])
projection_layer2 = layers.Dense(256, activation='relu')(concat_for_proj_layer2)
projection_output = layers.Dense(2*96+1, name='projection_output')(projection_layer)

# Build the full multitask model
model = Model(inputs=input_layer, outputs=[signed_distance_output, classification_output, projection_output])


def custom_loss(y_true, y_pred):

    # Extract patch/class labels and convert to one-hot to later reduce
    p_idx = y_true[:,-1]
    y_class = tf.keras.utils.to_categorical(p_idx, num_classes=96)

    # Regression targets for xi1, xi2, and gn
    xi1_true = y_true[:, 0:96]
    xi2_true = y_true[:, 96:192]
    xi1_pred = y_pred[:, 0:96]
    xi2_pred = y_pred[:, 96:192]

    # set_trace()

    # Weighted sum product using p_id_true as a mask to select the target patch
    xi1_loss = tf.reduce_sum(y_class * tf.square(xi1_true - xi1_pred), axis=1)
    xi2_loss = tf.reduce_sum(y_class * tf.square(xi2_true - xi2_pred), axis=1)

    # Total loss: Regression losses
    total_loss = tf.reduce_mean(xi1_loss + xi2_loss)
    return total_loss




# Compile the model with loss specified in order directly
model.compile(optimizer='adam',
            #   loss=['mse', 'sparse_categorical_crossentropy', 'mse'],
              loss=['mse', 'sparse_categorical_crossentropy', custom_loss],
              metrics={'classification_output': 'accuracy'})

# Verify that each output has the correct loss function
model.summary()

# Save the network architecture diagram
plot_model(model, to_file=os.path.join(output_dir, 'model_architecture.png'), show_shapes=True, show_layer_names=True)

# Save the model summary
with open(os.path.join(output_dir, 'model_summary.txt'), 'w') as f:
    model.summary(print_fn=lambda x: f.write(x + '\n'))




# Train the model
history = model.fit(
    x_train, 
    # [y_distance, y_classification, y_projection],
    [y_distance, y_classification, y_projection_one_hot],
    epochs=epochs, batch_size=32, validation_split=0.2
)



# Save the model in the .keras format
model.save(os.path.join(output_dir, model_name + '.keras'))

# Save the achieved accuracies and losses
history_df = pd.DataFrame(history.history)
history_df.to_csv(os.path.join(output_dir, 'training_history.csv'), index=False)





# Plot training and validation losses
plt.figure(figsize=(12, 4))
plt.subplot(1, 3, 1)
plt.plot(history.history['signed_distance_output_loss'], label='Train Signed Distance Loss')
plt.plot(history.history['val_signed_distance_output_loss'], label='Val Signed Distance Loss')
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.legend()
plt.title('Signed Distance Loss')

plt.subplot(1, 3, 2)
plt.plot(history.history['classification_output_loss'], label='Train Classification Loss')
plt.plot(history.history['val_classification_output_loss'], label='Val Classification Loss')
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.legend()
plt.title('Classification Loss')

plt.subplot(1, 3, 3)
plt.plot(history.history['projection_output_loss'], label='Train Projection Loss')
plt.plot(history.history['val_projection_output_loss'], label='Val Projection Loss')
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.legend()
plt.title('Projection Loss')

plt.tight_layout()

# Save the plot
plt.savefig(os.path.join(output_dir, 'training_losses.png'))

plt.show()



import numpy as np
import os

# Define the number of samples to check
num_samples_to_check = 20

# Randomly select samples from the dataset
sample_indices = np.random.choice(x_train.shape[0], num_samples_to_check, replace=False)
sample_data = x_train[sample_indices]
true_distance = y_distance[sample_indices]
true_classification = y_classification[sample_indices].flatten()  # Flatten to make it easier to match with predictions
true_projection = y_projection[sample_indices]

# Predict outputs
predicted_distance, predicted_classification, predicted_projection = model.predict(sample_data)




for i, idx in enumerate(sample_indices):
    # Get the true values for this sample
    true_p_id = true_classification[i]
    true_xi1, true_xi2 = true_projection[i]
    true_gn = true_distance[i][0]
    
    # Classification: get top 4 highest probabilities
    top4_indices = np.argsort(predicted_classification[i])[::-1][:4]
    top4_patches = [(patch, predicted_classification[i][patch]) for patch in top4_indices]
    
    # Predicted values for the top 2 patches for coordinates and signed distance
    top_patch_coords = [(patch, predicted_projection[i][patch], predicted_projection[i][96+patch]) for patch in top4_indices[:2]]
    top_patch_distances = [(patch, predicted_distance[i][0]) for patch in top4_indices[:2]]
    
    # Print the results
    print(f"Sample {idx}:")
    
    # True patch ID
    print(f"True p_id      -> {true_p_id}")
    
    # Predicted patch IDs with softmax values
    print("Predicted p_id -> ", end="")
    print(", ".join([f"{patch} ({value:.2e})" for patch, value in top4_patches]))
    
    # True (xi1, xi2)
    print(f"True (xi1, xi2)      ->           ({true_xi1:.2f}, {true_xi2:.2f})")
    
    # Predicted (xi1, xi2) for top 2 patches
    print("Predicted (xi1, xi2) -> ", end="")
    print(", ".join([f"patch {patch}: ({xi1:.2f}, {xi2:.2f})" for patch, xi1, xi2 in top_patch_coords]))
    
    # True signed distance
    print(f"True gn           -> {true_gn:.3f}")

    # Predicted signed distance for top 2 patches
    print(f"Predicted gn      -> {predicted_distance[i][0]:.3f}")
    print("-" * 30)
    

