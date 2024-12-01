import os
import sys
import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow.keras.models import load_model
from tensorflow.keras.utils import register_keras_serializable
from sklearn.metrics import classification_report, confusion_matrix
import seaborn as sns

import matplotlib.pyplot as plt

# Set the script's directory as the current working directory
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# Add the script's directory to the system path
sys.path.append(script_dir)


@register_keras_serializable()
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









# Load the Keras model
model_folder = 'OUTPUT_202412010019_multitask_DropOut0.0_50epochs_100percent_BatchNorm'
model_path = os.path.join(script_dir, model_folder+'/multitask_DropOut0.0_50epochs_100percent_BatchNorm.keras')
model = load_model(model_path)

# Load the data from the CSV file
data_path = os.path.join(script_dir, '..', 'sampled_data_1_percent.csv')
data = pd.read_csv(data_path)
data = data[(data['gn'] >= -0.5) & (data['gn'] <= 1.5)]

# Assuming 'data' is your DataFrame and 'X' is the input data for the model
# Select the first 3 columns if they are the relevant features
X = data.iloc[:, :3].values

# Now, X should have the shape (number_of_samples, 3)
predicted_distance, predicted_classification, predicted_projection = model.predict(X)

# Extract the predicted classes
y_true = data['p_id']
y_pred = np.argmax(predicted_classification, axis=1)

# Generate classification report
print(classification_report(y_true, y_pred))

# Compute confusion matrix
# Compute confusion matrix
cm = confusion_matrix(y_true, y_pred)
cm_log = np.log1p(cm)  # Apply logarithmic scale

classes = np.unique(y_true)
# Plot confusion matrix without numbers
plt.figure(figsize=(10, 8))
sns.heatmap(cm_log, annot=False, xticklabels=classes, yticklabels=classes, cmap='viridis')
plt.ylabel('Actual')
plt.xlabel('Predicted')
plt.show()



# cm = confusion_matrix(y_true, y_pred)

# # Plot confusion matrix
# plt.figure(figsize=(10, 8))
# sns.heatmap(cm, annot=True, fmt='.0f', xticklabels=classes, yticklabels=classes)
# plt.ylabel('Actual')
# plt.xlabel('Predicted')
# plt.show()
