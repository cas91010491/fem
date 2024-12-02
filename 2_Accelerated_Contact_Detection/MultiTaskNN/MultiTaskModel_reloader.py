import os
import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow.keras.models import load_model
from tensorflow.keras.callbacks import Callback
from sklearn.metrics import classification_report, confusion_matrix
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split

os.chdir(os.path.dirname(os.path.abspath(__file__)))




class SaveModelAndConfusionMatrix(Callback):
    def __init__(self,validation_data, save_dir, save_freq=50,initial_epoch=0):
        super(SaveModelAndConfusionMatrix, self).__init__()
        self.save_dir = save_dir
        self.save_freq = save_freq
        self.validation_data = validation_data
        self.initial_epoch = initial_epoch


    def on_train_begin(self, logs=None):
        
        os.makedirs(self.save_dir+"/ConfMtrx", exist_ok=True)

        if self.initial_epoch == 0:

            # Access the validation data
            X_val, y_true = self.validation_data
            y_true = y_true[1]  # Classification labels
            _, predicted_classification, _ = self.model.predict(X_val)
            y_pred = np.argmax(predicted_classification, axis=1)
            # Compute confusion matrix
            cm = confusion_matrix(y_true, y_pred)
            cm_log = np.log1p(cm)  # Apply logarithmic scale
            
            classes = np.unique(y_true)
            # Plot confusion matrix without numbers
            plt.figure(figsize=(10, 8))
            sns.heatmap(cm_log, annot=False, xticklabels=classes, yticklabels=classes, cmap='viridis')
            plt.ylabel('Actual')
            plt.xlabel('Predicted')
            colorbar = plt.gca().collections[0].colorbar
            colorbar.set_label(r'$\log(n_{samples})$')
            
            # Save confusion matrix plot

            cm_path = os.path.join(self.save_dir+"/ConfMtrx", f'confusion_matrix_epoch_0.png')
            plt.savefig(cm_path)
            plt.close()



    def on_epoch_end(self, epoch, logs=None):
        current_epoch = epoch + 1 + self.initial_epoch

        if current_epoch % self.save_freq == 0:
            # Save model
            model_path = os.path.join(self.save_dir, f'model_epoch_{current_epoch}.h5')
            self.model.save(model_path)
            
        if (current_epoch % 50 == 0) or (current_epoch in [1,2,5,10,20]):
            # Access the validation data
            X_val, y_true = self.validation_data
            y_true = y_true[1]  # Classification labels
            _, predicted_classification, _ = self.model.predict(X_val)
            y_pred = np.argmax(predicted_classification, axis=1)
            # Compute confusion matrix
            cm = confusion_matrix(y_true, y_pred)
            cm_log = np.log1p(cm)  # Apply logarithmic scale
            
            classes = np.unique(y_true)
            # Plot confusion matrix without numbers
            plt.figure(figsize=(10, 8))
            sns.heatmap(cm_log, annot=False, xticklabels=classes, yticklabels=classes, cmap='viridis')
            plt.ylabel('Actual')
            plt.xlabel('Predicted')
            
            # Save confusion matrix plot
            cm_path = os.path.join(self.save_dir+"/ConfMtrx", f'confusion_matrix_epoch_{current_epoch}.png')
            plt.savefig(cm_path)
            plt.close()

        # Save the training history after every epoch
        history_df = pd.DataFrame(self.model.history.history)
        history_df.to_csv(os.path.join(self.save_dir, 'training_history.csv'), index=False)


# Load the data
data = pd.read_csv('../sampled_data_50_percent.csv')
filtered_data = data[(data['gn'] >= -0.5) & (data['gn'] <= 1.5)]
n_data = filtered_data.shape[0]

# Separate features and labels
x_train = filtered_data[['x', 'y', 'z']].values
y_distance = filtered_data['gn'].values.reshape(-1, 1)
y_classification = filtered_data['p_id'].values.reshape(-1, 1)  # Integer labels for classification
y_projection = filtered_data[['xi1', 'xi2']].values

# Prepare y_projection_one_hot
y_projection_one_hot = np.zeros((n_data, 2*96+1))  # the '+1' is to store the patch number to be used in the custom loss function
y_projection_one_hot[:, 2*96] = y_classification.ravel()

for i in range(n_data):
    p_i = int(y_classification[i])
    y_projection_one_hot[i, [p_i, p_i+96]] = y_projection[i]

# Split the data into training and validation sets
x_train, x_test, y_distance_train, y_distance_test, y_classification_train, y_classification_test, y_projection_train, y_projection_test = train_test_split(
    x_train, y_distance, y_classification, y_projection_one_hot, test_size=0.2, random_state=42
)

# Pack outputs into tuples for Keras
y_train = (y_distance_train, y_classification_train, y_projection_train)
y_test = (y_distance_test, y_classification_test, y_projection_test)



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





# Load the trained model
model_dir = 'OUTPUT_202412011950_multitask_DropOut0.0_200epochs_50percent_BatchNorm_SUCCESS/'
model_name = 'multitask_DropOut0.0_200epochs_50percent_BatchNorm.keras'
model = load_model(model_dir+model_name, custom_objects={'custom_loss': custom_loss})

# Compile the model
opt = tf.keras.optimizers.Adadelta()
model.compile(optimizer=opt, loss=['mse', 'sparse_categorical_crossentropy', custom_loss], metrics={'classification_output': 'accuracy'})

# Create the callback
save_dir = model_dir[:-1]+'_ContinueTraining/'
save_callback = SaveModelAndConfusionMatrix(validation_data=(x_test, y_test), save_dir=save_dir, save_freq=50,initial_epoch=200)

# Continue training the model
history = model.fit(
    x_train, 
    [y_distance_train, y_classification_train, y_projection_train],
    epochs=200, batch_size=32, 
    validation_data=(x_test, [y_distance_test, y_classification_test, y_projection_test]),
    callbacks=[save_callback]
)

# Save the updated model in the .keras format
model.save(os.path.join(save_dir, 'updated_model_name.keras'))

# Save the achieved accuracies and losses
history_df = pd.DataFrame(history.history)
existing_history_df = pd.read_csv(os.path.join(model_dir, 'training_history.csv'))
updated_history_df = pd.concat([existing_history_df, history_df], ignore_index=True)
updated_history_df.to_csv(os.path.join(save_dir, 'updated_training_history.csv'), index=False)



