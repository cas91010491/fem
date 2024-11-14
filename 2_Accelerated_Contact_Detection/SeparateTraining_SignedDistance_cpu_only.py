import os
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
import pandas as pd
import tensorflow as tf
from tensorflow.keras import layers, Model
import matplotlib.pyplot as plt
import itertools

# Load data
data = pd.read_csv('sampled_data_1_percent.csv')
# Filter data based on the criterion -0.5 < gn < 1.5
filtered_data = data[(data['gn'] > -0.5) & (data['gn'] < 1.5)]
x_train = filtered_data[['x', 'y', 'z']].values
y_distance = filtered_data['gn'].values.reshape(-1, 1)

# Configuration parameters
layer_counts = [1, 2, 3]
neuron_options = [16, 32, 64]
epochs = 50

# Directory to save all models and results
output_dir = 'OUTPUT_signed_distance_models'
os.makedirs(output_dir, exist_ok=True)

# Initialize a list to collect performance metrics for each model
results = []

# Function to build, train, and evaluate each model configuration
def build_and_train_model(layer_config, x_train, y_distance, model_name):
    # Build the model
    inputs = tf.keras.Input(shape=(3,))
    x = inputs
    for neurons in layer_config:
        x = layers.Dense(neurons, activation='relu')(x)
    outputs = layers.Dense(1)(x)
    model = Model(inputs=inputs, outputs=outputs)
    
    # Compile the model with MSE loss and MAE metric
    model.compile(optimizer='adam', loss='mse', metrics=['mae'])

    # Prepare the model directory
    model_dir = os.path.join(output_dir, model_name)
    os.makedirs(model_dir, exist_ok=True)

    # Save model summary
    with open(os.path.join(model_dir, 'model_summary.txt'), 'w') as f:
        model.summary(print_fn=lambda x: f.write(x + '\n'))

    # Train the model
    history = model.fit(x_train, y_distance, epochs=epochs, batch_size=32, validation_split=0.2, verbose=0)

    # Save training history data
    history_df = pd.DataFrame(history.history)
    history_df.to_csv(os.path.join(model_dir, 'training_history.csv'), index=False)

    # Plot training and validation losses
    plt.figure()
    plt.plot(history.history['loss'], label='Training Loss')
    plt.plot(history.history['val_loss'], label='Validation Loss')
    plt.xlabel('Epochs')
    plt.ylabel('Loss (MSE)')
    plt.legend()
    plt.title(f'Training Loss for {model_name}')
    plt.savefig(os.path.join(model_dir, 'loss_plot.png'))
    plt.close()

    # Collect metrics for summary
    final_train_loss = history.history['loss'][-1]
    final_val_loss = history.history['val_loss'][-1]
    final_train_mae = history.history['mae'][-1]
    final_val_mae = history.history['val_mae'][-1]
    total_params = model.count_params()
    
    # Add results to the list
    results.append({
        'Model': model_name,
        'Layers': len(layer_config),
        'Neurons per Layer': layer_config,
        'Final Training Loss': final_train_loss,
        'Final Validation Loss': final_val_loss,
        'Final Training MAE': final_train_mae,
        'Final Validation MAE': final_val_mae,
        'Total Parameters': total_params
    })

# Loop through all configurations and train each model
for layer_count in layer_counts:
    for layer_config in itertools.product(neuron_options, repeat=layer_count):
        model_name = f"{layer_count}_layers_" + "_".join(map(str, layer_config)) + "_neurons"
        print(f"Training model: {model_name}")
        build_and_train_model(layer_config, x_train, y_distance, model_name)

# Save all results to a CSV file
results_df = pd.DataFrame(results)
results_df.to_csv(os.path.join(output_dir, 'results_summary.csv'), index=False)

print("Training completed for all configurations.")
print("Results summary saved to 'results_summary.csv'.")

# Load the results dataframe if it’s not already loaded
results_df = pd.read_csv(os.path.join(output_dir, 'results_summary.csv'))

# 1. **Minimum Validation Loss with Preference for Simpler Models**
min_val_loss = results_df['Final Validation Loss'].min()
optimal_models = results_df[results_df['Final Validation Loss'] <= 1.05 * min_val_loss]
best_model_by_val_loss = optimal_models.sort_values(by='Total Parameters').iloc[0]

print("Best Model by Validation Loss (within 5% of min loss):")
print(best_model_by_val_loss)

# 2. **Loss-to-Parameter Ratio for Efficiency**
results_df['Loss-to-Parameter Ratio'] = results_df['Final Validation Loss'] / results_df['Total Parameters']
best_model_by_ratio = results_df.loc[results_df['Loss-to-Parameter Ratio'].idxmin()]

print("\nBest Model by Loss-to-Parameter Ratio:")
print(best_model_by_ratio)

# 3. **Pareto Front for Balanced Trade-off**
# Sort by 'Total Parameters' and then by 'Final Validation Loss' to identify Pareto front models
sorted_df = results_df.sort_values(by=['Total Parameters', 'Final Validation Loss'])
pareto_front = sorted_df.drop_duplicates(subset='Total Parameters', keep='first')
best_model_pareto = pareto_front.iloc[0]

print("\nBest Model on Pareto Front (Balanced Trade-off):")
print(best_model_pareto)

# Save summary of best models to a text file
with open(os.path.join(output_dir, 'best_model_summary.txt'), 'w') as f:
    f.write("Best Model by Validation Loss (within 5% of min loss):\n")
    f.write(str(best_model_by_val_loss) + "\n\n")
    f.write("Best Model by Loss-to-Parameter Ratio:\n")
    f.write(str(best_model_by_ratio) + "\n\n")
    f.write("Best Model on Pareto Front (Balanced Trade-off):\n")
    f.write(str(best_model_pareto) + "\n")

print("\nSummary of best models saved to 'best_model_summary.txt'")

# Define acceptable error threshold
error_threshold = 0.01

# Function to calculate threshold accuracy (for MAE)
def calculate_threshold_accuracy(mae, threshold):
    return mae <= threshold

# Collect threshold accuracy for each model
threshold_accuracies = results_df['Final Validation MAE'].apply(lambda mae: calculate_threshold_accuracy(mae, error_threshold))
results_df['Threshold Accuracy'] = threshold_accuracies

# Select models that meet or exceed threshold accuracy
acceptable_models = results_df[results_df['Threshold Accuracy']]

print("\nModels with Validation MAE within the acceptable threshold (0.01):")
print(acceptable_models)

# Save threshold-acceptable models to a file
acceptable_models.to_csv(os.path.join(output_dir, 'acceptable_models.csv'), index=False)
