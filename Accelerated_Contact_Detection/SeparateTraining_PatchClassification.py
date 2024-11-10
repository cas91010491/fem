import os
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
y_classification = filtered_data['p_id'].values  # Patch IDs as integer labels (0 to 95)

# Configuration parameters
layer_counts = [1, 2, 3]
neuron_options = [64, 128, 256]
epochs = 50
num_classes = 96  # Number of patches for classification

# Directory to save all models and results
output_dir = 'OUTPUT_patch_classification_models'
os.makedirs(output_dir, exist_ok=True)

# Initialize a list to collect performance metrics for each model
results = []

def build_and_train_model(layer_config, x_train, y_classification, model_name):
    # Build the model
    inputs = tf.keras.Input(shape=(3,))
    x = inputs
    for neurons in layer_config:
        x = layers.Dense(neurons, activation='relu')(x)
    outputs = layers.Dense(num_classes, activation='softmax')(x)  # 96 units with softmax for classification
    model = Model(inputs=inputs, outputs=outputs)
    
    # Compile the model with sparse categorical crossentropy for classification and accuracy metric
    model.compile(optimizer='adam', loss='sparse_categorical_crossentropy', metrics=['accuracy'])

    # Prepare the model directory
    model_dir = os.path.join(output_dir, model_name)
    os.makedirs(model_dir, exist_ok=True)

    # Save model summary
    with open(os.path.join(model_dir, 'model_summary.txt'), 'w') as f:
        model.summary(print_fn=lambda x: f.write(x + '\n'))

    # Train the model
    history = model.fit(x_train, y_classification, epochs=epochs, batch_size=32, validation_split=0.2, verbose=0)

    # Save training history data
    history_df = pd.DataFrame(history.history)
    history_df.to_csv(os.path.join(model_dir, 'training_history.csv'), index=False)

    # Plot training and validation accuracy
    plt.figure()
    plt.plot(history.history['accuracy'], label='Training Accuracy')
    plt.plot(history.history['val_accuracy'], label='Validation Accuracy')
    plt.xlabel('Epochs')
    plt.ylabel('Accuracy')
    plt.legend()
    plt.title(f'Training and Validation Accuracy for {model_name}')
    plt.savefig(os.path.join(model_dir, 'accuracy_plot.png'))
    plt.close()

    # Collect metrics for summary
    final_train_loss = history.history['loss'][-1]
    final_val_loss = history.history['val_loss'][-1]
    final_train_acc = history.history['accuracy'][-1]
    final_val_acc = history.history['val_accuracy'][-1]
    total_params = model.count_params()
    
    # Add results to the list
    results.append({
        'Model': model_name,
        'Layers': len(layer_config),
        'Neurons per Layer': layer_config,
        'Final Training Loss': final_train_loss,
        'Final Validation Loss': final_val_loss,
        'Final Training Accuracy': final_train_acc,
        'Final Validation Accuracy': final_val_acc,
        'Total Parameters': total_params
    })

# Loop through all configurations and train each model
for layer_count in layer_counts:
    for layer_config in itertools.product(neuron_options, repeat=layer_count):
        model_name = f"{layer_count}_layers_" + "_".join(map(str, layer_config)) + "_neurons"
        print(f"Training model: {model_name}")
        build_and_train_model(layer_config, x_train, y_classification, model_name)

# Save all results to a CSV file
results_df = pd.DataFrame(results)
results_df.to_csv(os.path.join(output_dir, 'results_summary.csv'), index=False)

print("Training completed for all configurations.")
print("Results summary saved to 'results_summary.csv'.")

# Load the results dataframe if it’s not already loaded
results_df = pd.read_csv(os.path.join(output_dir, 'results_summary.csv'))

# 1. **Maximum Validation Accuracy with Preference for Simpler Models**
max_val_acc = results_df['Final Validation Accuracy'].max()
optimal_models = results_df[results_df['Final Validation Accuracy'] >= 0.95 * max_val_acc]
best_model_by_val_acc = optimal_models.sort_values(by='Total Parameters').iloc[0]

print("Best Model by Validation Accuracy (within 5% of max accuracy):")
print(best_model_by_val_acc)

# 2. **Accuracy-to-Parameter Ratio for Efficiency**
results_df['Accuracy-to-Parameter Ratio'] = results_df['Final Validation Accuracy'] / results_df['Total Parameters']
best_model_by_ratio = results_df.loc[results_df['Accuracy-to-Parameter Ratio'].idxmax()]

print("\nBest Model by Accuracy-to-Parameter Ratio:")
print(best_model_by_ratio)

# 3. **Pareto Front for Balanced Trade-off**
# Sort by 'Total Parameters' and then by 'Final Validation Accuracy' to identify Pareto front models
sorted_df = results_df.sort_values(by=['Total Parameters', 'Final Validation Accuracy'], ascending=[True, False])
pareto_front = sorted_df.drop_duplicates(subset='Total Parameters', keep='first')
best_model_pareto = pareto_front.iloc[0]

print("\nBest Model on Pareto Front (Balanced Trade-off):")
print(best_model_pareto)

# Save summary of best models to a text file
with open(os.path.join(output_dir, 'best_model_summary.txt'), 'w') as f:
    f.write("Best Model by Validation Accuracy (within 5% of max accuracy):\n")
    f.write(str(best_model_by_val_acc) + "\n\n")
    f.write("Best Model by Accuracy-to-Parameter Ratio:\n")
    f.write(str(best_model_by_ratio) + "\n\n")
    f.write("Best Model on Pareto Front (Balanced Trade-off):\n")
    f.write(str(best_model_pareto) + "\n")

print("\nSummary of best models saved to 'best_model_summary.txt'")
