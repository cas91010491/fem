import os
import pandas as pd

# Define the base directory and output file
base_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(base_dir)
output_file = os.path.join(base_dir, 'aggregated_results.csv')

# Initialize a list to hold the rows for the final CSV
data_rows = []

# Walk through the directories
for case in ['Ex1', 'Ex2', 'Ex3']:
    case_dir = os.path.join(base_dir, case)
    if not os.path.isdir(case_dir):
        continue

    for deformable in os.listdir(case_dir):
        deformable_dir = os.path.join(case_dir, deformable)
        if not os.path.isdir(deformable_dir):
            continue

        for experiment in os.listdir(deformable_dir):
            experiment_dir = os.path.join(deformable_dir, experiment)
            if not os.path.isdir(experiment_dir):
                continue

            # Parse the type from the folder name
            if "BFGS" in experiment:
                experiment_type = "BFGS-ANN" if "ANN" in experiment else "BFGS"
            elif "TR" in experiment:
                experiment_type = "TR-ANN" if "ANN" in experiment else "TR"
            else:
                continue

            # Paths to the required files
            counters_file = os.path.join(experiment_dir, 'COUNTERS.csv')
            timers_file = os.path.join(experiment_dir, 'TIMERS.csv')

            if not os.path.isfile(counters_file) or not os.path.isfile(timers_file):
                continue

            # Read the CSV files
            counters_data = pd.read_csv(counters_file)
            timers_data = pd.read_csv(timers_file)

            # Extract required data
            row = {
                'case': case,
                'deformable': deformable,
                'type': experiment_type,
                'Total': timers_data.iloc[0]['Total'],
                'Ktot+NRlinsolv': timers_data.iloc[0]['Ktot+NRlinsolv'],
                'Ktot+MINsolv': timers_data.iloc[0]['Ktot+MINsolv'],
                'fint': timers_data.iloc[0]['fint'],
                'fcon': timers_data.iloc[0]['fcon'],
                'ActSet_Updt': timers_data.iloc[0]['ActSet_Updt'],
                'NR_iters': counters_data.iloc[0]['NR_iters'],
                'mnmzn_iters': counters_data.iloc[0]['mnmzn_iters'],
                'mnmzn_fn_evals': counters_data.iloc[0]['mnmzn_fn_evals'],
            }

            data_rows.append(row)

# Create a DataFrame and save it to a CSV
output_df = pd.DataFrame(data_rows)


# import pdb; pdb.set_trace()

output_df.to_csv(output_file, index=False)

print(f"Aggregated data saved to {output_file}")
