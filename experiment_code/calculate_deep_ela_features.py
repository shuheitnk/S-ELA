import pandas as pd
from pflacco.deep_ela import load_medium_25d_v1, load_medium_50d_v1, load_large_25d_v1, load_large_50d_v1

# Load datasets
decision_data = pd.read_csv('C:/decision_200.csv')

# please make fitness_200.csv from decision_200.csv in R-function in advance
fitness_data = pd.read_csv('C:/fitness_200.csv')

# Load Deep ELA models
medium_25d = load_medium_25d_v1() 
medium_50d = load_medium_50d_v1() 
large_25d = load_large_25d_v1()
large_50d = load_large_50d_v1()

# Extract unique instances and roops
instance_set = list(set(fitness_data["instance"]))
roop_set = list(set(decision_data["roop"]))

# Helper function to extract values from dictionary
def extract_values_from_dict(data_dict):
    """Extract values from dictionary and return as a list."""
    return list(data_dict.values())

# Initialize lists to store feature sets
features_medium_25d, features_medium_50d = [], []
features_large_25d, features_large_50d = [], []

# Iterate through instances and roops to process the data
for instance in instance_set:
    for roop in roop_set:
        print(f"Processing instance {instance}, roop {roop}")
        
        # Filter data based on the current instance and roop
        filtered_decision_data = decision_data[decision_data["roop"] == roop].iloc[:, 1:]
        filtered_fitness_data = fitness_data[(fitness_data["instance"] == instance) & (fitness_data["roop"] == roop)].iloc[:, 1:]
        
        # Convert filtered data to numpy arrays
        X = filtered_decision_data.to_numpy().T
        Y = filtered_fitness_data.to_numpy().T
        
        # Compute Deep ELA features for the current data
        feat_medium_25d = extract_values_from_dict(medium_25d(X, Y, include_costs=True))
        feat_medium_50d = extract_values_from_dict(medium_50d(X, Y, include_costs=True))
        feat_large_25d = extract_values_from_dict(large_25d(X, Y, include_costs=True))
        feat_large_50d = extract_values_from_dict(large_50d(X, Y, include_costs=True))
        
        # Append the features along with instance and roop information
        features_medium_25d.append([instance, roop] + feat_medium_25d)
        features_medium_50d.append([instance, roop] + feat_medium_50d)
        features_large_25d.append([instance, roop] + feat_large_25d)
        features_large_50d.append([instance, roop] + feat_large_50d)

# Helper function to flatten nested lists
def flatten_list(nested_list):
    """Flatten a nested list into a single list."""
    return [item for sublist in nested_list for item in sublist]

# Flatten feature lists
flattened_medium_25d = [flatten_list(f) for f in features_medium_25d]
flattened_medium_50d = [flatten_list(f) for f in features_medium_50d]
flattened_large_25d = [flatten_list(f) for f in features_large_25d]
flattened_large_50d = [flatten_list(f) for f in features_large_50d]
