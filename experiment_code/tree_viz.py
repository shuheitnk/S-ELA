import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
import dtreeviz
import matplotlib.pyplot as plt

# Set font type for plotting
plt.rcParams['pdf.fonttype'] = plt.rcParams['ps.fonttype'] = 42

# Load the datasets containing target and feature data
target_df = pd.read_csv("wilcoxon_best_algo_1.csv")
features_df = pd.read_csv("s_ela_features.csv")


feature_columns = features_df.columns[1:]  # Skip the first column (instance identifier)
target_columns = target_df.columns[1:]  # Skip the first column (instance identifier)

# Prepare the feature matrix X and target vector Y
X, Y = [], []
for instance in target_df["instance"]:
    # Filter data by instance
    instance_target = target_df[target_df["instance"] == instance]
    instance_features = features_df[features_df["instance"] == instance]
    
    # Extract target values (ensure that the target is a single class)
    target_values = instance_target.loc[:, target_columns].values[0]
    if sum(target_values) == 1:  # Ensure that only one class is active
        Y.append(target_values)
        X.append(instance_features.loc[:, feature_columns].values[0])

# Convert to NumPy arrays for compatibility with scikit-learn
X, Y = np.array(X), np.array(Y)

# Y_classified represents a single class index
Y_classified = np.argmax(Y, axis=1)

# Initialize and train a Random Forest classifier
random_seed = np.random.randint(1, 100000)  # Generate a random seed for testing
model = RandomForestClassifier(max_depth=3, bootstrap=False, random_state=random_seed)
model.fit(X, Y_classified)

# Select a random tree from the trained RandomForest model
random_tree = model.estimators_[np.random.randint(0, len(model.estimators_))]

# Visualize the selected tree using dtreeviz
tree_viz = dtreeviz.model(
    model=random_tree,  # The selected random tree from the forest
    X_train=X,
    y_train=Y_classified,
    feature_names=feature_columns,
    target_name="Algorithm",
    class_names=["SMS-EMOA", "NSGA-II", "MOLE"]
)

# Display the tree visualization
tree_viz.view(label_fontsize=20, fontname='serif').show()
