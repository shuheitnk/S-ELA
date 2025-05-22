import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score

# ---- Graph settings for vector export ----
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# ---- Fix random seed ----
random.seed(100)

# ---- Load data ----
name = "s_ela" # "medium_25d", "medium_50d", "large_25d", "large_50d", "s_ela"
df_y = pd.read_csv("wilcoxon_best_algo_1.csv")
df_x = pd.read_csv(name+"_features.csv")

feature_cols = df_x.columns[1:]
target_cols = df_y.columns[1:]

# ---- Function to extract valid instances (only one-best algorithm) ----
def extract_valid_instances(df_x, df_y, feature_names, target_names):
    X, Y, instances = [], [], []
    for instance in df_y["instance"]:
        y_row = df_y[df_y["instance"] == instance].loc[:, target_names].values[0]
        if sum(y_row) == 1:
            x_row = df_x[df_x["instance"] == instance].loc[:, feature_names].values[0]
            X.append(x_row)
            Y.append(y_row)
            instances.append(instance)
    return np.array(X), np.array(Y), instances

# ---- Prepare data for initial feature selection ----
X_tmp, Y_tmp, _ = extract_valid_instances(df_x, df_y, feature_cols, target_cols)
Y_tmp = np.array([np.argmax(y) for y in Y_tmp])  

# ---- Feature selection (currently using all features) ----
selected_features = list(set(feature_cols))  
print(f"Selected Features: {selected_features}")

# ---- Create final dataset ----
X, Y, _ = extract_valid_instances(df_x, df_y, selected_features, target_cols)
Y = np.array([np.argmax(y) for y in Y]) 

# ---- Model evaluation settings ----
accuracy_scores = []
seeds = np.random.randint(1, 10000, size=31)

# ---- Cross-validation and training Random Forest model ----
for i, seed in enumerate(seeds, 1):
    print(f"Run {i}")
    kf = StratifiedKFold(n_splits=10, shuffle=True, random_state=seed)
    fold_scores = []

    for train_idx, test_idx in kf.split(X, Y):
        train_x, test_x = X[train_idx], X[test_idx]
        train_y, test_y = Y[train_idx], Y[test_idx]

        model = RandomForestClassifier(bootstrap=True, random_state=seed, class_weight='balanced')
        model.fit(train_x, train_y)
        pred = model.predict(test_x)

        fold_scores.append(accuracy_score(test_y, pred))

    accuracy_scores.append(np.mean(fold_scores))

# ---- Visualize accuracy scores ----
plt.figure(figsize=(10, 6))
plt.boxplot(accuracy_scores, vert=True, patch_artist=True, boxprops=dict(facecolor="lightblue"))
plt.title("Accuracy Scores Across 31 Repeated 10-Fold Cross-Validation", fontsize=14)
plt.ylabel("Accuracy", fontsize=12)
plt.grid(axis="y", linestyle="--", alpha=0.7)
plt.show()

# ---- Save accuracy scores to CSV ----
output_df = pd.DataFrame({
    "run": list(range(1, 32)),
     name: accuracy_scores
})
output_df.to_csv("accuracy_score_"+name+".csv", index=False)
