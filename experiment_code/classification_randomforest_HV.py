import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier

# ---- Plot settings for vector graphic export ----
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# ---- Fix random seed ----
random.seed(100)

# ---- Configuration ----
name = "s_ela"  # Options: "medium_25d", "medium_50d", "large_25d", "large_50d", "s_ela"

# ---- Load datasets ----
df_labels = pd.read_csv("wilcoxon_best_algo_1.csv")
df_features = pd.read_csv(f"{name}_features.csv")
df_hv = pd.read_csv("bbob_biobj_hv_data.csv")

# ---- Extract column names ----
feature_columns = df_features.columns[1:]
target_columns = df_labels.columns[1:]
hv_columns = df_hv.columns[5:]

# ---- Extract instances with exactly one best algorithm ----
def extract_valid_instances(features_df, labels_df, hv_df, feature_names, label_names, hv_names):
    X, Y, HV, valid_instances = [], [], [], []
    for instance_id in labels_df["instance"]:
        label_row = labels_df[labels_df["instance"] == instance_id][label_names].values[0]
        if sum(label_row) == 1:  # only one-best algorithm
            feature_row = features_df[features_df["instance"] == instance_id][feature_names].values[0]
            hv_row = hv_df[hv_df["instance"] == instance_id][hv_names].values[0]
            X.append(feature_row)
            Y.append(label_row)
            HV.append(hv_row)
            valid_instances.append(instance_id)
    return np.array(X), np.array(Y), np.array(HV), valid_instances

# ---- Compute mean HV based on predicted class per test instance ----
def calculate_mean_hv_per_prediction(hv_array, predicted_classes):
    hv_values = {0: [], 1: [], 2: []}
    for i, pred_class in enumerate(predicted_classes):
        hv_values[pred_class].append(hv_array[i][pred_class])
    all_hvs = sum(hv_values.values(), [])
    return (
        np.mean(all_hvs),
        np.mean(hv_values[0]) if hv_values[0] else None,
        np.mean(hv_values[1]) if hv_values[1] else None,
        np.mean(hv_values[2]) if hv_values[2] else None
    )

# ---- Prepare input data ----
X_tmp, Y_tmp, HV_tmp, _ = extract_valid_instances(df_features, df_labels, df_hv, feature_columns, target_columns, hv_columns)
Y_tmp = np.argmax(Y_tmp, axis=1)

# ---- Feature selection (currently: use all) ----
selected_features = list(feature_columns)
print(f"Selected Features: {selected_features}")

# ---- Final input data preparation ----
X, Y, HV, _ = extract_valid_instances(df_features, df_labels, df_hv, selected_features, target_columns, hv_columns)
Y = np.argmax(Y, axis=1)

# ---- Settings for cross-validation ----
num_trials = 31
seeds = np.random.randint(1, 10000, size=num_trials)
hv_all, hv_sms_all, hv_nsga_all, hv_mole_all = [], [], [], []

# ---- Run 31 trials of stratified 10-fold cross-validation ----
# ---- Cross-validation loop ----
for i, seed in enumerate(seeds, 1):
    print(f"Run {i}")
    kf = StratifiedKFold(n_splits=10, shuffle=True, random_state=seed)
    hv_per_fold, sms_per_fold, nsga_per_fold, mole_per_fold = [], [], [], []

    for train_idx, test_idx in kf.split(X, Y):
        train_x, test_x = X[train_idx], X[test_idx]
        train_y, test_y = Y[train_idx], Y[test_idx]
        hv_test = HV[test_idx]

        model = RandomForestClassifier(bootstrap=True, random_state=seed, class_weight='balanced')
        model.fit(train_x, train_y)
        pred = model.predict(test_x)

        # Accuracy is calculated but not stored—can be added if needed
        hv, sms_hv, nsga_hv, mole_hv = calculate_mean_hv_per_prediction(hv_test, pred)
        hv_per_fold.append(hv)
        if sms_hv is not None:
           sms_per_fold.append(sms_hv)
        if nsga_hv is not None:
            nsga_per_fold.append(nsga_hv)
        if mole_hv is not None:
            mole_per_fold.append(mole_hv)


    hv_all.append(np.mean(hv_per_fold))
    hv_sms_all.append(np.mean(sms_per_fold))
    hv_nsga_all.append(np.mean(nsga_per_fold))
    hv_mole_all.append(np.mean(mole_per_fold))




# ---- Save accuracy scores to CSV ----
output_df = pd.DataFrame({
    "run": list(range(1, num_trials + 1)),
    "hv_all": hv_all, 
    "sms": hv_sms_all,
    "naga": hv_nsga_all, 
    "mole": hv_mole_all
})
output_df.to_csv(f"hv_values_{name}2.csv", index=False)

