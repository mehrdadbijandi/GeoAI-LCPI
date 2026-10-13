import arcpy
import numpy as np
import math

def entropy_shannon_excel_style(norm_matrix):
    
    X = norm_matrix.astype(float)
    rows, cols = X.shape

    # Step 1: Column normalization (Excel-style)
    col_sums = X.sum(axis=0)
    P = X / col_sums

    # Step 2: Replace exact zeros for log(0) protection
    P_safe = np.where(P == 0, 1e-12, P)

    # Step 3: Excel-like precision
    P_excel = np.round(P_safe, 15)

    # Step 4: Compute P * ln(P) with rounding
    plnp = np.round(P_excel * np.log(P_excel), 15)

    # Step 5: Constant k
    k = 1.0 / np.log(rows)

    # Step 6: Entropy
    entropy = np.round(-k * plnp.sum(axis=0), 15)

    # Step 7: Redundancy
    d = np.round(1 - entropy, 15)

    # Step 8: Weights
    weights = np.round(d / d.sum(), 15)

    return weights


def main():
    # 1. Retrieve parameters from user input
    in_fc = arcpy.GetParameterAsText(0)           # Input Feature Class
    out_fc = arcpy.GetParameterAsText(1)          # Output Feature Class
    field_si = arcpy.GetParameterAsText(2)        # SI Field (Shape Irregularity)
    field_bc = arcpy.GetParameterAsText(3)        # BC Field (Boundary Complexity)
    field_tf = arcpy.GetParameterAsText(4)        # TF Field (Tenure Fragmentation)
    field_ac = arcpy.GetParameterAsText(5)        # AC Field (Accessibility)
    method = arcpy.GetParameterAsText(6)          # Weighting Method

    # Name the output score by the selected weighting method.
    lcpi_field_names = {
        "Equal Weights": "LCPI_EQ",
        "CRITIC": "LCPI_CR",
        "Shannon Entropy": "LCPI_SH",
        "PCA - 1st Component": "LCPI_PCA",
        "Standard Deviation": "LCPI_STD",
    }
    lcpi_field = lcpi_field_names[method]

    # List of input fields
    input_fields = [field_si, field_bc, field_tf, field_ac]
    
    arcpy.AddMessage(f"Starting calculation using method: {method}")
    
    # 2. Copy the input feature class to the output to preserve original data
    if arcpy.Exists(out_fc):
        arcpy.Delete_management(out_fc)
    arcpy.CopyFeatures_management(in_fc, out_fc)

    # 3. Read data into a NumPy Array with Float64 precision
    # Use 'handle null' logic to replace potential empty values with the column mean
    arr = arcpy.da.FeatureClassToNumPyArray(out_fc, input_fields)
    
    # Create raw matrix
    raw_matrix = np.zeros((len(arr), 4), dtype=np.float64)
    for i, f in enumerate(input_fields):
        col_data = arr[f]
        # Handle Null/NaN values: replace with the mean of the column
        if np.isnan(col_data).any():
            mean_val = np.nanmean(col_data)
            col_data = np.nan_to_num(col_data, nan=mean_val)
        raw_matrix[:, i] = col_data

    rows, cols = raw_matrix.shape

    # 4. Identical inverse min-max inputs for ALL five weighting methods.
    norm_matrix = np.zeros_like(raw_matrix)
    arcpy.AddMessage("Using common inverse Min-Max inputs for all methods")
    for j in range(cols):
        col_min = np.min(raw_matrix[:, j])
        col_max = np.max(raw_matrix[:, j])
        val_range = col_max - col_min
        if val_range == 0:
            norm_matrix[:, j] = 0.0
        else:
            norm_matrix[:, j] = (col_max - raw_matrix[:, j]) / val_range

    # 5. Calculate Weights based on the selected method
    weights = np.zeros(4) # [W_SI, W_BC, W_TF, W_AC]

    if method == "Equal Weights":
        weights[:] = 0.25

    elif method == "Standard Deviation":
        # Use ddof=1 for Sample Standard Deviation (matches Excel/SPSS behavior)
        stds = np.std(norm_matrix, axis=0, ddof=1)
        sum_std = np.sum(stds)
        
        if sum_std == 0:
            weights[:] = 0.25 # Fallback if all columns have zero variance
        else:
            weights = stds / sum_std

    elif method == "CRITIC":
        # Calculate Sample Standard Deviation
        stds = np.std(norm_matrix, axis=0, ddof=1)
        
        # Calculate Correlation Matrix
        # Note: If a column has zero variance, correlation results in NaN.
        with np.errstate(divide='ignore', invalid='ignore'):
            corr_matrix = np.corrcoef(norm_matrix, rowvar=False)
        
        # Replace NaN with 0 (handling constant columns)
        corr_matrix = np.nan_to_num(corr_matrix, nan=0.0)
        
        # Formula: C_j = sigma_j * sum(1 - r_jk)
        c_values = np.zeros(4)
        for j in range(4):
            sum_conflict = 0
            for k in range(4):
                sum_conflict += (1 - corr_matrix[j, k])
            c_values[j] = stds[j] * sum_conflict
            
        sum_c = np.sum(c_values)
        if sum_c == 0:
            weights[:] = 0.25
        else:
            weights = c_values / sum_c

    elif method == "Shannon Entropy":
        weights = entropy_shannon_excel_style(norm_matrix)

    elif method == "PCA - 1st Component":
        # ------------------------------------------------------------
        # STANDARD PCA WEIGHTING
        # ------------------------------------------------------------

        # Step 1: Correlation Matrix
        corr_matrix = np.corrcoef(norm_matrix, rowvar=False)

        # Step 2: Eigen Decomposition
        eigenvalues, eigenvectors = np.linalg.eig(corr_matrix)

        # Step 3: Select First Principal Component (largest eigenvalue)
        max_eig_idx = np.argmax(eigenvalues)
        pc1_vector = np.real(eigenvectors[:, max_eig_idx])

        # Step 4: Use Absolute Loadings
        abs_loadings = np.abs(pc1_vector)

        # Step 5: Normalize Loadings to Obtain Weights
        sum_loadings = np.sum(abs_loadings)

        if sum_loadings == 0:
            weights[:] = 0.25   # Fallback
        else:
            weights = abs_loadings / sum_loadings

    # 6. Display Calculated Weights
    arcpy.AddMessage("------------------------------------------------")
    arcpy.AddMessage(f"WEIGHTING METHOD: {method}")
    arcpy.AddMessage("------------------------------------------------")
    arcpy.AddMessage(f"Weight SI (Shape):  {weights[0]:.5f}")
    arcpy.AddMessage(f"Weight BC (Bound):  {weights[1]:.5f}")
    arcpy.AddMessage(f"Weight TF (Tenure): {weights[2]:.5f}")
    arcpy.AddMessage(f"Weight AC (Access): {weights[3]:.5f}")
    arcpy.AddMessage("------------------------------------------------")

    # 7. Calculate Final LCPI
    # LCPI = Sum(Weight_j * Normalized_Value_j)
    lcpi_values = np.dot(norm_matrix, weights)

    
    # 8. Save results to Output Feature Class
    # Define names for newly normalized fields and weight fields
    norm_field_names = ["SI", "BC", "TF", "AC"]
    weight_field_names = ["W_SI_Cal", "W_BC_Cal", "W_TF_Cal", "W_AC_Cal"]
    
    # Create fields for Normalized Values
    for fn in norm_field_names:
        arcpy.AddField_management(out_fc, fn, "DOUBLE")

    # Create fields for Weights
    for fn in weight_field_names:
        arcpy.AddField_management(out_fc, fn, "DOUBLE")
        
    arcpy.AddField_management(out_fc, lcpi_field, "DOUBLE")

    # Update list of fields to be written by Cursor
    fields_to_update = norm_field_names + weight_field_names + [lcpi_field]
    
    with arcpy.da.UpdateCursor(out_fc, fields_to_update) as cursor:
        i = 0
        for row in cursor:
            # Save Normalized Values
            row[0] = float(norm_matrix[i, 0])  # SI
            row[1] = float(norm_matrix[i, 1])  # BC
            row[2] = float(norm_matrix[i, 2])  # TF
            row[3] = float(norm_matrix[i, 3])  # AC
            
            # Save weights
            row[4] = float(weights[0])
            row[5] = float(weights[1])
            row[6] = float(weights[2])
            row[7] = float(weights[3])
            
            # Save calculated LCPI (Now guaranteed 0-1 for all methods)
            row[8] = float(lcpi_values[i])
            cursor.updateRow(row)
            i += 1

    arcpy.AddMessage(f"Calculation finished successfully. Output score field: {lcpi_field}")

if __name__ == "__main__":
    main()
