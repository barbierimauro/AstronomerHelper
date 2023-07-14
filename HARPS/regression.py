#!/scratch/mbarbier/miniconda3/bin/python3
import sys
import pandas as pd
import numpy as np
import scipy.stats as st
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, mean_absolute_error
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec

# Load the CSV data into a pandas DataFrame
df = pd.read_csv('catalog_data_for_regression.csv')

ranges = {
    'RV': [-100, 100],
    'vbroad': [5, 100],
    #'vbroad': [0.5, 2], #as log10
    'BPRP': [-1.2, 4],
    'logg': [3, 5],
    'FeH': [-2.5, 0.5],
    'bpabs': [-5, 15],
    'BPmag': [2.3, 14],
    'Teff': [2800, 10000],
    'RUWE' : [0.59,66],
    'Dist' : [0,  100]
}

# 'spectra_mag_b','spectra_mag_r','spectra_color_b_r',
# 'spectra_max_b','spectra_min_b','spectra_sum_b','spectra_mean_b',
# 'spectra_median_b','spectra_stdev_b','spectra_skew_b','spectra_kurt_b',
# 'spectra_iq99_b','spectra_iq01_b','spectra_lr_m_b','spectra_ks_v_b','spectra_rt_z_b',
# 'spectra_max_r','spectra_min_r','spectra_sum_r','spectra_mean_r',
# 'spectra_median_r','spectra_stdev_r','spectra_skew_r','spectra_kurt_r',
# 'spectra_iq99_r','spectra_iq01_r','spectra_lr_m_r','spectra_ks_v_r','spectra_rt_z_r',
# 'RUWE','BPmag','BPRP','RV','vbroad','Teff','logg','FeH','Dist','bpabs'


############################################################

# Define the columns of interest
#columns = ['spectra_skew_b', 'spectra_kurt_b', 'spectra_skew_r', 'spectra_kurt_r', 'spectra_color_b_r']
#columns = ['spectra_mag_b','spectra_mag_r']
#columns = ['spectra_mag_b','spectra_mag_r', 'spectra_skew_b']
#columns = ['spectra_skew_b','spectra_skew_r','spectra_kurt_b','spectra_kurt_r']
#columns = ['spectra_skew_b','spectra_kurt_r', 'spectra_color_b_r']
#columns = ['spectra_skew_b','spectra_kurt_b']
#columns = ['spectra_skew_r','spectra_kurt_r']
#columns = ['spectra_skew_b']
#columns = ['spectra_skew_r']
columns = ['spectra_mag_b','spectra_mag_r','spectra_color_b_r', 'spectra_max_b','spectra_min_b','spectra_sum_b','spectra_mean_b', 'spectra_median_b','spectra_stdev_b','spectra_skew_b','spectra_kurt_b', 'spectra_iq99_b','spectra_iq01_b','spectra_lr_m_b','spectra_ks_v_b','spectra_rt_z_b', 'spectra_max_r','spectra_min_r','spectra_sum_r','spectra_mean_r', 'spectra_median_r','spectra_stdev_r','spectra_skew_r','spectra_kurt_r', 'spectra_iq99_r','spectra_iq01_r','spectra_lr_m_r','spectra_ks_v_r','spectra_rt_z_r']
#columns = ['spectra_max_b','spectra_min_b','spectra_stdev_b','spectra_iq99_b','spectra_iq01_b','spectra_lr_m_b','spectra_ks_v_b','spectra_rt_z_b']
#columns = ['spectra_max_r','spectra_min_r','spectra_stdev_r','spectra_iq99_r','spectra_iq01_r','spectra_lr_m_r','spectra_ks_v_r','spectra_rt_z_r']
#columns = ['spectra_ks_v_b','spectra_ks_v_r']

transformed_columns = columns

#df[transformed_columns] = df[transformed_columns].apply(np.arcsinh)

target_column = 'Dist'

#df[target_column] = df[target_column].apply(np.arcsinh)
#df[target_column] = df[target_column].apply(np.log10)

study = target_column
polydeg = 1

############################################################


# Create a DataFrame that only includes rows where none of the columns of interest are NaN
df_filtered = df.dropna(subset=columns + [target_column])

# Split the data into features (X) and target variable (y)
X = df_filtered[columns]
y = df_filtered[target_column]
############################################################


############################################################


poly = PolynomialFeatures(degree=polydeg, interaction_only=False, include_bias=False)
X_poly = poly.fit_transform(X)
model = LinearRegression()
model.fit(X_poly, y)
y_pred = model.predict(X_poly)

# Print out the model's coefficients
#print("Intercept:    ", model.intercept_)
#print("Coefficients: ", model.coef_)
coef_df = pd.DataFrame({
    'feature': poly.get_feature_names_out(input_features=columns),
    'coefficient': model.coef_
})

print(coef_df)

# Calculate and print metrics
mse = mean_squared_error(y, y_pred)
rmse = np.sqrt(mse)
mae = mean_absolute_error(y, y_pred)
print("Mean Squared Error (MSE)      : ", mse)
print("Root Mean Squared Error (RMSE): ", rmse)
print("Mean Absolute Error (MAE)     : ", mae)




x_range_val = ranges[study]

residuals = y - y_pred
residuals_mean   = np.mean(residuals)
residuals_median = np.median(residuals)
residuals_std    = np.std(residuals)
residuals_sem    = st.sem(residuals)
residuals_skew   = st.skew(residuals)
residuals_kurt   = st.kurtosis(residuals)
nbins=500
residuals_hist,residuals_bins = np.histogram(residuals, bins=nbins*3)
residuals_hist = residuals_hist / residuals_hist.max()
residuals_bin_centers = 0.5*(residuals_bins[1:] + residuals_bins[:-1])  # Calculate bin centers
max_bin_index = np.argmax(residuals_hist)
residuals_mode = (residuals_bins[max_bin_index] + residuals_bins[max_bin_index + 1]) / 2

print("Residuals Mean                : ", residuals_mean  )
print("Residuals Median              : ", residuals_median)
print("Residuals Mode                : ", residuals_mode)
print("Residuals Standard Deviation  : ", residuals_std )
print("Residuals Standard Err Mean   : ", residuals_sem )
print("Residuals Skewness            : ", residuals_skew)
print("Residuals Kurtosis            : ", residuals_kurt)

ndev=5
y_range_val = [residuals_median-ndev*residuals_std,residuals_median+ndev*residuals_std]

# Create a gridspec figure
fig = plt.figure(figsize=(10, 8))
gs = gridspec.GridSpec(1, 4)


# Add a subplot for the rotated histogram of residuals (1/4 of the figure width)
subfig0 = plt.subplot(gs[0, :1])

residuals_pdf = st.norm.pdf(residuals_bin_centers, residuals_mean, residuals_std)
residuals_pdf = residuals_pdf / residuals_pdf.max()
subfig0.plot(residuals_hist, residuals_bin_centers, '-', drawstyle='steps-mid')  # Normalized histogram
subfig0.plot(residuals_pdf, residuals_bin_centers, '-')  # Normalized Gaussian
residuals_minus_gaussian = residuals_hist - residuals_pdf
subfig0.plot(residuals_minus_gaussian, residuals_bin_centers, '-', drawstyle='steps-mid')  # Histogram - Gaussian

subfig0.set_ylim(y_range_val)  # Use set_ylim for y-axis limit
subfig0.invert_xaxis()  # Invert the x-axis
subfig0.yaxis.set_ticks_position('left')  # Keep the y-axis labels on the left side
subfig0.set_ylabel('Residuals Values')
subfig0.grid(True)


# Add a subplot for the 2D histogram (3/4 of the figure width)
subfig1 = plt.subplot(gs[0, 1:])
subfig1.set_xlim(x_range_val)
subfig1.set_ylim(y_range_val)
subfig1.set_yticklabels([])

im = subfig1.hist2d(y, residuals, norm=LogNorm(), bins=(nbins,nbins), cmap=plt.cm.jet, range=[x_range_val, y_range_val])
subfig1.set_xlabel(target_column)
subfig1.grid(True)

fig.colorbar(im[3], ax=subfig1, label='Density')

plt.show()










sys.exit()
