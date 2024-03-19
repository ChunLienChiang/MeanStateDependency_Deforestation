import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

def read_dataset() -> tuple[dict, np.ndarray]:

    df = pd.read_csv(
        os.path.join(
            '../',
            'output/',
            'Output_Data/'
            'MSE_Budget_Summary/',
            'MSE_Budget_Summary.1970-2005.MC.csv'
        )
    )

    data = {
        key: df[key].values
        for key in df.columns if key not in ['En', 'Year', 'Group']
    }

    return data, df['Group'].values

def plot_barchart_mse_budget(data: dict, group: np.ndarray) -> None:

    # Set the output path
    output_path = os.path.join(
        '../',
        'output/',
        'Output_Figure/',
        'Plot_Barchart.MSE_Budget/',
    )
    output_file = 'MSE_Budget.MC.png'

    # Create the output directory if it does not exist
    if not os.path.exists(output_path):

        os.makedirs(output_path)

    # Create the figure
    fig, ax = plt.subplots(1, 1, figsize=(6, 5), dpi=300)

    for idx_variable, variable in enumerate(data.keys()):

        # Select the data for the variable
        data_variable_group1 = data[variable][group == 0]
        data_variable_group2 = data[variable][group == 3]

        ax.boxplot(
            [
                data_variable_group1,
            ],
            positions=[
                idx_variable * 3 + 0,
            ],
            whis=[5, 95],
            showmeans=True,
            meanline=True,
            showfliers=False,
            widths=0.6,
            patch_artist=True,
            boxprops=dict(facecolor='C0'),
            whiskerprops=dict(color='C0'),
            capprops=dict(color='C0'),
            medianprops=dict(color='black', linestyle=None),
            meanprops=dict(color='orange', linestyle=None),
        )
        
        ax.boxplot(
            [
                data_variable_group2,
            ],
            positions=[
                idx_variable * 3 + 1,
            ],
            whis=[5, 95],
            showmeans=True,
            meanline=True,
            showfliers=False,
            widths=0.6,
            patch_artist=True,
            boxprops=dict(facecolor='C3'),
            whiskerprops=dict(color='C3'),
            capprops=dict(color='C3'),
            medianprops=dict(color='black', linestyle=None),
            meanprops=dict(color='orange', linestyle=None),
        )

        ax.vlines(
            idx_variable * 3 + 2,
            -100,
            100,
            color='grey',
            linewidth=0.3,
        )
    
    ax.hlines(
        0,
        -1,
        len(data.keys()) * 3,
        color='black',
    )
    
    ax.set_xlim([-1, len(data.keys()) * 3 - 1])
    ax.set_xticks(np.arange(0, len(data.keys()) * 3, 3))
    ax.set_xticklabels(
        [
            'vertical mse advection',
            'horizontal mse advection',
            'latent heat',
            'sensible heat',
            'net longwave',
            'net shortwave',
            'residual'
        ],
        rotation=45,
    )
    ax.set_ylabel(r'MSE Budget ($W/m^{2}$)')
    ax.set_ylim([-20, 20])
    ax.set_title('MSE Budget')

    plt.tight_layout()
    plt.savefig(
        os.path.join(
            output_path,
            output_file,
        )
    )

    return

def plot_barchart_radiation_budget(data: dict) -> None:

    # Set the output path
    output_path = os.path.join(
        '../',
        'output/',
        'Output_Figure/',
        'Plot_Barchart.MSE_Budget/',
    )
    output_file = 'Radiation_Budget.MC.png'

    # Create the output directory if it does not exist
    if not os.path.exists(output_path):

        os.makedirs(output_path)

    # Create the figure
    fig, ax = plt.subplots(1, 1, figsize=(6, 5), dpi=300)

    for idx_variable, variable in enumerate(data.keys()):

        # Select the data for the variable
        data_variable_group1 = data[variable][group == 0]
        data_variable_group2 = data[variable][group == 3]

        ax.boxplot(
            [
                data_variable_group1,
            ],
            positions=[
                idx_variable * 3 + 0,
            ],
            whis=[5, 95],
            showmeans=True,
            meanline=True,
            showfliers=False,
            widths=0.6,
            patch_artist=True,
            boxprops=dict(facecolor='C0'),
            whiskerprops=dict(color='C0'),
            capprops=dict(color='C0'),
            medianprops=dict(color='black', linestyle=None),
            meanprops=dict(color='orange', linestyle=None),
        )
        
        ax.boxplot(
            [
                data_variable_group2,
            ],
            positions=[
                idx_variable * 3 + 1,
            ],
            whis=[5, 95],
            showmeans=True,
            meanline=True,
            showfliers=False,
            widths=0.6,
            patch_artist=True,
            boxprops=dict(facecolor='C3'),
            whiskerprops=dict(color='C3'),
            capprops=dict(color='C3'),
            medianprops=dict(color='black', linestyle=None),
            meanprops=dict(color='orange', linestyle=None),
        )
    
        ax.vlines(
            idx_variable * 3 + 2,
            -100,
            100,
            color='grey',
            linewidth=0.3,
        )
    
    ax.hlines(
        0,
        -1,
        len(data.keys()) * 3,
        color='black',
    )
    
    ax.set_xlim([-1, len(data.keys()) * 3 - 1])
    ax.set_xticks(np.arange(0, len(data.keys()) * 3, 3))
    ax.set_xticklabels(
        [
            'toa sw downward',
            'toa sw upward',
            'surface sw downward',
            'surface sw upward',
            'toa lw upward',
            'toa lw downward',
            'surface lw upward',
            'surface lw downward',
        ],
        rotation=45,
    )
    ax.set_ylabel(r'Radiation Budget ($W/m^{2}$)')
    ax.set_ylim([-10, 10])
    ax.set_title('Radiation Budget')

    plt.tight_layout()
    plt.savefig(
        os.path.join(
            output_path,
            output_file,
        )
    )

    return

if (__name__ == '__main__'):

    # Read dataset
    data, group = read_dataset()


    # Plot barchart for MSE budget
    plot_barchart_mse_budget(
        {
            key: data[key]
            for key in [
                'vertical_mse_advection_ano',
                'horizontal_mse_advection_ano',
                'latent_heat_flux_ano',
                'sensible_heat_flux_ano',
                'net_longwave_ano',
                'net_shortwave_ano',
                'residual_ano'
            ]
        },
        group,
    )

    # Plot barchart for radiation budget
    plot_barchart_radiation_budget(
        {
            key: data[key]
            for key in [
                'toa_sw_downward_ano',
                'toa_sw_upward_ano',
                'surface_sw_downward_ano',
                'surface_sw_upward_ano',
                'toa_lw_upward_ano',
                'toa_lw_downward_ano',
                'surface_lw_upward_ano',
                'surface_lw_downward_ano',
            ]
        }
    )