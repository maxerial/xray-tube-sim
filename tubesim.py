#!/usr/bin/env python3
"""
X-ray Spectrum Simulator with Stacked Filters

This script simulates X-ray spectra using the SpekPy library, allowing for
the definition of stacked filters via an external JSON configuration file.
It generates plots comparing the spectra for different filter combinations,
using PCHIP interpolation for visual smoothness while preserving sharp features.

Usage:
    python simulate_spectra.py

Configuration:
    Requires a 'simulation_config.json' file in the same directory.
    See accompanying example JSON for the required structure.

Dependencies:
    - spekpy
    - numpy
    - scipy
    - matplotlib

Output:
    - PNG (600 DPI) and PDF plots saved in a 'simulations' subfolder.
    - Filenames are generated based on the 'output_basename' and parameters
      defined in the JSON configuration.
"""

import json
import os
import sys
import spekpy as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy.interpolate import PchipInterpolator

# --- Configuration Constants ---
CONFIG_FILENAME = 'simulation_config.json'
OUTPUT_SUBFOLDER = 'simulations'
NUM_INTERP_POINTS = 3000 # Number of points for smooth interpolation plots

# --- Utility Functions ---

def load_config(filename):
    """Loads and validates the JSON configuration file."""
    print(f"Loading configuration from '{filename}'...")
    try:
        with open(filename, 'r') as f:
            config = json.load(f)
    except FileNotFoundError:
        print(f"Error: Configuration file '{filename}' not found.", file=sys.stderr)
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON from '{filename}': {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred reading '{filename}': {e}", file=sys.stderr)
        sys.exit(1)

    # Basic Validation
    if 'simulation_parameters' not in config or 'filter_combinations' not in config:
        print("Error: JSON must contain 'simulation_parameters' and 'filter_combinations'.", file=sys.stderr)
        sys.exit(1)
    if not isinstance(config['filter_combinations'], list):
         print("Error: 'filter_combinations' must be a list.", file=sys.stderr)
         sys.exit(1)

    # Check essential parameters
    required_params = ['kvp', 'anode_angle', 'min_energy_kev', 'max_energy_kev', 'output_basename']
    for param in required_params:
        if param not in config['simulation_parameters']:
            print(f"Error: Missing parameter '{param}' in 'simulation_parameters'.", file=sys.stderr)
            sys.exit(1)

    print("Configuration loaded successfully.")
    return config['simulation_parameters'], config['filter_combinations']

def create_output_dir(dirname):
    """Creates the output directory if it doesn't exist."""
    try:
        os.makedirs(dirname, exist_ok=True)
        print(f"Output directory: '{dirname}'")
    except OSError as e:
        print(f"Error creating output directory '{dirname}': {e}", file=sys.stderr)
        sys.exit(1)

def simulate_spectrum(kvp, angle, filter_list):
    """Generates a SpekPy spectrum, applies filters, and returns unique energy/fluence arrays."""
    try:
        # Create base spectrum (assumes default W target)
        s_combo = sp.Spek(kvp=kvp, th=angle)

        # Apply filters sequentially
        if not filter_list:
            print("  (No filters applied)")
        else:
            for filter_step in filter_list:
                 # Basic validation of filter format
                if len(filter_step) == 2 and isinstance(filter_step[0], str) and isinstance(filter_step[1], (int, float)):
                    filt_material, filt_thickness = filter_step
                    if filt_thickness < 0:
                         print(f"  Warning: Skipping filter with negative thickness: {filter_step}")
                         continue
                    print(f"  - Applying Filter: {filt_thickness:.2f} mm {filt_material}")
                    s_combo.filter(filt_material, filt_thickness)
                else:
                    print(f"  Warning: Skipping invalid filter format: {filter_step}")

        # Get spectrum data
        k_orig, phi_k_orig = s_combo.get_spectrum(edges=True)

        # Pre-process: Remove duplicate energy values for interpolation compatibility
        unique_indices = np.unique(k_orig, return_index=True)[1]
        if len(unique_indices) < len(k_orig):
            # print(f"  Info: Removed {len(k_orig) - len(unique_indices)} duplicate energy point(s).")
            k_unique = k_orig[unique_indices]
            phi_k_unique = phi_k_orig[unique_indices]
        else:
            k_unique = k_orig
            phi_k_unique = phi_k_orig

        return k_unique, phi_k_unique

    except Exception as e:
        print(f"  Error during SpekPy simulation: {e}", file=sys.stderr)
        return None, None


def interpolate_spectrum(k_unique, phi_k_unique, k_min_plot, k_max_plot, num_points):
    """Performs PCHIP interpolation on the spectrum data."""
    if k_unique is None or phi_k_unique is None or len(k_unique) <= 1:
        # print("  Info: Not enough points for interpolation.")
        return None, None # Cannot interpolate

    try:
        # Create the PCHIP interpolation function
        interp_func = PchipInterpolator(k_unique, phi_k_unique)

        # Create a fine energy grid within the plot limits and data range
        k_min_fine = max(np.min(k_unique), k_min_plot)
        k_max_fine = min(np.max(k_unique), k_max_plot)

        if k_max_fine <= k_min_fine:
            # print("  Info: Data range too narrow for interpolation within plot limits.")
            return None, None # Cannot create valid linspace

        k_fine = np.linspace(k_min_fine, k_max_fine, num_points)

        # Evaluate the interpolator
        phi_k_fine = interp_func(k_fine)

        # Ensure non-negativity (important for log plots)
        phi_k_fine = np.maximum(phi_k_fine, 1e-10)

        return k_fine, phi_k_fine

    except Exception as e:
        print(f"  Warning: Interpolation failed: {e}. Skipping interpolation.", file=sys.stderr)
        return None, None


def plot_and_save_spectra(results, params, output_basename, output_folder):
    """Generates, formats, and saves the comparison plot."""
    print("\nGenerating plot...")
    mpl.rcParams.update({'font.size': 16})
    plt.figure(figsize=(14, 9))

    plotted_something = False
    original_data_for_limits = []

    # Plot each result
    for res in results:
        label = res['label']
        k_orig, phi_orig = res['k_orig'], res['phi_orig']
        k_fine, phi_fine = res['k_fine'], res['phi_fine']

        # Store original data for accurate limit calculation
        if k_orig is not None and phi_orig is not None:
            original_data_for_limits.append({'k': k_orig, 'phi': phi_orig})

        # Plot interpolated data if available, otherwise plot original
        if k_fine is not None and phi_fine is not None:
            plt.plot(k_fine, phi_fine, label=label, linestyle='-', linewidth=1.5)
            plotted_something = True
        elif k_orig is not None and phi_orig is not None:
             print(f"  Info: Plotting original data for '{label}' (interpolation not available/failed).")
             plt.plot(k_orig, phi_orig, label=label, linestyle='-', linewidth=1.5)
             plotted_something = True
        else:
             print(f"  Warning: No valid data to plot for '{label}'.")


    if not plotted_something:
        print("Warning: No spectra were successfully plotted.", file=sys.stderr)
        plt.close() # Close empty figure
        return

    # --- Calculate Y Limits based on ORIGINAL data ---
    print("Calculating Y-axis limits...")
    global_y_min_in_range = np.inf
    global_y_max_in_range = 0
    first_spectrum = True
    min_k_plot = params['min_energy_kev']
    max_k_plot = params['max_energy_kev']

    for data in original_data_for_limits:
        k_data, phi_data = data['k'], data['phi']
        valid_indices = (k_data >= min_k_plot) & (k_data <= max_k_plot)
        if np.any(valid_indices):
             phi_in_range = phi_data[valid_indices]
             if phi_in_range.size > 0:
                 current_max = np.max(phi_in_range)
                 if first_spectrum:
                     global_y_max_in_range = current_max
                     first_spectrum = False
                 else:
                     global_y_max_in_range = max(global_y_max_in_range, current_max)

                 phi_positive_in_range = phi_in_range[phi_in_range > 1e-6]
                 if phi_positive_in_range.size > 0:
                      global_y_min_in_range = min(global_y_min_in_range, np.min(phi_positive_in_range))

    # Set Y limits
    if global_y_max_in_range > 0:
        reasonable_bottom = max(1e1, global_y_min_in_range * 0.5 if global_y_min_in_range != np.inf else 1e1)
        reasonable_top = global_y_max_in_range * 1.5
        plt.ylim(bottom=reasonable_bottom, top=reasonable_top)
        print(f"Y-limits set based on data range: [{reasonable_bottom:.2e}, {reasonable_top:.2e}]")
    else:
        print("Warning: Could not determine valid Y range from data. Using default Y limits.")

    # --- Final Plot Formatting ---
    print("Formatting final plot...")
    plt.xlabel('Energy [keV]')
    plt.ylabel('Differential fluence [photons cm$^{-2}$ keV$^{-1}$ per mAs at 100 cm]')
    title = f"{params['kvp']:.0f} kVp, {params['anode_angle']:.0f} deg Spectra - {output_basename}"
    plt.title(title)
    plt.yscale('log')
    plt.xlim(min_k_plot, max_k_plot)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)

    # Add legend if labels were generated
    handles, labels = plt.gca().get_legend_handles_labels()
    if handles:
        plt.legend(handles, labels, fontsize=10) # Adjust fontsize/ncol if needed
    else:
        print("Info: No labeled data to add to legend.")

    plt.tight_layout()

    # --- Save the plot (PNG and PDF) ---
    try:
        base_filename = f"{output_basename}_{params['kvp']:.0f}kv_{min_k_plot:.0f}-{max_k_plot:.0f}keV"
        png_filepath = os.path.join(output_folder, base_filename + ".png")
        pdf_filepath = os.path.join(output_folder, base_filename + ".pdf")

        plt.savefig(png_filepath, dpi=600)
        print(f"\nPlot saved as PNG: {png_filepath}")
        plt.savefig(pdf_filepath)
        print(f"Plot saved as PDF: {pdf_filepath}")
    except Exception as e:
        print(f"Error saving plot files: {e}", file=sys.stderr)

    # --- Show the plot ---
    plt.show()
    plt.close() # Close figure after showing/saving


# --- Main Execution ---
def main():
    """Main function to run the simulation and plotting workflow."""
    # Load configuration
    params, combinations = load_config(CONFIG_FILENAME)

    # Setup output
    create_output_dir(OUTPUT_SUBFOLDER)

    # Process each filter combination
    results_list = []
    print("\nProcessing filter combinations:")
    for combo in combinations:
        label = combo.get('label', 'Unnamed Combination') # Use .get for safety
        filter_list = combo.get('filters', [])

        print(f"\n- Processing: {label}")

        # Simulate spectrum
        k_unique, phi_k_unique = simulate_spectrum(
            params['kvp'], params['anode_angle'], filter_list
        )

        # Interpolate spectrum
        k_fine, phi_k_fine = interpolate_spectrum(
            k_unique, phi_k_unique,
            params['min_energy_kev'], params['max_energy_kev'],
            NUM_INTERP_POINTS
        )

        # Store results for plotting
        results_list.append({
            'label': label,
            'k_orig': k_unique,
            'phi_orig': phi_k_unique,
            'k_fine': k_fine,
            'phi_fine': phi_k_fine
        })

    # Generate plot
    if results_list:
        plot_and_save_spectra(
            results_list,
            params,
            params['output_basename'],
            OUTPUT_SUBFOLDER
        )
    else:
        print("\nNo results generated, skipping plot.", file=sys.stderr)

    print('\nFinished processing configuration!\n')


if __name__ == "__main__":
    main()