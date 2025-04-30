#!/usr/bin/env python3
"""
X-ray Spectrum Fixed Filter Thickness Comparison

This script simulates X-ray spectra using SpekPy to compare the effect
of fixed thicknesses (0.5, 1.0, 2.0 mm) for several common filter materials
(Sn, W, Au, Pb).

It generates a single plot showing all spectra, using color to group materials
and line style/opacity to distinguish thicknesses, highlighting the thickest
filter (2.0 mm) as solid/opaque.

Usage:
    python tubesim_filter_comparison.py

Configuration:
    Simulation parameters (kVp, angle, filters, thicknesses, plot range) are
    defined directly within the script constants.

Dependencies:
    - spekpy
    - numpy
    - matplotlib

Output:
    - PNG (600 DPI) and PDF plots saved in a 'simulations' subfolder.
    - Filename includes kVp and energy range settings.
"""

import os
import sys
import spekpy as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

# --- Configuration Constants ---
OUTPUT_SUBFOLDER = 'simulations'

# --- Simulation Parameters (Edit these directly) ---
KVP_SET = 450.0        # Set KVP [kV]
ANODE_ANGLE = 11.0     # Set anode angle [degrees]
MIN_K_PLOT = 100.0     # Min energy for plotting [keV]
MAX_K_PLOT = 450.0     # Max energy for plotting [keV]

# --- Filter Definitions (Edit these directly) ---
FILTERS_AND_THICKNESSES = {
    # Material: [List of thicknesses in mm]
    'Sn': [0.5, 1.0, 2.0],
    'W':  [0.5, 1.0, 2.0],
    'Au': [0.5, 1.0, 2.0],
    'Pb': [0.5, 1.0, 2.0]
}
# Thickness [mm] to plot with a solid line, full opacity (represents "best" hardening)
HIGHLIGHT_THICKNESS = 2.0

# --- Plot Styling Definitions ---
MATERIAL_COLORS = {
    'Sn': 'tab:blue',   # Blue for Tin
    'W':  'tab:orange', # Orange for Tungsten
    'Au': 'tab:green',  # Green for Gold
    'Pb': 'tab:red'     # Red for Lead
}
FIGURE_SIZE = (14, 9)   # Plot dimensions in inches
GLOBAL_FONT_SIZE = 16
LEGEND_FONT_SIZE = 10
LEGEND_COLUMNS = 2


def create_output_dir(dirname):
    """Creates the output directory if it doesn't exist."""
    try:
        os.makedirs(dirname, exist_ok=True)
        print(f"Output directory: '{dirname}'")
    except OSError as e:
        print(f"Error creating output directory '{dirname}': {e}", file=sys.stderr)
        sys.exit(1)

def main():
    """Runs the fixed filter thickness comparison simulation and plotting."""
    print('\nRunning script for fixed filter thickness comparison\n')

    # Setup output directory
    create_output_dir(OUTPUT_SUBFOLDER)

    # --- Plotting Setup ---
    mpl.rcParams.update({'font.size': GLOBAL_FONT_SIZE})
    plt.figure(figsize=FIGURE_SIZE)

    # --- Generate and Plot Unfiltered Spectrum ---
    print('Generating baseline unfiltered spectrum...')
    try:
        s_unfiltered = sp.Spek(kvp=KVP_SET, th=ANODE_ANGLE)
        k_unfiltered, phi_k_unfiltered = s_unfiltered.get_spectrum(edges=True)
        plt.plot(k_unfiltered, phi_k_unfiltered, label='Unfiltered',
                 linestyle='-', color='black', linewidth=1.5)
    except Exception as e:
        print(f"Error generating unfiltered spectrum: {e}", file=sys.stderr)
        # Continue execution if possible, but plotting might be affected
        k_unfiltered, phi_k_unfiltered = np.array([]), np.array([]) # Empty arrays

    # --- Variables for Y-limit calculation ---
    global_y_min_in_range = np.inf
    global_y_max_in_range = 0
    if k_unfiltered.size > 0: # Check if unfiltered spectrum was generated
        valid_indices_unfiltered = (k_unfiltered >= MIN_K_PLOT) & (k_unfiltered <= MAX_K_PLOT)
        if np.any(valid_indices_unfiltered):
            global_y_max_in_range = np.max(phi_k_unfiltered[valid_indices_unfiltered])
        else:
            global_y_max_in_range = 1e8 # Fallback
    else:
        global_y_max_in_range = 1e8 # Fallback

    # --- Generate and Plot Filtered Spectra ---
    print("\nGenerating filtered spectra:")
    plotted_filtered = False
    for material, thicknesses in FILTERS_AND_THICKNESSES.items():
        if material not in MATERIAL_COLORS:
            print(f"Warning: Color not defined for material '{material}'. Skipping.", file=sys.stderr)
            continue
        color = MATERIAL_COLORS[material]
        print(f"- Processing {material}")

        for thickness in thicknesses:
            print(f"  - Thickness: {thickness:.1f} mm")
            try:
                # Generate spectrum
                s_filtered = sp.Spek(kvp=KVP_SET, th=ANODE_ANGLE)
                s_filtered.filter(material, thickness)
                k_filtered, phi_k_filtered = s_filtered.get_spectrum(edges=True)

                # Determine line style, alpha, and linewidth based on thickness
                is_highlighted = np.isclose(thickness, HIGHLIGHT_THICKNESS)

                if is_highlighted:
                    linestyle = '-'
                    alpha = 1.0
                    linewidth = 1.5
                elif np.isclose(thickness, 1.0): # Hardcoding 1.0 as intermediate
                    linestyle = '--'
                    alpha = 0.6
                    linewidth = 1.2
                else: # Assuming the other is 0.5mm
                    linestyle = '-.'
                    alpha = 0.6
                    linewidth = 1.2

                # Plot on the main figure
                plt.plot(k_filtered, phi_k_filtered,
                         label=f'{material} {thickness:.1f} mm',
                         color=color,
                         linestyle=linestyle,
                         linewidth=linewidth,
                         alpha=alpha)
                plotted_filtered = True

                # Update global minimum y value found within the plotting range
                valid_indices = (k_filtered >= MIN_K_PLOT) & (k_filtered <= MAX_K_PLOT) & (phi_k_filtered * alpha > 1e-6)
                if np.any(valid_indices):
                     current_min = np.min(phi_k_filtered[valid_indices])
                     if current_min > 1e-6:
                         global_y_min_in_range = min(global_y_min_in_range, current_min)

            except Exception as e:
                 print(f"  Error processing {material} {thickness:.1f} mm: {e}", file=sys.stderr)

    # --- Final Plot Formatting ---
    if not plotted_filtered and k_unfiltered.size == 0:
         print("\nError: No spectra could be generated or plotted.", file=sys.stderr)
         plt.close()
         return # Exit if nothing was plotted

    print("\nFormatting final plot...")
    plt.xlabel('Energy [keV]')
    plt.ylabel('Differential fluence [photons cm$^{-2}$ keV$^{-1}$ per mAs at 100 cm]')
    plt.title(f'{KVP_SET:.0f} kVp, {ANODE_ANGLE:.0f} deg Anode Spectra - Filter Thickness Comparison')
    plt.yscale('log')
    plt.xlim(MIN_K_PLOT, MAX_K_PLOT)

    # Set Y limits
    if global_y_max_in_range > 0: # Check if reasonable max was found
        reasonable_bottom = max(1e1, global_y_min_in_range * 0.5 if global_y_min_in_range != np.inf else 1e1)
        reasonable_top = global_y_max_in_range * 1.5
        plt.ylim(bottom=reasonable_bottom, top=reasonable_top)
        print(f"Y-limits set based on data range: [{reasonable_bottom:.2e}, {reasonable_top:.2e}]")
    else:
        print("Warning: Could not determine valid Y range from data. Using default Y limits.")


    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    handles, labels = plt.gca().get_legend_handles_labels()
    if handles:
        plt.legend(handles, labels, fontsize=LEGEND_FONT_SIZE, ncol=LEGEND_COLUMNS)
    else:
        print("Info: No labeled data to add to legend.")

    plt.tight_layout()

    # --- Save the plot (PNG and PDF) ---
    try:
        # Construct filename dynamically
        base_filename = f'filtered_spectra_comparison_{KVP_SET:.0f}kv_{MIN_K_PLOT:.0f}-{MAX_K_PLOT:.0f}keV'
        # Use os.path.join to save into the subfolder (MODIFIED)
        png_filepath = os.path.join(OUTPUT_SUBFOLDER, base_filename + ".png")
        pdf_filepath = os.path.join(OUTPUT_SUBFOLDER, base_filename + ".pdf")

        plt.savefig(png_filepath, dpi=600)
        print(f"\nPlot saved as PNG: {png_filepath}")
        plt.savefig(pdf_filepath)
        print(f"Plot saved as PDF: {pdf_filepath}")
    except Exception as e:
        print(f"Error saving plot files: {e}", file=sys.stderr)

    # --- Show the plot ---
    plt.show()
    plt.close() # Close figure after showing/saving

    print('\nFinished fixed thickness comparison!\n')

# --- Main execution guard ---
if __name__ == "__main__":
    main()