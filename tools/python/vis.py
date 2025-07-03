#!/usr/bin/env python3
import argparse
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
import os

def parse_slice_string(slice_str):
    """Convert a string like ':,15,:' to a tuple of slice objects."""
    parts = slice_str.split(',')
    slice_list = []
    
    for part in parts:
        part = part.strip()
        if part == ':':
            slice_list.append(slice(None))
        elif ':' in part:
            slice_parts = part.split(':')
            start = int(slice_parts[0]) if slice_parts[0] else None
            stop = int(slice_parts[1]) if len(slice_parts) > 1 and slice_parts[1] else None
            step = int(slice_parts[2]) if len(slice_parts) > 2 and slice_parts[2] else None
            slice_list.append(slice(start, stop, step))
        else:
            slice_list.append(int(part))
    
    return tuple(slice_list)

def main():
    parser = argparse.ArgumentParser(description='Extract 2D slices from NIFTI files stack them horizontally and vis')
    parser.add_argument('input_files', nargs='+', help='Input NIFTI file(s)')
    parser.add_argument('--slice', '-s', required=True, 
                       help='Slice specification (e.g., ":,15,:" or "25,:,:" or ":,:,30")')
    parser.add_argument('--output', '-o', help='Output image file for the complete figure')
    
    args = parser.parse_args()
    
    # Parse slice once
    slice_tuple = parse_slice_string(args.slice)
    
    # Collect all slices and successful filenames
    slices = []
    successful_files = []
    
    # Process each NIFTI file
    for i, input_file in enumerate(args.input_files):
        print(f"\nProcessing file {i+1}/{len(args.input_files)}: {input_file}")
        
        try:
            # Load NIFTI file
            img = nib.load(input_file)
            data = img.get_fdata()
            
            print(f"  Data shape: {data.shape}")
            
            # Apply slice
            slice_2d = data[slice_tuple]
            print(f"  Slice {args.slice} -> shape: {slice_2d.shape}")
            
            # Ensure it's 2D
            if slice_2d.ndim != 2:
                print(f"  Warning: Result is {slice_2d.ndim}D, not 2D!")
                if slice_2d.ndim == 1:
                    slice_2d = slice_2d.reshape(-1, 1)
                elif slice_2d.ndim > 2:
                    print("  Squeezing extra dimensions...")
                    slice_2d = np.squeeze(slice_2d)
                    if slice_2d.ndim != 2:
                        raise ValueError(f"Could not reduce to 2D (shape: {slice_2d.shape})")
            
            slices.append(slice_2d)
            successful_files.append(input_file)
            
        except Exception as e:
            print(f"  Error processing {input_file}: {e}")
            continue
    
    if len(slices) == 0:
        print("No slices were successfully extracted!")
        return 1
    
    # Check if all slices have the same height
    heights = [s.shape[0] for s in slices]
    if len(set(heights)) > 1:
        raise ValueError(f"Slices have different heights: {heights}. All slices must have the same shape for stacking.")
    
    # Stack horizontally
    stacked = np.hstack(slices)
    print(f"\nStacked result shape: {stacked.shape}")
    print(f"Individual slice shapes: {[s.shape for s in slices]}")

    # Create title with filenames
    basenames = [os.path.basename(f) for f in successful_files]
    title = f'Slice "{args.slice}" from {", ".join(basenames)}'
    
    # Create and display the figure
    fig = plt.figure(figsize=(15, 5))
    plt.imshow(stacked, cmap='gray', aspect='auto')
    plt.title(title)
    
    # Add vertical lines to show boundaries
    x_pos = 0
    for i, s in enumerate(slices[:-1]):
        x_pos += s.shape[1]
        plt.axvline(x=x_pos, color='red', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    
    # Save if requested (before showing to ensure exact same figure is saved)
    if args.output:
        fig.savefig(args.output, dpi=300, bbox_inches='tight')
        print(f"Saved figure to: {args.output}")
    
    # Always show the figure
    plt.show()
    
    return 0

if __name__ == '__main__':
    exit(main())
