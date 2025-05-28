#!/usr/bin/env python
"""
Debug script to help understand PHACTS import issues.
This script checks the Python path and attempts to import phacts to
identify what's causing the ModuleNotFoundError.
"""

import os
import sys
import importlib.util
import traceback

def check_import_path(output_dir=None):
    """Check if PHACTS can be imported and diagnose any issues."""
    print("=" * 80)
    print("PHACTS Debugging Information")
    print("=" * 80)
    
    # Check Python version
    print(f"Python version: {sys.version}")
    print(f"Python executable: {sys.executable}")
    
    # Check current working directory
    print(f"Current working directory: {os.getcwd()}")
    
    # Check PYTHONPATH
    print("\nPYTHONPATH:")
    for path in sys.path:
        print(f"  - {path}")
    
    # Check if output_dir is provided
    if output_dir:
        phacts_path = os.path.join(output_dir, "db/phacts/PHACTS")
        phacts_parent = os.path.dirname(phacts_path)
        print(f"\nChecking provided output directory: {output_dir}")
        print(f"PHACTS expected path: {phacts_path}")
        
        # Check if directories exist
        print(f"Output directory exists: {os.path.exists(output_dir)}")
        print(f"PHACTS directory exists: {os.path.exists(phacts_path)}")
        
        # Check installation files
        phacts_py = os.path.join(phacts_path, "phacts.py") 
        print(f"phacts.py exists: {os.path.exists(phacts_py)}")
        print(f"__init__.py exists: {os.path.exists(os.path.join(phacts_path, '__init__.py'))}")
        
        # Add paths to Python path temporarily
        if os.path.exists(phacts_path):
            sys.path.insert(0, phacts_path)
            sys.path.insert(0, phacts_parent)
            print(f"\nAdded to sys.path: {phacts_path}")
            print(f"Added to sys.path: {phacts_parent}")
    
    # Try to import phacts
    print("\nAttempting to import 'phacts':")
    try:
        import phacts
        print(f"✅ Successfully imported phacts!")
        print(f"phacts location: {phacts.__file__}")
    except ImportError as e:
        print(f"❌ Failed to import phacts: {e}")
        print("\nTraceback:")
        traceback.print_exc()
        
    # Try an alternate import approach
    print("\nAttempting to find phacts module another way:")
    found_locations = []
    for path in sys.path:
        if not os.path.isdir(path):
            continue
        try:
            for item in os.listdir(path):
                if item == "phacts" or item == "PHACTS":
                    phacts_dir = os.path.join(path, item)
                    if os.path.isdir(phacts_dir):
                        found_locations.append(phacts_dir)
        except (PermissionError, FileNotFoundError):
            pass
    
    if found_locations:
        print(f"Found potential PHACTS directories:")
        for loc in found_locations:
            print(f"  - {loc}")
            # Check for initialization file
            if os.path.exists(os.path.join(loc, "__init__.py")):
                print(f"    Has __init__.py: Yes")
            else:
                print(f"    Has __init__.py: No")
    else:
        print("No 'phacts' directories found in Python path.")
    
    print("\nModules that are successfully imported:")
    for name, module in sys.modules.items():
        if "phacts" in name.lower():
            print(f"  - {name}: {getattr(module, '__file__', 'No file attribute')}")
            
    print("\n" + "=" * 80)
    print("End of Debugging Information")
    print("=" * 80)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Debug PHACTS import issues")
    parser.add_argument("--output-dir", help="The output directory where PHACTS is installed")
    args = parser.parse_args()
    
    check_import_path(args.output_dir)