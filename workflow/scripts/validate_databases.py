#!/usr/bin/env python3
"""
Database Validation Utility for Phage Analysis Pipeline

This script validates database configurations and provides helpful guidance
for setting up missing or incorrectly configured databases.

Usage:
    python workflow/scripts/validate_databases.py config/my_config.yaml
    python workflow/scripts/validate_databases.py config/my_config.yaml --entry-point viral_contigs
"""

import os
import sys
import argparse
import yaml
from pathlib import Path

def load_config(config_path):
    """Load and parse the YAML configuration file."""
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        return config
    except FileNotFoundError:
        print(f"❌ Error: Config file not found: {config_path}")
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"❌ Error: Invalid YAML in config file: {e}")
        sys.exit(1)

def is_placeholder_path(path):
    """Check if a path is likely a placeholder that should be ignored."""
    if not path:
        return True
    placeholder_indicators = ["/path/to/", "PLACEHOLDER", "YOUR_", "CHANGE_ME", "/your/path/"]
    return any(indicator in str(path) for indicator in placeholder_indicators)

def get_databases_for_entry_point(entry_point):
    """Return list of databases required for a given entry point."""
    databases = {
        "assembly": ["genomad", "checkv"],
        "reneo_output": ["genomad", "checkv"],
        "viral_contigs": ["genomad", "checkv"],
        "predicted_phages": [],
        "clustered_sequences": []
    }

    # Databases needed for analysis phase (all except clustered_sequences)
    if entry_point != "clustered_sequences":
        databases[entry_point].extend(["iphop", "phabox", "vcontact3"])

    return databases.get(entry_point, [])

def validate_database(db_name, db_config, required=True, auto_download_enabled=False):
    """
    Validate a single database configuration.

    Returns:
        tuple: (is_valid, status_message, recommendation)
    """
    if not db_config or not isinstance(db_config, dict):
        if required:
            return False, "❌ Missing configuration", f"Add {db_name} database configuration to your config file"
        else:
            return True, "⚠️  Optional (not configured)", "No action needed - this database is optional"

    db_path = db_config.get("db", "")
    auto_download = db_config.get("auto_download", False)

    # Check for placeholder paths
    if is_placeholder_path(db_path):
        if auto_download:
            return True, "⚠️  Placeholder path with auto-download enabled", "Database will be downloaded automatically when needed"
        else:
            return False, "❌ Placeholder path detected", f"Update path in config or enable auto_download for {db_name}"

    # Check if database exists
    if os.path.exists(db_path):
        # Additional checks for database completeness could go here
        return True, "✅ Found and accessible", f"Database ready at {db_path}"
    else:
        if auto_download:
            return True, "⚠️  Missing but auto-download enabled", "Database will be downloaded automatically when needed"
        else:
            download_commands = {
                "genomad": f"genomad download-database {db_path}",
                "checkv": f"checkv download-database {db_path}",
                "iphop": f"iphop download --out_dir {db_path}",
                "phabox": "Follow instructions at: https://github.com/KennthShang/PhaBox2",
                "vcontact3": "Follow instructions at: https://github.com/vcontact/vcontact3"
            }
            cmd = download_commands.get(db_name, f"Check documentation for {db_name} download instructions")
            return False, "❌ Database not found", f"Download with: {cmd}"

def main():
    parser = argparse.ArgumentParser(description="Validate database configurations for phage analysis pipeline")
    parser.add_argument("config", help="Path to configuration YAML file")
    parser.add_argument("--entry-point", default="assembly",
                       choices=["assembly", "reneo_output", "viral_contigs", "predicted_phages", "clustered_sequences"],
                       help="Pipeline entry point to validate databases for (default: assembly)")
    parser.add_argument("--check-all", action="store_true",
                       help="Check all databases regardless of entry point")

    args = parser.parse_args()

    # Load configuration
    config = load_config(args.config)
    databases_config = config.get("databases", {})

    if not databases_config:
        print("❌ Error: No databases section found in config file")
        sys.exit(1)

    print(f"🔍 Validating database configuration for entry point: {args.entry_point}")
    print(f"📁 Config file: {args.config}")
    print("=" * 80)

    # Determine which databases to check
    if args.check_all:
        databases_to_check = list(databases_config.keys())
    else:
        databases_to_check = get_databases_for_entry_point(args.entry_point)

        # Add databases that are configured and have features enabled
        if config.get("run_iphop", True) and "iphop" not in databases_to_check:
            databases_to_check.append("iphop")
        if config.get("run_consensus", True) and "phabox" not in databases_to_check:
            databases_to_check.append("phabox")
        if config.get("run_consensus", True) and "vcontact3" not in databases_to_check:
            databases_to_check.append("vcontact3")

    # Validate each database
    issues_found = False
    all_results = []

    for db_name in databases_to_check:
        db_config = databases_config.get(db_name)
        required = db_name in get_databases_for_entry_point(args.entry_point)

        is_valid, status, recommendation = validate_database(db_name, db_config, required)

        print(f"\n📂 {db_name.upper():12} {status}")
        if db_config:
            print(f"   Path: {db_config.get('db', 'Not specified')}")
            print(f"   Auto-download: {db_config.get('auto_download', False)}")
        print(f"   → {recommendation}")

        if not is_valid:
            issues_found = True

        all_results.append((db_name, is_valid, status, recommendation))

    # Summary
    print("\n" + "=" * 80)
    if issues_found:
        print("❌ VALIDATION FAILED: Issues found with database configuration")
        print("\n🔧 Recommended actions:")
        for db_name, is_valid, status, recommendation in all_results:
            if not is_valid:
                print(f"   • {db_name}: {recommendation}")
        print("\n📖 For detailed setup instructions, see README.md 'Database Setup' section")
        sys.exit(1)
    else:
        print("✅ VALIDATION PASSED: All required databases are properly configured")
        print(f"\n🚀 Your pipeline is ready to run with entry point '{args.entry_point}'!")
        sys.exit(0)

if __name__ == "__main__":
    main()