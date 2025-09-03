#!/usr/bin/env python3
import os
import re
import sys
from datetime import datetime

import tomli


def get_version_from_toml(file_path):
    """Extract version from pyproject.toml"""
    data = {}
    try:
        with open(file_path, "rb") as f:
            data = tomli.load(f)
        return data["tool"]["poetry"]["version"]
    except KeyError:
        # Try different possible locations for version
        try:
            return data["project"]["version"]
        except KeyError:
            print(f"Error: Could not find version in {file_path}")
            sys.exit(1)


def is_valid_yyyymmpp(version_str):
    """Check if version string follows YYYY.MM.PP format"""
    pattern = r"^(\d{4})\.(\d{1,2})\.(\d{1,2})$"
    match = re.match(pattern, version_str)
    if not match:
        return False

    year, month, patch = map(int, match.groups())

    # Validate date components
    if not (1900 <= year <= 2100):  # Reasonable year range
        return False
    if not (1 <= month <= 12):  # Valid month
        return False
    if not (0 <= patch <= 99):  # Reasonable patch range
        return False

    return True


def parse_yyyymmpp_version(version_str):
    """Parse YYYY.MM.PP version into components"""
    pattern = r"^(\d{4})\.(\d{1,2})\.(\d{1,2})$"
    match = re.match(pattern, version_str)
    if not match:
        raise ValueError(f"Invalid YYYY.MM.PP format: {version_str}")

    year, month, patch = map(int, match.groups())
    return year, month, patch


def compare_versions(base_version, pr_version):
    """Compare versions and check if PR version is an increment"""
    if not is_valid_yyyymmpp(pr_version):
        print(f"❌ PR version '{pr_version}' is not a valid YYYY.MM.PP version")
        print("Expected format: YYYY.MM.PP (e.g., 2025.8.1)")
        return False

    try:
        base_year, base_month, base_patch = parse_yyyymmpp_version(base_version)
        pr_year, pr_month, pr_patch = parse_yyyymmpp_version(pr_version)

        # Get current date for validation
        current_date = datetime.now()
        current_year = current_date.year
        current_month = current_date.month

        # Check if year is current
        if pr_year < current_year:
            print(
                f"❌ Error: PR version year {pr_year} is in the past compared to current year {current_year}"
            )
            return False

        # Check if month is current
        if pr_month < current_month:
            print(
                f"❌ Error: PR version month {pr_month} is in the past compared to current month {current_month}"
            )
            return False

        # Check if patch version is higher
        if pr_patch <= base_patch and pr_year == base_year and pr_month == base_month:
            print(
                f"❌ Error: new version {pr_patch} is not greater than base patch {base_patch}"
            )
            return False

    except Exception as e:
        print(f"❌ Error comparing versions: {e}")
        return False

    return True


def main():
    # Get the base branch version (main)
    os.system("git fetch origin main:main")

    # Get PR version
    pr_version = get_version_from_toml("pyproject.toml")
    print(f"PR version: {pr_version}")

    # Get base version from main branch
    os.system("git checkout main -- pyproject.toml 2>/dev/null || true")
    base_version = get_version_from_toml("pyproject.toml")
    print(f"Base version: {base_version}")

    # Validate base version format
    if not is_valid_yyyymmpp(base_version):
        print(
            f"⚠️  Warning: Base version '{base_version}' doesn't follow YYYY.MM.PP format"
        )

    # Check if version is incremented
    if compare_versions(base_version, pr_version):
        print("✅ Version increment check passed!")
        sys.exit(0)
    else:
        print("❌ Version increment check failed!")
        sys.exit(1)


if __name__ == "__main__":
    main()
