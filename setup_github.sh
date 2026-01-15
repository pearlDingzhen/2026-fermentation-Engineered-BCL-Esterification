#!/bin/bash

# GitHub Repository Setup Script
# Run this script to create and push to GitHub

# ========== Configuration ==========
# Replace these values with your information
GITHUB_USERNAME="your-github-username"
REPO_NAME="2026-fermentation-Engineered-BCL-Esterification"
REPO_DESCRIPTION="Computational and experimental data for engineered BCL esterification research"

# ========== Step 1: Create GitHub Repository ==========
# Option A: Using GitHub CLI (if installed)
# gh repo create "$REPO_NAME" --public --description "$REPO_DESCRIPTION"

# Option B: Manual creation via web browser
echo "Please create a new repository on GitHub:"
echo "1. Go to: https://github.com/new"
echo "2. Repository name: $REPO_NAME"
echo "3. Description: $REPO_DESCRIPTION"
echo "4. Select: Public"
echo "5. DO NOT initialize with README, .gitignore, or license"
echo "6. Click 'Create repository'"
echo ""
echo "After creating, copy the repository URL (e.g., https://github.com/$GITHUB_USERNAME/$REPO_NAME.git)"

# ========== Step 2: Add Remote and Push ==========
# Replace YOUR_GITHUB_USERNAME and YOUR_REPO_NAME with actual values
REMOTE_URL="https://github.com/YOUR_GITHUB_USERNAME/YOUR_REPO_NAME.git"

echo ""
echo "Adding remote origin..."
git remote add origin "$REMOTE_URL"

echo ""
echo "Pushing to GitHub..."
git branch -M main
git push -u origin main

echo ""
echo "Repository setup complete!"
echo "Your data is now available at: https://github.com/$GITHUB_USERNAME/$REPO_NAME"

# ========== Optional: Add DOI Badge ==========
echo ""
echo "To add a DOI badge to your README, consider registering your dataset with:"
echo "- Zenodo (https://zenodo.org)"
echo "- Figshare (https://figshare.com)"
echo "- OSF (https://osf.io)"
