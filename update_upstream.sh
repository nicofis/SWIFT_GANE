#!/bin/bash
set -e
# Run from the repository root, wherever it is
cd "$(dirname "$0")"

# Add the official SWIFT repository as 'upstream' the first time
git remote get-url upstream >/dev/null 2>&1 || \
  git remote add upstream https://github.com/SWIFTSIM/SWIFT.git

# Refuse to run with uncommitted changes
if ! git diff --quiet || ! git diff --cached --quiet; then
  echo "Commit or stash your changes first."; exit 1
fi

git fetch upstream
git checkout master
git merge upstream/master
git push origin master