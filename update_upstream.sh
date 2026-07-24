#!/bin/bash
#
cd ~/Desktop/SWIFT_GANE
git fetch upstream
git checkout main
git merge upstream/main
git push origin main
