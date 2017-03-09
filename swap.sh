#!/bin/bash
# Make swap

size=64G
path=/tmp/swapfile

# allocate file                  permissions             makeswap             enable swap
sudo fallocate -l $size $path && sudo chmod 600 $path && sudo mkswap $path && sudo swapon $path
