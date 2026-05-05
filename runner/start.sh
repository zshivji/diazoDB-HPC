#!/bin/bash
pkill -f runner.py 2>/dev/null
sleep 1
nohup python runner.py >> /tmp/runner.log 2>&1 &