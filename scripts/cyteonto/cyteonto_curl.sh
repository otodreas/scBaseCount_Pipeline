#!/usr/bin/env bash
curl -sS -X POST "https://cyteonto.nygen.io/compare" \
  -H 'Content-Type: application/json' \
  -d @"$(dirname "$0")/../../tmp/annotations.json" > "$(dirname "$0")/../../tmp/cyteonto_response.json"
