#!/bin/bash
set -e

# Run from project root
cd "$(dirname "$0")/.."

# Run all .vd tests
for testfile in tests/*.vd; do
    basename="${testfile%.vd}"
    echo "Running test: $basename"
    
    # Run the test and save output
    vd -b "$testfile"
    
    # Compare output with golden file
    if ! diff -u "tests/golden/$(basename "$basename").tsv" "$basename.tsv"; then
        echo "Test failed: $basename"
        exit 1
    fi
    
    echo "Test passed: $basename"
done

echo "All tests passed!"
