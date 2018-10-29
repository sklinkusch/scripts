#!/bin/bash

url="$1"

browser-sync start --directory --server --files "**/*.js, **/*.html, **/*.css" $url | ~/bin/bs-filter.pl
