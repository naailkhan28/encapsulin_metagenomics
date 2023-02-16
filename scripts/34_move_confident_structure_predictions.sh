#!/bin/bash

while read -r filename
  do
  mv "$filename" confident/
done < confident_predicted_structures.txt