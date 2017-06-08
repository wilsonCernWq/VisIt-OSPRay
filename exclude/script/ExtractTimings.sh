#!/bin/bash
echo "" > ./$1.txt
grep -r "Commit OSPRay patch" ./$1/* >> ./$1.txt
grep -r "Render OSPRay patch" ./$1/* >> ./$1.txt
grep -r "OSPRayRendering" ./$1/* >> ./$1.txt
grep -r "Sample point extraction" ./$1/* >> ./$1.txt
grep -r "Compositing" ./$1/* >> ./$1.txt
grep -r "Ray Tracing" ./$1/* >> ./$1.txt
grep -r "avtRayTracer" ./$1/* >> ./$1.txt
