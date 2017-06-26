#!/bin/bash
echo "" > ./$1.txt
grep -r "Commit OSPRay patch" ./$1/*.timings >> ./$1.txt
grep -r "Render OSPRay patch" ./$1/*.timings >> ./$1.txt
grep -r "AllPatchRendering" ./$1/*.timings >> ./$1.txt
grep -r "Sample point extraction" ./$1/*.timings >> ./$1.txt
grep -r "Compositing" ./$1/*.timings >> ./$1.txt
grep -r "Ray Tracing" ./$1/*.timings >> ./$1.txt
grep -r "avtRayTracer" ./$1/*.timings >> ./$1.txt
