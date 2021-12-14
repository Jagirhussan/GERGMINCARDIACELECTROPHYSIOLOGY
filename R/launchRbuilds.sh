#!/bin/bash
nohup  Rscript ./buildnetscript.R 0 2>&1 | tee trigNet0.txt &
nohup  Rscript ./buildnetscript.R 1 2>&1 | tee trigNet1.txt &
nohup  Rscript ./buildnetscript.R 2 2>&1 | tee trigNet2.txt &
nohup  Rscript ./buildnetscript.R 3 2>&1 | tee trigNet3.txt &
nohup  Rscript ./buildnetscript.R 4 2>&1 | tee trigNet4.txt &
nohup  Rscript ./buildnetscript.R 5 2>&1 | tee trigNet5.txt &
