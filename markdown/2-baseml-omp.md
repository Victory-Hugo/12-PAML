```sh
cd /mnt/f/OneDrive/文档（科研）/脚本/Download/12-PAML/src
gcc -O3 -std=c11 -fopenmp -o baseml_LLST baseml_LLST.c tools.c -lm
mv baseml_LLST ../bin/
```