#!/bin/bash
set -euo pipefail
#!使用前请一定配置好ctl文件!!!!
#!使用前请一定配置好ctl文件!!!!
#!使用前请一定配置好ctl文件!!!!
# 定义临时目录

BASEML_PATH="/mnt/f/OneDrive/文档（科研）/脚本/Download/12-PAML/bin/baseml_parra"
CTL_FILE="/mnt/f/OneDrive/文档（科研）/脚本/Download/12-PAML/conf/baseml_parra.ctl"
OUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/12-PAML/output"


# 拷贝所需ctl文件
cp "$CTL_FILE" \
    "$OUT_DIR/"

# 切换到输出目录
cd "$OUT_DIR"

# 运行PAML baseml
"${BASEML_PATH}" \
    "${CTL_FILE}"

rm "${OUT_DIR}/baseml2.ctl" \
   "${OUT_DIR}/lnf" \
   "${OUT_DIR}/rst" \
   "${OUT_DIR}/rst1" \
   "${OUT_DIR}/rst2" \
   "${OUT_DIR}/rub" \
   "${OUT_DIR}/2base.t"