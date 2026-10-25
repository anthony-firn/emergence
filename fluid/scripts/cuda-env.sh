# Source this to put the locally installed CUDA runtime libraries on the
# dynamic loader path:  source scripts/cuda-env.sh
_L="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/external/tools/cuda-libs/nvidia"
export LD_LIBRARY_PATH="$_L/cufft/lib:$_L/cuda_nvrtc/lib:$_L/cuda_runtime/lib:$_L/nvjitlink/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
unset _L
