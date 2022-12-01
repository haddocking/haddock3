local help_message = [[

HADDOCK3 environment based on miniconda with Python3.9

]]

help(help_message,"\n")


whatis("Name: haddock")
whatis("Version: 3.0.0")
whatis("Category: python conda haddock3")
whatis("Keywords: python conda haddock3")
whatis("Description: haddock3 miniconda3 python3.9 environment")

local conda_dir = "/trinity/home/software/miniconda3"
local funcs = "conda __conda_activate __conda_hashr __conda_reactivate __add_sys_prefix_to_path"

-- Specify where system and user environments should be created
setenv("CONDA_ENVS_PATH", conda_dir .. "/envs")
-- Directories are separated with a comma
setenv("CONDA_PKGS_DIRS", conda_dir .. "pkgs")
-- Initialize conda and activate environment
execute{cmd="source " .. conda_dir .. "/etc/profile.d/conda.sh; conda activate haddock3", modeA={"load"}}
-- Unload environments and clear conda from environment
execute{cmd="for i in $(seq ${CONDA_SHLVL:=0}); do conda deactivate; done; pre=" .. conda_dir .. "; \
        export LD_LIBRARY_PATH=$(echo ${LD_LIBRARY_PATH} | tr ':' '\\n' | grep . | grep -v $pre | tr '\\n' ':' | sed 's/:$//'); \
        export PATH=$(echo ${PATH} | tr ':' '\\n' | grep . | grep -v $pre | tr '\\n' ':' | sed 's/:$//'); \
        unset -f " .. funcs .. "; \
        unset $(env | grep -o \"[^=]*CONDA[^=]*\");", modeA={"unload"}}
-- Prevent from being loaded with another system python or conda environment
family("python")
