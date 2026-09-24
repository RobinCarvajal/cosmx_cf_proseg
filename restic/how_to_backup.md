


```bash

project_dir="/mnt/data/project0062/cosmx_cf_orig"
tag_name=$(basename "$project_dir")
host_name="mars"

restic "$backup project_dir" --tag "$tag_name" --host "$host_name"

```

```bash

project_dir="/mnt/data/project0062/cosmx_cf_proseg"
tag_name=$(basename "$project_dir")
host_name="mars"

restic "$backup project_dir" --tag "$tag_name" --host "$host_name"

```