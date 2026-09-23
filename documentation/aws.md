
# Backing up rationale

> The idea is to work mostly on the EC2 instance, and regularly update the S3 bucket. 

### EC2 to S3 (dry)

```bash

aws s3 sync \
    /data/proseg_data/ \
    s3://gray-lab-864138843198-eu-west-2-an/proseg/ \
    --exclude ".git/*" \
    --exclude "*/.git/*" \
    --dryrun

```

### EC2 to S3

```bash

aws s3 sync \
    /data/proseg_data/ \
    s3://gray-lab-864138843198-eu-west-2-an/proseg_data/ \
    --exclude ".git/*" \
    --exclude "*/.git/*" \


```

> If there are GPU intensive tasks, it's time to update the HPC disk and perform those tasks there, either scheduling a job or using jupyter. 

### S3 to HPC (dry)

```bash 

rclone sync \
    gray_s3:gray-lab-864138843198-eu-west-2-an/proseg_data \
    /mnt/data/project0062/proseg_data \
    -P \
    --dry-run \
    --exclude ".git/"

```

### S3 to HPC

```bash 

local_data="/mnt/data/project0062/proseg_data"
aws_data="gray_s3:gray-lab-864138843198-eu-west-2-an/proseg"

rclone copy "${aws_data}" "${local_data}" -P --dry-run --exclude ".git/"

```

> Once the job is finished, the S3 bucket should be updated. P.D. We exclude git files.

### HPC to S3 (dry)

```bash 

rclone sync \
    /mnt/data/project0062/proseg_data \
    gray_s3:gray-lab-864138843198-eu-west-2-an/proseg \
    -P \
    --dry-run \
    --exclude ".git/"

local_data="gray_s3:gray-lab-864138843198-eu-west-2-an/proseg_data"
aws_data="/mnt/data/project0062/proseg"

rclone copy "${local_data}" "${aws_data}" -P --dry-run --exclude ".git/"

```

### S3 to HPC

```bash 

rclone copy \
    /mnt/data/project0062/proseg_data \
    gray_s3:gray-lab-864138843198-eu-west-2-an/proseg \
    -P \
    --exclude ".git/"

```

> Lastly, we update the EC2 instance with all the changes changes saved on the S3 bucket.

### S3 to EC2 (dry)

```bash

aws s3 sync \
    s3://gray-lab-864138843198-eu-west-2-an/proseg_data/ \
    /data/proseg_data/ \
    --exclude ".git/*" \
    --exclude "*/.git/*" \
    --dryrun

```

### S3 to EC2 

```bash

aws s3 sync \
    s3://gray-lab-864138843198-eu-west-2-an/proseg_data/ \
    /data/proseg_data/ \
    --exclude ".git/*" \
    --exclude "*/.git/*" \

```

### Golden rule:

Always push your changes to S3 before switching machines, and always pull from S3 before starting work on the other machine.

```text

                    aws s3 sync
           ┌────────────────────────►┐
           │                         │
     ┌─────┴─────┐             ┌─────┴─────┐
     │           │             │           │
     │    EC2    │             │    S3     │
     │           │             │           │
     └─────▲─────┘             └─────┬─────┘
           │                         │
           │                         │ rclone sync
aws s3 sync│                         │
           │                         ▼
           │                   ┌───────────┐
           │                   │           │
           └────── S3 ◄────────│    HPC    │
                rclone copy    │           │
                               └───────────┘

```
