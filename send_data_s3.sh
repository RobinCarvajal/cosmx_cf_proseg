# frist argument is the source, second argument is the destination

rclone sync \
  /mnt/data/project0062/proseg_data \ 
  gray_s3:gray-lab-864138843198-eu-west-2-an/proseg \ 
  --dry-run -P


rclone copy /mnt/data/project0062/sc_pediatric \
  gray_s3:gray-lab-864138843198-eu-west-2-an/sc_pediatric \
  --dry-run -P
