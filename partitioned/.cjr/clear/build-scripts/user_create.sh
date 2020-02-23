#!/bin/bash

# Change Root Password - or podman cannot execute USER root after USER <user>
passwd root << EOF
$USER_PASSWORD
$USER_PASSWORD
EOF

# -- Install Sudo Package ------------------------------------------------------
swupd bundle-add sudo

# -- Add User-------------------------------------------------------------------
if [ -z "$USER_ID" ] || [ -z "GROUP_ID" ] ; then
  useradd -m -l -s /bin/bash $USER_NAME
else
  groupadd --gid $GROUP_ID $USER_NAME
  useradd -m -l -o -s /bin/bash --uid $USER_ID --gid $GROUP_ID $USER_NAME
  passwd $USER_NAME << EOL
$USER_PASSWORD
$USER_PASSWORD
EOL
fi

# -- Grant sudo ----------------------------------------------------------------
if [ "$GRANT_SUDO" = "TRUE" ] ; then
  (usermod -aG wheel $USER_NAME)
fi
