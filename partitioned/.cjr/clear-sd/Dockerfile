FROM clearlinux:latest

# -- install root dependencies -------------------------------------------------
RUN swupd bundle-add c-basic gdb openblas devpkg-fftw vim git sudo sysadmin-basic

# -- User arguments ------------------------------------------------------------
# USER_ID should be set to $(id -u) and GROUP_ID to $(id -g) for correct bind mounts permissions
ARG USER_NAME=user
ARG USER_PASSWORD=password
ARG GRANT_SUDO=false
ARG USER_ID
ARG GROUP_ID

# -- Generate restricted user account ------------------------------------------
RUN echo -e "$USER_PASSWORD\n$USER_PASSWORD" | passwd root
RUN if [ -z "$USER_ID" ] || [ -z "GROUP_ID" ] ; then \
useradd -m -l -s /bin/bash $USER_NAME ; \
else \
groupadd --gid $GROUP_ID $USER_NAME ; \
useradd -m -l -o -s /bin/bash --uid $USER_ID --gid $GROUP_ID $USER_NAME ; \
echo -e "$USER_PASSWORD\n$USER_PASSWORD" | passwd $USER_NAME ; \
fi
RUN if [ "$GRANT_SUDO" = "TRUE" ] ; then \
(usermod -aG wheel $USER_NAME); \
fi

# -- Set up user user dependancies ---------------------------------------------
WORKDIR /home/$USER_NAME
USER $USER_NAME

ENTRYPOINT ["bash", "-c"]
