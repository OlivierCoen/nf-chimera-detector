# Installation

## Nextflow

If your Java version is >= 17, you can install Nextflow in your `$HOME` folder using:

```bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow
mkdir -p $HOME/.local/bin/
mv nextflow $HOME/.local/bin/
echo "export PATH=$HOME/.local/bin/:$PATH" >> $HOME/.bashrc
```

> [!NOTE]
> If your Java version is < 17, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow.

> [!NOTE]
> You can also install Nextflow system-wide by moving the `nextflow` executable to a directory accessible to all users:
> ```bash
> sudo mv nextflow /usr/local/bin/
> ```

## Apptainer

When running the pipeline of multi-user server or on a cluster, the best practice is to use Apptainer (formerly Singularity).

You can install Apptainer by following these [instructions](https://apptainer.org/docs/admin/main/installation.html#).

In case you encounter the following error when running Apptainer:
```
ERROR  : Could not write info to setgroups: Permission denied
ERROR  : Error while waiting event for user namespace mappings: no event received
```
you may need to install the `apptainer-suid` package instead of `apptainer`:

```
# Debian / Ubuntu
sudo apt install apptainer-suid
# RHEL / CentOS
sudo yum install apptainer-suid
# Fedora
sudo dnf install apptainer-suid
```


## Docker / Podman / Shifter / Charliecloud

You can also use other container technologies to run the pipeline. See the following instructions to install:
* [Docker](https://docs.docker.com/engine/install/)
* [Podman](https://podman.io/docs/installation)
* [Shifter](https://shifter.readthedocs.io/en/latest/install_guides.html)
* [Charliecloud](https://charliecloud.io/latest/install.html)

> [!NOTE]
> For now, the only one that is tested is Docker. If you encounter any problem, please open an issue.


## Micromamba / Conda

This pipeline can be launched with Micromamba (or Conda) too and it has been tested successfully.

Although it is not the best practice, it should work perfectly and it can be very convenient when you cannot run containers (Apptainer, Docker, Podman, ...).

* To install Micromamba in your `$HOME` folder:

```bash
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```
(or see these [instructions](https://micromamba.readthedocs.io/en/latest/installation.html)).

* To install Miniconda in your `$HOME` folder, see these [instructions](https://www.anaconda.com/docs/getting-started/miniconda/install).
