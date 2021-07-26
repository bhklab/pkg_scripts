# Package Scripts

A version controller directory to store scripts related to various GitHub repos for R and Python packages.

## Usage

Pull this repo into the directory where you store GitHub repos.

Create a directory in this repo with the same name as the GitHub repo you want to store scripts for, if it doesn't already exist.

Change into the repo you want to store scripts for, then symbolic link the corresponding directory in this repo.
```sh
cd <pkg_name>
ln -s ../pkg_scripts/<pkg_name> `pwd`/scripts
```
