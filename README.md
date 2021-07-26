# Package Scripts

A version controlled directory to store scripts related to various GitHub repos for R and Python packages.

## Usage

Pull this repo into the directory where you store GitHub repos.

Create a directory in this repo with the same name as the GitHub repo you want to store scripts for, if it doesn't already exist.

Change into the repo you want to store scripts for, then symbolic link the corresponding directory in this repo.
```sh
cd <pkg_name>
ln -s ../pkg_scripts/<pkg_name> `pwd`/scripts
```
Now add the symbolc link to your `.gitignore` file, so it doesn't get pushed to the repo.

```sh
echo "scripts/*" >> .gitignore
```

To update the <pkg_name> repo, just commit as usual. To update the scripts directory, `cd scripts` then make your commits. This will commit all your scripts to the pkg_scripts directory.

If you want your scripts to be independent from another user, make a directory with your username before making the <pkg_name> directory.

For example, my personal scripts (i.e., scripts I don't want  others to use or modify) for AnnotationGx would live at:
```sh
pkg_scripts/ceeles/AnnotationGx
```

And I would link this into AnnotationGx like this:
```sh
mkdir -p pkg_scripts/ceeles/AnnotationGx
cd AnnotationGx
ln -s ../pkg_scripts/ceeles/AnnotationGx $(pwd)/scripts
```
