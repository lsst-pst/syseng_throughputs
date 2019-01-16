## This repository contains approved system engineering throughput components.

### Update process

Any updates to the repository must be approved by Steve Ritz (camera) or Chuck Claver (telescope).

Updates should be added in a branch and merged to master via a Pull Request.

Create a copy of the repo, if necessary.
```
git clone git@github.com:lsst-pst/syseng_throughputs syseng_throughputs
```

Create and checkout a new branch to add the information about the updated component <component_name>.
```
git checkout -b update/<component_name>
```
Note that the result of `git status` will now include the branch name. You can switch back to master by
```git checkout master``` and then back again by ```git checkout update/<component_name>```

Add your files or update the existing files in the repository, in this branch.

When ready to commit the files to git:
```
git add <my_updated_file>
```
and then commit these changes to git.
```git commit```
(you will be asked to include a git commit message .. make it descriptive of the changes, as this is a good way
to know what happened in this commit).

Ideally each commit includes a self-consistent "set" of files (i.e. all of the files related to a lens update,
along with any necessary documentation and python files), but this is not strictly necessary (we can
"squash" commits during a merge, assuming each PR & merge deals with a single change).

Please remember to update the top-level README, which summarizes changes to the throughputs repo for each
version of the files.

So these changes are local to your own copy of the repository -- you should now push them back to
github with
```git push```
to send to the remote branch. You can now go to github and start the Pull Request. Add additional documentation
in the PR, including plots of the previous and new throughput components, the source of the new throughput curves,
a short summary of the reason for the changes, and a summary of the effect of the changes. Then Steve and/or Chuck
can be added as the "Reviewers" -- at least one approval is required before the branch can be merged to master.

After it is merged to master, a new version # / tag should be assigned, to uniquely identify this "version" of the
throughputs repo. This will also be copied to the lsst/throughputs repo.


### Files in the repository

Tracking the inputs in the syseng_throughputs repository from the official files in Docushare is important.
For this reason, notebooks which CREATE the actual text throughput files from the Docushare input file will be
stored in the repository as well. These will be placed in the input_notebooks directory, together with a copy of
the Docushare file (or a link to the official docushare file?). These notebooks will be configured to write their
outputs directly to the appropriate directories in the components directory. If the format of the files changes as
a result, the corresponding code in the python/lsst/syseng/throughputs repository (which reads and creates the combined files)
should be updated in the same PR. There will be "integration tests" as part of the repository, which should be run
before a PR: these will exercise the python code to read the throughput files and present their outputs to
a user, as well as calculating <XXX> values to compare to previous/expected values.

