[Back](./Index.md) | [Home](https://github.com/acg-team/ProPIP/blob/master/ProPIP.wiki/ProPIP-Progressive-Multiple-Sequence-Alignment-with-Poisson-Indel-Process.md)

## How to commit to the repository

This repository follows the [GitFlow Workflow Guidelines](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow). We use the following branches:

1. `master` (latest stable version of the project: here we apply only patches until the next release)
2. `develop` (sandbox for developing process)
3. `release/*` (wrapping up branch before release: it branches out from develop)
4. `features/*` (implementation of independent features: it branches out from develop)
5. `hotfix/*` (implementation of hotfix to apply directly in the master branch)

Commit messages should be self-contained, self-explanatory and they should indicate the extent of the applied modifications. We use the following keywords:

`#ADD: Added file/module/method/routine which affects the following workflows, etc.`

`#REM: Removed file/module/method/routine due to ... `

`#FIX: Fixed bug affecting file/module/method/routine in ... `

`#REF: Refactored file/module/method/routine. `

`#TEST: Tested functionalities and outcomes. Test passed/not-passed`

`#ROLL: Rolled back to version/commit/... due to ... `


## Software versioning and build numbers

Every software release should enter the following process before being released:

        Master -> Develop -> features/* -> Develop -> release/* -> Master

Software version numbering must follow the schema: `major.minor[.build[.revision]]`, i.e.

- `[1]`.`[0]`.`[2]` build on `[commit]`


## Continuous Integration (Automatic Testing via Bitbucket pipelines)

The continuous integration (CI) system that we use to perform post-commit functionality tests is composed by the following steps:

1. Copying personalized environment docker image from `lorenzogatti89/minijati:latest`
2. Pulling a copy of the repository and checkout of the branch where the commit has taken place
3. Preparation of the `make files`
4. Compilation of the project with the `release` settings
5. Execution of the registered pipelines


The basic pipelines are run automatically on the following branches (git-flow ready):

- `master`
- `develop`
- `feature/*`
- `release/*`
- `hotfix/*`

Bitbucket Pipelines supports skipping builds.
If you don't want to run a build, for example you want to save build minutes, you can include `[skip ci]` or `[ci skip]` anywhere in your commit message of the HEAD commit. Any commits that include `[skip ci]` or `[ci skip]` in the message are ignored by Pipelines.


## How to read the log files

The log files produced by the GLOG library. In order to increase/decrease the verbosity, the following environment variables must be exported prior to the execution:

`export GLOG_v={0,1,2,3}`

`export GLOG_minloglevel={INFO,FATAL,WARN}`


For debugging logs, the project must be compiled with the DEBUG flag active.


## Where to find re-usable example
