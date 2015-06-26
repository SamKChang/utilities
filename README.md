This is the utility codes

The codes contains C, bash, perl, python 
at the moment of creation (20150626). 
Each code was written for specific project 
but it was also written such that 
it can be generally applied for future 
chemistry/physics applications.

All the codes are recollected from local machine
and uploaded to GitHub

The following commands were used:
create online repository at GitHub with repo_name.git

cd /path/to/local/repository

git init

git add *

git commit -m 'message for history'

git remote add upload_name https://github.com/user_name/repo_name.git


git push upload_name master

setup sshkey

git remote set-url upload_name git@github.com:user_name/repo_name.git

set up default push
git config --global push.default current

workflow: 
git add *
git commit -m "message"
git push or pull

NOTE: if the repo_name and local folder name are different, default "git push" command will require remote name to upload with "push.config current" setupt

NOTE: for default "git pull" command, .git/config needs to be modified to specify remote and branch by adding

[branch "master"]

    remote = origin
    
    merge = refs/heads/master
