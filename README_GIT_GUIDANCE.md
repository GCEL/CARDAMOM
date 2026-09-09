**DRAFT VERSION TO BE CHANGED FROM JPL FORMAT**
**This document was adapted for University of Edinbugh CARDAMOM Github based on the JPL-Stanford-UCSB CARDAMOM_GIT_MUST_READ.md, available at https://github.com/CARDAMOM-framework/**

# Two sections in this file:
1. CARDAMOM user GitHub instructions (New users are required to go through this section)
2. Instructions to correct an unintentional git commit and git push that bring files to an older version, or unintentionally stashed changes, committed and pushed to the master
3. A suggested protocol for collaborative development of the repository

***Original authors: Shuang Ma and others in the JPL team***
***University of Edinburgh modifiers: T. Luke Smallman, David T. Milodowski***

# 1. CARDAMOM user GitHub instructions (New users are required to go through this section)
## 1) MUST READ: 
* Always work in a separate branch, and merge back to master via pull / merge request (examples below). Never commit directly to master!
* Do NOT make ‘unnecessary’ changes to either the R or Fortran source codes, eg. change ‘paths’ in the code, as these changes do not apply to other users but will be pushed to the master. Project specific set up information, including paths, to different datasets are in the R control file (see ./example_files for examples). You should make a separate copy which is NOT tracked by github. If you have a new dataset of an existing data type please reformat to be consistent with an existing CARDAMOM compatible file format. Ask a team member for help if you are not sure. 
* Do NOT commit changes to the master branch unless you have just pulled from the master branch to merge any subsequent changes into whatever you have. If you attempt to pull, and are told there are conflicts, then DEFINITELY DO NOT commit at that point, because those conflicts will get propagated into the master branch. Instead, you should be committing to a separate branch, and then make a pull request to merge your branch back into master, which will let you resolve any conflicts on a line-by-line basis. 
## 2) Split out your own branch.
## 3) Edit your own branch until functional - regularly pull from master to keep your branch updated during the development.
* If you don't ensure you keep pace with the master you will have a much harder time merging when you have completed your developments
## 4) When you’re ready to merge your changes to the master:
1) Update your branch from the main line (resolve conflicts). 
2) And then pull request to the main line (This is where we could have someone else double check, still under discussion on details of action).
## 5) Check that example analyses and alternate models work correctly
* i.e. ensure that even though your developments work well for what you are doing that they also work for everyone else too!
## 6) Tips:

### To the CARDAMOM respository (NOTE: do not include the quotes)
1) git clone 'remote-repo-url'

### To clone a specific branch (NOTE: do not include the quotes)
1) git clone --branch 'branchname' --single-branch 'remote-repo-url'

### To create a new branch from your local and then push to the remote:
1) git status                 # there shouldn’t be any red or green files which are part of the framework before proceeding, if there are, commit them or resolve conflicts as needed
2) git checkout –b branchname # Create a new branch on your local
3) git push -u origin HEAD    # This step push your new branch to remote and track it
 
### If you changed your branch name from github (repo website), you need to manually change the name locally:
1) git branch -m oldbranch newbranch
2) git fetch origin
3) git branch -u origin/ newbranch newbranch  # This step is equivalent to git push -u origin HEAD

### To sync your local dev_branch to the latest master on remote:
1) git checkout dev_branch # Make sure you are on your local dev_branch
2) git pull origin master  # It will pull master from remote to your local master, then merge to your local dev_branch

### To push your local development branch to remote development branch:
1) git checkout dev_branch    # To stay on your local branch 
2) git push origin dev_branch # It pushes your local dev_branch to the remote dev_branch; ‘origin’ always means the remote repo

### To check which branch you are currently on (and those available)
1) git branch -a 

### To push (merge) your development branch changes to the master
1) git branch -a              # Check which branch you are on
2) git checkout dev_branch    # Change to the branch you want to push to the master
3) git pull origin master     # Ensure your branch is upto date with the master - check that everything works correctly
4) git push origin dev_branch # Update your remote branch with the local status
5) git checkout master        # Change to the master branch
6) git merge dev_branch       # Merge the master with the dev_branch

### To delete a disused branch
git branch -d dev_branch # Delete locally 
git branch -d -r origin/dev_branch # Delete remote copy

# 2. Instructions to correct an unintentional git commit and git push that bring files to an older version, or unintentionally stashed changes, committed and pushed to the master

0. If you need help to identify the problem. Let other CARDAMOM developers be aware of this issue and specify the problem. 
1. Find the commit ID that caused the problems (eg. the ID looks like ‘222fb39c3d39e70904371a63b14a346cfb29db08’)
2. To be safe, we encourage you to do this fix on a new branch, then pull request to merge after all looks good. Eg. create a new branch called ‘commit_revert_demo_April2021’
3. On the branch ‘commit_revert_demo_April2021’, and run: git revert 222fb39c3d39e70904371a63b14a346cfb29db08 
FYI: This is just an example, please remember to replace the commit ID 
Git revert basically creates a new commit at the head of your current branch that is the inverse of the reverted commit. Basically a selective ‘undo’. For users, when they go to push to master next time, they will have to pull this revert commit into their tree. Basically, reverting in git is seen like any other action, and it shows up as a commit in your history. You can use git revert to ‘undo’ either a historical commit or a latest commit that you want to call back. It will not revert changes made after a historical commit.
5. On the ‘commit_revert_demo_April2021’, double check 1) if the unintentionally changes files are changed back, and 2) if the later commits after ‘222fb39c3d39e70904371a63b14a346cfb29db08’ are kept; 3) if the later commits happen to have conflicts, manually resolve the conflicts. In theory step 1) and 2) should always be successful.
6. Now that you double checked all the files are correct on branch ‘commit_revert_demo_April2021’, you can do pull request and merge this branch in to master.
7. All users, whether or not you’ve made changes after the reverted commit, save your changes with ‘git add’, then pull from the master. 
More useful tips from Alex Patton: ‘If you have made changes to your local repository and are not up to date with the origin master branch, then doing git a git pull will attempt to auto-merge the changes from master into your local version when you run git pull. If git discovers that you have made changes that are incompatible with the changes on origin master branch, then it will alert you to a merge conflict, and you will have to go into your code and sort out the conflicts yourself. Git will let you know which files it can’t merge if that happens. (Note that this paragraph is how git always works with regard to commits, and is not unique to the git revert I did)
8. As noted in the user instructions, work on your own branch is the best way to avoid dealing with this kind of problems.


# 3. Suggested protocol for collaborative development of this repository

## How to contribute to CARDAMOM

You might want to:

-	Highlight a bug in the code
-	Add an new feature/dataset
-	Improve the existing code
-	Write some documentation

First of all, make sure you are working in your own development branch. If you are working on major revisions to the code base, or collaboratively working on a particular development, then taking branches of branches may be a useful approach.
Updates to the master branch should be undertaken via a pull request. Pull requests to the master branch should generally be reserved for changes that have been tested to ensure the code still functions as expected and represent a significant problem solved (i.e. a bug-fix, new feature), rather than part-completed components.
Together “Issues” and “Pull requests” provide a place for discussing changes to the source code and provide a documented history of the changes being made, and some of the reasoning behind why they were necessary.

## Our suggested workflow for collaborators is as follows:
1)	Raise an **issue** on the github page with a description of the problem.

    a.	For bugs, it is recommended that this is done as quickly as possible, as there may be an issue others are unaware of.

    b.	For enhancements, one could raise an issue at the beginning of the code development, or work on the change independently before opening the issue when ready to submit the changes. The former is advantageous as (i) it lets people know what you are planning to work on and avoid duplicating effort, and (ii) may lead to some helpful discussions/synergies.

    Raising “issues” is a useful way to track developments being made by different researchers across different branches. There are lots of options for keeping track of things, including check-lists and markdown support for easy formatting. They could also be used to flag for help if something isn’t working. The github page can then be used as a discussion platform that can aide the development and provide a log of the changes and underlying reasons for future records. We can use the weekly CARDAMOM group calls to look over open issues each week and ensure that development needs aren’t being ignored/provide guidance for tricky issues.

2)	**Pull the latest changes** from the master branch into your development branch (or create a new branch).

3)	Make the changes **in your own branch**.

4)	**Test** the changes.

5)	To merge back into the master branch (if required), **open a pull request**. On the github page, this can be **linked to specific issues** that are resolved by the changes being pulled into the master branch. Linked issues can then be automatically closed once the pull request is approved.

6)	**Resolve conflicts** between the development branch and the master. This is often straightforward, but a good way to avoid getting stuck in more complex conflict resolution is to make sure that you try to keep your development branch up-to-date with the master branch. 

7)	(Optional) Review of pull request – may be a good idea for substantive changes to the code base.

8)	When all conflicts are resolved, **close the pull request** to merge the branches.

This protocol is a working progress. If you find it useful, then great. If there are ways you can make it better, then any suggestions are very welcome.
