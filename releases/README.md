Release Making Instructions
===========================

To create a data release, run the appropriate script, e.g.:

    python DR1.py

It will install the released version of GAS (and appropriate required other
packages) then run the gridding scripts.


Building a new release
----------------------

Checklist for creating a new release:

  1. Create a new file DR#.py in this directory and populate it
  2. Change the version of GAS in `setup.py` to `DR#`
  3. Commit the changes (make sure to add `DR#.py`)
  4. `git tag DR#` to create a tag
  5. `git push --tags` to push the tags to the github repository
