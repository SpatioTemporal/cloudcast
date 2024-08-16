### Initial Setup

See setup.py to install the `cloudcast` package (or after adding any new modules).

```shell
    $ cd ~/myDev/cloudcast

    Mac:
        $ conda env export > environment.yml
        # Edit to remove references to things like pystare
        $ python -m pip install -e .
        $ pip install --editable .

    Linux:
        $ pip freeze > requirements.txt
        # Edit to remove references to things like pystare
        $ python setup.py build_ext --inplace
        $ pip install --editable .

    FlexFS:
        SSH into a Bayesics FlexFS server
        $ cd ~/.ssh; ssh-add <bayesics key>; cd ~/myDev/cloudcast
        $ ssh bayesics/bayesics1/bayesicsf/bayesics2
```
