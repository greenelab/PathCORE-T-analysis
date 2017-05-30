# Overview
This repository contains the scripts to run the analyses described in the
PathCORE paper. Running `./ANALYSIS.sh` is sufficient to reproduce
the results in the paper. To use PathCORE in your own analyses, please
review the sections from *The PathCORE analysis workflow* onwards in
this README.

## The `data` directory
A README is provided in the `./data` directory with details about the additional

## The `figures` directory
All figures in the PathCORE paper are also available here.

## The `jupyter-notebooks` directory
Scripts used to generate Figure 3 and Supplemental Figure 2 are provided
in notebook format. We have found that we can offer greater detail about
each of the figures in this format.

## The PathCORE analysis workflow
### Scripts:
- `run_network_creation.py`
- `run_permutation_test.py`

### Additional:
- `constants` directory
- `utils.py`

## Web application database setup
### Scripts:
- `web_initialize_db.py`
- `web_edge_page_data.py`
- `utils_setup_PAO1_example.py`

These scripts were used to set up the database that backs the
[PathCORE demo application](https://pathcore-demo.herokuapp.com/).
The demo application is built on the
[Flask microframework](http://flask.pocoo.org/) and deployed on
[Heroku](https://www.heroku.com/). The database is a MongoDB instance
hosted on [mLab](https://mlab.com/).

Both [Heroku](https://www.heroku.com/pricing) and
[mLab](https://mlab.com/plans/pricing/) provide free tier options for their
services.

### Step 1: mLab setup

Add a user to the new database that has write-access.
Create a credentials file (see example-mLab-credentials.yml)

### Step 2: Run `setup_pathcore_db.py`

### Step 3: Run `cache_pathcore_edge_pages.py`

### Step 4: Modify the PathCORE-demo-app source

