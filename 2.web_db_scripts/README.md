# Description
These scripts were used to set up the database that backs the
[PathCORE demo application](https://pathcore-demo.herokuapp.com/).
The demo application is built on the
[Flask microframework](http://flask.pocoo.org/) and deployed on
[Heroku](https://www.heroku.com/). The database is a MongoDB instance
hosted on [mLab](https://mlab.com/).

Both [Heroku](https://www.heroku.com/pricing) and
[mLab](https://mlab.com/plans/pricing/) provide free tier options for their
services.

## Step 1: mLab setup

Add a user to the new database that has write-access.
Create a credentials file (see example-mLab-credentials.yml)

## Step 2: Run `setup_pathcore_db.py`

## Step 3: Run `cache_pathcore_edge_pages.py`

## Step 4: Modify the PathCORE-demo-app source