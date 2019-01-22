#!/bin/bash

virtualenv ENV

source ENV/bin/activate

pip install -r requirements.txt

deactivate